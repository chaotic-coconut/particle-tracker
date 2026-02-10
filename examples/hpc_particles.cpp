#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>  // for std::cerr
#include <format>
#include <ctime>     // for std::time_t, std::ctime
#include <chrono>
#include <set>
#include <memory>
#include <random>
#include <regex>
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <unordered_map>
#include <map>
#include <iterator>
#include <filesystem>
#include "datetime_utils.hpp"
#include "data_prep_utils.hpp"       // Provides grid, time series, splines, etc.
#include "oneapi/tbb/parallel_for.h" // TBB parallel_for
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/global_control.h"

#define PK_IO_USE_ZLIB
#include "fixed_point_core.hpp"      // the fixed-point pack/unpack
#include "fixed_point_io_codec.hpp"  // the gzip helpers: write_records_gzip/read_records_gzip

using DataType=float;
using TimeType=double;

// -----------------------------------------------------------------------------
//  Globals                                                                     |
// -----------------------------------------------------------------------------

static const TimeType life_time_seconds=2*365*24*3600.;    // 2 yr
static const TimeType grace_seconds    =7*24*3600.;    // ~7 days
static PointCloud<DataType> initial_band_cloud;
static std::unique_ptr<KDTree<DataType>> kd_initial_band;
static DataType stay_initial_band_thresh_squared;

// --------------------------------------------------------------------
// Helper: Join two sets of splines for each variable
// --------------------------------------------------------------------
template<typename DT>
void joinSplines(
    std::map<std::string,std::vector<TimeSpline<DT,TimeType>>>&& splines_1,
    std::map<std::string,std::vector<TimeSpline<DT,TimeType>>>& splines_2,
    std::map<std::string,std::vector<TimeSpline<DT,TimeType>>>& joined_splines)
{
    joined_splines=std::move(splines_1); // take ownership

    for (auto &kv : splines_2)
    {
        auto &vec=joined_splines[kv.first];
        vec.insert(vec.end(),
        std::make_move_iterator(kv.second.begin()),
        std::make_move_iterator(kv.second.end())); // move, not copy
    }
    for (auto &kv : splines_2) kv.second.shrink_to_fit();
    splines_2.clear(); // optional: drop map nodes
}

namespace std            // open the namespace, add one specialization
{
    template<>
    struct hash<Date>
    {
        std::size_t operator()(const Date& d) const noexcept
        {
            /* pack YYYY, MM, DD into one 64-bit integer */
            return static_cast<std::size_t>(d.getYear())*32*13+static_cast<std::size_t>(d.getMonth())*32+static_cast<std::size_t>(d.getDay());
        }
    };
} // namespace std

// container keyed by release month
using DayKey        =Date;
using ParticlePacket=std::vector<Particle<DataType,TimeType>>;
static std::unordered_map<DayKey,ParticlePacket> all_particles;    // new hash specialization is needed here
static std::unordered_map<DayKey,std::uint64_t> released_total;
static std::unordered_map<DayKey,std::uint64_t> landed_total;
static std::unordered_map<DayKey,std::uint64_t> landed_east;
static std::unordered_map<DayKey,std::uint64_t> landed_west;
static std::unordered_map<DayKey,std::uint64_t> boring_total;      // never left origin band (killed at grace)
static std::unordered_map<DayKey,std::uint64_t> neverland_total;   // left origin but never landed

// Random number generator
static thread_local std::mt19937 gen(std::random_device{}());

// Helper: uniform random shift in [‑max,max]
template<class DataType>
inline DataType randomShift(DataType max)
{
    std::uniform_real_distribution<DataType> dis(-max,max);
    return dis(gen);
}

// ---------------------------------------------------------------------------
// helper: zero-pad an integer to 2 characters
// ---------------------------------------------------------------------------
inline std::string formatTwo(int v)
{
    std::ostringstream oss;
    oss<<std::setw(2)<<std::setfill('0')<<v;
    return oss.str();
}

// Helper function: Convert a Date into a string formatted as "YYYY_MM_DD".
inline std::string formatDate(Date const& d)
{
    std::ostringstream oss;
    oss<<std::setw(4)<<std::setfill('0')<<d.getYear()<<"_"
        <<std::setw(2)<<std::setfill('0')<<d.getMonth()<<"_"
        <<std::setw(2)<<std::setfill('0')<<d.getDay();
    return oss.str();
}

// -----------------------------------------------------------------------------
//  Utility: parse "YYYY-MM-DD"  or "YYYY_MM_DD" into Date                    |
// -----------------------------------------------------------------------------
inline Date parseDate(const std::string& s)
{
    std::regex re(R"((\d{4})[-_](\d{2})[-_](\d{2}))");
    std::smatch m;
    if (!std::regex_match(s,m,re))
    {
        std::cerr<<"Invalid date format: "<<s<<" (expected YYYY-MM-DD)\n";
        std::exit(EXIT_FAILURE);
    }
    int y=std::stoi(m[1]),mth=std::stoi(m[2]),d=std::stoi(m[3]);
    return Date(d,mth,y);
}

struct BBox
{
    DataType lon_min{},lon_max{},lat_min{},lat_max{};
    bool wraps{};
    static inline DataType norm_pi(DataType L) noexcept
    {
        while (L<=-pi) L+=2*pi;
        while (L>  pi) L-=2*pi;
        return L;
    };
    BBox(DataType lo_min,DataType lo_max,DataType la_min,DataType la_max)
    {
        // normalize tiny numeric noise
        const DataType eps=static_cast<DataType>(1e-12);
        lon_min=norm_pi(lo_min);lon_max=norm_pi(lo_max);
        lat_min=la_min;lat_max=la_max;
        if (std::abs(lon_min-lon_max)<eps) lon_max=lon_min; // degenerate ok
        wraps=(lon_min>lon_max);
    }
    bool contains(DataType lon,DataType lat) const noexcept
    {
        if (lat<lat_min || lat>lat_max) return false;
        lon=norm_pi(lon); // make robust to any input range
        if (!wraps) return lon>=lon_min && lon<=lon_max;
        return lon>=lon_min || lon<=lon_max;
    }
};

// ---------------------------------------------------------------------------
// helper: expand pattern for one date
// ---------------------------------------------------------------------------
std::string expandPattern(const std::string& pattern,const Date& d,int level)
{
    std::ostringstream lev_ss;
    lev_ss<<std::setw(4)<<std::setfill('0')<<level;
    std::string out=pattern;

    const std::pair<std::regex,std::string> subs[]=
    {
        {std::regex("\\{YYYY\\}"),std::to_string(d.getYear())},
        {std::regex("\\{MM\\}") , formatTwo(d.getMonth())    },
        {std::regex("\\{DD\\}") , formatTwo(d.getDay())      },
        {std::regex("\\{LEV\\}"), lev_ss.str()               }
    };
    for (auto const& s : subs)
    out=std::regex_replace(out,s.first,s.second);

    return out;
}

// ---------------------------------------------------------------------------
// build path vectors for both sub-domains
// ---------------------------------------------------------------------------
void buildPathVectors(const DateRange& dr,
                      int level,
                      const std::string& path_1,
                      const std::string& path_2,
                      std::vector<std::string>& out_1,
                      std::vector<std::string>& out_2)
{
    out_1.clear();out_2.clear();
    out_1.reserve(dr.length());
    out_2.reserve(dr.length());

    for (const Date& d : dr)
    {
        out_1.push_back(expandPattern(path_1,d,level));
        out_2.push_back(expandPattern(path_2,d,level));
    }
}

// comma-separated list of args to vector<string>
inline std::vector<std::string> splitArgs(const std::string& s)
{
    std::vector<std::string> v;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss,item,','))
    {
        if (!item.empty()) v.push_back(item);
    }
    return v;
}

// -----------------------------------------------------------------------------
//  Spline containers                                                           |
// -----------------------------------------------------------------------------

template<class DT> using spline_vec=std::vector<TimeSpline<DT,TimeType>>;

template<class DT> using spline_map=std::map<std::string,spline_vec<DT>>;

struct MonthData                 // lives one month, then discarded
{
    Date month{1,1,2000};        // default so MonthData() works
    spline_map<DataType> spl_1;    // sub-domain 1  (“left half”)
    spline_map<DataType> spl_2;    // sub-domain 2  (“right half”)
    spline_map<DataType> spl;
    Grid<DataType> grid_1,grid_2;
    std::vector<std::size_t> sea_idx_1,sea_idx_2;
    std::vector<DataType> sea_x,sea_y;
    PointCloud<DataType> sea_cloud,land_cloud;
    std::unique_ptr<KDTree<DataType>> kd_sea  {nullptr};
    std::unique_ptr<KDTree<DataType>> kd_land {nullptr};
};

inline void dumpCoastCSV(const MonthData& m,const std::string& path,int stride=8)
{
    std::ofstream out(path);
    if (!out) return;
    out<<std::setprecision(10);
    // header
    out<<"lon_deg,lat_deg,label\n";
    // land points (used for beaching)
    for (size_t i=0;i+1<m.land_cloud.pts.size();i+=2*stride)
    {
        double lon_deg=m.land_cloud.pts[i]  *180./pi;
        double lat_deg=m.land_cloud.pts[i+1]*180./pi;
        out<<lon_deg<< ","<<lat_deg<<",land\n";
    }
    // (optional) sea points
    //for (size_t i=0;i+1<m.sea_cloud.pts.size();i+=2*stride)
    //{
    //    double lon_deg=m.sea_cloud.pts[i]  *180./pi;
    //    double lat_deg=m.sea_cloud.pts[i+1]*180./pi;
    //    out<<lon_deg<<","<<lat_deg<<",sea\n";
    //}
}

inline void buildKdTrees(MonthData& m)
{
    auto kd_params=nanoflann::KDTreeSingleIndexAdaptorParams(10);
    m.kd_sea =std::make_unique<KDTree<DataType>>(2,m.sea_cloud ,kd_params);
    m.kd_land=std::make_unique<KDTree<DataType>>(2,m.land_cloud,kd_params);
    m.kd_sea ->buildIndex();
    m.kd_land->buildIndex();
    assert(m.kd_sea && m.kd_land);
    assert(m.sea_cloud.kdtree_get_point_count() >0);
    assert(m.land_cloud.kdtree_get_point_count()>0);
}

MonthData loadOneMonth(
    Date month,               // 2000-01-01, 2000-02-01 ...
    int level,
    const std::string& pattern_1,
    const std::string& pattern_2,
    const std::vector<std::string>& grid_vars,
    const std::vector<std::string>& data_vars,
    TimeType dt_seconds)
{
    // ---------- 1. daily file list for *this* month (+ padding) ------------
    Date d_0=month.addDays(-1);                       // day before 1st
    Date d_1=month.addMonths(1).addDays(1);          // day after last
    DateRange span(d_0,d_1);

    std::vector<std::string> fp_1,fp_2;
    buildPathVectors(span,level,pattern_1,pattern_2,fp_1,fp_2);

    // ---------- 2. minimal grid init (first file of each sub-domain) -------
    MonthData m;m.month=month;
    m.grid_1=initializeGrid<DataType>(fp_1.front(),grid_vars);
    m.grid_2=initializeGrid<DataType>(fp_2.front(),grid_vars);

    auto norm_pi_inplace=[](std::vector<DataType>& a)
    {
        for (auto& L : a)
        {
            while (L<=-pi) L+=2*pi;
            while (L > pi) L-=2*pi;
        }
    };

    norm_pi_inplace(m.grid_1.lon_rad);
    norm_pi_inplace(m.grid_2.lon_rad);

    m.sea_idx_1=getSeaPointIndices(m.grid_1);
    m.sea_idx_2=getSeaPointIndices(m.grid_2);

    // ---------- 3. Build sea & land clouds + KD-trees for this month */
    auto buildClouds=[&](MonthData& m)
    {
        const std::size_t n_1=m.grid_1.lon_rad.size();
        const std::size_t n_2=m.grid_2.lon_rad.size();
        const std::size_t total=n_1+n_2;

        // 1) SEA CLOUD  (order must match spline vectors: spl_1 then spl_2)
        m.sea_cloud.pts.reserve((m.sea_idx_1.size()+m.sea_idx_2.size())*2);

        // domain 1
        for (auto idx : m.sea_idx_1)
        {
            m.sea_cloud.pts.push_back(m.grid_1.lon_rad[idx]);
            m.sea_cloud.pts.push_back(m.grid_1.lat_rad[idx]);
        }
        // domain 2
        for (auto idx : m.sea_idx_2)
        {
            m.sea_cloud.pts.push_back(m.grid_2.lon_rad[idx]);
            m.sea_cloud.pts.push_back(m.grid_2.lat_rad[idx]);
        }

        m.sea_x.clear();
        m.sea_y.clear();
        m.sea_x.reserve(m.sea_cloud.kdtree_get_point_count());
        m.sea_y.reserve(m.sea_cloud.kdtree_get_point_count());

        for (size_t i=0;i<m.sea_cloud.kdtree_get_point_count();i++)
        {
            m.sea_x.push_back(m.sea_cloud.kdtree_get_pt(i,0));
            m.sea_y.push_back(m.sea_cloud.kdtree_get_pt(i,1));
        }

        // 2) LAND CLOUD (anything not in sea_idx_*; order doesn't matter)
        std::vector<char> is_sea(total,0);
        for (auto i : m.sea_idx_1) is_sea[i]=1;
        for (auto i : m.sea_idx_2) is_sea[i+n_1]=1;

        m.land_cloud.pts.reserve((total-m.sea_idx_1.size()-m.sea_idx_2.size())*2);

        for (std::size_t idx=0;idx<total;idx++)
        {
            if (is_sea[idx]) continue;

            DataType lon,lat;
            if (idx<n_1)
            {
                lon=m.grid_1.lon_rad[idx];
                lat=m.grid_1.lat_rad[idx];
            }
            else
            {
                std::size_t jdx=idx-n_1;
                lon=m.grid_2.lon_rad[jdx];
                lat=m.grid_2.lat_rad[jdx];
            }
            m.land_cloud.pts.push_back(lon);
            m.land_cloud.pts.push_back(lat);
        }

        const auto sea_N =m.sea_cloud.kdtree_get_point_count();
        const auto land_N=m.land_cloud.kdtree_get_point_count();

        if (sea_N==0)
        {
            throw std::runtime_error("No SEA points for month "+formatDate(m.month)
            +" (check masks/variables/region).");
        }
        if (land_N==0)
        {
            throw std::runtime_error("No LAND points for month "+formatDate(m.month)
            +" (check masks/variables/region).");
        }

        // 3) KD trees
        //auto kd_params=nanoflann::KDTreeSingleIndexAdaptorParams(10);
        //m.kd_sea =std::make_unique<KDTree<DataType>>(2,m.sea_cloud ,kd_params);
        //m.kd_land=std::make_unique<KDTree<DataType>>(2,m.land_cloud,kd_params);
        //m.kd_sea ->buildIndex();
        //m.kd_land->buildIndex();
    };

    buildClouds(m);

    dumpCoastCSV(m,"coast_"+formatDate(m.month)+".csv",8);

    // ---------- 3. read time series + build splines ------------------------
    std::vector<TimeType> time_steps_1,time_steps_2;
    //spline_map<DataType> tmp_spl_1,tmp_spl_2;
    std::map<std::string,std::vector<std::vector<DataType>>> raw_time_series_1,raw_time_series_2;

    collectTimeSeriesData(m.grid_1,m.sea_idx_1,data_vars,fp_1,raw_time_series_1,time_steps_1,dt_seconds);
    collectTimeSeriesData(m.grid_2,m.sea_idx_2,data_vars,fp_2,raw_time_series_2,time_steps_2,dt_seconds);

    if (time_steps_1.empty() || time_steps_2.empty())
    {
        throw std::runtime_error("No time steps found for month "+formatDate(m.month)+" (check file patterns / level / variable names).");
    }

    /* fill gaps, build splines */
    auto patch=[&](auto& ts,auto& steps)
    {
        for (auto& v : data_vars)
        for (auto& vec : ts[v]) interpolateMissingData<DataType>(vec,steps);
    };
    patch(raw_time_series_1,time_steps_1);
    patch(raw_time_series_2,time_steps_2);

    createSplinesForSeaPoints(std::move(raw_time_series_1),m.spl_1,static_cast<TimeType>(time_steps_1.front()),dt_seconds);
    createSplinesForSeaPoints(std::move(raw_time_series_2),m.spl_2,static_cast<TimeType>(time_steps_2.front()),dt_seconds);

    m.spl.clear();
    joinSplines<DataType>(std::move(m.spl_1),m.spl_2,m.spl);
    m.spl_2.clear(); // we moved spl_1; spl_2 was copied into m.spl (shared_ptr),
    m.spl_2={}; // drop its map/vector overhead

    // Optional sanity check (debug builds)
    for (auto &kv : m.spl)
        assert(kv.second.size()==m.sea_cloud.kdtree_get_point_count());

    return m;
}

MonthData buf[2];static int older=0,newer=1;

// seconds since 2000-01-01 00:00
inline TimeType secondsSince2000(const Date& d)
{
    return static_cast<TimeType>(d.secondsSince(Date(1,1,2000)));
}

inline Date dateFromSecondsSince2000(TimeType sec)
{
    // brute force
    Date ref(1,1,2000);
    long long days=static_cast<long long>(std::floor(static_cast<long double>(sec)/86400.L));
    return ref.addDays(static_cast<int>(days));
}

// ============================================================================
// ONE-TIME INITIAL POSITION GENERATOR + DAILY SPAWNER
// ============================================================================

// Global cache of initial positions (reused every day)
static std::vector<point<DataType>> init_positions;

// Condition that the point is in the narrow band near the coast
inline bool inNearCoastBand(DataType lon,DataType lat,const MonthData& md,DataType sea_tol_squared,DataType land_min_squared,DataType land_max_squared)
{
    if (!md.kd_sea || !md.kd_land) throw std::runtime_error("KD trees not built");

    DataType q[2]={lon,lat};

    // nearest SEA point check
    {
        size_t idx;DataType dist_squared;
        nanoflann::KNNResultSet<DataType> rs(1);
        rs.init(&idx,&dist_squared);
        md.kd_sea->findNeighbors(rs,q,nanoflann::SearchParameters());
        if (dist_squared>sea_tol_squared) return false;
    }

    // band relative to LAND
    {
        size_t idx; DataType dist_squared;
        nanoflann::KNNResultSet<DataType> rs(1);
        rs.init(&idx, &dist_squared);
        md.kd_land->findNeighbors(rs,q,nanoflann::SearchParameters());
        if (dist_squared<land_min_squared || dist_squared>land_max_squared) return false;
    }

    return true;
}

std::vector<point<DataType>>
buildCoastalBand(const MonthData& md,const BBox& box,DataType sea_tol_squared,DataType land_min_squared,DataType land_max_squared)
{
    std::vector<point<DataType>> band;

    const std::size_t N=md.sea_cloud.kdtree_get_point_count();

    band.reserve(N/20);

    for (size_t i=0;i<N;i++)
    {
        DataType lon=md.sea_cloud.kdtree_get_pt(i,0);
        DataType lat=md.sea_cloud.kdtree_get_pt(i,1);
        if (!box.contains(lon,lat)) continue;

        if (inNearCoastBand(lon,lat,md,sea_tol_squared,land_min_squared,land_max_squared))
            band.push_back({lon,lat});
    }

    if (band.empty())
        std::cerr<<"[WARN] buildCoastalBand(): empty band\n";

    return band;
}

void buildInitialPositions(std::size_t N,const MonthData& md_ref,const BBox& box,DataType sea_tol,DataType land_min,DataType land_max,DataType jitter)
{
    if (!init_positions.empty()) return;    // already built

    const DataType sea_tol_squared=sea_tol*sea_tol;
    const DataType land_min_squared=land_min*land_min;
    const DataType land_max_squared=land_max*land_max;

    auto band=buildCoastalBand(md_ref,box,sea_tol_squared,land_min_squared,land_max_squared);
    if (band.empty())
    {
        std::cerr<<"[ERROR] no valid band points for initial seeding\n";
        return;
    }

    // --------- build coastal-band KD-tree ----------
    initial_band_cloud.pts.clear();
    initial_band_cloud.pts.reserve(band.size()*2);
    for (auto& p : band)
    {
        initial_band_cloud.pts.push_back(p[0]);
        initial_band_cloud.pts.push_back(p[1]);
    }
    auto kd_params=nanoflann::KDTreeSingleIndexAdaptorParams(10);
    kd_initial_band=std::make_unique<KDTree<DataType>>(2,initial_band_cloud,kd_params);
    kd_initial_band->buildIndex();

    // Set “still-in-band” threshold once (example: 10 km)
    DataType stay_initial_band_rad=10000./6371000.;
    stay_initial_band_thresh_squared=stay_initial_band_rad*stay_initial_band_rad;

    std::uniform_int_distribution<std::size_t> pick(0,band.size()-1);

    init_positions.reserve(N);

    std::size_t attempts=0;
    const std::size_t max_attempts=N*20;  // heuristic

    for (std::size_t accepted=0;accepted<N && attempts<max_attempts;attempts++)
    {
        auto wrap_pi=[](DataType L)
        {
            // wrap to (-pi, pi]
            while (L<=-pi) L+=2*pi;
            while (L>  pi) L-=2*pi;
            return L;
        };
        const auto& s=band[pick(gen)];
        DataType lon=wrap_pi(s[0]+randomShift(jitter));
        DataType lat=s[1]+randomShift(jitter);

        if (inNearCoastBand(lon,lat,md_ref,sea_tol_squared,land_min_squared,land_max_squared))
        {
            init_positions.push_back({lon,lat});
            ++accepted;
        }
    }

    if (init_positions.size()!=N)
    {
        std::cerr<<"Failed to generate "<<N<<" points (got "<<init_positions.size()<< "). Loosen thresholds or enlarge jitter/band.\n";
    }
    if (init_positions.empty())
    {
        throw std::runtime_error("Failed to seed initial positions (band empty or thresholds too strict).");
    }
}

void spawn(const Date& day,ParticlePacket& packet,TimeType release_time)
{
    // count once per day
    released_total[day]+=static_cast<std::uint64_t>(init_positions.size());

    packet.reserve(packet.size()+init_positions.size());
    for (auto const& xy : init_positions)
    {
        packet.push_back(Particle<DataType,TimeType>{
            xy[0],xy[1],release_time,
            xy[0],xy[1],release_time,
            EndReason::none
        });
    }
}

// ---------------------------------------------------------------
// Build a NeighborsData object for one month (uses that month’s KD-sea)
// ---------------------------------------------------------------
inline NeighborsData<DataType,TimeType>
makeNeighborsData(const MonthData& m,DataType search_radius_rad,DataType shape_param)
{
    //const size_t N=m.sea_cloud.kdtree_get_point_count();

    DataType radius_squared=search_radius_rad*search_radius_rad;
    return NeighborsData<DataType,TimeType>(*m.kd_sea,m.sea_x,m.sea_y,radius_squared,shape_param);
}

// ---------------------------------------------------------------
// Distance-to-land check (uses squared threshold)
// ---------------------------------------------------------------
inline bool nearLand(const point<DataType>& p,const MonthData& md,DataType thresh_squared)
{
    if (!md.kd_land) throw std::runtime_error("kd_land not built");
    size_t idx;
    DataType d_squared;
    nanoflann::KNNResultSet<DataType> rs(1);
    rs.init(&idx,&d_squared);
    md.kd_land->findNeighbors(rs,p.data(),nanoflann::SearchParameters());
    return d_squared<=thresh_squared;
}

inline bool stillInitialBand(const point<DataType>& p)
{
    if (!kd_initial_band) return false;  // safety
    size_t idx;
    DataType d_squared;
    nanoflann::KNNResultSet<DataType> rs(1);
    rs.init(&idx,&d_squared);
    kd_initial_band->findNeighbors(rs,p.data(),nanoflann::SearchParameters());
    return d_squared<=stay_initial_band_thresh_squared;
}

// ---------------------------------------------------------------
// Pick which month’s data to use for absolute time t
// ---------------------------------------------------------------
inline const MonthData& pickMonth(TimeType t,const MonthData& older_m,const MonthData& newer_m)
{
    Date d=dateFromSecondsSince2000(t);
    return (d.beginOfMonth()==older_m.month) ? older_m : newer_m;
}

inline NeighborsData<DataType,TimeType>& pickNeighbors(TimeType t,
                                                       const MonthData& older_m,
                                                       const MonthData& newer_m,
                                                       NeighborsData<DataType,TimeType>& nb_old,
                                                       NeighborsData<DataType,TimeType>& nb_new)
{
    Date d=dateFromSecondsSince2000(t);
    return (d.beginOfMonth()==older_m.month) ? nb_old : nb_new;
}

// ---------------------------------------------------------------
// Interpolate U,V using NeighborsData
// ---------------------------------------------------------------
//inline void interpolateUV(TimeType t,
//                          DataType lon,DataType lat,
//                          const std::vector<std::string>& vars,   // {"water_u","water_v"}
//                          const MonthData& md,
//                          NeighborsData<DataType,TimeType>& nb,
//                          DataType& u,DataType& v)
//{
//    std::vector<Neighbor<DataType>> neigh;
//    nb.computeNeighbors(neigh,lon,lat);
//    if (neigh.empty()){u=v=0;return;}
//
//    auto vals=nb.interpolateVariables(t,neigh,md.spl,vars);
//    u=vals.at(vars[0]);
//    v=vals.at(vars[1]);
//}

inline void interpolateUV(TimeType t,
                          DataType lon, DataType lat,
                          const std::vector<std::string>& vars,   // {"water_u","water_v"}
                          const MonthData& md,
                          NeighborsData<DataType,TimeType>& nb,
                          std::vector<Neighbor<DataType>>& neigh, // <-- scratch
                          DataType& u, DataType& v)
{
    neigh.clear();                                     // reuse capacity
    nb.computeNeighbors(neigh,lon,lat);
    if (neigh.empty()){u=v=0;return;}

    auto vals=nb.interpolateVariables(t,neigh,md.spl,vars);
    u=vals.at(vars[0]);
    v=vals.at(vars[1]);
}

// ---------------------------------------------------------------------------
// propagateWindow(): back-propagate all alive particles over [win_beg, win_end]
// using explicit Euler, spline interpolation & kd-trees.
// Stops on: land hit, lifetime expired, window boundary.
// ---------------------------------------------------------------------------
void propagateWindow(const Date& win_beg,
                     const Date& win_end,
                     const MonthData& older_m,
                     const MonthData& newer_m,
                     const std::vector<std::string>& data_vars,
                     TimeType dt_seconds,
                     DataType land_thresh_rad,
                     TimeType life_time_seconds)
{
    const TimeType t_beg=secondsSince2000(win_beg);      // inclusive (earliest)
    //const TimeType t_end=secondsSince2000(win_end);      // inclusive (latest)
    const TimeType dt_sec=dt_seconds;

    DataType search_radius=static_cast<DataType>(.08*1.1/180.*pi);
    DataType shape_param  =search_radius/(1.1*1.1);

    // Build neighbor structures for both buffers (no static: rebuild each window)
    auto nb_old=makeNeighborsData(older_m,search_radius,shape_param);
    auto nb_new=makeNeighborsData(newer_m,search_radius,shape_param);

    const DataType land_thresh_squared=land_thresh_rad*land_thresh_rad;

    std::atomic<uint64_t> a_past_grace{0},c_allow{0},c_beach{0};

    // Iterate through *all* release-month buckets
    for (auto& bucket : all_particles)
    {
        auto& vec=bucket.second;

        if (!older_m.kd_sea || !newer_m.kd_sea) throw std::runtime_error("NeighborsData: kd_sea missing");

        oneapi::tbb::parallel_for(
            oneapi::tbb::blocked_range<std::size_t>(0,vec.size(),4096),        //32768
            [&](const oneapi::tbb::blocked_range<std::size_t>& r)
        {
            std::vector<Neighbor<DataType>> neigh;
            neigh.reserve(8); // typical neighbor count (4/8/16/32/64/128)

            for (std::size_t i=r.begin();i!=r.end();i++)
            {
                auto& p=vec[i];
                // skip dead particles
                if (p.reason!=EndReason::none) continue;

                bool grace_done=false;

                //while (p.reason==EndReason::none && p.time>=t_beg)
                //{
                //    TimeType age=p.init_time-p.time;

                //    if (age>=life_time_seconds){p.reason=EndReason::lifetime;break;}

                //    if (!grace_done && age>=grace_seconds)
                //    {
                //        if (stillInitialBand({p.lon,p.lat})){p.reason=EndReason::neverleft;break;}
                //        grace_done=true;
                //    }

                //    const MonthData& md=pickMonth(p.time,older_m,newer_m);
                //    auto& nb           =pickNeighbors(p.time,older_m,newer_m,nb_old,nb_new);

                //    // Only beach if (past grace) OR (we’ve left the initial band)
                //    bool allow_beach=(age>=grace_seconds) || !stillInitialBand({p.lon,p.lat});
                //    if (allow_beach && nearLand({p.lon,p.lat},md,land_thresh_squared))
                //    {
                //        p.reason=EndReason::landed;
                //        break;
                //    }

                //    DataType u,v;
                //    interpolateUV(p.time,p.lon,p.lat,data_vars,md,nb,u,v);
                //    DataType dx=-u*dt_sec,dy=-v*dt_sec;
                //    inverseTransform(dx,dy,p.lon,p.lat);

                //    p.time-=dt_sec;
                //    if (p.time<t_beg) break;
                //}

                const bool started_past_grace=(p.init_time-p.time)>=grace_seconds;
                bool grace_counted=started_past_grace;  // already past -> don't count a crossing

                while (p.reason==EndReason::none && p.time>=t_beg)
                {
                    TimeType age=p.init_time-p.time;

                    if (!grace_counted && age>=grace_seconds)
                    {
                        a_past_grace.fetch_add(1,std::memory_order_relaxed);
                        grace_counted=true;
                    }

                    if (age>=life_time_seconds){p.reason=EndReason::lifetime;break;}

                    const MonthData& md=pickMonth(p.time,older_m,newer_m);
                    auto& nb=pickNeighbors(p.time,older_m,newer_m,nb_old,nb_new);

                    //[[maybe_unused]] bool in_band  =stillInitialBand({p.lon,p.lat});
                    bool past_grace=(age>=grace_seconds);
                    //bool allow_beach=(past_grace || !in_band);
                    bool allow_beach=past_grace;
                    if (allow_beach) c_allow.fetch_add(1,std::memory_order_relaxed);

                    if (allow_beach && nearLand({p.lon,p.lat},md,land_thresh_squared))
                    {
                        c_beach.fetch_add(1,std::memory_order_relaxed);
                        p.reason=EndReason::landed;
                        break;
                    }

                    // 2) (Optional) Only *mark* never-left at grace,
                    //    but don’t break here — let it keep going further back in time.
                    if (!grace_done && past_grace)
                    {
                        //if (in_band)
                        //{
                        //    // p.reason = EndReason::neverleft; break;   // <- REMOVE this early break
                        //    // If you want to count later, track a flag:
                        //    // p.flags |= EVER_IN_BAND_AT_GRACE;  (or keep a side map if you can’t change Particle)
                        //}
                        grace_done=true;
                    }

                    DataType u,v;
                    interpolateUV(p.time,p.lon,p.lat,data_vars,md,nb,neigh,u,v);
                    DataType dx=-u*dt_sec,dy=-v*dt_sec;
                    inverseTransform(dx,dy,p.lon,p.lat);

                    auto wrap_lon_inplace=[](DataType& L)
                    {
                        while (L<=-pi) L+=2*pi;
                        while (L>  pi) L-=2*pi;
                    };
                    auto clamp_lat_inplace=[](DataType& L)
                    {
                        if (L<-pi/2) L=-pi/2;
                        if (L> pi/2) L= pi/2;
                    };

                    wrap_lon_inplace(p.lon);
                    clamp_lat_inplace(p.lat);

                    p.time-=dt_sec;
                    if (p.time<t_beg) break;
                }
            }
        });
    }

    std::cerr<<"window "<<formatDate(win_beg)<<": past_grace="<<a_past_grace.load()<<" allow_beach="<<c_allow.load()<<" beached="<<c_beach.load()<<"\n";
}

inline PackedParticle packBinary(const Particle<DataType,TimeType>& p)
{
    uint8_t flags=static_cast<uint8_t>(p.reason); // 0..3

    return encode_record<DataType,TimeType>(
    p.init_lon,  // was p.lon0
    p.init_lat,  // was p.lat0
    p.init_time, // was p.t0
    p.lon,       // was p.lon1
    p.lat,       // was p.lat1
    p.time,      // was p.t1
    flags
    );
}

inline std::string dayTag(const Date& d)
{
    std::string out="trajectories_"+formatDate(d)+".pkd.gz";
    return out;
}

inline std::filesystem::path outPath(const std::filesystem::path& dir,
                                     const Date& d,
                                     int level)
{
    const auto filename=std::format("trajectories_{}_{}m.pkd.gz",formatDate(d),level);
    return dir / filename;
}

// --------------------------------------------------------------------
// Main
// --------------------------------------------------------------------
int main(int argc,char* argv[]) try
{
    if (argc!=12 && argc!=13)
    {
        std::cerr<<"usage: "<<argv[0]
        <<" REL_START REL_END LON_MIN LON_MAX LAT_MIN LAT_MAX LEVEL TIMESTEP PATTERN1 PATTERN2 OUT_DIR"
        <<" [[GRID_LON,GRID_LAT,DATA1,DATA2,...]]\n";
        return 1;
    }

    //counter on evaluation time
    auto start=std::chrono::system_clock::now();
    std::time_t start_time=std::chrono::system_clock::to_time_t(start);
    std::cerr<<"started computation at "<<std::ctime(&start_time)<<'\n';

    //oneapi::tbb::global_control control(oneapi::tbb::global_control::max_allowed_parallelism,64); // Limit to XX threads

    /* release window ----------------------------------------------------- */
    Date start_release=parseDate(argv[1]);
    Date end_release  =parseDate(argv[2]);
    //DateRange sim_range(start_release,end_release);

    /* spatial filter ----------------------------------------------------- */
    DataType lon_min=std::stof(argv[3])*pi/180.;
    DataType lon_max=std::stof(argv[4])*pi/180.;
    DataType lat_min=std::stof(argv[5])*pi/180.;
    DataType lat_max=std::stof(argv[6])*pi/180.;

    int level=std::stoi(argv[7]);
    int timestep=std::stoi(argv[8]);
    std::string pattern_1=argv[9];
    std::string pattern_2=argv[10];
    std::filesystem::path out_dir=argv[11];
    std::error_code ec;
    std::filesystem::create_directories(out_dir,ec);
    if (ec)
    {
        std::cerr<<"Cannot create output dir '"<<out_dir<<"': "<< ec.message()<<"\n";
        return 1;
    }

    // ------------------------------------------------------------------
    // variable names (default to water_u / water_v)
    // ------------------------------------------------------------------
    std::vector<std::string> data_variable_names={"water_u","water_v"};
    std::vector<std::string> grid_variable_names={"lon","lat","water_u"};

    if (argc==13)
    {
        auto names=splitArgs(argv[12]);           // lon,lat,data1,...

        if (names.size()<3)
        {
            std::cerr<<"ERROR: variable list must contain lon_var,lat_var,at least one data_var\n";
            return 1;
        }

        // replace grids
        grid_variable_names[0]=names[0];         // lon variable
        grid_variable_names[1]=names[1];         // lat variable
        grid_variable_names[2]=names[2];         // first data-var (for NaN mask)

        // replace data list
        data_variable_names.assign(names.begin()+2,names.end());
    }

    //TimeType dt_hours=timestep;           // in hours
    TimeType dt_seconds=static_cast<TimeType>(timestep)*3600.;
    Date cur_month=end_release.beginOfMonth();          // e.g. 2004-01-01


    // Basic sanity on CLI angles (degrees expected, convert to rad above)
    auto in_range=[](double x,double a,double b){return x>=a && x<=b;};

    if (!in_range(std::stod(argv[3]),-180.,180.) ||
        !in_range(std::stod(argv[4]),-180.,180.) ||
        !in_range(std::stod(argv[5]), -90., 90.) ||
        !in_range(std::stod(argv[6]), -90., 90.))
    {
        std::cerr<<"Longitude/latitude degrees out of range. "
            << "Use lon in [-180,180], lat in [-90,90].\n";
        return 1;
    }

    if (timestep<=0)
    {
        std::cerr<<"TIMESTEP must be positive (hours).\n";
        return 1;
    }

    // load Month 0 and Month -1 (+padding handled inside)
    buf[older]=loadOneMonth(
        cur_month,
        level,
        pattern_1,pattern_2,
        grid_variable_names,
        data_variable_names,
        dt_seconds);
    buildKdTrees(buf[older]);

    buf[newer]=loadOneMonth(
        cur_month.addMonths(-1),
        level,
        pattern_1,pattern_2,
        grid_variable_names,
        data_variable_names,
        dt_seconds);
    buildKdTrees(buf[newer]);

    /* ------------------------------------------------------------------------- */
    /*  MAIN MONTH-ADVANCING LOOP                                                */
    /* ------------------------------------------------------------------------- */

    BBox bbox{lon_min,lon_max,lat_min,lat_max};

    DataType sea_tol =.2*pi/180.;
    DataType land_min=2000. /6371000.;
    DataType land_max=50000./6371000.;
    DataType jitter =.0005;
    std::size_t N_per_day=2'800'000;    //2'000'000

    buildInitialPositions(N_per_day,buf[older],bbox,sea_tol,land_min,land_max,jitter);

    //const Date last_month=end_release.beginOfMonth();           // e.g. 2004-12-01
    const Date first_month=start_release.addMonths(-24).beginOfMonth();  // earliest month

    while (cur_month>=first_month)
    {
        /* 1. one-month simulation window [cur_month , cur_month+1) ------------- */
        Date win_beg=cur_month;
        Date win_end=cur_month.addMonths(1).addDays(-1);     // inclusive last day

        for (Date day=win_beg;day<=win_end;day.increment())
        {
            if (day<start_release || day>end_release) continue;  // <-- guard
            TimeType rel_time=secondsSince2000(day);

            //MonthKey key=day.beginOfMonth();
            DayKey key=day;
            ParticlePacket& pack=all_particles[key];

            spawn(day,pack,rel_time);
        }

        DataType land_thresh=(0.5)*pi/180.;//5000./6371000.;
        propagateWindow(win_beg,win_end,buf[older],buf[newer],
            data_variable_names,
            dt_seconds,
            land_thresh,
            life_time_seconds);

        for (auto it=all_particles.begin();it!=all_particles.end();)
        {
                const DayKey day=it->first;        // copy key (iterator may change)
                auto&        vec=it->second;

            std::vector<PackedParticle> recs;
            recs.reserve(vec.size());

            vec.erase(std::remove_if(vec.begin(),vec.end(),[&](const auto& p)
            {
                if (p.reason==EndReason::none) return false;

                switch (p.reason)
                {
                    case EndReason::landed:
                    {
                        recs.push_back(packBinary(p));
                        landed_total[day]++;

                        // classify by longitude at landing
                        double lon_deg=p.lon*180./pi;
                        double lon360 =std::fmod(lon_deg+360.,360.); // [0,360)

                        if (lon360>=180.) landed_west[day]++;
                        else              landed_east[day]++;
                        break;
                    }

                    case EndReason::neverleft:
                        boring_total[day]++;
                        break;

                    case EndReason::lifetime:
                        neverland_total[day]++;
                        break;

                    default:
                        break;
                }
                return true; // drop ended from RAM
            }),
            vec.end());

            if (!recs.empty())
            {
                const auto path=outPath(out_dir,day,level);
                std::error_code ec2;
                const bool append=std::filesystem::exists(path,ec2) && !ec2;
                write_records_gzip(path.string(),recs,released_total[day],/*level=*/6,/*append=*/append);
                // optional: std::cerr<<"write "<<recs.size()<<" to "<<path<<(append ? " (append)\n" : " (new)\n");
            }

            if (vec.empty()) it=all_particles.erase(it);
            else             ++it;
        }

        std::uint64_t L=0,B=0,N=0,R=0,E=0,W=0;
        for (Date day=win_beg;day<=win_end;day.increment())
        {
            L+=landed_total[day];
            B+=boring_total[day];
            N+=neverland_total[day];
            R+=released_total[day];
            E+=landed_east[day];
            W+=landed_west[day];
        }
        std::cerr<<"["<<formatDate(win_beg)<<"–"<<formatDate(win_end)<<"] "<<"released="<<R<<" landed="<<L<<" landed east ="<<E<<" landed west="<<W<<
        '\n'<<"neverleft="<<B<<" lifetime="<<N<<"\n";


        /* 2. advance the ring buffer by one month ---------------------------- */
        cur_month=cur_month.addMonths(-1);     // Feb → Mar → Apr …

        older^=1;                           // swap indices (0 ↔ 1)
        newer^=1;

        /* load the *new* Month+1 (two months ahead of 'older') */
        if (cur_month>first_month)
        {
            buf[newer]=MonthData{}; // drop old newer before loading fresh
            buf[newer]=loadOneMonth(
                cur_month.addMonths(-1),
                level,
                pattern_1,pattern_2,
                grid_variable_names,
                data_variable_names,
                dt_seconds);
            buildKdTrees(buf[newer]);
        }
    }

    for (auto& [day,vec] : all_particles)
        for (auto& p : vec)
            if (p.reason==EndReason::none)
                neverland_total[day]++;

    all_particles.clear();

    auto end=std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds=end-start;
    std::time_t end_time=std::chrono::system_clock::to_time_t(end);
    std::cerr<<'\n'<<"finished computation at "<<std::ctime(&end_time)<<"elapsed time: "<<elapsed_seconds.count()<<"s\n";

    return 0;
}
catch (const std::exception& e)
{
    std::cerr<<"FATAL: "<< e.what()<<"\n";
    return 2;
}
catch (...)
{
    std::cerr<<"FATAL: unknown exception\n";
    return 3;
}
