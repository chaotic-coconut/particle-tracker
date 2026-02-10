// hpc_trajectories.cpp

#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <format>
#include <ctime>
#include <chrono>
#include <set>
#include <memory>
#include <random>
#include <limits>
#include <regex>
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <unordered_map>
#include <filesystem>
#include <string_view>
#include <map>
#include <vector>
#include <iterator>
#include <cctype>
#include <stdexcept>

#include "datetime_utils.hpp"
#include "data_prep_utils.hpp"
#include "seed_loader_v2.hpp"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/global_control.h"
#include "oneapi/tbb/enumerable_thread_specific.h"

#define PK_IO_USE_ZLIB
#include "fixed_point_core.hpp"
#include "fixed_point_io_codec.hpp"
#include "fixed_point_pkd2.hpp"

using DataType=float;
using TimeType=double;

static const TimeType life_limit_seconds=2*365*24*3600.;

// --------------------------------------------------------------------
// Helper: Join two sets of splines for each variable
// --------------------------------------------------------------------
template<typename DT>
void joinSplines(
    std::map<std::string,std::vector<TimeSpline<DT,TimeType>>>&& splines_1,
    std::map<std::string,std::vector<TimeSpline<DT,TimeType>>>& splines_2,
    std::map<std::string,std::vector<TimeSpline<DT,TimeType>>>& joined_splines)
{
    joined_splines=std::move(splines_1);
    for (auto &kv : splines_2)
    {
        auto &vec=joined_splines[kv.first];
        vec.insert(vec.end(),
            std::make_move_iterator(kv.second.begin()),
            std::make_move_iterator(kv.second.end()));
    }
    for (auto &kv : splines_2) kv.second.shrink_to_fit();
    splines_2.clear();
}

namespace std {
    template<> struct hash<Date>
    {
        std::size_t operator()(const Date& d) const noexcept
        {
            return static_cast<std::size_t>(d.getYear())*32*13
                 + static_cast<std::size_t>(d.getMonth())*32
                 + static_cast<std::size_t>(d.getDay());
        }
    };
}

// ---------------- small helpers ----------------
static inline std::string lower(std::string s){
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return s;
}
static inline bool ends_with_ci(const std::string& s, const std::string& suf)
{
    if (s.size() < suf.size()) return false;
    return lower(s.substr(s.size() - suf.size())) == lower(suf);
}
static inline SeedPoint opposite_seed_point(SeedPoint p)
{
    return (p==SeedPoint::First) ? SeedPoint::Last : SeedPoint::First;
}

// --------------------------------------------------------------------
// Trajectory storage
// --------------------------------------------------------------------
struct Trajectory
{
    uint64_t seed_id=0;
    std::vector<float> lon;
    std::vector<float> lat;
    std::vector<float> u;
    std::vector<float> v;
    int32_t first_idx=std::numeric_limits<int32_t>::max();

    inline void append(float lon_r,float lat_r,float u_ms,float v_ms,int32_t idx)
    {
        lon.push_back(lon_r);
        lat.push_back(lat_r);
        u.push_back(u_ms);
        v.push_back(v_ms);
        if (idx<first_idx) first_idx=idx;
    }

    inline void flush_to(pkd2::Writer& w_pos,pkd2::Writer& w_vel,uint32_t stride=1)
    {
        if (lon.empty()) return;
        std::reverse(lon.begin(),lon.end());
        std::reverse(lat.begin(),lat.end());
        std::reverse(u.begin(),u.end());
        std::reverse(v.begin(),v.end());
        // Positions
        w_pos.add_traj_rad<float>(
            seed_id,
            lon.data(),lat.data(),
            static_cast<uint32_t>(lon.size()),
            first_idx,stride
        );
        // Velocities
        // lon channel = u, lat channel = v (units: m/s).
        w_vel.add_traj_rad<float>(
            seed_id,
            u.data(),v.data(),
            static_cast<uint32_t>(u.size()),
            first_idx,stride
        );
        lon.clear();lat.clear();u.clear();v.clear();
        first_idx=std::numeric_limits<int32_t>::max();
    }
};

// --------------------------------------------------------------------
// Particle: integer time indexing + optional stop target
// --------------------------------------------------------------------
struct ParticleState
{
    uint64_t seed_id=0;
    DataType lon=0,lat=0;

    int64_t  init_idx=0;     // start index on dt grid
    int64_t  time_idx=0;     // current index

    bool     has_stop=false;
    int64_t  stop_idx=0;     // valid only if has_stop==true
    DataType stop_lon=0, stop_lat=0;

    EndReason reason=EndReason::none;
    Trajectory traj;
};

using DayKey=Date;
using ParticlePacket=std::vector<ParticleState>;
static std::unordered_map<DayKey,ParticlePacket> all_particles;
static std::unordered_map<DayKey,std::uint64_t> released_total;

// --------------------------------------------------------------------
// Formatting helpers
// --------------------------------------------------------------------
inline std::string formatTwo(int v)
{
    std::ostringstream oss;
    oss<<std::setw(2)<<std::setfill('0')<<v;
    return oss.str();
}

inline std::string formatDate(Date const& d)
{
    std::ostringstream oss;
    oss<<std::setw(4)<<std::setfill('0')<<d.getYear()<<"_"
       <<std::setw(2)<<std::setfill('0')<<d.getMonth()<<"_"
       <<std::setw(2)<<std::setfill('0')<<d.getDay();
    return oss.str();
}

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

inline std::vector<std::string> splitArgs(const std::string& s)
{
    std::vector<std::string> v;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss,item,',')) if (!item.empty()) v.push_back(item);
    return v;
}

template<class DT> using spline_vec=std::vector<TimeSpline<DT,TimeType>>;
template<class DT> using spline_map=std::map<std::string,spline_vec<DT>>;

struct MonthData
{
    Date month{1,1,2000};
    spline_map<DataType> spl_1;
    spline_map<DataType> spl_2;
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
    out<<"lon_deg,lat_deg,label\n";
    for (size_t i=0;i+1<m.land_cloud.pts.size();i+=2*stride)
    {
        double lon_deg=m.land_cloud.pts[i]  *180./pi;
        double lat_deg=m.land_cloud.pts[i+1]*180./pi;
        out<<lon_deg<< ","<<lat_deg<<",land\n";
    }
}

inline void buildKdTrees(MonthData& m)
{
    auto kd_params=nanoflann::KDTreeSingleIndexAdaptorParams(10);
    m.kd_sea =std::make_unique<KDTree<DataType>>(2,m.sea_cloud ,kd_params);
    m.kd_land=std::make_unique<KDTree<DataType>>(2,m.land_cloud,kd_params);
    m.kd_sea ->buildIndex();
    m.kd_land->buildIndex();
}

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

// seconds since 2000-01-01 00:00
inline TimeType secondsSince2000(const Date& d)
{
    return static_cast<TimeType>(d.secondsSince(Date(1,1,2000)));
}

inline Date dateFromSecondsSince2000(TimeType sec)
{
    Date ref(1,1,2000);
    long long days=static_cast<long long>(std::floor(static_cast<long double>(sec)/86400.L));
    return ref.addDays(static_cast<int>(days));
}

// --------------------------------------------------------------------
// Epoch + snapping utilities
// --------------------------------------------------------------------
enum class SeedTimeEpoch {Unix,Since2000};
constexpr int64_t UNIX_TO_2000=946684800;

static inline int64_t snap_to_grid(int64_t t_s,int32_t dt_s)
{
    const double x=static_cast<double>(t_s)/static_cast<double>(dt_s);
    return static_cast<int64_t>(std::llround(x))*static_cast<int64_t>(dt_s);
}

static inline int64_t sec_to_idx_exact(int64_t t_s, int32_t dt_s)
{
    if (dt_s<=0) throw std::runtime_error("dt_s <= 0");
    if (t_s % dt_s == 0) return t_s / dt_s;
    return static_cast<int64_t>(std::llround(static_cast<double>(t_s)/static_cast<double>(dt_s)));
}

inline void normalize_seed_times_to_t2000(std::vector<SeedInit>& seeds,SeedTimeEpoch epoch)
{
    if (epoch==SeedTimeEpoch::Since2000) return;
    for (auto& s : seeds)
    {
        s.t_start_s -= UNIX_TO_2000;
        if (s.has_stop) s.t_stop_s -= UNIX_TO_2000;
    }
}

// IMPORTANT: do NOT “fix” ordering by swapping only times.
// We only snap times here.
inline void snap_seed_times_to_grid(std::vector<SeedInit>& seeds,int32_t dt_s)
{
    for (auto& s : seeds)
    {
        s.t_start_s = snap_to_grid(s.t_start_s, dt_s);
        if (s.has_stop) s.t_stop_s = snap_to_grid(s.t_stop_s, dt_s);
    }
}

// --------------------------------------------------------------------
// Enqueue seeds (sets integer indices and stop target)
// --------------------------------------------------------------------
static inline void enqueue_seeds(const std::vector<SeedInit>& seeds, int32_t dt_s)
{
    for (const auto& s : seeds)
    {
        const Date d=dateFromSecondsSince2000(static_cast<TimeType>(s.t_start_s));
        auto& vec=all_particles[d];

        ParticleState ps;
        ps.seed_id   = s.seed_id;

        // Start point = lon_start_rad/lat_start_rad by our record-building logic in main()
        ps.lon       = static_cast<DataType>(s.lon_start_rad);
        ps.lat       = static_cast<DataType>(s.lat_start_rad);

        ps.init_idx  = sec_to_idx_exact(s.t_start_s, dt_s);
        ps.time_idx  = ps.init_idx;

        ps.has_stop  = s.has_stop;
        if (s.has_stop)
        {
            ps.stop_lon = static_cast<DataType>(s.lon_stop_rad);
            ps.stop_lat = static_cast<DataType>(s.lat_stop_rad);
            ps.stop_idx = sec_to_idx_exact(s.t_stop_s, dt_s);
        }
        else
        {
            ps.stop_idx = 0; // unused
        }

        ps.reason    = EndReason::none;
        ps.traj.seed_id = s.seed_id;

        vec.push_back(std::move(ps));
        released_total[d]++;
    }
}

// Map month-start Date -> seeds
using SeedsByMonth=std::unordered_map<Date,std::vector<SeedInit>>;

static inline Date month_start_from_t2000(int64_t t2000_s)
{
    return dateFromSecondsSince2000(static_cast<TimeType>(t2000_s)).beginOfMonth();
}

static inline void bucketize_seeds(const std::vector<SeedInit>& seeds,SeedsByMonth& bym)
{
    bym.reserve(seeds.size()/32+8);
    for (const auto& s : seeds)
        bym[month_start_from_t2000(s.t_start_s)].push_back(s);
}

static inline void activate_month(SeedsByMonth& bym,const Date& win_beg,const Date& win_end,int32_t dt_s)
{
    const Date key=win_beg.beginOfMonth();
    auto it=bym.find(key);
    if (it==bym.end()) return;

    const int64_t tmin=static_cast<int64_t>(secondsSince2000(win_beg));
    const int64_t tmax=static_cast<int64_t>(secondsSince2000(win_end.addDays(1)));

    std::vector<SeedInit> window_seeds;
    window_seeds.reserve(it->second.size());
    for (const auto& s : it->second)
        if (s.t_start_s>=tmin && s.t_start_s<tmax)
            window_seeds.push_back(s);

    if (!window_seeds.empty())
        enqueue_seeds(window_seeds, dt_s);

    bym.erase(it);
}

inline NeighborsData<DataType,TimeType>
makeNeighborsData(const MonthData& m,DataType search_radius_rad,DataType shape_param)
{
    DataType radius_squared=search_radius_rad*search_radius_rad;
    return NeighborsData<DataType,TimeType>(*m.kd_sea,m.sea_x,m.sea_y,radius_squared,shape_param);
}

struct ReconGuards
{
    DataType max_current_mps     =3.;
    int      max_empty_neighbors =10;
    DataType step_slack          =3.;
};

void propagateWindow(const Date& win_beg,
                     const Date& win_end,
                     const MonthData& older_m,
                     const MonthData& newer_m,
                     const std::vector<std::string>& data_vars,
                     TimeType dt_seconds,
                     TimeType life_time_seconds,
                     pkd2::Writer& pkd_pos,
                     pkd2::Writer& pkd_vel,
                     int32_t base_dt_s,
                     const ReconGuards& G=ReconGuards{})
{
    (void)win_end;
    (void)dt_seconds;

    if (data_vars.size()<2)
        throw std::runtime_error("propagateWindow: need at least 2 data_vars (u,v)");
    if (!older_m.kd_sea || !newer_m.kd_sea)
        throw std::runtime_error("propagateWindow: kd_sea missing");
    if (base_dt_s<=0)
        throw std::runtime_error("propagateWindow: base_dt_s must be positive");

    const TimeType dtT  =static_cast<TimeType>(base_dt_s);

    // Window lower bound in *index* (exact: win_beg is midnight, base_dt_s divides 86400 for 3h)
    const int64_t t_beg_s   = static_cast<int64_t>(secondsSince2000(win_beg));
    const int64_t t_beg_idx = t_beg_s / static_cast<int64_t>(base_dt_s);

    const DataType search_radius=static_cast<DataType>(.08*1.1/180.*pi);
    const DataType shape_param  =search_radius/(1.1*1.1);

    auto nb_old=makeNeighborsData(older_m,search_radius,shape_param);
    auto nb_new=makeNeighborsData(newer_m,search_radius,shape_param);

    std::atomic<uint64_t> stop_life{0},stop_no_nb{0},stop_fast{0},stop_step{0}, stop_seed0{0};

    oneapi::tbb::enumerable_thread_specific<std::vector<Trajectory>> done_tls;

    for (auto& bucket : all_particles)
    {
        auto& vec=bucket.second;

        oneapi::tbb::parallel_for(
            oneapi::tbb::blocked_range<std::size_t>(0,vec.size(),4096),
            [&](const oneapi::tbb::blocked_range<std::size_t>& r)
        {
            std::vector<Neighbor<DataType>> neigh;
            neigh.reserve(16);

            for (std::size_t i=r.begin();i!=r.end();i++)
            {
                auto& p=vec[i];
                if (p.reason!=EndReason::none) continue;

                p.traj.seed_id=p.seed_id;

                int empty_nb_streak=0;
                int fast_streak    =0;

                // Stop if we hit window bound OR stop target (seed0) OR lifetime
                while (p.reason==EndReason::none &&
                       p.time_idx>=t_beg_idx &&
                       (!p.has_stop || p.time_idx>=p.stop_idx))
                {
                    // If exactly at seed0 time, snap to seed0 position, record, and finish.
                    if (p.has_stop && p.time_idx==p.stop_idx)
                    {
                        p.lon = p.stop_lon;
                        p.lat = p.stop_lat;

                        if (p.time_idx<std::numeric_limits<int32_t>::min() || p.time_idx>std::numeric_limits<int32_t>::max())
                        {
                            p.reason=EndReason::lifetime;
                            ++stop_life;
                            break;
                        }

                        p.traj.append(static_cast<float>(p.lon),
                                      static_cast<float>(p.lat),
                                      std::numeric_limits<float>::quiet_NaN(),
                                      std::numeric_limits<float>::quiet_NaN(),
                                      static_cast<int32_t>(p.time_idx));

                        p.reason=EndReason::lifetime; // finished
                        ++stop_seed0;
                        break;
                    }

                    if (p.time_idx<std::numeric_limits<int32_t>::min() || p.time_idx>std::numeric_limits<int32_t>::max())
                    {
                        p.reason=EndReason::lifetime;
                        ++stop_life;
                        break;
                    }

                    // exact lifetime in index space
                    const int64_t age_idx = (p.init_idx - p.time_idx);
                    const TimeType age_s  = static_cast<TimeType>(age_idx) * dtT;
                    if (age_s>=life_time_seconds)
                    {
                        p.reason=EndReason::lifetime;
                        ++stop_life;
                        break;
                    }

                    // absolute time for interpolation
                    const TimeType t_s = static_cast<TimeType>(p.time_idx) * dtT;

                    const Date dcur=dateFromSecondsSince2000(t_s);
                    const bool use_older=(dcur.beginOfMonth()==older_m.month);
                    const MonthData& md =use_older ? older_m : newer_m;
                    auto& nb            =use_older ? nb_old  : nb_new;

                    neigh.clear();
                    nb.computeNeighbors(neigh,p.lon,p.lat);

                    if (neigh.empty())
                    {
                        // record position, but velocity unknown
                        p.traj.append(
                            static_cast<float>(p.lon),
                            static_cast<float>(p.lat),
                            std::numeric_limits<float>::quiet_NaN(),
                            std::numeric_limits<float>::quiet_NaN(),
                            static_cast<int32_t>(p.time_idx)
                        );
                    
                        if (++empty_nb_streak>G.max_empty_neighbors)
                        {
                            p.reason=EndReason::lifetime;
                            stop_no_nb++;
                            break;
                        }
                        --p.time_idx;
                        continue;
                    }
                    empty_nb_streak=0;

                    auto vals=nb.interpolateVariables(t_s,neigh,md.spl,data_vars);
                    DataType u=vals.at(data_vars[0]);
                    DataType v=vals.at(data_vars[1]);

                    p.traj.append(
                        static_cast<float>(p.lon),
                        static_cast<float>(p.lat),
                        static_cast<float>(u),
                        static_cast<float>(v),
                        static_cast<int32_t>(p.time_idx)
                    );

                    const double spd=std::sqrt(double(u)*double(u)+double(v)*double(v));
                    if (spd>G.max_current_mps)
                    {
                        if (++fast_streak>4)
                        {
                            p.reason=EndReason::lifetime;
                            ++stop_fast;
                            break;
                        }
                    }
                    else fast_streak=0;

                    const DataType dt=static_cast<DataType>(base_dt_s);
                    DataType dx=-u*dt;
                    DataType dy=-v*dt;

                    {
                        const double max_step_km=G.max_current_mps*double(base_dt_s)*G.step_slack/1000.;
                        const double dr_km=std::sqrt(double(dx)*double(dx)+double(dy)*double(dy))/1000.;
                        if (dr_km>max_step_km)
                        {
                            p.reason=EndReason::lifetime;
                            ++stop_step;
                            break;
                        }
                    }

                    inverseTransform(dx,dy,p.lon,p.lat);

                    while (p.lon<=-pi) p.lon+=2*pi;
                    while (p.lon>  pi) p.lon-=2*pi;
                    if (p.lat<-pi/2) p.lat=-pi/2;
                    if (p.lat> pi/2) p.lat= pi/2;

                    --p.time_idx;
                }

                if (p.reason!=EndReason::none)
                {
                    done_tls.local().push_back(std::move(p.traj));
                    p.traj=Trajectory{};
                    p.traj.seed_id=p.seed_id;
                }
            }
        });
    }

    // flush finished trajectories serially (Writer not thread-safe)
    for (auto& local : done_tls)
    {
        for (auto& tr : local)
            tr.flush_to(pkd_pos,pkd_vel,/*stride=*/1);
        local.clear();
    }

    std::cerr<<std::format(
        "window {} | lifetime={} no_nb={} fast={} step={} seed0={}\n",
        formatDate(win_beg),
        stop_life.load(),
        stop_no_nb.load(),
        stop_fast.load(),
        stop_step.load(),
        stop_seed0.load()
    );
}

// ---------------- loadOneMonth (unchanged logic) ----------------
MonthData loadOneMonth(
    Date month,
    int level,
    const std::string& pattern_1,
    const std::string& pattern_2,
    const std::vector<std::string>& grid_vars,
    const std::vector<std::string>& data_vars,
    TimeType dt_seconds)
{
    Date d_0=month.addDays(-1);
    Date d_1=month.addMonths(1).addDays(1);
    DateRange span(d_0,d_1);

    std::vector<std::string> fp_1,fp_2;
    buildPathVectors(span,level,pattern_1,pattern_2,fp_1,fp_2);

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

    auto buildClouds=[&](MonthData& m)
    {
        const std::size_t n_1=m.grid_1.lon_rad.size();
        const std::size_t n_2=m.grid_2.lon_rad.size();
        const std::size_t total=n_1+n_2;

        m.sea_cloud.pts.reserve((m.sea_idx_1.size()+m.sea_idx_2.size())*2);

        for (auto idx : m.sea_idx_1)
        {
            m.sea_cloud.pts.push_back(m.grid_1.lon_rad[idx]);
            m.sea_cloud.pts.push_back(m.grid_1.lat_rad[idx]);
        }
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

        if (sea_N==0) throw std::runtime_error("No SEA points for month "+formatDate(m.month));
        if (land_N==0) throw std::runtime_error("No LAND points for month "+formatDate(m.month));
    };

    buildClouds(m);
    dumpCoastCSV(m,"coast_"+formatDate(m.month)+".csv",8);

    std::vector<TimeType> time_steps_1,time_steps_2;
    std::map<std::string,std::vector<std::vector<DataType>>> raw_time_series_1,raw_time_series_2;

    collectTimeSeriesData(m.grid_1,m.sea_idx_1,data_vars,fp_1,raw_time_series_1,time_steps_1,dt_seconds);
    collectTimeSeriesData(m.grid_2,m.sea_idx_2,data_vars,fp_2,raw_time_series_2,time_steps_2,dt_seconds);

    if (time_steps_1.empty() || time_steps_2.empty())
        throw std::runtime_error("No time steps found for month "+formatDate(m.month));

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
    m.spl_2.clear();
    m.spl_2={};

    for (auto &kv : m.spl)
        assert(kv.second.size()==m.sea_cloud.kdtree_get_point_count());

    return m;
}

MonthData buf[2]; static int older=0,newer=1;

// --------------------------------------------------------------------
// Main
// --------------------------------------------------------------------
int main(int argc, char* argv[]) try
{
    const auto wall_start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(wall_start);
    std::cerr << "started computation at " << std::ctime(&start_time) << '\n';

    if (argc < 9)
    {
        std::cerr
          << "usage: " << argv[0] << "\n"
          << "  REL_START REL_END LEVEL TIMESTEP_H PATTERN1 PATTERN2 OUT_DIR SEEDS_PATH\n"
          << " [SEED_EPOCH] [VARLIST]\n"
          << " [--keep-p p] [--cap N] [--seed S] [--seed-point first|last]\n";
        return 1;
    }

    const Date start_release = parseDate(argv[1]);
    const Date end_release   = parseDate(argv[2]);
    const int  level         = std::stoi(argv[3]);
    const int  timestep_h    = std::stoi(argv[4]);
    const std::string pattern_1 = argv[5];
    const std::string pattern_2 = argv[6];
    const std::filesystem::path out_dir = argv[7];
    const std::string seeds_path = argv[8];

    if (timestep_h <= 0) { std::cerr << "TIMESTEP_H must be positive.\n"; return 1; }
    const TimeType dt_seconds = static_cast<TimeType>(timestep_h) * 3600.0;
    const int32_t base_dt_s   = static_cast<int32_t>(dt_seconds);

    int argi = 9;

    SeedTimeEpoch epoch = ends_with_ci(seeds_path, ".pkd2") ? SeedTimeEpoch::Since2000 : SeedTimeEpoch::Unix;

    std::vector<std::string> data_variable_names = {"water_u","water_v"};
    std::vector<std::string> grid_variable_names = {"lon","lat","water_u"};

    double      keep_p        = 1.0;
    std::size_t cap           = std::numeric_limits<std::size_t>::max();
    uint32_t    thin_seed     = 0;
    bool        use_thin_seed = false;

    SeedPoint seed_point = SeedPoint::Last; // backtracking start point

    auto is_flag = [](const char* s)->bool {
        return s && std::string_view(s).rfind("--",0)==0;
    };

    if (argi < argc && !is_flag(argv[argi]))
    {
        const std::string tok = lower(argv[argi]);
        if (tok=="unix" || tok=="t2000" || tok=="since2000")
        {
            epoch = (tok=="unix") ? SeedTimeEpoch::Unix : SeedTimeEpoch::Since2000;
            ++argi;
        }
    }

    if (argi < argc && !is_flag(argv[argi]))
    {
        const std::string tok = argv[argi];
        if (tok.find(',') != std::string::npos)
        {
            auto names = splitArgs(tok);
            if (names.size() < 3)
            {
                std::cerr << "ERROR: variable list must contain lon_var,lat_var,at least one data_var\n";
                return 1;
            }
            grid_variable_names[0] = names[0];
            grid_variable_names[1] = names[1];
            grid_variable_names[2] = names[2];
            data_variable_names.assign(names.begin()+2, names.end());
            ++argi;
        }
    }

    for (; argi < argc; ++argi)
    {
        const std::string a = argv[argi];
        if (a=="--keep-p")
        {
            if (argi+1>=argc) { std::cerr<<"ERROR: --keep-p requires a value\n"; return 1; }
            keep_p = std::stod(argv[++argi]);
        }
        else if (a=="--cap")
        {
            if (argi+1>=argc) { std::cerr<<"ERROR: --cap requires a value\n"; return 1; }
            cap = static_cast<std::size_t>(std::stoull(argv[++argi]));
        }
        else if (a=="--seed")
        {
            if (argi+1>=argc) { std::cerr<<"ERROR: --seed requires a value\n"; return 1; }
            thin_seed = static_cast<uint32_t>(std::stoul(argv[++argi]));
            use_thin_seed = true;
        }
        else if (a=="--seed-point")
        {
            if (argi+1>=argc) { std::cerr<<"ERROR: --seed-point requires first|last\n"; return 1; }
            const std::string v = lower(argv[++argi]);
            if (v=="first") seed_point = SeedPoint::First;
            else if (v=="last") seed_point = SeedPoint::Last;
            else { std::cerr<<"ERROR: --seed-point must be first|last\n"; return 1; }
        }
        else
        {
            std::cerr << "WARN: ignoring unknown argument '" << a << "'\n";
        }
    }

    if (!(keep_p >= 0.0 && keep_p <= 1.0))
    {
        std::cerr << "ERROR: --keep-p must be in [0,1], got " << keep_p << "\n";
        return 1;
    }

    // If PKD2: epoch must be Since2000 (t2000). If user passes unix, fail loudly.
    if (ends_with_ci(seeds_path, ".pkd2") && epoch==SeedTimeEpoch::Unix)
    {
        std::cerr << "FATAL: SEEDS_PATH is .pkd2; seed times are already t2000. Do NOT pass 'unix'.\n";
        return 2;
    }

    {
        std::error_code ec;
        std::filesystem::create_directories(out_dir, ec);
        if (ec)
        {
            std::cerr << "Cannot create output dir '" << out_dir << "': " << ec.message() << "\n";
            return 1;
        }
    }

    const int64_t T0_EPOCH_S = 0;

    const auto pkd_pos_path=(out_dir / std::format("trajectories_{}h.pkd2",timestep_h)).string();
    const auto pkd_vel_path=(out_dir / std::format("velocities_{}h.pkd2",  timestep_h)).string();

    pkd2::Writer pkd_pos(pkd_pos_path,T0_EPOCH_S,base_dt_s,/*pos_scale=*/1'000'000u,/*zlib=*/6);
    pkd2::Writer pkd_vel(pkd_vel_path,T0_EPOCH_S,base_dt_s,/*pos_scale=*/1'000'000u,/*zlib=*/6);

    // -------- LOAD SEEDS (build explicit start/stop records) --------
    std::vector<SeedInit> seeds;
    {
        auto starts = load_seeds_auto(seeds_path, seed_point);
        if (starts.empty())
        {
            std::cerr<<"No seeds loaded from "<<seeds_path<<"\n";
            return 1;
        }

        if (ends_with_ci(seeds_path, ".pkd2"))
        {
            auto stops = load_seeds_auto(seeds_path, opposite_seed_point(seed_point));
            if (stops.size() != starts.size())
                throw std::runtime_error("PKD2: starts.size() != stops.size() (corrupt TOC or mismatch)");

            seeds.reserve(starts.size());
            for (size_t i=0;i<starts.size();++i)
            {
                if (starts[i].seed_id != stops[i].seed_id)
                    throw std::runtime_error("PKD2: seed_id mismatch between start/stop streams");

                SeedInit r{};
                r.seed_id        = starts[i].seed_id;

                // start = user-chosen seed_point
                r.lon_start_rad  = starts[i].lon_rad;
                r.lat_start_rad  = starts[i].lat_rad;
                r.t_start_s      = starts[i].t_start_s;

                // stop = opposite point
                r.has_stop       = true;
                r.lon_stop_rad   = stops[i].lon_rad;
                r.lat_stop_rad   = stops[i].lat_rad;
                r.t_stop_s       = stops[i].t_start_s;

                seeds.push_back(r);
            }
        }
        else
        {
            // CSV: start only
            seeds = std::move(starts);
            for (auto& s : seeds)
            {
                s.lon_start_rad = s.lon_rad;
                s.lat_start_rad = s.lat_rad;
                s.has_stop      = false;
                s.t_stop_s      = 0;
                s.lon_stop_rad  = 0;
                s.lat_stop_rad  = 0;
            }
        }
    }

    normalize_seed_times_to_t2000(seeds, epoch);
    snap_seed_times_to_grid(seeds, base_dt_s);

    auto minmax = std::minmax_element(seeds.begin(), seeds.end(),
        [](auto const& a, auto const& b){ return a.t_start_s < b.t_start_s; });
    
    std::cerr << "SEEDS t_start_s min=" << minmax.first->t_start_s
              << " max=" << minmax.second->t_start_s << "\n";
    std::cerr << "SEEDS month min=" << formatDate(month_start_from_t2000(minmax.first->t_start_s))
              << " max=" << formatDate(month_start_from_t2000(minmax.second->t_start_s)) << "\n";


    // Critical invariant for BACKWARD integration:
    // start must be later-or-equal to stop.
    for (const auto& s : seeds)
    {
        if (s.has_stop && s.t_stop_s > s.t_start_s)
        {
            throw std::runtime_error(
                "Seed ordering invalid for backtracking: t_stop_s > t_start_s for seed_id="
                + std::to_string(s.seed_id)
                + ". Use --seed-point last for backtracking, or fix seed construction.");
        }
    }

    // Thin (keeps PKD2 start/stop paired because it's one record)
    if (keep_p<1. || cap!=std::numeric_limits<std::size_t>::max())
    {
        std::mt19937 rng;
        if (use_thin_seed) rng.seed(thin_seed);
        else rng.seed(std::random_device{}());

        std::bernoulli_distribution keep(keep_p);

        std::size_t w=0;
        const std::size_t N=seeds.size();
        for (std::size_t i=0;i<N && w<cap;i++)
        {
            if (keep_p>=1.0 || keep(rng))
                seeds[w++]=std::move(seeds[i]);
        }
        seeds.resize(w);

        std::cerr<<"After thinning: "<<seeds.size()<<"\n";
    }

    SeedsByMonth bym;
    bucketize_seeds(seeds, bym);

    Date cur_month=end_release.beginOfMonth();

    buf[older]=loadOneMonth(cur_month,level,pattern_1,pattern_2,grid_variable_names,data_variable_names,dt_seconds);
    buildKdTrees(buf[older]);

    buf[newer]=loadOneMonth(cur_month.addMonths(-1),level,pattern_1,pattern_2,grid_variable_names,data_variable_names,dt_seconds);
    buildKdTrees(buf[newer]);

    // 24 months back (e.g. 1997-01 -> 1995-01)
    const Date first_month=start_release.addMonths(-24).beginOfMonth();

    while (cur_month>=first_month)
    {
        const Date win_beg=cur_month;
        const Date win_end=cur_month.addMonths(1).addDays(-1);

        activate_month(bym,win_beg,win_end,base_dt_s);

        propagateWindow(win_beg,win_end,buf[older],buf[newer],data_variable_names,dt_seconds,life_limit_seconds,pkd_pos,pkd_vel,base_dt_s);

        for (auto it=all_particles.begin();it!=all_particles.end();)
        {
            auto& vec=it->second;
            vec.erase(
                std::remove_if(vec.begin(),vec.end(),[](const auto& p){return p.reason!=EndReason::none;}),
                vec.end());
            if (vec.empty()) it=all_particles.erase(it);
            else ++it;
        }

        cur_month=cur_month.addMonths(-1);
        older^=1; newer^=1;

        if (cur_month>first_month)
        {
            buf[newer]=MonthData{};
            buf[newer]=loadOneMonth(cur_month.addMonths(-1),level,pattern_1,pattern_2,grid_variable_names,data_variable_names,dt_seconds);
            buildKdTrees(buf[newer]);
        }
    }

    // flush survivors
    for (auto& bucket : all_particles)
        for (auto& p : bucket.second)
            p.traj.flush_to(pkd_pos,pkd_vel,/*stride=*/1);

    pkd_pos.close();
    pkd_vel.close();

    const auto wall_end=std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds=wall_end-wall_start;
    std::time_t end_time=std::chrono::system_clock::to_time_t(wall_end);
    std::cerr<<"\nfinished computation at "<<std::ctime(&end_time)
             <<"elapsed time: "<<elapsed_seconds.count()<<" s\n";

    return 0;
}
catch (const std::exception& e)
{
    std::cerr << "FATAL: " << e.what() << "\n";
    return 2;
}
catch (...)
{
    std::cerr << "FATAL: unknown exception\n";
    return 3;
}

