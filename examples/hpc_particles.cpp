#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/global_control.h"
#include "oneapi/tbb/parallel_for.h" // TBB parallel_for
#include "particle_tracker/data_prep_utils.hpp" // Provides grid, time series, splines, etc.
#include "particle_tracker/detail/datetime_utils.hpp"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime> // for std::time_t, std::ctime
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream> // for std::cerr
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <unordered_map>

#include "particle_tracker/fixed_point/fixed_point_core.hpp" // the fixed-point pack/unpack
#include "particle_tracker/fixed_point/fixed_point_io_codec.hpp" // the gzip helpers: write_records_gzip/read_records_gzip

using particle_tracker::DataType;
using particle_tracker::dateFromSecondsSince2000;
using particle_tracker::EndReason;
using particle_tracker::formatDate;
using particle_tracker::inverseTransform;
using particle_tracker::KDTree;
using particle_tracker::loadOneMonth;
using particle_tracker::MonthData;
using particle_tracker::Neighbor;
using particle_tracker::NeighborsData;
using particle_tracker::parseDate;
using particle_tracker::Particle;
using particle_tracker::point;
using particle_tracker::PointCloud;
using particle_tracker::secondsSince2000;
using particle_tracker::splitArgs;
using particle_tracker::TimeSpline;
using particle_tracker::TimeType;

// -----------------------------------------------------------------------------
//  Globals |
// -----------------------------------------------------------------------------

// inline constexpr double number_pi = std::numbers::pi_v<double>;
static const TimeType life_time_seconds = 2 * 365 * 24 * 3600.; // 2 yr
static const TimeType grace_seconds = 7 * 24 * 3600.;           // ~7 days
static PointCloud<DataType> initial_band_cloud;
static std::unique_ptr<KDTree<DataType>> kd_initial_band;
static DataType stay_initial_band_thresh_squared;

// container keyed by release month
using DayKey = Date;
using ParticlePacket = std::vector<Particle<DataType, TimeType>>;
static std::unordered_map<DayKey, ParticlePacket>
    all_particles; // new hash specialization is needed here
static std::unordered_map<DayKey, std::uint64_t> released_total;
static std::unordered_map<DayKey, std::uint64_t> landed_total;
static std::unordered_map<DayKey, std::uint64_t> landed_east;
static std::unordered_map<DayKey, std::uint64_t> landed_west;
static std::unordered_map<DayKey, std::uint64_t>
    boring_total; // never left origin band (killed at grace)
static std::unordered_map<DayKey, std::uint64_t>
    neverland_total; // left origin but never landed

// Random number generator
static thread_local std::mt19937 gen(std::random_device{}());

// Helper: uniform random shift in [-max,max]
template <class DataType> inline DataType randomShift(DataType max) {
  std::uniform_real_distribution<DataType> dis(-max, max);
  return dis(gen);
}

struct BBox {
  DataType lon_min{}, lon_max{}, lat_min{}, lat_max{};
  bool wraps{};
  static inline DataType norm_pi(DataType L) noexcept {
    while (L <= -number_pi)
      L += 2 * number_pi;
    while (L > number_pi)
      L -= 2 * number_pi;
    return L;
  };
  BBox(DataType lo_min, DataType lo_max, DataType la_min, DataType la_max) {
    // normalize tiny numeric noise
    const DataType eps = static_cast<DataType>(1e-12);
    lon_min = norm_pi(lo_min);
    lon_max = norm_pi(lo_max);
    lat_min = la_min;
    lat_max = la_max;
    if (std::abs(lon_min - lon_max) < eps)
      lon_max = lon_min; // degenerate ok
    wraps = (lon_min > lon_max);
  }
  bool contains(DataType lon, DataType lat) const noexcept {
    if (lat < lat_min || lat > lat_max)
      return false;
    lon = norm_pi(lon); // make robust to any input range
    if (!wraps)
      return lon >= lon_min && lon <= lon_max;
    return lon >= lon_min || lon <= lon_max;
  }
};

// ============================================================================
// ONE-TIME INITIAL POSITION GENERATOR + DAILY SPAWNER
// ============================================================================

// Global cache of initial positions (reused every day)
static std::vector<point<DataType>> init_positions;

// Condition that the point is in the narrow band near the coast
inline bool inNearCoastBand(DataType lon, DataType lat, const MonthData &md,
                            DataType sea_tol_squared, DataType land_min_squared,
                            DataType land_max_squared) {
  if (!md.kd_sea || !md.kd_land)
    throw std::runtime_error("KD trees not built");

  DataType q[2] = {lon, lat};

  // nearest SEA point check
  {
    size_t idx;
    DataType dist_squared;
    nanoflann::KNNResultSet<DataType> rs(1);
    rs.init(&idx, &dist_squared);
    md.kd_sea->findNeighbors(rs, q, nanoflann::SearchParameters());
    if (dist_squared > sea_tol_squared)
      return false;
  }

  // band relative to LAND
  {
    size_t idx;
    DataType dist_squared;
    nanoflann::KNNResultSet<DataType> rs(1);
    rs.init(&idx, &dist_squared);
    md.kd_land->findNeighbors(rs, q, nanoflann::SearchParameters());
    if (dist_squared < land_min_squared || dist_squared > land_max_squared)
      return false;
  }

  return true;
}

std::vector<point<DataType>>
buildCoastalBand(const MonthData &md, const BBox &box, DataType sea_tol_squared,
                 DataType land_min_squared, DataType land_max_squared) {
  std::vector<point<DataType>> band;

  const std::size_t N = md.sea_cloud.kdtree_get_point_count();

  band.reserve(N / 20);

  for (size_t i = 0; i < N; i++) {
    DataType lon = md.sea_cloud.kdtree_get_pt(i, 0);
    DataType lat = md.sea_cloud.kdtree_get_pt(i, 1);
    if (!box.contains(lon, lat))
      continue;

    if (inNearCoastBand(lon, lat, md, sea_tol_squared, land_min_squared,
                        land_max_squared))
      band.push_back({lon, lat});
  }

  if (band.empty())
    std::cerr << "[WARN] buildCoastalBand(): empty band\n";

  return band;
}

void buildInitialPositions(std::size_t N, const MonthData &md_ref,
                           const BBox &box, DataType sea_tol, DataType land_min,
                           DataType land_max, DataType jitter) {
  if (!init_positions.empty())
    return; // already built

  const DataType sea_tol_squared = sea_tol * sea_tol;
  const DataType land_min_squared = land_min * land_min;
  const DataType land_max_squared = land_max * land_max;

  auto band = buildCoastalBand(md_ref, box, sea_tol_squared, land_min_squared,
                               land_max_squared);
  if (band.empty()) {
    std::cerr << "[ERROR] no valid band points for initial seeding\n";
    return;
  }

  // --------- build coastal-band KD-tree ----------
  initial_band_cloud.pts.clear();
  initial_band_cloud.pts.reserve(band.size() * 2);
  for (auto &p : band) {
    initial_band_cloud.pts.push_back(p[0]);
    initial_band_cloud.pts.push_back(p[1]);
  }
  auto kd_params = nanoflann::KDTreeSingleIndexAdaptorParams(10);
  kd_initial_band =
      std::make_unique<KDTree<DataType>>(2, initial_band_cloud, kd_params);
  kd_initial_band->buildIndex();

  // Set "still-in-band" threshold once (example: 10 km)
  DataType stay_initial_band_rad = 10000. / 6371000.;
  stay_initial_band_thresh_squared =
      stay_initial_band_rad * stay_initial_band_rad;

  std::uniform_int_distribution<std::size_t> pick(0, band.size() - 1);

  init_positions.reserve(N);

  std::size_t attempts = 0;
  const std::size_t max_attempts = N * 20; // heuristic

  for (std::size_t accepted = 0; accepted < N && attempts < max_attempts;
       attempts++) {
    auto wrap_pi = [](DataType L) {
      // wrap to (-pi, pi]
      while (L <= -number_pi)
        L += 2 * number_pi;
      while (L > number_pi)
        L -= 2 * number_pi;
      return L;
    };
    const auto &s = band[pick(gen)];
    DataType lon = wrap_pi(s[0] + randomShift(jitter));
    DataType lat = s[1] + randomShift(jitter);

    if (inNearCoastBand(lon, lat, md_ref, sea_tol_squared, land_min_squared,
                        land_max_squared)) {
      init_positions.push_back({lon, lat});
      ++accepted;
    }
  }

  if (init_positions.size() != N) {
    std::cerr << "Failed to generate " << N << " points (got "
              << init_positions.size()
              << "). Loosen thresholds or enlarge jitter/band.\n";
  }
  if (init_positions.empty()) {
    throw std::runtime_error("Failed to seed initial positions (band empty or "
                             "thresholds too strict).");
  }
}

void spawn(const Date &day, ParticlePacket &packet, TimeType release_time) {
  // count once per day
  released_total[day] += static_cast<std::uint64_t>(init_positions.size());

  packet.reserve(packet.size() + init_positions.size());
  for (auto const &xy : init_positions) {
    packet.push_back(Particle<DataType, TimeType>{xy[0], xy[1], release_time,
                                                  xy[0], xy[1], release_time,
                                                  EndReason::none});
  }
}

// ---------------------------------------------------------------
// Distance-to-land check (uses squared threshold)
// ---------------------------------------------------------------
inline bool nearLand(const point<DataType> &p, const MonthData &md,
                     DataType thresh_squared) {
  if (!md.kd_land)
    throw std::runtime_error("kd_land not built");
  size_t idx;
  DataType d_squared;
  nanoflann::KNNResultSet<DataType> rs(1);
  rs.init(&idx, &d_squared);
  md.kd_land->findNeighbors(rs, p.data(), nanoflann::SearchParameters());
  return d_squared <= thresh_squared;
}

inline bool stillInitialBand(const point<DataType> &p) {
  if (!kd_initial_band)
    return false; // safety
  size_t idx;
  DataType d_squared;
  nanoflann::KNNResultSet<DataType> rs(1);
  rs.init(&idx, &d_squared);
  kd_initial_band->findNeighbors(rs, p.data(), nanoflann::SearchParameters());
  return d_squared <= stay_initial_band_thresh_squared;
}

// ---------------------------------------------------------------
// Pick which month's data to use for absolute time t
// ---------------------------------------------------------------
inline const MonthData &pickMonth(TimeType t, const MonthData &older_m,
                                  const MonthData &newer_m) {
  Date d = dateFromSecondsSince2000(t);
  return (d.beginOfMonth() == older_m.month) ? older_m : newer_m;
}

inline NeighborsData<DataType, TimeType> &
pickNeighbors(TimeType t, const MonthData &older_m, const MonthData &newer_m,
              NeighborsData<DataType, TimeType> &nb_old,
              NeighborsData<DataType, TimeType> &nb_new) {
  Date d = dateFromSecondsSince2000(t);
  return (d.beginOfMonth() == older_m.month) ? nb_old : nb_new;
}

// ---------------------------------------------------------------
// Interpolate U,V using NeighborsData
// ---------------------------------------------------------------
inline void
interpolateUV(TimeType t, DataType lon, DataType lat,
              const std::vector<std::string> &vars, // {"water_u","water_v"}
              const MonthData &md, NeighborsData<DataType, TimeType> &nb,
              std::vector<Neighbor<DataType>> &neigh, // <-- scratch
              DataType &u, DataType &v) {
  neigh.clear(); // reuse capacity
  nb.computeNeighbors(neigh, lon, lat);
  if (neigh.empty()) {
    u = v = 0;
    return;
  }

  auto vals = nb.interpolateVariables(t, neigh, md.spl, vars);
  u = vals.at(vars[0]);
  v = vals.at(vars[1]);
}

// ---------------------------------------------------------------------------
// propagateWindow(): back-propagate all alive particles over [win_beg, win_end]
// using explicit Euler, spline interpolation & kd-trees.
// Stops on: land hit, lifetime expired, window boundary.
// ---------------------------------------------------------------------------
void propagateWindow(const Date &win_beg, const Date &win_end,
                     const MonthData &older_m, const MonthData &newer_m,
                     const std::vector<std::string> &data_vars,
                     TimeType dt_seconds, DataType land_thresh_rad,
                     TimeType life_time_seconds) {
  const TimeType t_beg = secondsSince2000(win_beg); // inclusive (earliest)
  // const TimeType t_end=secondsSince2000(win_end);      // inclusive (latest)
  const TimeType dt_sec = dt_seconds;

  DataType search_radius = static_cast<DataType>(.08 * 1.1 / 180. * number_pi);
  DataType shape_param = search_radius * search_radius / (1.1 * 1.1);

  // Build neighbor structures for both buffers (no static: rebuild each window)
  auto nb_old = makeNeighborsData(older_m, search_radius, shape_param);
  auto nb_new = makeNeighborsData(newer_m, search_radius, shape_param);

  const DataType land_thresh_squared = land_thresh_rad * land_thresh_rad;

  std::atomic<uint64_t> a_past_grace{0}, c_allow{0}, c_beach{0};

  // Iterate through *all* release-month buckets
  for (auto &bucket : all_particles) {
    auto &vec = bucket.second;

    if (!older_m.kd_sea || !newer_m.kd_sea)
      throw std::runtime_error("NeighborsData: kd_sea missing");

    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<std::size_t>(0, vec.size(), 4096), // 32768
        [&](const oneapi::tbb::blocked_range<std::size_t> &r) {
          std::vector<Neighbor<DataType>> neigh;
          neigh.reserve(8); // typical neighbor count (4/8/16/32/64/128)

          for (std::size_t i = r.begin(); i != r.end(); i++) {
            auto &p = vec[i];
            // skip dead particles
            if (p.reason != EndReason::none)
              continue;

            bool grace_done = false;

            const bool started_past_grace =
                (p.init_time - p.time) >= grace_seconds;
            bool grace_counted =
                started_past_grace; // already past -> don't count a crossing

            while (p.reason == EndReason::none && p.time >= t_beg) {
              TimeType age = p.init_time - p.time;

              if (!grace_counted && age >= grace_seconds) {
                a_past_grace.fetch_add(1, std::memory_order_relaxed);
                grace_counted = true;
              }

              if (age >= life_time_seconds) {
                p.reason = EndReason::lifetime;
                break;
              }

              const MonthData &md = pickMonth(p.time, older_m, newer_m);
              auto &nb =
                  pickNeighbors(p.time, older_m, newer_m, nb_old, nb_new);

              //[[maybe_unused]] bool in_band  =stillInitialBand({p.lon,p.lat});
              bool past_grace = (age >= grace_seconds);
              // bool allow_beach=(past_grace || !in_band);
              bool allow_beach = past_grace;
              if (allow_beach)
                c_allow.fetch_add(1, std::memory_order_relaxed);

              if (allow_beach &&
                  nearLand({p.lon, p.lat}, md, land_thresh_squared)) {
                c_beach.fetch_add(1, std::memory_order_relaxed);
                p.reason = EndReason::landed;
                break;
              }

              // 2) (Optional) Only *mark* never-left at grace,
              //    but don't break here  -  let it keep going further back in
              //    time.
              if (!grace_done && past_grace) {
                grace_done = true;
              }

              DataType u, v;
              interpolateUV(p.time, p.lon, p.lat, data_vars, md, nb, neigh, u,
                            v);
              DataType dx = -u * dt_sec, dy = -v * dt_sec;
              inverseTransform(dx, dy, p.lon, p.lat);

              auto wrap_lon_inplace = [](DataType &L) {
                while (L <= -number_pi)
                  L += 2 * number_pi;
                while (L > number_pi)
                  L -= 2 * number_pi;
              };
              auto clamp_lat_inplace = [](DataType &L) {
                if (L < -number_pi / 2)
                  L = -number_pi / 2;
                if (L > number_pi / 2)
                  L = number_pi / 2;
              };

              wrap_lon_inplace(p.lon);
              clamp_lat_inplace(p.lat);

              p.time -= dt_sec;
              if (p.time < t_beg)
                break;
            }
          }
        });
  }

  std::cerr << "window " << formatDate(win_beg)
            << ": past_grace=" << a_past_grace.load()
            << " allow_beach=" << c_allow.load()
            << " beached=" << c_beach.load() << "\n";
}

inline PackedParticle packBinary(const Particle<DataType, TimeType> &p) {
  uint8_t flags = static_cast<uint8_t>(p.reason); // 0..3

  return encode_record<DataType, TimeType>(p.init_lon,  // was p.lon0
                                           p.init_lat,  // was p.lat0
                                           p.init_time, // was p.t0
                                           p.lon,       // was p.lon1
                                           p.lat,       // was p.lat1
                                           p.time,      // was p.t1
                                           flags);
}

inline std::string dayTag(const Date &d) {
  std::string out = "trajectories_" + formatDate(d) + ".pkd.gz";
  return out;
}

inline std::filesystem::path outPath(const std::filesystem::path &dir,
                                     const Date &d, int level) {
  const auto filename =
      std::format("trajectories_{}_{}m.pkd.gz", formatDate(d), level);
  return dir / filename;
}

// --------------------------------------------------------------------
// Main
// --------------------------------------------------------------------
int main(int argc, char *argv[]) try {
  if (argc != 12 && argc != 13) {
    std::cerr << "usage: " << argv[0]
              << " REL_START REL_END LON_MIN LON_MAX LAT_MIN LAT_MAX LEVEL "
                 "TIMESTEP PATTERN1 PATTERN2 OUT_DIR"
              << " [[GRID_LON,GRID_LAT,DATA1,DATA2,...]]\n";
    return 1;
  }

  // counter on evaluation time
  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  std::cerr << "started computation at " << std::ctime(&start_time) << '\n';

  /* release window ----------------------------------------------------- */
  Date start_release = parseDate(argv[1]);
  Date end_release = parseDate(argv[2]);
  // DateRange sim_range(start_release,end_release);

  /* spatial filter ----------------------------------------------------- */
  DataType lon_min = std::stof(argv[3]) * number_pi / 180.;
  DataType lon_max = std::stof(argv[4]) * number_pi / 180.;
  DataType lat_min = std::stof(argv[5]) * number_pi / 180.;
  DataType lat_max = std::stof(argv[6]) * number_pi / 180.;

  int level = std::stoi(argv[7]);
  int timestep = std::stoi(argv[8]);
  std::string pattern_1 = argv[9];
  std::string pattern_2 = argv[10];
  std::filesystem::path out_dir = argv[11];
  std::error_code ec;
  std::filesystem::create_directories(out_dir, ec);
  if (ec) {
    std::cerr << "Cannot create output dir '" << out_dir
              << "': " << ec.message() << "\n";
    return 1;
  }

  // ------------------------------------------------------------------
  // variable names (default to water_u / water_v)
  // ------------------------------------------------------------------
  std::vector<std::string> data_variable_names = {"water_u", "water_v"};
  std::vector<std::string> grid_variable_names = {"lon", "lat", "water_u"};

  if (argc == 13) {
    auto names = splitArgs(argv[12]); // lon,lat,data1,...

    if (names.size() < 3) {
      std::cerr << "ERROR: variable list must contain lon_var,lat_var,at least "
                   "one data_var\n";
      return 1;
    }

    // replace grids
    grid_variable_names[0] = names[0]; // lon variable
    grid_variable_names[1] = names[1]; // lat variable
    grid_variable_names[2] = names[2]; // first data-var (for NaN mask)

    // replace data list
    data_variable_names.assign(names.begin() + 2, names.end());
  }

  // TimeType dt_hours=timestep;           // in hours
  TimeType dt_seconds = static_cast<TimeType>(timestep) * 3600.;
  Date cur_month = end_release.beginOfMonth(); // e.g. 2004-01-01

  // Basic sanity on CLI angles (degrees expected, convert to rad above)
  auto in_range = [](double x, double a, double b) { return x >= a && x <= b; };

  if (!in_range(std::stod(argv[3]), -180., 180.) ||
      !in_range(std::stod(argv[4]), -180., 180.) ||
      !in_range(std::stod(argv[5]), -90., 90.) ||
      !in_range(std::stod(argv[6]), -90., 90.)) {
    std::cerr << "Longitude/latitude degrees out of range. "
              << "Use lon in [-180,180], lat in [-90,90].\n";
    return 1;
  }

  if (timestep <= 0) {
    std::cerr << "TIMESTEP must be positive (hours).\n";
    return 1;
  }

  MonthData buf[2];
  static int older = 0, newer = 1;

  // load Month 0 and Month -1 (+padding handled inside)
  buf[older] =
      loadOneMonth(cur_month, level, pattern_1, pattern_2, grid_variable_names,
                   data_variable_names, dt_seconds);
  buildKdTrees(buf[older]);

  buf[newer] =
      loadOneMonth(cur_month.addMonths(-1), level, pattern_1, pattern_2,
                   grid_variable_names, data_variable_names, dt_seconds);
  buildKdTrees(buf[newer]);

  /* -------------------------------------------------------------------------
   */
  /*  MAIN MONTH-ADVANCING LOOP */
  /* -------------------------------------------------------------------------
   */

  BBox bbox{lon_min, lon_max, lat_min, lat_max};

  DataType sea_tol = .2 * number_pi / 180.;
  DataType land_min = 2000. / 6371000.;
  DataType land_max = 50000. / 6371000.;
  DataType jitter = .0005;
  std::size_t N_per_day = 2'800'000; // 2'000'000

  buildInitialPositions(N_per_day, buf[older], bbox, sea_tol, land_min,
                        land_max, jitter);

  // const Date last_month=end_release.beginOfMonth();           // e.g.
  // 2004-12-01
  const Date first_month =
      start_release.addMonths(-24).beginOfMonth(); // earliest month

  while (cur_month >= first_month) {
    /* 1. one-month simulation window [cur_month , cur_month+1) ------------- */
    Date win_beg = cur_month;
    Date win_end = cur_month.addMonths(1).addDays(-1); // inclusive last day

    for (Date day = win_beg; day <= win_end; day.increment()) {
      if (day < start_release || day > end_release)
        continue; // <-- guard
      TimeType rel_time = secondsSince2000(day);

      // MonthKey key=day.beginOfMonth();
      DayKey key = day;
      ParticlePacket &pack = all_particles[key];

      spawn(day, pack, rel_time);
    }

    DataType land_thresh = (0.5) * number_pi / 180.; // 5000./6371000.;
    propagateWindow(win_beg, win_end, buf[older], buf[newer],
                    data_variable_names, dt_seconds, land_thresh,
                    life_time_seconds);

    for (auto it = all_particles.begin(); it != all_particles.end();) {
      const DayKey day = it->first; // copy key (iterator may change)
      auto &vec = it->second;

      std::vector<PackedParticle> recs;
      recs.reserve(vec.size());

      vec.erase(std::remove_if(vec.begin(), vec.end(),
                               [&](const auto &p) {
                                 if (p.reason == EndReason::none)
                                   return false;

                                 switch (p.reason) {
                                 case EndReason::landed: {
                                   recs.push_back(packBinary(p));
                                   landed_total[day]++;

                                   // classify by longitude at landing
                                   double lon_deg = p.lon * 180. / number_pi;
                                   double lon360 = std::fmod(lon_deg + 360.,
                                                             360.); // [0,360)

                                   if (lon360 >= 180.)
                                     landed_west[day]++;
                                   else
                                     landed_east[day]++;
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

      if (!recs.empty()) {
        const auto path = outPath(out_dir, day, level);
        std::error_code ec2;
        const bool append = std::filesystem::exists(path, ec2) && !ec2;
        write_records_gzip(path.string(), recs, released_total[day],
                           /*level=*/6, /*append=*/append);
        // optional: std::cerr<<"write "<<recs.size()<<" to "<<path<<(append ? "
        // (append)\n" : " (new)\n");
      }

      if (vec.empty())
        it = all_particles.erase(it);
      else
        ++it;
    }

    std::uint64_t L = 0, B = 0, N = 0, R = 0, E = 0, W = 0;
    for (Date day = win_beg; day <= win_end; day.increment()) {
      L += landed_total[day];
      B += boring_total[day];
      N += neverland_total[day];
      R += released_total[day];
      E += landed_east[day];
      W += landed_west[day];
    }
    std::cerr << "[" << formatDate(win_beg) << "–" << formatDate(win_end)
              << "] " << "released=" << R << " landed=" << L
              << " landed east =" << E << " landed west=" << W << '\n'
              << "neverleft=" << B << " lifetime=" << N << "\n";

    /* 2. advance the ring buffer by one month ---------------------------- */
    cur_month = cur_month.addMonths(-1); // Feb -> Mar -> Apr ...

    older ^= 1; // swap indices (0 <-> 1)
    newer ^= 1;

    /* load the *new* Month+1 (two months ahead of 'older') */
    if (cur_month > first_month) {
      buf[newer] = MonthData{}; // drop old newer before loading fresh
      buf[newer] =
          loadOneMonth(cur_month.addMonths(-1), level, pattern_1, pattern_2,
                       grid_variable_names, data_variable_names, dt_seconds);
      buildKdTrees(buf[newer]);
    }
  }

  for (auto &[day, vec] : all_particles)
    for (auto &p : vec)
      if (p.reason == EndReason::none)
        neverland_total[day]++;

  all_particles.clear();

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cerr << '\n'
            << "finished computation at " << std::ctime(&end_time)
            << "elapsed time: " << elapsed_seconds.count() << "s\n";

  return 0;
} catch (const std::exception &e) {
  std::cerr << "FATAL: " << e.what() << "\n";
  return 2;
} catch (...) {
  std::cerr << "FATAL: unknown exception\n";
  return 3;
}
