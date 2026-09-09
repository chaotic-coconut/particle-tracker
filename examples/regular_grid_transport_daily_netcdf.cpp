// regular_grid_transport_daily_netcdf.cpp
//
// Daily release-batch transport example with NetCDF output.
//
// For every release day in [REL_START, REL_END], this program releases the
// same regular lon/lat seed grid, integrates each particle forward for a fixed
// duration, and writes one self-contained NetCDF file for that release day.
//
// Important design choices:
//   - Simulation is parallelized only over particles.
//   - NetCDF writing is serial and happens after the parallel integration has
//     completed, so there is no race between simulation and I/O.
//   - Velocity fields are held in a small rolling month cache. For daily
//     release dates and a 30-day integration, end-of-month releases can cross
//     two month boundaries, so this code allows up to MAX_CACHED_MONTHS = 3
//     months in RAM. Two months are not always sufficient for a release on,
//     say, 1997-01-31 with a 30-day integration.
//   - Output contains numeric status codes, not strings.
//
// Requires the netCDF-C++4 and netCDF-C libraries in addition to the
// dependencies of the other examples.

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"

#include "particle_tracker/data_prep_utils.hpp"
#include "particle_tracker/detail/datetime_utils.hpp"
#include "particle_tracker/geography_utils.hpp"
#include "particle_tracker/interpolation_utils.hpp"

#include <netcdf>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using particle_tracker::buildKdTrees;
using particle_tracker::DataType;
using particle_tracker::dateFromSecondsSince2000;
using particle_tracker::formatDate;
using particle_tracker::GCD;
using particle_tracker::inverseTransform;
using particle_tracker::loadOneMonth;
using particle_tracker::makeNeighborsData;
using particle_tracker::MonthData;
using particle_tracker::Neighbor;
using particle_tracker::NeighborsData;
using particle_tracker::parseDate;
using particle_tracker::point;
using particle_tracker::secondsSince2000;
using particle_tracker::splitArgs;
using particle_tracker::TimeType;

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------

static constexpr double EARTH_RADIUS_M = particle_tracker::earth_mean_radius_m;
static constexpr double DEFAULT_SEED_STEP_DEG = 0.08;
static constexpr double INTERP_RADIUS_DEG = 0.08;
static constexpr double LAND_THRESH_DEG = DEFAULT_SEED_STEP_DEG;
static constexpr std::size_t MAX_CACHED_MONTHS = 3;

// -----------------------------------------------------------------------------
// Status codes stored in NetCDF
// -----------------------------------------------------------------------------

enum class StatusCode : std::int8_t {
  ok = 1,
  landed = 2,
  no_sea = 3,
  start_land = 4,
  start_no_sea = 5
};

// -----------------------------------------------------------------------------
// Small coordinate helpers
// -----------------------------------------------------------------------------

inline DataType degToRad(double deg) {
  return static_cast<DataType>(deg * number_pi / 180.);
}

inline double radToDeg(DataType rad) {
  return static_cast<double>(rad) * 180. / number_pi;
}

inline double wrapDeg(double lon) {
  while (lon <= -180.)
    lon += 360.;
  while (lon > 180.)
    lon -= 360.;
  return lon;
}

// -----------------------------------------------------------------------------
// Domain geometry and field sampling
// -----------------------------------------------------------------------------

struct BBox {
  DataType lon_min{}, lon_max{}, lat_min{}, lat_max{};
  bool wraps{};
  bool full_lon{};

  static inline DataType norm_pi(DataType L) noexcept {
    while (L <= -number_pi)
      L += 2 * number_pi;
    while (L > number_pi)
      L -= 2 * number_pi;
    return L;
  };

  BBox(DataType lo_min, DataType lo_max, DataType la_min, DataType la_max) {
    const DataType eps = static_cast<DataType>(1e-12);

    lon_min = norm_pi(lo_min);
    lon_max = norm_pi(lo_max);
    lat_min = la_min;
    lat_max = la_max;

    full_lon =
        std::abs(lo_max - lo_min) >= static_cast<DataType>(2 * number_pi) - eps;

    if (std::abs(lon_min - lon_max) < eps)
      lon_max = lon_min;

    wraps = (lon_min > lon_max);
  }

  bool contains(DataType lon, DataType lat) const noexcept {
    if (lat < lat_min || lat > lat_max)
      return false;

    if (full_lon)
      return true;

    lon = norm_pi(lon);

    if (!wraps)
      return lon >= lon_min && lon <= lon_max;

    return lon >= lon_min || lon <= lon_max;
  }
};

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

inline void interpolateUV(TimeType t, DataType lon, DataType lat,
                          const std::vector<std::string> &vars,
                          const MonthData &md,
                          NeighborsData<DataType, TimeType> &nb,
                          std::vector<Neighbor<DataType>> &neigh, DataType &u,
                          DataType &v) {
  neigh.clear();
  nb.computeNeighbors(neigh, lon, lat);

  if (neigh.empty()) {
    u = v = 0;
    return;
  }

  auto vals = nb.interpolateVariables(t, neigh, md.spl, vars);

  u = vals.at(vars[0]);
  v = vals.at(vars[1]);
}

// -----------------------------------------------------------------------------
// Data structures
// -----------------------------------------------------------------------------

struct Seed {
  std::size_t lon_idx{};
  std::size_t lat_idx{};
  DataType lon{};
  DataType lat{};
};

struct TrajectoryResult {
  std::size_t seed_idx{};
  std::size_t lon_idx{};
  std::size_t lat_idx{};

  DataType lon0{};
  DataType lat0{};
  DataType lon{};
  DataType lat{};

  float mean_speed{};
  float path_length{};
  float displacement{};

  std::int32_t valid_steps{};
  StatusCode status{StatusCode::ok};
};

struct MonthBundle {
  Date month;
  MonthData data;
  std::unique_ptr<NeighborsData<DataType, TimeType>> neighbors;

  MonthBundle(Date m, MonthData &&d) : month(m), data(std::move(d)) {}

  MonthBundle(const MonthBundle &) = delete;
  MonthBundle &operator=(const MonthBundle &) = delete;
  MonthBundle(MonthBundle &&) noexcept = default;
  MonthBundle &operator=(MonthBundle &&) noexcept = default;
};

class MonthCache {
public:
  MonthCache(int level, std::string pattern_1, std::string pattern_2,
             std::vector<std::string> grid_variable_names,
             std::vector<std::string> data_variable_names, TimeType dt_seconds,
             DataType interp_radius, DataType shape_param)
      : level_(level), pattern_1_(std::move(pattern_1)),
        pattern_2_(std::move(pattern_2)),
        grid_variable_names_(std::move(grid_variable_names)),
        data_variable_names_(std::move(data_variable_names)),
        dt_seconds_(dt_seconds), interp_radius_(interp_radius),
        shape_param_(shape_param) {}

  MonthBundle &get(const Date &month) {
    for (auto &b : bundles_) {
      if (b.month == month)
        return b;
    }

    std::cerr << "loading month " << formatDate(month) << '\n';

    MonthData md =
        loadOneMonth(month, level_, pattern_1_, pattern_2_,
                     grid_variable_names_, data_variable_names_, dt_seconds_);

    // Important: MonthData must first be placed into its final cache location.
    // KD-trees and NeighborsData can keep pointers/references into MonthData.
    // Building them before moving MonthData into the cache can leave dangling
    // internal pointers and cause a crash at the first nearLand/interpolation
    // call.
    bundles_.emplace_back(month, std::move(md));
    MonthBundle &stored = bundles_.back();

    buildKdTrees(stored.data);
    stored.neighbors = std::make_unique<NeighborsData<DataType, TimeType>>(
        makeNeighborsData(stored.data, interp_radius_, shape_param_));

    if (bundles_.size() > MAX_CACHED_MONTHS) {
      std::cerr << "evicting month " << formatDate(bundles_.front().month)
                << '\n';
      bundles_.pop_front();
    }

    return stored;
  }

  const MonthBundle &getLoadedForTime(TimeType t) const {
    Date d = dateFromSecondsSince2000(t);
    Date m = d.beginOfMonth();

    for (const auto &b : bundles_) {
      if (b.month == m)
        return b;
    }

    throw std::runtime_error("MonthCache::getLoadedForTime(): requested month "
                             "was not prefetched: " +
                             formatDate(m));
  }

  void prefetchForRelease(const Date &release, TimeType integration_seconds) {
    const Date m0 = release.beginOfMonth();
    get(m0);

    const Date end_date = dateFromSecondsSince2000(secondsSince2000(release) +
                                                   integration_seconds);
    const Date mend = end_date.beginOfMonth();

    Date m = m0;
    while (m < mend || m == mend) {
      get(m);
      m = m.addMonths(1);
    }
  }

private:
  int level_;
  std::string pattern_1_;
  std::string pattern_2_;
  std::vector<std::string> grid_variable_names_;
  std::vector<std::string> data_variable_names_;
  TimeType dt_seconds_;
  DataType interp_radius_;
  DataType shape_param_;
  std::deque<MonthBundle> bundles_;
};

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

static void printUsage(const char *prog) {
  std::cerr
      << "usage: " << prog
      << " REL_START REL_END LON_MIN LON_MAX LAT_MIN LAT_MAX LEVEL "
         "DT_HOURS INTEGRATION_DAYS SEED_STEP_DEG PATTERN1 PATTERN2 OUT_DIR"
      << " [GRID_LON,GRID_LAT,DATA1,DATA2,...]\n\n"
      << "  REL_START REL_END       daily release-date window, YYYY-MM-DD\n"
      << "  LON_MIN LON_MAX         degrees in [-180,180]; lon_min > lon_max "
         "means dateline wrap\n"
      << "  LAT_MIN LAT_MAX         degrees in [-90,90]\n"
      << "  LEVEL                   vertical level, e.g. 0\n"
      << "  DT_HOURS                Euler step in hours\n"
      << "  INTEGRATION_DAYS        forward integration duration in days\n"
      << "  SEED_STEP_DEG           regular seed-grid spacing; <=0 uses 0.08 "
         "degree\n"
      << "  PATTERN1 PATTERN2       NetCDF patterns for the two file "
         "subdomains\n"
      << "  OUT_DIR                 output directory for daily NetCDF files\n"
      << "  optional variable list  lon,lat,data1,data2,...; data1 is mask "
         "var\n";
}

static std::vector<Seed> buildSeeds(const BBox &bbox, double lon_min_deg,
                                    double lon_max_deg, double lat_min_deg,
                                    double lat_max_deg, double step_deg) {
  std::vector<double> lon_values;
  std::vector<double> lat_values;

  auto addLonSegment = [&](double lo, double hi) {
    if (hi < lo)
      return;

    const std::size_t n =
        static_cast<std::size_t>(std::floor((hi - lo) / step_deg)) + 1;

    lon_values.reserve(lon_values.size() + n);

    for (std::size_t i = 0; i < n; ++i) {
      const double lon = lo + static_cast<double>(i) * step_deg;
      const double wrapped = wrapDeg(lon);

      if (!lon_values.empty() && std::abs(lon_values.back() - wrapped) < 1e-10)
        continue;

      lon_values.push_back(wrapped);
    }
  };

  if (std::abs(lon_max_deg - lon_min_deg) >= 360. - 1e-12) {
    addLonSegment(-180. + step_deg, 180.);
  } else if (lon_min_deg <= lon_max_deg) {
    addLonSegment(lon_min_deg, lon_max_deg);
  } else {
    addLonSegment(lon_min_deg, 180.);
    addLonSegment(-180. + step_deg, lon_max_deg);
  }

  const std::size_t n_lat = static_cast<std::size_t>(std::floor(
                                (lat_max_deg - lat_min_deg) / step_deg)) +
                            1;

  lat_values.reserve(n_lat);

  for (std::size_t j = 0; j < n_lat; ++j) {
    lat_values.push_back(lat_min_deg + static_cast<double>(j) * step_deg);
  }

  std::vector<Seed> seeds;
  seeds.reserve(lon_values.size() * lat_values.size());

  for (std::size_t j = 0; j < lat_values.size(); ++j) {
    for (std::size_t i = 0; i < lon_values.size(); ++i) {
      const DataType lon = degToRad(lon_values[i]);
      const DataType lat = degToRad(lat_values[j]);

      if (!bbox.contains(lon, lat))
        continue;

      seeds.push_back(Seed{i, j, lon, lat});
    }
  }

  return seeds;
}

static TrajectoryResult integrateSeed(std::size_t seed_idx, const Seed &seed,
                                      const Date &release,
                                      const MonthCache &month_cache,
                                      const std::vector<std::string> &vars,
                                      TimeType dt_seconds, std::size_t n_steps,
                                      DataType land_thresh_squared) {
  TrajectoryResult r;

  r.seed_idx = seed_idx;
  r.lon_idx = seed.lon_idx;
  r.lat_idx = seed.lat_idx;
  r.lon0 = seed.lon;
  r.lat0 = seed.lat;
  r.lon = seed.lon;
  r.lat = seed.lat;

  std::vector<Neighbor<DataType>> neigh;
  neigh.reserve(32);

  TimeType t = secondsSince2000(release);

  const MonthBundle &start_bundle = month_cache.getLoadedForTime(t);

  if (nearLand({r.lon, r.lat}, start_bundle.data, land_thresh_squared)) {
    r.status = StatusCode::start_land;
    return r;
  }

  double speed_sum = 0.;
  double path_length = 0.;

  for (std::size_t step = 0; step < n_steps; ++step) {
    const MonthBundle &bundle = month_cache.getLoadedForTime(t);

    DataType u{}, v{};

    interpolateUV(t, r.lon, r.lat, vars, bundle.data, *bundle.neighbors, neigh,
                  u, v);

    if (neigh.empty()) {
      r.status =
          (r.valid_steps == 0) ? StatusCode::start_no_sea : StatusCode::no_sea;
      break;
    }

    const double ud = static_cast<double>(u);
    const double vd = static_cast<double>(v);
    const double speed = std::hypot(ud, vd);

    speed_sum += speed;

    const double dx_m = ud * static_cast<double>(dt_seconds);
    const double dy_m = vd * static_cast<double>(dt_seconds);

    path_length += std::hypot(dx_m, dy_m);

    DataType dx = static_cast<DataType>(dx_m);
    DataType dy = static_cast<DataType>(dy_m);

    const TimeType t_next = t + dt_seconds;

    inverseTransform(dx, dy, r.lon, r.lat);
    r.lon = BBox::norm_pi(r.lon);

    ++r.valid_steps;

    const MonthBundle &after_bundle = month_cache.getLoadedForTime(t_next);

    if (nearLand({r.lon, r.lat}, after_bundle.data, land_thresh_squared)) {
      r.status = StatusCode::landed;
      break;
    }

    t = t_next;
  }

  r.path_length = static_cast<float>(path_length);

  if (r.valid_steps > 0)
    r.mean_speed =
        static_cast<float>(speed_sum / static_cast<double>(r.valid_steps));

  const double central_angle =
      static_cast<double>(GCD(r.lon0, r.lat0, r.lon, r.lat));
  r.displacement = static_cast<float>(EARTH_RADIUS_M * central_angle);

  return r;
}

static std::filesystem::path outPath(const std::filesystem::path &dir,
                                     const Date &d, int level) {
  std::ostringstream name;
  name << "regular_grid_transport_" << formatDate(d) << "_" << level << "m.nc";
  return dir / name.str();
}

static void addUnits(netCDF::NcVar &var, const std::string &units) {
  var.putAtt("units", units);
}

static void writeDailyNetCDF(const std::filesystem::path &path,
                             const Date &release, int level, double dt_hours,
                             double integration_days, double seed_step_deg,
                             const std::vector<TrajectoryResult> &results,
                             bool omit_start_failures) {
  std::vector<std::int64_t> particle_id;
  std::vector<std::int32_t> lon_idx;
  std::vector<std::int32_t> lat_idx;
  std::vector<float> lon_init_deg;
  std::vector<float> lat_init_deg;
  std::vector<float> lon_final_deg;
  std::vector<float> lat_final_deg;
  std::vector<float> mean_speed;
  std::vector<float> path_length;
  std::vector<float> displacement;
  std::vector<std::int32_t> valid_steps;
  std::vector<signed char> status;

  particle_id.reserve(results.size());
  lon_idx.reserve(results.size());
  lat_idx.reserve(results.size());
  lon_init_deg.reserve(results.size());
  lat_init_deg.reserve(results.size());
  lon_final_deg.reserve(results.size());
  lat_final_deg.reserve(results.size());
  mean_speed.reserve(results.size());
  path_length.reserve(results.size());
  displacement.reserve(results.size());
  valid_steps.reserve(results.size());
  status.reserve(results.size());

  for (const auto &r : results) {
    if (omit_start_failures && (r.status == StatusCode::start_land ||
                                r.status == StatusCode::start_no_sea))
      continue;

    particle_id.push_back(static_cast<std::int64_t>(r.seed_idx));
    lon_idx.push_back(static_cast<std::int32_t>(r.lon_idx));
    lat_idx.push_back(static_cast<std::int32_t>(r.lat_idx));
    lon_init_deg.push_back(static_cast<float>(radToDeg(r.lon0)));
    lat_init_deg.push_back(static_cast<float>(radToDeg(r.lat0)));
    lon_final_deg.push_back(static_cast<float>(radToDeg(BBox::norm_pi(r.lon))));
    lat_final_deg.push_back(static_cast<float>(radToDeg(r.lat)));
    mean_speed.push_back(r.mean_speed);
    path_length.push_back(r.path_length);
    displacement.push_back(r.displacement);
    valid_steps.push_back(r.valid_steps);
    status.push_back(static_cast<signed char>(r.status));
  }

  netCDF::NcFile nc(path.string(), netCDF::NcFile::replace);

  nc.putAtt("title",
            "Daily regular-grid forward particle transport diagnostics");
  nc.putAtt("release_date", formatDate(release));
  nc.putAtt("release_time_seconds_since_2000", netCDF::NcType::nc_INT64,
            static_cast<long long>(secondsSince2000(release)));
  nc.putAtt("level_m", netCDF::NcType::nc_INT, level);
  nc.putAtt("dt_hours", netCDF::NcType::nc_DOUBLE, dt_hours);
  nc.putAtt("integration_days", netCDF::NcType::nc_DOUBLE, integration_days);
  nc.putAtt("seed_step_deg", netCDF::NcType::nc_DOUBLE, seed_step_deg);
  nc.putAtt("status_mapping",
            "1=ok, 2=landed, 3=no_sea, 4=start_land, 5=start_no_sea");
  nc.putAtt(
      "note",
      "straightness is not stored; calculate displacement_m/path_length_m");

  auto particle_dim = nc.addDim("particle", particle_id.size());

  auto v_particle_id =
      nc.addVar("particle_id", netCDF::NcType::nc_INT64, particle_dim);
  auto v_lon_idx = nc.addVar("lon_idx", netCDF::NcType::nc_INT, particle_dim);
  auto v_lat_idx = nc.addVar("lat_idx", netCDF::NcType::nc_INT, particle_dim);
  auto v_lon0 =
      nc.addVar("lon_init_deg", netCDF::NcType::nc_FLOAT, particle_dim);
  auto v_lat0 =
      nc.addVar("lat_init_deg", netCDF::NcType::nc_FLOAT, particle_dim);
  auto v_lon =
      nc.addVar("lon_final_deg", netCDF::NcType::nc_FLOAT, particle_dim);
  auto v_lat =
      nc.addVar("lat_final_deg", netCDF::NcType::nc_FLOAT, particle_dim);
  auto v_mean_speed =
      nc.addVar("mean_speed_m_per_s", netCDF::NcType::nc_FLOAT, particle_dim);
  auto v_path_length =
      nc.addVar("path_length_m", netCDF::NcType::nc_FLOAT, particle_dim);
  auto v_displacement =
      nc.addVar("displacement_m", netCDF::NcType::nc_FLOAT, particle_dim);
  auto v_valid_steps =
      nc.addVar("valid_steps", netCDF::NcType::nc_INT, particle_dim);
  auto v_status = nc.addVar("status", netCDF::NcType::nc_BYTE, particle_dim);

  addUnits(v_lon0, "degrees_east");
  addUnits(v_lat0, "degrees_north");
  addUnits(v_lon, "degrees_east");
  addUnits(v_lat, "degrees_north");
  addUnits(v_mean_speed, "m s-1");
  addUnits(v_path_length, "m");
  addUnits(v_displacement, "m");
  v_status.putAtt("flag_values", "1, 2, 3, 4, 5");
  v_status.putAtt("flag_meanings", "ok landed no_sea start_land start_no_sea");

  if (!particle_id.empty()) {
    v_particle_id.putVar(particle_id.data());
    v_lon_idx.putVar(lon_idx.data());
    v_lat_idx.putVar(lat_idx.data());
    v_lon0.putVar(lon_init_deg.data());
    v_lat0.putVar(lat_init_deg.data());
    v_lon.putVar(lon_final_deg.data());
    v_lat.putVar(lat_final_deg.data());
    v_mean_speed.putVar(mean_speed.data());
    v_path_length.putVar(path_length.data());
    v_displacement.putVar(displacement.data());
    v_valid_steps.putVar(valid_steps.data());
    v_status.putVar(status.data());
  }
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

int main(int argc, char *argv[]) try {
  if (argc != 14 && argc != 15) {
    printUsage(argv[0]);
    return 1;
  }

  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);

  std::cerr << "started computation at " << std::ctime(&start_time) << '\n';

  Date start_release = parseDate(argv[1]);
  Date end_release = parseDate(argv[2]);

  if (end_release < start_release) {
    std::cerr << "REL_END must not be before REL_START\n";
    return 1;
  }

  auto in_range = [](double x, double a, double b) { return x >= a && x <= b; };

  const double lon_min_deg = std::stod(argv[3]);
  const double lon_max_deg = std::stod(argv[4]);
  const double lat_min_deg = std::stod(argv[5]);
  const double lat_max_deg = std::stod(argv[6]);

  if (!in_range(lon_min_deg, -180., 180.) ||
      !in_range(lon_max_deg, -180., 180.) ||
      !in_range(lat_min_deg, -90., 90.) || !in_range(lat_max_deg, -90., 90.)) {
    std::cerr << "Longitude/latitude degrees out of range. "
              << "Use lon in [-180,180], lat in [-90,90].\n";
    return 1;
  }

  if (lat_min_deg > lat_max_deg) {
    std::cerr << "LAT_MIN must be <= LAT_MAX\n";
    return 1;
  }

  const int level = std::stoi(argv[7]);
  const double dt_hours = std::stod(argv[8]);
  const double integration_days = std::stod(argv[9]);
  const double seed_step_arg = std::stod(argv[10]);

  if (level < 0) {
    std::cerr << "LEVEL must be non-negative\n";
    return 1;
  }

  if (dt_hours <= 0. || integration_days <= 0.) {
    std::cerr << "DT_HOURS and INTEGRATION_DAYS must be positive\n";
    return 1;
  }

  // The integration length is bounded by what the rolling cache can cover: a
  // trajectory may not span more calendar months than MAX_CACHED_MONTHS.
  if (integration_days > 62.) {
    std::cerr
        << "This daily-output version is intended for <=62 days. "
        << "Increase MAX_CACHED_MONTHS if you need longer integrations.\n";
    return 1;
  }

  const double seed_step_deg =
      (seed_step_arg > 0.) ? seed_step_arg : DEFAULT_SEED_STEP_DEG;

  const TimeType dt_seconds = static_cast<TimeType>(dt_hours * 3600.);
  const TimeType integration_seconds =
      static_cast<TimeType>(integration_days * 24. * 3600.);

  const double step_ratio = static_cast<double>(integration_seconds) /
                            static_cast<double>(dt_seconds);
  const std::size_t n_steps = static_cast<std::size_t>(std::floor(step_ratio));

  const double remainder = step_ratio - static_cast<double>(n_steps);
  if (remainder > 1e-12) {
    std::cerr << "WARNING: integration duration is not an integer multiple of "
              << "DT_HOURS; final fractional step is ignored.\n";
  }

  if (n_steps == 0) {
    std::cerr << "Integration duration is shorter than one time step\n";
    return 1;
  }

  const std::string pattern_1 = argv[11];
  const std::string pattern_2 = argv[12];
  const std::filesystem::path out_dir = argv[13];

  std::error_code ec;
  std::filesystem::create_directories(out_dir, ec);

  if (ec) {
    std::cerr << "Cannot create output dir '" << out_dir
              << "': " << ec.message() << "\n";
    return 1;
  }

  std::vector<std::string> data_variable_names = {"water_u", "water_v"};
  std::vector<std::string> grid_variable_names = {"lon", "lat", "water_u"};

  if (argc == 15) {
    auto names = splitArgs(argv[14]);

    if (names.size() < 4) {
      std::cerr
          << "ERROR: variable list must contain lon_var,lat_var,u_var,v_var\n";
      return 1;
    }

    grid_variable_names[0] = names[0];
    grid_variable_names[1] = names[1];
    grid_variable_names[2] = names[2];

    data_variable_names.assign(names.begin() + 2, names.end());
  }

  const BBox bbox{degToRad(lon_min_deg), degToRad(lon_max_deg),
                  degToRad(lat_min_deg), degToRad(lat_max_deg)};

  const std::vector<Seed> seeds = buildSeeds(
      bbox, lon_min_deg, lon_max_deg, lat_min_deg, lat_max_deg, seed_step_deg);

  if (seeds.empty()) {
    std::cerr << "No seed points were generated\n";
    return 1;
  }

  std::cerr << "seed_count=" << seeds.size()
            << " seed_step_deg=" << seed_step_deg << " n_steps=" << n_steps
            << " daily releases from " << formatDate(start_release) << " to "
            << formatDate(end_release) << '\n';

  const DataType interp_radius =
      static_cast<DataType>(INTERP_RADIUS_DEG * number_pi / 180.);
  const DataType shape_param =
      interp_radius * interp_radius / static_cast<DataType>(1.1 * 1.1);

  const DataType land_thresh =
      static_cast<DataType>(LAND_THRESH_DEG * number_pi / 180.);
  const DataType land_thresh_squared = land_thresh * land_thresh;

  MonthCache month_cache(level, pattern_1, pattern_2, grid_variable_names,
                         data_variable_names, dt_seconds, interp_radius,
                         shape_param);

  Date release = start_release;

  while (release <= end_release) {
    std::cerr << "running daily release " << formatDate(release) << '\n';

    // Load required months before entering the parallel loop. This avoids any
    // parallel thread triggering NetCDF reads or cache mutation.
    month_cache.prefetchForRelease(release, integration_seconds);

    std::vector<TrajectoryResult> results(seeds.size());

    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<std::size_t>(0, seeds.size(), 4096),
        [&](const oneapi::tbb::blocked_range<std::size_t> &range) {
          for (std::size_t i = range.begin(); i != range.end(); ++i) {
            results[i] = integrateSeed(i, seeds[i], release, month_cache,
                                       data_variable_names, dt_seconds, n_steps,
                                       land_thresh_squared);
          }
        });

    std::size_t n_ok = 0;
    std::size_t n_landed = 0;
    std::size_t n_no_sea = 0;
    std::size_t n_start_land = 0;
    std::size_t n_start_no_sea = 0;

    for (const auto &r : results) {
      if (r.status == StatusCode::ok)
        ++n_ok;
      else if (r.status == StatusCode::landed)
        ++n_landed;
      else if (r.status == StatusCode::no_sea)
        ++n_no_sea;
      else if (r.status == StatusCode::start_land)
        ++n_start_land;
      else if (r.status == StatusCode::start_no_sea)
        ++n_start_no_sea;
    }

    const auto path = outPath(out_dir, release, level);

    // Seeds that never entered the water carry no trajectory, so they are left
    // out of the file. Set this to false to have every seed ID appear in every
    // daily file regardless of status.
    constexpr bool OMIT_START_FAILURES = true;

    writeDailyNetCDF(path, release, level, dt_hours, integration_days,
                     seed_step_deg, results, OMIT_START_FAILURES);

    std::cerr << "summary "
              << "ok=" << n_ok << " landed=" << n_landed
              << " no_sea=" << n_no_sea << " start_land=" << n_start_land
              << " start_no_sea=" << n_start_no_sea << '\n';

    std::cerr << "finished release " << formatDate(release) << " -> " << path
              << '\n';

    release = dateFromSecondsSince2000(secondsSince2000(release) + 24 * 3600);
  }

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
