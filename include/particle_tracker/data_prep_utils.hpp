#ifndef DATA_PREP_UTILS_HPP
#define DATA_PREP_UTILS_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numbers> // std::numbers::pi_v
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <oneapi/tbb/parallel_for.h>

#include "particle_tracker/detail/ncdf_utils.hpp"
#include "particle_tracker/fixed_point/fixed_point_core.hpp"
#include "particle_tracker/interpolation_utils.hpp"
#include <particle_tracker/detail/datetime_utils.hpp>

// Ensure C++17 or higher for filesystem support
#if __cplusplus < 201703L
#error "C++17 or higher is required"
#endif

inline constexpr double number_pi = std::numbers::pi_v<double>;

// ---------------------------------------------------------------------------
// Optional debug logging (compile-time gated)
// ---------------------------------------------------------------------------
// Define DATA_PREP_DEBUG_TIME_STEPS to print time-step diagnostics.
// #define DATA_PREP_DEBUG_TIME_STEPS

namespace std {
template <> struct hash<::Date> {
  std::size_t operator()(const ::Date &d) const noexcept {
    return static_cast<std::size_t>(d.getYear()) * 32 * 13 +
           static_cast<std::size_t>(d.getMonth()) * 32 +
           static_cast<std::size_t>(d.getDay());
  }
};
} // namespace std

namespace particle_tracker {

// Grid struct templated on DataType
template <typename DT> struct Grid {
  std::vector<DT> lon_rad;
  std::vector<DT> lat_rad;
  std::vector<bool> in_the_sea;
  std::size_t num_points;

  Grid() : num_points(0) {}
};

// TimeStepData struct templated on DataType and TimeType
template <typename DT, typename TT = double> struct TimeStepData {
  TT time_value;
  std::map<std::string, std::vector<DT>> variable_data;
};

// ---------------------------------------------------------------------------
// Function to initialize the grid
//
// variable_names is expected to contain:
//   [0] lon variable name
//   [1] lat variable name
//   [2] a "mask variable" used to detect sea points (NaN = land)
// ---------------------------------------------------------------------------
/**
 * @brief Initialize a lon/lat grid and a sea-mask from a netCDF file.
 *
 * The netCDF variables lon and lat are expected to be 1D arrays.
 * The mask variable is expected to match the flattened (lon x lat) grid size.
 *
 * @tparam DataType numeric type for stored grid coordinates.
 * @param file_path path to the netCDF file.
 * @param variable_names vector containing lon_var, lat_var, mask_var.
 * @return Grid<DataType> populated grid (lon_rad, lat_rad, in_the_sea).
 *
 * @throws std::runtime_error on any inconsistency or missing data.
 */
template <typename DT>
Grid<DT> initializeGrid(const std::string &file_path,
                        const std::vector<std::string> &variable_names) {
  if (variable_names.size() < 3) {
    throw std::runtime_error("initializeGrid(): variable_names must contain "
                             "lon_var, lat_var, mask_var.");
  }

  const std::string &lon_var = variable_names[0];
  const std::string &lat_var = variable_names[1];
  const std::string &mask_var = variable_names[2];

  Grid<DT> grid;

  // Create the processor for the grid initialization
  NetCDFProcessor<DT> processor(file_path, variable_names);

  // Retrieve variable shapes
  const auto &lon_shape = processor.getVariableShape(lon_var);
  const auto &lat_shape = processor.getVariableShape(lat_var);

  // Retrieve variable data (lon/lat are treated as static and stored under key
  // 0.0)
  const auto &lon_data = processor.getVariableData(lon_var).at(0.0);
  const auto &lat_data = processor.getVariableData(lat_var).at(0.0);

  // Use first available time slice for mask_var (works for both time-dependent
  // and not)
  const auto &mask_map = processor.getVariableData(mask_var);
  if (mask_map.empty()) {
    throw std::runtime_error("initializeGrid(): mask variable '" + mask_var +
                             "' has no data.");
  }
  const auto &mask_data = mask_map.begin()->second;

  // Check if 'lon' and 'lat' are 1D arrays
  if (!(lon_shape.size() == 1 && lat_shape.size() == 1)) {
    throw std::runtime_error("initializeGrid(): unsupported dimensions for '" +
                             lon_var + "' and '" + lat_var +
                             "'. Expected 1D lon/lat arrays.");
  }

  // Sizes of 'lon' and 'lat'
  std::size_t n_lon = lon_shape[0];
  std::size_t n_lat = lat_shape[0];

  grid.num_points = n_lon * n_lat;

  grid.lon_rad.resize(grid.num_points);
  grid.lat_rad.resize(grid.num_points);
  grid.in_the_sea.resize(grid.num_points);

  // Convert 'lon' and 'lat' to radians
  const DT grad_to_rad_factor = DT(number_pi / 180.);
  std::vector<DT> lon_rad_1d(n_lon);
  std::vector<DT> lat_rad_1d(n_lat);

  for (std::size_t i = 0; i < n_lon; i++) {
    lon_rad_1d[i] = lon_data[i] * grad_to_rad_factor;
  }

  for (std::size_t j = 0; j < n_lat; j++) {
    lat_rad_1d[j] = lat_data[j] * grad_to_rad_factor;
  }

  // Create grid points by combining 'lon' and 'lat'
  for (std::size_t j = 0; j < n_lat; j++) {
    for (std::size_t i = 0; i < n_lon; i++) {
      std::size_t idx = j * n_lon + i;
      grid.lon_rad[idx] = lon_rad_1d[i];
      grid.lat_rad[idx] = lat_rad_1d[j];
    }
  }

  // Ensure that mask_data size matches grid.num_points
  if (mask_data.size() != grid.num_points) {
    throw std::runtime_error("initializeGrid(): mask variable '" + mask_var +
                             "' size does not match grid size.");
  }

  // Determine in_the_sea vector using mask variable (NaN = land)
  for (std::size_t idx = 0; idx < grid.num_points; idx++) {
    const DT m_value = mask_data[idx];
    grid.in_the_sea[idx] = !std::isnan(m_value);
  }

  return grid;
}

// ---------------------------------------------------------------------------
// Function to populate data from multiple files
// ---------------------------------------------------------------------------
/**
 * @brief Populate time-step data by reading variables from multiple files.
 *
 * This function fills all_time_steps_data with one entry per time_value.
 * For points outside the sea-mask, it writes NaN.
 *
 * @tparam DataType numeric type.
 * @tparam TimeType time type.
 * @param file_paths list of netCDF files.
 * @param grid initialized by initializeGrid().
 * @param all_time_steps_data output vector of time-step payloads.
 * @param data_variable_names list of variables to load.
 * @param check_land_points unused flag (kept for compatibility).
 */
template <typename DT, typename TT = double>
void populateData(const std::vector<std::string> &file_paths,
                  const Grid<DT> &grid,
                  std::vector<TimeStepData<DT, TT>> &all_time_steps_data,
                  const std::vector<std::string> &data_variable_names,
                  bool check_land_points = false) {
  (void)check_land_points;

  for (const auto &file_path : file_paths) {
    try {
      NetCDFProcessor<DT, TT> processor(file_path, data_variable_names);

      // Processing function for data population
      auto processFunc = [&](const std::string &var_name, TT time_value,
                             const std::vector<DT> &data_vector) {
        // Ensure data size matches grid size
        if (data_vector.size() != grid.num_points) {
          throw std::runtime_error("Data size (" +
                                   std::to_string(data_vector.size()) +
                                   ") does not match grid size (" +
                                   std::to_string(grid.num_points) + ").");
        }

        // Find or create TimeStepData for this time_value
        auto it =
            std::find_if(all_time_steps_data.begin(), all_time_steps_data.end(),
                         [&](const TimeStepData<DT, TT> &tsd) {
                           return tsd.time_value == time_value;
                         });

        if (it == all_time_steps_data.end()) {
          // Create a new TimeStepData
          TimeStepData<DT, TT> tsd;
          tsd.time_value = time_value;

          // Initialize data vectors for each variable
          for (const auto &var : data_variable_names) {
            tsd.variable_data[var].resize(grid.num_points,
                                          std::numeric_limits<DT>::quiet_NaN());
          }

          all_time_steps_data.push_back(std::move(tsd));
          it = std::prev(all_time_steps_data.end());
        }

        // Populate data for the current variable
        auto &var_data = it->variable_data[var_name];

        for (std::size_t i = 0; i < grid.num_points; i++) {
          if (grid.in_the_sea[i]) {
            var_data[i] = data_vector[i];
          } else {
            var_data[i] = std::numeric_limits<DT>::quiet_NaN();
          }
        }
      };

      // Process variables and populate data
      processor.processVariables(processFunc);
    } catch (const std::exception &e) {
      std::cerr << "Error processing data from file '" << file_path
                << "': " << e.what() << std::endl;
      continue;
    }
  }
}

// Function to extract sea point indices
template <typename DT>
std::vector<std::size_t> getSeaPointIndices(const Grid<DT> &grid) {
  std::vector<std::size_t> sea_point_indices;
  sea_point_indices.reserve(grid.num_points / 2);

  for (std::size_t idx = 0; idx < grid.num_points; idx++) {
    if (grid.in_the_sea[idx]) {
      sea_point_indices.push_back(idx);
    }
  }
  return sea_point_indices;
}

template <typename TT>
std::vector<TT> generateExpectedTimeSteps(TT start_time, TT end_time,
                                          TT time_step_interval) {
  std::vector<TT> expected_time_steps;
  for (TT t = start_time; t <= end_time; t += time_step_interval) {
    expected_time_steps.push_back(t);
  }
  return expected_time_steps;
}

/**
 * @brief Collect time series for sea points across a set of files.
 *
 * Builds variable_time_series[var][sea_point_i][t_index].
 * The time axis is generated from earliest_time..latest_time with the provided
 * interval.
 *
 * @tparam DataType numeric type.
 * @tparam TimeType time type (should be integer-like seconds in the pipeline).
 *
 * @throws std::runtime_error if no time steps are found.
 */
template <typename DT, typename TT>
void collectTimeSeriesData(
    const Grid<DT> &grid, const std::vector<std::size_t> &sea_point_indices,
    const std::vector<std::string> &variable_names,
    const std::vector<std::string> &file_paths,
    std::map<std::string, std::vector<std::vector<DT>>> &variable_time_series,
    std::vector<TT> &expected_time_steps, TT time_step_interval) {
  (void)grid;

  // Determine the start and end times
  TT earliest_time = std::numeric_limits<TT>::max();
  TT latest_time = std::numeric_limits<TT>::lowest();
  std::map<TT, std::string> time_to_file_map;

  // Collect all available time steps from the files
  std::set<TT> actual_time_steps_set;

  for (const auto &file_path : file_paths) {
    try {
      NetCDFProcessor<DT, TT> processor(file_path, variable_names);

      // Get time values
      const auto &time_values = processor.getTimeValues();

      for (const auto &time_step : time_values) {
        actual_time_steps_set.insert(time_step);
        time_to_file_map[time_step] = file_path;
        earliest_time = std::min(earliest_time, time_step);
        latest_time = std::max(latest_time, time_step);
      }
    } catch (const std::exception &e) {
      std::cerr << "Error processing data from file '" << file_path
                << "': " << e.what() << std::endl;
      continue;
    }
  }

  if (actual_time_steps_set.empty()) {
    throw std::runtime_error("No time steps found in data files.");
  }

#ifdef DATA_PREP_DEBUG_TIME_STEPS
  std::cerr << "\n=== Actual time steps found in data files ===\n";
  std::cerr << "Earliest: " << earliest_time << "\n";
  std::cerr << "Latest:   " << latest_time << "\n";
  std::cerr << "All actual time steps:\n";
  for (auto t : actual_time_steps_set)
    std::cerr << t << " ";
  std::cerr << "\n\n";
#endif

  // Generate expected time steps
  expected_time_steps =
      generateExpectedTimeSteps(earliest_time, latest_time, time_step_interval);
  std::size_t num_time_steps = expected_time_steps.size();

#ifdef DATA_PREP_DEBUG_TIME_STEPS
  std::cerr << "=== Expected time steps (" << num_time_steps << ") ===\n";
#endif

  // Initialize time series with NaNs
  for (const auto &var_name : variable_names) {
    variable_time_series[var_name].resize(sea_point_indices.size());
    for (auto &ts_data : variable_time_series[var_name]) {
      ts_data.resize(num_time_steps, std::numeric_limits<DT>::quiet_NaN());
    }
  }

  // Map actual time steps to indices in expected time steps
  std::map<TT, std::size_t> time_step_indices;
  for (std::size_t idx = 0; idx < expected_time_steps.size(); idx++) {
    time_step_indices[expected_time_steps[idx]] = idx;
  }

  // Load data and place it in the correct position in time series
  for (const auto &[time_step, file_path] : time_to_file_map) {
    try {
      NetCDFProcessor<DT, TT> processor(file_path, variable_names);

      // Load data for this time step
      std::map<std::string, std::vector<DT>> data_at_time;

      for (const auto &var_name : variable_names) {
        const auto data_vector =
            processor.getVariableDataAtTime(var_name, time_step);

        // Collect data for sea points
        std::vector<DT> sea_point_data;
        sea_point_data.reserve(sea_point_indices.size());
        for (auto idx : sea_point_indices) {
          sea_point_data.push_back(data_vector[idx]);
        }
        data_at_time[var_name] = std::move(sea_point_data);
      }

      // Find the index in expected_time_steps
      auto it = time_step_indices.find(time_step);
      if (it != time_step_indices.end()) {
        std::size_t t_idx = it->second;

        // Place data into time series
        for (const auto &var_name : variable_names) {
          const auto &var_data = data_at_time.at(var_name);
          auto &time_series = variable_time_series[var_name];

          for (std::size_t i = 0; i < sea_point_indices.size(); i++) {
            time_series[i][t_idx] = var_data[i];
          }
        }
      } else {
        std::cerr << "Warning: Time step " << time_step
                  << " not in expected time steps." << std::endl;
      }
    } catch (const std::exception &e) {
      std::cerr << "Error processing data from file '" << file_path
                << "': " << e.what() << std::endl;
      continue;
    }
  }
}

template <typename DT, typename TT = double>
void interpolateMissingData(std::vector<DT> &time_series,
                            const std::vector<TT> &time_points) {
  std::size_t n = time_series.size();
  if (n != time_points.size()) {
    throw std::runtime_error("Time series and time points size mismatch.");
  }

  std::size_t start_idx = 0;
  while (start_idx < n) {
    // Skip over existing data
    while (start_idx < n && !std::isnan(time_series[start_idx])) {
      start_idx++;
    }

    if (start_idx == n)
      break;

    // Find the end of the NaN segment
    std::size_t end_idx = start_idx;
    while (end_idx < n && std::isnan(time_series[end_idx])) {
      end_idx++;
    }

    // Handle edge cases
    if (start_idx == 0 || end_idx == n) {
      // Cannot interpolate, decide how to handle
      DT fill_value;
      if (start_idx == 0 && end_idx < n) {
        // Missing data at the beginning, fill with first non-NaN value
        fill_value = time_series[end_idx];
      } else if (end_idx == n && start_idx > 0) {
        // Missing data at the end, fill with last non-NaN value
        fill_value = time_series[start_idx - 1];
      } else {
        // Entire time series is NaN
        fill_value = static_cast<DT>(0);
      }

      for (std::size_t idx = start_idx; idx < end_idx; idx++) {
        time_series[idx] = fill_value;
      }
    } else {
      // Interpolate between existing data
      DT y0 = time_series[start_idx - 1];
      DT y1 = time_series[end_idx];
      TT x0 = time_points[start_idx - 1];
      TT x1 = time_points[end_idx];

      for (std::size_t idx = start_idx; idx < end_idx; idx++) {
        TT xi = time_points[idx];
        DT yi = y0 + ((y1 - y0) * (xi - x0) / (x1 - x0));
        time_series[idx] = yi;
      }
    }

    start_idx = end_idx;
  }
}

template <typename DT, typename TT = double>
void createSplinesForSeaPoints(
    std::map<std::string, std::vector<std::vector<DT>>> &&variable_time_series,
    std::map<std::string, std::vector<TimeSpline<DT, TT>>> &variable_splines,
    TT start_time, TT time_step) {
  for (auto &[var_name, time_series_data] : variable_time_series) {
    std::vector<TimeSpline<DT, TT>> &splines = variable_splines[var_name];
    splines.resize(time_series_data.size());

    // each i touches a unique element -> safe to move in parallel
    oneapi::tbb::parallel_for(
        std::size_t(0), time_series_data.size(), [&](std::size_t i) {
          auto &data = time_series_data[i];
          TimeSpline<DT, TT> spline;
          spline.values = std::move(data); // <-- move, no copy
          spline.createSpline(start_time, time_step);
          splines[i] = std::move(spline);
        });

    // optional: release the now-empty storage for this variable
    time_series_data.clear();
    time_series_data.shrink_to_fit();
  }

  // optional: free the whole (now moved-from) map early
  variable_time_series.clear();
  variable_time_series = {};
}

// ---------------------------------------------------------------------------
// helper: zero-pad an integer to 2 characters
// ---------------------------------------------------------------------------
inline std::string formatTwo(int v) {
  std::ostringstream oss;
  oss << std::setw(2) << std::setfill('0') << v;
  return oss.str();
}

// Helper function: Convert a Date into a string formatted as "YYYY_MM_DD".
inline std::string formatDate(Date const &d) {
  std::ostringstream oss;
  oss << std::setw(4) << std::setfill('0') << d.getYear() << "_" << std::setw(2)
      << std::setfill('0') << d.getMonth() << "_" << std::setw(2)
      << std::setfill('0') << d.getDay();
  return oss.str();
}

// -----------------------------------------------------------------------------
//  Utility: parse "YYYY-MM-DD"  or "YYYY_MM_DD" into Date                    |
// -----------------------------------------------------------------------------
inline Date parseDate(const std::string &s) {
  std::regex re(R"((\d{4})[-_](\d{2})[-_](\d{2}))");
  std::smatch m;
  if (!std::regex_match(s, m, re)) {
    std::cerr << "Invalid date format: " << s << " (expected YYYY-MM-DD)\n";
    std::exit(EXIT_FAILURE);
  }
  int y = std::stoi(m[1]), mth = std::stoi(m[2]), d = std::stoi(m[3]);
  return Date(d, mth, y);
}

// comma-separated list of args to vector<string>
inline std::vector<std::string> splitArgs(const std::string &s) {
  std::vector<std::string> v;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (!item.empty())
      v.push_back(item);
  }
  return v;
}

// ---------------------------------------------------------------------------
// helper: expand pattern for one date
// ---------------------------------------------------------------------------
inline std::string expandPattern(const std::string &pattern, const Date &d,
                                 int level) {
  std::ostringstream lev_ss;
  lev_ss << std::setw(4) << std::setfill('0') << level;
  std::string out = pattern;

  const std::pair<std::regex, std::string> subs[] = {
      {std::regex("\\{YYYY\\}"), std::to_string(d.getYear())},
      {std::regex("\\{MM\\}"), formatTwo(d.getMonth())},
      {std::regex("\\{DD\\}"), formatTwo(d.getDay())},
      {std::regex("\\{LEV\\}"), lev_ss.str()}};
  for (auto const &s : subs)
    out = std::regex_replace(out, s.first, s.second);

  return out;
}

// ---------------------------------------------------------------------------
// build path vectors for both sub-domains
// ---------------------------------------------------------------------------
inline void buildPathVectors(const DateRange &dr, int level,
                             const std::string &path_1,
                             const std::string &path_2,
                             std::vector<std::string> &out_1,
                             std::vector<std::string> &out_2) {
  out_1.clear();
  out_2.clear();
  out_1.reserve(dr.length());
  out_2.reserve(dr.length());

  for (const Date &d : dr) {
    out_1.push_back(expandPattern(path_1, d, level));
    out_2.push_back(expandPattern(path_2, d, level));
  }
}

using DataType = float;
using TimeType = double;

// -----------------------------------------------------------------------------
//  Spline containers |
// -----------------------------------------------------------------------------

template <class DT> using spline_vec = std::vector<TimeSpline<DT, TimeType>>;

template <class DT> using spline_map = std::map<std::string, spline_vec<DT>>;

struct MonthData // lives one month, then discarded
{
  Date month{1, 1, 2000};     // default so MonthData() works
  spline_map<DataType> spl_1; // sub-domain 1  (“left half”)
  spline_map<DataType> spl_2; // sub-domain 2  (“right half”)
  spline_map<DataType> spl;
  Grid<DataType> grid_1, grid_2;
  std::vector<std::size_t> sea_idx_1, sea_idx_2;
  std::vector<DataType> sea_x, sea_y;
  PointCloud<DataType> sea_cloud, land_cloud;
  std::unique_ptr<KDTree<DataType>> kd_sea{nullptr};
  std::unique_ptr<KDTree<DataType>> kd_land{nullptr};
};

// --------------------------------------------------------------------
// Helper: Join two sets of splines for each variable
// --------------------------------------------------------------------
template <typename DT>
void joinSplines(
    std::map<std::string, std::vector<TimeSpline<DT, TimeType>>> &&splines_1,
    std::map<std::string, std::vector<TimeSpline<DT, TimeType>>> &splines_2,
    std::map<std::string, std::vector<TimeSpline<DT, TimeType>>>
        &joined_splines) {
  joined_splines = std::move(splines_1); // take ownership

  for (auto &kv : splines_2) {
    auto &vec = joined_splines[kv.first];
    vec.insert(vec.end(), std::make_move_iterator(kv.second.begin()),
               std::make_move_iterator(kv.second.end())); // move, not copy
  }
  for (auto &kv : splines_2)
    kv.second.shrink_to_fit();
  splines_2.clear(); // optional: drop map nodes
}

inline void dumpCoastCSV(const MonthData &m, const std::string &path,
                         int stride = 8) {
  std::ofstream out(path);
  if (!out)
    return;
  out << std::setprecision(10);
  // header
  out << "lon_deg,lat_deg,label\n";
  // land points (used for beaching)
  for (size_t i = 0; i + 1 < m.land_cloud.pts.size(); i += 2 * stride) {
    double lon_deg = m.land_cloud.pts[i] * 180. / number_pi;
    double lat_deg = m.land_cloud.pts[i + 1] * 180. / number_pi;
    out << lon_deg << "," << lat_deg << ",land\n";
  }
  // (optional) sea points
  // for (size_t i=0;i+1<m.sea_cloud.pts.size();i+=2*stride)
  //{
  //    double lon_deg=m.sea_cloud.pts[i]  *180./pi;
  //    double lat_deg=m.sea_cloud.pts[i+1]*180./pi;
  //    out<<lon_deg<<","<<lat_deg<<",sea\n";
  //}
}

inline void buildKdTrees(MonthData &m) {
  auto kd_params = nanoflann::KDTreeSingleIndexAdaptorParams(10);
  m.kd_sea = std::make_unique<KDTree<DataType>>(2, m.sea_cloud, kd_params);
  m.kd_land = std::make_unique<KDTree<DataType>>(2, m.land_cloud, kd_params);
  m.kd_sea->buildIndex();
  m.kd_land->buildIndex();
  assert(m.kd_sea && m.kd_land);
  assert(m.sea_cloud.kdtree_get_point_count() > 0);
  assert(m.land_cloud.kdtree_get_point_count() > 0);
}

inline MonthData loadOneMonth(Date month, // 2000-01-01, 2000-02-01 ...
                              int level, const std::string &pattern_1,
                              const std::string &pattern_2,
                              const std::vector<std::string> &grid_vars,
                              const std::vector<std::string> &data_vars,
                              TimeType dt_seconds) {
  // ---------- 1. daily file list for *this* month (+ padding) ------------
  Date d_0 = month.addDays(-1);             // day before 1st
  Date d_1 = month.addMonths(1).addDays(1); // day after last
  DateRange span(d_0, d_1);

  std::vector<std::string> fp_1, fp_2;
  buildPathVectors(span, level, pattern_1, pattern_2, fp_1, fp_2);

  // ---------- 2. minimal grid init (first file of each sub-domain) -------
  MonthData m;
  m.month = month;
  m.grid_1 = initializeGrid<DataType>(fp_1.front(), grid_vars);
  m.grid_2 = initializeGrid<DataType>(fp_2.front(), grid_vars);

  auto norm_pi_inplace = [](std::vector<DataType> &a) {
    for (auto &L : a) {
      while (L <= -number_pi)
        L += 2 * number_pi;
      while (L > number_pi)
        L -= 2 * number_pi;
    }
  };

  norm_pi_inplace(m.grid_1.lon_rad);
  norm_pi_inplace(m.grid_2.lon_rad);

  m.sea_idx_1 = getSeaPointIndices(m.grid_1);
  m.sea_idx_2 = getSeaPointIndices(m.grid_2);

  // ---------- 3. Build sea & land clouds + KD-trees for this month */
  auto buildClouds = [&](MonthData &m) {
    const std::size_t n_1 = m.grid_1.lon_rad.size();
    const std::size_t n_2 = m.grid_2.lon_rad.size();
    const std::size_t total = n_1 + n_2;

    // 1) SEA CLOUD  (order must match spline vectors: spl_1 then spl_2)
    m.sea_cloud.pts.reserve((m.sea_idx_1.size() + m.sea_idx_2.size()) * 2);

    // domain 1
    for (auto idx : m.sea_idx_1) {
      m.sea_cloud.pts.push_back(m.grid_1.lon_rad[idx]);
      m.sea_cloud.pts.push_back(m.grid_1.lat_rad[idx]);
    }
    // domain 2
    for (auto idx : m.sea_idx_2) {
      m.sea_cloud.pts.push_back(m.grid_2.lon_rad[idx]);
      m.sea_cloud.pts.push_back(m.grid_2.lat_rad[idx]);
    }

    m.sea_x.clear();
    m.sea_y.clear();
    m.sea_x.reserve(m.sea_cloud.kdtree_get_point_count());
    m.sea_y.reserve(m.sea_cloud.kdtree_get_point_count());

    for (size_t i = 0; i < m.sea_cloud.kdtree_get_point_count(); i++) {
      m.sea_x.push_back(m.sea_cloud.kdtree_get_pt(i, 0));
      m.sea_y.push_back(m.sea_cloud.kdtree_get_pt(i, 1));
    }

    // 2) LAND CLOUD (anything not in sea_idx_*; order doesn't matter)
    std::vector<char> is_sea(total, 0);
    for (auto i : m.sea_idx_1)
      is_sea[i] = 1;
    for (auto i : m.sea_idx_2)
      is_sea[i + n_1] = 1;

    m.land_cloud.pts.reserve((total - m.sea_idx_1.size() - m.sea_idx_2.size()) *
                             2);

    for (std::size_t idx = 0; idx < total; idx++) {
      if (is_sea[idx])
        continue;

      DataType lon, lat;
      if (idx < n_1) {
        lon = m.grid_1.lon_rad[idx];
        lat = m.grid_1.lat_rad[idx];
      } else {
        std::size_t jdx = idx - n_1;
        lon = m.grid_2.lon_rad[jdx];
        lat = m.grid_2.lat_rad[jdx];
      }
      m.land_cloud.pts.push_back(lon);
      m.land_cloud.pts.push_back(lat);
    }

    const auto sea_N = m.sea_cloud.kdtree_get_point_count();
    const auto land_N = m.land_cloud.kdtree_get_point_count();

    if (sea_N == 0) {
      throw std::runtime_error("No SEA points for month " +
                               formatDate(m.month) +
                               " (check masks/variables/region).");
    }
    if (land_N == 0) {
      throw std::runtime_error("No LAND points for month " +
                               formatDate(m.month) +
                               " (check masks/variables/region).");
    }

    // 3) KD trees
    // auto kd_params=nanoflann::KDTreeSingleIndexAdaptorParams(10);
    // m.kd_sea =std::make_unique<KDTree<DataType>>(2,m.sea_cloud ,kd_params);
    // m.kd_land=std::make_unique<KDTree<DataType>>(2,m.land_cloud,kd_params);
    // m.kd_sea ->buildIndex();
    // m.kd_land->buildIndex();
  };

  buildClouds(m);

  dumpCoastCSV(m, "coast_" + formatDate(m.month) + ".csv", 8);

  // ---------- 3. read time series + build splines ------------------------
  std::vector<TimeType> time_steps_1, time_steps_2;
  // spline_map<DataType> tmp_spl_1,tmp_spl_2;
  std::map<std::string, std::vector<std::vector<DataType>>> raw_time_series_1,
      raw_time_series_2;

  collectTimeSeriesData(m.grid_1, m.sea_idx_1, data_vars, fp_1,
                        raw_time_series_1, time_steps_1, dt_seconds);
  collectTimeSeriesData(m.grid_2, m.sea_idx_2, data_vars, fp_2,
                        raw_time_series_2, time_steps_2, dt_seconds);

  if (time_steps_1.empty() || time_steps_2.empty()) {
    throw std::runtime_error(
        "No time steps found for month " + formatDate(m.month) +
        " (check file patterns / level / variable names).");
  }

  /* fill gaps, build splines */
  auto patch = [&](auto &ts, auto &steps) {
    for (auto &v : data_vars)
      for (auto &vec : ts[v])
        interpolateMissingData<DataType>(vec, steps);
  };
  patch(raw_time_series_1, time_steps_1);
  patch(raw_time_series_2, time_steps_2);

  createSplinesForSeaPoints(std::move(raw_time_series_1), m.spl_1,
                            static_cast<TimeType>(time_steps_1.front()),
                            dt_seconds);
  createSplinesForSeaPoints(std::move(raw_time_series_2), m.spl_2,
                            static_cast<TimeType>(time_steps_2.front()),
                            dt_seconds);

  m.spl.clear();
  joinSplines<DataType>(std::move(m.spl_1), m.spl_2, m.spl);
  m.spl_2.clear(); // we moved spl_1; spl_2 was copied into m.spl
  m.spl_2 = {};    // drop its map/vector overhead

  // Optional sanity check (debug builds)
  for (auto &kv : m.spl)
    assert(kv.second.size() == m.sea_cloud.kdtree_get_point_count());

  return m;
}

// ---------------------------------------------------------------
// Build a NeighborsData object for one month (uses that month’s KD-sea)
// ---------------------------------------------------------------
inline NeighborsData<DataType, TimeType>
makeNeighborsData(const MonthData &m, DataType search_radius_rad,
                  DataType shape_param) {
  // const size_t N=m.sea_cloud.kdtree_get_point_count();

  DataType radius_squared = search_radius_rad * search_radius_rad;
  return NeighborsData<DataType, TimeType>(*m.kd_sea, m.sea_x, m.sea_y,
                                           radius_squared, shape_param);
}

// seconds since 2000-01-01 00:00
inline TimeType secondsSince2000(const Date &d) {
  return static_cast<TimeType>(d.secondsSince(Date(1, 1, 2000)));
}

inline Date dateFromSecondsSince2000(TimeType sec) {
  // brute force
  Date ref(1, 1, 2000);
  long long days = static_cast<long long>(
      std::floor(static_cast<long double>(sec) / 86400.L));
  return ref.addDays(static_cast<int>(days));
}

enum class EndReason : uint8_t {
  none = 0,
  landed = 1,
  neverleft = 2,
  lifetime = 3
};

// --------------------------------------------------------------------
// Structure to store final trajectory results for a particle
// --------------------------------------------------------------------
template <typename DT, typename TT> struct Particle {
  DT init_lon;  // initial longitude (radians)
  DT init_lat;  // initial latitude (radians)
  TT init_time; // initial time (seconds since 01.01.2000 0:00)
  DT lon;       // final longitude (radians)
  DT lat;       // final latitude (radians)
  TT time;      // final time (seconds since 01.01.2000 0:00)
  EndReason reason;
};

} // namespace particle_tracker

#endif // DATA_PREP_UTILS_HPP
