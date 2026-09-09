#pragma once
/**
 * @file interpolation_utils.hpp
 * @brief Time and spherical (lon/lat) interpolation utilities, including an
 * experimental "neighbors-of-neighbors" spatial derivative estimator.
 *
 * ## What this header provides
 * - **Time interpolation** of scalar fields using a cubic B-spline (Boost).
 * - **Spatial neighbor queries** on a sphere (lon/lat) via a KD-tree
 * (nanoflann) using a fast approximate metric suitable for KD-tree pruning.
 * - **Spatial interpolation** using Gaussian weights derived from *true*
 * great-circle distance.
 * - **Experimental** spatial derivative estimation using a two-level
 * neighborhood
 *   ("neighbors of neighbors") and distance-derivative information.
 *
 * ## Coordinate conventions
 * - Coordinates are in radians:
 *   - longitude `lon` in (-pi, pi]
 *   - latitude  `lat` in [-pi/2, pi/2]
 *
 * ## Dependencies
 * - Boost (either:
 *   - `<boost/math/interpolators/cardinal_cubic_b_spline.hpp>` (newer)
 *   - or `<boost/math/interpolators/cubic_b_spline.hpp>` (older / deprecated in
 * some newer Boost)
 * - nanoflann (header-only) for KD-tree support.
 *
 * ## Error policy
 * This header prefers **loud failures** over silent bad behavior:
 * - If enable nanoflann support and the header cannot be found: **compile-time
 * error**.
 * - If disable nanoflann support but try to instantiate KD-tree dependent
 * types: **compile-time error**.
 *
 * ## Threading
 * - No internal synchronization is provided.
 * - Reads are safe if underlying data structures are immutable.
 * - The caller must protect shared mutable structures.
 */

#include <array>
#include <cmath>   // std::cos, std::exp
#include <cstddef> // std::size_t
#include <map>
#include <numbers>   // std::numbers::pi_v
#include <stdexcept> // std::runtime_error, std::out_of_range
#include <string>
#include <type_traits> // std::false_type
#include <vector>

#include <particle_tracker/geography_utils.hpp>

//------------------------------------------------------------------------------
// User-configurable feature switches
//------------------------------------------------------------------------------

/**
 * @def PARTICLE_TRACKER_HAS_NANOFLANN
 * @brief Set to 1 to enable nanoflann-backed KD-tree functionality.
 *
 * If set to 0, KD-tree dependent classes/types in this header become
 * unavailable and will produce a compile-time error if instantiated.
 */
#ifndef PARTICLE_TRACKER_HAS_NANOFLANN
#define PARTICLE_TRACKER_HAS_NANOFLANN 1
#endif

//------------------------------------------------------------------------------
// Boost spline selection (new vs old header layouts)
//------------------------------------------------------------------------------
//
// Newer Boost provides interpolators::cardinal_cubic_b_spline.
// Older Boost provides cubic_b_spline.
// We detect availability by header presence and define consistent macros.
//------------------------------------------------------------------------------

#if __has_include(<boost/math/interpolators/cardinal_cubic_b_spline.hpp>)
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#define PARTICLE_TRACKER_HAS_BOOST_CARDINAL_CUBIC_B_SPLINE 1
#define PARTICLE_TRACKER_HAS_BOOST_CUBIC_B_SPLINE 0
#elif __has_include(<boost/math/interpolators/cubic_b_spline.hpp>)
#include <boost/math/interpolators/cubic_b_spline.hpp>
#define PARTICLE_TRACKER_HAS_BOOST_CARDINAL_CUBIC_B_SPLINE 0
#define PARTICLE_TRACKER_HAS_BOOST_CUBIC_B_SPLINE 1
#else
#error                                                                         \
    "particle_tracker: Boost spline header not found. Need either <boost/math/interpolators/cardinal_cubic_b_spline.hpp> or <boost/math/interpolators/cubic_b_spline.hpp>."
#endif

//------------------------------------------------------------------------------
// nanoflann detection and include (only if enabled)
//------------------------------------------------------------------------------
//
// Supported include layouts:
//  - <nanoflann.hpp>                                  (system installs)
//  - <nanoflann/nanoflann.hpp>                        (some packagers)
//  - <third_party/nanoflann/nanoflann.hpp> (vendored under include/third_party)
//------------------------------------------------------------------------------

#if PARTICLE_TRACKER_HAS_NANOFLANN
#if __has_include(<nanoflann.hpp>)
#include <nanoflann.hpp>
#elif __has_include(<nanoflann/nanoflann.hpp>)
#include <nanoflann/nanoflann.hpp>
#elif __has_include(<third_party/nanoflann/nanoflann.hpp>)
#include <third_party/nanoflann/nanoflann.hpp>
#else
#error                                                                         \
    "particle_tracker: nanoflann enabled (PARTICLE_TRACKER_HAS_NANOFLANN=1) but nanoflann header not found. Install nanoflann or vendor it under include/third_party/nanoflann and include it as <third_party/nanoflann/nanoflann.hpp>."
#endif
#endif

namespace particle_tracker {

#if PARTICLE_TRACKER_HAS_NANOFLANN
namespace nf = nanoflann;
#endif

//------------------------------------------------------------------------------
// Constants and helpers
//------------------------------------------------------------------------------

/** @brief pi in double precision (radians). */
inline constexpr double pt_pi = std::numbers::pi_v<double>;
/** @brief 2*pi in double precision (radians). */
inline constexpr double pt_two_pi = 2. * pt_pi;

/**
 * @brief Wrap longitude into (-pi, pi].
 * @tparam T Floating-point type.
 * @param lon Longitude in radians.
 * @return Wrapped longitude in (-pi, pi].
 *
 * This is used to avoid discontinuities across the dateline.
 */
template <typename T> inline T pt_wrap_lon(T lon) {
  // Use pt_pi constants but cast to T so comparisons behave correctly.
  const T pi = static_cast<T>(pt_pi);
  const T two_pi = static_cast<T>(pt_two_pi);

  while (lon <= -pi)
    lon += two_pi;
  while (lon > pi)
    lon -= two_pi;
  return lon;
}

/**
 * @brief 2D point type used for (lon, lat).
 * @tparam DataType numeric type
 */
template <typename DataType> using point = std::array<DataType, 2>;

/**
 * @brief Default time spacing used in some HYCOM-style datasets (3 hours).
 *
 * This constant is not required by the classes but is useful as a conventional
 * default.
 */
inline constexpr double dtt = 10800.0;

//------------------------------------------------------------------------------
// TimeSpline
//------------------------------------------------------------------------------

/**
 * @class TimeSpline
 * @brief Wrapper around a Boost cubic spline interpolator for uniformly spaced
 * samples.
 *
 * @tparam DataType Value type of the samples (e.g. float, double).
 * @tparam TimeType Time type used for evaluation (typically double).
 *
 * ### Usage
 * - Fill @ref values with uniformly spaced samples.
 * - Call @ref createSpline(start_time, time_step).
 * - Evaluate with `spline(t)` and time-derivative with `spline.prime(t)` if
 * supported.
 *
 * ### Notes
 * - This class stores the data and the spline object.
 * - The underlying Boost spline type depends on what Boost provides:
 *   - `boost::math::interpolators::cardinal_cubic_b_spline` (preferred when
 * available)
 *   - `boost::math::cubic_b_spline` (fallback)
 */
template <typename DataType, typename TimeType = double> class TimeSpline {
public:
  /// Sample values at uniform time spacing.
  std::vector<DataType> values;

#if PARTICLE_TRACKER_HAS_BOOST_CUBIC_B_SPLINE
  /// Spline implementation from older Boost.
  boost::math::cubic_b_spline<DataType> spline;
#elif PARTICLE_TRACKER_HAS_BOOST_CARDINAL_CUBIC_B_SPLINE
  /// Spline implementation from newer Boost (interpolators namespace).
  boost::math::interpolators::cardinal_cubic_b_spline<DataType> spline;
#endif

  TimeSpline() = default;

  /**
   * @brief Build the spline from the current @ref values.
   * @param start_time Time of the first sample.
   * @param time_step  Uniform time step between samples.
   *
   * @throws std::runtime_error if @ref values is empty.
   */
  void createSpline(TimeType start_time, TimeType time_step) {
    if (values.empty())
      throw std::runtime_error("particle_tracker::TimeSpline: no data to "
                               "create spline (values is empty).");

#if PARTICLE_TRACKER_HAS_BOOST_CUBIC_B_SPLINE
    // Old Boost API accepts pointer + size + (t0, dt)
    spline = boost::math::cubic_b_spline<DataType>(
        values.data(), values.size(), static_cast<DataType>(start_time),
        static_cast<DataType>(time_step));
#elif PARTICLE_TRACKER_HAS_BOOST_CARDINAL_CUBIC_B_SPLINE
    // This Boost version accepts pointer + size + (t0, dt)
    spline = boost::math::interpolators::cardinal_cubic_b_spline<DataType>(
        values.data(), values.size(), static_cast<DataType>(start_time),
        static_cast<DataType>(time_step));
#endif
  }
};

//------------------------------------------------------------------------------
// nanoflann dataset adaptor + distance adaptor + KDTree alias
//------------------------------------------------------------------------------

#if PARTICLE_TRACKER_HAS_NANOFLANN

/**
 * @struct PointCloud
 * @brief Minimal dataset adaptor for nanoflann storing (lon, lat) in a flat
 * array.
 *
 * @tparam DataType numeric type (float/double).
 *
 * Points are stored as:
 * ```
 * pts = [lon0, lat0, lon1, lat1, ..., lon(N-1), lat(N-1)]
 * ```
 *
 * @warning `pts.size()` must be even (2 values per point).
 */
template <typename DataType> struct PointCloud {
  /// Flat storage: [lon0, lat0, lon1, lat1, ...]. Must have even length.
  std::vector<DataType> pts;

  /// @return number of points stored (pts.size()/2).
  inline std::size_t kdtree_get_point_count() const { return pts.size() / 2; }

  /**
   * @brief Coordinate accessor for nanoflann.
   * @param idx Point index [0, point_count)
   * @param dim Coordinate index (0 -> lon, 1 -> lat)
   * @return coordinate value
   */
  inline DataType kdtree_get_pt(const std::size_t idx,
                                const std::size_t dim) const {
    return pts[2 * idx + dim];
  }

  /// Direct access to raw contiguous data.
  inline const DataType *data() const { return pts.data(); }

  /**
   * @brief Bounding box hook for nanoflann.
   * @return false (no precomputed bounding box is provided)
   */
  template <class BBOX> bool kdtree_get_bbox(BBOX &) const { return false; }
};

/**
 * @struct LonLatDistanceAdaptor
 * @brief Approximate metric for KD-tree pruning in (lon, lat).
 *
 * @tparam DataType numeric type.
 * @tparam DatasetAdaptor type providing `pts` in [lon,lat,lon,lat,...] layout.
 *
 * ### What this metric is (and is not)
 * - This is **not** a great-circle distance.
 * - It is a cheap local approximation used *only* for KD-tree pruning:
 *   ```
 *   d^2 ~= (dlat)^2 + (cos(lat_avg) * dlon_wrapped)^2
 *   ```
 * - After KD-tree returns candidate neighbors, **true** great-circle distance
 * is computed via `GCD()` / `GCD_deriv()` and used to compute weights.
 *
 * ### Dateline behavior
 * - Longitude differences are wrapped into (-pi, pi] for continuity across the
 * dateline.
 */
template <typename DataType, class DatasetAdaptor>
struct LonLatDistanceAdaptor {
  typedef DataType ElementType;
  typedef DataType DistanceType;

  const DatasetAdaptor &data;

  explicit LonLatDistanceAdaptor(const DatasetAdaptor &data_) : data(data_) {}

  /// @return pointer to (lon, lat) for dataset point idx.
  inline const DataType *getPoint(std::size_t idx) const {
    return &data.pts[2 * idx];
  }

  /**
   * @brief Squared approximate distance between two explicit points.
   * @param a pointer to {lon, lat}
   * @param b pointer to {lon, lat}
   */
  inline DistanceType evalMetric(const DataType *a, const DataType *b,
                                 std::size_t /*size*/) const {
    DataType dlon = a[0] - b[0];
    if (dlon < -static_cast<DataType>(pt_pi))
      dlon += static_cast<DataType>(pt_two_pi);
    else if (dlon > static_cast<DataType>(pt_pi))
      dlon -= static_cast<DataType>(pt_two_pi);

    DataType dlat = a[1] - b[1];
    DataType avgLat = (a[1] + b[1]) * static_cast<DataType>(0.5);

    DataType scaled_dlon = dlon * std::cos(avgLat);
    return scaled_dlon * scaled_dlon + dlat * dlat;
  }

  /**
   * @brief Squared approximate distance between query point a and dataset point
   * idx_b.
   */
  inline DistanceType evalMetric(const DataType *a, std::size_t idx_b,
                                 std::size_t size) const {
    const DataType *b = getPoint(idx_b);
    return evalMetric(a, b, size);
  }

  /**
   * @brief Incremental squared contribution along a single dimension (KD-tree
   * pruning).
   *
   * @param a coordinate of query
   * @param b coordinate of dataset
   * @param dim dimension index (0 -> lon, 1 -> lat)
   */
  template <typename T, typename IndexType>
  inline DistanceType accum_dist(const T a, const T b, IndexType dim) const {
    DistanceType diff = static_cast<DistanceType>(a - b);

    if (dim == 0) {
      if (diff < -static_cast<DistanceType>(pt_pi))
        diff += static_cast<DistanceType>(pt_two_pi);
      else if (diff > static_cast<DistanceType>(pt_pi))
        diff -= static_cast<DistanceType>(pt_two_pi);
    }

    return diff * diff;
  }
};

/**
 * @brief KD-tree alias used throughout this header when nanoflann is enabled.
 */
template <typename DataType>
using KDTree = nf::KDTreeSingleIndexAdaptor<
    LonLatDistanceAdaptor<DataType, PointCloud<DataType>>, PointCloud<DataType>,
    2>;

#else // PARTICLE_TRACKER_HAS_NANOFLANN == 0

// Provide "poison pill" templates that fail loudly if instantiated.

template <typename...> struct pt_dependent_false : std::false_type {};

/**
 * @brief Dummy KDTree type when nanoflann is disabled.
 *
 * Any attempt to instantiate NeighborsData / NeighborWithNeighborsData with
 * this will fail at compile time with a clear error.
 */
template <typename DataType> struct KDTree {
  static_assert(pt_dependent_false<DataType>::value,
                "particle_tracker: KDTree is unavailable because "
                "PARTICLE_TRACKER_HAS_NANOFLANN=0. "
                "Enable nanoflann support or avoid KD-tree dependent APIs.");
};

#endif // PARTICLE_TRACKER_HAS_NANOFLANN

//------------------------------------------------------------------------------
// Neighbor and NeighborsData
//------------------------------------------------------------------------------

/**
 * @struct Neighbor
 * @brief One spatial neighbor (grid node) used in interpolation.
 *
 * @tparam DataType numeric type.
 *
 * Fields:
 * - @ref ind: index into grid coordinate arrays (xx_sparse/yy_sparse) and
 * per-node data arrays.
 * - @ref dist: great-circle distance (from @c GCD), **not** KD-tree metric
 * distance.
 * - @ref weight: Gaussian weight, typically `exp(-(dist*dist)/shape)`.
 *
 * @note This is a lightweight, POD-like aggregate in practice (trivially
 * copyable if DataType is), but it has a user-provided constructor, so it is
 * **not** a "POD type" in the strict C++03 sense. In modern C++, the relevant
 * traits are:
 *       - standard-layout: yes (for typical DataType)
 *       - trivially copyable: yes (for typical DataType)
 */
template <typename DataType> struct Neighbor {
  std::size_t ind{};
  DataType dist{};
  DataType weight{};

  Neighbor(std::size_t index = 0, DataType distance = DataType{},
           DataType node_weight = DataType{})
      : ind(index), dist(distance), weight(node_weight) {}
};

#if PARTICLE_TRACKER_HAS_NANOFLANN

/**
 * @class NeighborsData
 * @brief Spatial neighbor query + Gaussian weighting + multi-variable
 * interpolation.
 *
 * @tparam DataType numeric type.
 * @tparam TimeType time type.
 *
 * This class holds references to:
 * - a pre-built @ref KDTree
 * - grid longitudes and latitudes (`xx_sparse`, `yy_sparse`)
 * and parameters:
 * - `rad2`   : squared KD-tree search radius (approximate metric)
 * - `shape`  : Gaussian length scale for weights in terms of **great-circle
 * distance**
 *
 * ### Dateline policy
 * - Query longitude is wrapped into (-pi, pi] before search.
 * - Grid longitudes are wrapped before calling GCD(). The formula is
 *   2*pi-periodic in the longitude difference, so this is not needed for
 *   correctness; it keeps the difference small so that cos() is not evaluated on
 *   a large argument, where it would lose precision.
 */
template <typename DataType, typename TimeType = double> class NeighborsData {
private:
  const KDTree<DataType> &index;
  const std::vector<DataType> &xx_sparse{};
  const std::vector<DataType> &yy_sparse{};
  const DataType rad2{};
  const DataType shape{};

public:
  NeighborsData(const KDTree<DataType> &idx, const std::vector<DataType> &xxs,
                const std::vector<DataType> &yys, DataType r2,
                DataType shape_param)
      : index(idx), xx_sparse(xxs), yy_sparse(yys), rad2(r2),
        shape(shape_param) {}

  /**
   * @brief Compute neighbor list around query point (x,y).
   *
   * - Uses KD-tree radius search with `rad2` (approx metric).
   * - Computes true great-circle distance via @c GCD for each candidate.
   * - Computes Gaussian weights from true distance.
   *
   * @param neighbors output vector cleared and filled.
   * @param x query longitude (radians)
   * @param y query latitude  (radians)
   */
  void computeNeighbors(std::vector<Neighbor<DataType>> &neighbors, DataType x,
                        DataType y) {
    x = pt_wrap_lon(x);

    nf::SearchParameters searchparams;
    searchparams.sorted = false;

    std::vector<nf::ResultItem<unsigned int, DataType>> indices_dists;

    neighbors.clear();

    point<DataType> query_point = {x, y};
    index.radiusSearch(query_point.data(), rad2, indices_dists, searchparams);

    neighbors.reserve(indices_dists.size());

    const DataType lon_q = pt_wrap_lon(query_point[0]);

    for (const auto &result : indices_dists) {
      Neighbor<DataType> nbr;
      nbr.ind = result.first;

      const DataType lon_i = pt_wrap_lon(xx_sparse.at(nbr.ind));

      nbr.dist = GCD(lon_q, query_point[1], lon_i, yy_sparse.at(nbr.ind));
      nbr.weight = std::exp(-nbr.dist * nbr.dist / shape);

      neighbors.push_back(nbr);
    }
  }

  /**
   * @brief Recompute dist/weight for an existing neighbor list without changing
   * indices.
   *
   * Use this only if the query point moved *slightly*
   * the neighbor identity may no longer be optimal. For larger moves, call
   * @ref computeNeighbors again.
   */
  void updateNeighborWeights(std::vector<Neighbor<DataType>> &neighbors,
                             DataType x, DataType y) {
    x = pt_wrap_lon(x);

    point<DataType> query_point = {x, y};
    const DataType lon_q = pt_wrap_lon(query_point[0]);

    for (auto &nbr : neighbors) {
      const DataType lon_i = pt_wrap_lon(xx_sparse.at(nbr.ind));
      nbr.dist = GCD(lon_q, query_point[1], lon_i, yy_sparse.at(nbr.ind));
      nbr.weight = std::exp(-nbr.dist * nbr.dist / shape);
    }
  }

  /**
   * @brief Weighted spatial interpolation for multiple variables at time t.
   *
   * @param t evaluation time
   * @param neighbors neighbor list produced by
   * computeNeighbors/updateNeighborWeights
   * @param variable_splines map(var_name -> spline-per-gridpoint)
   * @param variable_names variables to evaluate
   * @return map(var_name -> interpolated value)
   */
  std::map<std::string, DataType> interpolateVariables(
      TimeType t, const std::vector<Neighbor<DataType>> &neighbors,
      const std::map<std::string, std::vector<TimeSpline<DataType, TimeType>>>
          &variable_splines,
      const std::vector<std::string> &variable_names) {
    std::map<std::string, DataType> out;
    DataType den = static_cast<DataType>(0);

    for (const auto &var_name : variable_names)
      out[var_name] = static_cast<DataType>(0);

    for (const auto &nbr : neighbors) {
      for (const auto &var_name : variable_names) {
        const auto &splines = variable_splines.at(var_name);

        if (nbr.ind >= splines.size()) {
          throw std::out_of_range("particle_tracker::NeighborsData: index " +
                                  std::to_string(nbr.ind) +
                                  " out of bounds for variable '" + var_name +
                                  "'.");
        }

        const DataType v = splines[nbr.ind].spline(t);
        out[var_name] += nbr.weight * v;
      }
      den += nbr.weight;
    }

    if (den > static_cast<DataType>(1e-8)) {
      for (const auto &var_name : variable_names)
        out[var_name] /= den;
    } else {
      for (const auto &var_name : variable_names)
        out[var_name] = static_cast<DataType>(0);
    }

    return out;
  }

  /**
   * @brief Convenience: computeNeighbors + interpolateVariables.
   */
  std::map<std::string, DataType> calculateVariables(
      TimeType t, std::vector<Neighbor<DataType>> &neighbors,
      const std::map<std::string, std::vector<TimeSpline<DataType, TimeType>>>
          &variable_splines,
      const std::vector<std::string> &variable_names, DataType x, DataType y) {
    computeNeighbors(neighbors, x, y);
    return interpolateVariables(t, neighbors, variable_splines, variable_names);
  }

  /**
   * @brief Convenience: updateNeighborWeights + interpolateVariables.
   */
  std::map<std::string, DataType> calculateVariablesUpdate(
      TimeType t, std::vector<Neighbor<DataType>> &neighbors,
      const std::map<std::string, std::vector<TimeSpline<DataType, TimeType>>>
          &variable_splines,
      const std::vector<std::string> &variable_names, DataType x, DataType y) {
    updateNeighborWeights(neighbors, x, y);
    return interpolateVariables(t, neighbors, variable_splines, variable_names);
  }
};

#else // PARTICLE_TRACKER_HAS_NANOFLANN == 0

template <typename DataType, typename TimeType = double> class NeighborsData {
  static_assert(pt_dependent_false<DataType, TimeType>::value,
                "particle_tracker::NeighborsData is unavailable because "
                "PARTICLE_TRACKER_HAS_NANOFLANN=0. "
                "Enable nanoflann support or avoid KD-tree dependent APIs.");
};

#endif // PARTICLE_TRACKER_HAS_NANOFLANN

//------------------------------------------------------------------------------
// Experimental derivative estimator: NeighborPrime, NeighborWithNeighbors,
// VariableDerivatives
//------------------------------------------------------------------------------

/**
 * @struct NeighborPrime
 * @brief A "secondary neighbor" (neighbor of a primary neighbor).
 *
 * @tparam DataType numeric type.
 *
 * Fields:
 * - @ref ind_prime: index of the secondary grid node.
 * - @ref dist_prime: array [d, gx, gy] returned by @c GCD_deriv(primary,
 * secondary).
 *   - dist_prime[0] = d    : great-circle distance
 *   - dist_prime[1] = gx   : east component of `grad(-d^2/2)`
 *   - dist_prime[2] = gy   : north component of `grad(-d^2/2)`
 * - @ref weight_prime: Gaussian weight `exp(-(d*d)/shape)` based on
 * dist_prime[0].
 *
 */
template <typename DataType> struct NeighborPrime {
  std::size_t ind_prime{};
  std::array<DataType, 3> dist_prime{};
  DataType weight_prime{};

  NeighborPrime(std::size_t index_prime = 0,
                std::array<DataType, 3> distance_prime = {DataType{},
                                                          DataType{},
                                                          DataType{}},
                DataType node_weight_prime = DataType{})
      : ind_prime(index_prime), dist_prime(distance_prime),
        weight_prime(node_weight_prime) {}
};

/**
 * @struct NeighborWithNeighbors
 * @brief A primary Neighbor plus a list of its secondary neighbors.
 *
 * Used for the experimental derivative estimator.
 */
template <typename DataType> struct NeighborWithNeighbors : Neighbor<DataType> {
  std::vector<NeighborPrime<DataType>> neighbors_prime;

  NeighborWithNeighbors() = default;

  NeighborWithNeighbors(std::size_t index, DataType distance,
                        DataType node_weight,
                        const std::vector<NeighborPrime<DataType>> &nbr_prime)
      : Neighbor<DataType>(index, distance, node_weight),
        neighbors_prime(nbr_prime) {}
};

/**
 * @struct VariableDerivatives
 * @brief Output container: value + time derivative + (heuristic) spatial
 * derivatives.
 *
 * @tparam DataType numeric type.
 */
template <typename DataType> struct VariableDerivatives {
  DataType value{};
  DataType time_derivative{};
  DataType spatial_derivative_x{};
  DataType spatial_derivative_y{};
};

#if PARTICLE_TRACKER_HAS_NANOFLANN

/**
 * @class NeighborWithNeighborsData
 * @brief Two-level neighborhood builder + heuristic spatial derivative
 * estimator.
 *
 * @tparam DataType numeric type.
 * @tparam TimeType time type.
 *
 * ## Two-level neighborhood
 * 1) Find primary neighbors around the query point (x,y).
 * 2) For each primary neighbor i, find secondary neighbors around grid point i.
 *
 * ## Derivative estimator (heuristic)
 * - For each primary neighbor i:
 *   - For each secondary neighbor j:
 *     - dv = v_j - v_i
 *     - accumulate using weight_prime * (gx * dv) and (gy * dv)
 * - Normalize per-primary-neighbor by sum of weight_prime (den_prime)
 * - Average over primary neighbors with primary weights
 * - Apply final scaling factor (2/shape)
 *
 * ## Dateline policy
 * - Query longitude x is wrapped into (-pi, pi].
 * - Grid longitudes are wrapped before calling @c GCD and @c GCD_deriv.
 *
 * ## Loud failure behavior
 * - If secondary neighborhood search returns empty for a primary neighbor, the
 * primary neighbor is **skipped** (because derivatives cannot be formed). This
 * changes the primary set. This is intentional and explicit; it prevents silent
 * garbage derivatives.
 */
template <typename DataType, typename TimeType = double>
class NeighborWithNeighborsData {
private:
  const KDTree<DataType> &index;
  const std::vector<DataType> &xx_sparse{};
  const std::vector<DataType> &yy_sparse{};
  const DataType rad2{};
  const DataType shape{};

public:
  NeighborWithNeighborsData(const KDTree<DataType> &idx,
                            const std::vector<DataType> &xxs,
                            const std::vector<DataType> &yys, DataType r2,
                            DataType shape_param)
      : index(idx), xx_sparse(xxs), yy_sparse(yys), rad2(r2),
        shape(shape_param) {}

  /**
   * @brief Build primary+secondary neighborhoods around (x,y).
   *
   * @param neighbors output vector cleared and filled
   * @param x query longitude (radians)
   * @param y query latitude  (radians)
   */
  void computeNeighbors(std::vector<NeighborWithNeighbors<DataType>> &neighbors,
                        DataType x, DataType y) {
    x = pt_wrap_lon(x);

    point<DataType> query_point = {x, y};

    nf::SearchParameters searchparams;
    searchparams.sorted = false;

    std::vector<nf::ResultItem<unsigned int, DataType>> indices_dists;
    index.radiusSearch(query_point.data(), rad2, indices_dists, searchparams);

    neighbors.clear();
    neighbors.reserve(indices_dists.size());

    const DataType lon_q = pt_wrap_lon(query_point[0]);

    for (const auto &result : indices_dists) {
      NeighborWithNeighbors<DataType> nbr;
      nbr.ind = result.first;

      const DataType lon_i = pt_wrap_lon(xx_sparse.at(nbr.ind));
      const DataType lat_i = yy_sparse.at(nbr.ind);

      // True great-circle distance query -> primary neighbor (for primary
      // weights).
      nbr.dist = GCD(lon_q, query_point[1], lon_i, lat_i);
      nbr.weight = std::exp(-nbr.dist * nbr.dist / shape);

      // Secondary search is centered at the primary grid node.
      point<DataType> query_point_prime = {lon_i, lat_i};

      std::vector<nf::ResultItem<unsigned int, DataType>> indices_dists_prime;
      index.radiusSearch(query_point_prime.data(), rad2, indices_dists_prime,
                         searchparams);

      if (indices_dists_prime.empty()) {
        // Loud behavior: do not keep a primary neighbor with no secondary
        // neighborhood because derivative estimation would be undefined/noisy.
        continue;
      }

      nbr.neighbors_prime.clear();
      nbr.neighbors_prime.reserve(indices_dists_prime.size());

      for (const auto &result_prime : indices_dists_prime) {
        NeighborPrime<DataType> nbr_prime;
        nbr_prime.ind_prime = result_prime.first;

        if (nbr_prime.ind_prime >= xx_sparse.size() ||
            nbr_prime.ind_prime >= yy_sparse.size()) {
          throw std::out_of_range(
              "particle_tracker::NeighborWithNeighborsData: index " +
              std::to_string(nbr_prime.ind_prime) +
              " out of bounds in computeNeighbors (secondary).");
        }

        const DataType lon_j = pt_wrap_lon(xx_sparse.at(nbr_prime.ind_prime));
        const DataType lat_j = yy_sparse.at(nbr_prime.ind_prime);

        // The primary node is part of its own radius search. Avoid calling
        // GCD_deriv for that directionless self-pair.
        if (nbr_prime.ind_prime == nbr.ind) {
          nbr_prime.dist_prime = {static_cast<DataType>(0),
                                  static_cast<DataType>(0),
                                  static_cast<DataType>(0)};
        } else {
          nbr_prime.dist_prime = GCD_deriv(query_point_prime[0],
                                           query_point_prime[1], lon_j, lat_j);
        }
        nbr_prime.weight_prime = std::exp(-nbr_prime.dist_prime[0] *
                                          nbr_prime.dist_prime[0] / shape);

        nbr.neighbors_prime.push_back(nbr_prime);
      }

      neighbors.push_back(nbr);
    }
  }

  /**
   * @brief Update only primary (query->primary) distances and weights for an
   * existing neighbor list.
   *
   * Secondary neighborhoods and their weights are not recomputed.
   */
  void
  updateNeighborWeights(std::vector<NeighborWithNeighbors<DataType>> &neighbors,
                        DataType x, DataType y) {
    x = pt_wrap_lon(x);

    point<DataType> query_point = {x, y};
    const DataType lon_q = pt_wrap_lon(query_point[0]);

    for (auto &nbr : neighbors) {
      const DataType lon_i = pt_wrap_lon(xx_sparse.at(nbr.ind));
      nbr.dist = GCD(lon_q, query_point[1], lon_i, yy_sparse.at(nbr.ind));
      nbr.weight = std::exp(-nbr.dist * nbr.dist / shape);
    }
  }

  /**
   * @brief Interpolate value/time-derivative and estimate spatial derivatives
   * for requested variables.
   *
   * @param t evaluation time
   * @param neighbors two-level neighborhood list
   * @param variable_splines map(var_name -> spline-per-gridpoint)
   * @param variable_names which variables to compute
   * @param shape normalization length scale used in final scaling (2/shape)
   * @return map(var_name -> VariableDerivatives)
   *
   * @warning The spatial derivatives are heuristic and must be validated.
   */
  std::map<std::string, VariableDerivatives<DataType>> interpolateDerivatives(
      TimeType t, const std::vector<NeighborWithNeighbors<DataType>> &neighbors,
      const std::map<std::string, std::vector<TimeSpline<DataType, TimeType>>>
          &variable_splines,
      const std::vector<std::string> &variable_names, DataType shape) const {
    // Primary accumulators
    std::map<std::string, DataType> value_acc;
    std::map<std::string, DataType> dt_acc;
    std::map<std::string, DataType> dx_acc;
    std::map<std::string, DataType> dy_acc;

    DataType den = static_cast<DataType>(0);

    for (const auto &name : variable_names) {
      value_acc[name] = static_cast<DataType>(0);
      dt_acc[name] = static_cast<DataType>(0);
      dx_acc[name] = static_cast<DataType>(0);
      dy_acc[name] = static_cast<DataType>(0);
    }

    for (const auto &nbr : neighbors) {
      const DataType w = nbr.weight;
      den += w;

      // Value + time derivative from splines at primary neighbor grid nodes
      for (const auto &name : variable_names) {
        const auto &splines = variable_splines.at(name);

        if (nbr.ind >= splines.size()) {
          throw std::out_of_range(
              "particle_tracker::NeighborWithNeighborsData: index " +
              std::to_string(nbr.ind) + " out of bounds for variable '" + name +
              "'.");
        }

        const DataType v = splines[nbr.ind].spline(t);
        const DataType dv = splines[nbr.ind].spline.prime(t);

        value_acc[name] += w * v;
        dt_acc[name] += w * dv;
      }

      // Per-primary-neighbor spatial derivative accumulators (heuristic)
      DataType den_prime = static_cast<DataType>(0);

      std::map<std::string, DataType> dx_prime;
      std::map<std::string, DataType> dy_prime;

      for (const auto &name : variable_names) {
        dx_prime[name] = static_cast<DataType>(0);
        dy_prime[name] = static_cast<DataType>(0);
      }

      for (const auto &nbr_prime : nbr.neighbors_prime) {
        // Skip self-pairs (no geometric meaning for dv)
        if (nbr_prime.ind_prime == nbr.ind)
          continue;

        const DataType gx = nbr_prime.dist_prime[1];
        const DataType gy = nbr_prime.dist_prime[2];

        // IMPORTANT: den_prime must be independent of number of variables
        // requested.
        den_prime += nbr_prime.weight_prime;

        for (const auto &name : variable_names) {
          const auto &splines = variable_splines.at(name);

          if (nbr_prime.ind_prime >= splines.size()) {
            throw std::out_of_range(
                "particle_tracker::NeighborWithNeighborsData: index " +
                std::to_string(nbr_prime.ind_prime) +
                " out of bounds for variable '" + name + "' (secondary).");
          }

          const DataType v_i = splines[nbr.ind].spline(t);
          const DataType v_j = splines[nbr_prime.ind_prime].spline(t);
          const DataType dv = v_j - v_i;

          // GCD_deriv already returns grad(-d^2/2); multiplying by d again
          // would introduce an erroneous extra distance factor.
          dx_prime[name] += nbr_prime.weight_prime * (gx * dv);
          dy_prime[name] += nbr_prime.weight_prime * (gy * dv);
        }
      }

      if (den_prime > static_cast<DataType>(1e-16)) {
        for (const auto &name : variable_names) {
          dx_acc[name] += w * (dx_prime[name] / den_prime);
          dy_acc[name] += w * (dy_prime[name] / den_prime);
        }
      }
    }

    std::map<std::string, VariableDerivatives<DataType>> out;

    if (den > static_cast<DataType>(1e-16)) {
      for (const auto &name : variable_names) {
        VariableDerivatives<DataType> vd;

        vd.value = value_acc[name] / den;
        vd.time_derivative = dt_acc[name] / den;

        vd.spatial_derivative_x =
            (static_cast<DataType>(2) / shape) * (dx_acc[name] / den);
        vd.spatial_derivative_y =
            (static_cast<DataType>(2) / shape) * (dy_acc[name] / den);

        out[name] = vd;
      }
    } else {
      // No usable primary weights -> return zeros.
      for (const auto &name : variable_names)
        out[name] = VariableDerivatives<DataType>{};
    }

    return out;
  }

  /**
   * @brief Convenience: computeNeighbors + interpolateDerivatives.
   */
  std::map<std::string, VariableDerivatives<DataType>> calculateDerivatives(
      TimeType t, std::vector<NeighborWithNeighbors<DataType>> &neighbors,
      const std::map<std::string, std::vector<TimeSpline<DataType, TimeType>>>
          &variable_splines,
      const std::vector<std::string> &variable_names, DataType x, DataType y,
      DataType shape) {
    computeNeighbors(neighbors, x, y);
    return interpolateDerivatives(t, neighbors, variable_splines,
                                  variable_names, shape);
  }

  /**
   * @brief Convenience: updateNeighborWeights + interpolateDerivatives.
   */
  std::map<std::string, VariableDerivatives<DataType>>
  calculateDerivativesUpdate(
      TimeType t, std::vector<NeighborWithNeighbors<DataType>> &neighbors,
      const std::map<std::string, std::vector<TimeSpline<DataType, TimeType>>>
          &variable_splines,
      const std::vector<std::string> &variable_names, DataType x, DataType y,
      DataType shape) {
    updateNeighborWeights(neighbors, x, y);
    return interpolateDerivatives(t, neighbors, variable_splines,
                                  variable_names, shape);
  }
};

#else // PARTICLE_TRACKER_HAS_NANOFLANN == 0

template <typename DataType, typename TimeType = double>
class NeighborWithNeighborsData {
  static_assert(pt_dependent_false<DataType, TimeType>::value,
                "particle_tracker::NeighborWithNeighborsData is unavailable "
                "because PARTICLE_TRACKER_HAS_NANOFLANN=0. "
                "Enable nanoflann support or avoid KD-tree dependent APIs.");
};

#endif // PARTICLE_TRACKER_HAS_NANOFLANN

} // namespace particle_tracker
