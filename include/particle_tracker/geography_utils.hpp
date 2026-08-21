#pragma once
/**
 * @file geography_utils.hpp
 * @brief Small geographic utilities for spherical (lon/lat) computations in radians.
 *
 * This header provides:
 * - A generic clamp helper.
 * - Great-circle distance (central angle) on a sphere.
 * - Great-circle distance and a squared-distance kernel-gradient term.
 * - A simple "inverse transform" that applies a local tangent-plane offset (x,y in meters)
 *   to a reference lon/lat and returns the new lon/lat.
 *
 * ## Conventions
 * - All angles (lon/lat) are in **radians**.
 * - Great-circle distance returned by GCD / GCD_deriv is the **central angle** (radians),
 *   i.e. arc length is `R * angle`.
 *
 * ## Numerical robustness
 * - clamp() is used to keep floating-point roundoff from violating the domain constraints
 *   of trig inverse operations.
 * - GCD_deriv() does **not** throw. If the derivative is mathematically undefined or numerically
 *   unstable because the two points are coincident (or nearly coincident), it returns:
 *     { d, 0, 0 }
 *   and emits a very cheap warning to std::cerr (at most once).
 *
 * ## Optional warning control
 * - By default, one warning per process is printed the first time a degenerate derivative is detected.
 * - Define `PARTICLE_TRACKER_GCD_DERIV_WARN=0` to disable all warnings.
 * - Define `PARTICLE_TRACKER_GCD_DERIV_WARN_ALWAYS=1` to warn every time (not recommended on HPC).
 */

#include <algorithm>    // std::min, std::max
#include <array>
#include <cmath>        // std::sin, std::cos, std::sqrt, std::atan2, std::asin
#include <iostream>     // std::cerr
#include <atomic>       // std::atomic_flag

namespace particle_tracker
{

//------------------------------------------------------------------------------
// Compile-time configuration
//------------------------------------------------------------------------------

#ifndef PARTICLE_TRACKER_GCD_DERIV_WARN
#define PARTICLE_TRACKER_GCD_DERIV_WARN 1
#endif

#ifndef PARTICLE_TRACKER_GCD_DERIV_WARN_ALWAYS
#define PARTICLE_TRACKER_GCD_DERIV_WARN_ALWAYS 0
#endif

//------------------------------------------------------------------------------
// Utility helpers
//------------------------------------------------------------------------------

/**
 * @brief Clamp a value to the closed interval [lo, hi].
 *
 * @tparam DataType numeric type
 * @param val input value
 * @param lo lower bound
 * @param hi upper bound
 * @return clamped value
 */
template<typename DataType>
inline DataType clamp(DataType val, DataType lo, DataType hi)
{
    return std::max(lo, std::min(val, hi));
}

//------------------------------------------------------------------------------
// Great-circle distance
//------------------------------------------------------------------------------

/**
 * @brief Great-circle distance (central angle) between two points on a sphere.
 *
 * @tparam DataType floating-point type (float/double/long double)
 * @param lon1 longitude of point 1 (radians)
 * @param lat1 latitude  of point 1 (radians)
 * @param lon2 longitude of point 2 (radians)
 * @param lat2 latitude  of point 2 (radians)
 * @return central angle distance in radians (in [0, pi])
 *
 * Implementation:
 * - Computes cos(d) using spherical law of cosines.
 * - Uses atan2(sin(d), cos(d)) for numerical robustness.
 */
template<typename DataType>
inline DataType GCD(const DataType& lon1, const DataType& lat1,
                    const DataType& lon2, const DataType& lat2)
{
    DataType cn =
        std::sin(lat1) * std::sin(lat2) +
        std::cos(lat1) * std::cos(lat2) * std::cos(lon1 - lon2);

    cn = clamp(cn, static_cast<DataType>(-1), static_cast<DataType>(1));

    DataType sn = std::sqrt(static_cast<DataType>(1) - cn * cn);

    return std::atan2(sn, cn);
}

/**
 * @brief Great-circle distance and the local gradient of `-d^2/2`.
 *
 * @tparam DataType floating-point type
 * @param lon1 longitude of point 1 (radians)
 * @param lat1 latitude  of point 1 (radians)
 * @param lon2 longitude of point 2 (radians)
 * @param lat2 latitude  of point 2 (radians)
 * @return array { d, kernel_grad_east, kernel_grad_north }
 *
 * Where:
 * - d                 : central-angle distance (radians)
 * - kernel_grad_east  : `-d * grad(d)` in the local east direction
 * - kernel_grad_north : `-d * grad(d)` in the local north direction
 *
 * The east/north components use a local orthonormal tangent basis. In
 * particular, the east component includes the `1/cos(lat1)` metric factor;
 * it is not the coordinate derivative with respect to longitude.
 *
 * ### Degenerate case handling (d ~ 0)
 * The analytic expressions contain division by sin(d). When the two points coincide (or are
 * extremely close), the direction of grad(d) is undefined. The product
 * `-d*grad(d)` nevertheless tends to zero.
 *
 * Policy:
 * - Return {d, 0, 0} (safe "no directional information" fallback).
 * - Optionally print a warning to std::cerr:
 *     - by default: at most once per process (very cheap)
 *     - disable:   PARTICLE_TRACKER_GCD_DERIV_WARN=0
 *     - always:    PARTICLE_TRACKER_GCD_DERIV_WARN_ALWAYS=1 (not recommended for HPC)
 *
 * This quantity is the geometric factor in the derivative of a Gaussian
 * weight: `grad(exp(-d^2/shape)) = (2/shape) * weight * grad(-d^2/2)`.
 */
template<typename DataType>
inline std::array<DataType, 3> GCD_deriv(const DataType& lon1, const DataType& lat1,
                                        const DataType& lon2, const DataType& lat2)
{
    // cos(d)
    DataType cn =
        std::sin(lat1) * std::sin(lat2) +
        std::cos(lat1) * std::cos(lat2) * std::cos(lon1 - lon2);

    cn = clamp(cn, static_cast<DataType>(-1), static_cast<DataType>(1));

    // sin(d)
    DataType sn = std::sqrt(static_cast<DataType>(1) - cn * cn);

    // d
    DataType d = std::atan2(sn, cn);

    // Guard: derivatives contain division by sin(d)
    // Choose eps relative to DataType.
    // - For double, 1e-14 is conservative; for float, it still works.
    const DataType eps = static_cast<DataType>(1e-14);

    if (sn <= eps)
    {
#if PARTICLE_TRACKER_GCD_DERIV_WARN
#if PARTICLE_TRACKER_GCD_DERIV_WARN_ALWAYS
        std::cerr
            << "particle_tracker::GCD_deriv: sin(d) ~ 0 (coincident or nearly coincident points). "
            << "Returning {d,0,0} to avoid NaN/Inf.\n";
#else
        // Print at most once per process.
        static std::atomic_flag warned = ATOMIC_FLAG_INIT;
        if (!warned.test_and_set(std::memory_order_relaxed))
        {
            std::cerr
                << "particle_tracker::GCD_deriv: sin(d) ~ 0 (coincident or nearly coincident points). "
                << "Returning {d,0,0} to avoid NaN/Inf. (This warning is printed once.)\n";
        }
#endif
#endif
        return {d, static_cast<DataType>(0), static_cast<DataType>(0)};
    }

    // East component of grad(-d^2/2) in the local orthonormal tangent basis.
    DataType kernel_grad_east =
        d * std::cos(lat2) * std::sin(lon2 - lon1) / sn;

    // North component of grad(-d^2/2).
    DataType kernel_grad_north =
        d * (std::cos(lat1) * std::sin(lat2) -
             std::sin(lat1) * std::cos(lat2) * std::cos(lon2 - lon1)) / sn;

    return {d, kernel_grad_east, kernel_grad_north};
}

//------------------------------------------------------------------------------
// Local tangent-plane offset -> lon/lat update
//------------------------------------------------------------------------------

/**
 * @brief Apply a local tangent-plane offset (x,y in meters) to a reference lon/lat (radians).
 *
 * This implements the spherical **azimuthal equidistant (AEQD)** inverse transform
 * centered at the reference point (lon,lat). The offset (x,y) is interpreted in the
 * local tangent plane:
 *   - x: east (meters)
 *   - y: north (meters)
 *
 * @tparam DataType floating-point type
 * @param x eastward offset  (meters)
 * @param y northward offset (meters)
 * @param lon [in,out] reference longitude (radians), updated in-place
 * @param lat [in,out] reference latitude  (radians), updated in-place
 *
 * The returned point lies at great-circle distance |(x,y)| from the center, with
 * bearing given by (x,y).
 *
 * @note Spherical Earth model with radius EARTH_RADIUS (not WGS84 ellipsoid).
 * @note Accuracy degrades for very large offsets (hundreds–thousands of km).
 */
template<typename DataType>
inline void inverseTransform(const DataType& x, const DataType& y, DataType& lon, DataType& lat)
{
    constexpr DataType EARTH_RADIUS = static_cast<DataType>(6378000.0);

    DataType nrm = std::sqrt(x * x + y * y);

    const DataType eps = static_cast<DataType>(1e-10);
    if (nrm < eps)
        return;

    DataType gcd = nrm / EARTH_RADIUS;

    DataType cn = std::cos(gcd);
    DataType sn = std::sin(gcd);

    DataType cl = std::cos(lat);
    DataType sl = std::sin(lat);

    DataType arg = sl * cn + (y / nrm) * cl * sn;
    arg = clamp(arg, static_cast<DataType>(-1), static_cast<DataType>(1));

    DataType lat_prime = std::asin(arg);

    DataType lon_prime = lon + std::atan2(x * sn, nrm * cl * cn - y * sl * sn);

    lat = lat_prime;
    lon = lon_prime;
}

} // namespace particle_tracker
