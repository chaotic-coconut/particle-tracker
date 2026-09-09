#include "test_helpers.hpp"

#include <particle_tracker/geography_utils.hpp>

#include <array>
#include <cmath>
#include <numbers>

namespace {

constexpr double earth_radius_m = particle_tracker::earth_mean_radius_m;

std::array<double, 2> forward_aeqd(double center_lon, double center_lat,
                                   double lon, double lat) {
  const double delta_lon = lon - center_lon;
  const double sin_center_lat = std::sin(center_lat);
  const double cos_center_lat = std::cos(center_lat);
  const double sin_lat = std::sin(lat);
  const double cos_lat = std::cos(lat);

  const double cos_c = particle_tracker::clamp(
      sin_center_lat * sin_lat +
          cos_center_lat * cos_lat * std::cos(delta_lon),
      -1.0, 1.0);
  const double c = std::acos(cos_c);

  if (c < 1e-14)
    return {0.0, 0.0};

  const double scale = earth_radius_m * c / std::sin(c);
  return {scale * cos_lat * std::sin(delta_lon),
          scale * (cos_center_lat * sin_lat -
                   sin_center_lat * cos_lat * std::cos(delta_lon))};
}

void check_aeqd_round_trip(double center_lon, double center_lat, double x,
                           double y) {
  double lon = center_lon;
  double lat = center_lat;
  particle_tracker::inverseTransform(x, y, lon, lat);

  const auto recovered = forward_aeqd(center_lon, center_lat, lon, lat);
  require_near(recovered[0], x, 1e-5, "AEQD eastward round trip");
  require_near(recovered[1], y, 1e-5, "AEQD northward round trip");

  const double central_angle =
      particle_tracker::GCD(center_lon, center_lat, lon, lat);
  require_near(central_angle * earth_radius_m, std::hypot(x, y), 1e-5,
               "AEQD radial distance");
}

void check_kernel_gradient() {
  constexpr double lon1 = 0.4;
  constexpr double lat1 = 0.6;
  constexpr double lon2 = 0.43;
  constexpr double lat2 = 0.58;
  constexpr double epsilon = 1e-6;

  const auto gradient =
      particle_tracker::GCD_deriv(lon1, lat1, lon2, lat2);
  const double derivative_lon =
      (particle_tracker::GCD(lon1 + epsilon, lat1, lon2, lat2) -
       particle_tracker::GCD(lon1 - epsilon, lat1, lon2, lat2)) /
      (2.0 * epsilon);
  const double derivative_lat =
      (particle_tracker::GCD(lon1, lat1 + epsilon, lon2, lat2) -
       particle_tracker::GCD(lon1, lat1 - epsilon, lon2, lat2)) /
      (2.0 * epsilon);

  const double expected_east =
      -gradient[0] * derivative_lon / std::cos(lat1);
  const double expected_north = -gradient[0] * derivative_lat;
  require_near(gradient[1], expected_east, 2e-10,
               "GCD kernel-gradient east convention");
  require_near(gradient[2], expected_north, 2e-10,
               "GCD kernel-gradient north convention");
}

} // namespace

int main() {
  constexpr double pi = std::numbers::pi_v<double>;

  require_near(particle_tracker::GCD(0.0, 0.0, 0.0, 0.0), 0.0, 1e-15,
               "Coincident great-circle distance");
  require_near(particle_tracker::GCD(0.0, 0.0, pi / 2.0, 0.0), pi / 2.0,
               1e-15, "Quarter-circumference great-circle distance");
  require_near(particle_tracker::GCD(0.0, 0.0, pi, 0.0), pi, 1e-15,
               "Antipodal great-circle distance");

  check_aeqd_round_trip(0.2, 0.5, 1'500.0, -2'600.0);
  check_aeqd_round_trip(pi - 1e-4, 0.7, 4'000.0, 1'500.0);
  check_aeqd_round_trip(-2.4, 1.1, -1'200.0, 800.0);
  check_kernel_gradient();

  double lon = -1.2;
  double lat = 0.4;
  particle_tracker::inverseTransform(0.0, 0.0, lon, lat);
  require_near(lon, -1.2, 0.0, "Zero offset longitude");
  require_near(lat, 0.4, 0.0, "Zero offset latitude");
}
