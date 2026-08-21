#include "test_helpers.hpp"

#include <particle_tracker/fixed_point/fixed_point_core.hpp>

#include <cstdint>

int main() {
  constexpr double lon0 = 3.5;
  constexpr double lat0 = 0.45;
  constexpr double t0 = 1234.4;
  constexpr double lon1 = -3.4;
  constexpr double lat1 = -0.72;
  constexpr double t1 = 9876.6;
  constexpr std::uint8_t input_flags = 0xA5;

  const PackedParticle packed =
      encode_record(lon0, lat0, t0, lon1, lat1, t1, input_flags);

  double decoded_lon0 = 0.0;
  double decoded_lat0 = 0.0;
  double decoded_t0 = 0.0;
  double decoded_lon1 = 0.0;
  double decoded_lat1 = 0.0;
  double decoded_t1 = 0.0;
  std::uint8_t decoded_flags = 0;

  decode_record(packed, decoded_lon0, decoded_lat0, decoded_t0, decoded_lon1,
                decoded_lat1, decoded_t1, decoded_flags);

  constexpr double coordinate_tolerance = 0.5 / fp_ticks_per_rad + 1e-15;
  require_near(decoded_lon0, wrap_lon(lon0), coordinate_tolerance,
               "Packed lon0 round trip");
  require_near(decoded_lat0, lat0, coordinate_tolerance,
               "Packed lat0 round trip");
  require_near(decoded_lon1, wrap_lon(lon1), coordinate_tolerance,
               "Packed lon1 round trip");
  require_near(decoded_lat1, lat1, coordinate_tolerance,
               "Packed lat1 round trip");
  require_near(decoded_t0, 1234.0, 0.0, "Packed t0 rounding");
  require_near(decoded_t1, 9877.0, 0.0, "Packed t1 rounding");
  require(decoded_flags == input_flags, "Packed flags round trip");
}
