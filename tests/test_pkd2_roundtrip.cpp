#include "test_helpers.hpp"

#include <particle_tracker/fixed_point/fixed_point_pkd2.hpp>
#include <particle_tracker/seed_loader_v2.hpp>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace {

struct TemporaryFile {
  std::filesystem::path path;

  ~TemporaryFile() {
    std::error_code error;
    std::filesystem::remove(path, error);
  }
};

double unwrap_after(double previous, double value) {
  while (value <= -fp_pi)
    value += 2.0 * fp_pi;
  while (value > fp_pi)
    value -= 2.0 * fp_pi;

  const double delta = value - previous;
  if (delta > fp_pi)
    value -= 2.0 * fp_pi;
  else if (delta < -fp_pi)
    value += 2.0 * fp_pi;
  return value;
}

} // namespace

int main() {
  const auto unique_suffix =
      std::chrono::steady_clock::now().time_since_epoch().count();
  TemporaryFile file{std::filesystem::temp_directory_path() /
                     ("particle_tracker_pkd2_" +
                      std::to_string(unique_suffix) + ".pkd2")};

  constexpr std::int64_t epoch_s = 946'684'800;
  constexpr std::int32_t time_step_s = 10'800;
  constexpr std::uint64_t seed_id = 42;
  const std::vector<double> lon{3.10, 3.13, -3.13, -3.10};
  const std::vector<double> lat{0.20, 0.21, 0.22, 0.23};

  {
    pkd2::Writer writer(file.path.string(), epoch_s, time_step_s);
    writer.add_traj_rad(seed_id, lon.data(), lat.data(),
                        static_cast<std::uint32_t>(lon.size()), 7);
    writer.close();
  }

  pkd2::Reader reader(file.path.string());
  require(reader.header().version == 3, "PKD2 version remains 3");
  require(reader.header().n_traj == 1, "PKD2 trajectory count");
  require(reader.has(seed_id), "PKD2 seed lookup");

  std::vector<double> decoded_lon;
  std::vector<double> decoded_lat;
  std::vector<std::int64_t> decoded_time;
  reader.read_traj(seed_id, decoded_lon, decoded_lat, decoded_time);

  require(decoded_lon.size() == lon.size(), "PKD2 longitude sample count");
  require(decoded_lat.size() == lat.size(), "PKD2 latitude sample count");
  require(decoded_time.size() == lon.size(), "PKD2 time sample count");

  constexpr double coordinate_tolerance = 2.1e-6;
  double expected_lon = wrap_lon(lon.front());
  require_near(decoded_lon.front(), expected_lon, coordinate_tolerance,
               "PKD2 initial longitude");
  require_near(decoded_lat.front(), lat.front(), coordinate_tolerance,
               "PKD2 initial latitude");

  for (std::size_t i = 1; i < lon.size(); ++i) {
    expected_lon = unwrap_after(expected_lon, lon[i]);
    require_near(decoded_lon[i], expected_lon, coordinate_tolerance,
                 "PKD2 unwrapped longitude");
    require_near(decoded_lat[i], lat[i], coordinate_tolerance,
                 "PKD2 latitude");
  }

  for (std::size_t i = 0; i < decoded_time.size(); ++i) {
    const std::int64_t expected_time =
        epoch_s + static_cast<std::int64_t>(7 + i) * time_step_s;
    require(decoded_time[i] == expected_time, "PKD2 sample time");
  }

  const auto last_seed = load_seeds_auto(file.path.string(), SeedPoint::Last);
  require(last_seed.size() == 1, "PKD2 seed loader count");
  require(last_seed.front().seed_id == seed_id, "PKD2 seed loader id");
  require_near(last_seed.front().lon_rad, wrap_lon(lon.back()),
               coordinate_tolerance, "PKD2 last seed longitude");
  require_near(last_seed.front().lat_rad, lat.back(), coordinate_tolerance,
               "PKD2 last seed latitude");
  require(last_seed.front().t_start_s ==
              epoch_s + static_cast<std::int64_t>(10) * time_step_s,
          "PKD2 last seed time");
  require(last_seed.front().has_stop, "PKD2 last seed stop metadata");
  require_near(last_seed.front().lon_stop_rad, wrap_lon(lon.front()),
               coordinate_tolerance, "PKD2 stop longitude");

  bool missing_seed_threw = false;
  try {
    reader.read_traj(999, decoded_lon, decoded_lat, decoded_time);
  } catch (const std::runtime_error &) {
    missing_seed_threw = true;
  }
  require(missing_seed_threw, "PKD2 missing seed must throw");

}
