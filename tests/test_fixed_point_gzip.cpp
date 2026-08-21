#include "test_helpers.hpp"

#include <particle_tracker/fixed_point/fixed_point_io_codec.hpp>

#include <chrono>
#include <cstring>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

namespace {

struct TemporaryFile {
  std::filesystem::path path;

  ~TemporaryFile() {
    std::error_code error;
    std::filesystem::remove(path, error);
  }
};

bool records_equal(const PackedParticle &lhs, const PackedParticle &rhs) {
  return std::memcmp(&lhs, &rhs, sizeof(PackedParticle)) == 0;
}

} // namespace

int main() {
  const auto unique_suffix =
      std::chrono::steady_clock::now().time_since_epoch().count();
  TemporaryFile file{std::filesystem::temp_directory_path() /
                     ("particle_tracker_pkd1_" +
                      std::to_string(unique_suffix) + ".pkd.gz")};

  const std::vector<PackedParticle> initial{
      encode_record(0.1, 0.2, 1.0, 0.3, 0.4, 2.0, 1),
      encode_record(-2.8, -0.5, 3.0, 2.9, 0.6, 4.0, 2)};
  const std::vector<PackedParticle> appended{
      encode_record(1.1, -0.2, 5.0, 1.2, -0.1, 6.0, 3)};

  write_records_gzip(file.path.string(), initial, 12'345, 6, false);
  write_records_gzip(file.path.string(), appended, 12'345, 6, true);
  const auto decoded = read_records_gzip(file.path.string());

  require(decoded.size() == initial.size() + appended.size(),
          "Concatenated gzip record count");
  require(records_equal(decoded[0], initial[0]),
          "First gzip record round trip");
  require(records_equal(decoded[1], initial[1]),
          "Second gzip record round trip");
  require(records_equal(decoded[2], appended[0]),
          "Appended gzip record round trip");
}
