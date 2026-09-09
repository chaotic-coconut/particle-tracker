#include "test_helpers.hpp"

#include <particle_tracker/interpolation_utils.hpp>

#include <map>
#include <string>
#include <vector>

int main() {
  using particle_tracker::KDTree;
  using particle_tracker::Neighbor;
  using particle_tracker::NeighborWithNeighbors;
  using particle_tracker::NeighborWithNeighborsData;
  using particle_tracker::NeighborsData;
  using particle_tracker::PointCloud;
  using particle_tracker::TimeSpline;

  constexpr double h = 1e-3;
  const std::vector<double> lon{h, -h, 0.0, 0.0};
  const std::vector<double> lat{0.0, 0.0, h, -h};

  PointCloud<double> cloud;
  for (std::size_t i = 0; i < lon.size(); ++i) {
    cloud.pts.push_back(lon[i]);
    cloud.pts.push_back(lat[i]);
  }

  const nanoflann::KDTreeSingleIndexAdaptorParams params(2);
  KDTree<double> index(2, cloud, params);
  index.buildIndex();

  std::map<std::string, std::vector<TimeSpline<double>>> fields;
  auto &constant = fields["constant"];
  auto &linear = fields["linear"];
  constant.resize(lon.size());
  linear.resize(lon.size());

  for (std::size_t i = 0; i < lon.size(); ++i) {
    constant[i].values.assign(8, 7.25);
    linear[i].values.assign(8, 4.0 + 2.0 * lon[i] - 3.0 * lat[i]);
    constant[i].createSpline(0.0, 1.0);
    linear[i].createSpline(0.0, 1.0);
  }

  NeighborsData<double> interpolator(index, lon, lat, 1.1 * h * h, h * h);
  std::vector<Neighbor<double>> neighbors;
  const auto values = interpolator.calculateVariables(
      3.5, neighbors, fields, {"constant", "linear"}, 0.0, 0.0);

  require(neighbors.size() == 4,
          "Symmetric interpolation test must use all four neighbors");
  require_near(values.at("constant"), 7.25, 1e-12,
               "Normalized interpolation reproduces constants");
  require_near(values.at("linear"), 4.0, 1e-12,
               "Symmetric interpolation reproduces a linear field at center");

  const std::vector<double> derivative_lon{0.0, h, -h, 0.0, 0.0};
  const std::vector<double> derivative_lat{0.0, 0.0, 0.0, h, -h};
  PointCloud<double> derivative_cloud;
  for (std::size_t i = 0; i < derivative_lon.size(); ++i) {
    derivative_cloud.pts.push_back(derivative_lon[i]);
    derivative_cloud.pts.push_back(derivative_lat[i]);
  }
  KDTree<double> derivative_index(2, derivative_cloud, params);
  derivative_index.buildIndex();

  std::map<std::string, std::vector<TimeSpline<double>>> derivative_fields;
  auto &derivative_linear = derivative_fields["linear"];
  auto &derivative_constant = derivative_fields["constant"];
  derivative_linear.resize(derivative_lon.size());
  derivative_constant.resize(derivative_lon.size());
  for (std::size_t i = 0; i < derivative_lon.size(); ++i) {
    derivative_linear[i].values.assign(
        8, 4.0 + 2.0 * derivative_lon[i] - 3.0 * derivative_lat[i]);
    derivative_linear[i].createSpline(0.0, 1.0);
    derivative_constant[i].values.assign(8, 7.25);
    derivative_constant[i].createSpline(0.0, 1.0);
  }

  NeighborWithNeighbors<double> primary;
  primary.ind = 0;
  primary.weight = 1.0;
  for (std::size_t i = 1; i < derivative_lon.size(); ++i) {
    particle_tracker::NeighborPrime<double> secondary;
    secondary.ind_prime = i;
    secondary.dist_prime = particle_tracker::GCD_deriv(
        0.0, 0.0, derivative_lon[i], derivative_lat[i]);
    secondary.weight_prime =
        std::exp(-secondary.dist_prime[0] * secondary.dist_prime[0] /
                 (h * h));
    primary.neighbors_prime.push_back(secondary);
  }

  NeighborWithNeighborsData<double> derivative_interpolator(
      derivative_index, derivative_lon, derivative_lat, 4.1 * h * h, h * h);
  const auto derivatives = derivative_interpolator.interpolateDerivatives(
      3.5, {primary}, derivative_fields, {"constant", "linear"}, h * h);

  require_near(derivatives.at("linear").spatial_derivative_x, 2.0, 1e-5,
               "Analytic linear-field east derivative");
  require_near(derivatives.at("linear").spatial_derivative_y, -3.0, 1e-5,
               "Analytic linear-field north derivative");

  // Zeroth-order consistency: a constant field has no gradient. This is exact
  // rather than approximate, because every difference v_j - v_i is exactly zero
  // regardless of the stencil, so it is asserted at zero tolerance.
  require_near(derivatives.at("constant").spatial_derivative_x, 0.0, 0.0,
               "Constant field has zero east derivative");
  require_near(derivatives.at("constant").spatial_derivative_y, 0.0, 0.0,
               "Constant field has zero north derivative");
  require_near(derivatives.at("constant").value, 7.25, 1e-12,
               "Constant field value is reproduced by the derivative path");
}
