# Particle Transport Framework

[![ci](https://github.com/chaotic-coconut/particle-transport-framework/actions/workflows/ci.yml/badge.svg)](https://github.com/chaotic-coconut/particle-transport-framework/actions/workflows/ci.yml)

This repository contains the C++20 particle-transport code I developed for
simulations in time-dependent ocean velocity fields. It reads gridded NetCDF
data, interpolates the velocity field in space and time, advances particles on
a spherical Earth, and writes compact trajectory records.

The public drivers reflect the structure of the production calculations. They
are not standalone examples: running them requires the original external
NetCDF datasets and the corresponding variable names and file patterns.

## Main components

- NetCDF C++4 input for HYCOM-style gridded fields
- cubic spline interpolation in time
- 2D longitude/latitude kd-tree searches with nanoflann
- Gaussian spatial weighting based on great-circle distance
- explicit Euler particle updates in a local east/north tangent plane
- spherical azimuthal-equidistant conversion back to longitude/latitude
- shared-memory parallel particle propagation with Intel oneAPI TBB
- two-month rolling field buffers
- fixed-point trajectory records compressed with zlib

## Production workload represented by the code

`examples/hpc_particles.cpp` is configured for 2,800,000 released particles per
day and a maximum lifetime of 730 days. That is about 1.02 billion release
events per 365-day model year. Coastal seed positions are generated once and
reused for each release day. Completed particles are written and removed from
memory while two adjacent months of velocity fields remain loaded.

These numbers describe the workload for which the driver was written, not a
benchmark. The repository does not contain the input data, production logs,
batch scripts, or downstream ensemble-analysis workflow.

## Build

The code requires CMake 3.20 or newer, a C++20 compiler, Boost, NetCDF C and
C++4, Intel oneAPI TBB, and zlib. On Ubuntu:

```sh
sudo apt-get install cmake ninja-build g++ zlib1g-dev libtbb-dev \
  libnetcdf-dev libnetcdf-c++4-dev libboost-dev
```

Configure, build, and run the tests:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The build produces two executables:

- `hpc_particles` releases daily particle ensembles and stores compressed
  start/end records for landed particles.
- `hpc_trajectories` reconstructs position and velocity time series from seed
  records.

Running either executable without arguments prints its expected command-line
arguments.

## Tests

The CTest suite covers geographic transformations, great-circle calculations,
the `GCD_deriv` kernel-gradient convention and scaling, calendar edge cases,
fixed-point packing, gzip append/read behaviour, PKD2 position round trips, and
small deterministic interpolation cases.

These are regression and numerical sanity checks. They do not replace
validation with the original ocean fields or a production-scale run.

## License

BSD-3-Clause. See [LICENSE](LICENSE).
