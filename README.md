# particle-transport-framework

This repository shows the interfaces and overall structure of a C++ simulation framework
I use for long-running numerical experiments on particle dispersal in the ocean.
At the moment, the public version mainly contains header files that define the main components,
data flow, and validation logic in the code.

The implementation itself is currently being cleaned up and documented.

The headers provide C++ wrappers for reading HYCOM NetCDF data and reorganizing gridded flow fields.
They include utilities for geometric transformations, in particular the construction of local
tangent planes using an equidistant azimuthal projection.
The code supports kd-tree–based spatial searches (via nanoflann) at fixed depth levels,
nearest-neighbour queries, and Gaussian-weighted reconstruction of field values.
It also contains experimental routines for estimating first spatial derivatives and wrappers
for spline-based time interpolation using Boost.

