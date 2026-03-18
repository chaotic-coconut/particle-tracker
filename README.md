# particle-transport-framework

C++ framework for large-scale particle transport simulation and ocean flow data processing.

The framework is designed for long-running simulations based on gridded environmental velocity fields (e.g., HYCOM NetCDF data), with emphasis on performance, modular architecture, and reproducibility.

It provides building blocks for constructing simulation pipelines that operate on structured flow fields and track particle trajectories over time in complex spatial domains.

## Features

- Modular architecture for simulation workflows and data processing pipelines  
- Efficient handling and reorganization of gridded flow data (NetCDF, HYCOM)  
- Spatial indexing using kd-trees (nanoflann) for fast nearest-neighbor queries  
- Local geometric transformations and coordinate handling  
- Interpolation and reconstruction of physical fields (Gaussian weighting, spline-based methods)  
- Parallel processing using Intel TBB for efficient handling of large datasets and simulation workloads
- Support for long-duration simulations and repeated runs

## Tech Stack

- C++
- CMake
- Boost
- NetCDF
- Intel TBB
- Linux

## Status

The repository exposes the core structure and interfaces of the framework.  
Implementation is being actively cleaned, modularized, and documented.


