cmake_minimum_required(VERSION 3.20)
project(particle_tracker LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(PK_USE_ZLIB "Use zlib backend" ON)
option(PK_USE_ZSTD "Use zstd backend" OFF)

if(PK_USE_ZLIB AND PK_USE_ZSTD)
  message(FATAL_ERROR "Choose only one of PK_USE_ZLIB or PK_USE_ZSTD")
endif()

add_library(particle_tracker INTERFACE)
target_include_directories(particle_tracker INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include/particle_tracker
)

if(PK_USE_ZLIB)
  target_compile_definitions(particle_tracker INTERFACE PK_IO_USE_ZLIB)
elseif(PK_USE_ZSTD)
  target_compile_definitions(particle_tracker INTERFACE PK_IO_USE_ZSTD)
endif()

# ---- Dependencies (system-installed) ----
find_package(ZLIB REQUIRED)

# TBB: prefer config package if available, otherwise fall back to tbb library
find_package(TBB QUIET)
if(NOT TBB_FOUND)
  message(STATUS "TBB package not found; linking with -ltbb fallback")
endif()

# NetCDF: try config; many HPC systems provide netcdf libs without config
find_package(netCDF QUIET)

# ---- Examples ----
add_executable(hpc_particles examples/hpc_particles.cpp)
target_link_libraries(hpc_particles PRIVATE particle_tracker ZLIB::ZLIB)
if(TBB_FOUND)
  target_link_libraries(hpc_particles PRIVATE TBB::tbb)
else()
  target_link_libraries(hpc_particles PRIVATE tbb)
endif()
# netcdf fallback (keep minimal; you can refine later)
target_link_libraries(hpc_particles PRIVATE netcdf netcdf_c++4)

add_executable(hpc_trajectories examples/hpc_trajectories.cpp)
target_link_libraries(hpc_trajectories PRIVATE particle_tracker ZLIB::ZLIB)
if(TBB_FOUND)
  target_link_libraries(hpc_trajectories PRIVATE TBB::tbb)
else()
  target_link_libraries(hpc_trajectories PRIVATE tbb)
endif()
target_link_libraries(hpc_trajectories PRIVATE netcdf netcdf_c++4)
