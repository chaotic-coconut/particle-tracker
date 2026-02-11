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
)

if(PK_USE_ZLIB)
  target_compile_definitions(particle_tracker INTERFACE PK_IO_USE_ZLIB)
elseif(PK_USE_ZSTD)
  target_compile_definitions(particle_tracker INTERFACE PK_IO_USE_ZSTD)
endif()

# ---- Dependencies ----
find_package(ZLIB REQUIRED)

find_package(TBB QUIET)
if(NOT TBB_FOUND)
  message(STATUS "TBB package not found; linking with -ltbb fallback")
endif()

find_package(Boost REQUIRED)
target_link_libraries(particle_tracker INTERFACE Boost::boost)

find_package(netCDF QUIET)

# ---- Helper to link common deps ----
function(link_common_deps tgt)
  target_link_libraries(${tgt} PRIVATE particle_tracker ZLIB::ZLIB)

  if(TBB_FOUND)
    target_link_libraries(${tgt} PRIVATE TBB::tbb)
  else()
    target_link_libraries(${tgt} PRIVATE tbb)
  endif()

  if(netCDF_FOUND)
    # Different distros expose different target names; try the common ones.
    if(TARGET netCDF::netcdf)
      target_link_libraries(${tgt} PRIVATE netCDF::netcdf)
    endif()
    if(TARGET netCDF::netcdf-cxx4)
      target_link_libraries(${tgt} PRIVATE netCDF::netcdf-cxx4)
    elseif(TARGET netCDF::netcdf_c++4)
      target_link_libraries(${tgt} PRIVATE netCDF::netcdf_c++4)
    endif()
  else()
    # Fallback: may or may not exist on the machine.
    target_link_libraries(${tgt} PRIVATE netcdf netcdf_c++4)
  endif()
endfunction()

# ---- Examples ----
add_executable(hpc_particles examples/hpc_particles.cpp)
link_common_deps(hpc_particles)

add_executable(hpc_trajectories examples/hpc_trajectories.cpp)
link_common_deps(hpc_trajectories)
