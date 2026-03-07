cmake_minimum_required(VERSION 3.20)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

project(particle_tracker LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(PK_USE_ZLIB "Use zlib backend" ON)
option(PK_USE_ZSTD "Use zstd backend" OFF)

if(PK_USE_ZLIB AND PK_USE_ZSTD)
  message(FATAL_ERROR "Choose only one of PK_USE_ZLIB or PK_USE_ZSTD")
endif()

add_library(particle_tracker INTERFACE)

target_include_directories(particle_tracker
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/third_party/nanoflann>
)

target_compile_features(particle_tracker INTERFACE cxx_std_20)

if(PK_USE_ZLIB)
  target_compile_definitions(particle_tracker INTERFACE PK_IO_USE_ZLIB)
elseif(PK_USE_ZSTD)
  target_compile_definitions(particle_tracker INTERFACE PK_IO_USE_ZSTD)
endif()

# ---- Dependencies ----
find_package(ZLIB REQUIRED)

find_package(TBB QUIET)
if(NOT TBB_FOUND)
  message(STATUS "TBB package not found; will try raw -ltbb fallback")
endif()

find_package(Boost REQUIRED CONFIG)
target_link_libraries(particle_tracker INTERFACE Boost::boost)

find_package(PkgConfig QUIET)
find_package(netCDF QUIET)

set(PK_NETCDF_C_LIB "")
set(PK_NETCDF_CXX_LIB "")
set(PK_NETCDF_INCLUDE_DIRS "")

# ---- netCDF C / C++ detection ----
# Try imported CMake targets first
if(netCDF_FOUND)
  if(TARGET netCDF::netcdf)
    set(PK_NETCDF_C_LIB netCDF::netcdf)
  endif()

  if(TARGET netCDF::netcdf-cxx4)
    set(PK_NETCDF_CXX_LIB netCDF::netcdf-cxx4)
  elseif(TARGET netCDF::netcdf_c++4)
    set(PK_NETCDF_CXX_LIB netCDF::netcdf_c++4)
  endif()
endif()

# Fallback to pkg-config on systems like Arch
if(PkgConfig_FOUND)
  if(NOT PK_NETCDF_C_LIB)
    pkg_check_modules(NETCDF_C QUIET netcdf)
    if(NETCDF_C_FOUND)
      set(PK_NETCDF_C_LIB ${NETCDF_C_LINK_LIBRARIES})
      list(APPEND PK_NETCDF_INCLUDE_DIRS ${NETCDF_C_INCLUDE_DIRS})
    endif()
  endif()

  if(NOT PK_NETCDF_CXX_LIB)
    pkg_check_modules(NETCDF_CXX QUIET netcdf-cxx4)
    if(NETCDF_CXX_FOUND)
      set(PK_NETCDF_CXX_LIB ${NETCDF_CXX_LINK_LIBRARIES})
      list(APPEND PK_NETCDF_INCLUDE_DIRS ${NETCDF_CXX_INCLUDE_DIRS})
    endif()
  endif()
endif()

# Final fallback: raw library search
if(NOT PK_NETCDF_C_LIB)
  find_library(PK_NETCDF_C_LIB NAMES netcdf REQUIRED)
endif()

if(NOT PK_NETCDF_CXX_LIB)
  find_library(PK_NETCDF_CXX_LIB NAMES netcdf_c++4 netcdf-cxx4 REQUIRED)
endif()

list(REMOVE_DUPLICATES PK_NETCDF_INCLUDE_DIRS)

message(STATUS "netCDF C library   : ${PK_NETCDF_C_LIB}")
message(STATUS "netCDF C++ library : ${PK_NETCDF_CXX_LIB}")
if(PK_NETCDF_INCLUDE_DIRS)
  message(STATUS "netCDF include dirs: ${PK_NETCDF_INCLUDE_DIRS}")
endif()

# ---- Helper to link common deps ----
function(link_common_deps tgt)
  target_link_libraries(${tgt}
    PRIVATE
      particle_tracker
      ZLIB::ZLIB
      Boost::boost
      ${PK_NETCDF_C_LIB}
      ${PK_NETCDF_CXX_LIB}
  )

  if(PK_NETCDF_INCLUDE_DIRS)
    target_include_directories(${tgt} PRIVATE ${PK_NETCDF_INCLUDE_DIRS})
  endif()

  if(TBB_FOUND)
    target_link_libraries(${tgt} PRIVATE TBB::tbb)
  else()
    target_link_libraries(${tgt} PRIVATE tbb)
  endif()
endfunction()

# ---- Examples ----
add_executable(hpc_particles examples/hpc_particles.cpp)
link_common_deps(hpc_particles)

add_executable(hpc_trajectories examples/hpc_trajectories.cpp)
link_common_deps(hpc_trajectories)
