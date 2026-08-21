cmake_minimum_required(VERSION 3.20)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

project(particle_transport_framework LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel)
endif()

option(PARTICLE_TRACKER_BUILD_EXAMPLES "Build the simulation executables" ON)

add_library(particle_tracker INTERFACE)
add_library(particle_transport_framework::particle_tracker ALIAS particle_tracker)

target_include_directories(particle_tracker
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/third_party/nanoflann>
)

target_compile_features(particle_tracker INTERFACE cxx_std_20)
target_compile_definitions(particle_tracker INTERFACE PK_IO_USE_ZLIB)

# ---- Dependencies ----
find_package(ZLIB REQUIRED)
find_package(TBB REQUIRED)
if(NOT TARGET TBB::tbb)
  message(FATAL_ERROR "TBB was found, but the TBB::tbb target is unavailable")
endif()
find_package(Boost REQUIRED CONFIG)

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

target_link_libraries(particle_tracker
  INTERFACE
    ZLIB::ZLIB
    TBB::tbb
    Boost::boost
    ${PK_NETCDF_C_LIB}
    ${PK_NETCDF_CXX_LIB}
)

if(PK_NETCDF_INCLUDE_DIRS)
  target_include_directories(particle_tracker INTERFACE ${PK_NETCDF_INCLUDE_DIRS})
endif()

# ---- Examples ----
if(PARTICLE_TRACKER_BUILD_EXAMPLES)
  add_executable(hpc_particles examples/hpc_particles.cpp)
  target_link_libraries(hpc_particles
    PRIVATE particle_transport_framework::particle_tracker)

  add_executable(hpc_trajectories examples/hpc_trajectories.cpp)
  target_link_libraries(hpc_trajectories
    PRIVATE particle_transport_framework::particle_tracker)
endif()

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
