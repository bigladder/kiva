cmake_minimum_required (VERSION 2.8)
project (kiva)
enable_testing()

option(KIVA_COVERAGE "Add coverage reports" OFF)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if (KIVA_COVERAGE)
  find_package(codecov)
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

set(PROCESSOR_ARCHITECTURE "x86_64" CACHE STRING "Set architecture used for version string.")

# Set a default build type if none was specified
if(NOT PROCESSOR_ARCHITECTURE)
  message(STATUS "Setting processor architecture to 'x86_64' as none was specified.")
  set(PROCESSOR_ARCHITECTURE x86_64 CACHE STRING "Choose the processor architecture." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE PROCESSOR_ARCHITECTURE PROPERTY STRINGS "x86_64" "x86")
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(build_architecture "macosx-${CMAKE_CXX_COMPILER_ID}-${PROCESSOR_ARCHITECTURE}")
else()
  set(build_architecture "${CMAKE_SYSTEM_NAME}-${CMAKE_CXX_COMPILER_ID}-${PROCESSOR_ARCHITECTURE}")
endif()

string(TOLOWER ${build_architecture} build_architecture)

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )

option(KIVA_3D "Enable 3D Calculations" ON)
option(ENABLE_OPENMP "Use OpenMP" OFF)

set(KIVA_EXE_BUILD ON)

if(CMAKE_COMPILER_IS_GNUCXX OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  if (UNIX)
    add_definitions("-fPIC")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-register -Wno-deprecated-declarations")
endif()

add_definitions("-DBOOST_ALL_NO_LIB")

include_directories( ${CMAKE_SOURCE_DIR}/src/libkiva/)
include_directories( ${CMAKE_SOURCE_DIR}/src/libgroundplot/)
include_directories( SYSTEM ${CMAKE_SOURCE_DIR}/vendor/boost-1.61.0/)
include_directories( SYSTEM ${CMAKE_SOURCE_DIR}/vendor/mathgl-2.3.5.1/include/)
include_directories( SYSTEM ${CMAKE_BINARY_DIR}/vendor/mathgl-2.3.5.1/include/)
include_directories( SYSTEM ${CMAKE_SOURCE_DIR}/vendor/yaml-cpp-0.5.3/include/)
include_directories( SYSTEM ${CMAKE_SOURCE_DIR}/vendor/eigen-da9b4e14c255/)
include_directories( SYSTEM ${CMAKE_SOURCE_DIR}/vendor/gtest/include/)

set(KIVA_STATIC_LIB OFF CACHE BOOL "" FORCE)

add_subdirectory(vendor)
add_subdirectory(src/libkiva)
add_subdirectory(src/libgroundplot)
add_subdirectory(src/kiva)

# Testing
add_subdirectory(test)
if (KIVA_COVERAGE)
  coverage_evaluate()
endif()
