cmake_minimum_required (VERSION 2.8)
project (kiva)
enable_testing()

option(KIVA_COVERAGE "Add coverage reports" OFF)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if (KIVA_COVERAGE)
  set(ENABLE_COVERAGE ON CACHE BOOL "" FORCE)
  find_package(codecov)
endif()

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

add_subdirectory(vendor)
add_subdirectory(src/libkiva)
add_subdirectory(src/libgroundplot)
add_subdirectory(src/kiva)

# Testing
add_subdirectory(test)
if (KIVA_COVERAGE)
  coverage_evaluate()
endif()

include(cmake/Install.cmake)
