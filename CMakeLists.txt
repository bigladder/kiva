cmake_minimum_required (VERSION 2.8)
project (kiva)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )

option(ENABLE_OPENMP "Use OpenMP" OFF)
option( ENABLE_OPENGL "Use OpenGL to perform shading calculations." OFF )

if(CMAKE_COMPILER_IS_GNUCXX OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  if (UNIX)
    add_definitions("-fPIC")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-register -Wno-deprecated-declarations")
endif()

add_definitions("-DBOOST_ALL_NO_LIB")

include_directories( ${CMAKE_SOURCE_DIR}/src/libkiva/)
include_directories( ${CMAKE_SOURCE_DIR}/vendor/boost-1.61.0/)
include_directories( ${CMAKE_SOURCE_DIR}/vendor/lis-1.5.66/include/)
include_directories( ${CMAKE_SOURCE_DIR}/vendor/mathgl-2.3.5.1/include/)
include_directories( ${CMAKE_BINARY_DIR}/vendor/mathgl-2.3.5.1/include/)
include_directories( ${CMAKE_SOURCE_DIR}/vendor/yaml-cpp-0.5.3/include/)

add_subdirectory(vendor)
add_subdirectory(src/libkiva)
add_subdirectory(src/kiva)
