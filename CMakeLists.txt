cmake_minimum_required (VERSION 2.8)
project (kiva)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are:
		None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(NOT CMAKE_PREFIX_PATH)
  set(CMAKE_PREFIX_PATH "/usr/local" CACHE STRING "Primary prefix path for bin, lib, include directories." FORCE)
endif(NOT CMAKE_PREFIX_PATH)

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )

option( ENABLE_ITERATIVE_SOLVER "Use LIS as an iterative solver instead of UMFPACK" ON )
option( ENABLE_OPENGL "Use OpenGL to perform shading calculations." OFF )

include_directories( ${CMAKE_SOURCE_DIR}/src/core/ )

add_subdirectory(src/core)
add_subdirectory(src/kiva)
