# Lis: a Library of Iterative Solvers for Linear Systems (http://www.ssisc.org/lis/)

cmake_minimum_required(VERSION 2.8.7)
project( Lis )

if( NOT WIN32)
	enable_language(Fortran)
  INCLUDE(CMakeFortranInformation)
  INCLUDE(FortranCInterface)
endif()
INCLUDE(cmake/CheckTypeSizes.cmake)
INCLUDE(cmake/FindIncludeHeader.cmake)

######################################################################
# Options

OPTION(LIS_ENABLE_OPENMP "Use OpenMP" OFF)
OPTION(LIS_ENABLE_MPI "Use MPI" OFF)
OPTION(LIS_ENABLE_FORTRAN "build Fortran 77 interfaces" OFF)
OPTION(LIS_ENABLE_SAAMG "build SA-AMG preconditioner" OFF)
OPTION(LIS_ENABLE_QUAD "enable quadruple precision operations" OFF)
OPTION(LIS_ENABLE_SSE2 "use Intel Streaming SIMD Extensions" ON)

######################################################################
# General build setting

# Set build directories
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )
# Enable Visual Studio project folder grouping
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
# Set compiler helper variables
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(COMPILER_IS_CLANG TRUE)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(COMPILER_IS_GCC TRUE)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    SET(COMPILER_IS_INTEL TRUE)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    SET(COMPILER_IS_MSVC TRUE)
ENDIF () # CMAKE_CXX_COMPILER_ID
# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

######################################################################
# Lis build settting

ADD_DEFINITIONS(-DHAVE_CONFIG_H)
# Create Configure.h
CONFIGURE_FILE (include/Configure.h.in ${PROJECT_BINARY_DIR}/Configure.h)
INCLUDE_DIRECTORIES (${PROJECT_BINARY_DIR}/)

IF(LIS_ENABLE_OPENMP)
	FIND_PACKAGE( OpenMP REQUIRED )
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
	SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	IF(NOT COMPILER_IS_MSVC)
		SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp" )
	ENDIF()
ENDIF()

IF(LIS_ENABLE_MPI)
	FIND_PACKAGE(MPI)
	IF(MPI_FOUND)
		SET(CMAKE_C_COMPILER ${MPI_COMPILER})
		SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
		ADD_DEFINITIONS(-DUSE_MPI)
	ELSE(MPI_FOUND)
		MESSAGE (FATAL_ERROR "Aborting: MPI implementation is not found!")
	ENDIF(MPI_FOUND)
ENDIF()

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/include/malloc/")
endif()

######################################################################
# Add subdirectories

ADD_SUBDIRECTORY( src )
