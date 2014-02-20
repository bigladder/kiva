cmake_minimum_required (VERSION 2.8)
project (Kiva)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are:
		None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(NOT CMAKE_PREFIX_PATH)
  set(CMAKE_PREFIX_PATH "/usr/local" CACHE STRING "Primary prefix path for bin, lib, include directories." FORCE)
endif(NOT CMAKE_PREFIX_PATH)

message(${CMAKE_PREFIX_PATH})

set(kiva_src Main.c++
             Algorithms.c++
             Domain.c++
             Functions.c++
             Geometry.c++
             Ground.c++
             GroundPlot.c++
             Input.c++
             InputParser.c++
             Mesher.c++
             WeatherData.c++)

set(kiva_hdr Algorithms.h
             Domain.h
             Functions.h
             Geometry.h
             Ground.h
             GroundPlot.h
             Input.h
             InputParser.h
             Mesher.h
             WeatherData.h)

link_directories(${CMAKE_PREFIX_PATH}/lib)
include_directories(${CMAKE_PREFIX_PATH}/include ${CMAKE_PREFIX_PATH}/include/boost-numeric-bindings)

add_executable(Kiva ${kiva_src} ${kiva_hdr})


set(kiva_lib boost_date_time
             boost_thread
             boost_filesystem
             boost_system
             boost_program_options
             lis
             umfpack
             amd
             cholmod
             openblas
             colamd
             suitesparseconfig
             mgl
             yaml-cpp
             GLEW
             glfw3
             pthread
             Xrandr
             Xxf86vm
             Xi
             X11
             GL)

target_link_libraries(Kiva ${kiva_lib})
                    
         
                    
