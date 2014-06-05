set(iterative_solver 1 )

if (${iterative_solver})
  add_definitions(-D USE_LIS_SOLVER )
endif (${iterative_solver})

execute_process(COMMAND git describe --tags --long
                RESULT_VARIABLE git_exit_status
                OUTPUT_VARIABLE kiva_verion
                OUTPUT_STRIP_TRAILING_WHITESPACE)

if(${git_exit_status} MATCHES "0")
    SET(VERSION_NUMBER ${kiva_verion})
    message("Building Kiva version: ${VERSION_NUMBER}")
else()
    message("Error finding version number.")
endif(${git_exit_status} MATCHES "0")
                
configure_file (
    "${PROJECT_SOURCE_DIR}/source/Version.h.in"
    "${PROJECT_SOURCE_DIR}/source/Version.h"
)

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
             PixelCounter.c++
             WeatherData.c++ )

set(kiva_hdr Algorithms.h
             Domain.h
             Functions.h
             Geometry.h
             Ground.h
             GroundPlot.h
             Input.h
             InputParser.h
             Mesher.h
             PixelCounter.h
             WeatherData.h
             Version.h )

link_directories(${CMAKE_PREFIX_PATH}/lib)
include_directories(${CMAKE_PREFIX_PATH}/include ${CMAKE_PREFIX_PATH}/include/boost-numeric-bindings)

add_executable(kiva ${kiva_src} ${kiva_hdr})


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(kiva_link_flags "-fopenmp -framework OpenGL -framework Cocoa -framework IOKit")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(kiva_link_flags "-fopenmp")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${kiva_link_flags}")



set(kiva_lib boost_date_time
             boost_filesystem
             boost_system
             boost_program_options
             mgl
             yaml-cpp
             GLEW
             glfw3 )
             
if (${iterative_solver})
  set(kiva_lib_solver lis )
else (${iterative_solver})
  set(kiva_lib_solver umfpack
                      amd
                      cholmod
                      colamd
                      suitesparseconfig ) 
endif (${iterative_solver})


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(kiva_lib_os )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(kiva_lib_os pthread
                  Xrandr
                  Xxf86vm
                  Xi
                  X11
                  GL )
  if (NOT ${iterative_solver})
    set (kiva_lib_os ${kiva_lib_os} openblas)
  endif (NOT ${iterative_solver})
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")


target_link_libraries(kiva ${kiva_lib} ${kiva_lib_solver} ${kiva_lib_os})



