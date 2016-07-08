if (${ENABLE_ITERATIVE_SOLVER})
  add_definitions(-D USE_LIS_SOLVER )
endif (${ENABLE_ITERATIVE_SOLVER})

if (${ENABLE_OPENGL})
  add_definitions(-D ENABLE_OPENGL )
endif (${ENABLE_OPENGL})

set(kiva_src Algorithms.c++
             Algorithms.h
             Domain.c++
             Domain.h
             Functions.c++
             Functions.h
             Geometry.c++
             Geometry.h
             Ground.c++
             Ground.h
             GroundPlot.c++
             GroundPlot.h
             Input.c++
             Input.h
             Mesher.c++
             Mesher.h
             WeatherData.c++
             WeatherData.h )

if (${ENABLE_OPENGL})
  set(kiva_src ${kiva_src} PixelCounter.c++ PixelCounter.h )
endif (${ENABLE_OPENGL})

link_directories(${CMAKE_PREFIX_PATH}/lib)
include_directories(${CMAKE_PREFIX_PATH}/include ${CMAKE_PREFIX_PATH}/include/boost-numeric-bindings)

add_library(kivalib ${kiva_src})

set(kiva_link_flags "-fopenmp")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if (${ENABLE_OPENGL})
    set(kiva_link_flags "${kiva_link_flags} -framework OpenGL -framework Cocoa -framework IOKit" )
  endif (${ENABLE_OPENGL})
endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${kiva_link_flags}")

set(kivalib_lib boost_date_time
                boost_filesystem
                boost_system
                mgl )

if (${ENABLE_ITERATIVE_SOLVER})
  set(kivalib_lib ${kivalib_lib} lis )
else (${ENABLE_ITERATIVE_SOLVER})
  set(kivalib_lib ${kivalib_lib} umfpack
                                 amd
                                 cholmod
                                 colamd
                                 suitesparseconfig )
endif (${ENABLE_ITERATIVE_SOLVER})

if (${ENABLE_OPENGL})
  set(kivalib_lib ${kivalib_lib} GLEW glfw3 )
endif (${ENABLE_OPENGL})

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  if (${ENABLE_OPENGL})
    set(kivalib_lib ${kivalib_lib} pthread
                                   Xrandr
                                   Xxf86vm
                                   Xi
                                   X11
                                   GL )
  endif (${ENABLE_OPENGL})
  if (NOT ${ENABLE_ITERATIVE_SOLVER})
    set (kivalib_lib ${kivalib_lib} openblas)
  endif (NOT ${ENABLE_ITERATIVE_SOLVER})
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

target_link_libraries(kivalib ${kivalib_lib})
