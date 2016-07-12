if(${ENABLE_OPENGL})
  add_definitions(-D ENABLE_OPENGL )
endif()

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
endif ()

add_library(libkiva ${kiva_src})
set_target_properties(libkiva PROPERTIES OUTPUT_NAME kiva)

set(kiva_link_flags "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${ENABLE_OPENGL})
    set(kiva_link_flags "${kiva_link_flags} -framework OpenGL -framework Cocoa -framework IOKit" )
  endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(${ENABLE_OPENMP})
  set(kiva_link_flags "${kiva_link_flags} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${kiva_link_flags}")

set(links boost_date_time
          boost_filesystem
          boost_system
          mgl
          lis )

if(${ENABLE_OPENGL})
  set(links ${links} GLEW glfw3 )
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  if(${ENABLE_OPENGL})
    set(links ${links} pthread
                       Xrandr
                       Xxf86vm
                       Xi
                       X11
                       GL )
  endif()
endif()

target_link_libraries(libkiva ${links})
