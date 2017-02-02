project(groundplot)

set(kiva_src GroundPlot.cpp
             GroundPlot.hpp )

include_directories(${CMAKE_CURRENT_BINARY_DIR}/../libkiva/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libkiva/)

add_library(groundplot STATIC ${kiva_src})

set(links boost_filesystem
          boost_system
          mgl-static )

if( MSVC )
  set_target_properties(mgl-static PROPERTIES COMPILE_FLAGS "/W0")
else()
  set_target_properties(mgl-static PROPERTIES COMPILE_FLAGS "-w")
endif()

target_link_libraries(groundplot ${links} libkiva)

add_coverage(groundplot)
