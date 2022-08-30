project(groundplot)

set(groundplot_src  GroundPlot.cpp
                    GroundPlot.hpp )

add_library(groundplot STATIC ${groundplot_src})

target_link_libraries(groundplot PUBLIC fmt mgl-static libkiva PRIVATE kiva_common_interface)

target_include_directories(groundplot PUBLIC "${MathGL2_SOURCE_DIR}/include" "${MathGL2_BINARY_DIR}/include")

if (KIVA_COVERAGE)
  add_coverage(groundplot)
endif()
