project(groundplot)

set(kiva_src GroundPlot.cpp
             GroundPlot.hpp )

add_library(groundplot STATIC ${kiva_src})

target_link_libraries(groundplot fmt mgl-static libkiva)

target_include_directories(groundplot PUBLIC "${MathGL2_SOURCE_DIR}/include" "${MathGL2_BINARY_DIR}/include")

if (MSVC)
  # Remove "potential loss of data" warning in Boost
  set_source_files_properties(${kiva_src} PROPERTIES COMPILE_FLAGS "/wd4100 /wd4458 /wd4702")
endif()

target_compile_options(groundplot PUBLIC -W4)

if (KIVA_COVERAGE)
  add_coverage(groundplot)
endif()
