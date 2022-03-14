project(groundplot)

set(groundplot_src  GroundPlot.cpp
                    GroundPlot.hpp )

add_library(groundplot STATIC ${groundplot_src})

target_link_libraries(groundplot fmt mgl-static libkiva)

target_include_directories(groundplot PUBLIC "${MathGL2_SOURCE_DIR}/include" "${MathGL2_BINARY_DIR}/include")

if (MSVC)
  set_source_files_properties(${groundplot_src} PROPERTIES COMPILE_FLAGS "/wd4100 /wd4458 /wd4702")
else()
  set_source_files_properties(${groundplot_src} PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter")
endif()

target_compile_options(groundplot PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>                              # Eventually add /WX
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>  # Eventually add -Werror
)

if (KIVA_COVERAGE)
  add_coverage(groundplot)
endif()
