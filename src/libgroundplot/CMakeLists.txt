project(groundplot)

set(groundplot_src  GroundPlot.cpp
                    GroundPlot.hpp )

add_library(groundplot STATIC ${groundplot_src})

target_link_libraries(groundplot fmt mgl-static libkiva)

target_include_directories(groundplot PUBLIC "${MathGL2_SOURCE_DIR}/include" "${MathGL2_BINARY_DIR}/include")

target_compile_options(groundplot PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>    # Eventually add /WX
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>: # Checks for GCC and Clang
  -Wall -Wextra -Wpedantic -Wno-unknown-pragmas>         # Eventually add -Werror
)

if (KIVA_COVERAGE)
  add_coverage(groundplot)
endif()
