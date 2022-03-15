project(groundplot)

set(groundplot_src  GroundPlot.cpp
                    GroundPlot.hpp )

add_library(groundplot STATIC ${groundplot_src})

target_link_libraries(groundplot fmt mgl-static libkiva)

target_include_directories(groundplot PUBLIC "${MathGL2_SOURCE_DIR}/include" "${MathGL2_BINARY_DIR}/include")

if (MSVC)
  set_source_files_properties(GroundPlot.cpp PROPERTIES COMPILE_FLAGS "/wd4100 /wd4458 /wd4702")
elseif (IOS)
  set_source_files_properties(GroundPlot.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter")
elseif (UNIX)
  # unknown-pragmas comes from using pragma omp parallel. It is a knonw issue. Possible solution will be to add '-fopenmp'
  set_source_files_properties(GroundPlot.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter -Wno-misleading-indentation -Wno-maybe-uninitialized -Wno-unknown-pragmas")
endif()

target_compile_options(groundplot PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>    # Eventually add /WX
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>: # Checks for GCC and Clang
  -Wall -Wextra -Wpedantic>         # Eventually add -Werror
)

if (KIVA_COVERAGE)
  add_coverage(groundplot)
endif()
