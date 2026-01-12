project(groundplot)

set(groundplot_src  GroundPlot.cpp
                    GroundPlot.hpp )

find_package(fmt REQUIRED)

add_library(groundplot STATIC ${groundplot_src})

set(enable-lgpl ON)
find_package(MathGL2 CONFIG REQUIRED)

# if (MSVC)
#     target_compile_options(mgl PRIVATE /wd4477 /wd4334)
#     target_compile_options(mgl-static PRIVATE /wd4477 /wd4334)
# else ()
#     target_compile_options(mgl PRIVATE -Wno-deprecated-declarations -Wno-tautological-constant-out-of-range-compare -Wno-format-security -Wno-unused-result)
#     target_compile_options(mgl-static PRIVATE -Wno-deprecated-declarations -Wno-tautological-constant-out-of-range-compare -Wno-format-security -Wno-unused-result)
# endif ()

# set_target_properties(mgl PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

target_link_libraries(groundplot PUBLIC fmt::fmt libkiva PRIVATE kiva_common_interface)
target_link_libraries(groundplot PRIVATE $<IF:$<TARGET_EXISTS:mgl>, mgl, mgl-static>)

if (KIVA_COVERAGE)
  add_coverage(groundplot)
endif()
