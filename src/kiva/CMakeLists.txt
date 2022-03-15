set(kiva_src Main.cpp
             Input.cpp
             Input.hpp
             InputParser.cpp
             InputParser.hpp
             Simulator.cpp
             Simulator.hpp
             WeatherData.cpp
             WeatherData.hpp )

set(CMAKE_INSTALL_RPATH "@executable_path")

add_executable(kiva ${kiva_src})

target_link_libraries(kiva libkiva groundplot boost_program_options yaml-cpp)

target_include_directories(kiva PRIVATE "${groundplot_SOURCE_DIR}")

set(kiva_suppress_warnings Main.cpp
                           Input.cpp
                           InputParser.cpp
                           Simulator.cpp)
if (MSVC)
  set_source_files_properties(${kiva_suppress_warnings} PROPERTIES COMPILE_FLAGS /wd4458)
elseif (APPLE)
  set_source_files_properties(${kiva_suppress_warnings} PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter")
elseif (UNIX)
  # unknown-pragmas comes from using pragma omp parallel. It is a knonw issue. Possible solution will be to add '-fopenmp'
  #set_source_files_properties(${kiva_suppress_warnings} PROPERTIES COMPILE_FLAGS "-Wno-unknown-pragmas")
endif()

target_compile_options(kiva PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>    # Eventually add /WX
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>: # Checks for GCC and Clang
  -Wall -Wextra -Wpedantic>         # Eventually add -Werror
)

if (KIVA_COVERAGE)
  add_coverage(kiva)
endif()

install(TARGETS kiva DESTINATION bin/)
