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

target_compile_options(kiva PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
  -Wall -Wextra -Wpedantic -Werror>
)

if (KIVA_COVERAGE)
  add_coverage(kiva)
endif()

install(TARGETS kiva DESTINATION bin/)
