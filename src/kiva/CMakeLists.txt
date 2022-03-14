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

if (MSVC)
  # Remove "potential loss of data" warning in Boost
  set_source_files_properties(${kiva_src} PROPERTIES COMPILE_FLAGS /wd4458)
endif()

target_compile_options(kiva PUBLIC -W4)

if (KIVA_COVERAGE)
  add_coverage(kiva)
endif()

install(TARGETS kiva DESTINATION bin/)
