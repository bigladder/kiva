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

include_directories(${kiva_BINARY_DIR}/src/libkiva/)
add_executable(kiva ${kiva_src})

set(kiva_link_flags "")

if(${ENABLE_OPENMP})
  set(kiva_link_flags "${kiva_link_flags} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${kiva_link_flags}")

set(links boost_date_time
          boost_program_options
          boost_filesystem
          boost_system
          yaml-cpp )


if( MSVC )
  set_target_properties(yaml-cpp PROPERTIES COMPILE_FLAGS "/W0")
else()
  set_target_properties(yaml-cpp PROPERTIES COMPILE_FLAGS "-w")
endif()

target_link_libraries(kiva ${links} libkiva groundplot)

if (KIVA_COVERAGE)
  add_coverage(kiva)
endif()

install(TARGETS kiva DESTINATION bin/)
