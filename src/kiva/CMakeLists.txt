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

if (KIVA_EXE_BUILD)
    set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)
    # set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    # set(YAML_CPP_CLANG_FORMAT_EXE OFF CACHE BOOL "" FORCE)
    # set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
    # set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
    # set(YAML_MSVC_SHARED_RT OFF CACHE BOOL "" FORCE)
    find_package(yaml-cpp CONFIG REQUIRED)

    if (NOT MSVC)
        target_compile_options(yaml-cpp INTERFACE -Wno-shadow -Wno-effc++)
    endif ()
endif ()

target_link_libraries(kiva PRIVATE groundplot boost_program_options yaml-cpp::yaml-cpp kiva_common_interface)

target_include_directories(kiva PRIVATE "${groundplot_SOURCE_DIR}")

if (KIVA_COVERAGE)
  add_coverage(kiva)
endif()

install(TARGETS kiva DESTINATION bin/)
