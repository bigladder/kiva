execute_process(COMMAND git describe --tags --long
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                RESULT_VARIABLE git_exit_status
                OUTPUT_VARIABLE kiva_version
                OUTPUT_STRIP_TRAILING_WHITESPACE)

if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set (kiva_version "${kiva_version}-Debug")
endif()

if(${git_exit_status} MATCHES "0")
  set(VERSION_NUMBER ${kiva_version})
  message("Building Kiva version: ${VERSION_NUMBER}")
else()
  message("Error finding version number.")
endif()

string(TIMESTAMP YEAR %Y)

configure_file (
    "${PROJECT_SOURCE_DIR}/src/kiva/Version.h.in"
    "${PROJECT_SOURCE_DIR}/src/kiva/Version.h"
)

set(kiva_src Main.c++
             InputParser.c++
             InputParser.h
             Version.h )

add_executable(kiva ${kiva_src})

set(kiva_link_flags "")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(${ENABLE_OPENMP})
  set(kiva_link_flags "${kiva_link_flags} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${kiva_link_flags}")

set(links boost_date_time
             boost_program_options
             yaml-cpp
             lis )

target_link_libraries(kiva ${links} libkiva)
