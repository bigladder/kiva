
add_custom_target(results_directory ALL COMMAND ${CMAKE_COMMAND} -E make_directory results)

add_custom_target(results_subdirectory ALL COMMAND ${CMAKE_COMMAND} -E make_directory results/${build_architecture})

set(SCRIPT_DIR "${CMAKE_SOURCE_DIR}/scripts")
set(ENV{BUNDLE_GEMFILE} "${SCRIPT_DIR}/Gemfile")

function( add_integration_test )
  set(options)
  set(oneValueArgs IN_FILE EPW_FILE)
  set(multiValueArgs)
  cmake_parse_arguments(INT_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set(TEST_NAME "integration.${INT_TEST_IN_FILE}")

  add_custom_target(results_${TEST_NAME}_directory ALL COMMAND ${CMAKE_COMMAND} -E make_directory results/${build_architecture}/${TEST_NAME})

  if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    add_test(
      NAME ${TEST_NAME}
      WORKING_DIRECTORY ${SCRIPT_DIR}
      COMMAND C:\Ruby200\bin\bundle exec C:\Ruby200\bin\bundle\ruby run-dummy.rb)
  else()
    add_test(
      NAME ${TEST_NAME}
      WORKING_DIRECTORY ${SCRIPT_DIR}
      COMMAND bundle exec ruby run-kiva.rb
        ${CMAKE_BINARY_DIR}/kiva
        ${CMAKE_SOURCE_DIR}/examples/${INT_TEST_IN_FILE}.yaml
        ${CMAKE_SOURCE_DIR}/weather/${INT_TEST_EPW_FILE}.epw
        ${CMAKE_CURRENT_BINARY_DIR}/results/${build_architecture}/${TEST_NAME}/out.csv)
  endif()

endfunction()

add_integration_test( IN_FILE "slab" EPW_FILE "USA_DC_Washington")
add_integration_test( IN_FILE "basement" EPW_FILE "USA_IL_Chicago")
add_integration_test( IN_FILE "crawlspace" EPW_FILE "USA_FL_Tampa")
