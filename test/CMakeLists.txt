add_subdirectory(unit)

add_custom_target(results_directory ALL COMMAND ${CMAKE_COMMAND} -E make_directory results)

add_custom_target(rm_results_subdirectory ALL COMMAND ${CMAKE_COMMAND} -E remove_directory results/${build_architecture} DEPENDS results_directory)
add_custom_target(results_subdirectory ALL COMMAND ${CMAKE_COMMAND} -E make_directory results/${build_architecture} DEPENDS rm_results_subdirectory)

set(SCRIPT_DIR "${CMAKE_SOURCE_DIR}/scripts")

function( add_example_case )
  set(options)
  set(oneValueArgs IN_FILE EPW_FILE)
  set(multiValueArgs)
  cmake_parse_arguments(INT_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set(TEST_NAME "integration.example.${INT_TEST_IN_FILE}")

  add_custom_target(results_${INT_TEST_IN_FILE}_directory ALL COMMAND ${CMAKE_COMMAND} -E make_directory results/${build_architecture}/${TEST_NAME} DEPENDS results_subdirectory)

  add_test(NAME ${TEST_NAME} COMMAND ruby run-kiva.rb
    $<TARGET_FILE:kiva>
    ${CMAKE_SOURCE_DIR}/examples/${INT_TEST_IN_FILE}.yaml
    ${CMAKE_SOURCE_DIR}/weather/${INT_TEST_EPW_FILE}.epw
    ${CMAKE_CURRENT_BINARY_DIR}/results/${build_architecture}/${TEST_NAME}/out.csv
    WORKING_DIRECTORY ${SCRIPT_DIR})

endfunction()

function( add_bestest_case )
  set(options)
  set(oneValueArgs IN_FILE EPW_FILE)
  set(multiValueArgs)
  cmake_parse_arguments(INT_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set(TEST_NAME "integration.BESTEST.${INT_TEST_IN_FILE}")

  add_custom_target(results_${INT_TEST_IN_FILE}_directory ALL COMMAND ${CMAKE_COMMAND} -E make_directory results/${build_architecture}/${TEST_NAME} DEPENDS results_subdirectory)

  add_test(NAME ${TEST_NAME} COMMAND ruby run-kiva.rb
    $<TARGET_FILE:kiva>
    ${CMAKE_SOURCE_DIR}/test/bestest/${INT_TEST_IN_FILE}.yaml
    ${CMAKE_SOURCE_DIR}/test/bestest/weather/${INT_TEST_EPW_FILE}.epw
    ${CMAKE_CURRENT_BINARY_DIR}/results/${build_architecture}/${TEST_NAME}/out.csv
    WORKING_DIRECTORY ${SCRIPT_DIR})

endfunction()


add_example_case( IN_FILE "slab" EPW_FILE "USA_DC_Washington")
add_example_case( IN_FILE "basement" EPW_FILE "USA_IL_Chicago")
add_example_case( IN_FILE "crawlspace" EPW_FILE "USA_FL_Tampa")

add_bestest_case( IN_FILE "GC10a" EPW_FILE "GCSS-W40TM2")
add_bestest_case( IN_FILE "GC30a" EPW_FILE "GCSS-W40TM2")
add_bestest_case( IN_FILE "GC30b" EPW_FILE "GCSS-W20TM2")
add_bestest_case( IN_FILE "GC30c" EPW_FILE "GCSS-W40TM2")
add_bestest_case( IN_FILE "GC40a" EPW_FILE "GCSP-W40TM2")
add_bestest_case( IN_FILE "GC40b" EPW_FILE "GCSP-W20TM2")
add_bestest_case( IN_FILE "GC40c" EPW_FILE "GCSP-W40TM2")
add_bestest_case( IN_FILE "GC45b" EPW_FILE "GCSP-W20TM2")
add_bestest_case( IN_FILE "GC45c" EPW_FILE "GCSP-W40TM2")
add_bestest_case( IN_FILE "GC50b" EPW_FILE "GCSP-W20TM2")
add_bestest_case( IN_FILE "GC55b" EPW_FILE "GCSP-W20TM2")
add_bestest_case( IN_FILE "GC55c" EPW_FILE "GCSP-W40TM2")
add_bestest_case( IN_FILE "GC60b" EPW_FILE "GCSS-W20TM2")
add_bestest_case( IN_FILE "GC65b" EPW_FILE "GCSS-W01TM2")
add_bestest_case( IN_FILE "GC70b" EPW_FILE "GCSP-W01TM2")
add_bestest_case( IN_FILE "GC80b" EPW_FILE "GCSP-W20TM2")
add_bestest_case( IN_FILE "GC80c" EPW_FILE "GCSP-W40TM2")
