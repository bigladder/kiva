set(kiva_test_src
  bestest.unit.cpp
  main.cpp
)

include_directories(${CMAKE_BINARY_DIR}/src/libkiva/)

add_executable(kiva_tests ${kiva_test_src})

target_link_libraries(kiva_tests gtest libkiva)

add_test(NAME BESTEST.GC10a COMMAND kiva_tests "--gtest_filter=BESTEST.GC10a")
add_test(NAME BESTEST.GC30a COMMAND kiva_tests "--gtest_filter=BESTEST.GC30a")
