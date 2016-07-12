cmake_minimum_required(VERSION 2.8)

project(BoostDateTime)

include_directories(${CMAKE_SOURCE_DIR}/../..)



FILE(GLOB temp_files src/gregorian/*.hpp src/gregorian/*.cpp src/posix_time/*.hpp src/posix_time/*.cpp)
SOURCE_GROUP( dt FILES ${temp_files} )
LIST(APPEND src_files ${temp_files})

add_library(boost_date_time ${src_files})
