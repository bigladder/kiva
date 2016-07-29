cmake_minimum_required(VERSION 2.8)

project(BoostSystem)

include_directories(${CMAKE_SOURCE_DIR}/../..)



FILE(GLOB temp_files src/*.hpp src/*.cpp)
SOURCE_GROUP( s FILES ${temp_files} )
LIST(APPEND src_files ${temp_files})

add_library(boost_system ${src_files})
