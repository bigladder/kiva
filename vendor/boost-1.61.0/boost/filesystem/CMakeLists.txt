cmake_minimum_required(VERSION 2.8)

project(BoostFilesystem)

include_directories(${CMAKE_SOURCE_DIR}/../..)

FILE(GLOB temp_files src/*.hpp src/*.cpp)
SOURCE_GROUP( fs FILES ${temp_files} )
LIST(APPEND src_files ${temp_files})

add_library(boost_filesystem ${src_files})

if( MSVC )
  set_target_properties(boost_filesystem PROPERTIES COMPILE_FLAGS "/W0")
else()
  set_target_properties(boost_filesystem PROPERTIES COMPILE_FLAGS "-w")
endif()
