if ( MSVC AND NOT ( CMAKE_CXX_COMPILER_ID STREQUAL "Intel" ) )
  add_definitions("-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING")
  set(gtest_force_shared_crt ON)
endif()

add_subdirectory(fmt)

if (KIVA_EXE_BUILD OR KIVA_GROUND_PLOT)
  add_subdirectory(zlib)
  add_subdirectory(libpng)
  add_subdirectory(mathgl-2.4.4)
  if (KIVA_EXE_BUILD)
    add_subdirectory(boost-1.77.0)
    add_subdirectory(yaml-cpp)
  endif()
endif()

if (KIVA_TESTING)
  if (NOT TARGET gtest)
    # Prevent GoogleTest from overriding our compiler/linker options
    # when building with Visual Studio
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(BUILD_GTEST ON CACHE BOOL "" FORCE MARK)
    set(BUILD_GMOCK ON CACHE BOOL "" FORCE)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    mark_as_advanced(BUILD_GTEST BUILD_GMOCK INSTALL_GTEST)
    add_subdirectory(googletest)
  endif()
endif()
