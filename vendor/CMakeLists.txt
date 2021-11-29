include(initialize-submodules)
initialize_submodules()

add_subdirectory(fmt)

if (KIVA_EXE_BUILD OR KIVA_GROUND_PLOT)

  set(ZLIB_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
  add_subdirectory(zlib-ng)

  set(PNG_TESTS OFF CACHE BOOL "" FORCE)
  set(PNG_SHARED OFF CACHE BOOL "" FORCE)
  add_subdirectory(libpng)

  set(enable-openmp OFF CACHE BOOL "" FORCE)
  add_subdirectory(mathgl-2.4.4)

  if (KIVA_EXE_BUILD)
    add_subdirectory(boost-1.77.0)

    set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)
    set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(YAML_CPP_CLANG_FORMAT_EXE OFF CACHE BOOL "" FORCE)
    set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
    set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
    set(YAML_MSVC_SHARED_RT OFF CACHE BOOL "" FORCE)
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
