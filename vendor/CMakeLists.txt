include(initialize-submodules)

add_subdirectory(boost-1.88.0) #set by update-boost.py

if (KIVA_EXE_BUILD OR KIVA_GROUND_PLOT)

    set(SKIP_INSTALL_ALL ON)
    set(INSTALL_UTILS OFF)
    # if (NOT TARGET miniziplib)
    #     set(ZLIB_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    #     set(ZLIB_COMPAT ON CACHE BOOL "" FORCE)
    find_package(zlib-ng)
    #     set(ZLIB_LIBRARIES zlibstatic)
    #     set(ZLIB_INCLUDE_DIRS "${zlib_SOURCE_DIR}" "${zlib_BINARY_DIR}")
    #     set(ZLIB_INCLUDE_DIR ${ZLIB_INCLUDE_DIRS})
    #     set(ZLIB_LIBRARY ${ZLIB_LIBRARIES})
    #     if (MSVC)
    #         target_compile_options(zlibstatic PRIVATE /wd5105)
    #         target_compile_options(zlib PRIVATE /wd5105)
    #     endif ()
    #     set_target_properties(zlib PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
    # else ()
    #     set(ZLIB_LIBRARIES miniziplib)
    #     set(ZLIB_INCLUDE_DIRS "${miniziplib_SOURCE_DIR}" "${miniziplib_BINARY_DIR}")
    #     set(ZLIB_INCLUDE_DIR ${ZLIB_INCLUDE_DIRS})
    #     set(ZLIB_LIBRARY ${ZLIB_LIBRARIES})
    # endif ()
    set(PNG_BUILD_ZLIB ON CACHE BOOL "" FORCE)

    set(PNG_TESTS OFF CACHE BOOL "" FORCE)
    set(PNG_SHARED OFF CACHE BOOL "" FORCE)
    find_package(PNG CONFIG REQUIRED)

    set(enable-openmp OFF CACHE BOOL "" FORCE)
    set(PNG_LIBRARY png_static)
    set(PNG_INCLUDE_DIR "${libpng_SOURCE_DIR}" "${libpng_BINARY_DIR}")
    set(PNG_LIBRARIES ${PNG_LIBRARY})
    set(PNG_INCLUDE_DIRS ${PNG_INCLUDE_DIR})

endif ()

if (KIVA_TESTING)
    # if (NOT TARGET gtest)
    #     # Prevent GoogleTest from overriding our compiler/linker options
    #     # when building with Visual Studio
    #     set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    #     set(BUILD_GTEST ON CACHE BOOL "" FORCE MARK)
    #     set(BUILD_GMOCK ON CACHE BOOL "" FORCE)
    #     set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    #     mark_as_advanced(BUILD_GTEST BUILD_GMOCK INSTALL_GTEST)
    #     find_package(GTest)
    # endif ()
    find_package(GTest CONFIG REQUIRED)
endif ()
