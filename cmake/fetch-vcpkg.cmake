include(FetchContent)

MESSAGE(STATUS "Fetching vcpkg...")

# FetchContent_Populate, with SOURCE_DIR specified, is being deprecated.
FetchContent_Declare(
    vcpkg
    GIT_REPOSITORY https://github.com/microsoft/vcpkg.git
    GIT_TAG        "2025.12.12"
    OVERRIDE_FIND_PACKAGE
)
FetchContent_MakeAvailable(vcpkg)
MESSAGE(STATUS "vcpkg source downloaded to: ${vcpkg_SOURCE_DIR}")

if(WIN32)
    set(VCPKG_BOOTSTRAP_SCRIPT "bootstrap-vcpkg.bat" CACHE STRING "")
else()
    set(VCPKG_BOOTSTRAP_SCRIPT "bootstrap-vcpkg.sh" CACHE STRING "")
endif()

if(NOT EXISTS "${vcpkg_SOURCE_DIR}/vcpkg${CMAKE_EXECUTABLE_SUFFIX}")
    execute_process(
        COMMAND "${VCPKG_BOOTSTRAP_SCRIPT}"
        WORKING_DIRECTORY "${vcpkg_SOURCE_DIR}"
        RESULT_VARIABLE bootstrap_result
        OUTPUT_VARIABLE bootstrap_output
        ERROR_VARIABLE bootstrap_error
    )
    if(NOT bootstrap_result EQUAL 0)
        message(FATAL_ERROR "vcpkg bootstrap failed: ${bootstrap_error}")
    endif()
endif()

# Could use the following; instead we have a chain of toolchain files in Presets
#list(APPEND CMAKE_PROJECT_TOP_LEVEL_INCLUDES "${vcpkg_SOURCE_DIR}/scripts/buildsystems/vcpkg.cmake")
#list(APPEND CMAKE_TRY_COMPILE_PLATFORM_VARIABLES CMAKE_PROJECT_TOP_LEVEL_INCLUDES)

