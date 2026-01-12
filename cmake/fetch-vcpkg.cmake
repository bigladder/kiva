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

