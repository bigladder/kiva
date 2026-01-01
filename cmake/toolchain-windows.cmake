# SPDX-FileCopyrightText: © 2025 Big Ladder Software <info@bigladdersoftware.com>
# SPDX-License-Identifier: BSD-3-Clause

if(CMAKE_SYSTEM_NAME STREQUAL Windows)

    # Set the Windows SDK version (e.g., 10.0.22000.0)
    # CMAKE_SYSTEM_VERSION (defaults to CMAKE_HOST_SYSTEM_VERSION) is ignored in CMake 3.27+, because
    # it actually corresponds to the Windows build version, which may or may not exactly match an SDK version

    # CMAKE_GENERATOR_PLATFORM is strict; it will issue an error if it cannot be matched
    set(CMAKE_GENERATOR_PLATFORM version=10.0.26100)

    # Set the platform toolset (e.g., v143 for Visual Studio 2022)
    # Note: v142 supports a "version" keyword; v143 no longer does
    set(CMAKE_GENERATOR_TOOLSET "v143")
    set(REQUESTED_CPP_COMPLIER_VERSION "19.44.35213.0" CACHE INTERNAL "")

endif()
