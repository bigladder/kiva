# SPDX-FileCopyrightText: © 2025 Big Ladder Software <info@bigladdersoftware.com>
# SPDX-License-Identifier: BSD-3-Clause

macro(git_command args output_variable)
    execute_process(
            COMMAND ${GIT_EXECUTABLE} ${args}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            RESULT_VARIABLE exit_status
            ERROR_VARIABLE error_output
            OUTPUT_VARIABLE temp_variable
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (NOT ${exit_status} MATCHES "0")
        message(SEND_ERROR "Unable to retrieve ${output_variable}.")
        message(FATAL_ERROR ${error_output})
    endif ()
    set(${output_variable} "${temp_variable}")
endmacro()

macro(process_git_version)

    find_package(Git)

    if (GIT_FOUND)

        # Branch
        git_command("rev-parse;--abbrev-ref;HEAD" ${PROJECT_NAME}_git_branch)

        # SHA
        git_command("rev-parse;--verify;--short;HEAD" ${PROJECT_NAME}_git_sha)

        # Tag
        git_command("tag" tag_list)

        if (tag_list MATCHES "^$")
            # Use v0.0.0 if there aren't tags yet
            set(${PROJECT_NAME}_git_tag "v0.0.0")
            git_command("rev-list;--count;HEAD" ${PROJECT_NAME}_git_build_number)
        else ()
            git_command("describe;--tags;--abbrev=0" ${PROJECT_NAME}_git_tag)

            # Build number
            git_command("rev-list;--count;HEAD;^${${PROJECT_NAME}_git_tag}" ${PROJECT_NAME}_git_build_number)

        endif ()

        # Status (dirty?)
        git_command("diff;--shortstat" ${PROJECT_NAME}_git_status)

        if (${PROJECT_NAME}_git_status MATCHES "changed")
            set(${PROJECT_NAME}_git_status "dirty")
        else ()
            unset(${PROJECT_NAME}_git_status)
        endif ()

        # Get version components
        set(semver_normal_regex "v(0|[1-9][0-9]*)\\.(0|[1-9][0-9]*)\\.(0|[1-9][0-9]*)")
        set(semver_prerelease_regex "(-((0|[1-9][0-9]*|[0-9]*[a-zA-Z-][0-9a-zA-Z-]*)(\\.(0|[1-9][0-9]*|[0-9]*[a-zA-Z-][0-9a-zA-Z-]*))*))")
        set(semver_regex "^${semver_normal_regex}${semver_prerelease_regex}?$")
        if (NOT ${PROJECT_NAME}_git_tag MATCHES ${semver_regex})
            message(FATAL_ERROR "Tag, \"${${PROJECT_NAME}_git_tag}\", does not conform to SemVer \"${${PROJECT_NAME}_semver_regex}\".")
        endif ()

        string(REGEX REPLACE ${semver_regex} "\\1" ${PROJECT_NAME}_version_major "${${PROJECT_NAME}_git_tag}")
        string(REGEX REPLACE ${semver_regex} "\\2" ${PROJECT_NAME}_version_minor "${${PROJECT_NAME}_git_tag}")
        string(REGEX REPLACE ${semver_regex} "\\3" ${PROJECT_NAME}_version_patch "${${PROJECT_NAME}_git_tag}")
        if (${PROJECT_NAME}_git_tag MATCHES "^${semver_normal_regex}${semver_prerelease_regex}$")
            string(REGEX REPLACE ${semver_regex} "\\4" ${PROJECT_NAME}_version_prerelease "${${PROJECT_NAME}_git_tag}")
        else ()
            set(${PROJECT_NAME}_version_prerelease "")
        endif ()

        # Build meta info
        if (NOT ${PROJECT_NAME}_git_build_number MATCHES "^0$")
            set(${PROJECT_NAME}_version_meta "+${${PROJECT_NAME}_git_branch}.${${PROJECT_NAME}_git_sha}.${${PROJECT_NAME}_git_build_number}")
            if (DEFINED ${PROJECT_NAME}_git_status)
                set(${PROJECT_NAME}_version_meta "${${PROJECT_NAME}_version_meta}.${${PROJECT_NAME}_git_status}")
            endif ()
        else ()
            if (DEFINED ${PROJECT_NAME}_git_status)
                set(${PROJECT_NAME}_version_meta "+${${PROJECT_NAME}_git_status}")
            else ()
                set(${PROJECT_NAME}_version_meta "")
            endif ()
        endif ()


    else ()

        message(WARNING "Failed to find git executable. Unable to establish version information for ${PROJECT_NAME}.")
        set(${PROJECT_NAME}_version_major "0")
        set(${PROJECT_NAME}_version_minor "0")
        set(${PROJECT_NAME}_version_patch "0")
        set(${PROJECT_NAME}_version_meta "+git.not.found")

    endif ()

    set(${PROJECT_NAME}_version "v${${PROJECT_NAME}_version_major}.${${PROJECT_NAME}_version_minor}.${${PROJECT_NAME}_version_patch}${${PROJECT_NAME}_version_prerelease}${${PROJECT_NAME}_version_meta}")

endmacro()

macro(add_version_variable version_variable)
    string(TOUPPER ${version_variable} temp_upper)
    string(APPEND content "#define ${temp_upper} ${${version_variable}}\n")
endmacro()

macro(make_version_header)

    set(content "#pragma once\n\n")

    add_version_variable(${PROJECT_NAME}_git_branch)
    add_version_variable(${PROJECT_NAME}_git_sha)
    add_version_variable(${PROJECT_NAME}_git_tag)
    add_version_variable(${PROJECT_NAME}_git_build_number)
    add_version_variable(${PROJECT_NAME}_git_status)
    add_version_variable(${PROJECT_NAME}_version_major)
    add_version_variable(${PROJECT_NAME}_version_minor)
    add_version_variable(${PROJECT_NAME}_version_patch)
    add_version_variable(${PROJECT_NAME}_version_prerelease)
    add_version_variable(${PROJECT_NAME}_version_meta)
    add_version_variable(${PROJECT_NAME}_version)

    file(CONFIGURE
            OUTPUT "${PROJECT_BINARY_DIR}/include/${PROJECT_NAME}/version.h"
            CONTENT ${content}
    )
endmacro()