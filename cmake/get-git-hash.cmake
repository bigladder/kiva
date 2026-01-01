# See the LICENSE.txt file for additional terms and conditions.
# Adapted from https://jonathanhamberg.com/post/cmake-embedding-git-hash/
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE PROJECT_GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
