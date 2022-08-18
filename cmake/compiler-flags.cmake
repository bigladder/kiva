# Empty default flags
# set(CMAKE_C_FLAGS "")
# set(CMAKE_CXX_FLAGS "")
# set(CMAKE_CXX_FLAGS_RELEASE "")
# set(CMAKE_CXX_FLAGS_DEBUG "")
# set(CMAKE_EXE_LINKER_FLAGS "")
# set(CMAKE_EXE_LINKER_FLAGS_RELEASE "")
# set(CMAKE_EXE_LINKER_FLAGS_DEBUG "")

add_library(kiva_common_interface INTERFACE)

  #================#
  # Compiler flags #
  #================#

target_compile_options(kiva_common_interface INTERFACE
$<$<CXX_COMPILER_ID:MSVC>:
  /W4
  $<$<BOOL:${KIVA_WERROR}>:
    /WX
  >
>
$<$<CXX_COMPILER_ID:GNU>:
  -Wall 
  -Wextra 
  -Wpedantic 
  $<$<BOOL:${KIVA_WERROR}>:
    -Werror
  >
>

$<$<CXX_COMPILER_ID:Clang>:
  -Wall 
  -Wextra 
  -Wpedantic
  $<$<BOOL:${KIVA_WERROR}>:
    -Werror
  >
>
)