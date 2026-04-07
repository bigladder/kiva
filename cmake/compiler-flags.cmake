# Interface library contains the flags the library will use for each compiler
# Interface library is used instead of CMAKE_CXX flags to prevent it from propagating it to undesired places
add_library(kiva_common_interface INTERFACE)

  #================#
  # Compiler flags #
  #================#

target_compile_options(kiva_common_interface INTERFACE
  $<$<CXX_COMPILER_ID:MSVC>:
    /W4 # Warning level (default is W3)
    $<$<BOOL:${KIVA_WERROR}>:
      /WX # Turn warnings into errors
    >
  >
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:
    -Wall
    -Wextra
    -Wpedantic
    $<$<BOOL:${KIVA_WERROR}>:
      -Werror # Turn warnings into errors
    >
  >
)

# This macro will encapsulate the CMAKE_CXX flags that should only be set for executables
macro(SET_CXX_FLAGS)
  # Remove unwanted CMake defaults from global flags
  if (MSVC)
  # See https://gitlab.kitware.com/cmake/cmake/-/blob/master/Modules/Platform/Windows-MSVC.cmake
  set(CMAKE_CXX_FLAGS
    /EHsc         #*Specifies the model of exception handling (sc options).
    /DWIN32       #*Windows Platform (regardless of architecture)
    /D_WINDOWS    #*
  )
  set(CMAKE_CXX_FLAGS_RELEASE
    /O2           #*Creates fast code (Og+Oi+Ot+Oy+Ob2+GF+Gy).
    # /Ob2        #*Controls inline expansion (level 2). (part of O2)
    /DNDEBUG      #*Enables or disables compilation of assertions (CSE traditionally did not define this. See https://github.com/cse-sim/cse/issues/285)
  )
  set(CMAKE_CXX_FLAGS_DEBUG
    /Ob0          #*Controls inline expansion (level 0 -- disabled).
    /Od           #*Disables optimization.
    /Zi           #*Generates complete debugging information.
    /RTC1         #*Enables run-time error checking.
  )
  else () # GCC or Clang or AppleClang
  # See https://gitlab.kitware.com/cmake/cmake/-/blob/master/Modules/Compiler/GNU.cmake
  set(CMAKE_CXX_FLAGS "")
  set(CMAKE_CXX_FLAGS_RELEASE
    -O3           #*Maximum optimization (see https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html#Optimize-Options).
    -DNDEBUG      #*Enables or disables compilation of assertions (CSE traditionally did not define this. See https://github.com/cse-sim/cse/issues/285)
  )
  set(CMAKE_CXX_FLAGS_DEBUG
    -g            #*Produce debugging information in the operating system’s native format.
  )
  endif()

  # Convert lists to space-separated strings
  list(JOIN CMAKE_CXX_FLAGS " " CMAKE_CXX_FLAGS)
  list(JOIN CMAKE_CXX_FLAGS_RELEASE " " CMAKE_CXX_FLAGS_RELEASE)
  list(JOIN CMAKE_CXX_FLAGS_DEBUG " " CMAKE_CXX_FLAGS_DEBUG)
endmacro()