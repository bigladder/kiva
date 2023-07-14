# - Config file for the MathGL package
# It defines the following variables
 
set(MathGL2_INCLUDE_DIRS "/usr/local/include")
set(MathGL2_LIBRARIES_DIRS "/usr/local/lib")
set(MathGL2_HAVE_QT5 "")
set(MathGL2_HAVE_QT4 "")
set(MathGL2_HAVE_WX "OFF")
set(MathGL2_HAVE_FLTK "OFF")
set(MathGL2_HAVE_GLUT "OFF")
set(MathGL2_HAVE_PTHREAD "1")
set(MathGL2_HAVE_OPENMP "FALSE")

include(CMakeFindDependencyMacro)
# Adding dependency for Threads imported target
if (MathGL2_HAVE_PTHREAD STRGREATER "")
  find_dependency(Threads)
endif()

# Adding dependency for OpenMP imported target
if (MathGL2_HAVE_OPENMP STRGREATER "")
  find_dependency(OpenMP)
endif()

# Compute paths
get_filename_component(MathGL2_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Our library dependencies (contains definitions for IMPORTED targets)
if(NOT TARGET mgl AND NOT MathGL_BINARY_DIR)
  include("${MathGL2_CMAKE_DIR}/MathGLTargets.cmake")
endif()
