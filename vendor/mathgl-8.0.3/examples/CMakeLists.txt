
#file(COPY ${CMAKE_SOURCE_DIR}/udav_new.png DESTINATION ${CMAKE_BINARY_DIR}/examples)
file(COPY ${CMAKE_SOURCE_DIR}/examples/iris.dat DESTINATION ${CMAKE_BINARY_DIR}/examples)
file(COPY ${CMAKE_SOURCE_DIR}/examples/Equirectangular-projection.jpg DESTINATION ${CMAKE_BINARY_DIR}/examples)
file(COPY ${CMAKE_SOURCE_DIR}/examples/samples.cpp DESTINATION ${CMAKE_BINARY_DIR}/examples)
file(COPY ${CMAKE_SOURCE_DIR}/examples/wnd_samples.cpp DESTINATION ${CMAKE_BINARY_DIR}/examples)
add_executable(mgl_example wnd_samples.cpp full_test.cpp samples.cpp)
target_link_libraries(mgl_example mgl-static ${getopt_lib-static})
if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND CMAKE_COMPILER_IS_GNUCXX AND enable-lcov)
	setup_target_for_coverage(
	    NAME mgl_coverage
	    EXECUTABLE $<TARGET_FILE:mgl_example> -mgl
	    DEPENDENCIES mgl_example
	)
endif(CMAKE_BUILD_TYPE STREQUAL "Debug" AND CMAKE_COMPILER_IS_GNUCXX AND enable-lcov)

if(MGL_HAVE_FLTK)
	include_directories(${FLTK_INCLUDE_DIR})
	add_executable(mgl_fltk_example wnd_samples.cpp fltk_example.cpp)
	target_link_libraries(mgl_fltk_example mgl-fltk)
endif(MGL_HAVE_FLTK)

if(MGL_HAVE_GLUT)
	add_executable(mgl_glut_example wnd_samples.cpp glut_example.cpp)
	target_link_libraries(mgl_glut_example mgl-glut)
endif(MGL_HAVE_GLUT)

if(MGL_HAVE_WX)
	include(${wxWidgets_USE_FILE})
	add_executable(mgl_wx_example wnd_samples.cpp wx_example.cpp)
	target_link_libraries(mgl_wx_example mgl-wx)
endif(MGL_HAVE_WX)

if(TARGET mgl-qt)
        include(../scripts/qt.cmake)
        find_qt_libs(${mgl_qt_major} ${mgl_qt_version} ON ON ${mgl_qt_components})

	add_executable(mgl_qt_example wnd_samples.cpp qt_example.cpp)
	target_link_libraries(mgl_qt_example mgl-qt)

	if(MGL_HAVE_OPENGL)
		add_executable(mgl_qgl_example wnd_samples.cpp qgl_example.cpp)
		target_link_libraries(mgl_qgl_example mgl ${MGL_QT_LIBS})
	endif(MGL_HAVE_OPENGL)
endif(TARGET mgl-qt)

if(MGL_HAVE_LTDL)
	add_library(mgl_module MODULE mgl_module.cpp)
	target_link_libraries(mgl_module mgl)	# for compatibility with win32
endif(MGL_HAVE_LTDL)

if(MGL_HAVE_MPI)
	add_executable(mgl_mpi_example mpi_test.cpp)
	target_link_libraries(mgl_mpi_example mgl-mpi mgl ${MPI_LIBRARIES})
	target_include_directories(mgl_mpi_example SYSTEM PUBLIC ${MPI_CXX_INCLUDE_PATH})
endif(MGL_HAVE_MPI)
