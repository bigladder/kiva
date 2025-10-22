add_executable(make_pas make_pas.cpp)

if(MSVC)
set(link_type -static)
else(MSVC)
set(link_type)
endif(MSVC)

add_executable(mgltask mgltask.cpp)
target_link_libraries(mgltask mgl${link_type} ${getopt_lib-static})
install(
	TARGETS mgltask
	EXPORT MathGLTargets
	RUNTIME DESTINATION ${MathGL_INSTALL_BIN_DIR}
)

add_executable(mglconv mglconv.cpp)
target_link_libraries(mglconv mgl${link_type} ${getopt_lib-static})
install(
	TARGETS mglconv
	EXPORT MathGLTargets
	RUNTIME DESTINATION ${MathGL_INSTALL_BIN_DIR}
)

add_executable(mgl.cgi mglcgi.cpp)
target_link_libraries(mgl.cgi mgl-static)
install(
	TARGETS mgl.cgi
	EXPORT MathGLTargets
# should be /usr/lib/cgi-bin/
	RUNTIME DESTINATION ${MGL_CGI_PATH}
)

mgl_po_src(mglconv.cpp mglview.cpp mglcgi.cpp mgltask.cpp)

if(MGL_HAVE_FLTK)
	add_definitions(-DUSE_FLTK)
	add_executable(mglview mglview.cpp)
	target_link_libraries(mglview mgl-fltk ${getopt_lib-static} ${FLTK_LIBRARIES})
	install(
		TARGETS mglview
		EXPORT MathGLTargets
		RUNTIME DESTINATION ${MathGL_INSTALL_BIN_DIR}
	)

else(MGL_HAVE_FLTK)
	if(TARGET mgl-qt)
		include(../scripts/qt.cmake)
                find_qt_libs(${mgl_qt_major} ${mgl_qt_version} ON ON ${mgl_qt_components})

		add_executable(mglview mglview.cpp)
		target_link_libraries(mglview mgl-qt${link_type} ${getopt_lib-static} ${MGL_QT_LIBS})

		install(
			TARGETS mglview
			EXPORT MathGLTargets
			RUNTIME DESTINATION ${MathGL_INSTALL_BIN_DIR}
		)
	endif(TARGET mgl-qt)
endif(MGL_HAVE_FLTK)
