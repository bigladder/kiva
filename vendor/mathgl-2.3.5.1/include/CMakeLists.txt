configure_file(${MathGL_SOURCE_DIR}/include/config.h.in ${MathGL_BINARY_DIR}/include/mgl2/config.h)

install(DIRECTORY mgl2/ DESTINATION ${MGL_INCLUDE_PATH}
                        PATTERN ".svn" EXCLUDE
                        PATTERN "fltk.h" EXCLUDE
                        PATTERN "glut.h" EXCLUDE
                        PATTERN "wx.h" EXCLUDE
                        PATTERN "qt.h" EXCLUDE
                        PATTERN "qmathgl.h" EXCLUDE
                        PATTERN "window.h" EXCLUDE
                        PATTERN "*old.h" EXCLUDE
                        PATTERN "*~" EXCLUDE
                        PATTERN "*.fs"
                        PATTERN "*.pas"
                        PATTERN "*.h")
