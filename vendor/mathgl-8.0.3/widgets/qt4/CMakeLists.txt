include(GenerateExportHeader)
include(../../scripts/qt.cmake)

if(enable-qt4 OR (enable-qt AND NOT TARGET mgl-qt))

	set(MGL_QT_FILES ../qt.cpp ../../include/mgl2/qt.h ../../include/mgl2/qmathgl.h)

        if(enable-qt4)
                find_qt_libs(4 ${mgl_qt4_version} ON ON ${mgl_qt4_components})
        else(enable-qt4)
                find_qt_libs(4 ${mgl_qt4_version} OFF ON ${mgl_qt4_components})
        endif(enable-qt4)

        if(MGL_HAVE_QT)
                if(enable-qt4asqt OR enable-qt)
                        mgl_add_qt_lib(4 ON)
                else(enable-qt4asqt OR enable-qt)
                        mgl_add_qt_lib(4 OFF)
                endif(enable-qt4asqt OR enable-qt)
        endif(MGL_HAVE_QT)

endif(enable-qt4 OR (enable-qt AND NOT TARGET mgl-qt))
