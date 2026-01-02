include(GenerateExportHeader)
include(../../scripts/qt.cmake)

if(enable-qt5 OR (enable-qt AND NOT TARGET mgl-qt))
	set(MGL_QT_FILES ../qt.cpp ../../include/mgl2/qt.h ../../include/mgl2/qmathgl.h)

        if(enable-qt5)
                find_qt_libs(5 ${mgl_qt5_version} ON ON ${mgl_qt5_components})
        else(enable-qt5)
                find_qt_libs(5 ${mgl_qt5_version} OFF ON ${mgl_qt5_components})
        endif(enable-qt5)

        if(MGL_HAVE_QT)
                if(enable-qt5asqt OR enable-qt)
                        mgl_add_qt_lib(5 ON)
                else(enable-qt5asqt OR enable-qt)
                        mgl_add_qt_lib(5 OFF)
                endif(enable-qt5asqt OR enable-qt)
        endif(MGL_HAVE_QT)

endif(enable-qt5 OR (enable-qt AND NOT TARGET mgl-qt))
