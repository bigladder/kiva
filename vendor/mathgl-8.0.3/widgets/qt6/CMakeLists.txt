include(GenerateExportHeader)
include(../../scripts/qt.cmake)

if(enable-qt6 OR (enable-qt AND NOT TARGET mgl-qt))

	set(MGL_QT_FILES ../qt.cpp ../../include/mgl2/qt.h ../../include/mgl2/qmathgl.h)

        if(enable-qt6)
                find_qt_libs(6 ${mgl_qt6_version} ON ON ${mgl_qt6_components})
        else(enable-qt6)
                find_qt_libs(6 ${mgl_qt6_version} OFF ON ${mgl_qt6_components})
        endif(enable-qt6)

        if(MGL_HAVE_QT)
                if(enable-qt6asqt OR enable-qt)
                        mgl_add_qt_lib(6 ON)
                else(enable-qt6asqt OR enable-qt)
                        mgl_add_qt_lib(6 OFF)
                endif(enable-qt6asqt OR enable-qt)
        endif(MGL_HAVE_QT)

endif(enable-qt6 OR (enable-qt AND NOT TARGET mgl-qt))
