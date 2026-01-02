set(src_imp_dep
../include/mgl2/type.h
../include/mgl2/data.h
../include/mgl2/mgl.h
mgl.i
data.i
type.i
)

set(mgl_clean_files "")

SET_SOURCE_FILES_PROPERTIES(mathgl.i PROPERTIES CPLUSPLUS ON)

FIND_PACKAGE(SWIG)
if(NOT SWIG_FOUND)
	message(SEND_ERROR "Couldn't find swig needed for interfaces compiling.")
endif(NOT SWIG_FOUND)
INCLUDE(${SWIG_USE_FILE})

if(enable-python)

	find_package (Python COMPONENTS Interpreter Development)
	if(NOT Python_FOUND)
		message(SEND_ERROR "Couldn't find python interpreter or development.")
	endif(NOT Python_FOUND)

	execute_process(
		COMMAND ${Python_EXECUTABLE} -c "import numpy; print (numpy.get_include())"
		OUTPUT_VARIABLE NUMPY_INCLUDE_PATH
		RESULT_VARIABLE NUMPY_ERR
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	if(NOT NUMPY_INCLUDE_PATH)
		message(SEND_ERROR "Couldn't find numpy. ${NUMPY_ERR}")
	endif(NOT NUMPY_INCLUDE_PATH)
	include_directories(${Python_INCLUDE_DIRS} ${NUMPY_INCLUDE_PATH})
	if(NOT numpy_i_once)
		find_file(numpy_i NAMES numpy.i PATHS ${NUMPY_INCLUDE_PATH} NO_CACHE PATH_SUFFIXES numpy)
		if(numpy_i)
			message(STATUS "Using system numpy.i")
			file(COPY ${numpy_i} DESTINATION ${MathGL2_BINARY_DIR}/include/numpy)
		else(numpy_i)
			message(STATUS "Your system has no numpy.i. Using numpy.i from MathGL")
			file(COPY ${MathGL2_SOURCE_DIR}/lang/numpy.i DESTINATION ${MathGL2_BINARY_DIR}/include/numpy)
		endif(numpy_i)
		set(numpy_i_once TRUE CACHE INTERNAL "Do once")
	endif(NOT numpy_i_once)

	if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
		SWIG_ADD_MODULE(mathgl python mathgl.i)
	else(${CMAKE_VERSION} VERSION_LESS "3.8.0")
		SWIG_ADD_LIBRARY(mathgl LANGUAGE python SOURCES mathgl.i)
	endif(${CMAKE_VERSION} VERSION_LESS "3.8.0")

	SWIG_LINK_LIBRARIES(mathgl mgl ${Python_LIBRARIES})
	set_target_properties(mathgl PROPERTIES BUILD_WITH_INSTALL_RPATH ON)

	add_custom_command(OUTPUT mathgl.pyc
		COMMAND ${Python_EXECUTABLE}
		ARGS -c \"from py_compile import compile\; compile('mathgl.py') \"
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lang
		DEPENDS mathgl
	)
	add_custom_target(mgl_python_module ALL DEPENDS mathgl mathgl.pyc)

	install(FILES ${MathGL2_BINARY_DIR}/lang/mathgl.py DESTINATION ${Python_SITELIB})
	install(FILES ${MathGL2_BINARY_DIR}/lang/__pycache__/mathgl.cpython-${Python_VERSION_MAJOR}${Python_VERSION_MINOR}.pyc DESTINATION ${Python_SITELIB}/__pycache__)
	install (TARGETS mathgl LIBRARY DESTINATION ${Python_SITEARCH})
	set(mgl_clean_files ${mgl_clean_files} mathgl.py)

endif(enable-python)

if(enable-lua)
	INCLUDE(FindLua51)
	if(NOT LUA51_FOUND)
		message(SEND_ERROR "Couldn't find Lua 5.1 library.")
	endif(NOT LUA51_FOUND)

	include_directories(${LUA_INCLUDE_DIR})
	if(${CMAKE_VERSION} VERSION_LESS "3.8.0") 
		SWIG_ADD_MODULE(mgl-lua lua mathgl.i)
	else(${CMAKE_VERSION} VERSION_LESS "3.8.0")
		SWIG_ADD_LIBRARY(mgl-lua LANGUAGE lua SOURCES mathgl.i)
	endif(${CMAKE_VERSION} VERSION_LESS "3.8.0")
	SWIG_LINK_LIBRARIES(mgl-lua mgl ${LUA_LIBRARIES})
	set_target_properties(mgl-lua PROPERTIES PREFIX "" BUILD_WITH_INSTALL_RPATH ON)

	install (TARGETS mgl-lua LIBRARY DESTINATION ${MathGL_INSTALL_LIB_DIR})
endif(enable-lua)

if(enable-octave)
	find_program(oct_prog octave-config)
	if(NOT oct_prog)
		message(SEND_ERROR "Couldn't find octave-config needed for octave interfaces compiling.")
	endif(NOT oct_prog)
	find_program(oct_exec octave)
	if(NOT oct_exec)
		message(SEND_ERROR "Couldn't find octave needed for octave interfaces compiling.")
	endif(NOT oct_exec)
	find_program(oct_mk mkoctfile)
	if(NOT oct_mk)
		message(SEND_ERROR "Couldn't find mkoctfile needed for octave interfaces compiling.")
	endif(NOT oct_mk)
	find_program(oct_tar tar)
	if(NOT oct_tar)
		message(SEND_ERROR "Couldn't find tar needed for octave interfaces creation.")
	endif(NOT oct_tar)

	execute_process(COMMAND ${oct_prog} -p CANONICAL_HOST_TYPE
	OUTPUT_VARIABLE oct_host
	OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND ${oct_prog} -p API_VERSION
	OUTPUT_VARIABLE oct_api
	OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND ${oct_mk} -p INCFLAGS
	OUTPUT_VARIABLE oct_inc
	OUTPUT_STRIP_TRAILING_WHITESPACE)
	separate_arguments(oct_inc)
message(STATUS "${oct_prog} ${oct_host} ${oct_api}")
	set(oct_arch ${oct_host}-${oct_api})
	set(pkg_name mathgl)
	set(mgl_pkg_dir ${MathGL2_BINARY_DIR}/lang/${pkg_name}/inst/${oct_arch})

	add_compile_options(${oct_inc})
	if(${CMAKE_VERSION} VERSION_LESS "3.8.0") 
		SWIG_ADD_MODULE(mgl-oct octave mathgl.i)
	else(${CMAKE_VERSION} VERSION_LESS "3.8.0")
		SWIG_ADD_LIBRARY(mgl-oct LANGUAGE octave SOURCES mathgl.i)
	endif(${CMAKE_VERSION} VERSION_LESS "3.8.0")
	SWIG_LINK_LIBRARIES(mgl-oct mgl)
	set_target_properties(mgl-oct PROPERTIES OUTPUT_NAME mathgl PREFIX "" SUFFIX ".oct" BUILD_WITH_INSTALL_RPATH ON)
	add_custom_command(OUTPUT mathgl.tar.gz
		COMMAND ${CMAKE_COMMAND} -E make_directory ${mgl_pkg_dir}
		COMMAND ${CMAKE_COMMAND} -E copy ${MathGL2_SOURCE_DIR}/COPYING ${MathGL2_BINARY_DIR}/lang/${pkg_name}
		COMMAND ${CMAKE_COMMAND} -E copy ${MathGL2_SOURCE_DIR}/lang/DESCRIPTION ${MathGL2_BINARY_DIR}/lang/${pkg_name}
		COMMAND ${CMAKE_COMMAND} -E copy ${MathGL2_SOURCE_DIR}/lang/INDEX ${MathGL2_BINARY_DIR}/lang/${pkg_name}
		COMMAND ${CMAKE_COMMAND} -E copy ${MathGL2_SOURCE_DIR}/lang/PKG_ADD_template ${MathGL2_BINARY_DIR}/lang/${pkg_name}
		COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:mgl-oct> ${mgl_pkg_dir}
		COMMAND ${oct_tar} cpzf mathgl.tar.gz ${pkg_name}
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lang
		DEPENDS mgl-oct
	)
	add_custom_target(mgl_octave_module ALL DEPENDS mathgl.tar.gz)

	install(FILES ${CMAKE_BINARY_DIR}/lang/mathgl.tar.gz DESTINATION ${MGL_DATA_PATH})
	if(enable-octave-install)
		install(CODE "execute_process(
			COMMAND ${oct_exec} -q --eval \"pkg install mathgl.tar.gz\"
			WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lang )")
	endif(enable-octave-install)
	set(mgl_clean_files ${mgl_clean_files} mathgl)
endif(enable-octave)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${mgl_clean_files}")
