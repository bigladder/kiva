# src/
# Get source files
FILE(GLOB temp_Files array/*.h array/*.c)
SOURCE_GROUP( lis\\src\\array FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files})
FILE(GLOB temp_Files esolver/*.h esolver/*.c)
SOURCE_GROUP( lis\\src\\esolver FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files})
IF (LIS_ENABLE_SAAMG OR LIS_ENABLE_FORTRAN)
	FILE(GLOB temp_Files fortran/*.h fortran/*.c fortran/*.F)
	SOURCE_GROUP( lis\\src\\fortran_c FILES ${temp_Files} )
	LIST(APPEND lis_src_Files ${temp_Files})
ENDIF()
IF (LIS_ENABLE_SAAMG)
	IF (LIS_ENABLE_MPI)
		SET(temp_Files
			fortran/amg/lis_m_aggregate_mod.F90
      fortran/amg/lis_m_count_time_mod.F90
      fortran/amg/lis_m_data_creation_AMGCG.F90
      fortran/amg/lis_m_data_structure_for_AMG.F90
      fortran/amg/lis_m_finit.F90
      fortran/amg/lis_m_hash_mod.F90
			fortran/amg/lis_m_isort.F90
			fortran/amg/lis_m_queue_mod.F90
      fortran/amg/lis_m_solver_AMGCG.F90
      fortran/amg/lis_m_solver_Gnumbering.F90
			fortran/amg/lis_m_solver_SR2.F90
		)
	ELSE()
		SET(temp_Files
      fortran/amg/lis_s_aggregate_mod.F90
      fortran/amg/lis_s_count_time_mod.F90
      fortran/amg/lis_s_data_creation_AMGCG.F90
      fortran/amg/lis_s_data_structure_for_AMG.F90
      fortran/amg/lis_s_finit.F90
      fortran/amg/lis_s_hash_mod.F90
			fortran/amg/lis_s_isort.F90
			fortran/amg/lis_s_queue_mod.F90
      fortran/amg/lis_s_solver_AMGCG.F90
      fortran/amg/lis_s_solver_Gnumbering.F90
			fortran/amg/lis_s_solver_SR2.F90
		)
	ENDIF()
	SOURCE_GROUP( lis\\src\\fortran\\amg FILES ${temp_Files} )
	LIST(APPEND lis_src_Files ${temp_Files})
ENDIF()
FILE(GLOB temp_Files matrix/*.h matrix/*.c)
SOURCE_GROUP( lis\\src\\matrix FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files})
FILE(GLOB temp_Files matvec/*.h matvec/*.c)
SOURCE_GROUP( lis\\src\\matvec FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files})
FILE(GLOB temp_Files precision/*.h precision/*.c)
SOURCE_GROUP( lis\\src\\precision FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files})
FILE(GLOB temp_Files precon/*.h precon/*.c)
SOURCE_GROUP( lis\\src\\precon FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files})
FILE(GLOB temp_Files solver/*.h solver/*.c)
SOURCE_GROUP( lis\\src\\solver FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files})
FILE(GLOB temp_Files system/*.h system/*.c)
SOURCE_GROUP( lis\\src\\system FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files})
FILE(GLOB temp_Files vector/*.h vector/*.c)
SOURCE_GROUP( lis\\src\\vector FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files})

FILE(GLOB lis_include_Files ../include/*.h)
SOURCE_GROUP( lis\\include FILES ${lis_include_Files} )


# Create the executable
INCLUDE_DIRECTORIES(
	.
	../include
)

ADD_LIBRARY( lis
        ${lis_src_Files}
)
#SET_TARGET_PROPERTIES(lis PROPERTIES LINKER_LANGUAGE C)

# Link the standard C math library
find_library(M_LIB m)
TARGET_LINK_LIBRARIES( lis ${M_LIB} )

# Link MPI
IF (LIS_ENABLE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	TARGET_LINK_LIBRARIES( lis ${MPI_C_LIBRARIES} )
ENDIF ()
