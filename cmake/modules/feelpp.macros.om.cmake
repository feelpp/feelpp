#create the directory of the lib

# generate the library using the custom makefile
execute_process( COMMAND ${OMC_COMPILER} ${FMU_SCRIPT_NAME} WORKING_DIRECTORY ${OMWRAPPER_LIBDIR} )
