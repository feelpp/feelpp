#create the directory of the lib

# generate the library using the custom makefile
execute_process( COMMAND make -f ${OMWRAPPER_LIBDIR}/${OMWRAPPER_MAKEFILE} WORKING_DIRECTORY ${OMWRAPPER_LIBDIR} )
