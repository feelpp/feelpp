IF ( ENABLE_INSTANTIATION_MODE )
  SET( LIFEFILTERS_INST_SRCS 
    gmshtensorizeddomain.cpp gmshsimplexdomain.cpp 
    exporterensight.cpp exportergmsh.cpp exportergnuplot.cpp
    )
  endif()

add_library( lifefilters
  exporter.cpp gmsh.cpp importergambit.cpp ${LIFEFILTERS_INST_SRCS}
  )

target_link_libraries (lifefilters lifediscr)
