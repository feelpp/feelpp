IF ( ENABLE_INSTANTIATION_MODE )
  SET( LIFEFILTERS_INST_SRCS
    gmshtensorizeddomain.cpp gmshsimplexdomain.cpp gmshellipsoiddomain.cpp
   # exporterensight.cpp exportergmsh.cpp # exportergnuplot.cpp
    )
  endif()

add_library( lifefilters SHARED
  #exporterfactorysimplex.cpp exporterfactorysimplexproduct.cpp
  exporter.cpp gmsh.cpp importergambit.cpp ${LIFEFILTERS_INST_SRCS}
  )
set_target_properties( lifefilters PROPERTIES VERSION ${LIFE_SHARED_VERSION} SOVERSION ${LIFE_SHARED_SOVERSION})
target_link_libraries (lifefilters lifediscr)
