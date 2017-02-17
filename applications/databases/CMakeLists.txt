
if (FEELPP_HAS_HDF5 )
  if ( FEELPP_HAS_ACUSIM )
    feelpp_add_application( converter_acusim SRCS converter_database_acusim.cpp INCLUDE_IN_ALL PROJECT databases INSTALL)
  endif()
  feelpp_add_application( export SRCS export_database.cpp INCLUDE_IN_ALL PROJECT databases INSTALL)
  feelpp_add_application( pod SRCS pod_database.cpp INCLUDE_IN_ALL PROJECT databases INSTALL)
endif()
