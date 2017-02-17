
if (FEELPP_HAS_HDF5 AND FEELPP_HAS_ACUSIM )
  feelpp_add_application( converter_acusim SRCS converter_database_acusim.cpp INCLUDE_IN_ALL PROJECT databases INSTALL)
endif ()
