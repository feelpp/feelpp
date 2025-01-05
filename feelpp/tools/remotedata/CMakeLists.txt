message( STATUS "--------------------------------------------")
message( STATUS "[tools/remotedata] Enable feelpp_remotedata")
feelpp_add_application( remotedata SRCS remotedata.cpp INSTALL MAN remotedata )

