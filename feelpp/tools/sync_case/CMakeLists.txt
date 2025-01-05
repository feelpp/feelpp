message( STATUS "--------------------------------------------")
message( STATUS "[tools/sync_case] Enable feelpp_sync_case")
feelpp_add_application( sync_case SRCS sync_case.cpp  INSTALL MAN sync_case )
