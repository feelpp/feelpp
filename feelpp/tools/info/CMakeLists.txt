message( STATUS "--------------------------------------------")
message( STATUS "[tools/info] Enable feelpp_info")
feelpp_add_application(info SRCS info.cpp INSTALL MAN info)
