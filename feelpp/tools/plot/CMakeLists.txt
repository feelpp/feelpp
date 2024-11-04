message( STATUS "[feelpp_plot] Configuring plot: FEELPP_HAS_PLOT:${FEELPP_HAS_MATPLOT}" )
if ( FEELPP_HAS_MATPLOT )
    message( STATUS "[feelpp_plot] Matplot++ found" )
    feelpp_add_application( plot SRCS plot.cpp  INSTALL )
#    target_link_libraries( feelpp_plot PRIVATE Matplot++::matplot )

endif()