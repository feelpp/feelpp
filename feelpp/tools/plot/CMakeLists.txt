if ( FEELPP_ENABLE_MATPLOT )

    feelpp_add_application( plot SRCS plot.cpp  INSTALL )
    target_link_libraries( feelpp_plot matplot )

endif()