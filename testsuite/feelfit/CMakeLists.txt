set_directory_properties(PROPERTIES LABEL testfit )

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )


feelpp_add_test( nofeel_interpolator CFG data.txt)
target_link_libraries( feelpp_test_nofeel_interpolator gsl )

feelpp_add_test( interpolator CFG data.txt test_fit.cfg )

