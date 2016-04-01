set_directory_properties(PROPERTIES LABEL testfit )

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

feelpp_add_test( interpolator CFG data.txt data_interp.txt )
target_link_libraries( feelpp_test_interpolator gsl )

