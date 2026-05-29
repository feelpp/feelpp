set_directory_properties(PROPERTIES LABELS "toolboxes-tests;cfpdes")

feelpp_add_test( coefficientformpdes LINK_LIBRARIES Feelpp::feelpp_toolbox_coefficientformpdes_lib )
