find_package(Ruby 3.0)
find_program(BUNDLE_EXECUTABLE bundle NAMES bundle3.0)
message(STATUS "[feelpp] ruby is ${Ruby_FOUND}: ${Ruby_EXECUTABLE}, ${BUNDLE_EXECUTABLE}" )
if ( Ruby_EXECUTABLE )

    configure_file(doc.rb.in ${CMAKE_CURRENT_SOURCE_DIR}/doc.rb @ONLY)
    add_custom_target(mantoolboxes
        COMMAND ${BUNDLE_EXECUTABLE} config set config set --local path 'vendor/bundle'
        COMMAND ${BUNDLE_EXECUTABLE} install
        COMMAND ${BUNDLE_EXECUTABLE} exec ${RUBY} doc.rb
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDS doc.rb.in
        COMMENT "Generating manual for toolbox applications"
        VERBATIM)
else()
    
endif()
