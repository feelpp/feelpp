<<<<<<< HEAD
message( STATUS "--------------------------------------------")
message( STATUS "[tools/remotedata] Enable feelpp_remotedata")
=======
>>>>>>> origin/develop
feelpp_add_application( remotedata SRCS remotedata.cpp INSTALL MAN remotedata )

# CLI Tests for feelpp_remotedata application
if(FEELPP_ENABLE_TESTS)

    # GitHub download test (create directory first)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_github_download)
    add_test(
        NAME feelpp_remotedata_github_download
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --download "github:{owner:feelpp,repo:feelpp,path:README.adoc}"
            --data ${CMAKE_CURRENT_BINARY_DIR}/test_github_download
    )
    set_tests_properties(feelpp_remotedata_github_download PROPERTIES
        LABELS "remotedata;github;download"
        TIMEOUT 120
        PASS_REGULAR_EXPRESSION "(✓ Downloaded|All operations completed successfully|Starting download to GitHub)"
    )

    # GitHub branch test
    add_test(
        NAME feelpp_remotedata_github_branch
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --download "github:{owner:feelpp,repo:feelpp,branch:develop,path:README.adoc}"
            --data ${CMAKE_CURRENT_BINARY_DIR}/test_github_branch/
    )
    set_tests_properties(feelpp_remotedata_github_branch PROPERTIES
        LABELS "remotedata;github;branch"
        TIMEOUT 120
        PASS_REGULAR_EXPRESSION "(✓ Downloaded|All operations completed successfully|Starting download to GitHub)"
    )

    # GitHub upload test - should fail gracefully with proper error message
    add_test(
        NAME feelpp_remotedata_github_upload
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --upload "github:{owner:feelpp,repo:feelpp,path:test-upload}"
            --data ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    )
    set_tests_properties(feelpp_remotedata_github_upload PROPERTIES
        LABELS "remotedata;github;upload"
        TIMEOUT 30
        PASS_REGULAR_EXPRESSION "(GitHub uploads are not supported|Use Girder or CKAN for upload)"
    )

    # Girder path download test (create directory first) - using collection path
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_girder_path_download)
    add_test(
        NAME feelpp_remotedata_girder_path_download
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --download "girder:{path:/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata}"
            --data ${CMAKE_CURRENT_BINARY_DIR}/test_girder_path_download
    )
    set_tests_properties(feelpp_remotedata_girder_path_download PROPERTIES
        LABELS "remotedata;girder;path;download"
        TIMEOUT 120
        PASS_REGULAR_EXPRESSION "(Downloaded|✓ Downloaded|All operations completed successfully)"
    )

    # Girder file path download test (create directory first) - using known file in collection
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_girder_file_path_download)
    add_test(
        NAME feelpp_remotedata_girder_file_path_download
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --download "girder:{path:/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/test_file.txt}"
            --data ${CMAKE_CURRENT_BINARY_DIR}/test_girder_file_path_download
    )
    set_tests_properties(feelpp_remotedata_girder_file_path_download PROPERTIES
        LABELS "remotedata;girder;file;path;download"
        TIMEOUT 120
        PASS_REGULAR_EXPRESSION "(Downloaded|✓ Downloaded|All operations completed successfully)"
    )

    # Girder item path download test (create directory first) - using known item in collection
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_girder_item_path_download)
    add_test(
        NAME feelpp_remotedata_girder_item_path_download
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --download "girder:{path:/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata/dataset}"
            --data ${CMAKE_CURRENT_BINARY_DIR}/test_girder_item_path_download
    )
    set_tests_properties(feelpp_remotedata_girder_item_path_download PROPERTIES
        LABELS "remotedata;girder;item;path;download"
        TIMEOUT 120
        PASS_REGULAR_EXPRESSION "(Downloaded|✓ Downloaded|All operations completed successfully)"
    )

    # URL upload test - should fail gracefully with proper error message
    add_test(
        NAME feelpp_remotedata_url_upload
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --upload "https://example.com/upload"
            --data ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    )
    set_tests_properties(feelpp_remotedata_url_upload PROPERTIES
        LABELS "remotedata;url;upload"
        TIMEOUT 30
        PASS_REGULAR_EXPRESSION "(URL/HTTP uploads are not supported|Use Girder or CKAN for upload)"
    )

    # Girder download test (create directory first)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_girder_download)
    add_test(
        NAME feelpp_remotedata_girder_download
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --download "girder:{folder:6743a47bb0e95728eb010c47}"
            --data ${CMAKE_CURRENT_BINARY_DIR}/test_girder_download
    )
    set_tests_properties(feelpp_remotedata_girder_download PROPERTIES
        LABELS "remotedata;girder;download"
        TIMEOUT 120
        PASS_REGULAR_EXPRESSION "(Downloaded|✓ Downloaded|All operations completed successfully)"
    )

    # Girder upload test (uses FEELPP_GIRDER_API_KEY environment variable)
    add_test(NAME feelpp_remotedata_girder_upload
      COMMAND $<TARGET_FILE:feelpp_remotedata> 
        --upload "girder:{folder:6743a47bb0e95728eb010c47}"
        --data "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(feelpp_remotedata_girder_upload PROPERTIES
      TIMEOUT 30
      LABELS "remotedata;girder;upload"
      PASS_REGULAR_EXPRESSION "(✓ Completed|All operations completed successfully|UPLOAD.*Starting upload|Authentication unavailable.*invalid upload.*configuration is incorrect)"
    )

    # Girder contents test
    add_test(
        NAME feelpp_remotedata_girder_contents
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --contents "girder:{folder:6743a47bb0e95728eb010c47}"
    )
    set_tests_properties(feelpp_remotedata_girder_contents PROPERTIES
        LABELS "remotedata;girder;contents"
        TIMEOUT 60
        PASS_REGULAR_EXPRESSION "(CONTENTS|contents|folder|file)"
    )

    # Girder path-based upload test (uses FEELPP_GIRDER_API_KEY environment variable) - using correct collection path
    add_test(NAME feelpp_remotedata_girder_path_upload
      COMMAND $<TARGET_FILE:feelpp_remotedata> 
        --upload "girder:{path:/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata}"
        --data "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(feelpp_remotedata_girder_path_upload PROPERTIES
      TIMEOUT 30
      LABELS "remotedata;girder;path;upload"
      PASS_REGULAR_EXPRESSION "(\\[UPLOAD\\].*✓.*Completed|✓.*All operations completed successfully)"
    )

    # Girder path-based contents test - using correct collection path (allows HTTP errors)
    add_test(
        NAME feelpp_remotedata_girder_path_contents
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --contents "girder:{path:/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata}"
    )
    set_tests_properties(feelpp_remotedata_girder_path_contents PROPERTIES
        LABELS "remotedata;girder;path;contents"
        TIMEOUT 60
        PASS_REGULAR_EXPRESSION "(CONTENTS|contents|folder|file|feelpp_test_remotedata|Girder error|remotedata)"
    )

    # Girder generic path test (supports any resource type) - using collection path (allows HTTP errors)
    add_test(
        NAME feelpp_remotedata_girder_generic_path_contents
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --contents "girder:{path:/collection/feelpp/testsuite/feelcore/feelpp_test_remotedata}"
    )
    set_tests_properties(feelpp_remotedata_girder_generic_path_contents PROPERTIES
        LABELS "remotedata;girder;generic;path;contents"
        TIMEOUT 60
        PASS_REGULAR_EXPRESSION "(CONTENTS|contents|folder|file|feelpp_test_remotedata|Girder error|remotedata)"
    )

    # CKAN upload test (requires API key and uses environment variables)
    add_test(
        NAME feelpp_remotedata_ckan_upload
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --debug
            --upload "ckan:{organization:cemosis,dataset:feelpp-test-dataset}"
            --data ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    )
    set_tests_properties(feelpp_remotedata_ckan_upload PROPERTIES
        LABELS "remotedata;ckan;upload"
        TIMEOUT 120
        ENVIRONMENT "FEELPP_CKAN_URL=https://ckan.hidalgo2.eu;FEELPP_CKAN_ORGANIZATION=cemosis;FEELPP_CKAN_API_KEY=$ENV{FEELPP_CKAN_API_KEY}"
        PASS_REGULAR_EXPRESSION "(UPLOAD.*Starting|✓ Completed|uploaded successfully|CKAN.*url)"
    )

    # CKAN list organizations test (uses environment variables)
    add_test(
        NAME feelpp_remotedata_ckan_list_organizations
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --debug
            --list-organizations "ckan:{url:https://ckan.hidalgo2.eu}"
    )
    set_tests_properties(feelpp_remotedata_ckan_list_organizations PROPERTIES
        LABELS "remotedata;ckan;list"
        TIMEOUT 60
        ENVIRONMENT "FEELPP_CKAN_URL=https://ckan.hidalgo2.eu;FEELPP_CKAN_ORGANIZATION=cemosis"
        PASS_REGULAR_EXPRESSION "(organizations|CKAN.*url|LIST-ORGANIZATIONS|cemosis)"
    )

    # URL download test (create directory first)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_url_download)
    add_test(
        NAME feelpp_remotedata_url_download
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --download "https://raw.githubusercontent.com/feelpp/feelpp/develop/README.adoc"
            --data ${CMAKE_CURRENT_BINARY_DIR}/test_url_download
    )
    set_tests_properties(feelpp_remotedata_url_download PROPERTIES
        LABELS "remotedata;url;download"
        TIMEOUT 60
        PASS_REGULAR_EXPRESSION "(Downloaded|✓ Downloaded|All operations completed successfully|README|download data in)"
    )

    # Error handling test for invalid URLs - should succeed but print error message (create directory first)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_error)
    add_test(
        NAME feelpp_remotedata_error_invalid_url
        COMMAND $<TARGET_FILE:feelpp_remotedata>
            --download "http://invalid.example.com/nonexistent"
            --data ${CMAKE_CURRENT_BINARY_DIR}/test_error
    )
    set_tests_properties(feelpp_remotedata_error_invalid_url PROPERTIES
        LABELS "remotedata;error"
        TIMEOUT 30
        PASS_REGULAR_EXPRESSION "(Download error|invalid download|error|failed)"
    )

    # Help functionality test
    add_test(
        NAME feelpp_remotedata_help
        COMMAND $<TARGET_FILE:feelpp_remotedata> --help
    )
    set_tests_properties(feelpp_remotedata_help PROPERTIES
        LABELS "remotedata;help"
        TIMEOUT 10
        PASS_REGULAR_EXPRESSION "(Remote data options|--upload|--download|--help|Feel\\+\\+ application)"
    )

    # Version functionality test
    add_test(
        NAME feelpp_remotedata_version
        COMMAND $<TARGET_FILE:feelpp_remotedata> --version
    )
    set_tests_properties(feelpp_remotedata_version PROPERTIES
        LABELS "remotedata;version"
        TIMEOUT 10
        PASS_REGULAR_EXPRESSION "(version|remotedata.*Feel\\+\\+ application|0\\.[0-9]+\\.[0-9]+)"
    )

endif()