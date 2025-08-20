# --- Ruby 3 interpreter only (no headers/libs) -------------------------------
# Try specific ruby3 names first, then generic "ruby"
find_program(Ruby_EXECUTABLE NAMES ruby3.3 ruby3.2 ruby3.1 ruby3 ruby)
if (NOT Ruby_EXECUTABLE)
  message(FATAL_ERROR "Ruby 3 interpreter not found (tried ruby3*, ruby).")
endif()

# Query the actual version from the interpreter
execute_process(
  COMMAND "${Ruby_EXECUTABLE}" -e "print RUBY_VERSION"
  OUTPUT_VARIABLE Ruby_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "[feelpp] Ruby interpreter: ${Ruby_EXECUTABLE} (version ${Ruby_VERSION})")

if (NOT Ruby_VERSION MATCHES "^3\\.")
  message(FATAL_ERROR "Found Ruby ${Ruby_VERSION}, but Ruby 3.x is required.")
endif()

# --- Bundler ---------------------------------------------------------------
find_program(BUNDLE_EXECUTABLE NAMES bundle bundle3 bundle3.0)
if (NOT BUNDLE_EXECUTABLE)
  message(FATAL_ERROR "Bundler not found. Install it (e.g. `gem install bundler`).")
endif()
message(STATUS "[feelpp] Bundler: ${BUNDLE_EXECUTABLE}")

if ( Ruby_EXECUTABLE AND BUNDLE_EXECUTABLE)
  # --- Generate and build the manuals ----------------------------------------
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc.rb.in
                 ${CMAKE_CURRENT_BINARY_DIR}/doc.rb @ONLY)

  add_custom_target(mantoolboxes
    COMMAND "${BUNDLE_EXECUTABLE}" config set --local path vendor/bundle
    COMMAND "${BUNDLE_EXECUTABLE}" install
    COMMAND "${Ruby_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/doc.rb"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/doc.rb.in"
    COMMENT "Generating manual for toolbox applications"
    VERBATIM)
endif()