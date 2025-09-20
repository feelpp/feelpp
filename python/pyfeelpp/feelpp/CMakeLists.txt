# Add the core subdirectory which contains all the Python bindings
add_subdirectory(core)

# feelpp is now an implicit namespace package (PEP 420)
# No __init__.py file needed - this allows multiple distributions
# (pyfeelpp, pyfeelpp-mor, pyfeelpp-toolboxes) to contribute to the feelpp namespace