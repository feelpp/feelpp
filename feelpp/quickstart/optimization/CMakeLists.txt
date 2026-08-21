feelpp_add_application(
  optimization_field
  SRCS qs_optimization_field.cpp
  EXEC QS_OPTIMIZATION_FIELD
  INSTALL
  TESTS )

add_dependencies(quickstart ${QS_OPTIMIZATION_FIELD})
