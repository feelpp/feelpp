# Actually include the module
include(metabench)

# Add new datasets
metabench_add_dataset(bench_poly_dataset1 "lagrange2d.cpp.erb" "[1,2,3,4]")
target_link_libraries(bench_poly_dataset1 Feelpp::feelpp)
metabench_add_dataset(bench_poly_dataset2 "lagrange3d.cpp.erb" "[1,2,3,4]")
target_link_libraries(bench_poly_dataset2 Feelpp::feelpp)

# Add a new chart
metabench_add_chart(poly_chart DATASETS bench_poly_dataset1 bench_poly_dataset2)

add_custom_target(bench_comp_poly)
#feelpp_add_application( compile_element SRCS  element.cpp EXEC B_C_ELEMENT)
#add_dependencies(bench_comp_poly ${B_C_ELEMENT})
