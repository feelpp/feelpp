# Actually include the module
include(metabench)

# Add new datasets
metabench_add_dataset(bench_discr_pch "pch.cpp.erb" "[1,2,3]")
target_link_libraries(bench_discr_pch Feelpp::feelpp)
metabench_add_dataset(bench_discr_pch2 "pch2.cpp.erb" "[1,2,3]")
target_link_libraries(bench_discr_pch2 Feelpp::feelpp)
metabench_add_dataset(bench_discr_pch3 "pch3.cpp.erb" "[1,2,3]")
target_link_libraries(bench_discr_pch3 Feelpp::feelpp)
metabench_add_dataset(bench_discr_pch4 "pch4.cpp.erb" "[1,2,3]")
target_link_libraries(bench_discr_pch4 Feelpp::feelpp)
#metabench_add_dataset(dataset3 "path/to/dataset3.cpp.erb" "(1...20).step(5)")

# Add a new chart
metabench_add_chart(bench_discr_chart DATASETS bench_discr_pch bench_discr_pch2 bench_discr_pch3 bench_discr_pch4)

add_custom_target(bench_comp_discr)
#feelpp_add_application( compile_element SRCS  element.cpp EXEC B_C_ELEMENT)
#add_dependencies(bench_comp_mesh ${B_C_ELEMENT})
