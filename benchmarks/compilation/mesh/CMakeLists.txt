# Actually include the module
include(metabench)

# Add new datasets
metabench_add_dataset(bench_mesh_geoentity "geoentity.cpp.erb" "[1,2,3]")
target_link_libraries(bench_mesh_geoentity Feelpp::feelpp)

metabench_add_dataset(bench_mesh_geond "geond.cpp.erb" "[1,2,3]")
target_link_libraries(bench_mesh_geond Feelpp::feelpp)

metabench_add_dataset(bench_mesh_geoelement "geoelement.cpp.erb" "[1,2,3]")
target_link_libraries(bench_mesh_geoelement Feelpp::feelpp)


metabench_add_dataset(bench_mesh_geomap "geomap.cpp.erb" "[1,2,3]")
target_link_libraries(bench_mesh_geomap Feelpp::feelpp)

metabench_add_dataset(bench_mesh_mesh "mesh.cpp.erb" "[1,2,3]")
target_link_libraries(bench_mesh_mesh Feelpp::feelpp) 

#metabench_add_dataset(dataset2 "path/to/dataset2.cpp.erb" "(1...15)")
#metabench_add_dataset(dataset3 "path/to/dataset3.cpp.erb" "(1...20).step(5)")

# Add a new chart
metabench_add_chart(bench_mesh_chart DATASETS bench_mesh_geoentity bench_mesh_geond bench_mesh_geoelement bench_mesh_geomap bench_mesh_mesh)

add_custom_target(bench_comp_mesh)
#feelpp_add_application( compile_element SRCS  element.cpp EXEC B_C_ELEMENT)
#add_dependencies(bench_comp_mesh ${B_C_ELEMENT})
