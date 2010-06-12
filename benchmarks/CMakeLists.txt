project( life-benchmarks )

# set(petsc_SOLVER_OPTIONS "-ksp_monitor -pc_type asm -pc_asm_blocks 2" )
#set(petsc_SOLVER_OPTIONS "-pc_type lu -pc_factor_mat_solver_package umfpack" )
set(petsc_SOLVER_OPTIONS "-pc_type lu" )
#set(petsc_LU_SOLVER_OPTIONS "-ksp_monitor -pc_type lu -pc_factor_mat_solver_package umfpack -ksp_view" )
set(petsc_LU_SOLVER_OPTIONS " -pc_type lu" )
set(petsc_ILU_SOLVER_OPTIONS " -pc_type ilu" )

macro(life_add_bench examplename)

  set(targetname life_bench_${examplename})
  set(filename ${examplename}.cpp)
  add_executable(${targetname} ${filename})
  target_link_libraries(${targetname} ${LIFE_LIBRARIES} )
endmacro(life_add_bench)

foreach( DIR convergence navierstokes/kovasznay perf wrapper/python )
  add_subdirectory( ${DIR} )
endforeach()
