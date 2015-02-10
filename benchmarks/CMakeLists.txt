set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

project( feel-benchmarks )

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_custom_target(benchmarks)
# set(petsc_SOLVER_OPTIONS "-ksp_monitor -pc_type asm -pc_asm_blocks 2" )
#set(petsc_SOLVER_OPTIONS "-pc_type lu -pc_factor_mat_solver_package umfpack" )
set(petsc_SOLVER_OPTIONS "-pc_type lu" )

# When PETSc is compiled with external direct solvers (UMFPACK, SuperLU, etc.),
# they can be selected like this:
set(petsc_LU_SOLVER_OPTIONS " -pc_type lu " )
set(petsc_LU_UMFPACK_SOLVER_OPTIONS " -pc_type lu -pc_factor_mat_solver_package umfpack -ksp_monitor " )

set(petsc_JACOBI_SOLVER_OPTIONS " -pc_type jacobi" )

#
# By default, we use the following options (GMRES with ILU(6) preconditioner and
# RCMK renumbering):
set(petsc_ILU_SOLVER_OPTIONS " -pc_type ilu -pc_factor_levels 6 -pc_factor_mat_ordering rcm " )

macro(feel_add_bench examplename)

  set(targetname feel_bench_${examplename})
  set(filename ${examplename}.cpp)
  add_executable(${targetname} ${filename})
  target_link_libraries(${targetname} feelpp ${FEELPP_LIBRARIES} )
  set_property(TARGET ${targetname} PROPERTY LABELS benchmarks)
  add_dependencies( benchmarks ${targetname} )
endmacro(feel_add_bench)

foreach( DIR navierstokes darcy curl perf wrapper/python python compilation )
  add_subdirectory( ${DIR} )
endforeach()
