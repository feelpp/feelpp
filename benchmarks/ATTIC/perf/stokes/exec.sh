#!/bin/.sh

for h in 0.2 0.1 0.05 0.025 0.0125 0.01
do

  ./feel_perf_stokes_p2p1 --hsize=$h --bctype=1 -ksp_constant_null_space true -pc_type fieldsplit -pc_fieldsplit_type schur -fieldsplit_0_pc_type lu -fieldsplit_0_ksp_type gmres -fieldsplit_1_pc_type none -fieldsplit_1_ksp_type gmres   -fieldsplit_1_ksp_max_it 1000 -fieldsplit_1_ksp_converged_use_initial_residual_norm true > log-$h 2>&1

echo "=============================="
echo "h=$h"
cat $HOME/feel/perf/stokes/stokes/Simplex_2_1_2/P2/h_$h/toto-1.0 | grep "_error||_2"
cat $HOME/feel/perf/stokes/stokes/Simplex_2_1_2/P2/h_$h/toto-1.0 | grep time
cat $HOME/feel/perf/stokes/stokes/Simplex_2_1_2/P2/h_$h/toto-1.0 | grep elements
cat $HOME/feel/perf/stokes/stokes/Simplex_2_1_2/P2/h_$h/toto-1.0 | grep dof





done
