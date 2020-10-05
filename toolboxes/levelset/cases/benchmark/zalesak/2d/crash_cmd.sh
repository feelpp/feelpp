BIN_DIR="/ssd/derhovsepian/feelpp/Debug_install_develop_rpath_custom_ls_toolbox/bin/feelpp_toolbox_levelset_2d"
BIN_DIR="/ssd/derhovsepian/feelpp/Release_build_develop_rpath_custom_ls_toolbox/toolboxes/levelset/feelpp_toolbox_levelset_2d"
CFG_DIR="/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/zalesak/2d"
echo ${BIN_DIR}
mpirun -np 1 ${BIN_DIR} \
       --config-file=${CFG_DIR}/slotteddisk.cfg \
       --levelset.mesh.filename=${CFG_DIR}/domain.geo \
       --directory=toolboxes/levelset/zalesak2d/test_exportsspace_convergence/GALS_reinitevery_-1/h_0.32/quadorder_2 \
       --ts.time-step=36.57 \
       --ts.time-final=1.0 \
       --levelset.reinit-every=-1 \
       --levelset.ts.order=2 \
       --levelset.gmsh.hsize=0.32 \
       --levelset.stabilization.method=GALS \
       --levelset.quad.order=2 
