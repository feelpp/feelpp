#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --exclusive

module purge
module load feelpp.profile
module load boost/1.72.0_gcc830_openmpi402
# module load feelpp-toolboxes/feature/fix_machine_atlas_gcc830_openmpi402

np=$1
# porder=$2
# gorder=$3
cfgfile=$2
pc=$3 # typically /data/scratch/ricka/src-testfluid/feelpp/toolbox/fluid/pmm.cfg

declare -A dict
dict[1]=a
dict[2]=b

pathtoerrors='/home/u2/ricka/feel/'$(grep 'directory=.*' $cfgfile | awk -F"=" '{print $2}')'/np_'$np'/fluid.measures.csv'
if [ -f "$pathtoerrors" ]; then
   rm $pathtoerrors
fi
pathtoinfo='/home/u2/ricka/feel/'$(grep 'directory=.*' $cfgfile | awk -F"=" '{print $2}')'/np_'$np'/fluid.informations.txt'
pathtoadoc='/home/u2/ricka/feel/'$(grep 'directory=.*' $cfgfile | awk -F"=" '{print $2}')'/np_'$np'/fluid.informations.adoc'
pathtojournal='/home/u2/ricka/feel/'$(grep 'directory=.*' $cfgfile | awk -F"=" '{print $2}')'/np_'$np'/journal.json'
pathtotime='/home/u2/ricka/feel/'$(grep 'directory=.*' $cfgfile | awk -F"=" '{print $2}')'/fluid.scalibility.FluidMechanicsSolve.data'

caseoutfile='/data/scratch/ricka/results/'$4'.csv'
echo 'h solvetime nelements dofu dofp porder gorder perr uerrh1 uerrl2' > $caseoutfile

for porder in 1 2; do
    for gorder in 1 2; do
        outfile='/data/scratch/ricka/results/'$4${dict[$gorder]}${dict[$porder]}'.csv'
        echo 'h solvetime nelements dofu dofp porder gorder perr uerrh1 uerrl2' > $outfile

        for h in 0.04 0.02 0.01; do
            # mpiexec -n $np -bind-to core feelpp_toolbox_fluid --config-files $cfgfile $pc --fluid.gmsh.straighten=0 --fluid.gmsh.hsize=$h --case.discretization='P'$(($porder+1))'P'$porder'G'$gorder --fluid.scalability-save=1
            # mpiexec -n $np -bind-to core feelpp_toolbox_fluid --config-file $cfgfile --fluid.gmsh.straighten=0 --fluid.gmsh.hsize=$h --case.discretization='P'$(($porder+1))'P'$porder'G'$gorder --fluid.scalability-save=1
            # mpiexec -n $np /data/scratch/ricka/src-testfluid/build/toolboxes/fluid/feelpp_toolbox_fluid --config-files $cfgfile $pc --fluid.gmsh.straighten=0 --fluid.gmsh.hsize=$h --case.discretization='P'$(($porder+1))'P'$porder'G'$gorder --fluid.scalability-save=1
            mpiexec -n $np /data/scratch/ricka/src-testfluid/build/toolboxes/fluid/feelpp_toolbox_fluid --config-file $cfgfile --fluid.gmsh.straighten=0 --fluid.gmsh.hsize=$h --case.discretization='P'$(($porder+1))'P'$porder'G'$gorder --fluid.scalability-save=1
            raw=$(grep 'nb dof (u)      *' $pathtoinfo | awk -F":" '{print $2}')
            dofu=$(echo $raw| cut -d'(' -f 1)
            raw=$(grep 'nb dof (p)      *' $pathtoinfo | awk -F":" '{print $2}')
            dofp=$(echo $raw| cut -d'(' -f 1)
            raw=$(grep -A1 'n_elements' $pathtoadoc | grep -v "n_elements" | grep -zoP '[0-9]+')
            nel=$(echo $raw)
            raw=$(grep -o '[0-9].[0-9]*[e][+-][0-9]*' $pathtotime | tail -1)
            tsolve=$(echo $raw)
            output=$h' '$tsolve' '$nel' '$dofu' '$dofp' '$porder' '$gorder' '$(grep -o '[0-9].[0-9]*e-[0-9]*' $pathtoerrors | tr "\n" " ")
            echo ${output::-1} >> $outfile
            echo ${output::-1} >> $caseoutfile
        done
    done
done

cat $caseoutfile
