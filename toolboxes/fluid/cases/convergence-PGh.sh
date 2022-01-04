#!/bin/bash

#SBATCH --nodes=24
#SBATCH --time=12:00:00
#SBATCH --exclusive

#module load feelpp-toolboxes/develop_gcc830_openmpi402 
module load boost/1.72.0_gcc830_openmpi402 

np=24

outfile='/data/scratch/ricka/'$2'.csv'
echo 'h,perr,uerrh1,uerrl2' > $outfile

pathtoerrors='/home/u2/ricka/feel/'$(grep 'directory=.*' $1 | awk -F"=" '{print $2}')'/np_'$np'/fluid.measures.csv'
rm $pathtoerrors

for gorder in 1 ; do
    for porder in 1 2 ; do
        echo 'P'$(($porder+1))'P'$porder'G'$gorder >> $outfile
        for h in 0.16 0.08 0.04 0.02 ; do
            mpirun -np $np /data/scratch/ricka/src-testfluid/build/toolboxes/fluid/feelpp_toolbox_fluid --config-files $1 /data/scratch/ricka/src-testfluid/feelpp/toolboxes/fluid/cases/lsc.cfg --fluid.gmsh.hsize=$h --case.discretization='P'$(($porder+1))'P'$porder'G'$gorder
            output=$h','$(grep -o '[0-9].[0-9]*e-[0-9]*' $pathtoerrors | tr "\n" ",")
            echo $output >> $outfile
        done
        echo '' >> $outfile
    done
    echo '' >> $outfile
done

cat $outfile