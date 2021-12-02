#!/bin/bash

#SBATCH --nodes=16
#SBATCH --time=04:00:00
#SBATCH --exclusive

module load feelpp-toolboxes/develop_gcc830_openmpi402 

np=16

outfile='/data/scratch/ricka/'$2'.csv'
echo '|h|perr|uerrh1|uerrl2|' > $outfile
echo '|-|-|-|-|' >> $outfile

pathtoerrors='/home/u2/ricka/feel/'$(grep 'directory=.*' $1 | awk -F"=" '{print $2}')'/np_'$np'/fluid.measures.csv'
rm $pathtoerrors

for h in 0.08 0.04 0.02 0.01 ; do
    mpirun -np $np feelpp_toolbox_fluid --config-file $1 --fluid.gmsh.hsize=$h #/data/scratch/ricka/src-testfluid/build/toolboxes/fluid/
    output='|'$h'|'$(grep -o '[0-9].[0-9]*e-[0-9]*' $pathtoerrors | tr "\n" "|")
    echo $output >> $outfile
done

cat $outfile