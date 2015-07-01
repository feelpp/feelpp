#!/bin/bash

# This script is to test the precAFP behavior

# $1 - nbProcs
# $2 - Dim
# $3 - mu
# $4 - hsize
function simu {
      echo "mpirun -np $1 ./feelpp_test_precAFD${2}D                "           
      echo "  --config-files precAFC${2}D.cfg backend.cfg           "           
      echo "  --mu ${3}                                             "   
      echo "  --gmsh.hsize ${4}                                     "      
      echo "  --ms.backend.pc-type=$5                               "      
      echo "  --ms.backend.ksp-type=$6                              "      
      echo "  --ms.blockms.11.backend.pc-type=$7                    "        
      echo "  --ms.blockms.11.backend.ksp-type=$8                   "   
      echo "  --ms.blockms.22.backend.pc-type=$9                    "   
      echo "  --ms.blockms.22.backend.ksp-type=${10}                  "    
      echo "  --title ${2}D-h-${4}-mu-${3}-$5-${6}_$7-${8}_${9}-${10} | grep RES >> ${2}_${11}"
}

NPROCS=6
OUTFILE=res.txt

for D in `seq 2 3`;
do
  for mu in `perl -le'for my $i (0..7) { print 10**-$i }'`;
  do
    for h in `perl -le'for my $i (1..7) { print 1/(2**$i) }'`; 
    do
      simu $NPROCS $D $mu $h gmres lu gmres lu gmres lu $OUTFILE
      simu $NPROCS $D $mu $h gmres blockms gmres lu gmres lu $OUTFILE
      simu $NPROCS $D $mu $h gmres blockms gmres gamg gmres gamg $OUTFILE
    done
  done
done
