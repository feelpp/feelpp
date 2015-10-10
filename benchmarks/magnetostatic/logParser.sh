#!/bin/bash
# set -x
# This script is to gather results from test_prec.sh

# $1 - nbProcs
# $2 - Dim
# $3 - mu
# $4 - hsize
# $5-11 - ksp/pc config
# $12 - string - appDir
# $13 - test type : lin or sin
function res() {
  # List of files we are looking in
  title=${2}D-h-*-mu-${3}-$5-${6}_$7-${8}_${9}-${10}-${13}_${11}
  # Saved file
  out=OUT_${2}D-mu-${3}-$5-${6}_$7-${8}_${9}-${10}-${13}_${11}
  echo -e "RES\thSize\tnDof\tnNz\tmu\terr\terr_rel" > $out
  grep -h RES $title 1>> $out 2>&-
}

function nbIter() {
  # List of files we are looking in
  title=${2}D-h-*-mu-${3}-$5-${6}_$7-${8}_${9}-${10}-${13}_${11}
  # Saved file
  out=ITER_${2}D-mu-${3}-$5-${6}_$7-${8}_${9}-${10}-${13}_${11}
  echo -e "h\tTimer\tCount\tTotal\tMax\tMin\tMean\tStdDev" > $out
  for i in $title; do
    h=$(echo $i| sed "s/-/ /g" | awk '{print $3}')
    r=$(grep -h -m 1 "Inverse" $i)
    if [ ! -z "$r" ]; then
      echo -e "$h\t$r" >> $out
    fi
  done
}
function plot() {
  out=OUT_${2}D-mu-${3}-*-${13}_${11}
  echo "set logscale"
  echo "set title \"Res from $out\""
  echo "set datafile missing \"-\""
  echo "set xtics nomirror rotate by -45"
  echo "set key noenhanced"
  echo -n "plot ";
  for i in `ls $out`;
    do 
      echo -n "\"$i\" using \"hSize\":\"err_rel\" w lp,"; 
  done
  echo ""
}

NPROCS=6
OUTFILE=res.txt
appDir=`pwd`
h=1
for D in `seq 2 3`;
do
  for poly in `echo poly sin`;
  do
    for mu in `perl -le'for my $i (0..7) { print 10**-$i }'`;
    do
      # LU
      res $NPROCS $D $mu $h cg lu cg lu cg lu $OUTFILE $appDir $poly
      nbIter $NPROCS $D $mu $h cg lu cg lu cg lu $OUTFILE $appDir $poly
      # Block : LU LU
      res $NPROCS $D $mu $h cg blockms cg lu cg lu $OUTFILE $appDir $poly
      nbIter $NPROCS $D $mu $h cg blockms cg lu cg lu $OUTFILE $appDir $poly
      # Block : Gamg Gamg
      res $NPROCS $D $mu $h cg blockms cg gamg cg gamg $OUTFILE $appDir $poly
      nbIter $NPROCS $D $mu $h cg blockms cg gamg cg gamg $OUTFILE $appDir $poly
      # Block : hypre hypre
      res $NPROCS $D $mu $h cg blockms cg hypre cg hypre $OUTFILE $appDir $poly
      nbIter $NPROCS $D $mu $h cg blockms cg hypre cg hypre $OUTFILE $appDir $poly
      # Generate Gnuplot Command
      # plot $NPROCS $D $mu $h gmres blockms gmres gamg gmres gamg $OUTFILE $appDir $poly
    done
  done
done

