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
  echo "*************"
  grep RES $title > $out
  echo "*************"
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
      res $NPROCS $D $mu $h gmres lu gmres lu gmres lu $OUTFILE $appDir $poly
      # Block : LU LU
      res $NPROCS $D $mu $h gmres blockms gmres lu gmres lu $OUTFILE $appDir $poly
      # Block : Gamg Gamg
      res $NPROCS $D $mu $h gmres blockms gmres gamg gmres gamg $OUTFILE $appDir $poly
    done
  done
done

