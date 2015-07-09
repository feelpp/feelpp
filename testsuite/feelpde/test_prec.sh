#!/bin/bash
# set -x
# This script is to test the precAFP behavior

# $1 - nbProcs
# $2 - Dim
# $3 - mu
# $4 - hsize
# $5-11 - ksp/pc config
# $12 - string
function simu() {
  title=${2}D-h-${4}-mu-${3}-$5-${6}_$7-${8}_${9}-${10}
  out=${title}.batch
  rm $out; touch $out
  echo "#!/bin/bash" >>$out
  echo "# Lines with SBATCH starting with ## are comments and starting with # are actual commands for sbatch">>$out
  echo "#SBATCH --job-name  ${2}D-h-${4}-mu-${3}-$5-${6}_$7-${8}_${9}-${10}">>$out
  echo "source \$HOME/.zsh/my/atlas">>$out
  echo "unset LC_CTYPE">>$out
  echo "##SBATCH -p public">>$out
  echo "# number of cores">>$out
  echo "##SBATCH -n $1">>$out
  echo "# min-max number of nodes">>$out
  echo "##SBATCH -N 1-4">>$out
  echo "# max time of exec (will be killed afterwards)">>$out
  echo "##SBATCH -t 12:00:00">>$out
  echo "# number of tasks per node">>$out
  echo "##SBATCH --tasks-per-node 1">>$out
  echo "# specify execution constraitns">>$out
  echo "##SBATCH --constraint \"intel\"">>$out
  echo "# min mem size">>$out
  echo "##SBATCH --mem=16684">>$out
  echo "# display info about cpu binding">>$out
  echo "##SBATCH --cpu_bind=verbose">>$out
  echo "# send a mail at the end of the exec">>$out
  echo "#SBATCH --mail-type=END">>$out
  echo "#SBATCH --mail-user=vincent.huber@cemosis.fr">>$out
  echo "">>$out
  echo "# If you want to have access to Feel++ logs">>$out
  echo "# export the FEELPP_SCRATCHDIR variable to an NFS mounted directory">>$out
  echo "export FEELPP_SCRATCHDIR=/scratch/job.\${SLURM_JOB_ID}/log">>$out
  echo "export FEELPP_WORKDIR=/scratch/job.\${SLURM_JOB_ID}/export">>$out
  echo "">>$out
  echo "#################### OPTIONAL: ">>$out
  echo "module load gcc490.profile" >>$out
  echo "#################### OPTIONAL:">>$out
  echo "">>$out
  echo "# Finally launch the job">>$out
  echo "# mpirun of openmpi is natively interfaced with Slurm">>$out
  echo "# No need to precise the number of processors to use">>$out
  echo "cd ${12}">>$out
  echo "mpirun --bind-to core -x LD_LIBRARY_PATH ./feelpp_test_precAFP${2}D --config-files precAFC${2}D.cfg backend.cfg --functions.m ${3} --gmsh.hsize ${4} --ms.pc-type=$5 --ms.ksp-type=$6 --ms.blockms.11.pc-type=$7 --ms.blockms.11.ksp-type=$8 --ms.blockms.22.pc-type=$9 --ms.blockms.22.ksp-type=${10} --title $title --generateMD true > ${title}_${11}">>$out

  echo "mkdir -p /data/`whoami`/prec_behavior/${title} ">>$out
  echo "cp -r /scratch/job.\${SLURM_JOB_ID}/* /data/`whoami`/prec_behavior/${title} ">>$out
  echo "cp ${title}_${11} /data/`whoami`/prec_behavior/${title} ">>$out
  echo "cp *md /data/`whoami`/prec_behavior/${title} ">>$out
}

NPROCS=6
OUTFILE=res.txt
appDir=`pwd`
for D in `seq 2 3`;
do
  for mu in `perl -le'for my $i (0..7) { print 10**-$i }'`;
  do
    for h in `perl -le'for my $i (1..7) { print 1/(2**$i) }'`; 
    do
      # LU
      simu $NPROCS $D $mu $h gmres lu gmres lu gmres lu $OUTFILE $appDir
      # Block : LU LU
      simu $NPROCS $D $mu $h gmres blockms gmres lu gmres lu $OUTFILE $appDir
      # Block : Gamg Gamg
      simu $NPROCS $D $mu $h gmres blockms gmres gamg gmres gamg $OUTFILE $appDir
    done
  done
done

