#!/bin/bash
# set -x
# This script is to test the precAFP behavior

declare -a ksp_pc_map

# $1 : kind of simu
function setKspPc(){
ksp_pc_map[0]=nan 
ksp_pc_map[1]=nan 
ksp_pc_map[2]=nan 
ksp_pc_map[3]=nan 
ksp_pc_map[4]=nan 
ksp_pc_map[5]=nan 
ksp_pc_map[6]=nan 
ksp_pc_map[7]=nan 
ksp_pc_map[8]=nan 
ksp_pc_map[9]=nan 
case $1 in 
  0)
  ksp_pc_map[0]=minres
  ksp_pc_map[1]=lu 
  ;;
  1)
  ksp_pc_map[0]=minres
  ksp_pc_map[1]=blockms
  ksp_pc_map[2]=cg
  ksp_pc_map[3]=lu
  ksp_pc_map[8]=cg
  ksp_pc_map[9]=lu
  ;;
  2)
  ksp_pc_map[0]=minres
  ksp_pc_map[1]=blockms
  ksp_pc_map[2]=cg
  ksp_pc_map[3]=lu
  ksp_pc_map[8]=cg
  ksp_pc_map[9]=boomeramg
  ;;
  3)
  ksp_pc_map[0]=minres
  ksp_pc_map[1]=blockms
  ksp_pc_map[2]=cg
  ksp_pc_map[3]=AS
  ksp_pc_map[4]=cg
  ksp_pc_map[5]=boomeramg
  ksp_pc_map[6]=cg
  ksp_pc_map[7]=boomeramg
  ksp_pc_map[8]=cg
  ksp_pc_map[9]=boomeramg
  ;;
  4)
  ksp_pc_map[0]=minres
  ksp_pc_map[1]=blockms
  ksp_pc_map[2]=cg
  ksp_pc_map[3]=ams
  ksp_pc_map[4]=nan
  ksp_pc_map[5]=nan
  ksp_pc_map[6]=nan
  ksp_pc_map[7]=nan
  ksp_pc_map[8]=cg
  ksp_pc_map[9]=boomeramg
  ;;
esac
}

# $1 - kind of simu
function getConf()
{
  setKspPc $1
  echo "${ksp_pc_map[0]}-${ksp_pc_map[1]}_${ksp_pc_map[2]}-${ksp_pc_map[3]}_${ksp_pc_map[4]}-${ksp_pc_map[5]}_${ksp_pc_map[6]}-${ksp_pc_map[7]}_${ksp_pc_map[8]}-${ksp_pc_map[9]}"
}

function simu() {
  title=${2}D-h-${4}-mu-${3}-$5-${6}_$7-${8}_${9}-${10}-${13}
  echo "mpirun -np $1 ./feelpp_test_precAFP${2}D --config-files precAFP${2}D_${13}.cfg backend.cfg --functions.m ${3} --gmsh.hsize ${4} --ms.ksp-type=$5 --ms.pc-type=$6 --ms.blockms.11.ksp-type=$7 --ms.blockms.11.pc-type=$8 --ms.blockms.22.ksp-type=$9 --ms.blockms.22.pc-type=${10} --title $title --generateMD true > ${title}_${11}"
}

# $1 - nbProcs
# $2 - Dim
# $3 - mu
# $4 - hsize
# $5 - OUTFILE
# $6 - appDir
# $7 - test type : lin or sin
# $8 - ksp/pc config
function simuBatch() {
  setKspPc $8
  ksp_pc=$(getConf $8)
  title=${2}D-h-${4}-mu-${3}-${ksp_pc}-${7}
  out=${title}.batch
  echo $out
  rm $out; touch $out
  echo "#!/bin/bash" >>$out
  echo "# Lines with SBATCH starting with ## are comments and starting with # are actual commands for sbatch">>$out
  echo "#SBATCH --job-name  $title">>$out
  echo "source \$HOME/.zsh/my/atlas">>$out
  echo "unset LC_CTYPE">>$out
  echo "##SBATCH -p public">>$out
  echo "# number of cores">>$out
  echo "#SBATCH -n $1">>$out
  echo "# min-max number of nodes">>$out
  echo "##SBATCH -N 1-4">>$out
  echo "# max time of exec (will be killed afterwards)">>$out
  echo "#SBATCH -t 01:00:00">>$out
  echo "# number of tasks per node">>$out
  echo "##SBATCH --tasks-per-node 4">>$out
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
  echo "module load feelpp.profile" >>$out
  echo "#################### OPTIONAL:">>$out
  echo "">>$out
  echo "# Finally launch the job">>$out
  echo "# mpirun of openmpi is natively interfaced with Slurm">>$out
  echo "# No need to precise the number of processors to use">>$out
  echo "cd ${12}">>$out
  echo "mpirun --bind-to core -x LD_LIBRARY_PATH ./feelpp_test_precAFP${2}D --config-files precAFP${2}D_${7}.cfg backend.cfg --functions.m ${3} --gmsh.hsize ${4} --ms.ksp-type=${ksp_pc_map[0]} --ms.pc-type=${ksp_pc_map[1]} --ms.blockms.11.ksp-type=${ksp_pc_map[2]} --ms.blockms.11.pc-type=${ksp_pc_map[3]} --ms.blockms.11.1.ksp-type=${ksp_pc_map[4]} --ms.blockms.11.1.pc-type=${ksp_pc_map[5]} --ms.blockms.11.2.ksp-type=${ksp_pc_map[6]} --ms.blockms.11.2.pc-type=${ksp_pc_map[7]} --ms.blockms.22.ksp-type=${ksp_pc_map[8]} --ms.blockms.22.pc-type=${ksp_pc_map[9]} --title $title --generateMD true --saveTimers true > ${title}_${5}" >>$out
  echo "mkdir -p /data/`whoami`/prec_behavior/${title} ">>$out
  echo "cp -r /scratch/job.\${SLURM_JOB_ID}/* /data/`whoami`/prec_behavior/${title} ">>$out
  echo "cp ${title}_${11} /data/`whoami`/prec_behavior/${title} ">>$out
  echo "cp *md /data/`whoami`/prec_behavior/${title} ">>$out
}

NPROCS=1
OUTFILE=res.txt
appDir=`pwd`
# Dimension
for D in `seq 2 3`;
do
  # test case 
  for poly in `echo poly sin`;
  do
    # Permeability
    for mu in `perl -le'for my $i (0..7) { print 10**-$i }'`;
    do
      # hsize
      for h in `perl -le'for my $i (1..7) { print 1/(2**$i) }'`; 
      do
        # ksp/pc conf (see setKspPc)
        for kind in `seq 0 4`;
        do
          simuBatch $NPROCS $D $mu $h $OUTFILE $appDir $poly $kind
        done
      done
    done
  done
done

