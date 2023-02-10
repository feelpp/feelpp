#!/bin/bash

# dataFolder="$~/feelppdb/nirb/heat"
# dataFolder="${HOME}/Documents/FEELTEST/nirb/heat/np_1"

# removeOldDatas="rm -r ${dataFolder}"

# ${removeOldDatas}


Nproc="1 2 4 8 16 32"
for n in $Nproc;
do


echo " ------------------------------------------ "
echo "  Restarting program with nb procs = : $n "
echo " ------------------------------------------ "

test_para="mpiexec -n $n -bind-to core python3 test_perf_nirb.py --config-file model/square/square.cfg --savetime 1 --Ntest 50"

# test_para="mpiexec -n $n -bind-to core python3 test_perf_nirb.py --config-file model/thermal-fin-3d/thermal-fin.cfg --savetime 1"

${test_para}


done
