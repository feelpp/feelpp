#!/bin/bash

test="rectification"
test=$1

# if [ $test=="parallel" ];
# then

# echo " ------------------------------------------ "  
# echo "          Run parallel test : $test         "  
# echo " ------------------------------------------ "

# ### test parallel computing 
# Nproc="1 2 4 8 16 32"
# for n in $Nproc;
# do


# echo " ------------------------------------------ "  
# echo "  Restarting program with nb procs = : $n "  
# echo " ------------------------------------------ "

# run="mpiexec -n $n -bind-to core python3 test_perf_nirb.py --config-file model/square/square.cfg --timeexec 1 --convergence 0 --Ntest 50"

# # test_para="mpiexec -n $n -bind-to core python3 test_perf_nirb.py --config-file model/thermal-fin-3d/thermal-fin.cfg --savetime 1"

# ${run}
# done 
# fi 

if [ $test=='rectification' ];
then 

echo " ------------------------------------------ "  
echo "          Run rectificatio test $test       "  
echo " ------------------------------------------ "

#### test rectification parameter 
Nl="0 1 3 5 7 9 10 11 12"  # lambda = 1.e^{Nl}
Nl="3"
model="/data/home/elarif/nirbDatas/model/square/square4/heat-cube.cfg"

for n in $Nl;
do


echo " ------------------------------------------ "  
echo "  Restarting program with lambda = : 1.E-$n "  
echo " ------------------------------------------ "

run="python3 test_perf_nirb.py --config-file $model --timeexec 0 --idmodel s4 --Ntest 30 --lmdexp $n"

# test_para="mpiexec -n $n -bind-to core python3 test_perf_nirb.py --config-file model/thermal-fin-3d/thermal-fin.cfg --savetime 1"

${run}


done 
fi 
