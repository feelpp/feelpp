#!/bin/bash

dataFolder="$~/feelppdb/nirb/heat"
# dataFolder="${HOME}/Documents/FEELTEST/nirb/heat/np_1"

removeOldDatas="rm -r ${dataFolder}"

# ${removeOldDatas}

Nsnap="1 2 4 6 10 12 14 16 20 25 30 35 40 45 50 60 70 100 150 200"

Rectification=0
biorthonormal=0

echo "N timeToolbox timeNirb" > "${dataFolder}/nirbOnline_time_exec.dat"
# echo "N l2_min linf_min l2_mean linf_mean l2_max linf_max" > "${dataFolder}/nirb_error.dat"
echo "N l2 linf" > "${dataFolder}/nirb_error.dat"

for n in $Nsnap;
do

Ns=$n

offline="python3 nirbOffline.py ${Ns}"
online="python3 nirbOnline.py"

${offline}
${online}

echo " ------------------------------------------ "  
echo "  Restarting program with Ns = : $Ns+1 "  
echo " ------------------------------------------ "


done 
exit 
