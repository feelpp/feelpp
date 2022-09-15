#!/bin/bash

dataFolder="$~/feelppdb/nirb/heat"

removeOldDatas="rm -r ${dataFolder}"

# ${removeOldDatas}

Nsnap="2 4 8 16 32 64 128"

Rectification=0
biorthonormal=0


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



