#!/bin/bash

dataFolder="$~/feelppdb/nirb/heat"

removeOldDatas="rm -r ${dataFolder}"

# ${removeOldDatas}

Nsnap="1 2 4 6 10 12 14 16 20 25 30 35 40 45 50 60 70"

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



