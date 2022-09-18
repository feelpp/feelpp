#!/bin/bash

dataFolder="$~/feelppdb/nirb/heat"

removeOldDatas="rm -r ${dataFolder}"

# ${removeOldDatas}

Nsnap="1 2 4 6 10 12 14 16 20 25 30 35 40 45 50 70 80 100"
Rectification=0
biorthonormal=0


for n in $Nsnap;
do

Ns=$n

echo " ------------------------------------------ "  
echo "  Restarting program with Ns = : $Ns "  
echo " ------------------------------------------ "

offline="python3 nirbOffline.py ${Ns}"
online="python3 nirbOnline.py"

${offline}
${online}



done 
exit 



