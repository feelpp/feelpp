# # On "allume" les dépots (les deux premiers ne fonctionnent pas)
# sudo rm /etc/yum.repos.d/rpm*
# sudo rm /etc/yum.repos.d/slc6-cernonly*
# sudo sed -i "s/enabled=0/enabled=1/g" -e /etc/yum.repos.d
# sudo yum upgrade
# 
# # Install de quelques paquets nécessaires.
# sudo yum -y update
# sudo yum -y update yum
# sudo yum -y install cmake28 texinfo glibc-devel.i686 icu.x86_64 libicu-devel.x86_64 libmpc.x86_64 mpfr.x86_64 gmp-devel.x86_64 mpfr-devel.x86_64 libmpc-devel.x86_64

#Pour installer ailleurs: changer la variable $workdir (s'assurer d'avoir les droits)
#espace nécessaire: 10Go (possibilité de minimiser beaucoup en supprimant les repertoires de compilation)
export cmake=cmake28
export workdir=/usr/local/feelpp/work
mkdir $workdir

export gccVersion=4.7.2
export gccDir=$work/gcc-$gccVersion
export openmpiDir=$work/openmpi/gcc-$gccVersion/
export boostDir=$work/boost/gcc-$gccVersion/
export gmshDir=$work/gmsh/gcc-$gccVersion/
export petscDir=$work/petscDir/gcc-$gccVersion/
export feelppDir=$work/feelppDir/gcc-$gccVersion/

export nbProc=`cat /proc/cpuinfo |grep processor | tail -n 1 | awk '{print $3}'|awk '{print "("$1"+1)/3"}' | bc`

#télécharge, compile et install gcc dans $workdir/gcc
sh gcc47.sh
#télécharge, compile et install openmpi dans $workdir/openmpi
so openmpi.sh
#télécharge, compile et install boost dans $workdir/boost
sh boost.sh
#télécharge, compile et install gmsh dans $workdir/gmsh
sh gmsh.sh
#télécharge, compile et install petsc dans $workdir/petsc-version
sh petsc.sh
#télécharge, compile et install feelpp dans $workdir/feelppLib
sh feelpp.sh
