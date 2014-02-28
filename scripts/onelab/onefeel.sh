#!/bin/bash
# Script executing a Feel++ app with Gmsh through OneLab
# Initial contributors: Carolina Diaz, Jérôme Boeglin, Sébastien Landre
# Authors: Alexandre Ancel

PCONFIG=""

#	Argument check
if [[ $# -lt 2 ]]
then
	echo "Usage : $0 <Nprocs> <Application Path> [options]"
	#echo "Usage : $0 <Application Path> [options]"
	exit
fi

# Store desired number of processors
NPROCS=$1
shift

if [[ $NPROCS -lt 1 ]]; then
    echo "Invalid number of MPI processes ($NPROCS)"
    exit
fi

#	Check that the app exists
#if [[ ! -f "$1" ]]; then
#    echo "$1 does not exist, please check the path you have entered"
#    exit
#fi

# if a onelab.cfg file exists
if [[ -f "$1.onelab.cfg" ]]; then
	# we do not restore the previous configuration for the moment
    #echo "Restoring previous configuration from \"$1.onelab.cfg\""
    #PCONFIG="--config-file $1.onelab.cfg"
    PCONFIG=""
fi

# handle MPI
if [[ $NPROCS -eq 1 ]]; then
    # single process code
    eval "$* --generate-ol ${PCONFIG}"
else
    # Launch the feel++ application to generate the OneLab files
    eval "mpirun -np ${NPROCS} $* --generate-ol ${PCONFIG}"
fi

if [[ $? -ne 0 ]]; then
	exit 1
fi

# Launch Gmsh with the OneLab files
~/svn/gmsh/build/gmsh $1.ol
