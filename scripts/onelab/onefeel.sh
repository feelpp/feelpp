#!/bin/bash
# Script executing a Feel++ app with Gmsh through OneLab
# Initial contributors: Carolina Diaz, Jérôme Boeglin, Sébastien Lan
# Authors: Alexandre Ancel

PCONFIG=""

#	Argument check
if [[ $# -lt 1 ]]
then
	echo "Usage : $0 <Application Path> [options]"
	exit
fi

#	Check that the app exists
if [[ ! -f "$1" ]]; then
    echo "$1 does not exist, please check the path you have entered"
    exit
fi

if [[ -f "$1.onelab.cfg" ]]; then
		echo "Restoring previous configuration from \"$1.onelab.cfg\""
		PCONFIG="--config-file $1.onelab.cfg"
fi

# Launch the feel++ application to generate the OneLab files
# echo "$* --generate-ol ${PCONFIG}"
eval "$* --generate-ol ${PCONFIG}"

# echo $?

# Launch Gmsh with the OneLab files
#Lancement de Gmsh avec les fichiers crées 
gmsh $1.ol
