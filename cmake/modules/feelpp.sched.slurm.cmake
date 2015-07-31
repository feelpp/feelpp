###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2013-04-26
#
#  Copyright (C) 2013-2015 Feel++ Consortium
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
if ( FEELPP_ENABLE_SCHED_SLURM )

    SET(SLURM_NCORES "16")
    # attempt to get the number of cores from slurm on atlas 
    if(FEELPP_MACHINE_NAME MATCHES "irma-atlas")
        execute_process(COMMAND sinfo -h -o "%C" OUTPUT_VARIABLE SINFO_CORES0)
        # remove end of lines
        STRING(REGEX REPLACE "\n" "" SINFO_CORES ${SINFO_CORES0})
        # get the maximum number of cores available if we were able to obtain it
        if(SINFO_CORES)
            string(REPLACE "/" ";" SINFO_CORES_LIST ${SINFO_CORES})
            list(GET SINFO_CORES_LIST 3 SLURM_NCORES) 
        endif(SINFO_CORES)
    endif()

    # write header of file
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm 
"#! /bin/bash
# type 'sbatch ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm' to submit the job,
# sbatch will pich the number of cores
# given in the script by '#SBATCH -n xxx', you can override this value by
# typing 'sbatch -n xxxx ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm'

#SBATCH -n ${SLURM_NCORES} #need ${SLURM_NCORES} cores (one thread by core)
source $HOME/.bash_profile
unset LC_CTYPE

#export IMPORTANT_VAR=important_value
")

    if(FEELPP_MACHINE_NAME MATCHES "irma-atlas")
        file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm 
"
# In case you want to use modules.
# You first have to activate the module command
# by uncommenting the following line
# source /etc/profile.d/modules.sh

# Source the configuration for Feel++ or your custom configuration
# By uncommenting the following lines
# PREVPATH=`pwd`
# cd /data/software/config/etc
# source feelpprc.sh
# cd ${PREVPATH}

# Load the latest feelpp dependencies modules here
module load feelpp.profile

#cd /workdir/math/whoami
cd ${CMAKE_CURRENT_BINARY_DIR}/
")

    endif()

    # If we use a cfg specify it
    if ( FEELPP_APP_CFG )
        foreach( cfg ${FEELPP_APP_CFG} )
            get_filename_component( CFG_NAME ${cfg} NAME )
            file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm 
            "mpirun --bind-to core ${CMAKE_CURRENT_BINARY_DIR}/${execname} --config-file=${cfg}  
# "
            )
        endforeach()
    # Otherwise do not provide the config-file option
    else( FEELPP_APP_CFG )
        file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm
            "mpirun --bind-to core ${CMAKE_CURRENT_BINARY_DIR}/${execname} 
# "
        )
    endif( FEELPP_APP_CFG )
endif( FEELPP_ENABLE_SCHED_SLURM )
