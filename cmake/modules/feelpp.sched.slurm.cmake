###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
#       Date: 2013-04-26
#
#  Copyright (C) 2013 Université de Strasbourg
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
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm "#! /bin/bash
# type 'sbatch ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm' to submit the job,
# sbatch will pich the number of cores
# given in the script by '#SBATCH -n xxx', you can override this value by
# typing 'sbatch -n xxxx ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm'

#SBATCH -n 32 #need 32 cores (one thread by core)
source $HOME/.bash_profile
unset LC_CTYPE

#export IMPORTANT_VAR=important_value

#cd /workdir/math/whoami
cd ${CMAKE_CURRENT_BINARY_DIR}/
")
  if ( FEELPP_APP_CFG )
    foreach(  cfg ${FEELPP_APP_CFG} )
      get_filename_component( CFG_NAME ${cfg} NAME )
      file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm "
mpirun --bind-to-core -x LD_LIBRARY_PATH ${CMAKE_CURRENT_BINARY_DIR}/${execname} --config-file=${NAME}  # add other fel++  options here")
    endforeach()
  else( FEELPP_APP_CFG )
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${execname}.slurm
"mpirun --bind-to-core -x LD_LIBRARY_PATH ${CMAKE_CURRENT_BINARY_DIR}/${execname} # add other feel++ options here")
  endif( FEELPP_APP_CFG )
endif( FEELPP_ENABLE_SCHED_SLURM )
