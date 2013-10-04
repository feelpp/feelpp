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
if ( FEELPP_ENABLE_SCHED_LOADLEVELER )
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${execname}.ll "#@ shell = /usr/bin/ksh
#@ job_type = MPICH
#@ environment = COPY_ALL
#@ island_count=1
#@ energy_policy_tag = NONE
#@ job_name = ${execname}
#@ class = test
#@ wall_clock_limit = 00:40:00
#@ network.MPI = sn_all,not_shared,us
#@ notification = never
#@ initialdir = ${CMAKE_CURRENT_BINARY_DIR}
#@ output = ${execname}.$(jobid).out
#@ error =  ${execname}.$(jobid).err
#@ node = 4
#@ total_tasks =  64
#@ queue

")
  if ( FEELPP_APP_CFG )
    foreach(  cfg ${FEELPP_APP_CFG} )
      get_filename_component( CFG_NAME ${cfg} NAME )
      file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${execname}.ll "
mpirun --bind-to-core ${execname} --config-file=${cfg}  # add other feel++ options here ")
    endforeach()
  else( FEELPP_APP_CFG )
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${execname}.ll "
mpirun --bind-to-core ${execname} # add other feel++ options here")
  endif( FEELPP_APP_CFG )
endif( FEELPP_ENABLE_SCHED_LOADLEVELER )
