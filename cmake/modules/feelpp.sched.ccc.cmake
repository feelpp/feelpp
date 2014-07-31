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
  if (FEELPP_ENABLE_SCHED_CCC )
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${execname}.msub "#! /bin/bash
#MSUB -r ${execname}         # Request name
#MSUB -n 64                  # Number of tasks to use
#MSUB -T 1800                # Elapsed time limit in seconds of the job (default: 1800)
#MSUB -o ${execname}_%I.o    # Standard output. %I is the job id
#MSUB -e ${execname}_%I.e    # Error output. %I is the job id
#MSUB -A gen7334  # Hifimagnet
#MSUB -A gen7335  # Bloodflow
#MSUB -A ra0840              # Project ID
#MSUB -q standard            # Choosing large nodes
##MSUB -@ noreply@cea.fr:end # Uncomment this line for being notified at the end of the job by sending a mail at the given address

#set -x
cd \${BRIDGE_MSUB_PWD}        # BRIDGE_MSUB_PWD is a environment variable which contains the directory where the script was submitted
unset LC_CTYPE
ccc_mprun ${execname}  # you can add Feel++ options here
")
  endif()
