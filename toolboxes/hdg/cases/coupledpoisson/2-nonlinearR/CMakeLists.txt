###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 30 Mar 2016
#
#  Copyright (C) 2016 Feel++ Consortium
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

feelpp_add_fmu( test3d0d_nonlinearR CLASS test3d0d_nonlinearR SRCS test3d0d_nonlinearR.mo nonlinearResistor.mo voltage_nonlinearTest3d0d.mo voltageInterfaceSource.mo currentInterfaceSource.mo CATEGORY toolboxes/hdg/cases/coupledpoisson/2-nonlinearR )
