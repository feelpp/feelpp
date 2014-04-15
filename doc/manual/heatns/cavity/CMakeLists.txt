###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2012-07-27
#
#  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)
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

feelpp_add_application( heatns_natural_convection_cavity_2d SRCS
  convection_main.cpp convection_crb_specifics.cpp ${lib_convection} CFG
  convection.cfg LABELS doc DEFS FEELPP_USE_LM CONVECTION_DIM=2 NO_TEST )

feelpp_add_application( heatns_natural_convection_cavity_2d_crb SRCS
  convection_main.cpp convection_crb_specifics.cpp ${lib_convection} CFG
  convection.cfg LABELS doc DEFS FEELPP_USE_LM CONVECTION_DIM=2 CRB_SOLVER=1
  NO_TEST)

feelpp_add_application(  heatns_natural_convection_cavity_3d SRCS convection_main.cpp
  convection_crb_specifics.cpp ${lib_convection} CFG convection.cfg
  LABELS doc DEFS FEELPP_USE_LM CONVECTION_DIM=3 CRB_SOLVER=0
  NO_TEST )

feelpp_add_application( heatns_natural_convection_cavity_3d_crb SRCS
  convection_main.cpp convection_crb_specifics.cpp ${lib_convection} CFG
  convection.cfg LABELS doc DEFS FEELPP_USE_LM CONVECTION_DIM=3  CRB_SOLVER=1
  NO_TEST )

feelpp_add_application( P3P2P3_heatns_natural_convection_cavity_2d_crb
  SRCS convection_main.cpp convection_crb_specifics.cpp
  ${lib_convection} CFG convection.cfg LABELS doc DEFS FEELPP_USE_LM
  CONVECTION_DIM=2 CRB_SOLVER=1 CONVECTION_ORDER_U=3  CONVECTION_ORDER_P=2 CONVECTION_ORDER_T=3
  NO_TEST )

feelpp_add_application( P3P2P3_heatns_natural_convection_cavity_3d
  SRCS convection_main.cpp convection_crb_specifics.cpp
  ${lib_convection} CFG convection.cfg LABELS doc DEFS FEELPP_USE_LM
  CONVECTION_DIM=3 CRB_SOLVER=0 CONVECTION_ORDER_U=3   CONVECTION_ORDER_P=2 CONVECTION_ORDER_T=3
  NO_TEST )

feelpp_add_application( P3P2P3_heatns_natural_convection_cavity_3d_crb
  SRCS convection_main.cpp convection_crb_specifics.cpp
  ${lib_convection} CFG convection.cfg LABELS doc DEFS FEELPP_USE_LM
  CONVECTION_DIM=3 CRB_SOLVER=1 CONVECTION_ORDER_U=3  CONVECTION_ORDER_P=2 CONVECTION_ORDER_T=3
  NO_TEST )
