###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2011-04-23
#
#  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
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
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

#
# opus_models: a library of models
#
add_library(feel++_opus_models SHARED
  eads.cpp

  opusdata.cpp opuscomponent.cpp opusmodelbase.cpp
  opusmodelfactory.cpp opusmodelthermal.cpp
  opusmodelfluid.cpp
#  opusmodelfluidoseen.cpp
  opusmodel.cpp

#  opusmodel_211.cpp
  opusmodel_212.cpp
#  opusmodel_213.cpp
#  opusmodel_214.cpp

  opusmodelrb_212.cpp

  eadsmfemapp.cpp
  eadsscmapp.cpp
)
target_link_libraries(feel++_opus_models ${FEEL_LIBRARIES} ${MPI_LIBRARIES} )
set_target_properties( feel++_opus_models PROPERTIES VERSION ${FEELPP_OPUS_SHARED_VERSION} SOVERSION ${FEELPP_OPUS_SHARED_SOVERSION})
install(TARGETS feel++_opus_models DESTINATION lib COMPONENT Bin)
set_property(TARGET feel++_opus_models PROPERTY LABELS opus)
add_dependencies(opus feel++_opus_models)

#
# executables
#
# mfem
opus_add_executable( eadsmfem eadsmfemmain.cpp LINK_LIBRARIES feel++_opus_models ${FEEL_LIBRARIES} CFG eadsmfem.cfg )
# pfem
#opus_add_executable( eadspfem eadspfemmain.cpp LINK_LIBRARIES feel++_opus_models ${FEEL_LIBRARIES} CFG eadspfem.cfg )
opus_add_executable( eadspfem eadspfemmain.cpp  LINK_LIBRARIES feel++_opus_models ${FEEL_LIBRARIES} CFG eadspfem.cfg )
# SCM
opus_add_executable( eadsscm eadsscmmain.cpp LINK_LIBRARIES feel++_opus_models   ${FEEL_LIBRARIES} CFG eadsscm.cfg )
# CRB
opus_add_executable( eadscrb eadscrbmain.cpp LINK_LIBRARIES  feel++_opus_models ${FEEL_LIBRARIES} CFG eadscrb.cfg )

add_subdirectory(python)
add_subdirectory(octave)
