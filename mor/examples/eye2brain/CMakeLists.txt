# ##  CMakeLists.txt; coding: utf-8 ---

# Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#            Thomas Saigre <saigre@math.unistra.fr>
# Date: 01 April 2023
#
# Copyright (C) 2023 Feel++ Consortium
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

mor_add_library(eye2brain SRCS eye2brain.cpp EXEC CRB_LIB_EYE2BRAIN PLUGIN )
target_compile_definitions(${CRB_LIB_EYE2BRAIN} PRIVATE FEELPP_MOR_PLUGIN_NAME=eye2brain FEELPP_MOR_PLUGIN_LIBNAME=${CRB_LIB_EYE2BRAIN})
# install(FILES eye2brain.hpp DESTINATION include/feelpp/crb/eye2brain COMPONENT Devel)

set(FEELPP_CRB_EYE2BRAIN_LIBRARIES ${CRB_LIB_EYE2BRAIN} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_EYE2BRAIN_LIBRARIES PARENT_SCOPE)

# mor_add_model(eye2brain Eye2Brain HDRS eye2brain.hpp  LINK_LIBRARIES ${CRB_LIB_EYE2BRAIN} )
add_custom_target(feelpp_mor_eye2brain ALL)

foreach( ORDER 1 2 )
foreach( DIM 2 3 )
    # configure_file( eye2brain_${D}d.cfg.in eye2brain_${D}d.cfg @ONLY )
    mor_add_model(
        eye2brain_${DIM}dP${ORDER}
        Eye2Brain
        CLASS Eye2Brain<${ORDER}, ${DIM}>
        HDRS eye2brain.hpp
        OPT eye2brain
        LINK_LIBRARIES ${CRB_LIB_EYE2BRAIN}
    )
    add_dependencies(feelpp_mor_eye2brain feelpp_mor_eye2brain_${DIM}dP${ORDER}app)
endforeach()
endforeach()


feelpp_add_testcase(eye2brain PREFIX feelpp_crb CATEGORY CRB)

if (FEELPP_ENABLE_OPENTURNS)
     feelpp_add_application( sensitivity_analysis SRCS sensitivity_analysis.cpp
         PROJECT mor
         LINK_LIBRARIES Feelpp::feelpp_mor Feelpp::feelpp_openturns tbb
     )
endif()

#[==[
## test : disabled for the runner : meshed are only present locally
foreach( ORDER 1 2 )

    set (TEST_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/test_2dP${ORDER}.sh")
    file (WRITE ${TEST_SCRIPT}
        "#!/bin/bash\n"
        "./feelpp_mor_eye2brain_2dP${ORDER}app --config-file ${CMAKE_CURRENT_SOURCE_DIR}/eye2brain/eye2brain-2d.cfg\
        --eye2brain.model-name eye2brain_2dP${ORDER}_test --eye2brain.run.mode 2 --eye2brain.run.sampling.size 4\
        --crb.dimension-max 3\n"
        "../../apps/feelpp_mor_onlinerun --crbmodel.name eye2brain_2dP${ORDER}_test --plugin.dir ${CMAKE_BINARY_DIR}/install/lib\n"
    )

    add_test( NAME feelpp_mor_eye2brain_2dP${ORDER}_test
        COMMAND bash ${TEST_SCRIPT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    set (TEST_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/test_3dP${ORDER}.sh")
    file (WRITE ${TEST_SCRIPT}
        "#!/bin/bash\n"
        "mpirun -np 12 ./feelpp_mor_eye2brain_3dP${ORDER}app --config-file ${CMAKE_CURRENT_SOURCE_DIR}/eye2brain/eye2brain-3d.cfg\
        --eye2brain.model-name eye2brain_3dP${ORDER}_test --eye2brain.run.mode 2 --eye2brain.run.sampling.size 4\
        --crb.dimension-max 3\n"
        "../../apps/feelpp_mor_onlinerun --crbmodel.name eye2brain_3dP${ORDER}_test --plugin.dir ${CMAKE_BINARY_DIR}/install/lib\n"
    )

    add_test( NAME feelpp_mor_eye2brain_3dP${ORDER}_test
        COMMAND bash ${TEST_SCRIPT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

endforeach()
#]==]