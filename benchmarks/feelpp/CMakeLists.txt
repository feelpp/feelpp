# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

add_custom_target( benchs )

add_subdirectory( functionspace_manager )

feelpp_add_application( env SRCS 01-env.cpp LINK_LIBRARIES  benchmark EXEC B_ENV )
feelpp_add_application( mesh SRCS 02-mesh.cpp LINK_LIBRARIES  benchmark EXEC B_MESH )
feelpp_add_application( space SRCS 03-space.cpp LINK_LIBRARIES  benchmark EXEC B_SPACE )
feelpp_add_application( integrate SRCS 04-integrate.cpp LINK_LIBRARIES  benchmark EXEC B_INTEGRATE )

feelpp_add_application( basis SRCS 05-basis.cpp LINK_LIBRARIES benchmark EXEC B_BASIS )

add_subdirectory( feelpoly )

if ( FEELPP_ENABLE_SPECX )
    feelpp_add_application( specx SRCS 99-specx.cpp LINK_LIBRARIES Feelpp::feelpp_specx )
endif()
add_dependencies( benchs ${B_ENV} ${B_MESH} ${B_SPACE} ${B_INTEGRATE} ${B_BASIS}
                  ${B_FSPACE_MANAGER_BASELINE} ${B_FSPACE_MANAGER} benchs-feelpoly )
