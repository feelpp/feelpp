# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

add_custom_target( benchs-feelpoly )

feelpp_add_application( geomap
    SRCS 00_geomap.cpp
    LINK_LIBRARIES benchmark
    EXEC B_FEELPOLY_GEOMAP )

feelpp_add_application( geomap_inverse
    SRCS 01_geomap_inverse.cpp
    LINK_LIBRARIES benchmark
    EXEC B_FEELPOLY_GEOMAP_INVERSE )

add_dependencies( benchs-feelpoly ${B_FEELPOLY_GEOMAP} )

add_dependencies( benchs-feelpoly
    ${B_FEELPOLY_GEOMAP}
    ${B_FEELPOLY_GEOMAP_INVERSE}
    ${B_FEELPOLY_GEOMESH_STORAGE} )
