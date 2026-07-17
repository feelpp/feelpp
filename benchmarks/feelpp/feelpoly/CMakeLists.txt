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

add_dependencies( benchs-feelpoly ${B_FEELPOLY_GEOMAP} )
