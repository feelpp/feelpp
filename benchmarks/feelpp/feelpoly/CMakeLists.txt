# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

add_custom_target( benchs-feelpoly )

# Baseline measurements are opt-in and carry no timing assertions in CTest.
feelpp_add_application( order_baseline
  SRCS order_baseline.cpp allocation_probe.cpp
  LINK_LIBRARIES benchmark
  NO_TEST EXCLUDE_FROM_ALL EXEC B_ORDER_BASELINE )

feelpp_add_application( geomap
    SRCS 00_geomap.cpp
    LINK_LIBRARIES benchmark
    EXEC B_FEELPOLY_GEOMAP )

feelpp_add_application( geomap_inverse
    SRCS 01_geomap_inverse.cpp
    LINK_LIBRARIES benchmark
    EXEC B_FEELPOLY_GEOMAP_INVERSE )

feelpp_add_application( geomesh_storage
    SRCS 02_geomesh_storage.cpp allocation_probe.cpp
    LINK_LIBRARIES benchmark
    EXEC B_FEELPOLY_GEOMESH_STORAGE )

add_dependencies( benchs-feelpoly
    ${B_ORDER_BASELINE}
    ${B_FEELPOLY_GEOMAP}
    ${B_FEELPOLY_GEOMAP_INVERSE}
    ${B_FEELPOLY_GEOMESH_STORAGE} )
