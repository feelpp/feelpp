# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2018-2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
feel_append_src(feelsystem "${HDR}")
