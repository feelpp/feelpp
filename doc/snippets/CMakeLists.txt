###  TEMPLATE.txt.tpl; coding: utf-8 ---

# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileContributor: Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>
#
# SPDX-FileCopyrightText: 2011 Université Joseph Fourier (Grenoble I)
# SPDX-FileCopyrightText: 2011-2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
feelpp_add_application( dofpoints SRCS dofpoints.cpp LABELS doc)
feelpp_add_application( meshho SRCS meshho.cpp LABELS doc)
feelpp_add_application( bdfpod SRCS bdfpod.cpp LABELS doc)
feelpp_add_application( pmesh SRCS pmesh.cpp CFG pmesh.cfg LABELS doc)
feelpp_add_application( dofneighbors SRCS dofneighbors.cpp LABELS doc)
feelpp_add_application( dofboundary SRCS dofboundary.cpp LABELS doc)
feelpp_add_application( forms SRCS forms.cpp LABELS doc)
