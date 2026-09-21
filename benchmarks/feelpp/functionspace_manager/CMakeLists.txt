# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

feelpp_add_application(
  functionspace_manager_baseline
  SRCS functionspace_manager_baseline.cpp
  LINK_LIBRARIES benchmark
  EXEC B_FSPACE_MANAGER_BASELINE
)

feelpp_add_application(
  functionspace_manager
  SRCS functionspace_manager.cpp
  LINK_LIBRARIES benchmark
  EXEC B_FSPACE_MANAGER
)
