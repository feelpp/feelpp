# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later

option(FEELPP_ENABLE_KOKKOS "Enable optional Kokkos integration targets" OFF)

set(FEELPP_HAS_KOKKOS 0)
set(FEELPP_KOKKOS_VERSION "")
set(FEELPP_KOKKOS_KERNELS_VERSION "")
if(FEELPP_ENABLE_KOKKOS)
  find_package(Kokkos 5.1 CONFIG REQUIRED)
  find_package(KokkosKernels 5.1 CONFIG REQUIRED)
  foreach(_feelpp_kokkos_target Kokkos::kokkos Kokkos::kokkoskernels)
    if(NOT TARGET ${_feelpp_kokkos_target})
      message(FATAL_ERROR "Missing Kokkos integration target ${_feelpp_kokkos_target}")
    endif()
  endforeach()
  set(FEELPP_HAS_KOKKOS 1)
  set(FEELPP_KOKKOS_VERSION "${Kokkos_VERSION}")
  set(FEELPP_KOKKOS_KERNELS_VERSION "${KokkosKernels_VERSION}")
  list(APPEND FEELPP_ENABLED_OPTIONS Kokkos)
endif()
