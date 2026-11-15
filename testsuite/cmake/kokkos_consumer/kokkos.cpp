/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include <feel/feelcore/info.hpp>
#include <feel/feelconfig.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_config.h>
#include <iostream>

#if !defined(FEELPP_HAS_KOKKOS) || !defined(FEELPP_KOKKOS_DEPENDENCY_SENTINEL)
#error "Kokkos integration usage requirements are missing"
#endif

int main()
{
    std::cout << Feel::Info::versionString() << " Kokkos " << KOKKOS_VERSION
              << " Kokkos Kernels " << KOKKOSKERNELS_VERSION << '\n';
    return Kokkos::is_initialized() ? 1 : 0;
}
