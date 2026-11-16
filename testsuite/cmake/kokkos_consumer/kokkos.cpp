/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <feel/feelconfig.h>
#include <feel/feelcore/info.hpp>
#include <feel/feelcore/kokkosruntime.hpp>
#include <iostream>

#if !defined( FEELPP_HAS_KOKKOS ) || !defined( FEELPP_KOKKOS_DEPENDENCY_SENTINEL )
#error "Kokkos integration usage requirements are missing"
#endif

/** Exercise the installed runtime API without initializing MPI or Kokkos. */
int main()
{
    std::cout << Feel::Info::versionString() << " Kokkos " << KOKKOS_VERSION
              << " Kokkos Kernels " << KOKKOSKERNELS_VERSION << '\n';
    Feel::KokkosRuntime runtime( Feel::KokkosInitPolicy::off );
    return runtime.status().initialized || !Feel::KokkosRuntime::capabilities().available ? 1 : 0;
}
