/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include <feel/feelcore/info.hpp>
#include <iostream>

#ifdef FEELPP_KOKKOS_DEPENDENCY_SENTINEL
#error "Kokkos usage requirements reached the ordinary CPU target"
#endif

int main()
{
    std::cout << Feel::Info::versionString() << '\n';
}
