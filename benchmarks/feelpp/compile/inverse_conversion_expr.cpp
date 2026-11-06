/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include <numbers>

#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

template <int Dim>
auto
makeMesh()
{
    if constexpr ( Dim == 2 )
        return unitSquare();
    else
        return unitCube();
}

template <int Dim>
int
instantiate()
{
    auto mesh = makeMesh<Dim>();

    if constexpr ( Dim == 2 )
    {
        auto v = vec( cst( 1.0 ), Px(), cst( 2.0 ) + Py() );
        auto m = vec( cst( 1.0 ), cst( std::numbers::sqrt2_v<double> )*Px(), cst( 2.0 ) + Py() );
        auto integral = integrate( _range=elements( mesh ),
                                   _expr=inner( unvoigt( v ), unvoigt( v ) ) +
                                         inner( unmandel( m ), unmandel( m ) ) );
        return sizeof( integral ) > 0 ? 0 : 0;
    }
    else
    {
        auto v = vec( cst( 1.0 ), Px(), Py(),
                      cst( 2.0 ) + Px(), cst( 3.0 ) + Pz(), cst( 4.0 ) + Px() );
        auto m = vec( cst( 1.0 ), cst( std::numbers::sqrt2_v<double> )*Px(), cst( std::numbers::sqrt2_v<double> )*Py(),
                      cst( 2.0 ) + Px(), cst( std::numbers::sqrt2_v<double> )*( cst( 3.0 ) + Pz() ), cst( 4.0 ) + Px() );
        auto integral = integrate( _range=elements( mesh ),
                                   _expr=inner( unvoigt( v ), unvoigt( v ) ) +
                                         inner( unmandel( m ), unmandel( m ) ) );
        return sizeof( integral ) > 0 ? 0 : 0;
    }
}

} // namespace

int
main()
{
    return instantiate<2>() + instantiate<3>();
}
