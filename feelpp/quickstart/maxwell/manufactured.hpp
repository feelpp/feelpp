/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_QUICKSTART_MAXWELL_MANUFACTURED_HPP
#define FEELPP_QUICKSTART_MAXWELL_MANUFACTURED_HPP

#include <feel/feelvf/ginac.hpp>

namespace Feel::quickstart::maxwell
{
// Shared datasets for the form-based and Problem API Maxwell quickstarts.
template <int Dim>
auto exactField( bool polynomial )
{
    static_assert( Dim == 2 || Dim == 3 );
    if constexpr ( Dim == 2 )
        return vf::expr<2, 1>( polynomial ? "{-y,x}:x:y" :
                                      "{sin(pi*y),sin(pi*x)}:x:y" );
    else
        return vf::expr<3, 1>( polynomial ? "{-y,x,0}:x:y:z" :
                                      "{sin(pi*y),sin(pi*z),sin(pi*x)}:x:y:z" );
}

template <int Dim>
auto exactCurl( bool polynomial )
{
    static_assert( Dim == 2 || Dim == 3 );
    if constexpr ( Dim == 2 )
        return vf::expr( polynomial ? "2" : "pi*(cos(pi*x)-cos(pi*y)):x:y" );
    else
        return vf::expr<3, 1>( polynomial ? "{0,0,2}" :
                                      "{-pi*cos(pi*z),-pi*cos(pi*x),-pi*cos(pi*y)}:x:y:z" );
}
} // namespace Feel::quickstart::maxwell

#endif
