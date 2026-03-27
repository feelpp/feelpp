/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme
       Date: 2026-03-25

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#define BOOST_TEST_MODULE trial_test_basis testsuite
#include <feel/feelcore/testsuite.hpp>

#include <cmath>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

template <typename FormType, typename TestElement, typename TrialElement>
double
formEnergy( FormType& form,
            TestElement const& testElement,
            TrialElement const& trialElement )
{
    form.close();
    return form.matrixPtr()->energy( testElement, trialElement );
}

} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( trial_test_basis_suite )

BOOST_AUTO_TEST_CASE( scalar_basis_proxy_matches_legacy_syntax )
{
    auto mesh = unitSquare();
    auto Xh = Pch<2>( mesh );

    auto uLegacy = Xh->element( "u_legacy" );
    auto vLegacy = Xh->element( "v_legacy" );
    auto u = trial( Xh, "u_phase1" );
    auto v = test( Xh, "v_phase1" );

    auto legacy = form2( _test=Xh, _trial=Xh );
    legacy = integrate( _range=elements( mesh ),
                        _expr=inner( gradt( uLegacy ), grad( vLegacy ) ) + idt( uLegacy )*id( vLegacy ) );

    auto modern = form2( _test=Xh, _trial=Xh );
    modern = integrate( _range=elements( mesh ),
                        _expr=inner( grad( u ), grad( v ) ) + u*v );

    auto trialProbe = Xh->element( Px() + 2.0*Py(), "trial_probe" );
    auto testProbe = Xh->element( cst( 1.0 ) + Px()*Py(), "test_probe" );

    double legacyEnergy = formEnergy( legacy, testProbe, trialProbe );
    double modernEnergy = formEnergy( modern, testProbe, trialProbe );

    BOOST_CHECK_SMALL( std::abs( legacyEnergy-modernEnergy ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( vector_basis_proxy_matches_legacy_symm_grad )
{
    auto mesh = unitSquare();
    auto Xh = Pchv<2>( mesh );

    auto uLegacy = Xh->element( "u_legacy" );
    auto vLegacy = Xh->element( "v_legacy" );
    auto u = trial( uLegacy );
    auto v = test( vLegacy );

    auto legacy = form2( _test=Xh, _trial=Xh );
    legacy = integrate( _range=elements( mesh ),
                        _expr=inner( sym( gradt( uLegacy ) ), sym( grad( vLegacy ) ) ) + inner( idt( uLegacy ), id( vLegacy ) ) );

    auto modern = form2( _test=Xh, _trial=Xh );
    modern = integrate( _range=elements( mesh ),
                        _expr=inner( symm_grad( u ), symm_grad( v ) ) + inner( u, v ) );

    auto trialProbe = Xh->element( vec( Px()+Py(), 2.0*Px()-Py() ), "trial_probe" );
    auto testProbe = Xh->element( vec( cst( 1.0 )+Px(), Px()*Py() ), "test_probe" );

    double legacyEnergy = formEnergy( legacy, testProbe, trialProbe );
    double modernEnergy = formEnergy( modern, testProbe, trialProbe );

    BOOST_CHECK_SMALL( std::abs( legacyEnergy-modernEnergy ), 1e-12 );
}

BOOST_AUTO_TEST_SUITE_END()
