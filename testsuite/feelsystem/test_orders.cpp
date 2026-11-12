/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE feelsystem_orders
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelsystem/system.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template <typename SpacePtr>
void checkProblem( SpacePtr const& Xh )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto const mesh = Xh->mesh();
    auto uh = Xh->element();
    auto exact = cst( 1.0 ) + Px() + Py();

    // Store expressions built from temporary named bases. Assembly happens
    // after every basis temporary has been destroyed.
    auto problem = ellipticProblem( _space = Xh )
        .lhs( integrate( _range = elements( mesh ),
                         _expr = inner( grad( trial( Xh, "u" ) ), grad( test( Xh, "v" ) ) )
                               + trial( Xh ) * test( Xh ) ) )
        .rhs( integrate( _range = elements( mesh ), _expr = exact * test( Xh ) ) )
        .dirichlet( _range = boundaryfaces( mesh ), _element = uh, _expr = exact );

    static_assert( systems::LinearProblemLike<decltype( problem )> );
    static_assert( Feel::vf::detail::BasisProxyType<decltype( trial( Xh ) )> );
    static_assert( Feel::vf::detail::basis_proxy_type_t<decltype( trial( Xh ) )>::role == __TRIAL );

    auto system = problem.assemble();
    auto result = system.solve( _solution = uh );
    BOOST_REQUIRE( result.converged );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ), _expr = idv( uh ) - exact ), 1e-8 );

    // Replay the descriptor after the first system has been solved.
    auto replay = problem.assemble();
    auto matrixDelta = replay.matrixPtr()->clone();
    matrixDelta->addMatrix( -1.0, system.matrixPtr() );
    matrixDelta->close();
    BOOST_CHECK_SMALL( matrixDelta->linftyNorm(), 1e-10 );
}
} // namespace

BOOST_AUTO_TEST_CASE( static_and_dynamic_orders_support_stored_problems )
{
    using namespace Feel;
    auto mesh = unitSquare( 0.4 );
    checkProblem( Pch<1>( mesh ) );
    checkProblem( Pch<2>( mesh ) );
    for ( uint16_type order : { 1, 2 } )
        checkProblem( Pch<Dynamic>( mesh, RuntimeOrder{ order } ) );
}
