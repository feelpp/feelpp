/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#pragma once

#include <algorithm>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel::quickstart
{
// A scalar DG0 mass solve computes physical cell means, including the varying
// geometric Jacobian on curved cells. Each local system has one unknown.
template<typename Element, typename ConstantSpacePtr>
auto cellMean( Element const& field, ConstantSpacePtr const& constants )
{
    using namespace boost::hana::literals;
    CHECK_EQ( constants->fe()->order(), 0 ) << "cellMean requires a DG0 space";
    auto space = product( constants );
    auto result = space.element();
    auto mean = result( 0_c );
    auto mass = blockform2( space, solve::strategy::local, backend() );
    auto rhs = blockform1( space, solve::strategy::local, backend() );
    auto range = elements( constants->mesh() );
    // A full-dimensional simplex has a polynomial Jacobian determinant of
    // degree dim*(g-1). Use the same rule for both physical cell integrals.
    int const quadratureOrder = std::max( 1, static_cast<int>( field.functionSpace()->fe()->order() ) +
        constants->mesh()->nDim * ( static_cast<int>( constants->mesh()->order() ) - 1 ) );
    mass( 0_c, 0_c ) = integrate( _range=range, _expr=idt(mean)*id(mean), _quad=quadratureOrder );
    rhs( 0_c ) = integrate( _range=range, _expr=idv(field)*id(mean), _quad=quadratureOrder );
    mass.solve( _solution=result, _rhs=rhs, _local=true );
    return result( 0_c );
}

template<typename Element, typename ReferenceElement, typename ConstantSpacePtr>
void matchCellMean( Element& field, ReferenceElement const& reference, ConstantSpacePtr const& constants )
{
    auto fieldMean = cellMean( field, constants );
    auto referenceMean = cellMean( reference, constants );
    // Different runtime orders share an Element type: coefficient arithmetic
    // would not lift the DG0 correction to every local interpolation node.
    auto corrected = field.functionSpace()->element();
    corrected.on( _range=elements( constants->mesh() ),
                  _expr=idv(field) - idv(fieldMean) + idv(referenceMean) );
    field = corrected;
}
}
