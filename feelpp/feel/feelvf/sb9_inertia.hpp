/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_VF_SB9_INERTIA_HPP
#define FEELPP_VF_SB9_INERTIA_HPP 1

#include <feel/feelvf/basis.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/shellgeometric.hpp>

namespace Feel::vf
{
/**
 * Cell-local displacement enrichment per unit internal SB9 amplitude.
 *
 * The normal belongs to the reference shell cell. On an affine, orthogonal
 * thickness map, differentiating (1-zeta^2)*alpha along the normal gives the
 * existing W9 strain -4*zeta*alpha/h. Alpha is an enrichment amplitude, not
 * the absolute displacement of an additional mesh node.
 */
[[nodiscard]] inline auto
sb9DisplacementBubble()
{
    auto zt = zeta();
    return ( cst( 1.0 ) - zt * zt ) * shellNormal();
}

/** Assemble the Newmark inertia contribution to the effective mixed operator.
 *
 * The consistent model uses u* = u + (1-zeta^2)*alpha*n. Its mass matrix has
 * four blocks. The massless model retains only the physical Q1 mass.
 * Full volume quadrature is required: vertex Lobatto quadrature annihilates
 * the internal bubble. Order six also integrates the quartic bubble product
 * on affine cells without reusing the reduced stiffness rule.
 */
template <typename BilinearForm, RangeConcept Range,
          detail::BasisProxyType U, detail::BasisProxyType V,
          detail::BasisProxyType Alpha, detail::BasisProxyType Beta>
void
assembleSb9Inertia( BilinearForm& a, Range const& range,
                    U const& u, V const& v, Alpha const& alpha, Beta const& beta,
                    double density, double accelerationCoefficient,
                    bool consistentInternalInertia )
{
    double const scale = density * accelerationCoefficient;
    auto bubble = sb9DisplacementBubble();
    a( 0_c, 0_c ) += integrate( _range=range, _quad=_Q<6>(), _quad1=_Q<6>(),
                                _expr=scale * inner( id( u ), id( v ) ) );
    if ( consistentInternalInertia )
    {
        a( 0_c, 1_c ) += integrate( _range=range, _quad=_Q<6>(), _quad1=_Q<6>(),
                                    _expr=scale * inner( bubble * id( alpha ), id( v ) ) );
        a( 1_c, 0_c ) += integrate( _range=range, _quad=_Q<6>(), _quad1=_Q<6>(),
                                    _expr=scale * inner( id( u ), bubble * id( beta ) ) );
        a( 1_c, 1_c ) += integrate( _range=range, _quad=_Q<6>(), _quad1=_Q<6>(),
                                    _expr=scale * inner( bubble * id( alpha ), bubble * id( beta ) ) );
    }
}

/** Assemble M times the Newmark acceleration-history polynomial in both rows.
 * The caller clears the complete mixed RHS before assembling a new step.
 */
template <typename LinearForm, RangeConcept Range,
          detail::BasisProxyType V, detail::BasisProxyType Beta,
          typename DisplacementHistory, typename InternalHistory>
void
assembleSb9InertiaHistory( LinearForm& l, Range const& range,
                           V const& v, Beta const& beta,
                           DisplacementHistory const& displacementHistory,
                           InternalHistory const& internalHistory,
                           double density, bool consistentInternalInertia )
{
    auto bubble = sb9DisplacementBubble();
    l( 0_c ) += integrate( _range=range, _quad=_Q<6>(), _quad1=_Q<6>(),
                            _expr=density * inner( idv( displacementHistory ), id( v ) ) );
    if ( consistentInternalInertia )
    {
        auto history = idv( displacementHistory ) + bubble * idv( internalHistory );
        l( 0_c ) += integrate( _range=range, _quad=_Q<6>(), _quad1=_Q<6>(),
                                _expr=density * inner( bubble * idv( internalHistory ), id( v ) ) );
        l( 1_c ) += integrate( _range=range, _quad=_Q<6>(), _quad1=_Q<6>(),
                                _expr=density * inner( history, bubble * id( beta ) ) );
    }
}
} // namespace Feel::vf

#endif
