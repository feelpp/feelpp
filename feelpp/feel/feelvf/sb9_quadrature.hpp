/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
 * @file sb9_quadrature.hpp
 * @brief SB9-specific quadrature rules used by the shell element formulation.
 */

#ifndef FEELPP_VF_SB9_QUADRATURE_HPP
#define FEELPP_VF_SB9_QUADRATURE_HPP 1

#include <array>
#include <cmath>

#include <feel/feelpoly/im.hpp>

namespace Feel
{
namespace vf
{
/**
 * @brief Reduced Gauss-Lobatto quadrature through the SB9 shell thickness.
 *
 * This rule is tailored to the SB9 stiffness terms. It places all points at
 * the in-plane reference center, \f$(\xi,\eta)=(0,0)\f$, and uses the
 * five-point Gauss-Lobatto rule in the thickness coordinate \f$\zeta\f$:
 * \f[
 *   \zeta = \left(-1,-\sqrt{3/7},0,\sqrt{3/7},1\right).
 * \f]
 *
 * The one-dimensional weights are scaled by the in-plane reference measure
 * \f$4\f$, so the total weight is \f$8\f$, the measure of the reference
 * hexahedron \f$[-1,1]^3\f$.
 *
 * @tparam T Floating-point value type used for quadrature points and weights.
 */
template <typename T = double>
class SB9ThroughThicknessLobatto5
    : public IMGeneral<3, T, Hypercube>
{
public:
    /// Base integration method type.
    using super = IMGeneral<3, T, Hypercube>;
    /// Point-set quadrature type associated with a trilinear reference hexahedron.
    using pointset_super = PointSetQuadrature<Hypercube<3, 1, 3>, T>;
    /// Scalar value type for coordinates and weights.
    using value_type = T;
    /// Point-set return type used by Feel++ quadrature infrastructure.
    using return_type = typename pointset_super::return_type;
    /// Matrix type storing quadrature point coordinates.
    using nodes_type = typename super::nodes_type;
    /// Vector type storing quadrature weights.
    using weights_type = typename super::weights_type;

    /// Number of quadrature points in the rule.
    static inline const uint32_type Npoints = 5;

    /**
     * @brief Build the five-point through-thickness Lobatto rule.
     *
     * The points are ordered from the lower thickness face to the upper
     * thickness face. Each point has \f$\xi=\eta=0\f$.
     */
    SB9ThroughThicknessLobatto5()
        : super()
    {
        this->M_order = 5;
        this->M_name = "sb9-through-thickness-lobatto5";
        this->M_w.resize( Npoints );
        this->M_prod.resize( Npoints );
        this->M_exprq.resize( Npoints );
        this->M_w_sum = value_type( 0 );

        setWeight( 0, value_type( 4 ) * value_type( 1 ) / value_type( 10 ) );
        setWeight( 1, value_type( 4 ) * value_type( 49 ) / value_type( 90 ) );
        setWeight( 2, value_type( 4 ) * value_type( 32 ) / value_type( 45 ) );
        setWeight( 3, value_type( 4 ) * value_type( 49 ) / value_type( 90 ) );
        setWeight( 4, value_type( 4 ) * value_type( 1 ) / value_type( 10 ) );

        value_type const a = std::sqrt( value_type( 3 ) / value_type( 7 ) );
        nodes_type points( 3, Npoints );

        points( 0, 0 ) = value_type( 0 );
        points( 1, 0 ) = value_type( 0 );
        points( 2, 0 ) = value_type( -1 );

        points( 0, 1 ) = value_type( 0 );
        points( 1, 1 ) = value_type( 0 );
        points( 2, 1 ) = -a;

        points( 0, 2 ) = value_type( 0 );
        points( 1, 2 ) = value_type( 0 );
        points( 2, 2 ) = value_type( 0 );

        points( 0, 3 ) = value_type( 0 );
        points( 1, 3 ) = value_type( 0 );
        points( 2, 3 ) = a;

        points( 0, 4 ) = value_type( 0 );
        points( 1, 4 ) = value_type( 0 );
        points( 2, 4 ) = value_type( 1 );

        this->setPoints( points );
    }

    FEELPP_DEFINE_VISITABLE();

private:
    /**
     * @brief Store a quadrature weight and update the cached total weight.
     *
     * @param q Quadrature point index.
     * @param weight Weight associated with point @p q.
     */
    void setWeight( uint16_type q, value_type weight )
    {
        this->M_w( q ) = weight;
        this->M_w_sum += weight;
    }
};

/**
 * @brief Create the reduced five-point SB9 stiffness quadrature rule.
 *
 * Use this rule for the SB9 membrane-bending, shear, pinching and internal
 * mode stiffness blocks that follow the MATLAB/documentation reduced
 * through-thickness integration.
 *
 * @par Example
 * @code{.cpp}
 * auto sb9Quad = Feel::vf::sb9ThroughThicknessLobatto5();
 *
 * a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
 *                             _quad=sb9Quad,
 *                             _quad1=sb9Quad,
 *                             _expr=ddot( C, epsShellTrial, epsShellTest ) );
 * @endcode
 *
 * @tparam T Floating-point value type used for quadrature points and weights.
 * @return A ready-to-use SB9ThroughThicknessLobatto5 quadrature object.
 */
template <typename T = double>
[[nodiscard]] inline SB9ThroughThicknessLobatto5<T>
sb9ThroughThicknessLobatto5()
{
    return SB9ThroughThicknessLobatto5<T>{};
}

/**
 * @brief Full vertex Gauss-Lobatto quadrature for lumped SB9 displacement mass.
 *
 * This rule is intended for the physical Q1 displacement mass matrix only. It
 * uses the tensor-product two-point Lobatto rule in each reference direction,
 * i.e. the eight vertices of \f$[-1,1]^3\f$, each with weight \f$1\f$.
 *
 * For Q1 Lagrange displacement basis functions on a hexahedron, evaluating the
 * mass matrix at the interpolation vertices diagonalizes the consistent mass
 * contribution. The internal SB9 scalar mode is not a physical displacement
 * degree of freedom and should not receive an inertia contribution from this
 * rule.
 *
 * @tparam T Floating-point value type used for quadrature points and weights.
 */
template <typename T = double>
class SB9LumpedMassLobatto
    : public IMGeneral<3, T, Hypercube>
{
public:
    /// Base integration method type.
    using super = IMGeneral<3, T, Hypercube>;
    /// Point-set quadrature type associated with a trilinear reference hexahedron.
    using pointset_super = PointSetQuadrature<Hypercube<3, 1, 3>, T>;
    /// Scalar value type for coordinates and weights.
    using value_type = T;
    /// Point-set return type used by Feel++ quadrature infrastructure.
    using return_type = typename pointset_super::return_type;
    /// Matrix type storing quadrature point coordinates.
    using nodes_type = typename super::nodes_type;

    /// Number of quadrature points in the rule.
    static inline const uint32_type Npoints = 8;

    /**
     * @brief Build the eight-vertex tensor-product Lobatto rule.
     *
     * The rule visits all sign combinations
     * \f$(\xi,\eta,\zeta)\in\{-1,1\}^3\f$ and assigns weight \f$1\f$ to each
     * point, yielding total reference weight \f$8\f$.
     */
    SB9LumpedMassLobatto()
        : super()
    {
        this->M_order = 1;
        this->M_name = "sb9-lumped-mass-lobatto";
        this->M_w.resize( Npoints );
        this->M_prod.resize( Npoints );
        this->M_exprq.resize( Npoints );
        this->M_w_sum = value_type( 0 );

        nodes_type points( 3, Npoints );
        std::array<value_type, 2> const lobattoPoints = {{ value_type( -1 ), value_type( 1 ) }};
        uint16_type q = 0;
        for ( value_type z : lobattoPoints )
        {
            for ( value_type y : lobattoPoints )
            {
                for ( value_type x : lobattoPoints )
                {
                    points( 0, q ) = x;
                    points( 1, q ) = y;
                    points( 2, q ) = z;
                    setWeight( q, value_type( 1 ) );
                    ++q;
                }
            }
        }

        this->setPoints( points );
    }

    FEELPP_DEFINE_VISITABLE();

private:
    /**
     * @brief Store a quadrature weight and update the cached total weight.
     *
     * @param q Quadrature point index.
     * @param weight Weight associated with point @p q.
     */
    void setWeight( uint16_type q, value_type weight )
    {
        this->M_w( q ) = weight;
        this->M_w_sum += weight;
    }
};

/**
 * @brief Create the eight-vertex SB9 lumped mass quadrature rule.
 *
 * Use this rule to assemble the mass matrix on the Q1 displacement space. Do
 * not use it for SB9 stiffness terms, and do not assemble it on the internal
 * SB9 scalar mode.
 *
 * @par Example
 * @code{.cpp}
 * auto Uh = Pchv<1>( mesh );
 * auto u = trial( Uh, "u" );
 * auto v = test( Uh, "v" );
 * auto massQuad = Feel::vf::sb9LumpedMassLobatto();
 *
 * auto m = form2( _trial=Uh, _test=Uh );
 * m = integrate( _range=elements( mesh ),
 *                _quad=massQuad,
 *                _quad1=massQuad,
 *                _expr=cst( rho ) * inner( u, v ) );
 * @endcode
 *
 * @tparam T Floating-point value type used for quadrature points and weights.
 * @return A ready-to-use SB9LumpedMassLobatto quadrature object.
 */
template <typename T = double>
[[nodiscard]] inline SB9LumpedMassLobatto<T>
sb9LumpedMassLobatto()
{
    return SB9LumpedMassLobatto<T>{};
}
} // namespace vf
} // namespace Feel

#endif // FEELPP_VF_SB9_QUADRATURE_HPP
