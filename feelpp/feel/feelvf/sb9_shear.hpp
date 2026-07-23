/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
   \file sb9_shear.hpp
   \brief SB9 shear shell kinematic operators
 */
#ifndef FEELPP_VF_SB9_SHEAR_HPP
#define FEELPP_VF_SB9_SHEAR_HPP 1

#include <feel/feelvf/sb9_common.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{
/**
 * \brief Selects an SB9 transverse-shear coefficient block.
 *
 * `Bc0` is the physical transverse-shear operator used in the SB9 elastic
 * strain. `Bc1` and `Bc2` are the two Hallquist transverse-shear stabilization
 * blocks used by the stabilized formulation.
 */
enum class SB9ShearKind
{
    /// Physical transverse-shear block contributing to eps_xz and eps_yz.
    Bc0,
    /// First transverse-shear stabilization block.
    Bc1,
    /// Second transverse-shear stabilization block.
    Bc2
};

/**
 * \brief Element-local cache for SB9 transverse-shear coefficients.
 *
 * The cache computes the geometry-dependent auxiliary rows for the three SB9
 * transverse-shear blocks once per element and reuses them while Feel++
 * evaluates the trial/test basis coefficients. Public users normally access
 * this through the free functions \ref sb9Bc0, \ref sb9Bc1, \ref sb9Bc2 and
 * \ref sb9Shear.
 *
 * \tparam GeometryDataType Shell cell geometry cache type provided by
 *         ShellCellGeometryTensorBase.
 */
template <typename GeometryDataType>
class SB9ShearKernelCache : public SB9KernelBase<GeometryDataType>
{
public:
    /// Shared SB9 geometry cache base.
    using base_type = SB9KernelBase<GeometryDataType>;
    /// Scalar type used by the geometry data and generated coefficients.
    using value_type = typename base_type::value_type;
    /// Three-component local frame row type.
    using row_type = typename base_type::row_type;
    /// Pair of auxiliary rows defining one transverse-shear block.
    using shear_aux_type = std::array<row_type, 2>;
    /// Number of geometric nodes in the SB9 Q1 hexahedral element.
    static constexpr uint16_type node_count = base_type::node_count;
    /// Number of displacement components handled by the vector operator.
    static constexpr uint16_type component_count = base_type::component_count;

    /**
     * \brief Build the element-local transverse-shear coefficient cache.
     *
     * \param data Element-local shell geometry data.
     */
    explicit SB9ShearKernelCache( GeometryDataType const& data )
        :
        base_type( data )
    {
        for ( uint16_type node = 0; node < node_count; ++node )
        {
            M_shear[0][node] = this->computeShearAux( node, SB9ShearKind::Bc0 );
            M_shear[1][node] = this->computeShearAux( node, SB9ShearKind::Bc1 );
            M_shear[2][node] = this->computeShearAux( node, SB9ShearKind::Bc2 );
        }
    }

    /**
     * \brief Fill one SB9 shear coefficient column.
     *
     * The raw `Bc*` operator uses the local SB9 vector-operator convention:
     * only slots 4 and 5 are written. The public \ref sb9Shear helper maps
     * those two slots into the compact symmetric Mandel strain vector used by
     * the assembled elastic formulation.
     *
     * \tparam Kind Compile-time selector for `Bc0`, `Bc1`, or `Bc2`.
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector.
     * \param node Local geometric node index.
     * \param component Displacement component index.
     */
    template <SB9ShearKind Kind, typename VectorType>
    void fillVectorCoefficients( VectorType& coeff, uint16_type node, uint16_type component ) const
    {
        static_assert( Kind == SB9ShearKind::Bc0 || Kind == SB9ShearKind::Bc1 || Kind == SB9ShearKind::Bc2,
                       "unsupported SB9 shear vector kind" );

        auto const& frame = this->M_frame[component];
        auto const& aux = M_shear[shearKindIndex<Kind>()][node];

        coeff( 4 ) = this->mandelShearScale() * this->dot( aux[0], frame );
        coeff( 5 ) = this->mandelShearScale() * this->dot( aux[1], frame );
    }

private:
    /**
     * \brief Map the compile-time shear-kind tag to the cache array index.
     *
     * \tparam Kind Compile-time selector for `Bc0`, `Bc1`, or `Bc2`.
     * \return Cache index in `M_shear`.
     */
    template <SB9ShearKind Kind>
    static consteval int shearKindIndex()
    {
        if constexpr ( Kind == SB9ShearKind::Bc0 )
            return 0;
        else if constexpr ( Kind == SB9ShearKind::Bc1 )
            return 1;
        else
            return 2;
    }

    /**
     * \brief Compute the two auxiliary rows for one node and one shear block.
     *
     * The returned rows are later projected on the local displacement frame.
     * The formulas follow the SB9 documentation/MATLAB decomposition of the
     * transverse-shear blocks and keep the geometry-dependent terms separated
     * from the trial/test component projection.
     *
     * \param node Local geometric node index.
     * \param kind Runtime selector used while filling the element cache.
     * \return Two auxiliary rows associated with the selected shear block.
     */
    shear_aux_type computeShearAux( uint16_type node, SB9ShearKind kind ) const
    {
        auto const& Ja = this->M_data.Ja;
        auto const& Jb = this->M_data.Jb;
        auto const& Jc = this->M_data.Jc;
        auto const& Jd = this->M_data.Jd;
        auto const& invJa = this->M_data.invJa;
        auto const& invJb = this->M_data.invJb;
        auto const& invJc = this->M_data.invJc;
        auto const& invJd = this->M_data.invJd;
        auto const& Bx = this->M_data.bx;
        auto const& By = this->M_data.by;
        auto const& Bz = this->M_data.bz;
        auto const& Vgamma = this->M_data.vgamma;

        value_type f11bz = value_type( 0 );
        value_type f11g1 = value_type( 0 );
        value_type f12bz = value_type( 0 );
        value_type f12g1 = value_type( 0 );
        value_type f13bx = value_type( 0 );
        value_type f13by = value_type( 0 );
        value_type f13g3 = value_type( 0 );
        value_type f21bz = value_type( 0 );
        value_type f21g1 = value_type( 0 );
        value_type f21g2 = value_type( 0 );
        value_type f22bz = value_type( 0 );
        value_type f22g1 = value_type( 0 );
        value_type f22g2 = value_type( 0 );
        value_type f23bx = value_type( 0 );
        value_type f23by = value_type( 0 );
        value_type f23g3 = value_type( 0 );

        if ( kind == SB9ShearKind::Bc0 )
        {
            f11bz = this->M_invJ0( 0, 0 )*( Ja( 0, 0 ) + Jc( 0, 0 ) )/value_type( 2 );
            f11g1 = this->M_invJ0( 0, 0 )*( invJc( 2, 2 )*Jc( 0, 0 ) - invJa( 2, 2 )*Ja( 0, 0 ) )/value_type( 2 );

            f12bz = this->M_invJ0( 0, 0 )*( Ja( 0, 1 ) + Jc( 0, 1 ) )/value_type( 2 );
            f12g1 = this->M_invJ0( 0, 0 )*( invJc( 2, 2 )*Jc( 0, 1 ) - invJa( 2, 2 )*Ja( 0, 1 ) )/value_type( 2 );

            f13bx = f11bz;
            f13by = f12bz;
            f13g3 = this->M_invJ0( 0, 0 )*( invJc( 0, 0 )*Jc( 0, 0 ) + invJc( 1, 0 )*Jc( 0, 1 ) -
                                            invJa( 0, 0 )*Ja( 0, 0 ) - invJa( 1, 0 )*Ja( 0, 1 ) )/value_type( 2 );

            f21bz = this->M_invJ0( 1, 0 )*( Ja( 0, 0 ) + Jc( 0, 0 ) )/value_type( 2 ) +
                    this->M_invJ0( 1, 1 )*( Jb( 1, 0 ) + Jd( 1, 0 ) )/value_type( 2 );
            f21g1 = this->M_invJ0( 1, 0 )*( invJc( 2, 2 )*Jc( 0, 0 ) - invJa( 2, 2 )*Ja( 0, 0 ) )/value_type( 2 );
            f21g2 = this->M_invJ0( 1, 1 )*( invJb( 2, 2 )*Jb( 1, 0 ) - invJd( 2, 2 )*Jd( 1, 0 ) )/value_type( 2 );

            f22bz = this->M_invJ0( 1, 0 )*( Ja( 0, 1 ) + Jc( 0, 1 ) )/value_type( 2 ) +
                    this->M_invJ0( 1, 1 )*( Jb( 1, 1 ) + Jd( 1, 1 ) )/value_type( 2 );
            f22g1 = this->M_invJ0( 1, 0 )*( invJc( 2, 2 )*Jc( 0, 1 ) - invJa( 2, 2 )*Ja( 0, 1 ) )/value_type( 2 );
            f22g2 = this->M_invJ0( 1, 1 )*( invJb( 2, 2 )*Jb( 1, 1 ) - invJd( 2, 2 )*Jd( 1, 1 ) )/value_type( 2 );

            f23bx = f21bz;
            f23by = f22bz;
            f23g3 = this->M_invJ0( 1, 0 )*( invJc( 0, 0 )*Jc( 0, 0 ) + invJc( 1, 0 )*Jc( 0, 1 ) -
                                            invJa( 0, 0 )*Ja( 0, 0 ) - invJa( 1, 0 )*Ja( 0, 1 ) )/value_type( 2 );
            f23g3 += this->M_invJ0( 1, 1 )*( invJb( 0, 1 )*Jb( 1, 0 ) + invJb( 1, 1 )*Jb( 1, 1 ) -
                                             invJd( 0, 1 )*Jd( 1, 0 ) - invJd( 1, 1 )*Jd( 1, 1 ) )/value_type( 2 );
        }
        else if ( kind == SB9ShearKind::Bc1 )
        {
            f11bz = this->M_invJ0( 0, 0 )*( Jc( 0, 0 ) - Ja( 0, 0 ) )/value_type( 2 );
            f11g1 = this->M_invJ0( 0, 0 )*( invJc( 2, 2 )*Jc( 0, 0 ) + invJa( 2, 2 )*Ja( 0, 0 ) )/value_type( 2 );

            f12bz = this->M_invJ0( 0, 0 )*( Jc( 0, 1 ) - Ja( 0, 1 ) )/value_type( 2 );
            f12g1 = this->M_invJ0( 0, 0 )*( invJc( 2, 2 )*Jc( 0, 1 ) + invJa( 2, 2 )*Ja( 0, 1 ) )/value_type( 2 );

            f13bx = f11bz;
            f13by = f12bz;
            f13g3 = this->M_invJ0( 0, 0 )*( invJc( 0, 0 )*Jc( 0, 0 ) + invJc( 1, 0 )*Jc( 0, 1 ) +
                                            invJa( 0, 0 )*Ja( 0, 0 ) + invJa( 1, 0 )*Ja( 0, 1 ) )/value_type( 2 );

            f21bz = this->M_invJ0( 1, 0 )*( Jc( 0, 0 ) - Ja( 0, 0 ) )/value_type( 2 );
            f21g1 = this->M_invJ0( 1, 0 )*( invJc( 2, 2 )*Jc( 0, 0 ) + invJa( 2, 2 )*Ja( 0, 0 ) )/value_type( 2 );

            f22bz = this->M_invJ0( 1, 0 )*( Jc( 0, 1 ) - Ja( 0, 1 ) )/value_type( 2 );
            f22g1 = this->M_invJ0( 1, 0 )*( invJc( 2, 2 )*Jc( 0, 1 ) + invJa( 2, 2 )*Ja( 0, 1 ) )/value_type( 2 );

            f23bx = f21bz;
            f23by = f22bz;
            f23g3 = this->M_invJ0( 1, 0 )*( invJc( 0, 0 )*Jc( 0, 0 ) + invJc( 1, 0 )*Jc( 0, 1 ) +
                                            invJa( 0, 0 )*Ja( 0, 0 ) + invJa( 1, 0 )*Ja( 0, 1 ) )/value_type( 2 );
        }
        else
        {
            f21bz = this->M_invJ0( 1, 1 )*( Jb( 1, 0 ) - Jd( 1, 0 ) )/value_type( 2 );
            f21g2 = this->M_invJ0( 1, 1 )*( invJb( 2, 2 )*Jb( 1, 0 ) + invJd( 2, 2 )*Jd( 1, 0 ) )/value_type( 2 );

            f22bz = this->M_invJ0( 1, 1 )*( Jb( 1, 1 ) - Jd( 1, 1 ) )/value_type( 2 );
            f22g2 = this->M_invJ0( 1, 1 )*( invJb( 2, 2 )*Jb( 1, 1 ) + invJd( 2, 2 )*Jd( 1, 1 ) )/value_type( 2 );

            f23bx = f21bz;
            f23by = f22bz;
            f23g3 = this->M_invJ0( 1, 1 )*( invJb( 0, 1 )*Jb( 1, 0 ) + invJb( 1, 1 )*Jb( 1, 1 ) +
                                            invJd( 0, 1 )*Jd( 1, 0 ) + invJd( 1, 1 )*Jd( 1, 1 ) )/value_type( 2 );
        }

        value_type const aux11 = f11bz*Bz( node ) + f11g1*Vgamma( node, 0 );
        value_type const aux12 = f12bz*Bz( node ) + f12g1*Vgamma( node, 0 );
        value_type const aux13 = f13bx*Bx( node ) + f13g3*Vgamma( node, 2 ) + f13by*By( node );

        value_type const aux21 = f21bz*Bz( node ) + f21g1*Vgamma( node, 0 ) + f21g2*Vgamma( node, 1 );
        value_type const aux22 = f22bz*Bz( node ) + f22g1*Vgamma( node, 0 ) + f22g2*Vgamma( node, 1 );
        value_type const aux23 = f23bx*Bx( node ) + f23g3*Vgamma( node, 2 ) + f23by*By( node );

        return { row_type{ aux11, aux12, aux13 },
                 row_type{ aux21, aux22, aux23 } };
    }

    /// Cached auxiliary rows indexed by shear block and local geometric node.
    std::array<std::array<shear_aux_type, node_count>, 3> M_shear{};
};
} // namespace detail

/**
 * \brief Build the physical SB9 transverse-shear coefficient expression.
 *
 * The returned expression is the `Bc0` block used to assemble the physical
 * transverse-shear strain contribution. It is a six-component compact
 * symmetric vector expression in the same storage convention as the other SB9
 * strain blocks.
 *
 * \par Example
 * \code{.cpp}
 * auto shearWeight = cst( 5.0 / 4.0 ) * ( cst( 1.0 ) - zeta()*zeta() );
 * auto epsS = sb9Shear( u, shearWeight );
 *
 * a += integrate( _range=elements( mesh ),
 *                 _quad=sb9ThroughThicknessLobatto5(),
 *                 _expr=ddot( C, epsS, sb9Shear( v, shearWeight ) ) );
 * \endcode
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the `Bc0` coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bc0( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9ShearKernelCache,
                                                detail::SB9ShearKind::Bc0>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the first SB9 transverse-shear stabilization block.
 *
 * `Bc1` is not part of the physical strain vector. It is used in the
 * Hallquist transverse-shear stabilization term, typically together with
 * \ref sb9Bc2.
 *
 * \par Example
 * \code{.cpp}
 * a += integrate( _range=elements( mesh ),
 *                 _quad=sb9ThroughThicknessLobatto5(),
 *                 _expr=inner( sb9Bc1( u ), sb9Bc1( v ) ) );
 * \endcode
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the `Bc1` coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bc1( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9ShearKernelCache,
                                                detail::SB9ShearKind::Bc1>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the second SB9 transverse-shear stabilization block.
 *
 * `Bc2` complements \ref sb9Bc1 in the Hallquist transverse-shear
 * stabilization term. It is kept in this shear header, while the mode
 * stabilization blocks `Bs1..Bs4` live in `sb9_stabilization.hpp`.
 *
 * \par Example
 * \code{.cpp}
 * a += integrate( _range=elements( mesh ),
 *                 _quad=sb9ThroughThicknessLobatto5(),
 *                 _expr=inner( sb9Bc2( u ), sb9Bc2( v ) ) );
 * \endcode
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the `Bc2` coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bc2( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9ShearKernelCache,
                                                detail::SB9ShearKind::Bc2>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the weighted physical SB9 transverse-shear strain expression.
 *
 * The helper extracts the `xz` and `yz` slots from \ref sb9Bc0, applies the
 * supplied through-thickness shear weight, and returns a compact symmetric
 * Mandel vector with only the transverse-shear components populated.
 *
 * \par Example
 * \code{.cpp}
 * auto shearWeight = cst( 5.0 / 4.0 ) * ( cst( 1.0 ) - zeta()*zeta() );
 * auto epsTrial = sb9Shear( u, shearWeight );
 * auto epsTest = sb9Shear( v, shearWeight );
 *
 * a += integrate( _range=elements( mesh ),
 *                 _quad=sb9ThroughThicknessLobatto5(),
 *                 _expr=ddot( C, epsTrial, epsTest ) );
 * \endcode
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \tparam ShearWeightExprT Feel++ scalar expression type for the shear weight.
 * \param proxy Trial or test basis proxy.
 * \param shearWeight Scalar weight applied to the transverse-shear slots.
 * \return Compact symmetric Mandel vector expression for SB9 shear strain.
 */
template <detail::BasisProxyType ProxyType, typename ShearWeightExprT>
[[nodiscard]] inline auto
sb9Shear( ProxyType const& proxy, ShearWeightExprT const& shearWeight )
{
    auto bc0 = sb9Bc0( proxy );
    return mandel_vec<3>( cst( 0.0 ),
                          cst( 0.0 ),
                          cst( 0.0 ),
                          cst( 0.0 ),
                          shearWeight * component<4, 0>( bc0 ),
                          shearWeight * component<5, 0>( bc0 ) );
}
} // namespace vf
} // namespace Feel

#endif
