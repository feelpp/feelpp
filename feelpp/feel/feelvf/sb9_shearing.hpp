/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
   \file sb9_shearing.hpp
   \brief SB9 transverse-shearing shell kinematic operators
 */
#ifndef FEELPP_VF_SB9_SHEARING_HPP
#define FEELPP_VF_SB9_SHEARING_HPP 1

#include <feel/feelvf/sb9_common.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{
/**
 * \brief Selects the SB9 transverse-shearing coefficient family.
 *
 * The values are used as compile-time tags by \ref SB9ShearingKernelCache and
 * \ref SB9VectorOperator.
 * 
 * This enumeration is kept to preserve the same architecture as the SB9 membrane, 
 * bending and pinching operators.
 */
enum class SB9ShearingKind
{
    /// Mid-surface transverse-shearing coefficient matrix, usually denoted Bts0.
    Bts0
};

/**
 * \brief Element-local cache for SB9 transverse-shearing coefficients.
 *
 * The cache reuses \ref SB9KernelBase for frame and Jacobian data and computes the
 * transverse-shearing derivative coefficients. The resulting coefficients are later
 * projected into Mandel symmetric storage by \ref fillVectorCoefficients.
 *
 * \tparam GeometryDataType Shell cell geometry cache type provided by
 *         ShellCellGeometryTensorBase.
 */
template <typename GeometryDataType>
class SB9ShearingKernelCache : public SB9KernelBase<GeometryDataType>
{
public:
    /// Shared SB9 geometry cache base.
    using base_type = SB9KernelBase<GeometryDataType>;
    /// Scalar type used by the geometry and generated coefficients.
    using value_type = typename base_type::value_type;
    /// Number of geometric nodes in the current SB9 Q1 hexahedral element.
    static constexpr uint16_type node_count = base_type::node_count;
    /// Number of displacement components handled by SB9 transverse-shearing operator.
    static constexpr uint16_type component_count = base_type::component_count;

    /**
     * \brief Build transverse-shearing coefficient data for one element.
     *
     * This constructor assembles the intermediate SB9 shearing matrices and
     * computes the element-local mid-surface transverse-shearing operator
     * from the shell geometry data. The resulting coefficients are stored in 
     * a compact per-node representation for later reuse during assembly.
     *
     * \param data Element-local shell geometry data.
     */
    explicit SB9ShearingKernelCache( GeometryDataType const& data )
        :
        base_type( data )
    {
        Eigen::Matrix<value_type, 2, 24> Bts = Eigen::Matrix<value_type, 2, 24>::Zero();
        for ( uint16_type node = 0; node < node_count; ++node )
        {
            Bts( 0, node ) = this->M_data.bz( node );
            Bts( 0, 2*node_count + node ) = this->M_data.bx( node );

            Bts( 1, node_count + node ) = this->M_data.bz( node );
            Bts( 1, 2*node_count + node ) = this->M_data.by( node );
        }

        Eigen::Matrix<value_type, 2, 24> Bats = Eigen::Matrix<value_type, 2, 24>::Zero();
        for ( uint16_type node = 0; node < node_count; ++node )
        {
            Bats( 0, node ) = - this->M_data.invJa( 2, 2 ) * this->M_data.vgamma( node, 0 );
            Bats( 0, 2*node_count + node ) = - this->M_data.invJa( 0, 0 ) * this->M_data.vgamma( node, 2 );

            Bats( 1, node_count + node ) = - this->M_data.invJa( 2, 2 ) * this->M_data.vgamma( node, 0 );
            Bats( 1, 2*node_count + node ) = - this->M_data.invJa( 1, 0 ) * this->M_data.vgamma( node, 2 );
        }
        Eigen::Matrix<value_type, 2, 24> Bbts = Eigen::Matrix<value_type, 2, 24>::Zero();
        for ( uint16_type node = 0; node < node_count; ++node )
        {
            Bbts( 0, node ) = this->M_data.invJb( 2, 2 ) * this->M_data.vgamma( node, 1 );
            Bbts( 0, 2*node_count + node ) = this->M_data.invJb( 0, 1 ) * this->M_data.vgamma( node, 2 );

            Bbts( 1, node_count + node ) = this->M_data.invJb( 2, 2 ) * this->M_data.vgamma( node, 1 );
            Bbts( 1, 2*node_count + node ) = this->M_data.invJb( 1, 1 ) * this->M_data.vgamma( node, 2 );
        }
        Eigen::Matrix<value_type, 2, 24> Bcts = Eigen::Matrix<value_type, 2, 24>::Zero();
        for ( uint16_type node = 0; node < node_count; ++node )
        {
            Bcts( 0, node ) = this->M_data.invJc( 2, 2 ) * this->M_data.vgamma( node, 0 );
            Bcts( 0, 2*node_count + node ) = this->M_data.invJc( 0, 0 ) * this->M_data.vgamma( node, 2 );

            Bcts( 1, node_count + node ) = this->M_data.invJc( 2, 2 ) * this->M_data.vgamma( node, 0 );
            Bcts( 1, 2*node_count + node ) = this->M_data.invJc( 1, 0 ) * this->M_data.vgamma( node, 2 );
        }
        Eigen::Matrix<value_type, 2, 24> Bdts = Eigen::Matrix<value_type, 2, 24>::Zero();
        for ( uint16_type node = 0; node < node_count; ++node )
        {
            Bdts( 0, node ) = - this->M_data.invJd( 2, 2 ) * this->M_data.vgamma( node, 1 );
            Bdts( 0, 2*node_count + node ) = - this->M_data.invJd( 0, 1 ) * this->M_data.vgamma( node, 2 );

            Bdts( 1, node_count + node ) = - this->M_data.invJd( 2, 2 ) * this->M_data.vgamma( node, 1 );
            Bdts( 1, 2*node_count + node ) = - this->M_data.invJd( 1, 1 ) * this->M_data.vgamma( node, 2 );
        }
        Eigen::Matrix<value_type, 2, 24> Ba = Bts + Bats;
        Eigen::Matrix<value_type, 2, 24> Bb = Bts + Bbts;
        Eigen::Matrix<value_type, 2, 24> Bc = Bts + Bcts;
        Eigen::Matrix<value_type, 2, 24> Bd = Bts + Bdts;

        Eigen::Matrix<value_type, 8, 24> Ce = Eigen::Matrix<value_type, 8, 24>::Zero();
        Eigen::Matrix<value_type, 2, 2> Ja;
        Eigen::Matrix<value_type, 2, 2> Jb;
        Eigen::Matrix<value_type, 2, 2> Jc;
        Eigen::Matrix<value_type, 2, 2> Jd;
        for ( uint16_type i = 0; i < 2; i++ )
        {
            for ( uint16_type j = 0; j < 2; j++ )
            {
                Ja(i,j) = this->M_data.Ja(i,j);
                Jb(i,j) = this->M_data.Jb(i,j);
                Jc(i,j) = this->M_data.Jc(i,j);
                Jd(i,j) = this->M_data.Jd(i,j);
            }
        }
        Ce.block(0, 0, 2, 24) = Ja * Ba;
        Ce.block(2, 0, 2, 24) = Jb * Bb;
        Ce.block(4, 0, 2, 24) = Jc * Bc;
        Ce.block(6, 0, 2, 24) = Jd * Bd;

        Eigen::Matrix<value_type, 2, 8> W = Eigen::Matrix<value_type, 2, 8>::Zero();
        value_type w = value_type( 0.5 );
        W( 0, 0 ) = w * this->M_invJ0( 0, 0 );
        W( 0, 4 ) = w * this->M_invJ0( 0, 0 );
        W( 1, 0 ) = w * this->M_invJ0( 1, 0 );
        W( 1, 3 ) = w * this->M_invJ0( 1, 1 );
        W( 1, 4 ) = w * this->M_invJ0( 1, 0 );
        W( 1, 7 ) = w * this->M_invJ0( 1, 1 );
        Eigen::Matrix<value_type, 2, 24> Bts0 = W * Ce;
        
        for ( uint16_type node = 0; node < node_count; ++node )
        {
            M_shearing[node][0][0] = Bts0( 0, node );
            M_shearing[node][1][0] = Bts0( 0, node_count + node );
            M_shearing[node][2][0] = Bts0( 0, node_count*2 + node );
            M_shearing[node][0][1] = Bts0( 1, node );
            M_shearing[node][1][1] = Bts0( 1, node_count + node );
            M_shearing[node][2][1] = Bts0( 1, node_count*2 + node );
        }
    }

    /**
     * \brief Fill one SB9 transverse-shearing-family coefficient vector.
     *
     * \tparam Kind Compile-time selector, \ref SB9BendingKind::Bts0.
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector in Feel++ symmetric storage order.
     * \param node Local geometric node index.
     * \param component Displacement component index.
     */
    template <SB9ShearingKind Kind, typename VectorType>
    void fillVectorCoefficients( VectorType& coeff, uint16_type node, uint16_type component ) const
    {
        static_assert( Kind == SB9ShearingKind::Bts0, "unsupported SB9 transverse-shearing vector kind" );

        this->fillShearingCoefficients( coeff, component, M_shearing[node][0][0], M_shearing[node][1][0], M_shearing[node][2][0], 
                                                          M_shearing[node][0][1], M_shearing[node][1][1], M_shearing[node][2][1] );
    }

private:
    /// Precomputed transverse-shearing derivative coefficients per node and shear component.
    std::array<std::array<std::array<value_type, 2>, 3>, node_count> M_shearing{};
};
} // namespace detail

/**
 * \brief Build the SB9 mid-surface transverse-shearing coefficient expression.
 *
 * The returned expression evaluates the `Bts0` contribution for a trial or test
 * basis proxy created with `trial(Vh, ...)` or `test(Vh, ...)`.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the SB9 `Bts0` coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bts0( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9ShearingKernelCache,
                                                detail::SB9ShearingKind::Bts0>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the SB9 transverse-shearing Mandel strain expression.
 *
 * This helper combines the transverse-shearing coefficient expressions and the 
 * Reissner function as `5/4 ( 1 - zeta² ) Bts0` in the shearing symmetric-storage 
 * entries and returns a six-component Mandel vector with all other components set
 * to zero.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \tparam ZetaExprT Feel++ expression type used for the through-thickness
 *         coordinate, typically `zeta()`.
 * \param proxy Trial or test basis proxy.
 * \param zetaExpr Through-thickness coordinate or scaling expression.
 * \return Feel++ Mandel-vector expression for the SB9 shearing strain.
 */
template <detail::BasisProxyType ProxyType, typename ZetaExprT>
[[nodiscard]] inline auto
sb9Shearing( ProxyType const& proxy, ZetaExprT const& zetaExpr )
{
    auto bts0 = sb9Bts0( proxy );
    auto reissner_func = cst( 1.25 ) * ( cst( 1 ) - zetaExpr * zetaExpr );

    return mandel_vec<3>( cst( 0.0 ),
                          cst( 0.0 ),
                          cst( 0.0 ),
                          cst( 0.0 ),
                          component<4, 0>( bts0 ) * reissner_func,
                          component<5, 0>( bts0 ) * reissner_func );
}
} // namespace vf
} // namespace Feel

#endif
