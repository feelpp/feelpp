/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
    SPDX-FileContributor: Hanna Chetouane

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
   \file sb9_stabilization.hpp
   \brief SB9 stabilization shell kinematic operators
 */
#ifndef FEELPP_VF_SB9_STABILIZATION_HPP
#define FEELPP_VF_SB9_STABILIZATION_HPP 1

#include <feel/feelvf/sb9_common.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{
/**
 * \brief Selects the SB9 stabilization coefficient blocks.
 *
 * The values are used as compile-time tags to select the corresponding
 * stabilization contribution in \ref SB9StabilizationKernelCache and
 * \ref SB9StabilizationOperator.
 */
enum class SB9StabilizationKind
{
    /// Mode stabilization contributions denoted Bs, obtained by splitting the
    /// coefficient matrix into four blocks.
    Bs1, Bs2, Bs3, Bs4,
    /// Transverse shear stabilization contributions denoted Bc, obtained by
    /// splitting the coefficient matrix into two blocks.
    Bc1, Bc2
};

/**
 * \brief Element-local cache for SB9 stabilization coefficients.
 *
 * The cache reuses \ref SB9KernelBase for the shell geometry data and computes
 * the stabilization coefficients. The coefficients are stored in a compact
 * representation and used during matrix assembly by
 * \ref fillStabilizationCoefficients.
 *
 * \tparam GeometryDataType Shell cell geometry cache type provided by
 *         ShellCellGeometryTensorBase.
 */
template <typename GeometryDataType>
class SB9StabilizationKernelCache : public SB9KernelBase<GeometryDataType>
{
public:
    /// Shared SB9 geometry cache base.
    using base_type = SB9KernelBase<GeometryDataType>;
    /// Scalar type used by the geometry and generated coefficients.
    using value_type = typename base_type::value_type;
    /// Number of geometric nodes in the current SB9 Q1 hexahedral element.
    static constexpr uint16_type node_count = base_type::node_count;
    /// Number of displacement components handled by SB9 bending operators.
    static constexpr uint16_type component_count = base_type::component_count;

    /**
     * \brief Build stabilization coefficients data for one element.
     *
     * This constructor assembles the intermediate stabilization coefficients
     * from the shell geometry data. The resulting coefficients are stored in a
     * compact per-node representation for later reuse during assembly.
     *
     * \param data Element-local shell geometry data.
     */
    explicit SB9StabilizationKernelCache( GeometryDataType const& data )
        :
        base_type( data )
    {
        auto Bc1f11bz = 0.5 * this->M_invJ0(0,0) * ( this->M_data.Jc(0,0) - this->M_data.Ja(0,0) );
        auto Bc1f11g1 = 0.5 * this->M_invJ0(0,0) * ( this->M_data.invJc(2,2) * this->M_data.Jc(0,0) + 
                                                     this->M_data.invJa(2,2) * this->M_data.Ja(0,0) );

        auto Bc1f12bz = 0.5 * this->M_invJ0(0,0) * ( this->M_data.Jc(0,1) - this->M_data.Ja(0,1) );
        auto Bc1f12g1 = 0.5 * this->M_invJ0(0,0) * ( this->M_data.invJc(2,2) * this->M_data.Jc(0,1) + 
                                                     this->M_data.invJa(2,2) * this->M_data.Ja(0,1) );

        auto Bc1f13g3 = 0.5 * this->M_invJ0(0,0) * ( this->M_data.invJc(0,0) * this->M_data.Jc(0,0) + 
                                                     this->M_data.invJc(1,0) * this->M_data.Jc(0,1) +
                                                     this->M_data.invJa(0,0) * this->M_data.Ja(0,0) +
                                                     this->M_data.invJa(1,0) * this->M_data.Ja(0,1) );

        auto Bc1f21bz = 0.5 * this->M_invJ0(1,0) * ( this->M_data.Jc(0,0) - this->M_data.Ja(0,0) );
        auto Bc1f21g1 = 0.5 * this->M_invJ0(1,0) * ( this->M_data.invJc(2,2) * this->M_data.Jc(0,0) + 
                                                     this->M_data.invJa(2,2) * this->M_data.Ja(0,0) );

        auto Bc1f22bz = 0.5 * this->M_invJ0(1,0) * ( this->M_data.Jc(0,1) - this->M_data.Ja(0,1) );
        auto Bc1f22g1 = 0.5 * this->M_invJ0(1,0) * ( this->M_data.invJc(2,2) * this->M_data.Jc(0,1) + 
                                                     this->M_data.invJa(2,2) * this->M_data.Ja(0,1) );

        auto Bc1f23g3 = 0.5 * this->M_invJ0(1,0) * ( this->M_data.invJc(0,0) * this->M_data.Jc(0,0) + 
                                                     this->M_data.invJc(1,0) * this->M_data.Jc(0,1) +
                                                     this->M_data.invJa(0,0) * this->M_data.Ja(0,0) +
                                                     this->M_data.invJa(1,0) * this->M_data.Ja(0,1) );


        auto Bc2f21bz = 0.5 * this->M_invJ0(1,1) * ( this->M_data.Jb(1,0) - this->M_data.Jd(1,0) );
        auto Bc2f21g2 = 0.5 * this->M_invJ0(1,1) * ( this->M_data.invJb(2,2) * this->M_data.Jb(1,0) + 
                                                     this->M_data.invJd(2,2) * this->M_data.Jd(1,0) );

        auto Bc2f22bz = 0.5 * this->M_invJ0(1,1) * ( this->M_data.Jb(1,1) - this->M_data.Jd(1,1) );
        auto Bc2f22g2 = 0.5 * this->M_invJ0(1,1) * ( this->M_data.invJb(2,2) * this->M_data.Jb(1,1) + 
                                                     this->M_data.invJd(2,2) * this->M_data.Jd(1,1) );

        auto Bc2f23g3 = 0.5 * this->M_invJ0(1,1) * ( this->M_data.invJb(0,1) * this->M_data.Jb(1,0) + 
                                                     this->M_data.invJb(1,1) * this->M_data.Jb(1,1) +
                                                     this->M_data.invJd(0,1) * this->M_data.Jd(1,0) +
                                                     this->M_data.invJd(1,1) * this->M_data.Jd(1,1) );

        for ( uint16_type node = 0; node < node_count; ++node )
        {
            M_Bc1[node][0][0] = Bc1f11bz * this->M_data.bz( node ) + Bc1f11g1 * this->M_data.vgamma( node, 0 );
            M_Bc1[node][1][0] = Bc1f12bz * this->M_data.bz( node ) + Bc1f12g1 * this->M_data.vgamma( node, 0 );
            M_Bc1[node][2][0] = Bc1f11bz * this->M_data.bx( node ) + Bc1f13g3 * this->M_data.vgamma( node, 2 ) + 
                                Bc1f12bz * this->M_data.by( node );

            M_Bc1[node][0][1] = Bc1f21bz * this->M_data.bz( node ) + Bc1f21g1 * this->M_data.vgamma( node, 0 );
            M_Bc1[node][1][1] = Bc1f22bz * this->M_data.bz( node ) + Bc1f22g1 * this->M_data.vgamma( node, 0 );
            M_Bc1[node][2][1] = Bc1f21bz * this->M_data.bx( node ) + Bc1f23g3 * this->M_data.vgamma( node, 2 ) + 
                                Bc1f22bz * this->M_data.by( node );

            M_Bc2[node][0] = Bc2f21bz * this->M_data.bz( node ) + Bc2f21g2 * this->M_data.vgamma( node, 1 );
            M_Bc2[node][1] = Bc2f22bz * this->M_data.bz( node ) + Bc2f22g2 * this->M_data.vgamma( node, 1 );
            M_Bc2[node][2] = Bc2f21bz * this->M_data.bx( node ) + Bc2f23g3 * this->M_data.vgamma( node, 2 ) + 
                             Bc2f22bz * this->M_data.by( node );

            for ( uint16_type i = 0; i < 4; ++i )
            {
                M_Bsgamma[node][i] = this->M_data.vgamma( node, i );
            }
        }
    }

    /**
     * \brief Fill the coefficient vector associated with one SB9 stabilization
     *        contribution.
     *
     * \tparam Kind Compile-time selector identifying the stabilization contribution.
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector.
     * \param node Local geometric node index.
     * \param component Displacement component index.
     */
    template <SB9StabilizationKind Kind, typename VectorType>
    void fillStabilizationCoefficients( VectorType& coeff, uint16_type node, uint16_type component ) const
    {
        static_assert( Kind == SB9StabilizationKind::Bs1 || Kind == SB9StabilizationKind::Bs2 || 
                       Kind == SB9StabilizationKind::Bs3 || Kind == SB9StabilizationKind::Bs4 ||
                       Kind == SB9StabilizationKind::Bc1 || Kind == SB9StabilizationKind::Bc2,
                       "unsupported SB9 stabilization vector kind" );  

        if constexpr ( Kind == SB9StabilizationKind::Bs1 )
            this->fillBs1Coefficients( coeff, component, M_Bsgamma[node][0] );
        else if constexpr ( Kind == SB9StabilizationKind::Bs2 )
            this->fillBs1Coefficients( coeff, component, M_Bsgamma[node][1] );
        else if constexpr ( Kind == SB9StabilizationKind::Bs3 )
            this->fillBs3Coefficients( coeff, component, M_Bsgamma[node][2] );
        else if constexpr ( Kind == SB9StabilizationKind::Bs4 )
            this->fillBs4Coefficients( coeff, component, M_Bsgamma[node][3] );
        else if constexpr ( Kind == SB9StabilizationKind::Bc1 )
            this->fillBc1Coefficients( coeff, component, M_Bc1[node][0][0], M_Bc1[node][1][0], M_Bc1[node][2][0],
                                                         M_Bc1[node][0][1], M_Bc1[node][1][1], M_Bc1[node][2][1] );
        else
            this->fillBc2Coefficients( coeff, component, M_Bc2[node][0], M_Bc2[node][1], M_Bc2[node][2] );
    }

private:
    /// Precomputed mode stabilization coefficients per node.
    std::array<std::array<value_type, 4>, node_count> M_Bsgamma{};
    /// Precomputed first transverse shear stabilization coefficients per node.
    std::array<std::array< std::array<value_type, 2>, 3>, node_count> M_Bc1{};
    /// Precomputed second transverse shear stabilization coefficients per node.
    std::array<std::array<value_type, 3>, node_count> M_Bc2{};
};
} // namespace detail

/**
 * \brief Build an SB9 mode stabilization coefficient expression.
 *
 * The returned expression evaluates one of the four mode stabilization
 * contribution (`Bs1`, `Bs2`, `Bs3` or `Bs4`) for a trial or test basis proxy
 * created with `trial(Vh, ...)` or `test(Vh, ...)`. Each expression represents
 * a three-component coefficient vector.
 * 
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the SB9 coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bs1( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9StabilizationOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9StabilizationKernelCache,
                                                detail::SB9StabilizationKind::Bs1,
                                                3>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bs2( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9StabilizationOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9StabilizationKernelCache,
                                                detail::SB9StabilizationKind::Bs2,
                                                3>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bs3( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9StabilizationOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9StabilizationKernelCache,
                                                detail::SB9StabilizationKind::Bs3,
                                                3>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bs4( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9StabilizationOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9StabilizationKernelCache,
                                                detail::SB9StabilizationKind::Bs4,
                                                3>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build a transverse shear stabilization coefficient expression.
 *
 * The returned expression evaluates one of the two transverse shear
 * stabilization coefficient vectors (`Bc1` or `Bc2`) for a trial or test basis
 * proxy created with `trial(Vh, ...)` or `test(Vh, ...)`. Each expression
 * represents a two-component coefficient vector.
 * 
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the SB9 coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bc1( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9StabilizationOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9StabilizationKernelCache,
                                                detail::SB9StabilizationKind::Bc1,
                                                2>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bc2( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9StabilizationOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9StabilizationKernelCache,
                                                detail::SB9StabilizationKind::Bc2,
                                                2>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the SB9 mode stabilization expression.
 *
 * This helper assembles the SB9 mode stabilization contribution from the four 
 * coefficient expressions (`Bs1`, `Bs2`, `Bs3` and `Bs4`). The returned
 * expression represents the element stiffness contribution associated with the
 * mode stabilization.
 * 
 * \tparam ProxyUType Feel++ trial basis proxy type.
 * \tparam ProxyVType Feel++ test basis proxy type.
 * \tparam XiExprT Feel++ expression type used for the reference `xi()`.
 * \tparam EtaExprT Feel++ expression type used for the reference `eta()`.
 * \tparam ZetaExprT Feel++ expression type used for the reference `zeta()`.
 * \tparam LambdaExprT Expression type of the first Lamé coefficient.
 * \tparam MuExprT Expression type of the second Lamé coefficient.
 * \param proxy_u Trial basis proxy.
 * \param proxy_v Test basis proxy.
 * \param xiExpr Expression representing the reference coordinate xi.
 * \param etaExpr Expression representing the reference coordinate eta.
 * \param zetaExpr Expression representing the reference coordinate zeta.
 * \param lambda First Lamé coefficient.
 * \param mu Second Lamé coefficient.
 * \param coefStabMembrane Membrane stabilization scaling factor, set to 1.0 by default.
 * \param coefStabBending Bending stabilization scaling factor, set to 1.0 by default.
 * \param coefStabPinching Pinching stabilization scaling factor, set to 1.0 by default.
 * \return Feel++ expression representing the element stiffness contribution
 *         associated with the mode stabilization.
 */
template <detail::BasisProxyType ProxyUType, detail::BasisProxyType ProxyVType, typename XiExprT, typename EtaExprT, typename ZetaExprT, typename LambdaExprT, typename MuExprT>
[[nodiscard]] inline auto
sb9ModeStabilization( ProxyUType const& proxy_u, ProxyVType const& proxy_v, XiExprT const& xiExpr, EtaExprT const& etaExpr, ZetaExprT const& zetaExpr, LambdaExprT const& lambda, MuExprT const& mu,
                       double const& coefStabMembrane = 1.0, double const& coefStabBending = 1.0, double const& coefStabPinching = 1.0 )
{
    auto Cdiag =  lambda + 2.0 * mu;
    auto Volume = shellThickness() * shellArea0(); 

    auto Ds1 = shellInvJ0_00() * shellInvJ0_00() * Cdiag * Volume / 3.0;
    auto Ds2 = ( shellInvJ0_10() * shellInvJ0_10() + shellInvJ0_11() * shellInvJ0_11() ) * Cdiag * Volume / 3.0;
    auto Ds3 = coefStabPinching * shellInvJ0_22() * shellInvJ0_22() * Cdiag * Volume / 3.0;

    auto Bs1_u = sb9Bs1( proxy_u );
    auto Bs1_v = sb9Bs1( proxy_v );
    auto Bs2_u = sb9Bs2( proxy_u );
    auto Bs2_v = sb9Bs2( proxy_v );
    auto Bs3_u = sb9Bs3( proxy_u );
    auto Bs3_v = sb9Bs3( proxy_v );
    auto Bs4_u = sb9Bs4( proxy_u );
    auto Bs4_v = sb9Bs4( proxy_v );

    auto Ks1 = Ds1 * ( coefStabMembrane * trans( component<0,0>( Bs3_u ) ) * component<0,0>( Bs3_v ) + 
                       coefStabBending/3 * trans( component<0,0>( Bs4_u ) ) * component<0,0>( Bs4_v ) );

    auto Ks2 = Ds2 * ( coefStabMembrane * trans( component<1,0>( Bs3_u ) ) * component<1,0>( Bs3_v ) + 
                       coefStabBending / 3.0 * trans( component<1,0>( Bs4_u ) ) * component<1,0>( Bs4_v ) );

    auto Ks3 = Ds3 * ( trans( component<2,0>( Bs1_u ) ) * component<2,0>( Bs1_v ) + 
                       trans( component<2,0>( Bs2_u ) ) * component<2,0>( Bs2_v ) + 
                       1.0 / 3.0 * trans( component<2,0>( Bs4_u ) ) * component<2,0>( Bs4_v ) );

    return Ks1 + Ks2 + Ks3;
}

/**
 * \brief Build the SB9 mode stabilization expression using the default
 * reference coordinate expressions \ref xi(), \ref eta() and \ref zeta().
 */
template <detail::BasisProxyType ProxyUType, detail::BasisProxyType ProxyVType, typename LambdaExprT, typename MuExprT>
[[nodiscard]] inline auto
sb9ModeStabilization( ProxyUType const& proxy_u, ProxyVType const& proxy_v, LambdaExprT const& lambda, MuExprT const& mu,
                 double const& coefStabMembrane = 1.0, double const& coefStabBending = 1.0, double const& coefStabPinching = 1.0 )
{
    return sb9ModeStabilization( proxy_u, proxy_v, xi(), eta(), zeta(), lambda, mu, coefStabMembrane, coefStabBending, coefStabPinching );
}

/**
 * \brief Build the SB9 transverse shear stabilization expression.
 *
 * This helper assembles the SB9 transverse shear stabilization contribution
 * from the two coefficient expressions (`Bc1` and `Bc2`). The returned
 * expression represents the element stiffness contribution associated with the
 * transverse shear stabilization.
 * 
 * \tparam ProxyUType Feel++ trial basis proxy type.
 * \tparam ProxyVType Feel++ test basis proxy type.
 * \tparam MuExprT Expression type of the second Lamé coefficient.
 * \param proxy_u Trial basis proxy.
 * \param proxy_v Test basis proxy.
 * \param mu Second Lamé coefficient.
 * \return Feel++ expression representing the element stiffness contribution
 *         associated with the transverse shear stabilization.
 */
template <detail::BasisProxyType ProxyUType, detail::BasisProxyType ProxyVType, typename MuExprT>
[[nodiscard]] inline auto
sb9ShearingStabilization( ProxyUType const& proxy_u, ProxyVType const& proxy_v, MuExprT const& mu )
{
    auto Bc1_u = sb9Bc1( proxy_u );
    auto Bc1_v = sb9Bc1( proxy_v );
    auto Bc2_u = sb9Bc2( proxy_u );
    auto Bc2_v = sb9Bc2( proxy_v );

    auto Cstab = 5.0/6.0* mu* shellThickness() * shellArea0() / 3.0;

    auto Kc1 = Cstab * ( trans( component<0,0>( Bc1_u ) ) * component<0,0>( Bc1_v ) );
    auto Kc2 = Cstab * ( trans( component<1,0>( Bc1_u ) ) * component<1,0>( Bc1_v ) + 
                         trans( component<1,0>( Bc2_u ) ) * component<1,0>( Bc2_v ) );

    return Kc1 + Kc2;
}

/**
 * \brief Build the complete SB9 stabilization expression.
 *
 * This helper combines the SB9 mode stabilization contribution and the
 * transverse shear stabilization contribution to build the complete
 * stabilization term of the SB9 formulation.
 */
template <detail::BasisProxyType ProxyUType, detail::BasisProxyType ProxyVType, typename LambdaExprT, typename MuExprT>
[[nodiscard]] inline auto
sb9Stabilization( ProxyUType const& proxy_u, ProxyVType const& proxy_v, LambdaExprT const& lambda, MuExprT const& mu,
                 double const& coefStabMembrane = 1.0, double const& coefStabBending = 1.0, double const& coefStabPinching = 1.0 )
{
    return sb9ModeStabilization( proxy_u, proxy_v, xi(), eta(), zeta(), lambda, mu, coefStabMembrane, coefStabBending, coefStabPinching ) + 
           sb9ShearingStabilization( proxy_u, proxy_v, mu);
}
} // namespace vf
} // namespace Feel

#endif
