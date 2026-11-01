/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file integrator_localform_dispatch.hpp
   \brief Small lowering dispatch helpers for Integrator::assemble.
 */
#ifndef FEELPP_VF_DETAIL_INTEGRATOR_LOCALFORM_DISPATCH_HPP
#define FEELPP_VF_DETAIL_INTEGRATOR_LOCALFORM_DISPATCH_HPP 1

template<RangeConcept Elements, typename Im, VfExprConcept Expr, typename Im2>
    requires QuadOrderConcept<Im> && QuadOrderConcept<Im2>
template<typename Elem1, typename Elem2, typename FormType>
bool
Integrator<Elements, Im, Expr, Im2>::tryAssembleLoweredScalarLocalform( std::shared_ptr<Elem1> const& __u,
                                                                        std::shared_ptr<Elem2> const& __v,
                                                                        FormType& __form ) const
{
    if ( M_local_form_lowering == local_form_lowering::disabled )
        return false;
    if constexpr ( !detail::can_lower_scalar_bilinear_localform_for_form_v<
                       expression_type, FormType> )
        if ( M_local_form_lowering != local_form_lowering::enabled )
            return false;

    typedef typename boost::is_same<typename eval::gmc_type::element_type,typename Elem1::mesh_type::element_type>::type same1_mesh_type;
    typedef typename boost::is_same<typename eval::gmc_type::element_type,typename Elem2::mesh_type::element_type>::type same2_mesh_type;
    typedef typename boost::mpl::and_< same1_mesh_type,same2_mesh_type>::type same_mesh_type;

    if constexpr (
        iDim == MESH_ELEMENTS && same_mesh_type::value &&
        ( detail::can_lower_scalar_bilinear_localform_for_form_v<
              expression_type, FormType> ||
          detail::can_lower_scalar_h1_stiffness_localform_for_form_v<
              expression_type, FormType> ||
          detail::can_lower_anisotropic_scalar_h1_stiffness_localform_for_form_v<
              expression_type, FormType> ||
          detail::can_lower_vector_mass_localform_for_form_v<
              expression_type, FormType> ||
          detail::can_lower_curl_curl_localform_for_form_v<
              expression_type, FormType> ||
          detail::can_lower_div_div_localform_for_form_v<
              expression_type, FormType> ||
          detail::can_lower_value_div_localform_for_form_v<
              expression_type, FormType> ||
          detail::can_lower_div_value_localform_for_form_v<
              expression_type, FormType> ) )
    {
        auto lowered_expr = [&]()
        {
            if constexpr (
                detail::can_lower_scalar_bilinear_localform_for_form_v<
                    expression_type, FormType> )
                return detail::lower_scalar_bilinear_localform(
                    this->expression() );
            else if constexpr (
                detail::can_lower_scalar_h1_stiffness_localform_for_form_v<
                    expression_type, FormType> )
                return detail::lower_scalar_h1_stiffness_localform(
                    this->expression() );
            else if constexpr (
                detail::can_lower_anisotropic_scalar_h1_stiffness_localform_for_form_v<
                    expression_type, FormType> )
                return detail::lower_anisotropic_scalar_h1_stiffness_localform(
                    this->expression() );
            else if constexpr (
                detail::can_lower_vector_mass_localform_for_form_v<
                    expression_type, FormType> )
                return detail::lower_canonical_bilinear_localform<
                    detail::vector_mass_contraction>( this->expression() );
            else if constexpr (
                detail::can_lower_curl_curl_localform_for_form_v<
                    expression_type, FormType> )
                return detail::lower_canonical_bilinear_localform<
                    detail::curl_curl_contraction>( this->expression() );
            else if constexpr (
                detail::can_lower_div_div_localform_for_form_v<
                    expression_type, FormType> )
                return detail::lower_canonical_bilinear_localform<
                    detail::divergence_divergence_contraction>(
                        this->expression() );
            else if constexpr (
                detail::can_lower_value_div_localform_for_form_v<
                    expression_type, FormType> )
                return detail::lower_canonical_bilinear_localform<
                    detail::value_divergence_contraction>(
                        this->expression() );
            else
                return detail::lower_canonical_bilinear_localform<
                    detail::divergence_value_contraction>(
                        this->expression() );
        }();
        using lowered_expr_type = decltype( lowered_expr );
        using lowered_quadptloc_ptrtype = std::shared_ptr<QuadPtLocalization<Elements, Im, lowered_expr_type> >;
        auto lowered_integrator =
            Integrator<Elements, Im, lowered_expr_type, Im2>(
                M_elts, M_im, lowered_expr, M_gt,
                M_geometry_evaluation, M_im2,
                M_use_harts,
                lowered_quadptloc_ptrtype{}, M_local_form_lowering );
        lowered_integrator.setBeginElement( M_eltbegin );
        lowered_integrator.setEndElement( M_eltend );
        lowered_integrator.assemble( __u, __v, __form );
        return true;
    }

    return false;
}

template<RangeConcept Elements, typename Im, VfExprConcept Expr, typename Im2>
    requires QuadOrderConcept<Im> && QuadOrderConcept<Im2>
template<typename Elem1, typename FormType>
bool
Integrator<Elements, Im, Expr, Im2>::tryAssembleLoweredScalarLocalform( std::shared_ptr<Elem1> const& __v,
                                                                        FormType& __form ) const
{
    if ( M_local_form_lowering == local_form_lowering::disabled )
        return false;

    typedef typename boost::is_same<typename eval::gmc_type::element_type,typename Elem1::mesh_type::element_type>::type same_mesh_type;

    if constexpr ( iDim == MESH_ELEMENTS && same_mesh_type::value &&
                   detail::can_lower_scalar_linear_localform_for_form_v<expression_type, FormType> )
    {
        auto lowered_expr = detail::lower_scalar_linear_localform( this->expression() );
        using lowered_expr_type = decltype( lowered_expr );
        using lowered_quadptloc_ptrtype = std::shared_ptr<QuadPtLocalization<Elements, Im, lowered_expr_type> >;
        auto lowered_integrator =
            Integrator<Elements, Im, lowered_expr_type, Im2>(
                M_elts, M_im, lowered_expr, M_gt,
                M_geometry_evaluation, M_im2,
                M_use_harts,
                lowered_quadptloc_ptrtype{}, M_local_form_lowering );
        lowered_integrator.setBeginElement( M_eltbegin );
        lowered_integrator.setEndElement( M_eltend );
        lowered_integrator.assemble( __v, __form );
        return true;
    }

    return false;
}
#endif
