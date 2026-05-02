/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
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
    typedef typename boost::is_same<typename eval::gmc_type::element_type,typename Elem1::mesh_type::element_type>::type same1_mesh_type;
    typedef typename boost::is_same<typename eval::gmc_type::element_type,typename Elem2::mesh_type::element_type>::type same2_mesh_type;
    typedef typename boost::mpl::and_< same1_mesh_type,same2_mesh_type>::type same_mesh_type;

    if constexpr ( iDim == MESH_ELEMENTS && same_mesh_type::value &&
                   detail::can_lower_scalar_bilinear_localform_for_form_v<expression_type, FormType> )
    {
        auto lowered_expr = detail::lower_scalar_bilinear_localform( this->expression() );
        using lowered_expr_type = decltype( lowered_expr );
        using lowered_quadptloc_ptrtype = std::shared_ptr<QuadPtLocalization<Elements, Im, lowered_expr_type> >;
        auto lowered_integrator =
            Integrator<Elements, Im, lowered_expr_type, Im2>( M_elts, M_im, lowered_expr, M_gt, M_im2,
                                                              M_use_tbb, M_use_harts, M_grainsize, M_partitioner,
                                                              lowered_quadptloc_ptrtype{} );
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
    typedef typename boost::is_same<typename eval::gmc_type::element_type,typename Elem1::mesh_type::element_type>::type same_mesh_type;

    if constexpr ( iDim == MESH_ELEMENTS && same_mesh_type::value &&
                   detail::can_lower_scalar_linear_localform_for_form_v<expression_type, FormType> )
    {
        auto lowered_expr = detail::lower_scalar_linear_localform( this->expression() );
        using lowered_expr_type = decltype( lowered_expr );
        using lowered_quadptloc_ptrtype = std::shared_ptr<QuadPtLocalization<Elements, Im, lowered_expr_type> >;
        auto lowered_integrator =
            Integrator<Elements, Im, lowered_expr_type, Im2>( M_elts, M_im, lowered_expr, M_gt, M_im2,
                                                              M_use_tbb, M_use_harts, M_grainsize, M_partitioner,
                                                              lowered_quadptloc_ptrtype{} );
        lowered_integrator.setBeginElement( M_eltbegin );
        lowered_integrator.setEndElement( M_eltend );
        lowered_integrator.assemble( __v, __form );
        return true;
    }

    return false;
}
#endif
