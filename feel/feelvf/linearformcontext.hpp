/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-04-27

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file linearformcontext.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-04-27
 */
#ifndef __LinearFormContext_H
#define __LinearFormContext_H 1

namespace Feel
{
namespace vf
{
namespace detail
{
//
// Context class for linear forms
//
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::Context( form_type& __form,
                                                                                                                                                   map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                                   map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                                   map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                                   ExprT const& expr,
                                                                                                                                                   IM const& im )
    : //super(),
      M_form( __form ),
      M_test_dof( __form.functionSpace()->dof().get() ),
      M_lb( __form.blockList() ),

      M_test_pc( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), fusion::at_key<gmc<0>>( _gmcTest )->pc()->nodes() ) ),
      M_test_pc_face( precomputeTestBasisAtPoints( im ) ),

      M_gmc( _gmcTest ),
      M_gmc_left( fusion::at_key<gmc<0>>( _gmcTest ) ),
      M_left_map( fusion::make_map<gmc<0>>( M_gmc_left ) ),
      M_test_fec( fusion::transform( M_gmc, vf::detail::FEContextInit<0, form_context_type>( __form.testFiniteElement<UseMortar>(), *this ) ) ),
      M_test_fec0( fusion::make_map<gmc<0>>( fusion::at_key<gmc<0>>( M_test_fec ) ) ),
      M_rep(),
      M_rep_2(),
      M_rep_mortar(),
      M_eval0_expr( new eval0_expr_type( expr, _gmcExpr, M_test_fec0 ) ),
      M_eval1_expr(),
      M_integrator( im )
{
    M_eval0_expr->init( im );
}

template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
template <typename IM2>
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::Context( form_type& __form,
                                                                                                                                                   map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                                   map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                                   map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                                   ExprT const& expr,
                                                                                                                                                   IM const& im,
                                                                                                                                                   IM2 const& im2 )
    : //super(),
      M_form( __form ),
      M_test_dof( __form.functionSpace()->dof().get() ),
      M_lb( __form.blockList() ),

      M_test_pc( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), im2.points() ) ),
      M_test_pc_face( precomputeTestBasisAtPoints( im2 ) ),

      M_gmc( _gmcTest ),
      M_gmc_left( fusion::at_key<gmc<0>>( _gmcTest ) ),
      M_left_map( fusion::make_map<gmc<0>>( M_gmc_left ) ),
      M_test_fec( fusion::transform( M_gmc, vf::detail::FEContextInit<0, form_context_type>( __form.testFiniteElement<UseMortar>(), *this ) ) ),
      M_test_fec0( fusion::make_map<gmc<0>>( fusion::at_key<gmc<0>>( M_test_fec ) ) ),
      M_rep(),
      M_rep_2(),
      M_rep_mortar(),
      M_eval0_expr( new eval0_expr_type( expr, _gmcExpr, M_test_fec0 ) ),
      M_eval1_expr(),
      M_integrator( im )
{
    M_eval0_expr->init( im2 );
}
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
template <typename IMExpr, typename IMTest>
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::Context( form_type& __form,
                                                                                                                                                   map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                                   map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                                   ExprT const& expr,
                                                                                                                                                   IM const& im,
                                                                                                                                                   IMExpr const& imExpr, IMTest const& imTest )
    : //super(),
      M_form( __form ),
      M_test_dof( __form.functionSpace()->dof().get() ),
      M_lb( __form.blockList() ),

      M_test_pc( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), fusion::at_key<gmc<0>>( _gmcTest )->pc()->nodes() ) ),
      M_test_pc_face( precomputeTestBasisAtPoints( imTest ) ),

      M_gmc( _gmcTest ),
      M_gmc_left( fusion::at_key<gmc<0>>( _gmcTest ) ),
      M_left_map( fusion::make_map<gmc<0>>( M_gmc_left ) ),
      M_test_fec( fusion::transform( M_gmc, vf::detail::FEContextInit<0, form_context_type>( __form.testFiniteElement<UseMortar>(), *this ) ) ),
      M_test_fec0( fusion::make_map<gmc<0>>( fusion::at_key<gmc<0>>( M_test_fec ) ) ),
      M_rep(),
      M_rep_2(),
      M_rep_mortar(),
      M_eval0_expr( new eval0_expr_type( expr, _gmcExpr, M_test_fec0 ) ),
      M_eval1_expr(),
      M_integrator( im )
{
    M_eval0_expr->init( imExpr );
}
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
template <typename IM2>
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::Context( form_type& __form,
                                                                                                                                                   map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                                   map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                                   map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                                   ExprT const& expr,
                                                                                                                                                   IM const& im,
                                                                                                                                                   IM2 const& im2,
                                                                                                                                                   mpl::int_<2> )
    : //super(),
      M_form( __form ),
      M_test_dof( __form.functionSpace()->dof().get() ),
      M_lb( __form.blockList() ),

      M_test_pc( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), im2.points() ) ),
      M_test_pc_face( precomputeTestBasisAtPoints( im2 ) ),

      M_gmc( _gmcTest ),
      M_gmc_left( fusion::at_key<gmc<0>>( _gmcTest ) ),
      M_gmc_right( fusion::at_key<gmc1>( _gmcTest ) ),
      M_left_map( fusion::make_map<gmc<0>>( M_gmc_left ) ),
      M_right_map( fusion::make_map<gmc<0>>( M_gmc_right ) ),
      M_test_fec( fusion::transform( M_gmc, vf::detail::FEContextInit<0, form_context_type>( __form.testFiniteElement<UseMortar>(), *this ) ) ),
      M_test_fec0( fusion::make_map<gmc<0>>( fusion::at_key<gmc<0>>( M_test_fec ) ) ),
      M_test_fec1( fusion::make_pair<gmc1>( fusion::at_key<gmc1>( M_test_fec ) ) ),
      M_rep(),
      M_rep_2(),
      M_rep_mortar(),
      M_eval0_expr( new eval0_expr_type( expr, _gmcExpr, M_test_fec0 ) ),
      M_eval1_expr( new eval1_expr_type( expr, _gmcExpr, M_test_fec1 ) ),
      M_integrator( im )

{
    M_eval0_expr->init( im2 );
    M_eval1_expr->init( im2 );
}
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
void LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                                       map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                                       map_geometric_mapping_expr_context_type const& _gmcExpr )
{
    M_gmc = _gmcTest;
    M_gmc_left = fusion::at_key<gmc<0>>( _gmcTest );
    M_left_map = fusion::make_map<gmc<0>>( M_gmc_left );
    fusion::for_each( M_test_fec, vf::detail::FEContextUpdate<0, form_context_type>( _gmcTest, *this ) );
    M_test_fec0 = fusion::make_map<gmc<0>>( fusion::at_key<gmc<0>>( M_test_fec ) );
    M_eval0_expr->update( _gmcExpr, M_test_fec0 );

    M_integrator.update( *fusion::at_key<gmc<0>>( _gmcExpr ) );
}
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
void LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                                                          map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                                                          map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                                                          std::vector<boost::tuple<size_type, size_type>> const& indexLocalToQuad )
{
    M_gmc = _gmcTest;
    M_gmc_left = fusion::at_key<gmc<0>>( _gmcTest );
    M_left_map = fusion::make_map<gmc<0>>( M_gmc_left );
    precomputeBasisAtPoints( fusion::at_key<gmc<0>>( _gmcTest )->xRefs() );
    fusion::for_each( M_test_fec, vf::detail::FEContextUpdate<0, form_context_type>( _gmcTest, *this ) );
    M_test_fec0 = fusion::make_map<gmc<0>>( fusion::at_key<gmc<0>>( M_test_fec ) );
    M_eval0_expr->update( _gmcExpr, M_test_fec0 );

    M_integrator.update( *fusion::at_key<gmc<0>>( _gmcExpr ), indexLocalToQuad );
}

template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
void LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                                       map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                                       map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                                       mpl::int_<2> )
{
//    typedef mpl::int_<fusion::result_of::template size<map_geometric_mapping_context_type>::type::value> map_size;
/*BOOST_MPL_ASSERT_MSG( (mpl::equal_to<mpl::int_<map_size::value>,mpl::int_<2> >::value),
                          INVALID_GEOMAP,
                          (map_size,map_geometric_mapping_context_type ));*/
//M_gmc = _gmc;
#if 0
    M_gmc_left = fusion::at_key<gmc<0> >( _gmcTest );
    M_gmc_right =  fusion::at_key<gmc1 >( _gmcTest );
    M_left_map = fusion::make_map<gmc<0> >( M_gmc_left );
    M_right_map = fusion::make_map<gmc<0> >( M_gmc_right );
#endif
    fusion::for_each( M_test_fec, vf::detail::FEContextUpdate<0, form_context_type>( _gmcTest, *this ) );
    M_test_fec0 = fusion::make_map<gmc<0>>( fusion::at_key<gmc<0>>( M_test_fec ) );
    M_test_fec1 = fusion::make_map<gmc1>( fusion::at_key<gmc1>( M_test_fec ) );
    M_eval0_expr->update( _gmcExpr, M_test_fec0 );
    M_eval1_expr->update( _gmcExpr, M_test_fec1 );

    M_integrator.update( *fusion::at_key<gmc<0>>( _gmcExpr ) );
}
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
void LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                                       map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                                       map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                                       IM const& im )
{
    M_integrator = im;
    this->update( _gmcTest, _gmcTrial, _gmcExpr );
}
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
void LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                                       map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                                       map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                                       IM const& im, mpl::int_<2> )
{
    M_integrator = im;
    this->update( _gmcTest, _gmcTrial, _gmcExpr, mpl::int_<2>() );
}

template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
void
    LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::integrate( mpl::int_<1> )
{
    typedef typename eval0_expr_type::shape shape;
    BOOST_MPL_ASSERT_MSG( ( shape::M == 1 && shape::N == 1 ),
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N>));

    if ( !UseMortar )
    {
        for ( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
        {
            M_rep( i ) = M_integrator( *M_eval0_expr, i, 0, 0 );
        }
    }
    else
    {
#if !defined( NDEBUG )
        CHECK( M_test_dof->mesh()->isBoundaryElement( M_gmc_left->id() ) ) << "element in context must be on boundary";
#endif
        for ( uint16_type i = 0; i < test_dof_type::nDofPerElement - 1; ++i )
        {
            M_rep_mortar( i ) = M_integrator( *M_eval0_expr, i, 0, 0 );
        }
    }
}
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
void
    LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::integrate( mpl::int_<2> )
{
    typedef mpl::int_<fusion::result_of::template size<map_test_geometric_mapping_context_type>::type::value> map_size;
    BOOST_MPL_ASSERT_MSG( map_size::value == 2, INVALID_GEOMAP, ( map_size, map_test_geometric_mapping_context_type ) );

    typedef typename eval0_expr_type::shape shape;
    BOOST_MPL_ASSERT_MSG( ( shape::M == 1 && shape::N == 1 ),
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N>));

    for ( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
    {
        uint16_type ii = i;
        // test dof element 0
        M_rep_2( ii ) = M_integrator( *M_eval0_expr, i, 0, 0 );

        ii = i + test_dof_type::nDofPerElement;
        // test dof element 1
        M_rep_2( ii ) = M_integrator( *M_eval1_expr, i, 0, 0 );
    }
}
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
void
    LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::integrateInCaseOfInterpolate( mpl::int_<1>,
                                                                                                                                                                            std::vector<boost::tuple<size_type, size_type>> const& indexLocalToQuad,
                                                                                                                                                                            bool isFirstExperience )
{
    typedef typename eval0_expr_type::shape shape;
    BOOST_MPL_ASSERT_MSG( ( shape::M == 1 && shape::N == 1 ),
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N>));

    if ( !UseMortar )
    {
        if ( isFirstExperience )
            for ( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
            {
                M_rep( i ) = M_integrator( *M_eval0_expr, i, 0, 0, indexLocalToQuad );
            }

        else
            for ( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
            {
                M_rep( i ) += M_integrator( *M_eval0_expr, i, 0, 0, indexLocalToQuad );
            }
    }
    else
    {
#if !defined( NDEBUG )
        CHECK( M_test_dof->mesh()->isBoundaryElement( M_gmc_left->id() ) ) << "element in context must be on boundary";
#endif
        if ( isFirstExperience )
            for ( uint16_type i = 0; i < test_dof_type::nDofPerElement - 1; ++i )
            {
                M_rep_mortar( i ) = M_integrator( *M_eval0_expr, i, 0, 0, indexLocalToQuad );
            }

        else
            for ( uint16_type i = 0; i < test_dof_type::nDofPerElement - 1; ++i )
            {
                M_rep_mortar( i ) += M_integrator( *M_eval0_expr, i, 0, 0, indexLocalToQuad );
            }
    }
}
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
void LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::assemble( size_type elt_0 )
{

    if ( !UseMortar )
    {
        M_local_rows = M_test_dof->localToGlobalIndices( elt_0, M_form.dofIdToContainerId() ).array();

#if !defined( NDEBUG )
        DVLOG( 2 ) << "lf dof (elt=" << elt_0 << ") = " << M_local_rows << "\n";
#endif
        if ( test_dof_type::is_modal || is_hdiv_conforming<test_fe_type>::value || is_hcurl_conforming<test_fe_type>::value )
        {
            M_local_rowsigns = M_test_dof->localToGlobalSigns( elt_0 );
            DVLOG( 2 ) << "rep = " << M_rep;
            M_rep.array() *= M_local_rowsigns.array().template cast<value_type>();
            DVLOG( 2 ) << "rep after sign change = " << M_rep;
        }
        M_form.addVector( M_local_rows.data(), M_local_rows.size(),
                          M_rep.data() );
    }
    else
    {
#if !defined( NDEBUG )
        CHECK( M_test_dof->mesh()->isBoundaryElement( elt_0 ) ) << "element in context must be on boundary";
#endif
        M_mortar_local_rows = M_test_dof->localToGlobalIndices( elt_0, M_form.dofIdToContainerId() ).array();
        if ( test_dof_type::is_modal || is_hdiv_conforming<test_fe_type>::value || is_hcurl_conforming<test_fe_type>::value )
        {
            CHECK( false ) << "TODO";
        }
        M_form.addVector( M_mortar_local_rows.data(), M_mortar_local_rows.size(),
                          M_rep_mortar.data() );
    }
}
template <typename SpaceType, typename VectorType, typename ElemContType>
template <typename GeomapContext, typename ExprT, typename IM, typename GeomapExprContext, typename GeomapTrialContext, int UseMortarType>
void LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext, ExprT, IM, GeomapExprContext, GeomapTrialContext, UseMortarType>::assemble( size_type elt_0, size_type elt_1 )
{
    M_local_rows_2.template head<test_dof_type::nDofPerElement>() = M_test_dof->localToGlobalIndices( elt_0, M_form.dofIdToContainerId() ).array();
    M_local_rows_2.template tail<test_dof_type::nDofPerElement>() = M_test_dof->localToGlobalIndices( elt_1, M_form.dofIdToContainerId() ).array();

    if ( test_dof_type::is_modal )
    {
        M_local_rowsigns_2.template head<test_dof_type::nDofPerElement>() = M_test_dof->localToGlobalSigns( elt_0 );
        M_local_rowsigns_2.template tail<test_dof_type::nDofPerElement>() = M_test_dof->localToGlobalSigns( elt_1 );

        M_rep_2.array() *= M_local_rowsigns_2.array().template cast<value_type>();
    }

    M_form.addVector( M_local_rows_2.data(), M_local_rows_2.size(),
                      M_rep_2.data() );
}
}
}
}
#endif /* __LinearFormContext_H */
