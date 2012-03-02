/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \file bilinearformcontext.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-04-27
 */
#ifndef __BilinearFormContext_H
#define __BilinearFormContext_H 1

namespace Feel
{
namespace vf
{
namespace detail
{
//
// Context
//
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::Context( form_type& __form,
                                                                                                                       map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                       map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                       map_geometric_mapping_expr_context_type const & gmcExpr,
                                                                                                                       ExprT const& expr,
                                                                                                                       IM const& im )
    :
    _M_form( __form ),
    _M_lb( __form.blockList() ),
    _M_test_dof( __form.testSpace()->dof().get() ),
    _M_trial_dof( __form.trialSpace()->dof().get() ),
    _M_test_gmc( _gmcTest ),
    _M_trial_gmc( _gmcTrial ),

    _M_test_pc( new test_precompute_type( _M_form.testSpace()->fe(), fusion::at_key<gmc<0> >( _gmcTest )->pc()->nodes() ) ),
    _M_trial_pc( new trial_precompute_type( _M_form.trialSpace()->fe(), fusion::at_key<gmc<0> >( _gmcTrial )->pc()->nodes() ) ),
    _M_test_pc_face( precomputeTestBasisAtPoints( im ) ),
    _M_trial_pc_face( precomputeTrialBasisAtPoints( im ) ),

    _M_test_fec( fusion::transform( _gmcTest,
                                    detail::FEContextInit<0,form_context_type>(__form.testSpace()->fe(),
                                                                               *this ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_trial_fec( getMap( _M_test_fec, fusion::transform( _gmcTrial,
                                                          detail::FEContextInit<1,form_context_type>( __form.trialSpace()->fe(),
                                                                                                      *this ) ) ) ),
    _M_trial_fec0( getMapL( _M_test_fec0, fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) ) ) ),

    _M_rep(),
    _M_rep_2(),
    _M_eval_expr00( new eval00_expr_type( expr, gmcExpr, _M_test_fec0, _M_trial_fec0 ) ),

    _M_eval_expr01(),
    _M_eval_expr10(),
    _M_eval_expr11(),

    M_integrator( im )
{
    _M_eval_expr00->init( im );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
template<typename IM2>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::Context( form_type& __form,
                                                                                                                       map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                       map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                       map_geometric_mapping_expr_context_type const & _gmcExpr,
                                                                                                                       ExprT const& expr,
                                                                                                                       IM const& im,
                                                                                                                       IM2 const& im2 )
    :
    _M_form( __form ),
    _M_lb( __form.blockList() ),
    _M_test_dof( __form.testSpace()->dof().get() ),
    _M_trial_dof( __form.trialSpace()->dof().get() ),
    _M_test_gmc( _gmcTest ),
    _M_trial_gmc( _gmcTrial ),

    _M_test_pc( new test_precompute_type( _M_form.testSpace()->fe(), im2.points() ) ),
    _M_trial_pc( new trial_precompute_type( _M_form.trialSpace()->fe(), im2.points() ) ),
    _M_test_pc_face( precomputeTestBasisAtPoints( im2 ) ),
    _M_trial_pc_face( precomputeTrialBasisAtPoints( im2 ) ),

    _M_test_fec( fusion::transform( _gmcTest, detail::FEContextInit<0,form_context_type>(__form.testSpace()->fe(), *this ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_trial_fec( getMap( _M_test_fec, fusion::transform( _gmcTrial, detail::FEContextInit<1,form_context_type>( __form.trialSpace()->fe(), *this ) ) ) ),
    _M_trial_fec0( getMapL( _M_test_fec0, fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) ) ) ),
    _M_rep(),
    _M_rep_2(),
    _M_eval_expr00( new eval00_expr_type( expr, _gmcExpr, _M_test_fec0, _M_trial_fec0 ) ),
    _M_eval_expr01(),
    _M_eval_expr10(),
    _M_eval_expr11(),
    M_integrator( im )
{
    // faces
    _M_eval_expr00->init( im2 );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
template<typename IM2>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::Context( form_type& __form,
                                                                                                                       map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                       map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                       map_geometric_mapping_expr_context_type const & _gmcExpr,
                                                                                                                       ExprT const& expr,
                                                                                                                       IM const& im,
                                                                                                                       IM2 const& im2,
                                                                                                                       mpl::int_<2> )
    :
    _M_form( __form ),
    _M_lb( __form.blockList() ),
    _M_test_dof( __form.testSpace()->dof().get() ),
    _M_trial_dof( __form.trialSpace()->dof().get() ),
    _M_test_gmc( _gmcTest ),
    _M_trial_gmc( _gmcTrial ),

    _M_test_pc( new test_precompute_type( _M_form.testSpace()->fe(), im2.points() ) ),
    _M_trial_pc( new trial_precompute_type( _M_form.trialSpace()->fe(), im2.points() ) ),
    _M_test_pc_face( precomputeTestBasisAtPoints( im2 ) ),
    _M_trial_pc_face( precomputeTrialBasisAtPoints( im2 ) ),

    _M_test_fec( fusion::transform( _gmcTest, detail::FEContextInit<0,form_context_type>(__form.testSpace()->fe(), *this ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_test_fec1( fusion::make_map<test_gmc1 >( fusion::at_key<test_gmc1 >( _M_test_fec ) ) ),
    _M_trial_fec( fusion::transform( _gmcTrial, detail::FEContextInit<1,form_context_type>( __form.trialSpace()->fe(), *this ) ) ),
    _M_trial_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) ) ),
    _M_trial_fec1( fusion::make_map<trial_gmc1 >( fusion::at_key<trial_gmc1 >( _M_trial_fec ) ) ),
    _M_rep(),
    _M_rep_2(),
    _M_eval_expr00( new eval00_expr_type( expr, _gmcExpr, _M_test_fec0, _M_trial_fec0 ) ),
    _M_eval_expr01( new eval01_expr_type( expr, _gmcExpr, _M_test_fec0, _M_trial_fec1 ) ),
    _M_eval_expr10( new eval10_expr_type( expr, _gmcExpr, _M_test_fec1, _M_trial_fec0 ) ),
    _M_eval_expr11( new eval11_expr_type( expr, _gmcExpr, _M_test_fec1, _M_trial_fec1 ) ),
    M_integrator( im )
{
    FEELPP_ASSERT( fusion::at_key<gmc<0> >( _M_test_fec0 ).get() != 0 ).error( "invalid test_fec");
    FEELPP_ASSERT( fusion::at_key<test_gmc1 >( _M_test_fec1 ).get() != 0 ).error( "invalid test_fec");
    FEELPP_ASSERT( fusion::at_key<gmc<0> >( _M_trial_fec0 ).get() != 0 ).error( "invalid trial_fec");
    FEELPP_ASSERT( fusion::at_key<trial_gmc1 >( _M_trial_fec1 ).get() != 0 ).error( "invalid trial_fec");

    _M_eval_expr00->init( im2 );
    _M_eval_expr01->init( im2 );
    _M_eval_expr10->init( im2 );
    _M_eval_expr11->init( im2 );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                      map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                      map_geometric_mapping_expr_context_type const& _gmcExpr)
{
    update( _gmcTest, _gmcTrial, _gmcExpr, boost::is_same<map_test_fecontext_type, map_trial_fecontext_type>() );
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmcExpr ) );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                      map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                      map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                      mpl::bool_<false> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmcTest, *this ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    fusion::for_each( _M_trial_fec, detail::FEContextUpdate<1,form_context_type>( _gmcTrial, *this ) );
    _M_trial_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) );
    _M_eval_expr00->update( _gmcExpr, _M_test_fec0, _M_trial_fec0 );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                      map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                      map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                      mpl::bool_<true> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmcTest, *this ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_eval_expr00->update( _gmcExpr, _M_test_fec0, _M_test_fec0 );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                      map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                      map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                      mpl::int_<2> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmcTest, *this ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_test_fec1 = fusion::make_map<test_gmc1 >( fusion::at_key<test_gmc1 >( _M_test_fec ) );
    fusion::for_each( _M_trial_fec, detail::FEContextUpdate<1,form_context_type>( _gmcTrial, *this ) );
    _M_trial_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) );
    _M_trial_fec1 = fusion::make_map<trial_gmc1 >( fusion::at_key<trial_gmc1 >( _M_trial_fec ) );

    FEELPP_ASSERT( fusion::at_key<gmc<0> >( _M_test_fec0 ).get() != 0 )
        ( 0 ).error( "invalid test_fec0" );
    FEELPP_ASSERT( fusion::at_key<gmc<1> >( _M_test_fec1 ).get() != 0 )
        ( 1 ).error( "invalid test_fec1" );
    FEELPP_ASSERT( fusion::at_key<gmc<0> >( _M_trial_fec0 ).get() != 0 )
        ( 0 ).error( "invalid trial_fec0" );
    FEELPP_ASSERT( fusion::at_key<gmc<1> >( _M_trial_fec1 ).get() != 0 )
        ( 0 ).error( "invalid trial_fec1" );

    _M_eval_expr00->update( _gmcExpr, _M_test_fec0, _M_trial_fec0 );
    _M_eval_expr01->update( _gmcExpr, _M_test_fec0, _M_trial_fec1 );
    _M_eval_expr10->update( _gmcExpr, _M_test_fec1, _M_trial_fec0 );
    _M_eval_expr11->update( _gmcExpr, _M_test_fec1, _M_trial_fec1 );

    M_integrator.update( *fusion::at_key<gmc<0> >( _gmcTest ) );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                         map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                         map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                         std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad )
{
    precomputeBasisAtPoints( fusion::at_key<gmc<0> >( _gmcTest )->xRefs(),
                             fusion::at_key<gmc<0> >( _gmcTrial )->xRefs() );///!!!!!!!
    //updateInCaseOfInterpolate( _gmcTest, _gmcTrial, _gmcExpr, boost::is_same<map_test_fecontext_type, map_trial_fecontext_type>() );
    updateInCaseOfInterpolate( _gmcTest, _gmcTrial, _gmcExpr, mpl::bool_<false>() ); // forcage!
    M_integrator.update(*fusion::at_key<gmc<0> >( _gmcExpr ), indexLocalToQuad );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                     map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                     map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                     mpl::bool_<false> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdateInCaseOfInterpolate<0,form_context_type>( _gmcTest, *this ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    fusion::for_each( _M_trial_fec, detail::FEContextUpdateInCaseOfInterpolate<1,form_context_type>( _gmcTrial, *this ) );
    _M_trial_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) );
    _M_eval_expr00->update( _gmcExpr, _M_test_fec0, _M_trial_fec0 );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                                         map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                                         map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                                         mpl::bool_<true> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdateInCaseOfInterpolate<0,form_context_type>( _gmcTest, *this ) ); //!!!!!
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_eval_expr00->update( _gmcExpr, _M_test_fec0, _M_test_fec0 );
}



template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::integrate( mpl::int_<1> )
{

    typedef test_geometric_mapping_context_type gmc_type;
    typedef typename eval00_expr_type::shape shape;
    static const bool cond = (shape::M == 1 && shape::N == 1);
    BOOST_MPL_ASSERT_MSG( cond,
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N> ) );

#if !defined(NDEBUG)
    test_geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_test_gmc );
    Debug( 5050 ) << "[BilinearForm::integrate] local assembly in element " << _gmc.id() << "\n";
#endif /* NDEBUG */

    if (_M_form.isPatternDefault() && boost::is_same<trial_dof_type,test_dof_type>::value &&
        trial_dof_type::is_product)
    {
        _M_rep = local_matrix_type::Zero();
        if ( _M_form.isPatternSymmetric() )
        {
            for( uint16_type c = 0; c < trial_dof_type::nComponents1; ++c )
                for( uint16_type j = 0; j < trial_dof_type::fe_type::nLocalDof; ++j )
                    for( uint16_type i = 0; i <= j; ++i )
                    {
                        uint16_type testLocalDofIndex = i+c*test_dof_type::fe_type::nLocalDof;
                        uint16_type trialLocalDofIndex = j+c*trial_dof_type::fe_type::nLocalDof;
                        _M_rep(testLocalDofIndex, trialLocalDofIndex ) = M_integrator( *_M_eval_expr00, testLocalDofIndex, trialLocalDofIndex, 0, 0 );
                        _M_rep(trialLocalDofIndex, testLocalDofIndex ) = _M_rep(testLocalDofIndex, trialLocalDofIndex );
                    }
        }
        else
        {
            for( uint16_type c = 0; c < trial_dof_type::nComponents1; ++c )
                for( uint16_type j = 0; j < trial_dof_type::fe_type::nLocalDof; ++j )
                    for( uint16_type i = 0; i < test_dof_type::fe_type::nLocalDof; ++i )
                    {
                        uint16_type testLocalDofIndex = i+c*test_dof_type::fe_type::nLocalDof;
                        uint16_type trialLocalDofIndex = j+c*trial_dof_type::fe_type::nLocalDof;
                        _M_rep(testLocalDofIndex, trialLocalDofIndex ) = M_integrator( *_M_eval_expr00, testLocalDofIndex, trialLocalDofIndex, 0, 0 );
                    }
        }
    }
    else
    {
        if ( boost::is_same<trial_dof_type,test_dof_type>::value && _M_form.isPatternSymmetric()  )
        {
            for( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
                for( uint16_type i = 0; i <= j; ++i )
                {
                    _M_rep(i, j ) = M_integrator( *_M_eval_expr00, i, j, 0, 0 );
                    _M_rep(j,i)=_M_rep(i,j);
                }
        }
        else
        {
            for( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
                for( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
                {
                    _M_rep(i, j ) = M_integrator( *_M_eval_expr00, i, j, 0, 0 );
                }
        }
    }
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::integrate( mpl::int_<2> )
{
    //geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_gmc );
    typedef test_geometric_mapping_context_type gmc_type;
    typedef typename eval00_expr_type::shape shape;
    BOOST_MPL_ASSERT_MSG( (mpl::and_<mpl::equal_to<mpl::int_<shape::M>,mpl::int_<1> >,
                                     mpl::equal_to<mpl::int_<shape::N>,mpl::int_<1> > >::value),
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N> ) );


    for( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
        for( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
        {
            uint16_type ii = i;
            uint16_type jj = j;
            // test dof element 0 - trial dof element 0
            _M_rep_2(i, j ) = M_integrator( *_M_eval_expr00, i, j, 0, 0 );

            ii = i;
            jj = j + trial_dof_type::nDofPerElement;
            // test dof element 0 - trial dof element 1
            _M_rep_2(ii,jj) = M_integrator( *_M_eval_expr01, i, j, 0, 0 );

            ii = i + test_dof_type::nDofPerElement;
            jj = j;
            // test dof element 1 - trial dof element 0
            _M_rep_2(ii,jj) = M_integrator( *_M_eval_expr10, i, j, 0, 0 );

            ii = i + test_dof_type::nDofPerElement;
            jj = j + trial_dof_type::nDofPerElement;
            // test dof element 1 - trial dof element 1
            _M_rep_2(ii,jj) = M_integrator( *_M_eval_expr11, i, j, 0, 0 );
        }
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::integrateInCaseOfInterpolate( mpl::int_<1>,
                                                                                                                                            std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad,
                                                                                                                                            bool isFirstExperience )
{

    typedef test_geometric_mapping_context_type gmc_type;
    typedef typename eval00_expr_type::shape shape;
    static const bool cond = (shape::M == 1 && shape::N == 1);
    BOOST_MPL_ASSERT_MSG( cond,
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N> ) );

#if !defined(NDEBUG)
    test_geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_test_gmc );
    Debug( 5050 ) << "[BilinearForm::integrate] local assembly in element " << _gmc.id() << "\n";
#endif /* NDEBUG */
    if (isFirstExperience)
        for( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
            for( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
            {
                _M_rep(i, j ) = M_integrator( *_M_eval_expr00, i, j, 0, 0, indexLocalToQuad );
            }
    else
        for( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
            for( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
            {
                _M_rep(i, j ) += M_integrator( *_M_eval_expr00, i, j, 0, 0, indexLocalToQuad );
            }
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::assemble( size_type elt_0 )
{
    size_type row_start = _M_lb.front().globalRowStart();
    size_type col_start = _M_lb.front().globalColumnStart();

#if !defined(NDEBUG)
    Debug( 5050 ) << "[BilinearForm::assemble] global assembly in element " << elt_0 << "\n";
    Debug( 5050 ) << "[BilinearForm::assemble] row start " << row_start << "\n";
    Debug( 5050 ) << "[BilinearForm::assemble] col start " << col_start << "\n";
#endif /* NDEBUG */
    bool do_less = ( ( _M_form.isPatternDefault() &&
                       ( _M_test_dof->nComponents == _M_trial_dof->nComponents ) ) &&
                     !_M_form.isPatternCoupled() );

    if ( do_less )
        //if ( 0 )
    {
        for( uint16_type c = 0; c < trial_dof_type::nComponents1; ++c )
        {
            M_c_rep = _M_rep.block( c*test_dof_type::fe_type::nLocalDof, c*trial_dof_type::fe_type::nLocalDof,
                                    test_dof_type::fe_type::nLocalDof, trial_dof_type::fe_type::nLocalDof );
            M_c_local_rows.array() = _M_test_dof->localToGlobalIndices(elt_0).array().segment(c*test_dof_type::fe_type::nLocalDof,
                                                                                              test_dof_type::fe_type::nLocalDof) + row_start;
            M_c_local_cols.array() = _M_trial_dof->localToGlobalIndices(elt_0).array().segment(c*trial_dof_type::fe_type::nLocalDof,
                                                                                               trial_dof_type::fe_type::nLocalDof) + col_start;

            if ( test_dof_type::is_modal || trial_dof_type::is_modal )
            {
                M_c_local_rowsigns = _M_test_dof->localToGlobalSigns(elt_0).segment(c*test_dof_type::fe_type::nLocalDof,test_dof_type::fe_type::nLocalDof);
                M_c_local_colsigns = _M_trial_dof->localToGlobalSigns(elt_0).segment(c*trial_dof_type::fe_type::nLocalDof,trial_dof_type::fe_type::nLocalDof);
                M_c_rep.array() *= (M_c_local_rowsigns*M_c_local_colsigns.transpose()).array().template cast<value_type>();
            }

            _M_form.addMatrix( M_c_local_rows.data(), M_c_local_rows.size(),
                               M_c_local_cols.data(), M_c_local_cols.size(),
                               M_c_rep.data() );
        }
    }
    else
    {
        M_local_rows.array() = _M_test_dof->localToGlobalIndices(elt_0).array() + row_start;
        M_local_cols.array() = _M_trial_dof->localToGlobalIndices(elt_0).array() + col_start;


        if ( test_dof_type::is_modal || trial_dof_type::is_modal )
        {
            M_local_rowsigns = _M_test_dof->localToGlobalSigns(elt_0);
            M_local_colsigns = _M_trial_dof->localToGlobalSigns(elt_0);
            _M_rep.array() *= (M_local_rowsigns*M_local_colsigns.transpose()).array().template cast<value_type>();
        }

        _M_form.addMatrix( M_local_rows.data(), M_local_rows.size(),
                           M_local_cols.data(), M_local_cols.size(),
                           _M_rep.data() );
    }
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::assemble( size_type elt_0, size_type elt_1  )
{
    size_type row_start = _M_lb.front().globalRowStart();
    size_type col_start = _M_lb.front().globalColumnStart();

    auto local_rows_0 = _M_test_dof->localToGlobalIndices(elt_0).array() + row_start;
    auto local_rows_1 = _M_test_dof->localToGlobalIndices(elt_1).array() + row_start;
    M_local_rows_2.template head<test_dof_type::nDofPerElement>().array() = local_rows_0;
    M_local_rows_2.template tail<test_dof_type::nDofPerElement>().array() = local_rows_1;

    auto local_cols_0 = _M_trial_dof->localToGlobalIndices(elt_0).array() + col_start;
    auto local_cols_1 = _M_trial_dof->localToGlobalIndices(elt_1).array() + col_start;
    M_local_cols_2.template head<trial_dof_type::nDofPerElement>().array() = local_cols_0;
    M_local_cols_2.template tail<trial_dof_type::nDofPerElement>().array() = local_cols_1;


    if ( test_dof_type::is_modal || trial_dof_type::is_modal )
    {
        auto local_rowsigns_0 = _M_test_dof->localToGlobalSigns(elt_0);
        auto local_rowsigns_1 = _M_test_dof->localToGlobalSigns(elt_1);
        M_local_rowsigns_2.template head<test_dof_type::nDofPerElement>() = local_rowsigns_0;
        M_local_rowsigns_2.template tail<test_dof_type::nDofPerElement>() = local_rowsigns_1;

        auto local_colsigns_0 = _M_trial_dof->localToGlobalSigns(elt_0);
        auto local_colsigns_1 = _M_trial_dof->localToGlobalSigns(elt_1);
        M_local_colsigns_2.template head<trial_dof_type::nDofPerElement>() = local_colsigns_0;
        M_local_colsigns_2.template tail<trial_dof_type::nDofPerElement>() = local_colsigns_1;

        _M_rep_2.array() *= (M_local_rowsigns_2*M_local_colsigns_2.transpose()).array().template cast<value_type>();
    }
    _M_form.addMatrix( M_local_rows_2.data(), M_local_rows_2.size(),
                       M_local_cols_2.data(), M_local_cols_2.size(),
                       _M_rep_2.data() );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext>::assembleInCaseOfInterpolate()
{
    size_type row_start = _M_lb.front().globalRowStart();
    size_type col_start = _M_lb.front().globalColumnStart();

    auto eltTest = fusion::at_key<gmc<0> >( _M_test_gmc )->id();
    auto eltTrial = fusion::at_key<gmc<0> >( _M_trial_gmc )->id();

    M_local_rows.array() = _M_test_dof->localToGlobalIndices(eltTest).array() + row_start;
    M_local_cols.array() = _M_trial_dof->localToGlobalIndices(eltTrial).array() + col_start;

    bool do_less = ( ( _M_form.isPatternDefault() &&
                       ( _M_test_dof->nComponents == _M_trial_dof->nComponents ) ) &&
                     !_M_form.isPatternCoupled() );

    if ( test_dof_type::is_modal || trial_dof_type::is_modal )
    {
        M_local_rowsigns = _M_test_dof->localToGlobalSigns(eltTest);
        M_local_colsigns = _M_trial_dof->localToGlobalSigns(eltTrial);
        _M_rep.array() *= (M_local_rowsigns*M_local_colsigns.transpose()).array().template cast<value_type>();
    }

    _M_form.addMatrix( M_local_rows.data(), M_local_rows.size(),
                       M_local_cols.data(), M_local_cols.size(),
                       _M_rep.data() );

}


}
}
    }
#endif /* __BilinearFormContext_H */
