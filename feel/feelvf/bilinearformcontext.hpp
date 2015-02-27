/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
   \file bilinearformcontext.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-04-27
 */
#ifndef FEELPP_BILINEARFORMCONTEXT_HPP
#define FEELPP_BILINEARFORMCONTEXT_HPP 1

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
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::Context( form_type& __form,
        map_test_geometric_mapping_context_type const& _gmcTest,
        map_trial_geometric_mapping_context_type const& _gmcTrial,
        map_geometric_mapping_expr_context_type const & gmcExpr,
        ExprT const& expr,
        IM const& im )
    :
    M_form( __form ),
    M_lb( __form.blockList() ),
    M_test_dof( __form.testSpace()->dof().get() ),
    M_trial_dof( __form.trialSpace()->dof().get() ),


    M_test_pc( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), fusion::at_key<gmc<0> >( _gmcTest )->pc()->nodes() ) ),
    M_test_pc_face( precomputeTestBasisAtPoints( im ) ),
    M_trial_pc( new trial_precompute_type( M_form.trialFiniteElement<UseMortar>(), fusion::at_key<gmc<0> >( _gmcTrial )->pc()->nodes() ) ),
    M_trial_pc_face( precomputeTrialBasisAtPoints( im ) ),

    M_test_gmc( _gmcTest ),
    M_trial_gmc( _gmcTrial ),

    M_test_fec( fusion::transform( _gmcTest,
                                    vf::detail::FEContextInit<0,form_context_type>( __form.testFiniteElement<UseMortar>(),
                                            *this ) ) ),
    M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_test_fec ) ) ),
    M_trial_fec( getMap( M_test_fec,
                         fusion::transform( _gmcTrial,
                                            vf::detail::FEContextInit<1,form_context_type>( __form.trialFiniteElement<UseMortar>(),
                                                                                            *this ) ),
                         __form.testSpace()->mesh()->isRelatedTo( __form.trialSpace()->mesh() )) ),
    M_trial_fec0( getMapL( M_test_fec0, fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_trial_fec ) ),
                           __form.testSpace()->mesh()->isRelatedTo( __form.trialSpace()->mesh() )) ),

    M_rep(),
    M_rep_2(),
    M_eval_expr00( new eval00_expr_type( expr, gmcExpr, M_test_fec0, M_trial_fec0 ) ),

    M_eval_expr01(),
    M_eval_expr10(),
    M_eval_expr11(),

    M_integrator( im )
{
    this->initDynamicEigenMatrix();

    if ( UseMortar )
    {

        LOG(INFO) << "mortar phi: ndof " << M_test_pc->fePtr()->nbDof();
        M_test_pc->print();
        LOG(INFO) << "mortar Phi context ";
        fusion::at_key<gmc<0>>( M_test_fec0 ).get()->print();
    }
    M_eval_expr00->init( im );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
template<typename IM2>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::Context( form_type& __form,
        map_test_geometric_mapping_context_type const& _gmcTest,
        map_trial_geometric_mapping_context_type const& _gmcTrial,
        map_geometric_mapping_expr_context_type const & _gmcExpr,
        ExprT const& expr,
        IM const& im,
        IM2 const& im2 )
    :
    M_form( __form ),
    M_lb( __form.blockList() ),
    M_test_dof( __form.testSpace()->dof().get() ),
    M_trial_dof( __form.trialSpace()->dof().get() ),

    M_test_pc( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), im2.points() ) ),
    M_test_pc_face( precomputeTestBasisAtPoints( im2 ) ),
    M_trial_pc( new trial_precompute_type( M_form.trialFiniteElement<UseMortar>(), im2.points() ) ),
    M_trial_pc_face( precomputeTrialBasisAtPoints( im2 ) ),

    M_test_gmc( _gmcTest ),
    M_trial_gmc( _gmcTrial ),

    M_test_fec( fusion::transform( _gmcTest, vf::detail::FEContextInit<0,form_context_type>( __form.testFiniteElement<UseMortar>(), *this ) ) ),
    M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_test_fec ) ) ),
    M_trial_fec( getMap( M_test_fec,
                         fusion::transform( _gmcTrial, vf::detail::FEContextInit<1,form_context_type>( __form.trialFiniteElement<UseMortar>(), *this ) ),
                         __form.testSpace()->mesh()->isRelatedTo( __form.trialSpace()->mesh() )) ),
    M_trial_fec0( getMapL( M_test_fec0, fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_trial_fec ) ),
                           __form.testSpace()->mesh()->isRelatedTo( __form.trialSpace()->mesh() )) ),
    M_rep(),
    M_rep_2(),
    M_eval_expr00( new eval00_expr_type( expr, _gmcExpr, M_test_fec0, M_trial_fec0 ) ),
    M_eval_expr01(),
    M_eval_expr10(),
    M_eval_expr11(),
    M_integrator( im )
{
    this->initDynamicEigenMatrix();
    // faces
    M_eval_expr00->init( im2 );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
template<typename IM2,typename IMTest,typename IMTrial>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::Context( form_type& __form,
        map_test_geometric_mapping_context_type const& _gmcTest,
        map_trial_geometric_mapping_context_type const& _gmcTrial,
        map_geometric_mapping_expr_context_type const & _gmcExpr,
        ExprT const& expr,
        IM const& im, IM2 const& im2, IMTest const& imTest, IMTrial const& imTrial )
    :
    M_form( __form ),
    M_lb( __form.blockList() ),
    M_test_dof( __form.testSpace()->dof().get() ),
    M_trial_dof( __form.trialSpace()->dof().get() ),

    M_test_pc( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), fusion::at_key<gmc<0> >( _gmcTest )->pc()->nodes() ) ),
    M_test_pc_face( precomputeTestBasisAtPoints( imTest ) ),
    M_trial_pc( new trial_precompute_type( M_form.trialFiniteElement<UseMortar>(), fusion::at_key<gmc<0> >( _gmcTrial )->pc()->nodes() ) ),
    M_trial_pc_face( precomputeTrialBasisAtPoints( imTrial ) ),

    M_test_gmc( _gmcTest ),
    M_trial_gmc( _gmcTrial ),

    M_test_fec( fusion::transform( _gmcTest, vf::detail::FEContextInit<0,form_context_type>( __form.testFiniteElement<UseMortar>(), *this ) ) ),
    M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_test_fec ) ) ),
    M_trial_fec( getMap( M_test_fec, fusion::transform( _gmcTrial, vf::detail::FEContextInit<1,form_context_type>( __form.trialFiniteElement<UseMortar>(), *this ) ),
                         __form.testSpace()->mesh()->isRelatedTo( __form.trialSpace()->mesh() )) ),
    M_trial_fec0( getMapL( M_test_fec0, fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_trial_fec ) ),
                           __form.testSpace()->mesh()->isRelatedTo( __form.trialSpace()->mesh() )) ),
    M_rep(),
    M_rep_2(),
    M_eval_expr00( new eval00_expr_type( expr, _gmcExpr, M_test_fec0, M_trial_fec0 ) ),
    M_eval_expr01(),
    M_eval_expr10(),
    M_eval_expr11(),
    M_integrator( im )
{
    this->initDynamicEigenMatrix();
    // faces
    M_eval_expr00->init( im2 );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
template</*typename IM2, */typename IMTest,typename IMTrial>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::Context( form_type& __form,
        map_test_geometric_mapping_context_type const& _gmcTest,
        map_trial_geometric_mapping_context_type const& _gmcTrial,
        map_geometric_mapping_expr_context_type const & _gmcExpr,
        ExprT const& expr,
        IM const& im, IMTest const& imTest, IMTrial const& imTrial )
    :
    M_form( __form ),
    M_lb( __form.blockList() ),
    M_test_dof( __form.testSpace()->dof().get() ),
    M_trial_dof( __form.trialSpace()->dof().get() ),

    M_test_pc( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), fusion::at_key<gmc<0> >( _gmcTest )->pc()->nodes() ) ),
    M_test_pc_face( precomputeTestBasisAtPoints( imTest ) ),
    M_trial_pc( new trial_precompute_type( M_form.trialFiniteElement<UseMortar>(), fusion::at_key<gmc<0> >( _gmcTrial )->pc()->nodes() ) ),
    M_trial_pc_face( precomputeTrialBasisAtPoints( imTrial ) ),

    /*M_test_pc( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), im2.points() ) ),
    M_test_pc_face( precomputeTestBasisAtPoints( im2 ) ),
    M_trial_pc( new trial_precompute_type( M_form.trialFiniteElement<UseMortar>(), im2.points() ) ),
    M_trial_pc_face( precomputeTrialBasisAtPoints( im2 ) ),*/

    M_test_gmc( _gmcTest ),
    M_trial_gmc( _gmcTrial ),

    M_test_fec( fusion::transform( _gmcTest, vf::detail::FEContextInit<0,form_context_type>( __form.testFiniteElement<UseMortar>(), *this ) ) ),
    M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_test_fec ) ) ),
    M_trial_fec( getMap( M_test_fec, fusion::transform( _gmcTrial, vf::detail::FEContextInit<1,form_context_type>( __form.trialFiniteElement<UseMortar>(), *this ) ),
                         __form.testSpace()->mesh()->isRelatedTo( __form.trialSpace()->mesh() )) ),
    M_trial_fec0( getMapL( M_test_fec0, fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_trial_fec ) ),
                           __form.testSpace()->mesh()->isRelatedTo( __form.trialSpace()->mesh() )) ),
    M_rep(),
    M_rep_2(),
    M_eval_expr00( new eval00_expr_type( expr, _gmcExpr, M_test_fec0, M_trial_fec0 ) ),
    M_eval_expr01(),
    M_eval_expr10(),
    M_eval_expr11(),
    M_integrator( im )
{
    this->initDynamicEigenMatrix();
    // faces
    M_eval_expr00->init( im );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
template<typename IM2>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::Context( form_type& __form,
        map_test_geometric_mapping_context_type const& _gmcTest,
        map_trial_geometric_mapping_context_type const& _gmcTrial,
        map_geometric_mapping_expr_context_type const & _gmcExpr,
        ExprT const& expr,
        IM const& im,
        IM2 const& im2,
        mpl::int_<2> )
    :
    M_form( __form ),
    M_lb( __form.blockList() ),
    M_test_dof( __form.testSpace()->dof().get() ),
    M_trial_dof( __form.trialSpace()->dof().get() ),

    M_test_pc( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), im2.points() ) ),
    M_test_pc_face( precomputeTestBasisAtPoints( im2 ) ),
    M_trial_pc( new trial_precompute_type( M_form.trialFiniteElement<UseMortar>(), im2.points() ) ),
    M_trial_pc_face( precomputeTrialBasisAtPoints( im2 ) ),

    M_test_gmc( _gmcTest ),
    M_trial_gmc( _gmcTrial ),

    M_test_fec( fusion::transform( _gmcTest, vf::detail::FEContextInit<0,form_context_type>( __form.testFiniteElement<UseMortar>(), *this ) ) ),
    M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_test_fec ) ) ),
    M_test_fec1( fusion::make_map<test_gmc1 >( fusion::at_key<test_gmc1 >( M_test_fec ) ) ),
    M_trial_fec( fusion::transform( _gmcTrial, vf::detail::FEContextInit<1,form_context_type>( __form.trialFiniteElement<UseMortar>(), *this ) ) ),
    M_trial_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_trial_fec ) ) ),
    M_trial_fec1( fusion::make_map<trial_gmc1 >( fusion::at_key<trial_gmc1 >( M_trial_fec ) ) ),
    M_rep(),
    M_rep_2(),
    M_eval_expr00( new eval00_expr_type( expr, _gmcExpr, M_test_fec0, M_trial_fec0 ) ),
    M_eval_expr01( new eval01_expr_type( expr, _gmcExpr, M_test_fec0, M_trial_fec1 ) ),
    M_eval_expr10( new eval10_expr_type( expr, _gmcExpr, M_test_fec1, M_trial_fec0 ) ),
    M_eval_expr11( new eval11_expr_type( expr, _gmcExpr, M_test_fec1, M_trial_fec1 ) ),
    M_integrator( im )
{
    FEELPP_ASSERT( fusion::at_key<gmc<0> >( M_test_fec0 ).get() != 0 ).error( "invalid test_fec" );
    FEELPP_ASSERT( fusion::at_key<test_gmc1 >( M_test_fec1 ).get() != 0 ).error( "invalid test_fec" );
    FEELPP_ASSERT( fusion::at_key<gmc<0> >( M_trial_fec0 ).get() != 0 ).error( "invalid trial_fec" );
    FEELPP_ASSERT( fusion::at_key<trial_gmc1 >( M_trial_fec1 ).get() != 0 ).error( "invalid trial_fec" );

    this->initDynamicEigenMatrix();

    M_eval_expr00->init( im2 );
    M_eval_expr01->init( im2 );
    M_eval_expr10->init( im2 );
    M_eval_expr11->init( im2 );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                      map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                      map_geometric_mapping_expr_context_type const& _gmcExpr )
{
    if ( UseMortar )
        LOG(INFO) << "bilin context same fe: " << boost::is_same<map_test_fecontext_type, map_trial_fecontext_type>::value;
    update( _gmcTest, _gmcTrial, _gmcExpr, boost::is_same<map_test_fecontext_type, map_trial_fecontext_type>() );
    // if we know that the result will be zero, don't update the integrator and return immediately
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmcExpr ) );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                      map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                      map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                      mpl::bool_<false> )
{
    fusion::for_each( M_test_fec, vf::detail::FEContextUpdate<0,form_context_type>( _gmcTest, *this ) );
    M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_test_fec ) );

    if ( UseMortar )
    {
        LOG(INFO) << "update Fe mortar Phi context element : ";
        fusion::at_key<gmc<0>>( M_test_fec0 ).get()->print();
    }


    fusion::for_each( M_trial_fec, vf::detail::FEContextUpdate<1,form_context_type>( _gmcTrial, *this ) );
    M_trial_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_trial_fec ) );
    M_eval_expr00->update( _gmcExpr, M_test_fec0, M_trial_fec0 );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::update( map_test_geometric_mapping_context_type const& _gmcTest,
                                                                                                                      map_trial_geometric_mapping_context_type const& _gmcTrial,
                                                                                                                      map_geometric_mapping_expr_context_type const& _gmcExpr,
                                                                                                                      mpl::bool_<true> )
{
    fusion::for_each( M_test_fec, vf::detail::FEContextUpdate<0,form_context_type>( _gmcTest, *this ) );
    M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_test_fec ) );
    M_eval_expr00->update( _gmcExpr, M_test_fec0, M_test_fec0 );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::update( map_test_geometric_mapping_context_type const& _gmcTest,
        map_trial_geometric_mapping_context_type const& _gmcTrial,
        map_geometric_mapping_expr_context_type const& _gmcExpr,
        mpl::int_<2> )
{
    fusion::for_each( M_test_fec, vf::detail::FEContextUpdate<0,form_context_type>( _gmcTest, *this ) );
    M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_test_fec ) );
    M_test_fec1 = fusion::make_map<test_gmc1 >( fusion::at_key<test_gmc1 >( M_test_fec ) );
    fusion::for_each( M_trial_fec, vf::detail::FEContextUpdate<1,form_context_type>( _gmcTrial, *this ) );
    M_trial_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_trial_fec ) );
    M_trial_fec1 = fusion::make_map<trial_gmc1 >( fusion::at_key<trial_gmc1 >( M_trial_fec ) );

    FEELPP_ASSERT( fusion::at_key<gmc<0> >( M_test_fec0 ).get() != 0 )
    ( 0 ).error( "invalid test_fec0" );
    FEELPP_ASSERT( fusion::at_key<gmc<1> >( M_test_fec1 ).get() != 0 )
    ( 1 ).error( "invalid test_fec1" );
    FEELPP_ASSERT( fusion::at_key<gmc<0> >( M_trial_fec0 ).get() != 0 )
    ( 0 ).error( "invalid trial_fec0" );
    FEELPP_ASSERT( fusion::at_key<gmc<1> >( M_trial_fec1 ).get() != 0 )
    ( 0 ).error( "invalid trial_fec1" );

    M_eval_expr00->update( _gmcExpr, M_test_fec0, M_trial_fec0 );
    M_eval_expr01->update( _gmcExpr, M_test_fec0, M_trial_fec1 );
    M_eval_expr10->update( _gmcExpr, M_test_fec1, M_trial_fec0 );
    M_eval_expr11->update( _gmcExpr, M_test_fec1, M_trial_fec1 );

    M_integrator.update( *fusion::at_key<gmc<0> >( _gmcExpr ) );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& _gmcTest,
        map_trial_geometric_mapping_context_type const& _gmcTrial,
        map_geometric_mapping_expr_context_type const& _gmcExpr,
        std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad )
{
    M_test_gmc = _gmcTest;
    M_trial_gmc = _gmcTrial;
    precomputeBasisAtPoints( fusion::at_key<gmc<0> >( _gmcTest )->xRefs(),
                             fusion::at_key<gmc<0> >( _gmcTrial )->xRefs() );///!!!!!!!
    //updateInCaseOfInterpolate( _gmcTest, _gmcTrial, _gmcExpr, boost::is_same<map_test_fecontext_type, map_trial_fecontext_type>() );
    updateInCaseOfInterpolate( _gmcTest, _gmcTrial, _gmcExpr, mpl::bool_<false>() ); // forcage!
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmcExpr ), indexLocalToQuad );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& _gmcTest,
        map_trial_geometric_mapping_context_type const& _gmcTrial,
        map_geometric_mapping_expr_context_type const& _gmcExpr,
        mpl::bool_<false> )
{
    fusion::for_each( M_test_fec, vf::detail::FEContextUpdate<0,form_context_type>( _gmcTest, *this ) );
    M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_test_fec ) );
    fusion::for_each( M_trial_fec, vf::detail::FEContextUpdate<1,form_context_type>( _gmcTrial, *this ) );
    M_trial_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_trial_fec ) );
    M_eval_expr00->update( _gmcExpr, M_test_fec0, M_trial_fec0 );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& _gmcTest,
        map_trial_geometric_mapping_context_type const& _gmcTrial,
        map_geometric_mapping_expr_context_type const& _gmcExpr,
        mpl::bool_<true> )
{
    fusion::for_each( M_test_fec, vf::detail::FEContextUpdate<0,form_context_type>( _gmcTest, *this ) ); //!!!!!
    M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( M_test_fec ) );
    M_eval_expr00->update( _gmcExpr, M_test_fec0, M_test_fec0 );
}



template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::integrate( mpl::int_<1> )
{

    typedef test_geometric_mapping_context_type gmc_type;
    typedef typename eval00_expr_type::shape shape;
    static const bool cond = ( shape::M == 1 && shape::N == 1 );
    BOOST_MPL_ASSERT_MSG( cond,
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          ( mpl::int_<shape::M>, mpl::int_<shape::N> ) );

    test_geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( M_test_gmc );
#if !defined(NDEBUG)
    DVLOG(2) << "[BilinearForm::integrate] local assembly in element " << _gmc.id() << "\n";
#endif /* NDEBUG */

    if ( M_form.isPatternDefault() && boost::is_same<trial_dof_type,test_dof_type>::value &&
            trial_dof_type::is_product )
    {
        //if ( useEigenDynamicAlloc )
        M_rep = local_matrix_type::Zero(nDofPerElementTest, nDofPerElementTrial);
        //else
        //M_rep = local_matrix_type::Zero();

        if ( M_form.isPatternSymmetric() )
        {
            for ( uint16_type c = 0; c < trial_dof_type::nComponents1; ++c )
                for ( uint16_type j = 0; j < trial_dof_type::fe_type::nLocalDof; ++j )
                    for ( uint16_type i = 0; i <= j; ++i )
                    {
                        uint16_type testLocalDofIndex = i+c*test_dof_type::fe_type::nLocalDof;
                        uint16_type trialLocalDofIndex = j+c*trial_dof_type::fe_type::nLocalDof;
                        M_rep( testLocalDofIndex, trialLocalDofIndex ) = M_integrator( *M_eval_expr00, testLocalDofIndex, trialLocalDofIndex, 0, 0 );
                        M_rep( trialLocalDofIndex, testLocalDofIndex ) = M_rep( testLocalDofIndex, trialLocalDofIndex );
                    }
        }

        else
        {
            for ( uint16_type c = 0; c < trial_dof_type::nComponents1; ++c )
                for ( uint16_type j = 0; j < trial_dof_type::fe_type::nLocalDof; ++j )
                    for ( uint16_type i = 0; i < test_dof_type::fe_type::nLocalDof; ++i )
                    {
                        uint16_type testLocalDofIndex = i+c*test_dof_type::fe_type::nLocalDof;
                        uint16_type trialLocalDofIndex = j+c*trial_dof_type::fe_type::nLocalDof;
                        M_rep( testLocalDofIndex, trialLocalDofIndex ) = M_integrator( *M_eval_expr00, testLocalDofIndex, trialLocalDofIndex, 0, 0 );
                    }
        }
    }

    else
    {
        if ( boost::is_same<trial_dof_type,test_dof_type>::value && M_form.isPatternSymmetric()  )
        {
            for ( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
                for ( uint16_type i = 0; i <= j; ++i )
                {
                    M_rep( i, j ) = M_integrator( *M_eval_expr00, i, j, 0, 0 );
                    M_rep( j,i )=M_rep( i,j );
                }
        }

        else
        {
            DVLOG(2) << "local Assembly for element " << _gmc.id()
                     << " UseMortar=" << UseMortar << " bdy: " << M_test_dof->mesh()->isBoundaryElement( _gmc.id() );
            if ( !UseMortar || !M_test_dof->mesh()->isBoundaryElement( _gmc.id() ) )
                for ( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
                    for ( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
                    {
                        M_rep( i, j ) = M_integrator( *M_eval_expr00, i, j, 0, 0 );
                    }
            else
            {
                DVLOG(2) << "local Assembly for element " << _gmc.id()
                         << "ntestdof : " << test_dof_type::nDofPerElement-1;
                for ( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
                    for ( uint16_type i = 0; i < uint16_type(test_dof_type::nDofPerElement-1); ++i )
                    {
                        M_mortar_rep( i, j ) = M_integrator( *M_eval_expr00, i, j, 0, 0 );
                        DVLOG(2) << "mortar_rep(" << i << "," << j << ")=" << M_mortar_rep( i, j );

                    }
                DVLOG(2) << "local Assembly for element " << _gmc.id()
                         << "matrix = " << M_mortar_rep;
            }
        }
    }
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::integrate( mpl::int_<2> )
{
    //geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( M_gmc );
    typedef test_geometric_mapping_context_type gmc_type;
    typedef typename eval00_expr_type::shape shape;
    BOOST_MPL_ASSERT_MSG( ( mpl::and_<mpl::equal_to<mpl::int_<shape::M>,mpl::int_<1> >,
                            mpl::equal_to<mpl::int_<shape::N>,mpl::int_<1> > >::value ),
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          ( mpl::int_<shape::M>, mpl::int_<shape::N> ) );


    for ( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
        for ( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
        {
            uint16_type ii = i;
            uint16_type jj = j;
            // test dof element 0 - trial dof element 0
            M_rep_2( i, j ) = M_integrator( *M_eval_expr00, i, j, 0, 0 );

            ii = i;
            jj = j + trial_dof_type::nDofPerElement;
            // test dof element 0 - trial dof element 1
            M_rep_2( ii,jj ) = M_integrator( *M_eval_expr01, i, j, 0, 0 );

            ii = i + test_dof_type::nDofPerElement;
            jj = j;
            // test dof element 1 - trial dof element 0
            M_rep_2( ii,jj ) = M_integrator( *M_eval_expr10, i, j, 0, 0 );

            ii = i + test_dof_type::nDofPerElement;
            jj = j + trial_dof_type::nDofPerElement;
            // test dof element 1 - trial dof element 1
            M_rep_2( ii,jj ) = M_integrator( *M_eval_expr11, i, j, 0, 0 );
        }
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::integrateInCaseOfInterpolate( mpl::int_<1>,
        std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad,
        bool isFirstExperience )
{

    typedef test_geometric_mapping_context_type gmc_type;
    typedef typename eval00_expr_type::shape shape;
    static const bool cond = ( shape::M == 1 && shape::N == 1 );
    BOOST_MPL_ASSERT_MSG( cond,
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          ( mpl::int_<shape::M>, mpl::int_<shape::N> ) );


    test_geometric_mapping_context_type const& _gmcTest = *fusion::at_key<gmc<0> >( M_test_gmc );
    trial_geometric_mapping_context_type const& _gmcTrial = *fusion::at_key<gmc<0> >( M_trial_gmc );
    DVLOG(2) << "[BilinearForm::integrate] local assembly in element test " << _gmcTest.id() << " trial : " << _gmcTrial.id();

    if ( !UseMortar || !M_test_dof->mesh()->isBoundaryElement( _gmcTest.id() ) )
    {
        if ( isFirstExperience )
            for ( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
                for ( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
                {
                    M_rep( i, j ) = M_integrator( *M_eval_expr00, i, j, 0, 0, indexLocalToQuad );
                }

        else
            for ( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
                for ( uint16_type i = 0; i < test_dof_type::nDofPerElement; ++i )
                {
                    M_rep( i, j ) += M_integrator( *M_eval_expr00, i, j, 0, 0, indexLocalToQuad );
                }
        DVLOG(2) << "M_rep = " << M_rep*6/0.25;
    }
    else
    {
        if ( isFirstExperience )
            for ( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
                for ( uint16_type i = 0; i < uint16_type(test_dof_type::nDofPerElement-1); ++i )
                {
                    M_mortar_rep( i, j ) = M_integrator( *M_eval_expr00, i, j, 0, 0, indexLocalToQuad );
                }

        else
            for ( uint16_type j = 0; j < trial_dof_type::nDofPerElement; ++j )
                for ( uint16_type i = 0; i < uint16_type(test_dof_type::nDofPerElement-1); ++i )
                {
                    M_mortar_rep( i, j ) += M_integrator( *M_eval_expr00, i, j, 0, 0, indexLocalToQuad );
                }
        DVLOG(2) << "M_mortar_rep = " << M_mortar_rep*6/0.25;
    }

}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::assemble( size_type elt_0 )
{
    size_type row_start = M_lb.front().globalRowStart();
    size_type col_start = M_lb.front().globalColumnStart();

#if !defined(NDEBUG)
    DVLOG(2) << "[BilinearForm::assemble] global assembly in element " << elt_0 << "\n";
    DVLOG(2) << "[BilinearForm::assemble] row start " << row_start << "\n";
    DVLOG(2) << "[BilinearForm::assemble] col start " << col_start << "\n";
#endif /* NDEBUG */
    bool do_less = ( ( M_form.isPatternDefault() &&
                       ( M_test_dof->nComponents == M_trial_dof->nComponents ) ) &&
                     !M_form.isPatternCoupled() );

    if ( do_less )
    {
        for ( uint16_type c = 0; c < trial_dof_type::nComponents1; ++c )
        {
            M_c_rep = M_rep.block( c*test_dof_type::fe_type::nLocalDof, c*trial_dof_type::fe_type::nLocalDof,
                                    test_dof_type::fe_type::nLocalDof, trial_dof_type::fe_type::nLocalDof );
            M_c_local_rows.array() = M_test_dof->localToGlobalIndices( elt_0 ).array().segment( c*test_dof_type::fe_type::nLocalDof,
                                     test_dof_type::fe_type::nLocalDof ) + row_start;
            M_c_local_cols.array() = M_trial_dof->localToGlobalIndices( elt_0 ).array().segment( c*trial_dof_type::fe_type::nLocalDof,
                                     trial_dof_type::fe_type::nLocalDof ) + col_start;

            if ( test_dof_type::is_modal || trial_dof_type::is_modal )
            {
                M_c_local_rowsigns = M_test_dof->localToGlobalSigns( elt_0 ).segment( c*test_dof_type::fe_type::nLocalDof,test_dof_type::fe_type::nLocalDof );
                M_c_local_colsigns = M_trial_dof->localToGlobalSigns( elt_0 ).segment( c*trial_dof_type::fe_type::nLocalDof,trial_dof_type::fe_type::nLocalDof );
                M_c_rep.array() *= ( M_c_local_rowsigns*M_c_local_colsigns.transpose() ).array().template cast<value_type>();
            }

            M_form.addMatrix( M_c_local_rows.data(), M_c_local_rows.size(),
                               M_c_local_cols.data(), M_c_local_cols.size(),
                               M_c_rep.data() );
        }
    }

    else
    {
        size_type trial_eid= this->trialElementId( elt_0 );
        //
        DCHECK( trial_eid != invalid_size_type_value )
            << "this case should have been taken care of earlier before the assembly process\n";

        DVLOG(2) << "local Assembly for element " << elt_0
                 << " UseMortar=" << UseMortar << " bdy: " << M_test_dof->mesh()->isBoundaryElement( elt_0 );
        if ( !UseMortar || !M_test_dof->mesh()->isBoundaryElement( elt_0 ) )
        {
            M_local_rows.array() = M_test_dof->localToGlobalIndices( elt_0 ).array() + row_start;
            M_local_cols.array() = M_trial_dof->localToGlobalIndices( trial_eid ).array() + col_start;

            if ( test_dof_type::is_modal || trial_dof_type::is_modal ||
                 is_hdiv_conforming<trial_fe_type>::value || is_hdiv_conforming<test_fe_type>::value ||
                 is_hcurl_conforming<trial_fe_type>::value || is_hcurl_conforming<test_fe_type>::value )
            {
                M_local_rowsigns = M_test_dof->localToGlobalSigns( elt_0 );
                M_local_colsigns = M_trial_dof->localToGlobalSigns( trial_eid );
                M_rep.array() *= ( M_local_rowsigns*M_local_colsigns.transpose() ).array().template cast<value_type>();
            }

            M_form.addMatrix( M_local_rows.data(), M_local_rows.size(),
                              M_local_cols.data(), M_local_cols.size(),
                              M_rep.data() );
        }
        else
        {
            M_mortar_local_rows.array() = M_test_dof->localToGlobalIndices( elt_0 ).array() + row_start;
            M_local_cols.array() = M_trial_dof->localToGlobalIndices( trial_eid ).array() + col_start;

            if ( test_dof_type::is_modal || trial_dof_type::is_modal ||
                 is_hdiv_conforming<trial_fe_type>::value || is_hdiv_conforming<test_fe_type>::value ||
                 is_hcurl_conforming<trial_fe_type>::value || is_hcurl_conforming<test_fe_type>::value )
            {
                M_local_rowsigns = M_test_dof->localToGlobalSigns( elt_0 );
                M_local_colsigns = M_trial_dof->localToGlobalSigns( trial_eid );
                M_rep.array() *= ( M_local_rowsigns*M_local_colsigns.transpose() ).array().template cast<value_type>();
            }

            M_form.addMatrix( M_mortar_local_rows.data(), M_mortar_local_rows.size(),
                              M_local_cols.data(), M_local_cols.size(),
                              M_mortar_rep.data() );
        }
    }
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::assemble( size_type elt_0, size_type elt_1  )
{
    size_type row_start = M_lb.front().globalRowStart();
    size_type col_start = M_lb.front().globalColumnStart();

    size_type trial_e0id= this->trialElementId( elt_0 );
    size_type trial_e1id= this->trialElementId( elt_1 );
    DCHECK( trial_e0id != invalid_size_type_value && trial_e1id != invalid_size_type_value )
        << "this case should have been taken care of earlier before the assembly process\n";

    M_local_rows_2.template head<test_dof_type::nDofPerElement>().array() = M_test_dof->localToGlobalIndices( elt_0 ).array() + row_start;
    M_local_rows_2.template tail<test_dof_type::nDofPerElement>().array() = M_test_dof->localToGlobalIndices( elt_1 ).array() + row_start;

    M_local_cols_2.template head<trial_dof_type::nDofPerElement>().array() = M_trial_dof->localToGlobalIndices( trial_e0id ).array() + col_start;
    M_local_cols_2.template tail<trial_dof_type::nDofPerElement>().array() = M_trial_dof->localToGlobalIndices( trial_e1id ).array() + col_start;

    if ( test_dof_type::is_modal || trial_dof_type::is_modal )
    {
        M_local_rowsigns_2.template head<test_dof_type::nDofPerElement>() = M_test_dof->localToGlobalSigns( elt_0 );
        M_local_rowsigns_2.template tail<test_dof_type::nDofPerElement>() = M_test_dof->localToGlobalSigns( elt_1 );

        M_local_colsigns_2.template head<trial_dof_type::nDofPerElement>() = M_trial_dof->localToGlobalSigns( trial_e0id );
        M_local_colsigns_2.template tail<trial_dof_type::nDofPerElement>() = M_trial_dof->localToGlobalSigns( trial_e1id );

        M_rep_2.array() *= ( M_local_rowsigns_2*M_local_colsigns_2.transpose() ).array().template cast<value_type>();
    }

    M_form.addMatrix( M_local_rows_2.data(), M_local_rows_2.size(),
                       M_local_cols_2.data(), M_local_cols_2.size(),
                       M_rep_2.data() );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapTestContext,typename ExprT,typename IM,typename GeomapExprContext,typename GeomapTrialContext,bool UseMortar>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortar>::assembleInCaseOfInterpolate()
{
    size_type row_start = M_lb.front().globalRowStart();
    size_type col_start = M_lb.front().globalColumnStart();

    auto eltTest = fusion::at_key<gmc<0> >( M_test_gmc )->id();
    auto eltTrial = fusion::at_key<gmc<0> >( M_trial_gmc )->id();

    if ( !UseMortar || !M_test_dof->mesh()->isBoundaryElement( eltTest ) )
    {
        M_local_rows.array() = M_test_dof->localToGlobalIndices( eltTest ).array() + row_start;
        M_local_cols.array() = M_trial_dof->localToGlobalIndices( eltTrial ).array() + col_start;
        DVLOG(2) << "M_local_rows: " << M_local_rows;
        DVLOG(2) << "M_local_cols: " << M_local_cols;
#if 0
        bool do_less = ( ( M_form.isPatternDefault() &&
                           ( M_test_dof->nComponents == M_trial_dof->nComponents ) ) &&
                         !M_form.isPatternCoupled() );
#endif

        if ( test_dof_type::is_modal || trial_dof_type::is_modal )
        {
            M_local_rowsigns = M_test_dof->localToGlobalSigns( eltTest );
            M_local_colsigns = M_trial_dof->localToGlobalSigns( eltTrial );
            M_rep.array() *= ( M_local_rowsigns*M_local_colsigns.transpose() ).array().template cast<value_type>();
        }
        DVLOG(2) << "add rep : " << M_rep;
        M_form.addMatrix( M_local_rows.data(), M_local_rows.size(),
                          M_local_cols.data(), M_local_cols.size(),
                          M_rep.data() );

    }
    else
    {
        M_mortar_local_rows.array() = M_test_dof->localToGlobalIndices( eltTest ).array() + row_start;

        M_local_cols.array() = M_trial_dof->localToGlobalIndices( eltTrial ).array() + col_start;
        DVLOG(2) << "M_mortar_local_rows: " << M_mortar_local_rows;
        DVLOG(2) << "M_local_cols: " << M_local_cols;
#if 0
        bool do_less = ( ( M_form.isPatternDefault() &&
                           ( M_test_dof->nComponents == M_trial_dof->nComponents ) ) &&
                         !M_form.isPatternCoupled() );
#endif

        if ( test_dof_type::is_modal || trial_dof_type::is_modal )
        {
            M_local_rowsigns = M_test_dof->localToGlobalSigns( eltTest );
            M_local_colsigns = M_trial_dof->localToGlobalSigns( eltTrial );
            M_rep.array() *= ( M_local_rowsigns*M_local_colsigns.transpose() ).array().template cast<value_type>();
        }
        DVLOG(2) << "add mortar rep : " << M_mortar_rep;
        M_form.addMatrix( M_mortar_local_rows.data(), M_mortar_local_rows.size(),
                          M_local_cols.data(), M_local_cols.size(),
                          M_mortar_rep.data() );

    }
}


}
}
}
#endif /* __BilinearFormContext_H */
