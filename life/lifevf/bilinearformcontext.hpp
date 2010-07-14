/* -*- mode: c++ -*-

  This file is part of the Life library

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

namespace Life
{
namespace vf
{
namespace detail
{
//
// Context
//
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                              map_geometric_mapping_context_type const& _gmc,
                                                                              ExprT const& expr,
                                                                              IM const& im )
    :
    _M_form( __form ),
    _M_lb( __form.blockList() ),
    _M_test_dof( __form.testSpace()->dof().get() ),
    _M_trial_dof( __form.trialSpace()->dof().get() ),
    _M_gmc( _gmc ),
    _M_test_fec( fusion::transform( _gmc,
                                    detail::FEContextInit<0,form_context_type>(__form.testSpace()->fe(),
                                                                       _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_trial_fec( getMap( _M_test_fec, fusion::transform( _gmc,
                                                          detail::FEContextInit<1,form_context_type>( __form.trialSpace()->fe(),
                                                                                              _M_form ) ) ) ),
    _M_trial_fec0( getMapL( _M_test_fec0, fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) ) ) ),
    _M_rep( boost::extents[fusion::size(_gmc)*test_fecontext_type::nDof][test_fecontext_type::nComponents1]
            [fusion::size(_gmc)*trial_fecontext_type::nDof][trial_fecontext_type::nComponents1] ),
    _M_eval_expr00( new eval00_expr_type( expr, _gmc, _M_test_fec0, _M_trial_fec0 ) ),
    _M_eval_expr01(),
    _M_eval_expr10(),
    _M_eval_expr11(),
    M_integrator( im )
{
    _M_eval_expr00->init( im );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
template<typename IM2>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                              map_geometric_mapping_context_type const& _gmc,
                                                                              ExprT const& expr,
                                                                              IM const& im,
                                                                              IM2 const& im2 )
    :
    _M_form( __form ),
    _M_lb( __form.blockList() ),
    _M_test_dof( __form.testSpace()->dof().get() ),
    _M_trial_dof( __form.trialSpace()->dof().get() ),
    _M_gmc( _gmc ),
    _M_test_fec( fusion::transform( _gmc, detail::FEContextInit<0,form_context_type>(__form.testSpace()->fe(), _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_trial_fec( getMap( _M_test_fec, fusion::transform( _gmc, detail::FEContextInit<1,form_context_type>( __form.trialSpace()->fe(), _M_form ) ) ) ),
    _M_trial_fec0( getMapL( _M_test_fec0, fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) ) ) ),
    _M_rep( boost::extents[fusion::size(_gmc)*test_fecontext_type::nDof][test_fecontext_type::nComponents1]
            [fusion::size(_gmc)*trial_fecontext_type::nDof][trial_fecontext_type::nComponents1] ),
    _M_eval_expr00( new eval00_expr_type( expr, _gmc, _M_test_fec0, _M_trial_fec0 ) ),
    _M_eval_expr01(),
    _M_eval_expr10(),
    _M_eval_expr11(),
    M_integrator( im )
{
    // faces
    _M_eval_expr00->init( im2 );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
template<typename IM2>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                                         map_geometric_mapping_context_type const& _gmc,
                                                                                         ExprT const& expr,
                                                                                         IM const& im,
                                                                                         IM2 const& im2,
                                                                                         mpl::int_<2> )
    :
    _M_form( __form ),
    _M_lb( __form.blockList() ),
    _M_test_dof( __form.testSpace()->dof().get() ),
    _M_trial_dof( __form.trialSpace()->dof().get() ),
    _M_gmc( _gmc ),
    _M_test_fec( fusion::transform( _gmc, detail::FEContextInit<0,form_context_type>(__form.testSpace()->fe(), _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_test_fec1( fusion::make_map<gmc1 >( fusion::at_key<gmc1 >( _M_test_fec ) ) ),
    _M_trial_fec( fusion::transform( _gmc, detail::FEContextInit<1,form_context_type>( __form.trialSpace()->fe(), _M_form ) ) ),
    _M_trial_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) ) ),
    _M_trial_fec1( fusion::make_map<gmc1 >( fusion::at_key<gmc1 >( _M_trial_fec ) ) ),
    _M_rep( boost::extents[fusion::size(_gmc)*test_fecontext_type::nDof][test_fecontext_type::nComponents1]
            [fusion::size(_gmc)*trial_fecontext_type::nDof][trial_fecontext_type::nComponents1]),
    _M_eval_expr00( new eval00_expr_type( expr, _gmc, _M_test_fec0, _M_trial_fec0 ) ),
    _M_eval_expr01( new eval01_expr_type( expr, _gmc, _M_test_fec0, _M_trial_fec1 ) ),
    _M_eval_expr10( new eval10_expr_type( expr, _gmc, _M_test_fec1, _M_trial_fec0 ) ),
    _M_eval_expr11( new eval11_expr_type( expr, _gmc, _M_test_fec1, _M_trial_fec1 ) ),
    M_integrator( im )
{
    LIFE_ASSERT( fusion::at_key<gmc<0> >( _M_test_fec0 ).get() != 0 ).error( "invalid test_fec");
    LIFE_ASSERT( fusion::at_key<gmc1 >( _M_test_fec1 ).get() != 0 ).error( "invalid test_fec");
    LIFE_ASSERT( fusion::at_key<gmc<0> >( _M_trial_fec0 ).get() != 0 ).error( "invalid trial_fec");
    LIFE_ASSERT( fusion::at_key<gmc1 >( _M_trial_fec1 ).get() != 0 ).error( "invalid trial_fec");

    _M_eval_expr00->init( im2 );
    _M_eval_expr01->init( im2 );
    _M_eval_expr10->init( im2 );
    _M_eval_expr11->init( im2 );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc )
{
    update( _gmc,  boost::is_same<map_test_fecontext_type, map_trial_fecontext_type>() );
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmc ) );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, mpl::bool_<false> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmc, _M_form ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    fusion::for_each( _M_trial_fec, detail::FEContextUpdate<1,form_context_type>( _gmc, _M_form ) );
    _M_trial_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) );
    _M_eval_expr00->update( _gmc, _M_test_fec0, _M_trial_fec0 );
    std::fill( _M_rep.data(), _M_rep.data()+_M_rep.num_elements(), value_type(0) );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, mpl::bool_<true> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmc, _M_form ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_eval_expr00->update( _gmc, _M_test_fec0, _M_test_fec0 );
    std::fill( _M_rep.data(), _M_rep.data()+_M_rep.num_elements(), value_type(0) );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, mpl::int_<2> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmc, _M_form ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_test_fec1 = fusion::make_map<gmc1 >( fusion::at_key<gmc1 >( _M_test_fec ) );
    fusion::for_each( _M_trial_fec, detail::FEContextUpdate<1,form_context_type>( _gmc, _M_form ) );
    _M_trial_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) );
    _M_trial_fec1 = fusion::make_map<gmc1 >( fusion::at_key<gmc1 >( _M_trial_fec ) );

    LIFE_ASSERT( fusion::at_key<gmc<0> >( _M_test_fec0 ).get() != 0 )
        ( 0 ).error( "invalid test_fec0" );
    LIFE_ASSERT( fusion::at_key<gmc1 >( _M_test_fec1 ).get() != 0 )
        ( 1 ).error( "invalid test_fec1" );
    LIFE_ASSERT( fusion::at_key<gmc<0> >( _M_trial_fec0 ).get() != 0 )
        ( 0 ).error( "invalid trial_fec0" );
    LIFE_ASSERT( fusion::at_key<gmc1 >( _M_trial_fec1 ).get() != 0 )
        ( 0 ).error( "invalid trial_fec1" );

    _M_eval_expr00->update( _gmc, _M_test_fec0, _M_trial_fec0 );
    _M_eval_expr01->update( _gmc, _M_test_fec0, _M_trial_fec1 );
    _M_eval_expr10->update( _gmc, _M_test_fec1, _M_trial_fec0 );
    _M_eval_expr11->update( _gmc, _M_test_fec1, _M_trial_fec1 );

    std::fill( _M_rep.data(), _M_rep.data()+_M_rep.num_elements(), value_type(0) );
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmc ) );
}



template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::integrate( mpl::int_<1> )
{

    typedef geometric_mapping_context_type gmc_type;
    typedef typename eval00_expr_type::shape shape;
    static const bool cond = (shape::M == 1 && shape::N == 1);
    BOOST_MPL_ASSERT_MSG( cond,
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N> ) );

#if !defined(NDEBUG)
    geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_gmc );
    Debug( 5050 ) << "[BilinearForm::integrate] local assembly in element " << _gmc.id() << "\n";
#endif /* NDEBUG */

    if ( ( _M_form.isPatternDefault() &&
           ( _M_test_dof->nComponents == _M_trial_dof->nComponents ) ) &&
         !_M_form.isPatternCoupled() )
    {
        for ( uint16_type i = 0; i < test_fecontext_type::nDof; ++i )
        {
            int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
            for ( uint16_type c1 = 0; c1 < ncomp; ++c1 )
            {
                indi.setIndex( boost::make_tuple( i, c1, 0 ) );
                LIFE_ASSERT( size_type(indi.index()) ==
                             size_type(c1*test_fecontext_type::nDof+i) )
                    ( i )( c1 )
                    ( indi.index() )
                    ( c1*test_fecontext_type::nDof+i ).error( "invalid test index" );

                for ( uint16_type j = 0; j < trial_fecontext_type::nDof; ++j )
                {
                        indj.setIndex( boost::make_tuple( j, c1, 0 ) );
                        LIFE_ASSERT( size_type(indj.index()) ==
                                     size_type(c1*trial_fecontext_type::nDof+j) )
                            ( j )( c1 )
                            ( indj.index() )
                            ( c1*trial_fecontext_type::nDof+j ).error( "invalid trial index" );

                        _M_rep[i][c1][j][c1] = M_integrator( *_M_eval_expr00,
                                                             indi,
                                                             indj,
                                                             0, 0 );
                    } // j, c2
            }
        }
    }
    else
    {
        for ( uint16_type i = 0; i < test_fecontext_type::nDof; ++i )
        {
            int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
            for ( uint16_type c1 = 0; c1 < ncomp; ++c1 )
            {
                indi.setIndex( boost::make_tuple( i, c1, 0 ) );
                LIFE_ASSERT( size_type(indi.index()) ==
                             size_type(c1*test_fecontext_type::nDof+i) )
                    ( i )( c1 )
                    ( indi.index() )
                    ( c1*test_fecontext_type::nDof+i ).error( "invalid test index" );

                for ( uint16_type j = 0; j < trial_fecontext_type::nDof; ++j )
                {
                    int ncomp2= (trial_fecontext_type::is_product?trial_fecontext_type::nComponents1:1);
                    for ( uint16_type c2 = 0; c2 < ncomp2; ++c2 )
                    {
                        indj.setIndex( boost::make_tuple( j, c2, 0 ) );
                        LIFE_ASSERT( size_type(indj.index()) ==
                                     size_type(c2*trial_fecontext_type::nDof+j) )
                            ( j )( c2 )
                            ( indj.index() )
                            ( c2*trial_fecontext_type::nDof+j ).error( "invalid trial index" );

                        _M_rep[i][c1][j][c2] = M_integrator( *_M_eval_expr00,
                                                             indi,
                                                             indj,
                                                             0, 0 );
                    } // j, c2
                }
            }
        }
    }
}
    template<typename FE1,  typename FE2, typename ElemContType>
        template<typename GeomapContext,typename ExprT,typename IM>
        void
    BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::integrate( mpl::int_<2> )
        {
            //geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_gmc );
            typedef geometric_mapping_context_type gmc_type;
            typedef typename eval00_expr_type::shape shape;
            BOOST_MPL_ASSERT_MSG( (mpl::and_<mpl::equal_to<mpl::int_<shape::M>,mpl::int_<1> >,
                                   mpl::equal_to<mpl::int_<shape::N>,mpl::int_<1> > >::value),


                                  INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                                  (mpl::int_<shape::M>, mpl::int_<shape::N> ) );

            //BOOST_MPL_ASSERT_MSG( eval_expr_type::shape::M == 1, INVALID_TENSOR_EXTENT_M, mpl::int_<eval_expr_type::shape::M> );
            //BOOST_MPL_ASSERT_MSG( eval_expr_type::shape::N == 1, INVALID_TENSOR_EXTENT_N, mpl::int_<eval_expr_type::shape::N> );

            if ( ( _M_form.isPatternDefault() &&
                   ( _M_test_dof->nComponents == _M_trial_dof->nComponents ) ) &&
                 !_M_form.isPatternCoupled() )
            {
                for ( uint16_type i = 0; i < test_fecontext_type::nDof; ++i )
                {
                    int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
                    for ( uint16_type c1 = 0; c1 < ncomp; ++c1 )
                    {
                        indi.setIndex( boost::make_tuple( i, c1, 0 ) );
                        for ( uint16_type j = 0; j < trial_fecontext_type::nDof; ++j )
                        {
                            indj.setIndex( boost::make_tuple( j, c1, 0 ) );
                            _M_rep[i][c1][j][c1] += M_integrator( *_M_eval_expr00, indi, indj, 0, 0 );

                            uint16_type ii = i;
                            uint16_type jj = j + trial_fecontext_type::nDof;

                            _M_rep[ii][c1][jj][c1] += M_integrator( *_M_eval_expr01, indi, indj, 0, 0 );

                            ii = i+ test_fecontext_type::nDof;
                            jj = j ;

                            _M_rep[ii][c1][jj][c1] += M_integrator( *_M_eval_expr10, indi, indj, 0, 0 );

                            ii = i+ test_fecontext_type::nDof;
                            jj = j + trial_fecontext_type::nDof;

                            _M_rep[ii][c1][jj][c1] += M_integrator( *_M_eval_expr11, indi, indj, 0, 0 );
                        }
                    }
                }
            }
            else
            {
                for ( uint16_type i = 0; i < test_fecontext_type::nDof; ++i )
                {
                    int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
                    for ( uint16_type c1 = 0; c1 < ncomp; ++c1 )
                    {
                        indi.setIndex( boost::make_tuple( i, c1, 0 ) );
                        for ( uint16_type j = 0; j < trial_fecontext_type::nDof; ++j )
                        {
                            int ncomp2= (trial_fecontext_type::is_product?trial_fecontext_type::nComponents1:1);
                            for ( uint16_type c2 = 0; c2 < ncomp2; ++c2 )
                            {
                                    indj.setIndex( boost::make_tuple( j, c2, 0 ) );
                                    _M_rep[i][c1][j][c2] += M_integrator( *_M_eval_expr00, indi, indj, 0, 0 );

                                    uint16_type ii = i;
                                    uint16_type jj = j + trial_fecontext_type::nDof;

                                    _M_rep[ii][c1][jj][c2] += M_integrator( *_M_eval_expr01, indi, indj, 0, 0 );

                                    ii = i+ test_fecontext_type::nDof;
                                    jj = j ;

                                    _M_rep[ii][c1][jj][c2] += M_integrator( *_M_eval_expr10, indi, indj, 0, 0 );

                                    ii = i+ test_fecontext_type::nDof;
                                    jj = j + trial_fecontext_type::nDof;

                                    _M_rep[ii][c1][jj][c2] += M_integrator( *_M_eval_expr11, indi, indj, 0, 0 );
                                }
                        }
                    }
                }
            }
        }
    template<typename FE1,  typename FE2, typename ElemContType>
        template<typename GeomapContext,typename ExprT,typename IM>
        void
        BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::assemble( size_type elt_0 )
        {

            size_type ig,jg;
            size_type row_start = _M_lb.front().globalRowStart();
            size_type col_start = _M_lb.front().globalColumnStart();
            int isign, jsign;

#if !defined(NDEBUG)
            Debug( 5050 ) << "[BilinearForm::assemble] global assembly in element " << elt_0 << "\n";
            Debug( 5050 ) << "[BilinearForm::assemble] row start " << row_start << "\n";
            Debug( 5050 ) << "[BilinearForm::assemble] col start " << col_start << "\n";
#endif /* NDEBUG */

            bool do_less = ( ( _M_form.isPatternDefault() &&
                               ( _M_test_dof->nComponents == _M_trial_dof->nComponents ) ) &&
                             !_M_form.isPatternCoupled() );
            if ( do_less )
            {
                for ( uint16_type i = 0 ; i < test_fecontext_type::nDof; ++i )
                {
                    int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
                    for ( uint16_type c1 = 0 ; c1 < ncomp; ++c1 )
                    {
                        boost::tie( ig, isign, boost::tuples::ignore) = _M_test_dof->localToGlobal( elt_0, i, c1 );
                        ig += row_start;
                        for ( uint16_type j = 0 ; j < trial_fecontext_type::nDof; j++ )
                        {
                            boost::tie( jg, jsign, boost::tuples::ignore ) = _M_trial_dof->localToGlobal( elt_0, j, c1 );
                            jg += col_start;
                            _M_form.add( ig, jg, value_type(isign*jsign)*_M_rep[i][c1][j][c1] );
                        }
                    }
                } // element loop
            }
            else // couple the components
            {
                for ( uint16_type i = 0 ; i < test_fecontext_type::nDof; ++i )
                {
                    int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
                    for ( uint16_type c1 = 0 ; c1 < ncomp; ++c1 )
                    {
                        boost::tie( ig, isign, boost::tuples::ignore) = _M_test_dof->localToGlobal( elt_0, i, c1 );
                        ig += row_start;
                        for ( uint16_type j = 0 ; j < trial_fecontext_type::nDof; j++ )
                        {
                            int ncomp2= (trial_fecontext_type::is_product?trial_fecontext_type::nComponents1:1);
                            for ( uint16_type c2 = 0 ; c2 < ncomp2; ++c2 )
                            {
                                boost::tie( jg, jsign, boost::tuples::ignore ) = _M_trial_dof->localToGlobal( elt_0, j, c2 );
                                jg += col_start;
                                _M_form.add( ig, jg, value_type(isign*jsign)*_M_rep[i][c1][j][c2] );
                            }
                        }
                    }
                } // element loop
            }
        }

    template<typename FE1,  typename FE2, typename ElemContType>
        template<typename GeomapContext,typename ExprT,typename IM>
        void
        BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::assemble( size_type elt_0, size_type elt_1  )
        {
            test_index_type indi;
            trial_index_type indj;
            size_type ig0,ig1,jg0,jg1;
            size_type row_start = _M_lb.front().globalRowStart();
            size_type col_start = _M_lb.front().globalColumnStart();
            int isign0, isign1, jsign0, jsign1;
            const uint16_type test_ndof = test_fecontext_type::nDof;
            const uint16_type trial_ndof = trial_fecontext_type::nDof;
            //std::cout << "rep_" << elt_0 << "=" << _M_rep << "\n";

            if ( ( _M_form.isPatternDefault() &&
                   ( _M_test_dof->nComponents == _M_trial_dof->nComponents ) ) &&
                 !_M_form.isPatternCoupled() )
            {
                for ( uint16_type i = 0 ; i < test_fecontext_type::nDof; ++i )
                {
                    int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
                    for ( uint16_type c1 = 0 ; c1 < ncomp; ++c1 )
                    {
                        boost::tie( ig0, isign0, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_0, i, c1 );
                        boost::tie( ig1, isign1, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_1, i, c1 );
                        ig0 += row_start;
                        ig1 += row_start;
                        for ( uint16_type j = 0 ; j < trial_fecontext_type::nDof; j++ )
                        {
                            boost::tie( jg0, jsign0, boost::tuples::ignore ) = _M_trial_dof->localToGlobal( elt_0, j, c1 );
                            boost::tie( jg1, jsign1, boost::tuples::ignore ) = _M_trial_dof->localToGlobal( elt_1, j, c1 );
                            jg0 += col_start;
                            jg1 += col_start;

                            _M_form.add( ig0, jg0, value_type(isign0*jsign0)*_M_rep[i][c1][j][c1] );
                            _M_form.add( ig0, jg1, value_type(isign0*jsign1)*_M_rep[i][c1][j+trial_ndof][c1] );
                            _M_form.add( ig1, jg0, value_type(isign0*jsign0)*_M_rep[i+test_ndof][c1][j][c1] );
                            _M_form.add( ig1, jg1, value_type(isign0*jsign1)*_M_rep[i+test_ndof][c1][j+trial_ndof][c1] );
                        }
                    }
                }
            }
            else
            {
                for ( uint16_type i = 0 ; i < test_fecontext_type::nDof; ++i )
                {
                    int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
                    for ( uint16_type c1 = 0 ; c1 < test_fecontext_type::nComponents1; ++c1 )
                    {
                        boost::tie( ig0, isign0, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_0, i, c1 );
                        boost::tie( ig1, isign1, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_1, i, c1 );
                        ig0 += row_start;
                        ig1 += row_start;
                        for ( uint16_type j = 0 ; j < trial_fecontext_type::nDof; j++ )
                        {
                            int ncomp2= (trial_fecontext_type::is_product?trial_fecontext_type::nComponents1:1);
                            for ( uint16_type c2 = 0 ; c2 < ncomp2; ++c2 )
                                {
                                    boost::tie( jg0, jsign0, boost::tuples::ignore ) = _M_trial_dof->localToGlobal( elt_0, j, c2 );
                                    boost::tie( jg1, jsign1, boost::tuples::ignore ) = _M_trial_dof->localToGlobal( elt_1, j, c2 );
                                    jg0 += col_start;
                                    jg1 += col_start;

                                    _M_form.add( ig0, jg0, value_type(isign0*jsign0)*_M_rep[i][c1][j][c2] );
                                    _M_form.add( ig0, jg1, value_type(isign0*jsign1)*_M_rep[i][c1][j+trial_ndof][c2] );
                                    _M_form.add( ig1, jg0, value_type(isign0*jsign0)*_M_rep[i+test_ndof][c1][j][c2] );
                                    _M_form.add( ig1, jg1, value_type(isign0*jsign1)*_M_rep[i+test_ndof][c1][j+trial_ndof][c2] );
                                }
                        }
                    }
                }
            }

        }
}
}
}
#endif /* __BilinearFormContext_H */
