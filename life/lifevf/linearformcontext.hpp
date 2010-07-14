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
   \file linearformcontext.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-04-27
 */
#ifndef __LinearFormContext_H
#define __LinearFormContext_H 1


namespace Life
{
namespace vf
{
namespace detail
{
//
// Context class for linear forms
//
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT, typename IM>
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                                           map_geometric_mapping_context_type const& _gmc,
                                                                                           ExprT const& expr,
                                                                                           IM const& im )
    :
    //super(),
    _M_form( __form ),
    _M_test_dof( __form.functionSpace()->dof().get() ),
    _M_lb( __form.blockList() ),
    _M_gmc( _gmc ),
    _M_gmc_left( fusion::at_key<gmc<0> >( _gmc ) ),
    _M_left_map( fusion::make_map<gmc<0> >( _M_gmc_left ) ),
    _M_test_fec( fusion::transform( _M_gmc, detail::FEContextInit<0,form_context_type>(__form.functionSpace()->fe(), _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_rep(boost::extents[fusion::size(_gmc)*test_fecontext_type::nDof]
           [test_fecontext_type::nComponents1]),
    _M_eval0_expr( new eval0_expr_type( expr, _gmc, _M_test_fec0 ) ),
    _M_eval1_expr(),
    M_integrator( im )
{
    _M_eval0_expr->init( im );
}

template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT, typename IM>
template<typename IM2>
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                                           map_geometric_mapping_context_type const& _gmc,
                                                                                           ExprT const& expr,
                                                                                           IM const& im,
                                                                                           IM2 const& im2 )
    :
    //super(),
    _M_form( __form ),
    _M_test_dof( __form.functionSpace()->dof().get() ),
    _M_lb( __form.blockList() ),
    _M_gmc( _gmc ),
    _M_gmc_left( fusion::at_key<gmc<0> >( _gmc ) ),
    _M_left_map( fusion::make_map<gmc<0> >( _M_gmc_left ) ),
    _M_test_fec( fusion::transform( _M_gmc, detail::FEContextInit<0,form_context_type>(__form.functionSpace()->fe(), _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_rep(boost::extents[fusion::size(_gmc)*test_fecontext_type::nDof]
           [test_fecontext_type::nComponents1]),
    _M_eval0_expr( new eval0_expr_type( expr, _gmc, _M_test_fec0 ) ),
    _M_eval1_expr(),
    M_integrator( im )
{
    _M_eval0_expr->init( im2 );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT, typename IM>
template<typename IM2>
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                                           map_geometric_mapping_context_type const& _gmc,
                                                                                           ExprT const& expr,
                                                                                           IM const& im,
                                                                                           IM2 const& im2,
                                                                                           mpl::int_<2> )
    :
    //super(),
    _M_form( __form ),
    _M_test_dof( __form.functionSpace()->dof().get() ),
    _M_lb( __form.blockList() ),
    _M_gmc( _gmc ),
    _M_gmc_left( fusion::at_key<gmc<0> >( _gmc ) ),
    _M_gmc_right( fusion::at_key<gmc1 >( _gmc ) ),
    _M_left_map( fusion::make_map<gmc<0> >( _M_gmc_left ) ),
    _M_right_map( fusion::make_map<gmc<0> >( _M_gmc_right ) ),
    _M_test_fec( fusion::transform( _M_gmc, detail::FEContextInit<0,form_context_type>(__form.functionSpace()->fe(), _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_test_fec1( fusion::make_pair<gmc1 >( fusion::at_key<gmc1 >( _M_test_fec ) ) ),
    _M_rep(boost::extents[fusion::size(_gmc)*test_fecontext_type::nDof]
           [test_fecontext_type::nComponents1]),
    _M_eval0_expr( new eval0_expr_type( expr, _gmc, _M_test_fec0 ) ),
    _M_eval1_expr( new eval1_expr_type( expr, _gmc, _M_test_fec1 ) ),
    M_integrator( im )

{
    _M_eval0_expr->init( im2 );
    _M_eval1_expr->init( im2 );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc )
{
    _M_gmc = _gmc;
    _M_gmc_left = fusion::at_key<gmc<0> >( _gmc );
    _M_left_map = fusion::make_map<gmc<0> >( _M_gmc_left );
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmc, _M_form ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_eval0_expr->update( _gmc, _M_test_fec0 );
    std::fill( _M_rep.data(), _M_rep.data()+_M_rep.num_elements(), value_type(0) );
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmc ) );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, mpl::int_<2> )
{
    //    typedef mpl::int_<fusion::result_of::template size<map_geometric_mapping_context_type>::type::value> map_size;
    /*BOOST_MPL_ASSERT_MSG( (mpl::equal_to<mpl::int_<map_size::value>,mpl::int_<2> >::value),
                          INVALID_GEOMAP,
                          (map_size,map_geometric_mapping_context_type ));*/
    //_M_gmc = _gmc;
#if 0
    _M_gmc_left = fusion::at_key<gmc<0> >( _gmc );
    _M_gmc_right =  fusion::at_key<gmc1 >( _gmc );
    _M_left_map = fusion::make_map<gmc<0> >( _M_gmc_left );
    _M_right_map = fusion::make_map<gmc<0> >( _M_gmc_right );
#endif
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmc, _M_form ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_test_fec1 = fusion::make_map<gmc1 >( fusion::at_key<gmc1 >( _M_test_fec ) );
    _M_eval0_expr->update( _gmc, _M_test_fec0 );
    _M_eval1_expr->update( _gmc, _M_test_fec1 );
    std::fill( _M_rep.data(), _M_rep.data()+_M_rep.num_elements(), value_type(0) );
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmc ) );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, IM const& im )
{
    M_integrator = im;
    this->update( _gmc );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, IM const& im, mpl::int_<2> )
{
    M_integrator = im;
    this->update( _gmc, mpl::int_<2>() );
}


template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::integrate( mpl::int_<1> )
{
    //geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_gmc );
    typedef geometric_mapping_context_type gmc_type;

    typedef typename eval0_expr_type::shape shape;
    BOOST_MPL_ASSERT_MSG( (shape::M == 1 && shape::N == 1),
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N> ) );

    test_index_type indi;
    int i, c1;
    value_type res;
    //#pragma omp parallel
    //{
        //#pragma omp single
        //Debug() << "[linearform::integrate] num threads: " << OMP_GET_NUM_THREADS << "\n";

    //#pragma omp for private(i,c1,indi, res)
        for ( i = 0; i < test_fecontext_type::nDof; ++i )
        {
            int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
            for( c1 = 0;c1 < ncomp; ++ c1 )
                {
                    indi.setIndex( boost::make_tuple( i, c1, 0 ) );
                    res = M_integrator( *_M_eval0_expr, indi, 0, 0 );
                    _M_rep[i][c1] = res;
                }
        }
        //}
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::integrate( mpl::int_<2> )
{
    //geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_gmc );
    typedef geometric_mapping_context_type gmc_type;

    typedef mpl::int_<fusion::result_of::template size<map_geometric_mapping_context_type>::type::value> map_size;
    BOOST_MPL_ASSERT_MSG( map_size::value == 2, INVALID_GEOMAP, (map_size,map_geometric_mapping_context_type ));

    typedef typename eval0_expr_type::shape shape;
    BOOST_MPL_ASSERT_MSG( (shape::M == 1 && shape::N == 1),
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N> ) );
    test_index_type indi;

    for ( uint16_type i = 0; i < test_fecontext_type::nDof; ++i )
    {
        int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
        for( uint16_type c1 = 0;c1 < ncomp; ++ c1 )
            {
                indi.setIndex( boost::make_tuple( i, c1, 0 ) );

                LIFE_ASSERT( M_integrator.isFaceIm() )
                    ( M_integrator ).error( "invalid face integrator" );
                _M_rep[i][c1] = M_integrator( *_M_eval0_expr, indi, 0, 0 );
                uint16_type ii = i + test_fecontext_type::nDof;
                _M_rep[ii][c1] = M_integrator( *_M_eval1_expr, indi, 0, 0 );
            }
    }
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::assemble( size_type elt_0 )
{
    test_index_type indi;
    size_type ig;
    int isign;
    size_type row_start = _M_lb.front().globalRowStart();
    for ( uint16_type k = 0 ; k < test_fecontext_type::nDof; k++ )
    {
        int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
        for( uint16_type c1 = 0;c1 < ncomp; ++ c1 )
            {
                //uint16_type c = test_fecontext_type::nComponents2*c1+c2;
                boost::tie( ig, isign, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_0, k, c1 );
                ig += row_start;
                _M_form.add( ig, value_type(isign)*_M_rep[k][c1] );
            }
    }
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::assemble( size_type elt_0, size_type elt_1 )
{
    test_index_type indi;
    size_type ig0,ig1;
    int isign0,isign1;
    size_type row_start = _M_lb.front().globalRowStart();
    for ( uint16_type k = 0 ; k < test_fecontext_type::nDof; k++ )
    {
        int ncomp= (test_fecontext_type::is_product?test_fecontext_type::nComponents1:1);
        for( uint16_type c1 = 0;c1 < ncomp; ++ c1 )
            //for( uint16_type c2 = 0;c2 < test_fecontext_type::nComponents2; ++ c2 )
            {
                //uint16_type c = test_fecontext_type::nComponents2*c1+c2;
                boost::tie( ig0, isign0, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_0, k, c1 );
                ig0 += row_start;
                _M_form.add( ig0, value_type(isign0)*_M_rep[k][c1] );

                boost::tie( ig1, isign1, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_1, k, c1 );
                ig1 += row_start;
                size_type l = test_fecontext_type::nDof + k;
                _M_form.add( ig1, value_type(isign1)*_M_rep[l][c1] );
            }
    }
}

}
}
}
 #endif /* __LinearFormContext_H */
