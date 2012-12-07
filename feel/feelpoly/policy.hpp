/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-12-03

  Copyright (C) 2005,2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file policy.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-12-03
 */
#ifndef __policy_H
#define __policy_H 1



#include <boost/mpl/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/mpl/if.hpp>
#include <boost/mpl/find.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
namespace Feel
{
namespace ublas = boost::numeric::ublas;

namespace fem
{
enum transformation_type { LINEAR, BILINEAR,  NONLINEAR };
}

/**
 * Policy for \c Scalar polynomials or polynomial set of dimension
 * \p Dim
 * \note \c Scalar can be seen as rank 0 tensor polynomials
 */
template<uint16_type Dim>
struct Scalar
{
    static const uint16_type rank = 0;
    static const uint16_type nDim = Dim;

    static const bool is_scalar = true;
    static const bool is_vectorial = false;
    static const bool is_tensor2 = false;
    static const bool is_tensor3 = false;

    static const uint16_type nComponents = 1;
    static const uint16_type nComponents1 = 1;
    static const uint16_type nComponents2 = 1;
    static const uint16_type nComponents3 = 1;
    static const uint16_type nComponentsLast = 1;

    template<typename T>
    static
    inline ublas::matrix<T> const&
    toMatrix( ublas::matrix<T> const&  __c )
    {
        return __c;
    }
    template<typename T>
    static
    inline ublas::matrix<T> const&
    toType( ublas::matrix<T> const&  __c )
    {
        return __c;
    }
};


/**
 * Policy for \c Vectorial polynomials or polynomial sets of dimension
 * \p Dim
 * \note \c Vectorial can be seen as rank 1 Tensor polynomials
 */
template<uint16_type Dim>
struct Vectorial
{
    static const uint16_type rank = 1;
    static const uint16_type nDim = Dim;

    static const bool is_scalar = false;
    static const bool is_vectorial = true;
    static const bool is_tensor2 = false;
    static const bool is_tensor3 = false;

    static const uint16_type nComponents = nDim;
    static const uint16_type nComponents1 = nDim;
    static const uint16_type nComponents2 = 1;
    static const uint16_type nComponents3 = 1;
    static const uint16_type nComponentsLast = nComponents1;
    template<typename T>
    static
    ublas::matrix<T>
    toMatrix( ublas::matrix<T> const&  __c )
    {
        typedef T value_type;
        // reshape the coefficients in the vectorial case
        const size_type nRows = __c.size1()/nComponents;
        const size_type nCols = __c.size2();
        ublas::matrix<T> __c_reshaped( ublas::zero_matrix<value_type>( nRows, nComponents*nCols ) );

        for ( int c = 0; c < nComponents; ++c )
        {
            ublas::project( __c_reshaped,
                            ublas::range( 0, nRows ),
                            ublas::range( c*nCols, ( c+1 )*nCols ) ) = ublas::project( __c,
                                    ublas::slice( c, nComponents, nRows ),
                                    ublas::slice( 0, 1, nCols ) );
        }

        return __c_reshaped;
    }

    template<typename AE>
    static
    ublas::matrix<typename ublas::matrix_expression<AE>::value_type>
    toType( ublas::matrix_expression<AE> const&  __c )
    {
        typedef typename ublas::matrix_expression<AE>::value_type value_type;

        // reshape the coefficients in the vectorial case
        const size_type nRows = __c().size1()*nComponents;
        const size_type nCols = __c().size2()/nComponents;
        ublas::matrix<value_type> __c_reshaped( ublas::zero_matrix<value_type>( nRows, nCols ) );

        for ( int c = 0; c < nComponents; ++c )
        {
            ublas::project( __c_reshaped,
                            ublas::slice( c, nComponents, nRows/nComponents ),
                            ublas::slice( 0, 1, nCols ) ) = ublas::project( __c(),
                                    ublas::range( 0, nRows/nComponents ),
                                    ublas::range( c*nCols, ( c+1 )*nCols ) );
        }

        return __c_reshaped;
    }

    template<typename T>
    static
    ublas::matrix<T>
    toType( ublas::matrix<T> const&  __c )
    {
        typedef T value_type;

        // reshape the coefficients in the vectorial case
        const size_type nRows = __c.size1()*nComponents;
        const size_type nCols = __c.size2()/nComponents;
        ublas::matrix<value_type> __c_reshaped( ublas::zero_matrix<value_type>( nRows, nCols ) );

        for ( int c = 0; c < nComponents; ++c )
        {
            ublas::project( __c_reshaped,
                            ublas::slice( c, nComponents, nRows/nComponents ),
                            ublas::slice( 0, 1, nCols ) ) = ublas::project( __c,
                                    ublas::range( 0, nRows/nComponents ),
                                    ublas::range( c*nCols, ( c+1 )*nCols ) );
        }

        return __c_reshaped;
    }
};

/**
 * Policy for \c Scalar polynomials or polynomial set of dimension
 * \p Dim
 */
namespace detail
{
template<uint16_type N, uint16_type M>
struct Field
{
    static const uint16_type rank = ( M > 1 );
    static const uint16_type nDim = N;
    static const uint16_type nVariables = N;

    static const bool is_scalar = ( M==1 );
    static const bool is_vectorial = ( N==M );
    static const bool is_tensor2 = false;
    static const bool is_tensor3 = false;

    static const uint16_type nComponents = M;
    static const uint16_type nComponents1 = M;
    static const uint16_type nComponents2 = 1;
    static const uint16_type nComponents3 = 1;
    static const uint16_type nComponentsLast = 1;

    template<typename T>
    static
    inline ublas::matrix<T> const&
    toMatrix( ublas::matrix<T> const&  __c )
    {
        return __c;
    }
    template<typename T>
    static
    inline ublas::matrix<T> const&
    toType( ublas::matrix<T> const&  __c )
    {
        return __c;
    }
};
}

template<uint16_type M>
struct Field
{
    template <uint16_type Nvar>
    struct apply
    {
        typedef detail::Field<Nvar,M> type;
    };
};

/**
 * Policy for rank 2 tensor polynomials or polynomial sets of
 * dimension \p Dim
 *
 */
template<uint16_type Dim>
struct Tensor2
{
    static const uint16_type rank = 2;
    static const uint16_type nDim = Dim;

    static const bool is_scalar = false;
    static const bool is_vectorial = false;
    static const bool is_tensor2 = true;
    static const bool is_tensor3 = false;

    static const uint16_type nComponents = nDim*nDim;
    static const uint16_type nComponents1 = nDim;
    static const uint16_type nComponents2 = nDim;
    static const uint16_type nComponents3 = 1;
    static const uint16_type nComponentsLast = nComponents2;

    template<typename T>
    static
    ublas::matrix<T>
    toMatrix( ublas::matrix<T> const&  __c )
    {
        typedef T value_type;
        // reshape the coefficients in the vectorial case
        const size_type nRows = __c.size1();
        const size_type nRows1= __c.size2()*nComponents;
        const size_type nCols = __c.size2();
        ublas::matrix<T> __c_reshaped( nRows/nComponents, nCols*nComponents );

        //__c_reshaped = ublas::scalar_matrix<value_type>( nRows/nComponents, nCols*nComponents, -1 );
        for ( int c1 = 0; c1 < nComponents; ++c1 )
        {
            uint16_type i1 = nRows1*c1;

            for ( int c2 = 0; c2 < nComponents; ++c2 )
            {
                ublas::project( __c_reshaped,
                                ublas::range( i1/nComponents, ( i1+nRows1 )/nComponents ),
                                ublas::range( c2*nCols, ( c2+1 )*nCols ) ) = ublas::project( __c,
                                        ublas::slice( i1+c2, nComponents, nRows1/nComponents ),
                                        ublas::slice( 0, 1, nCols ) );
            }
        }

        return __c_reshaped;
    }

    template<typename AE>
    static
    ublas::matrix<typename ublas::matrix_expression<AE>::value_type>
    toType( ublas::matrix_expression<AE> const&  __c )
    {

    }

    template<typename T>
    static
    ublas::matrix<T>
    toType( ublas::matrix<T> const&  __c )
    {
        typedef T value_type;
        // reshape the coefficients in the vectorial case
        const size_type nRows = __c.size1()*nComponents;
        const size_type nRows1= __c.size2();
        const size_type nCols = __c.size2()/nComponents;
        ublas::matrix<T> __c_reshaped( nRows, nCols );

        //__c_reshaped = ublas::scalar_matrix<value_type>( nRows, nCols, -1 );
        for ( int c1 = 0; c1 < nComponents; ++c1 )
        {
            uint16_type i1 = nRows1*c1;

            for ( int c2 = 0; c2 < nComponents; ++c2 )
            {
                ublas::project( __c_reshaped,
                                ublas::slice( i1+c2, nComponents, nRows1/nComponents ),
                                ublas::slice( 0, 1, nCols ) ) = ublas::project( __c,
                                        ublas::range( i1/nComponents, ( i1+nRows1 )/nComponents ),
                                        ublas::range( c2*nCols, ( c2+1 )*nCols ) );
            }
        }

        return __c_reshaped;
    }
};

/**
 * Policy for rank 2 tensor polynomials or polynomial sets of
 * dimension \p Dim
 *
 */
template<uint16_type Dim>
struct Tensor3
{
    static const uint16_type rank = 3;
    static const uint16_type nDim = Dim;

    static const bool is_scalar = false;
    static const bool is_vectorial = false;
    static const bool is_tensor2 = false;
    static const bool is_tensor3 = true;

    static const uint16_type nComponents = nDim*nDim*nDim;
    static const uint16_type nComponents1 = nDim;
    static const uint16_type nComponents2 = nDim;
    static const uint16_type nComponents3 = nDim;
    static const uint16_type nComponentsLast = nComponents3;
};

template<int Dim>
struct ListReturnTypes
{
    typedef mpl::vector<Scalar<Dim>, Vectorial<Dim>, Tensor2<Dim>, Tensor3<Dim> > return_types;
};
template<typename T1, typename T2>
struct ReturnSelect
{
    typedef typename mpl::if_<boost::is_same<T1, T2>,
            mpl::identity<T1>,
            typename mpl::if_<mpl::greater<mpl::int_<T1::rank>, mpl::int_<T2::rank> >,
            mpl::identity<T1>,
            mpl::identity<T2> >::type>::type::type type;
};

enum EnumIndex
{
    GLOBAL_COMPONENT,
    COMPONENT_IN_COMPONENT,
    GLOBAL_FUNCTION_INDEX,
    PER_COMPONENT_FUNCTION_INDEX,
    COMPONENT_IN_COMPONENT_FUNCTION_INDEX,
    FUNCTION_INDEX
};

const mpl::int_<GLOBAL_COMPONENT> INDEX_GLOBAL_COMPONENT = mpl::int_<GLOBAL_COMPONENT>();
const mpl::int_<COMPONENT_IN_COMPONENT> INDEX_COMPONENT_IN_COMPONENT = mpl::int_<COMPONENT_IN_COMPONENT>() ;
const mpl::int_<GLOBAL_FUNCTION_INDEX> INDEX_GLOBAL_FUNCTION_INDEX = mpl::int_<GLOBAL_FUNCTION_INDEX>();
const mpl::int_<PER_COMPONENT_FUNCTION_INDEX> INDEX_PER_COMPONENT_FUNCTION_INDEX = mpl::int_<PER_COMPONENT_FUNCTION_INDEX>();
const mpl::int_<COMPONENT_IN_COMPONENT_FUNCTION_INDEX> INDEX_COMPONENT_IN_COMPONENT_FUNCTION_INDEX = mpl::int_<COMPONENT_IN_COMPONENT_FUNCTION_INDEX>();
const mpl::int_<FUNCTION_INDEX> INDEX_FUNCTION_INDEX = mpl::int_<FUNCTION_INDEX>();

/**
 * Get the component type out the available types
 * \code
 * typedef typename Component<Vectorial<3> >::type component_type;
 * // component_type should be of type \c Scalar<3>
 * \endcode
 */
template<typename T>
struct Component
{
    static const uint16_type nDim = T::nDim;
    typedef mpl::vector<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim> > types;
#if 0
    typedef typename mpl::find<types, T>::type iter;
    typedef typename mpl::if_<boost::is_same<T, Scalar<nDim> >,
            mpl::identity<Scalar<nDim> >,
            mpl::identity<typename mpl::deref<typename mpl::prior<iter>::type>::type> >::type::type type;
#else
    typedef typename mpl::if_<boost::is_same<T, Scalar<nDim> >,
            mpl::identity<Scalar<nDim> >,
            typename mpl::if_<boost::is_same<T, Vectorial<nDim> >,
            mpl::identity<Vectorial<nDim> >,
            typename mpl::if_<boost::is_same<T, Tensor2<nDim> >,
            mpl::identity<Tensor2<nDim> > >::type>::type>::type::type type;
#endif
};

/**
 * Get the next rank
 * \code
 * typedef typename RankUp<Scalar<3> >::type rank_type;
 * // rank_type should be of type \c Vectorial<3>
 * \endcode
 */
template<typename T>
struct RankUp
{
    static const uint16_type nDim = T::nDim;
    typedef mpl::vector<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim> > types;
    typedef typename mpl::find<types, T>::type iter;
    typedef typename mpl::deref<typename mpl::next<iter>::type>::type type;
};

template<typename T>
struct RankUp2
{
    static const uint16_type nDim = T::nDim;
    typedef mpl::vector<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim> > types;
    typedef typename mpl::find<types, T>::type iter;
    typedef typename mpl::deref<typename mpl::next<typename mpl::next<iter>::type>::type>::type type;
};

template<typename T>
struct RankSame
{
    static const uint16_type nDim = T::nDim;
    typedef T type;
};

template<typename T>
struct RankDown
{
    static const uint16_type nDim = T::nDim;
    typedef mpl::vector<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim> > types;
    typedef typename mpl::find<types, T>::type iter;
    typedef typename mpl::deref<typename mpl::prior<iter>::type>::type _type;
    typedef typename mpl::if_<boost::is_same<_type,mpl::void_>,
            mpl::identity<Scalar<nDim> >,
            mpl::identity<_type> >::type::type type;
};
/**
 * Policy for orthogonal polynomial set which can be \p normalized or
 * not.  This can be used as a policy for orthogonal polynomial sets
 * in order to select either the normalized or un-normalized version
 */
template<bool normalized>
struct Normalized
{
    static const bool is_normalized = normalized;
};

/**
 * Storage Policy using ublas for numerical type \p T
 *
 *
 */
template<typename T>
struct StorageUBlas
{
    typedef T value_type;
    typedef ublas::vector<value_type> vector_type;
    typedef ublas::matrix<value_type> matrix_type;
    typedef ublas::vector<matrix_type> vector_matrix_type;
    typedef ublas::vector<vector_matrix_type> vector_vector_matrix_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;
    typedef typename matrix_node<value_type>::type points_type;
    typedef typename node<value_type>::type node_type;

};

} // Feel

#endif /* __policy_H */
