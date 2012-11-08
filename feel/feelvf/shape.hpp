/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-11-13

  Copyright (C) 2006,2007 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2012 Universite de Strasbourg

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
   \file shape.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-11-13
 */
#ifndef __Shape_H
#define __Shape_H 1

#include <boost/numeric/ublas/storage.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class Shape
 * \brief Shape class
 *
 * @author Christophe Prud'homme
 * @see
 */
template<uint16_type Dim, template<uint16_type D> class Type, bool is_transposed, bool diag = false>
class Shape
{
};

/**
 * \class Shape
 * \brief Shape class for scalar functions
 *
 * @author Christophe Prud'homme
 * @see
 */
template<uint16_type Dim, bool transpose, bool diag>
class Shape<Dim,Scalar,transpose,diag>
{
public:

    /** @name Constants
     */
    //@{

    static const uint16_type nDim = Dim;
    static const uint16_type M = 1;
    static const uint16_type N = 1;
    static const uint16_type O = 1;
    static const bool is_transposed = transpose;
    static const bool is_diagonalized = diag;
    static const uint16_type rank = Scalar<nDim>::rank;

    static const bool is_scalar = Scalar<nDim>::is_scalar;
    static const bool is_vectorial = Scalar<nDim>::is_vectorial;
    static const bool is_tensor2 = Scalar<nDim>::is_tensor2;
    static const bool is_tensor3 = Scalar<nDim>::is_tensor3;
    //@}

    /** @name Typedefs
     */
    //@{

    typedef Shape<Dim,Scalar,transpose, is_diagonalized> shape_type;

    //@}
};

/**
 * \class Shape
 * \brief Shape class for vector functions
 *
 * @author Christophe Prud'homme
 * @see
 */
template<uint16_type Dim, bool transpose, bool diag>
class Shape<Dim,Vectorial,transpose, diag>
{
public:

    /** @name Constants
     */
    //@{

    static const uint16_type nDim = Dim;
    static const uint16_type M = mpl::if_<mpl::equal_to<mpl::bool_<transpose>, mpl::bool_<true> >,
                             mpl::int_<1>, mpl::int_<nDim> >::type::value;
    static const uint16_type N = mpl::if_<mpl::equal_to<mpl::bool_<transpose>, mpl::bool_<true> >, mpl::int_<nDim>, mpl::int_<1> >::type::value;
    static const uint16_type O = 1;
    static const bool is_transposed = transpose;
    static const bool is_diagonalized = diag;
    static const uint16_type rank = Vectorial<nDim>::rank;

    static const bool is_scalar = Vectorial<nDim>::is_scalar;
    static const bool is_vectorial = Vectorial<nDim>::is_vectorial;
    static const bool is_tensor2 = Vectorial<nDim>::is_tensor2;
    static const bool is_tensor3 = Vectorial<nDim>::is_tensor3;
    //@}

    /** @name Typedefs
     */
    //@{

    typedef Shape<Dim,Vectorial,transpose, is_diagonalized> shape_type;
    //@}
};

/**
 * \class Shape
 * \brief Shape class for matricial functions
 *
 * @author Christophe Prud'homme
 * @see
 */
template<uint16_type Dim, bool transpose, bool diag>
class Shape<Dim,Tensor2,transpose,diag>
{
public:

    /** @name Constants
     */
    //@{

    static const uint16_type nDim = Dim;
    static const uint16_type M = nDim;
    static const uint16_type N = nDim;
    static const uint16_type O = 1;
    static const bool is_transposed = transpose;
    static const bool is_diagonalized = diag;
    static const uint16_type rank = Tensor2<nDim>::rank;

    static const bool is_scalar = Tensor2<nDim>::is_scalar;
    static const bool is_vectorial = Tensor2<nDim>::is_vectorial;
    static const bool is_tensor2 = Tensor2<nDim>::is_tensor2;
    static const bool is_tensor3 = Tensor2<nDim>::is_tensor3;
    //@}

    /** @name Typedefs
     */
    //@{

    typedef Shape<Dim,Tensor2,transpose,diag> shape_type;
    //@}
};

/**
 * \class Shape
 * \brief Shape class for matricial functions
 *
 * @author Christophe Prud'homme
 * @see
 */
template<uint16_type Dim, bool transpose, bool diag>
class Shape<Dim,Tensor3,transpose,diag>
{
public:

    /** @name Constants
     */
    //@{

    static const uint16_type nDim = Dim;
    static const uint16_type M = nDim;
    static const uint16_type N = nDim;
    static const uint16_type O = nDim;
    static const bool is_transposed = transpose;
    static const bool is_diagonalized = diag;
    static const uint16_type rank = Tensor3<nDim>::rank;

    static const bool is_scalar = Tensor3<nDim>::is_scalar;
    static const bool is_vectorial = Tensor3<nDim>::is_vectorial;
    static const bool is_tensor2 = Tensor3<nDim>::is_tensor2;
    static const bool is_tensor3 = Tensor3<nDim>::is_tensor3;
    //@}

    /** @name Typedefs
     */
    //@{

    typedef Shape<Dim,Tensor3,transpose,diag> shape_type;
    //@}
};

//template<uint16_type Dim, template<uint16_type D> class Type, bool transpose>
template<typename TheShape>
class Transpose
{
public:
    /** @name Typedefs
     */
    //@{
    //typedef Shape<Dim,Type,transpose> original_shape_type;
    typedef TheShape original_shape_type;

    typedef typename mpl::if_<mpl::bool_<original_shape_type::is_scalar>,
            mpl::identity<Shape<original_shape_type::nDim,Scalar,!original_shape_type::is_transposed,original_shape_type::is_diagonalized> >,
            typename mpl::if_<mpl::bool_<original_shape_type::is_vectorial>,
            mpl::identity<Shape<original_shape_type::nDim,Vectorial,!original_shape_type::is_transposed,original_shape_type::is_diagonalized> >,
            typename mpl::if_<mpl::bool_<original_shape_type::is_tensor2>,
            mpl::identity<Shape<original_shape_type::nDim,Tensor2,!original_shape_type::is_transposed,original_shape_type::is_diagonalized> >,
            mpl::identity<Shape<original_shape_type::nDim,Tensor3,!original_shape_type::is_transposed,original_shape_type::is_diagonalized> >
            >::type >::type >::type::type type;

    //@}
};

//template<uint16_type Dim, template<uint16_type D> class Type, bool transpose>
template<typename TheShape>
class Diag
{
public:
    /** @name Typedefs
     */
    //@{
    //typedef Shape<Dim,Type,transpose> original_shape_type;
    typedef TheShape original_shape_type;

    typedef typename mpl::if_<mpl::bool_<original_shape_type::is_scalar>,
            mpl::identity<Shape<original_shape_type::nDim,Scalar,original_shape_type::is_transposed,original_shape_type::is_diagonalized> >,
            typename mpl::if_<mpl::bool_<original_shape_type::is_vectorial>,
            mpl::identity<Shape<original_shape_type::nDim,Tensor2,original_shape_type::is_transposed,original_shape_type::is_diagonalized> >,
            mpl::identity<Shape<original_shape_type::nDim,Vectorial,original_shape_type::is_transposed,original_shape_type::is_diagonalized> >
            >::type >::type::type type;

    //@}
};


template<typename Left, typename Right>
struct shape_op_samerank
{
    static const bool is_rank_ok = ( ( Left::M == Right::M ) && ( Left::N == Right::N ) );
    BOOST_MPL_ASSERT_MSG( mpl::bool_<is_rank_ok>::value,
                          INVALID_TENSOR_OPERATION_SHOULD_BE_OF_SAME_RANK,
                          ( mpl::int_<Left::M>, mpl::int_<Left::N>,
                            mpl::int_<Right::M>, mpl::int_<Right::N>,
                            Left, Right ) );


    static const uint16_type nDim = Left::nDim;
    static const bool is_diagonalized = Left::is_diagonalized && Right::is_diagonalized;
    static const bool is_transposed = Left::is_transposed && Right::is_transposed;

    typedef typename mpl::if_<mpl::bool_<Left::is_scalar>,
            mpl::identity<Shape<nDim,Scalar,is_transposed,is_diagonalized> >,
            typename mpl::if_<mpl::bool_<Left::is_vectorial>,
            mpl::identity<Shape<nDim,Vectorial,is_transposed,is_diagonalized> >,
            typename mpl::if_<mpl::bool_<Left::is_tensor2>,
            mpl::identity<Shape<nDim,Tensor2,is_transposed,is_diagonalized> >,
            mpl::identity<Shape<nDim,Tensor3,is_transposed,is_diagonalized> >
            >::type >::type >::type::type type;

    static const int op = 0;

    template<bool left_is_zero, bool right_is_zero>
    struct is_zero
    {
        static const bool value = ( left_is_zero && right_is_zero );
        static const bool update_and_eval_left = !left_is_zero;
        static const bool update_and_eval_right = !right_is_zero;
    };
};

template<typename Left, typename Right>
struct shape_op_id
{
    static const uint16_type nDim = Left::nDim;
    static const bool is_diagonalized = Left::is_diagonalized && Right::is_diagonalized;
    static const bool is_transposed = Left::is_transposed && Right::is_transposed;

    typedef typename mpl::if_<mpl::bool_<Left::is_scalar>,
            mpl::identity<Shape<nDim,Scalar,is_transposed,is_diagonalized> >,
            typename mpl::if_<mpl::bool_<Left::is_vectorial>,
            mpl::identity<Shape<nDim,Vectorial,is_transposed,is_diagonalized> >,
            typename mpl::if_<mpl::bool_<Left::is_tensor2>,
            mpl::identity<Shape<nDim,Tensor2,is_transposed,is_diagonalized> >,
            mpl::identity<Shape<nDim,Tensor3,is_transposed,is_diagonalized> >
            >::type >::type >::type::type type;

    static const int op = 0;

    template<bool left_is_zero, bool right_is_zero>
    struct is_zero
    {
        static const bool value = ( left_is_zero||right_is_zero );
        static const bool update_and_eval_left = !value;
        static const bool update_and_eval_right = !value;
    };
};

struct INVALID_MULTIPLICATION {};
template<int D, uint16_type M, uint16_type N>
struct mn_to_shape
{
    typedef typename mpl::if_<mpl::greater<mpl::int_<M>,
            mpl::int_<N> >,
            mpl::identity<Shape<D, Vectorial, false, false> >,
            typename mpl::if_<mpl::greater<mpl::int_<N>,
            mpl::int_<M> >,
            mpl::identity<Shape<D, Vectorial, true, false> >,
            typename mpl::if_<mpl::equal_to<mpl::int_<N>,
            mpl::int_<1> >,
            mpl::identity<Shape<D, Scalar, false,false> >,
            mpl::identity<Shape<D, Tensor2, false,false> > >::type>::type>::type::type type;


};

template<typename Left, typename Right>
struct shape_op_mul
{

    typedef typename mpl::if_<mpl::bool_<Left::is_scalar>,
            mpl::identity<Right>,
            typename mpl::if_<mpl::bool_<Right::is_scalar>,
            mpl::identity<Left>,

            // check that Left::N == Right::M
            typename mpl::if_<mpl::equal_to<mpl::int_<Left::N>,
            mpl::int_<Right::M> >,
            mpl::identity<typename mn_to_shape<Left::nDim,
            Left::M,
            Right::N >::type >,
            mpl::identity<typename shape_op_id<Left, Right>::type> >::type >::type>::type::type type;

    static const int op = mpl::if_<mpl::or_<mpl::bool_<Left::is_scalar>,
                     mpl::bool_<Right::is_scalar> >,
                     mpl::int_<0>,
                     mpl::int_<1> >::type::value;
    template<bool left_is_zero, bool right_is_zero>
    struct is_zero
    {
        static const bool value = ( left_is_zero||right_is_zero );
        static const bool update_and_eval_left = !value;
        static const bool update_and_eval_right = !value;
    };
};


template<typename Left, typename Right>
struct shape_op_div
{

    BOOST_MPL_ASSERT_MSG( mpl::bool_<Right::is_scalar>::value,
                          INVALID_TENSOR_RANK_FOR_DIVISION_SHOULD_BE_0,
                          ( Left, Right ) );

    typedef typename mpl::if_<mpl::bool_<Right::is_scalar>,
            mpl::identity<mpl::identity<Left> >,
            mpl::identity<shape_op_id<Left, Right> > >::type::type::type type;

    static const int op = 0;

    template<bool left_is_zero, bool right_is_zero>
    struct is_zero
    {
        //BOOST_MPL_ASSERT_MSG( (!right_is_zero), (INVALID_OPERATION_DIV_BY_ZERO), (left_is_zero,right_is_zero));
        static const bool value = left_is_zero;
        static const bool update_and_eval_left = !value;
        static const bool update_and_eval_right = !value;
    };
};

template<uint16_type Dim, template<uint16_type D> class Type, bool is_transposed, bool diag>
std::ostream& operator<<( std::ostream& os, Shape<Dim, Type, is_transposed, diag> const& shape )
{
    os << "shape = [ndim = " << shape.nDim << ", rank = " << shape.rank << ", M = " << shape.M << ", N = " << shape.N << "]";

    return os;
}
/// \endcond
} // vf
} // Feel
#endif /* __Shape_H */
