/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-12-03

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2011-2016 Feel++ Consortium

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
#ifndef FEELPP_FEELPOLY_POLICY_HPP
#define FEELPP_FEELPOLY_POLICY_HPP 1


#include <Eigen/Core>

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/traits.hpp>
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/concepts.hpp>
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
struct Scalar : public ScalarBase
{
    static constexpr uint16_type rank = 0;
    static constexpr uint16_type nDim = Dim;

    static constexpr bool is_scalar = true;
    static constexpr bool is_vectorial = false;
    static constexpr bool is_tensor2 = false;
    static constexpr bool is_tensor3 = false;

    static constexpr uint16_type nComponents = 1;
    static constexpr uint16_type nComponents1 = 1;
    static constexpr uint16_type nComponents2 = 1;
    static constexpr uint16_type nComponents3 = 1;
    static constexpr uint16_type nComponentsLast = 1;

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
    template<typename Derived>
    requires EigenMatrix<Derived>
    static
    inline Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    toMatrix( Eigen::MatrixBase<Derived> const& __c )
    {
        return __c.eval();
    }
    template<typename Derived>
    requires EigenMatrix<Derived>
    static
    inline Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    toType( Eigen::MatrixBase<Derived> const& __c )
    {
        return __c.eval();
    }
};


/**
 * Policy for \c Vectorial polynomials or polynomial sets of dimension
 * \p Dim
 * \note \c Vectorial can be seen as rank 1 Tensor polynomials
 */
template<uint16_type Dim>
struct Vectorial : public VectorialBase
{
    static constexpr uint16_type rank = 1;
    static constexpr uint16_type nDim = Dim;

    static constexpr bool is_scalar = false;
    static constexpr bool is_vectorial = true;
    static constexpr bool is_tensor2 = false;
    static constexpr bool is_tensor3 = false;

    static constexpr uint16_type nComponents = nDim;
    static constexpr uint16_type nComponents1 = nDim;
    static constexpr uint16_type nComponents2 = 1;
    static constexpr uint16_type nComponents3 = 1;
    static constexpr uint16_type nComponentsLast = nComponents1;
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
    template<typename Derived>
    requires EigenMatrix<Derived>
    static
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    toMatrix( Eigen::MatrixBase<Derived> const& __c )
    {
        using value_type = typename Derived::Scalar;
        const Eigen::Index comp = static_cast<Eigen::Index>( nComponents );
        const Eigen::Index nRows = __c.rows() / comp;
        const Eigen::Index nCols = __c.cols();
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> __c_reshaped(
            nRows, comp * nCols );

        for ( Eigen::Index c = 0; c < comp; ++c )
        {
            for ( Eigen::Index r = 0; r < nRows; ++r )
            {
                __c_reshaped.block( r, c * nCols, 1, nCols ) =
                    __c.row( r * comp + c );
            }
        }

        return __c_reshaped;
    }

    template<typename Derived>
    requires EigenMatrix<Derived>
    static
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    toType( Eigen::MatrixBase<Derived> const& __c )
    {
        using value_type = typename Derived::Scalar;
        const Eigen::Index comp = static_cast<Eigen::Index>( nComponents );
        const Eigen::Index inRows = __c.rows();
        const Eigen::Index inCols = __c.cols();
        const Eigen::Index outRows = inRows * comp;
        const Eigen::Index outCols = inCols / comp;
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> __c_reshaped(
            outRows, outCols );

        for ( Eigen::Index c = 0; c < comp; ++c )
        {
            for ( Eigen::Index r = 0; r < inRows; ++r )
            {
                __c_reshaped.block( r * comp + c, 0, 1, outCols ) =
                    __c.block( r, c * outCols, 1, outCols );
            }
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
    static constexpr uint16_type rank = ( M > 1 );
    static constexpr uint16_type nDim = N;
    static constexpr uint16_type nVariables = N;

    static constexpr bool is_scalar = ( M==1 );
    static constexpr bool is_vectorial = ( N==M );
    static constexpr bool is_tensor2 = false;
    static constexpr bool is_tensor3 = false;

    static constexpr uint16_type nComponents = M;
    static constexpr uint16_type nComponents1 = M;
    static constexpr uint16_type nComponents2 = 1;
    static constexpr uint16_type nComponents3 = 1;
    static constexpr uint16_type nComponentsLast = 1;

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
    template<typename Derived>
    requires EigenMatrix<Derived>
    static
    inline Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    toMatrix( Eigen::MatrixBase<Derived> const& __c )
    {
        return __c.eval();
    }
    template<typename Derived>
    requires EigenMatrix<Derived>
    static
    inline Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    toType( Eigen::MatrixBase<Derived> const& __c )
    {
        return __c.eval();
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
struct Tensor2 : public Tensor2Base
{
    static constexpr uint16_type rank = 2;
    static constexpr uint16_type nDim = Dim;

    static constexpr bool is_scalar = false;
    static constexpr bool is_vectorial = false;
    static constexpr bool is_tensor2 = true;
    static constexpr bool is_tensor3 = false;

    static constexpr uint16_type nComponents = nDim*nDim;
    static constexpr uint16_type nComponents1 = nDim;
    static constexpr uint16_type nComponents2 = nDim;
    static constexpr uint16_type nComponents3 = 1;
    static constexpr uint16_type nComponentsLast = nComponents2;

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
                                ublas::range( c2*nCols, ( c2+1 )*nCols ) ) =
                    ublas::project( __c,
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
    template<typename Derived>
    requires EigenMatrix<Derived>
    static
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    toMatrix( Eigen::MatrixBase<Derived> const& __c )
    {
        using value_type = typename Derived::Scalar;
        const Eigen::Index comp = static_cast<Eigen::Index>( nComponents );
        const Eigen::Index nRows = __c.rows();
        const Eigen::Index nCols = __c.cols();
        const Eigen::Index nRows1 = nCols * comp;
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> __c_reshaped(
            nRows / comp, nCols * comp );

        for ( Eigen::Index c1 = 0; c1 < comp; ++c1 )
        {
            const Eigen::Index i1 = nRows1 * c1;
            for ( Eigen::Index c2 = 0; c2 < comp; ++c2 )
            {
                for ( Eigen::Index r = 0; r < nRows1 / comp; ++r )
                {
                    const Eigen::Index srcRow = i1 + c2 + r * comp;
                    const Eigen::Index dstRow = i1 / comp + r;
                    __c_reshaped.block( dstRow, c2 * nCols, 1, nCols ) =
                        __c.row( srcRow );
                }
            }
        }

        return __c_reshaped;
    }

    template<typename Derived>
    requires EigenMatrix<Derived>
    static
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    toType( Eigen::MatrixBase<Derived> const& __c )
    {
        using value_type = typename Derived::Scalar;
        const Eigen::Index comp = static_cast<Eigen::Index>( nComponents );
        const Eigen::Index inRows = __c.rows();
        const Eigen::Index inCols = __c.cols();
        const Eigen::Index outRows = inRows * comp;
        const Eigen::Index nRows1 = inCols;
        const Eigen::Index outCols = inCols / comp;
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> __c_reshaped(
            outRows, outCols );

        for ( Eigen::Index c1 = 0; c1 < comp; ++c1 )
        {
            const Eigen::Index i1 = nRows1 * c1;
            for ( Eigen::Index c2 = 0; c2 < comp; ++c2 )
            {
                for ( Eigen::Index r = 0; r < nRows1 / comp; ++r )
                {
                    const Eigen::Index srcRow = i1 / comp + r;
                    const Eigen::Index dstRow = i1 + c2 + r * comp;
                    __c_reshaped.block( dstRow, 0, 1, outCols ) =
                        __c.block( srcRow, c2 * outCols, 1, outCols );
                }
            }
        }

        return __c_reshaped;
    }
};

struct Tensor2SymmBase : Tensor2Base {};

/**
 * Policy for symmetric rank 2 tensor polynomials or polynomial sets of
 * dimension \p Dim
 *
 */
template<uint16_type Dim>
struct Tensor2Symm : public Tensor2SymmBase
{
    static constexpr uint16_type rank = 2;
    static constexpr uint16_type nDim = Dim;

    static constexpr bool is_scalar = false;
    static constexpr bool is_vectorial = false;
    static constexpr bool is_tensor2 = true;
    static constexpr bool is_tensor3 = false;

    static constexpr uint16_type nComponents = nDim*nDim;
    static constexpr uint16_type nComponents1 = nDim;
    static constexpr uint16_type nComponents2 = nDim;
    static constexpr uint16_type nComponents3 = 1;
    static constexpr uint16_type nComponentsLast = nComponents2;

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
                                ublas::range( c2*nCols, ( c2+1 )*nCols ) ) =
                    ublas::project( __c,
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
    template<typename Derived>
    requires EigenMatrix<Derived>
    static
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    toMatrix( Eigen::MatrixBase<Derived> const& __c )
    {
        using value_type = typename Derived::Scalar;
        const Eigen::Index comp = static_cast<Eigen::Index>( nComponents );
        const Eigen::Index nRows = __c.rows();
        const Eigen::Index nCols = __c.cols();
        const Eigen::Index nRows1 = nCols * comp;
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> __c_reshaped(
            nRows / comp, nCols * comp );

        for ( Eigen::Index c1 = 0; c1 < comp; ++c1 )
        {
            const Eigen::Index i1 = nRows1 * c1;
            for ( Eigen::Index c2 = 0; c2 < comp; ++c2 )
            {
                for ( Eigen::Index r = 0; r < nRows1 / comp; ++r )
                {
                    const Eigen::Index srcRow = i1 + c2 + r * comp;
                    const Eigen::Index dstRow = i1 / comp + r;
                    __c_reshaped.block( dstRow, c2 * nCols, 1, nCols ) =
                        __c.row( srcRow );
                }
            }
        }

        return __c_reshaped;
    }

    template<typename Derived>
    requires EigenMatrix<Derived>
    static
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
    toType( Eigen::MatrixBase<Derived> const& __c )
    {
        using value_type = typename Derived::Scalar;
        const Eigen::Index comp = static_cast<Eigen::Index>( nComponents );
        const Eigen::Index inRows = __c.rows();
        const Eigen::Index inCols = __c.cols();
        const Eigen::Index outRows = inRows * comp;
        const Eigen::Index nRows1 = inCols;
        const Eigen::Index outCols = inCols / comp;
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> __c_reshaped(
            outRows, outCols );

        for ( Eigen::Index c1 = 0; c1 < comp; ++c1 )
        {
            const Eigen::Index i1 = nRows1 * c1;
            for ( Eigen::Index c2 = 0; c2 < comp; ++c2 )
            {
                for ( Eigen::Index r = 0; r < nRows1 / comp; ++r )
                {
                    const Eigen::Index srcRow = i1 / comp + r;
                    const Eigen::Index dstRow = i1 + c2 + r * comp;
                    __c_reshaped.block( dstRow, 0, 1, outCols ) =
                        __c.block( srcRow, c2 * outCols, 1, outCols );
                }
            }
        }

        return __c_reshaped;
    }
};

/**
 * if T derives from a Tensor2SymmBase then is_symm is mpl::bool_<true>,
 * mpl::bool_<false> otherwise
 */
template<typename T>
struct is_symm
    : bool_c<std::is_base_of_v<Tensor2SymmBase,T>>
{
};

/**
 * helper constant that is true if the T derives from Tensor2SymmBase, false
 * orherwise
 */
template<typename T>
constexpr bool is_symm_v = std::is_base_of_v<Tensor2SymmBase,T>;

template<typename T>
using is_scalar_t = bool_c<std::is_base_of_v<ScalarBase,T>>;
template<typename T>
constexpr bool is_scalar_v = std::is_base_of_v<ScalarBase,T>;
template<typename T>
using is_vectorial_t = bool_c<std::is_base_of_v<VectorialBase,T>>;
template<typename T>
constexpr bool is_vectorial_v = std::is_base_of_v<VectorialBase,T>;
template<typename T>
using is_tensor2_t = bool_c<std::is_base_of_v<Tensor2Base,T>>;
template<typename T>
constexpr bool is_tensor2_v = std::is_base_of_v<Tensor2Base,T>;

/**
 * Policy for rank 2 tensor polynomials or polynomial sets of
 * dimension \p Dim
 *
 */
template<uint16_type Dim>
struct Tensor3
{
    static constexpr uint16_type rank = 3;
    static constexpr uint16_type nDim = Dim;

    static constexpr bool is_scalar = false;
    static constexpr bool is_vectorial = false;
    static constexpr bool is_tensor2 = false;
    static constexpr bool is_tensor3 = true;

    static constexpr uint16_type nComponents = nDim*nDim*nDim;
    static constexpr uint16_type nComponents1 = nDim;
    static constexpr uint16_type nComponents2 = nDim;
    static constexpr uint16_type nComponents3 = nDim;
    static constexpr uint16_type nComponentsLast = nComponents3;
};

template<int Dim>
struct ListReturnTypes
{
    using return_types = type_list<Scalar<Dim>, Vectorial<Dim>, Tensor2<Dim>, Tensor3<Dim> >;
};
template<typename T1, typename T2>
struct ReturnSelect
{
    using type = if_t<std::is_same_v<T1, T2>,
                      T1,
                      if_t<(T1::rank > T2::rank), T1, T2>>;
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

inline constexpr int_c<GLOBAL_COMPONENT> INDEX_GLOBAL_COMPONENT{};
inline constexpr int_c<COMPONENT_IN_COMPONENT> INDEX_COMPONENT_IN_COMPONENT{};
inline constexpr int_c<GLOBAL_FUNCTION_INDEX> INDEX_GLOBAL_FUNCTION_INDEX{};
inline constexpr int_c<PER_COMPONENT_FUNCTION_INDEX> INDEX_PER_COMPONENT_FUNCTION_INDEX{};
inline constexpr int_c<COMPONENT_IN_COMPONENT_FUNCTION_INDEX> INDEX_COMPONENT_IN_COMPONENT_FUNCTION_INDEX{};
inline constexpr int_c<FUNCTION_INDEX> INDEX_FUNCTION_INDEX{};


/**
 * Get the component type out the available types
 * \code
 * typedef typename Component<Vectorial<3> >::type component_type;
 * // component_type should be of type \c Scalar<3>
 * \endcode
 */
template<typename T>
struct GetComponent
{
    static constexpr uint16_type nDim = T::nDim;
    using type = if_t<std::is_same_v<T, Scalar<nDim> >,
                      Scalar<nDim>,
                      if_t<std::is_same_v<T, Vectorial<nDim> >,
                           Vectorial<nDim>,
                           if_t<std::is_same_v<T, Tensor2<nDim> >,
                                Tensor2<nDim>,
                                if_t<std::is_same_v<T, Tensor2Symm<nDim> >,
                                     Tensor2Symm<nDim>,
                                     mp::mp_void<>>>>>;
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
    static constexpr uint16_type nDim = T::nDim;
    using types = type_list<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim> >;
    static constexpr std::size_t index = mp::mp_find<types, T>::value;
    using type = mp::mp_at_c<types, index + 1>;
};

template<typename T>
struct RankUp2
{
    static constexpr uint16_type nDim = T::nDim;
    using types = type_list<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim> >;
    static constexpr std::size_t index = mp::mp_find<types, T>::value;
    using type = mp::mp_at_c<types, index + 2>;
};

template<typename T>
struct RankDown2
{
    static constexpr uint16_type nDim = T::nDim;
    using types = type_list<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim> >;
    static constexpr std::size_t index = mp::mp_find<types, T>::value;
    using type = mp::mp_at_c<types, index - 2>;
};

template<typename T>
struct RankSame
{
    static constexpr uint16_type nDim = T::nDim;
    typedef T type;
};
template<typename T>
struct Rank0
{
    static constexpr uint16_type nDim = T::nDim;
    typedef Scalar<nDim> type;
};
template<typename T>
struct Rank1
{
    static constexpr uint16_type nDim = T::nDim;
    typedef Vectorial<nDim> type;
};

template<typename T>
struct RankDown
{
    static constexpr uint16_type nDim = T::nDim;
    using types = if_t<std::is_base_of_v<Tensor2SymmBase,T>,
                       type_list<Scalar<nDim>, Vectorial<nDim>, Tensor2Symm<nDim>, Tensor3<nDim>>,
                       type_list<Scalar<nDim>, Vectorial<nDim>, Tensor2<nDim>, Tensor3<nDim>>>;
    static constexpr std::size_t index = mp::mp_find<types, T>::value;
    using type = mp::mp_eval_if_c<(index == 0),
                                  Scalar<nDim>,
                                  mp::mp_at,
                                  types, mp::mp_size_t<index - 1>>;
};
template<typename T>
using rankdown_t = typename RankDown<T>::type;

template<typename T>
struct RankCurl
{
    static constexpr uint16_type nDim = T::nDim;
    typedef typename if_t<( nDim == 3 ), RankSame<T>, RankDown<T> >::type type;
    static constexpr uint16_type value = ( nDim == 3 ) ? 3 : 1;
};

/**
 * Policy for orthogonal polynomial set which can be \p normalized or
 * not.  This can be used as a policy for orthogonal polynomial sets
 * in order to select either the normalized or un-normalized version
 */
template<bool normalized>
struct Normalized
{
    static constexpr bool is_normalized = normalized;
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

/**
 * Storage Policy using Eigen for numerical type \p T
 */
template<typename T>
struct StorageEigen
{
    using value_type = T;
    using vector_type = Eigen::Matrix<value_type, Eigen::Dynamic, 1>;
    using matrix_type = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_matrix_type = ublas::vector<matrix_type>;
    using vector_vector_matrix_type = ublas::vector<vector_matrix_type>;
    using matrix_node_type = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;
    using points_type = matrix_node_type;
    using node_type = Eigen::Matrix<value_type, Eigen::Dynamic, 1>;
};

} // Feel

#endif /* FEELPP_FEELPOLY_POLICY_HPP */
