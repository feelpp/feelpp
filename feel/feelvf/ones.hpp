/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-04-11

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file ones.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-04-11
 */
#ifndef __Ones_H
#define __Ones_H 1

//#include <blitz/array.h>
#include <boost/multi_array.hpp>


namespace Feel
{
namespace vf
{
/// \cond detail
namespace detail
{
/**
 * \class Ones
 * \brief Return a matrix expression or N-dimensional array whose elements are all 1
 *
 *If you need to create a matrix whose values are all the same, you
 * should use an expression like
 *
 *\code
 * val * ones<n,m>()
 *\endcode
 *
 * @author Christophe Prud'homme
 */
template<int M, int N>
class Ones
{
public:


    /** @name Typedefs
     */
    //@{
    static const size_type context = 0;

    static const uint16_type imorder = 0;
    static const bool imIsPoly = true;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = false;
    };

    typedef Ones<M,N> this_type;
    typedef double value_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Ones()
        :
        _M_ones( boost::extents[M][N] )
    {

        //_M_ones = 1;
        for ( int i = 0; i < M; ++i )
            for ( int j = 0; j < N; ++j )
                _M_ones[i][j] = 1;
    }

    Ones( Ones const & eig )
        :
        _M_ones( eig._M_ones )
    {}
    ~Ones()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //blitz::Array<value_type,2> ones() const { return _M_ones; }
    boost::multi_array<value_type,2> const& ones() const
    {
        return _M_ones;
    }

    //@}
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_type::value_type value_type;
        typedef value_type return_value_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >,
                mpl::identity<detail::gmc<0> >,
                mpl::identity<detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        struct INVALID_SHAPE {};
        static const bool eq11 = ( M==1 )&&( N==1 );
        static const bool eqD1 = ( M==gmc_type::nDim )&&( N==1 );
        static const bool eq1D = ( M==1 )&&( N==gmc_type::nDim );
        static const bool eqDD = ( M==gmc_type::nDim )&&( N==gmc_type::nDim );
        typedef typename mpl::if_< mpl::bool_<eq11>,
                mpl::identity<Shape<gmc_type::nDim, Scalar, false, false> >,
                typename mpl::if_< mpl::bool_<eqD1>,
                mpl::identity<Shape<gmc_type::nDim, Vectorial, false, false> >,
                typename mpl::if_< mpl::bool_<eq1D>,
                mpl::identity<Shape<gmc_type::nDim, Vectorial, true, false> >,
                typename mpl::if_< mpl::bool_<eqDD>,
                mpl::identity<Shape<gmc_type::nDim, Tensor2, false, false> >,
                mpl::identity<INVALID_SHAPE> >::type>::type>::type>::type::type shape;


        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,Geo_t const&, Basis_i_t const&, Basis_j_t const& )
            :
            _M_expr( expr )
        {
            //std::cout << "tensor::ones = " << _M_expr.ones() << "\n";
        }

        tensor( this_type const& expr,Geo_t const&, Basis_i_t const& )
            :
            _M_expr( expr )
        {
        }

        tensor( this_type const& expr, Geo_t const&  )
            :
            _M_expr( expr )
        {
        }

        void update( Geo_t const&, Basis_i_t const&, Basis_j_t const& )
        {
        }
        void update( Geo_t const&, Basis_i_t const& )
        {
        }
        void update( Geo_t const& )
        {
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            Feel::detail::ignore_unused_variable_warning( q );
            return eval( c1, c2, mpl::int_<shape::rank>() );
        }

        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            Feel::detail::ignore_unused_variable_warning( q );
            return eval( c1, c2, mpl::int_<shape::rank>() );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( q );
            return eval( c1, c2, mpl::int_<shape::rank>() );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( q );
            return eval( c1, c2, mpl::int_<shape::rank>() );
        }
    private:
        value_type
        eval( int c1, int c2, mpl::int_<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return _M_expr.ones()[0][0];
        }
        value_type
        eval( int c1, int c2, mpl::int_<1> ) const
        {
            if ( shape::is_transposed )
                return _M_expr.ones()[0][c2];

            return _M_expr.ones()[c1][0];
        }
        value_type
        eval( int c1, int c2, mpl::int_<2> ) const
        {
            return _M_expr.ones()[c1][c2];
        }
        this_type _M_expr;
    };
private:
    //blitz::Array<value_type,2> _M_ones;
    boost::multi_array<value_type,2> _M_ones;

};
} // detail
/// \endcond

/**
 * \class Ones
 * \brief Return a matrix expression or N-dimensional array whose elements are all 1
 *
 *If you need to create a matrix whose values are all the same, you
 * should use an expression like
 *
 *\code
 * val * ones<n,m>()
 *\endcode
 *
 * @author Christophe
 */
template<int M, int N>
inline
Expr<detail::Ones<M,N> >
ones()
{
    return Expr< detail::Ones<M,N> >(  detail::Ones<M, N>() );
}
} // vf
} // Feel
#endif /* __Ones_H */
