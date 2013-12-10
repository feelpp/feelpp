/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-08-04

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file products.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-08-04
 */
#ifndef __Products_H
#define __Products_H 1

#include <Eigen/Core>

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class Product
 * \brief inner and outer products
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename ExprL, typename ExprR, int Type = 1>
class Product
{
public:

    static const size_type context = ExprL::context |  ExprR::context;
    static const bool is_terminal = false;
    static const int product_type = Type;

    static const uint16_type imorder = ExprL::imorder+ExprR::imorder;
    static const bool imIsPoly = ExprL::imIsPoly && ExprR::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result =
            ExprL::template HasTestFunction<Func>::result ||
        ExprR::template HasTestFunction<Func>::result ;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result =
            ExprL::template HasTrialFunction<Func>::result||
        ExprR::template HasTrialFunction<Func>::result ;
    };


    /** @name Typedefs
     */
    //@{

    typedef ExprL left_expression_type;
    typedef ExprR right_expression_type;
    typedef typename left_expression_type::value_type value_type;
    typedef Product<ExprL,ExprR,Type> this_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Product( left_expression_type const & left_expr,
                      right_expression_type const & right_expr )
        :
        M_left_expr( left_expr ),
        M_right_expr( right_expr )
    {}
    Product( Product const & te )
        :
        M_left_expr( te.M_left_expr ),
        M_right_expr( te.M_right_expr )
    {}
    ~Product()
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

    left_expression_type const& left() const
    {
        return M_left_expr;
    }
    right_expression_type const& right() const
    {
        return M_right_expr;
    }

    //@}

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        typedef typename left_expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_tensor_expr_type;
        typedef typename right_expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_tensor_expr_type;
        typedef typename l_tensor_expr_type::value_type value_type;

        typedef typename l_tensor_expr_type::shape left_shape;
        typedef typename r_tensor_expr_type::shape right_shape;
        typedef Shape<left_shape::nDim,Scalar,false,false> shape;
        static const bool l_is_terminal = left_expression_type::is_terminal;
        static const bool r_is_terminal = right_expression_type::is_terminal;
        template <class Args> struct sig
        {
            typedef value_type type;
        };

        BOOST_MPL_ASSERT_MSG( (left_shape::M == right_shape::M) && (left_shape::N == right_shape::N) ,
                              INVALID_RANK_LEFT_AND_RIGHT_SHOULD_BE_THE_SAME,
                              (mpl::int_<left_shape::M>,mpl::int_<right_shape::M>,
                               mpl::int_<left_shape::N>,mpl::int_<right_shape::N>));

        struct is_zero
        {
            static const bool value = l_tensor_expr_type::is_zero::value || r_tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_l_tensor_expr( expr.left(), geom, fev, feu ),
            M_r_tensor_expr( expr.right(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_l_tensor_expr( expr.left(), geom, fev ),
            M_r_tensor_expr( expr.right(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_l_tensor_expr( expr.left(), geom ),
            M_r_tensor_expr( expr.right(), geom )
        {
        }

        template<typename IM>
        void init( IM const& im )
        {
            M_l_tensor_expr.init( im );
            M_r_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_l_tensor_expr.update( geom, fev, feu );
            M_r_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_l_tensor_expr.update( geom, fev );
            M_r_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_l_tensor_expr.update( geom );
            M_r_tensor_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_l_tensor_expr.update( geom, face );
            M_r_tensor_expr.update( geom, face );
        }
        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
            M_l_tensor_expr.updateContext( ctx );
            M_r_tensor_expr.updateContext( ctx );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            return evalijq( i, j, cc1, cc2, q, typename mpl::and_<mpl::bool_<l_is_terminal>,mpl::bool_<l_is_terminal> >::type() );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::bool_<0> ) const
        {
            double res = 0;

            if ( Type == 1 )
            {
                for ( int c2 = 0; c2 < left_shape::N; ++ c2 )
                    for ( int c1 = 0; c1 < left_shape::M; ++ c1 )
                    {
                        res += M_l_tensor_expr.evalijq( i, j, c1, c2, q )*M_r_tensor_expr.evalijq( i, j, c1, c2, q );
                    }
            }

            return res;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::bool_<1> ) const
        {
            //if ( Type == 1 )
            {
                return ( M_l_tensor_expr.evalijq( i,j,q ).adjoint()*M_r_tensor_expr.evalijq( i,j,q ) ).trace();
            }
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            double res = 0;

            if ( Type == 1 )
            {
                for ( int c2 = 0; c2 < left_shape::N; ++ c2 )
                    for ( int c1 = 0; c1 < left_shape::M; ++ c1 )
                    {
                        res += M_l_tensor_expr.evaliq( i, c1, c2, q )*M_r_tensor_expr.evaliq( i, c1, c2, q );
                    }
            }

            return res;
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            double res = 0;

            if ( Type == 1 )
            {
                for ( int c2 = 0; c2 < left_shape::N; ++ c2 )
                    for ( int c1 = 0; c1 < left_shape::M; ++ c1 )
                    {
                        res += M_l_tensor_expr.evalq( c1, c2, q )*M_r_tensor_expr.evalq( c1, c2, q );
                    }
            }

            return res;
        }

    private:
        l_tensor_expr_type M_l_tensor_expr;
        r_tensor_expr_type M_r_tensor_expr;
        mutable Eigen::Matrix<value_type,left_shape::M,left_shape::N> M1,M2;
    };

private:
    mutable left_expression_type  M_left_expr;
    mutable right_expression_type  M_right_expr;
};
/// \endcond

/**
 * \brief symetric part of a matricial expression
 */
template<typename ExprL, typename ExprR>
inline
Expr< Product<ExprL, ExprR,1> >
inner( ExprL l, ExprR r )
{
    typedef Product<ExprL, ExprR,1> product_t;
    return Expr< product_t >(  product_t( l, r ) );
}

/// \cond detail
/**
 * \class CrossProduct
 * \brief cross products
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename ExprL, typename ExprR>
class CrossProduct
{
public:

    static const size_type context = ExprL::context |  ExprR::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ExprL::imorder+ExprR::imorder;
    static const bool imIsPoly = ExprL::imIsPoly && ExprR::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result =
            ExprL::template HasTestFunction<Func>::result ||
        ExprR::template HasTestFunction<Func>::result ;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result =
            ExprL::template HasTrialFunction<Func>::result||
        ExprR::template HasTrialFunction<Func>::result ;
    };


    /** @name Typedefs
     */
    //@{

    typedef ExprL left_expression_type;
    typedef ExprR right_expression_type;
    typedef typename left_expression_type::value_type value_type;
    typedef CrossProduct<ExprL,ExprR> this_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit CrossProduct( left_expression_type const & left_expr,
                      right_expression_type const & right_expr )
        :
        M_left_expr( left_expr ),
        M_right_expr( right_expr )
    {}
    CrossProduct( CrossProduct const & te )
        :
        M_left_expr( te.M_left_expr ),
        M_right_expr( te.M_right_expr )
    {}
    ~CrossProduct()
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

    left_expression_type const& left() const
    {
        return M_left_expr;
    }
    right_expression_type const& right() const
    {
        return M_right_expr;
    }

    //@}

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        typedef typename left_expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_tensor_expr_type;
        typedef typename right_expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_tensor_expr_type;
        typedef typename l_tensor_expr_type::value_type value_type;

        typedef typename l_tensor_expr_type::shape left_shape;
        typedef typename r_tensor_expr_type::shape right_shape;
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<left_shape::nDim>,mpl::int_<2> >,
                                  mpl::identity<Shape<2,Scalar,false,false> >,
                                  mpl::identity<Shape<3,Vectorial,false,false> > >::type::type shape;
        static const bool l_is_terminal = left_expression_type::is_terminal;
        static const bool r_is_terminal = right_expression_type::is_terminal;

        BOOST_MPL_ASSERT_MSG( ( left_shape::nDim > 1 ),
                              CROSS_INVALID_DIMENSION,
                              (mpl::int_<left_shape::nDim>,mpl::int_<right_shape::nDim>));
        BOOST_MPL_ASSERT_MSG( left_shape::nDim == right_shape::nDim,
                              CROSS_INVALID_DIMENSION_LEFT_AND_RIGHT_SHOULD_BE_THE_SAME,
                              (mpl::int_<left_shape::nDim>,mpl::int_<right_shape::nDim>));
        BOOST_MPL_ASSERT_MSG( (left_shape::M == right_shape::M) && (left_shape::N == right_shape::N) ,
                              CROSS_INVALID_RANK_LEFT_AND_RIGHT_SHOULD_BE_THE_SAME,
                              (mpl::int_<left_shape::M>,mpl::int_<right_shape::M>,
                               mpl::int_<left_shape::N>,mpl::int_<right_shape::N>));


        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = l_tensor_expr_type::is_zero::value || r_tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_l_tensor_expr( expr.left(), geom, fev, feu ),
            M_r_tensor_expr( expr.right(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_l_tensor_expr( expr.left(), geom, fev ),
            M_r_tensor_expr( expr.right(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_l_tensor_expr( expr.left(), geom ),
            M_r_tensor_expr( expr.right(), geom )
        {
        }

        template<typename IM>
        void init( IM const& im )
        {
            M_l_tensor_expr.init( im );
            M_r_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_l_tensor_expr.update( geom, fev, feu );
            M_r_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_l_tensor_expr.update( geom, fev );
            M_r_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_l_tensor_expr.update( geom );
            M_r_tensor_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_l_tensor_expr.update( geom, face );
            M_r_tensor_expr.update( geom, face );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            return evalijq( i, j, cc1, cc2, q, mpl::int_<left_shape::nDim>() );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<2> ) const
        {
            double res =  M_l_tensor_expr.evalijq( i, j, 0, 0, q )*M_r_tensor_expr.evalijq( i, j, 1, 0, q );
            res -= M_l_tensor_expr.evalijq( i, j, 1, 0, q )*M_r_tensor_expr.evalijq( i, j, 0, 0, q );
            return res;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<3> ) const
        {
            double res =  M_l_tensor_expr.evalijq( i, j, (cc1+1)%3, 0, q )*M_r_tensor_expr.evalijq( i, j, (cc1+2)%3, 0, q );
            res -= M_l_tensor_expr.evalijq( i, j, (cc1+2)%3, 0, q )*M_r_tensor_expr.evalijq( i, j, (cc1+1)%3, 0, q );
            return res;
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            return evaliq( i, cc1, cc2, q, mpl::int_<left_shape::nDim>() );
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<2> ) const
        {
            double res =  M_l_tensor_expr.evaliq( i, 0, 0, q )*M_r_tensor_expr.evaliq( i, 1, 0, q );
            res -= M_l_tensor_expr.evaliq( i, 1, 0, q )*M_r_tensor_expr.evaliq( i, 0, 0, q );
            return res;
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<3> ) const
        {
            double res =  M_l_tensor_expr.evaliq( i, (cc1+1)%3, 0, q )*M_r_tensor_expr.evaliq( i, (cc1+2)%3, 0, q );
            res -= M_l_tensor_expr.evaliq( i, (cc1+2)%3, 0, q )*M_r_tensor_expr.evaliq( i, (cc1+1)%3, 0, q );
            return res;
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            return evalq( cc1, cc2, q, mpl::int_<left_shape::nDim>() );
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<2> ) const
        {
            double res =  M_l_tensor_expr.evalq( 0, 0, q )*M_r_tensor_expr.evalq( 1, 0, q );
            res -= M_l_tensor_expr.evalq( 1, 0, q )*M_r_tensor_expr.evalq( 0, 0, q );
            return res;
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<3> ) const
        {
            double res =  M_l_tensor_expr.evalq( (cc1+1)%3, 0, q )*M_r_tensor_expr.evalq( (cc1+2)%3, 0, q );
            res -= M_l_tensor_expr.evalq( (cc1+2)%3, 0, q )*M_r_tensor_expr.evalq( (cc1+1)%3, 0, q );
            return res;
        }

    private:
        l_tensor_expr_type M_l_tensor_expr;
        r_tensor_expr_type M_r_tensor_expr;
    };

private:
    mutable left_expression_type  M_left_expr;
    mutable right_expression_type  M_right_expr;
};
/// \endcond

/**
 * \brief symetric part of a matricial expression
 */
template<typename ExprL, typename ExprR>
inline
Expr< CrossProduct<ExprL, ExprR> >
cross( ExprL l, ExprR r )
{
    typedef CrossProduct<ExprL, ExprR> product_t;
    return Expr< product_t >(  product_t( l, r ) );
}

} // vf


} // Feel
#endif /* __Products_H */
