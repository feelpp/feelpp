/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-13

  Copyright (C) 2014 Feel++ Consortium

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
   \file minmax.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_MINMAX_HPP
#define FEELPP_VF_MINMAX_HPP 1

namespace Feel
{
namespace vf
{
template < typename ExprT1, typename ExprT2 >
class OpMax
{
public:
    static const size_type context = ExprT1::context | ExprT2::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ( ExprT1::imorder<ExprT2::imorder )*ExprT2::imorder + ( ExprT1::imorder>=ExprT2::imorder )*ExprT1::imorder;
    static const bool imIsPoly = ExprT1::imIsPoly && ExprT2::imIsPoly;

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

    typedef OpMax<ExprT1, ExprT2> this_type;
    typedef ExprT1 expression_1_type;
    typedef ExprT2 expression_2_type;
    typedef typename strongest_numeric_type<typename expression_1_type::value_type,
            typename expression_2_type::value_type>::type value_type;
    typedef value_type evaluate_type;
    explicit OpMax( expression_1_type const& __expr1, expression_2_type const& __expr2  )
        :
        M_expr_1( __expr1 ),
        M_expr_2( __expr2 )
    {
        DVLOG(2) << "OpMax::OpMax default constructor\n";
    }

    OpMax( OpMax const& __vfp  )
        :
        M_expr_1( __vfp.M_expr_1 ),
        M_expr_2( __vfp.M_expr_2 )
    {
        DVLOG(2) << "OpMax::OpMax copy constructor\n";
    }

    expression_1_type const& left() const
    {
        return M_expr_1;
    }
    expression_2_type const& right() const
    {
        return M_expr_2;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_type;
        typedef typename expression_2_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_type;

        typedef typename strongest_numeric_type<typename l_type::value_type,
                typename r_type::value_type>::type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        BOOST_MPL_ASSERT_MSG( ( boost::is_same<typename l_type::shape, typename r_type::shape>::value ),
                              INVALID_SHAPES_FOR_MIN,
                              ( mpl::int_<l_type::shape::M>,mpl::int_<l_type::shape::N>,mpl::int_<r_type::shape::M>,mpl::int_<r_type::shape::N> ) );
        typedef typename l_type::shape shape;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev, feu ),
            M_right( expr.right(), geom, fev, feu )
        {
        }
        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev ),
            M_right( expr.right(), geom, fev )
        {
        }
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom ),
            M_right( expr.right(), geom )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            M_left.init( im );
            M_right.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
            M_left.update( geom );
            M_right.update( geom );

        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
            M_left.update( geom, face );
            M_right.update( geom, face );

        }

        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            return evalq( c1, c2, q );
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );

        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type left = M_left.evalq( c1, c2, q );
            value_type right = M_right.evalq( c1, c2, q );
            return std::max( left, right );
        }

    private:
        gmc_ptrtype M_gmc;
        l_type M_left;
        r_type M_right;
    };

protected:
    OpMax() {}

    expression_1_type M_expr_1;
    expression_2_type M_expr_2;
};

template<typename ExprT1, typename ExprT2>
inline
Expr< OpMax<typename mpl::if_<boost::is_arithmetic<ExprT1>,
      mpl::identity<Cst<ExprT1> >,
      mpl::identity<ExprT1> >::type::type,
      typename mpl::if_<boost::is_arithmetic<ExprT2>,
      mpl::identity<Cst<ExprT2> >,
      mpl::identity<ExprT2> >::type::type
      > >
      max( ExprT1 const& __e1, ExprT2 const& __e2 )
{
    typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,
            mpl::identity<Cst<ExprT1> >,
            mpl::identity<ExprT1> >::type::type t1;
    typedef typename mpl::if_<boost::is_arithmetic<ExprT2>,
            mpl::identity<Cst<ExprT2> >,
            mpl::identity<ExprT2> >::type::type t2;
    typedef OpMax<t1, t2> expr_t;
    return Expr< expr_t >(  expr_t( t1( __e1 ), t2( __e2 ) ) );
}

template < typename ExprT1, typename ExprT2 >
class OpMin
{
public:
    static const size_type context = ExprT1::context | ExprT2::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ( ExprT1::imorder<ExprT2::imorder )*ExprT2::imorder + ( ExprT1::imorder>=ExprT2::imorder )*ExprT1::imorder;
    static const bool imIsPoly = ExprT1::imIsPoly && ExprT2::imIsPoly;

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

    typedef OpMin<ExprT1, ExprT2> this_type;
    typedef ExprT1 expression_1_type;
    typedef ExprT2 expression_2_type;
    typedef typename strongest_numeric_type<typename expression_1_type::value_type,
            typename expression_2_type::value_type>::type value_type;
    typedef value_type evaluate_type;
    explicit OpMin( expression_1_type const& __expr1, expression_2_type const& __expr2  )
        :
        M_expr_1( __expr1 ),
        M_expr_2( __expr2 )
    {
        DVLOG(2) << "OpMin::OpMin default constructor\n";
    }

    OpMin( OpMin const& __vfp  )
        :
        M_expr_1( __vfp.M_expr_1 ),
        M_expr_2( __vfp.M_expr_2 )
    {
        DVLOG(2) << "OpMin::OpMin copy constructor\n";
    }

    expression_1_type const& left() const
    {
        return M_expr_1;
    }
    expression_2_type const& right() const
    {
        return M_expr_2;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_type;
        typedef typename expression_2_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_type;

        typedef typename strongest_numeric_type<typename l_type::value_type,
                typename r_type::value_type>::type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        BOOST_MPL_ASSERT_MSG( ( boost::is_same<typename l_type::shape, typename r_type::shape>::value ),
                              INVALID_SHAPES_FOR_MIN,
                              ( mpl::int_<l_type::shape::M>,mpl::int_<l_type::shape::N>,mpl::int_<r_type::shape::M>,mpl::int_<r_type::shape::N> ) );
        typedef typename l_type::shape shape;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev, feu ),
            M_right( expr.right(), geom, fev, feu )
        {
        }
        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev ),
            M_right( expr.right(), geom, fev )
        {
        }
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom ),
            M_right( expr.right(), geom )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            M_left.init( im );
            M_right.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
            M_left.update( geom );
            M_right.update( geom );

        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
            M_left.update( geom, face );
            M_right.update( geom, face );

        }

        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return evalq( c1, c2, q );
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );

        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type left = M_left.evalq( c1, c2, q );
            value_type right = M_right.evalq( c1, c2, q );
            return std::min( left, right );
        }

    private:
        gmc_ptrtype M_gmc;
        l_type M_left;
        r_type M_right;
    };

protected:
    OpMin() {}

    expression_1_type M_expr_1;
    expression_2_type M_expr_2;
};

template<typename ExprT1, typename ExprT2>
inline
Expr< OpMin<typename mpl::if_<boost::is_arithmetic<ExprT1>,
      mpl::identity<Cst<ExprT1> >,
      mpl::identity<ExprT1> >::type::type,
      typename mpl::if_<boost::is_arithmetic<ExprT2>,
      mpl::identity<Cst<ExprT2> >,
      mpl::identity<ExprT2> >::type::type> >
      min( ExprT1 const& __e1, ExprT2 const& __e2 )
{
    typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,
            mpl::identity<Cst<ExprT1> >,
            mpl::identity<ExprT1> >::type::type t1;
    typedef typename mpl::if_<boost::is_arithmetic<ExprT2>,
            mpl::identity<Cst<ExprT2> >,
            mpl::identity<ExprT2> >::type::type t2;
    typedef OpMin<t1, t2> expr_t;
    return Expr< expr_t >(  expr_t( t1( __e1 ), t2( __e2 ) ) );
}

} // vf
} // Feel
#endif /* FEELPP_VF_MINMAX_HPP */
