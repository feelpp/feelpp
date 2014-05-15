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
   \file pow.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_POW_HPP
#define FEELPP_VF_POW_HPP 1

namespace Feel
{
namespace vf
{
template < typename ExprT1, typename ExprT2 >
class Pow
{
public:

    static const size_type context = ExprT1::context|ExprT2::context;
    static const bool is_terminal = false;

    /**
     * \warning the Pow order computation is wrong here, we actually need the
     * ExprT2 value (and not imorder) to multiply by ExprT1::imorder.
     */
    static const uint16_type imorder = ExprT1::imorder;
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

    typedef Pow<ExprT1, ExprT2> this_type;
    typedef ExprT1 expression_1_type;
    typedef ExprT2 expression_2_type;
    typedef typename expression_1_type::value_type value_1_type;
    typedef typename expression_2_type::value_type value_2_type;
    typedef value_1_type value_type;
    typedef value_type evaluate_type;

    // verify that all returning types are integral or floating types
    BOOST_STATIC_ASSERT( ::boost::is_arithmetic<value_1_type>::value  &&
                         ::boost::is_arithmetic<value_2_type>::value );

    explicit Pow( expression_1_type const& __expr1, expression_2_type const& __expr2  )
        :
        M_expr_1( __expr1 ),
        M_expr_2( __expr2 )
    {
        DVLOG(2) << "Pow::Pow default constructor\n";
    }

    Pow( Pow const& __vfp  )
        :
        M_expr_1( __vfp.M_expr_1 ),
        M_expr_2( __vfp.M_expr_2 )
    {
        DVLOG(2) << "Pow::Pow copy constructor\n";
    }

    bool isSymetric() const
    {
        return false;
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
        typedef typename l_type::shape shape;

        typedef typename Eigen::Matrix<value_type,shape::M,shape::N> loc_type;

        struct is_zero
        {
            static const bool value = l_type::is_zero::value;
        };

        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev, feu ),
            M_right( expr.right(), geom, fev, feu ),
            M_loc( boost::extents[M_gmc->nPoints()]  )
        {
            update( geom );
        }
        tensor( expression_type const& expr, Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom, fev ),
            M_right( expr.right(), geom, fev ),
            M_loc( boost::extents[M_gmc->nPoints()] )
        {
            update( geom );
        }
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_left( expr.left(),  geom ),
            M_right( expr.right(), geom ),
            M_loc(  boost::extents[M_gmc->nPoints()] )
        {
            update( geom );
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
            const int npts =  fusion::at_key<key_type>( geom ).get()->nPoints();
            M_left.update( geom );
            M_right.update( geom );

            for ( int q = 0; q < npts; ++q )
                for ( int c1 = 0; c1 < shape::M; ++c1 )
                    for ( int c2 = 0; c2 < shape::N; ++c2 )
                    {
                        value_type left = M_left.evalq( c1, c2, q );
                        value_type right = M_right.evalq( c1, c2, q );
                        M_loc[q]( c1,c2 ) = std::pow( left, right );
                    }
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
            M_left.update( geom, face );
            M_right.update( geom, face );

            for ( int q = 0; q < M_gmc->nPoints(); ++q )
                for ( int c1 = 0; c1 < shape::M; ++c1 )
                    for ( int c2 = 0; c2 < shape::N; ++c2 )
                    {
                        value_type left = M_left.evalq( c1, c2, q );
                        value_type right = M_right.evalq( c1, c2, q );
                        M_loc[q]( c1,c2 ) = std::pow( left, right );
                    }
        }

        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q  ) const
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
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q  ) const
        {
            return evalq( c1, c2, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q, mpl::int_<shape::rank>() );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return M_loc[q]( 0,0 );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1> ) const
        {
            if ( shape::M > shape::N )
                return M_loc[q]( c1,0 );

            return M_loc[q]( 0,c2 );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> ) const
        {
            return M_loc[q]( c1,c2 );
        }

    private:
        gmc_ptrtype M_gmc;
        l_type M_left;
        r_type M_right;
        //ublas::vector<double> M_loc;
        boost::multi_array<loc_type,1> M_loc;
    };

protected:
    Pow() {}

    expression_1_type M_expr_1;
    expression_2_type M_expr_2;
};

template<typename ExprT1,  typename ExprT2>
inline
Expr< Pow<typename mpl::if_<boost::is_arithmetic<ExprT1>,
      mpl::identity<Cst<ExprT1> >,
      mpl::identity<ExprT1> >::type::type,
      typename mpl::if_<boost::is_arithmetic<ExprT2>,
      mpl::identity<Cst<ExprT2> >,
      mpl::identity<ExprT2> >::type::type> >
      pow( ExprT1 const& __e1, ExprT2 const& __e2 )
{
    typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,
            mpl::identity<Cst<ExprT1> >,
            mpl::identity<ExprT1> >::type::type t1;
    typedef typename mpl::if_<boost::is_arithmetic<ExprT2>,
            mpl::identity<Cst<ExprT2> >,
            mpl::identity<ExprT2> >::type::type t2;
    typedef Pow<t1, t2> expr_t;
    return Expr< expr_t >(  expr_t( t1( __e1 ), t2( __e2 ) ) );
}


template<typename ExprT1,  typename ExprT2>
inline
Expr< Pow<typename mpl::if_<boost::is_arithmetic<ExprT1>,
                            mpl::identity<Cst<ExprT1> >,
                            mpl::identity<Expr<ExprT1> > >::type::type,
          typename mpl::if_<boost::is_arithmetic<ExprT2>,
                            mpl::identity<Cst<ExprT2> >,
                            mpl::identity<Expr<ExprT2> > >::type::type> >
operator^( typename mpl::if_<boost::is_arithmetic<ExprT1>,
                             mpl::identity<ExprT1>,
                             mpl::identity<Expr<ExprT1> > >::type::type const& __e1,
           typename mpl::if_<boost::is_arithmetic<ExprT2>,
                             mpl::identity<ExprT2>,
                             mpl::identity<Expr<ExprT2> > >::type::type const& __e2 )
{
    typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,
            mpl::identity<Cst<ExprT1> >,
            mpl::identity<ExprT1> >::type::type t1;
    typedef typename mpl::if_<boost::is_arithmetic<ExprT2>,
            mpl::identity<Cst<ExprT2> >,
            mpl::identity<ExprT2> >::type::type t2;
    typedef Pow<t1, t2> expr_t;
    return Expr< expr_t >(  expr_t( t1( __e1 ), t2( __e2 ) ) );
}
} // vf
} //Feel
#endif /* FEELPP_VF_POW_HPP */
