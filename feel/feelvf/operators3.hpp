/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-17

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Universite Joseph Fourier (Grenoble I)

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
   \file operators3.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-17
 */
#if !defined( __FEELPP_VF_OPERATORS3_HPP )
#define __FEELPP_VF_OPERATORS3_HPP 1

# include <boost/preprocessor/stringize.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
template <typename LeftExprType, typename RightExprType>
class OpDot
{
public:

    static const size_type context = LeftExprType::context|RightExprType::context;

    static const uint16_type imorder = LeftExprType::imorder + RightExprType::imorder;
    static const bool imIsPoly = LeftExprType::imIsPoly && RightExprType::imIsPoly;

    typedef LeftExprType left_expression_type;
    typedef RightExprType right_expression_type;
    typedef OpDot<LeftExprType, RightExprType> this_type;
    typedef this_type self_type;
    typedef typename mpl::if_<mpl::greater<mpl::sizeof_<typename left_expression_type::value_type>,
            mpl::sizeof_<typename right_expression_type::value_type> >,
            mpl::identity<typename left_expression_type::value_type>,
            mpl::identity<typename right_expression_type::value_type> >::type::type value_type;

    OpDot ( left_expression_type const& left,
            right_expression_type const& right )
        :
        M_left ( left ),
        M_right ( right )
    {
        DVLOG(2) << "[" BOOST_PP_STRINGIZE( OpDot ) "] default constructor\n";

    }
    OpDot( OpDot const& op )
        :
        M_left ( op.M_left ),
        M_right ( op.M_right )
    {
        DVLOG(2) << "[" BOOST_PP_STRINGIZE( OpDot ) "] copy constructor\n";

    }

    left_expression_type const& left() const
    {
        return M_left;
    }
    right_expression_type const& right() const
    {
        return M_right;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename LeftExprType::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_type;
        typedef typename RightExprType::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_type;
        typedef typename strongest_numeric_type<typename l_type::value_type,
                typename r_type::value_type>::type value_type;
        BOOST_STATIC_ASSERT( ( boost::is_same<typename l_type::return_value_type,
                               typename r_type::return_value_type>::value ) );

        typedef typename l_type::return_value_type return_value_type;

        static const uint16_type rank = return_value_type::rank+1;
        static const uint16_type nComponents = return_value_type::nComponents;
        static const bool do_reduction =  ( !boost::is_same<return_value_type,typename Basis_i_t::polyset_type>::value||
                                            !boost::is_same<return_value_type,typename Basis_j_t::polyset_type>::value );
        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_left( expr.left(), geom, fev, feu ),
            M_right( expr.right(), geom, fev, feu )
        {}

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_left.update( geom, fev, feu );
            M_right.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_left.update( geom, fev );
            M_right.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_left.update( geom );
            M_right.update( geom );
        }


        value_type
        operator()( uint16_type i,
                    uint16_type j,
                    int q ) const
        {
            return this->operator()( i, j, q, mpl::bool_<do_reduction>() );
        }

    private:

        value_type
        operator()( uint16_type i,
                    uint16_type j,
                    int q,
                    mpl::bool_<true> ) const
        {
            value_type res( 0 );
#if 0
            std::cout << "[dot] rank = "<< rank << " nc=" << return_value_type::nComponentsLast << "\n"
                      << "[dot] ranki= "<< uint16_type::rank << " rankj = "<< uint16_type::rank << "\n";
#endif
            typename Basis_i_t::template Index<rank> i_up ( i );
            typename Basis_i_t::template Index<rank> j_up ( j );

            // fail for g++ 4.2
#if 0
            BOOST_MPL_ASSERT_MSG( rank == uint16_type::rank+1, INVALID_INDEX_RANK, ( mpl::int_<uint16_type::rank>, mpl::int_<rank> ) );
            BOOST_MPL_ASSERT_MSG( rank == uint16_type::rank+1, INVALID_INDEX_RANK, ( mpl::int_<uint16_type::rank>, mpl::int_<rank> ) );
#endif

            for ( int c = 0; c < return_value_type::nComponentsLast; ++c )
            {
                i_up.setIndex( rank-1, c );
                j_up.setIndex( rank-1, c );
                res += M_left( i_up, j_up, q ) * M_right( i_up, j_up, q );
            }

            return res;
        }

        value_type
        operator()( uint16_type i,
                    uint16_type j,
                    int q,
                    mpl::bool_<false> ) const
        {
            return M_left( i, j, q ) * M_right( i, j, q );
        }
    private:
        l_type M_left;
        r_type M_right;

    };

protected:
    OpDot () {}


private:
    left_expression_type const& M_left;
    right_expression_type const& M_right;
};
template <class LeftExprType, typename RightExprType>
inline Expr< OpDot< LeftExprType, RightExprType> >
dot( LeftExprType const& left, RightExprType const& right  )
{
    typedef OpDot<LeftExprType, RightExprType> expr_t;
    return Expr< expr_t >(  expr_t( left, right ) );
}
// template <class LeftExprType, typename RightExprType>
// inline Expr< OpDot< LeftExprType, RightExprType> >
// operator,( LeftExprType const& left, RightExprType const& right  )
// {
//     typedef OpDot<LeftExprType, RightExprType> expr_t;
//     return Expr< expr_t >(  expr_t(left, right) );
// }
/// \endcond
} // vf
} // feel

#endif /* __FEELPP_VF_OPERATORS2_HPP */
