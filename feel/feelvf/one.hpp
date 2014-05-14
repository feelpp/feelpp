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
   \file one.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_ONE_HPP
#define FEELPP_VF_ONE_HPP 1

namespace Feel
{
namespace vf
{
template<int CType>
class One
{
public:
    static const size_type context = 0;
    static const bool is_terminal = false;

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


    typedef One<CType> this_type;

    typedef double value_type;
    typedef value_type evaluate_type;

    One() {}
    One( One const& /*__vff*/ ) {}

    template<typename... TheExpr>
    struct Lambda
    {
        typedef this_type type;
    };
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  ) { return this_type(); }


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef Shape<gmc_type::nDim, Vectorial, false, false> shape;
        static const bool theshape = ( shape::M == gmc_type::nDim && shape::N == 1 );
        BOOST_MPL_ASSERT_MSG( theshape,
                              INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_1,
                              ( mpl::int_<shape::M>, mpl::int_<shape::N> ) );

        typedef typename expression_type::value_type value_type;

        static const uint16_type nComponents = gmc_type::nDim;
        static const int16_type vector_comp = ( CType==-1 )?1:CType;

        typedef typename mpl::if_<mpl::equal_to<mpl::int_<CType>,mpl::int_<-1> >,
                mpl::identity<ublas::scalar_vector<scalar_type> >,
                mpl::identity<ublas::unit_vector<scalar_type> > >::type::type vector_type;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& /*expr*/,
                Geo_t const& /*geom*/,
                Basis_i_t const& /*fev*/,
                Basis_j_t const& /*feu*/ )
            :
            M_one( nComponents, vector_comp )
        {
            //std::cout << "one = " << M_one << "\n";
        }
        tensor( expression_type const& /*expr*/,
                Geo_t const& /*geom*/,
                Basis_i_t const& /*fev*/ )
            :
            M_one( nComponents, vector_comp )
        {
        }
        tensor( expression_type const& /*expr*/,
                Geo_t const& /*geom*/ )
            :
            M_one( nComponents, vector_comp )
        {
            //                 std::cout << "one = " << M_one << "\n"
            //                           << "M=" << shape::M << "\n"
            //                           << "N=" << shape::N << "\n";
        }
        template<typename IM>
        void init( IM const& /*im*/ )
        {
        }
        void update( Geo_t const& /*geom*/, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
        }
        void update( Geo_t const& /*geom*/, Basis_i_t const& /*fev*/ )
        {
        }
        void update( Geo_t const& /*geom*/ )
        {
        }
        void update( Geo_t const& /*geom*/, uint16_type /*face*/ )
        {
        }

        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return ( gmc_type::nDim>=c1 )&&( ( c1==(uint16_type)CType ) || ( CType==-1 ) );
            //return M_one[c1];
        }
        template<int PatternContext>
        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/,
                 mpl::int_<PatternContext> ) const
        {
            return ( gmc_type::nDim>=c1 )&&( ( c1==(uint16_type)CType ) || ( CType==-1 ) );
            //return M_one[c1];
        }

        FEELPP_STRONG_INLINE value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return ( gmc_type::nDim>=c1 )&&( ( c1==(uint16_type)CType ) || ( CType==-1 ) );
            //return M_one[c1];
        }
        FEELPP_STRONG_INLINE value_type
        evalq( uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return ( gmc_type::nDim>=c1 )&&( ( c1==(uint16_type)CType ) || ( CType==-1 ) );
            //return M_one[c1];
        }
        vector_type M_one;
    };

};

inline
Expr<One<-1> >
one()
{
    return Expr< One<-1> >(  One<-1>() );
}

inline
Expr<One<0> >
oneX()
{
    return Expr< One<0> >(  One<0>() );
}

inline
Expr<One<1> >
oneY()
{
    return Expr< One<1> >(  One<1>() );
}

inline
Expr<One<2> >
oneZ()
{
    return Expr< One<2> >(  One<2>() );
}
inline
Expr<One<0> >
unitX()
{
    return Expr< One<0> >(  One<0>() );
}

inline
Expr<One<1> >
unitY()
{
    return Expr< One<1> >(  One<1>() );
}

inline
Expr<One<2> >
unitZ()
{
    return Expr< One<2> >(  One<2>() );
}

} // vf
} // Feel
#endif /* FEELPP_VF_ONE_HPP */
