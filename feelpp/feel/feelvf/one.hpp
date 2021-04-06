/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-13

  Copyright (C) 2014-2016 Feel++ Consortium

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

#include <feel/feeldiscr/enums.hpp>

namespace Feel
{
namespace vf
{

/**
 * CType = -3 : dynamic case, the one compoent(s) is given in argument
 * CType = -2 : (0 0 0)
 * CType = -1 : (1 1 1)
 * CType =  0 : (1 0 0)
 * CType =  1 : (0 1 0)
 * CType =  2 : (0 0 1)
 */
template<int CType>
class One
{
public:
    static const size_type context = 0;
    static const bool is_terminal = true;

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
    template<typename Func>
    static const bool has_test_basis = false;
    template<typename Func>
    static const bool has_trial_basis = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;



    typedef One<CType> this_type;

    typedef double value_type;

    One() = default;
    One( One const& /*__vff*/ ) = default;
    One( ComponentType c )
        :
        M_oneComponentsDynamic( { c } )
        {}
    One( std::set<ComponentType> const& c)
        :
        M_oneComponentsDynamic( c )
        {}

    std::set<ComponentType> const& oneComponentsDynamic() const { return M_oneComponentsDynamic; }

    template<typename... TheExpr>
    struct Lambda
    {
        typedef this_type type;
    };
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  ) { return this_type(); }

    //! polynomial order
    constexpr uint16_type polynomialOrder() const { return 0; }

    //! expression is polynomial?
    constexpr bool isPolynomial() const { return true; }

    template <typename SymbolsExprType>
    this_type applySymbolsExpr( SymbolsExprType const& se ) const
        {
            return *this;
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        return One<-2>();
    }


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        using key_type = key_t<Geo_t>;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef Shape<gmc_type::nDim, Vectorial, false, false> shape;
        static const bool theshape = ( shape::M == gmc_type::nDim && shape::N == 1 );
        BOOST_MPL_ASSERT_MSG( theshape,
                              INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_1,
                              ( mpl::int_<shape::M>, mpl::int_<shape::N> ) );

        typedef typename expression_type::value_type value_type;

        static const uint16_type nComponents = gmc_type::nDim;
        //static const int16_type vector_comp = ( CType==-1 )?1:CType;

        using vector_type = Eigen::Matrix<value_type,nComponents,1>;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/,
                Basis_j_t const& /*feu*/ )
            :
            tensor( expr, geom )
            {}

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/ )
            :
            tensor( expr, geom )
            {}
        tensor( expression_type const& expr,
                Geo_t const& /*geom*/ )
            :
            M_one( vector_type::Zero() )
            {
                if constexpr ( CType == -1 )
                    M_one = vector_type::Ones();
                else if constexpr ( CType >= 0 && CType < gmc_type::nDim )
                    M_one(CType) = 1;
                else if constexpr ( CType == -3 )
                {
                    auto const& theDynComp = expr.oneComponentsDynamic();
                    if ( nComponents > 0 && theDynComp.find( ComponentType::X ) != theDynComp.end() )
                        M_one( 0 ) = 1;
                    if ( nComponents > 1 && theDynComp.find( ComponentType::Y ) != theDynComp.end() )
                        M_one( 1 ) = 1;
                    if ( nComponents > 2 && theDynComp.find( ComponentType::Z ) != theDynComp.end() )
                        M_one( 2 ) = 1;
                }
            }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                expression_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            tensor( expr, geom, theInitArgs... )
            {}
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
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {}

        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
        {
        }

        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evalq( c1, c2, q );
        }
        template<int PatternContext>
        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return this->evalq( c1, c2, q );
        }

        FEELPP_STRONG_INLINE value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evalq( c1, c2, q );
        }
        FEELPP_STRONG_INLINE value_type
        evalq( uint16_type c1, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            if ( c1 >= gmc_type::nDim )
                return value_type(0);
            else
                return M_one( c1 );
        }
        FEELPP_STRONG_INLINE Eigen::Map<const vector_type>
        evalijq( uint16_type i, uint16_type j, uint16_type /*q*/ ) const
        {
            return Eigen::Map<const vector_type>(M_one.data());
        }
        FEELPP_STRONG_INLINE Eigen::Map<const vector_type>
        evaliq( uint16_type i, uint16_type /*q*/ ) const
        {
            return Eigen::Map<const vector_type>(M_one.data());
        }
        FEELPP_STRONG_INLINE Eigen::Map<const vector_type>
        evalq( uint16_type /*q*/ ) const
        {
            return Eigen::Map<const vector_type>(M_one.data());
        }
    private:
        vector_type M_one;
    };
private :
    std::set<ComponentType> M_oneComponentsDynamic;
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

inline
Expr<One<-2> >
vector_zero()
{
    return Expr< One<-2> >(  One<-2>() );
}

inline
Expr<One<-3> >
one( ComponentType c)
{
    return Expr< One<-3> >(  One<-3>( c ) );
}
inline
Expr<One<-3> >
one( std::set<ComponentType> const& c)
{
    return Expr< One<-3> >(  One<-3>( c ) );
}
inline
Expr<One<-3> >
one( std::initializer_list<ComponentType> const& c)
{
    return Expr< One<-3> >(  One<-3>( std::set<ComponentType>( c ) ) );
}

} // vf
} // Feel
#endif /* FEELPP_VF_ONE_HPP */
