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
   \file cst.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_CST_HPP
#define FEELPP_VF_CST_HPP 1

#include <feel/feelvf/expr.hpp>

namespace Feel
{
namespace vf
{
class CstBase {};

template < class T>
class Cst : public CstBase
{
public:

    //BOOST_STATIC_ASSERT( ::boost::is_arithmetic<T>::value );

    static const size_type context = 0;
    static const bool is_terminal = false;

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

    typedef typename mpl::if_<boost::is_reference_wrapper<T>,
                              mpl::identity<T>,
                              mpl::identity<mpl::identity<T> > >::type::type::type _value_type;
    using value_type = std::decay_t<_value_type>;
    using evaluate_type = Eigen::Matrix<value_type,1,1>;

    typedef Cst<T> expression_type;

    Cst() : M_constant( 0 )
    {
        //DVLOG(2) << "Cst::Cst( default ) : constant value: " << M_constant << "\n";
    }
    constexpr explicit Cst( const T& value )
        :
        M_constant( value )
    {
    }
    constexpr explicit Cst( T&& value )
        :
        M_constant( std::forward<T>(value) )
        {
        }

    Cst( Cst && c ) = default;
    Cst( Cst const& c ) = default;
    Cst& operator=( Cst const& c ) = default;

    //! polynomial order
    constexpr uint16_type polynomialOrder() const { return 0; }

    //! expression is polynomial?
    constexpr bool isPolynomial() const { return true; }

    T& value_ref() 
        {
            return M_constant;
        }
    T const& value_ref()  const
        {
            return M_constant;
        }
    constexpr value_type value() const
    {
        return M_constant;
    }
#if 0
    constexpr value_type evaluate() const
    {
        return M_constant;
    }

    constexpr value_type evaluate( bool ) const
    {
        return M_constant;
    }
#endif
    constexpr evaluate_type evaluate( bool, worldcomm_ptr_t const& ) const
    {
        return evaluate_type::Constant( M_constant );
    }

    template<typename... TheExpr>
    struct Lambda
    {
        typedef expression_type type;
    };
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  ) { return typename Lambda<TheExpr...>::type(M_constant); }

    // void setParameterValues( std::map<std::string,value_type> const& mp ) {}

    template <typename SymbolsExprType>
    expression_type applySymbolsExpr( SymbolsExprType const& se ) const
        {
            return *this;
        }

    // template <typename TheSymbolExprType>
    // bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const { return false; }
    // template <typename TheSymbolExprType>
    // void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const {}

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        return expression_type(0);
    }




    template<typename Geo_t, typename Basis_i_t=mpl::void_, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename Cst<T>::expression_type expression_type;
        typedef typename Cst<T>::value_type value_type;

        using key_type = key_t<Geo_t>;
        using gmc_ptrtype = gmc_ptr_t<Geo_t>;
        using gmc_type = gmc_t<Geo_t>;
        typedef Shape<gmc_type::nDim, Scalar, false, false> shape;


        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
            :
            M_constant( expr.value_ref() )
        {
        }
        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/ )
            :
            M_constant( expr.value_ref() )
        {
        }
        tensor( expression_type const& expr, Geo_t const& /*geom*/ )
            :
            M_constant( expr.value_ref() )
        {
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
        void update( Geo_t const&, Basis_i_t const& , Basis_j_t const&  )
        {
        }
        void update( Geo_t const& , Basis_i_t const&  )
        {
        }
        void update( Geo_t const& )
        {
        }
        void update( Geo_t const&, uint16_type )
        {
        }
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
        {
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {}


        constexpr value_type
        evalij( uint16_type /*i*/, uint16_type /*j*/ ) const
        {
            return M_constant;
        }


        constexpr value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return M_constant;
        }
        template<int PatternContext>
        constexpr value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/,
                 mpl::int_<PatternContext> ) const
        {
            return M_constant;
        }

        constexpr value_type
        evaliq( uint16_type /*i*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return M_constant;
        }
        constexpr value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return M_constant;
        }
        T M_constant;
    };

protected:

    const T M_constant;
};

template<typename T>
inline
Expr< Cst<T> >
constant( T v )
{
    typedef Cst<T> cst_t;
    return Expr< cst_t >(  cst_t( v ) );
}

template<typename T>
inline
Expr< Cst<T> >
cst( T v )
{
    typedef Cst<T> cst_t;
    return Expr< cst_t >(  cst_t( v ) );
}
template<typename T>
inline
Expr< Cst<boost::reference_wrapper<T> > >
cst_ref( T& v )
{
    typedef Cst<boost::reference_wrapper<T> > cst_t;
    return Expr< cst_t >(  cst_t( boost::ref( v ) ) );
}
template<typename T>
inline
Expr< Cst<boost::reference_wrapper<T> > >
constant_ref( T& v )
{
    return cst_ref( v );
}

} // vf
} // Feel

#endif /* FEELPP_VF_CST_HPP */
