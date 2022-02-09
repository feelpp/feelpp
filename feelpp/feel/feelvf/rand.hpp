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
   \file rand.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_RAND_HPP
#define FEELPP_VF_RAND_HPP 1

#include <random>
#include <boost/math/special_functions/round.hpp>
#include <feel/feelvf/expr.hpp>

namespace Feel
{
namespace vf
{
namespace vfdetails
{
template<typename T>
class Rand_d
{
public:
    
    using value_type = T;
    using distribution_t = typename mpl::if_<std::is_integral<T>,
                                             mpl::identity<std::uniform_int_distribution<T>>,
                                             mpl::identity<std::uniform_real_distribution<T>>>::type::type;
                                             
    Rand_d () : Rand_d( T(0), T(1) ) {}
    Rand_d( value_type lo, value_type hi )
        :
        r( std::bind( distribution_t(lo,hi), std::default_random_engine() ) )
        {}
    value_type operator()() const { return r(); }
private:
    std::function<value_type()> r;
    
};

}
class RandBase {};

template < class T>
class Rand : public RandBase
{
public:

    //BOOST_STATIC_ASSERT( ::boost::is_arithmetic<T>::value );

    static const size_type context = vm::JACOBIAN;
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
            mpl::identity<mpl::identity<T> > >::type::type::type value_type;
    using evaluate_type = Eigen::Matrix<value_type,1,1>;

    typedef Rand<T> expression_type;

    constexpr explicit Rand()  
        : 
        M_low( value_type(0) ), 
        M_high( value_type(1) ), 
        M_r( M_low, M_high ) 
        {}
    Rand( value_type lo, value_type hi )  : M_low(lo), M_high(hi), M_r( lo, hi )  {}
    Rand( Rand && c ) = default;
    Rand( Rand const& c ) = default;
    Rand& operator=( Rand const& c ) = default;

    //! polynomial order
    uint16_type polynomialOrder() const { return 0; }

    //! expression is polynomial?
    constexpr bool isPolynomial() const { return true; }

    value_type low() const { return M_low; }
    value_type high() const { return M_high; }
    constexpr value_type value() const
    {
        return M_r();
    }

    constexpr value_type evaluate() const
    {
        return M_r();
    }

    constexpr value_type evaluate( bool ) const
    {
        return M_r();
    }

    constexpr evaluate_type evaluate( bool, WorldComm const& ) const
    {
        return evaluate_type::Constant( M_r() );
    }

    template<typename... TheExpr>
    struct Lambda
    {
        typedef expression_type type;
    };
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  ) { return typename Lambda<TheExpr...>::type(M_r()); }

    template<typename Geo_t, typename Basis_i_t=mpl::void_, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename Rand<T>::expression_type expression_type;
        typedef typename Rand<T>::value_type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >, mpl::identity<vf::detail::gmc<0> >, mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
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
            M_r(expr.low(), expr.high())
        {
        }
        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/ )
            :
            M_r(expr.low(),expr.high())
        {
        }
        tensor( expression_type const& expr, Geo_t const& /*geom*/ )
            :
            M_r(expr.low(), expr.high())
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
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
        {
        }

        constexpr value_type
        evalij( uint16_type /*i*/, uint16_type /*j*/ ) const
        {
            return M_r();
        }


        constexpr value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return M_r();
        }
        template<int PatternContext>
        constexpr value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/,
                 mpl::int_<PatternContext> ) const
        {
            return M_r();
        }

        constexpr value_type
        evaliq( uint16_type /*i*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return M_r();
        }
        constexpr value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return M_r();
        }
        vfdetails::Rand_d<value_type> M_r;
    };

private:
    value_type M_low;
    value_type M_high;
    vfdetails::Rand_d<value_type> M_r;
};


/**
 * provide random expression in ]lo,hi[
 * @ingroup DSEL-Variational-Formulation
 * @param lo lower end of the interval
 * @param hi upper of the interval
 */
template<typename T>
inline
Expr< Rand<T> >
rand()
{
    typedef Rand<T> rand_t;
    return Expr< rand_t >(  rand_t( T(0), T(1) ) );
}

/**
 * provide random expression in ]lo,hi[
 * @ingroup DSEL-Variational-Formulation
 * @param lo lower end of the interval
 * @param hi upper of the interval
 */
template<typename T>
inline
Expr< Rand<T> >
rand( T lo, T hi )
{
    typedef Rand<T> rand_t;
    return Expr< rand_t >(  rand_t( lo, hi ) );
}

} // vf
} // Feel

#endif /* FEELPP_VF_RAND_HPP */
