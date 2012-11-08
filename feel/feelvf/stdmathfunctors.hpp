/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-06-07

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2012 Universite Joseph Fourier (Grenoble I)

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
   \file functions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-06-07
 */
#if !defined( STD_MATH_UNARY_FUNCTORS_HPP )
#define STD_MATH_UNARY_FUNCTORS_HPP 1

# include <boost/preprocessor/comparison/less.hpp>
# include <boost/preprocessor/logical/and.hpp>
# include <boost/preprocessor/control/if.hpp>
# include <boost/preprocessor/list/at.hpp>
# include <boost/preprocessor/list/cat.hpp>
# include <boost/preprocessor/list/for_each_product.hpp>
# include <boost/preprocessor/logical/or.hpp>
# include <boost/preprocessor/tuple/to_list.hpp>
# include <boost/preprocessor/tuple/eat.hpp>
# include <boost/preprocessor/facilities/empty.hpp>
# include <boost/preprocessor/punctuation/comma.hpp>
# include <boost/preprocessor/facilities/identity.hpp>

#include <boost/utility/enable_if.hpp>

/// \cond detail
#include <feel/feelcore/traits.hpp>
#include <feel/feelvf/unaryfunctor.hpp>

namespace Feel
{
namespace vf
{
namespace details
{
/**
   sign function
 */
template<typename T>
T
sign( T const& x )
{
    return ( x > T( 0.0 ) )?T( 1.0 ):( ( x < T( 0.0 ) )?T( -1.0 ):0.0 );
}

} // details

}
}

# /* Information about functions  */
#
# /* Accessors for the operator datatype. */
# define VF_FUNC_SYMBOL(O)        BOOST_PP_TUPLE_ELEM(8, 0, O)
# define VF_FUNC_NAME(O)          BOOST_PP_TUPLE_ELEM(8, 1, O)
# define VF_FUNC_IMPL(O)          BOOST_PP_TUPLE_ELEM(8, 2, O)
# define VF_FUNC_NAME_STRING(O)   BOOST_PP_TUPLE_ELEM(8, 3, O)
# define VF_FUNC_DOMAIN(O)        BOOST_PP_TUPLE_ELEM(8, 4, O)
# define VF_FUNC_IS_FLOATING(O )  BOOST_PP_TUPLE_ELEM(8, 5, O)
# define VF_FUNC_IS_LOGICAL(O)    BOOST_PP_TUPLE_ELEM(8, 6, O)
# define VF_FUNC_IS_POLYNOMIAL(O) BOOST_PP_TUPLE_ELEM(8, 7, O)
#
# /* List of applicative unary functions. */
# define VF_APPLICATIVE_UNARY_FUNCS \
   BOOST_PP_TUPLE_TO_LIST( \
      15, \
      ( \
         ( abs  , __Abs__ , Feel::math::abs     ,"absolute value"      , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( cos  , __Cos__ , Feel::math::cos     ,"cosine"              , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( sin  , __Sin__ , Feel::math::sin     ,"sine"                , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( tan  , __Tan__ , Feel::math::tan     ,"tangent"             , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( acos , __ACos__, Feel::math::acos    ,"inverse cosine"      , BoundedDomain<value_type>(-1.0,1.0), 1, 0, 2), \
         ( asin , __ASin__, Feel::math::asin    ,"inverse sine"        , BoundedDomain<value_type>(-1.0,1.0), 1, 0, 2), \
         ( atan , __ATan__, Feel::math::atan    ,"inverse tangent"     , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( cosh , __Cosh__, Feel::math::cosh    ,"hyperbolic cosine"   , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( sinh , __Sinh__, Feel::math::sinh    ,"hyperbolic sine"     , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( tanh , __Tanh__, Feel::math::tanh    ,"hyperbolic tangent"  , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( exp  , __Exp__ , Feel::math::exp     ,"exponential"         , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( log  , __Log__ , Feel::math::log     ,"logarithm"           , PositiveDomain<value_type>()       , 1, 0, 2), \
         ( sqrt , __Sqrt__, Feel::math::sqrt    ,"square root"         , PositiveDomain<value_type>()       , 1, 0, 2), \
         ( sign , __Sign__, details::sign       ,"sign"                , UnboundedDomain<value_type>()      , 1, 0, 0), \
         ( chi  , __Chi__ ,                     ,"chi"                 , UnboundedDomain<value_type>()      , 0, 1, 0) \
      ) \
   ) \
   /**/
#
#
# /* Generates code for all binary operators and integral type pairs. */
# define VF_UNARY_FUNCTIONS(_, O) \
      VF_UNARY_FUNCTIONS_CODE O \
   /**/

#if defined( FEELPP_HAS_QD_H ) && defined(FEELPP_HAS_MPFR)
# define VF_CHECK_ARITHMETIC_TYPE()                                        \
   BOOST_STATIC_ASSERT( (::boost::is_arithmetic<value_1_type>::value ||    \
                         ::boost::is_same<value_1_type, std::complex<float> >::value || \
                         ::boost::is_same<value_1_type, std::complex<double> >::value || \
                         ::boost::is_same<value_1_type,mp_type>::value ||  \
                         ::boost::is_same<value_1_type,dd_real>::value ||  \
                         ::boost::is_same<value_1_type,qd_real>::value) ); \
   /**/
#elif defined( FEELPP_HAS_QD_H )
# define VF_CHECK_ARITHMETIC_TYPE()                                        \
   BOOST_STATIC_ASSERT( (::boost::is_arithmetic<value_1_type>::value ||    \
                         ::boost::is_same<value_1_type, std::complex<float> >::value || \
                         ::boost::is_same<value_1_type, std::complex<double> >::value || \
                         ::boost::is_same<value_1_type,dd_real>::value ||  \
                         ::boost::is_same<value_1_type,qd_real>::value) ); \
   /**/
#elif defined( FEELPP_HAS_MPFR )
# define VF_CHECK_ARITHMETIC_TYPE()                                        \
   BOOST_STATIC_ASSERT( (::boost::is_arithmetic<value_1_type>::value ||    \
                         ::boost::is_same<value_1_type, std::complex<float> >::value || \
                         ::boost::is_same<value_1_type, std::complex<double> >::value || \
                         ::boost::is_same<value_1_type,mp_type>::value) ); \
   /**/
#else
# define VF_CHECK_ARITHMETIC_TYPE()                                     \
    BOOST_STATIC_ASSERT( ( ::boost::is_arithmetic<value_1_type>::value || \
                           ::boost::is_same<value_1_type, std::complex<float> >::value || \
                           ::boost::is_same<value_1_type, std::complex<double> >::value ) \
                         );                                             \
   /**/
#endif

# define VF_IM_IS_POLY(O)                                               \
    BOOST_PP_IF( BOOST_PP_EQUAL( VF_FUNC_IS_POLYNOMIAL(O) , 0) ,        \
                 0 ,                                                    \
                 1   )                                                  \
    /**/




#define VF_UNARY_FUNCTIONS_CODE(O)                                      \
    template < typename ExprT1 >                                        \
class VF_FUNC_NAME( O ) : public UnaryFunctor<typename ExprT1::value_type>              \
    {                                                                   \
    public:                                                             \
                                                                        \
        static const size_type context = ExprT1::context;               \
        static const bool is_terminal = false;                          \
                                                                        \
        static const uint16_type imorder = VF_FUNC_IS_POLYNOMIAL(O)*ExprT1::imorder;             \
        static const bool imIsPoly = (VF_IM_IS_POLY(O) || (ExprT1::imorder==0)); \
                                                                        \
        template<typename Func>                                         \
            struct HasTestFunction                                      \
        {                                                               \
            static const bool result = false;                           \
        };                                                              \
                                                                        \
        template<typename Func>                                         \
            struct HasTrialFunction                                     \
        {                                                               \
            static const bool result = false;                           \
        };                                                              \
                                                                        \
        typedef UnaryFunctor<typename ExprT1::value_type> super;        \
        typedef typename super::functordomain_type functordomain_type;  \
        typedef typename super::functordomain_ptrtype functordomain_ptrtype; \
        typedef ExprT1 expression_1_type;                               \
        typedef VF_FUNC_NAME(O)<ExprT1> this_type;                      \
        typedef typename expression_1_type::value_type value_1_type;    \
        typedef value_1_type value_type;                                \
                                                                        \
        VF_CHECK_ARITHMETIC_TYPE()                                      \
                                                                        \
            explicit VF_FUNC_NAME(O)( expression_1_type const& __expr1  ) \
            :                                                           \
            super( VF_FUNC_NAME_STRING(O), functordomain_ptrtype(new VF_FUNC_DOMAIN(O) )), \
            _M_expr_1( __expr1 )                                        \
            {                                                           \
                Debug( 5051 ) << "VF_FUNC_NAME(O)::VF_FUNC_NAME(O) default constructor\n"; \
            }                                                           \
                                                                        \
        VF_FUNC_NAME(O)( VF_FUNC_NAME(O) const& __vfp  )                \
            :                                                           \
            super( VF_FUNC_NAME_STRING(O), functordomain_ptrtype(new VF_FUNC_DOMAIN(O) )), \
            _M_expr_1( __vfp._M_expr_1 )                                \
                {                                                       \
                    Debug( 5051 ) << "VF_FUNC_NAME(O)::VF_FUNC_NAME(O) copy constructor\n"; \
                }                                                       \
                                                                        \
        bool isSymetric() const { return false; }                       \
                                                                        \
        void eval( int nx, value_type const* x, value_type* f ) const   \
        {                                                               \
            for( int i = 0; i < nx; ++i )                               \
                f[i] = VF_FUNC_IMPL(O)( x[i] );                         \
        }                                                               \
                                                                        \
        expression_1_type const& expression() const { return _M_expr_1; } \
                                                                        \
        template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t> \
            struct tensor                                               \
        {                                                               \
            typedef this_type expression_type;                          \
            typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t,Basis_j_t> tensor2_expr_type; \
            typedef typename tensor2_expr_type::value_type value_type;  \
            typedef typename vf::detail::ExtractGm<Geo_t>::gmc_ptrtype gmc_ptrtype; \
            typedef typename vf::detail::ExtractGm<Geo_t>::gmc_type gmc_type; \
            typedef typename tensor2_expr_type::shape shape;            \
                                                                        \
            struct is_zero { static const bool value = tensor2_expr_type::is_zero::value; }; \
                                                                        \
            tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) \
                :                                                       \
                _M_expr( expr.expression(), geom ),                     \
                _M_gmc(vf::detail::ExtractGm<Geo_t>::get( geom ) )         \
                    {                                                   \
                        update( geom );                                 \
                    }                                                   \
            tensor( this_type const& expr,Geo_t const& geom, Basis_i_t const& /*fev*/ ) \
                :                                                       \
                _M_expr( expr.expression(), geom ),                     \
                _M_gmc(vf::detail::ExtractGm<Geo_t>::get( geom ) )         \
                    {                                                   \
                        update( geom );                                 \
                    }                                                   \
            tensor( this_type const& expr, Geo_t const& geom )             \
                :                                                       \
                _M_expr( expr.expression(), geom ),                     \
                _M_gmc(vf::detail::ExtractGm<Geo_t>::get( geom ) )         \
                    {                                                   \
                        update( geom );                                 \
                    }                                                   \
            template<typename IM>                                       \
                void init( IM const& im )                               \
            {                                                           \
                _M_expr.init( im );                                     \
            }                                                           \
            void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) \
            {                                                           \
                update( geom );                                         \
            }                                                           \
            void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )  \
            {                                                           \
                update( geom );                                         \
            }                                                           \
            void update( Geo_t const& geom )                            \
            {                                                           \
                _M_expr.update( geom );                                 \
            }                                                           \
            void update( Geo_t const& geom, uint16_type face )          \
            {                                                           \
                _M_expr.update( geom, face );                           \
            }                                                           \
                              \
                value_type                                              \
                evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
            template<int PatternContext> \
                value_type                                              \
                evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q, \
                         mpl::int_<PatternContext> ) const              \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
                                               \
                value_type                                              \
                evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
            value_type                                                  \
                evalq( uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return VF_FUNC_IMPL(O)( _M_expr.evalq( c1, c2, q ) );   \
            }                                                           \
        private:                                                        \
            tensor2_expr_type _M_expr;                                  \
            gmc_ptrtype _M_gmc;                                         \
        };                                                              \
                                                                        \
    protected:                                                          \
        VF_FUNC_NAME(O)() {}                                            \
                                                                        \
        expression_1_type _M_expr_1;                                    \
    };                                                                  \
                                                                        \
    template<typename ExprT1>                                           \
    inline                                                              \
    Expr< VF_FUNC_NAME( O )<typename mpl::if_<boost::is_arithmetic<ExprT1>, \
                                              mpl::identity<Cst<ExprT1> >, \
                                              mpl::identity<Expr<ExprT1> > >::type::type > > \
    VF_FUNC_SYMBOL( O )( Expr<ExprT1> const& __e1 )                     \
    {                                                                   \
        typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,         \
            mpl::identity<Cst<ExprT1> >,                                \
            mpl::identity<Expr<ExprT1> > >::type::type t1;               \
        typedef VF_FUNC_NAME(O)<t1> expr_t;                             \
        return Expr< expr_t >(  expr_t( t1( __e1 ) ) );                 \
    }                                                                   \
    template<typename ExprT1>                                           \
    inline                                                              \
    Expr< VF_FUNC_NAME( O )<Cst<ExprT1> > >                             \
    VF_FUNC_SYMBOL( O )( ExprT1 const& __e1, typename boost::enable_if<boost::is_arithmetic<ExprT1> >::type* dummy = 0 ) \
    {                                                                   \
        typedef Cst<ExprT1> t1;                                         \
        typedef VF_FUNC_NAME(O)<t1> expr_t;                             \
        return Expr< expr_t >(  expr_t( t1( __e1 ) ) );                 \
    }                                                                   \
    /**/
#

namespace Feel
{
namespace vf
{

BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_UNARY_FUNCTIONS, 1, ( VF_APPLICATIVE_UNARY_FUNCS ) )
}
}
/// \endcond
#endif /* STD_MATH_UNARY_FUNCTORS_HPP */
