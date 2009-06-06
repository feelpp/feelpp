/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-06-07

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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

/// \cond detail
#include <life/lifecore/traits.hpp>
#include <life/lifevf/unaryfunctor.hpp>

namespace Life
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
# define VF_FUNC_SYMBOL(O)      BOOST_PP_TUPLE_ELEM(7, 0, O)
# define VF_FUNC_NAME(O)        BOOST_PP_TUPLE_ELEM(7, 1, O)
# define VF_FUNC_IMPL(O)        BOOST_PP_TUPLE_ELEM(7, 2, O)
# define VF_FUNC_NAME_STRING(O) BOOST_PP_TUPLE_ELEM(7, 3, O)
# define VF_FUNC_DOMAIN(O)      BOOST_PP_TUPLE_ELEM(7, 4, O)
# define VF_FUNC_IS_FLOATING(O) BOOST_PP_TUPLE_ELEM(7, 5, O)
# define VF_FUNC_IS_LOGICAL(O)  BOOST_PP_TUPLE_ELEM(7, 6, O)
#
# /* List of applicative unary functions. */
# define VF_APPLICATIVE_UNARY_FUNCS \
   BOOST_PP_TUPLE_TO_LIST( \
      15, \
      ( \
         ( abs  , __Abs__ , Life::math::abs     ,"absolute value"      , UnboundedDomain<value_type>()      , 1, 0), \
         ( cos  , __Cos__ , Life::math::cos     ,"cosine"              , UnboundedDomain<value_type>()      , 1, 0), \
         ( sin  , __Sin__ , Life::math::sin     ,"sine"                , UnboundedDomain<value_type>()      , 1, 0), \
         ( tan  , __Tan__ , Life::math::tan     ,"tangent"             , UnboundedDomain<value_type>()      , 1, 0), \
         ( acos , __ACos__, Life::math::acos    ,"inverse cosine"      , BoundedDomain<value_type>(-1.0,1.0), 1, 0), \
         ( asin , __ASin__, Life::math::asin    ,"inverse sine"        , BoundedDomain<value_type>(-1.0,1.0), 1, 0), \
         ( atan , __ATan__, Life::math::atan    ,"inverse tangent"     , UnboundedDomain<value_type>()      , 1, 0), \
         ( cosh , __Cosh__, Life::math::cosh    ,"hyperbolic cosine"   , UnboundedDomain<value_type>()      , 1, 0), \
         ( sinh , __Sinh__, Life::math::sinh    ,"hyperbolic sine"     , UnboundedDomain<value_type>()      , 1, 0), \
         ( tanh , __Tanh__, Life::math::tanh    ,"hyperbolic tangent"  , UnboundedDomain<value_type>()      , 1, 0), \
         ( exp  , __Exp__ , Life::math::exp     ,"exponential"         , UnboundedDomain<value_type>()      , 1, 0), \
         ( log  , __Log__ , Life::math::log     ,"logarithm"           , PositiveDomain<value_type>()       , 1, 0), \
         ( sqrt , __Sqrt__, Life::math::sqrt    ,"square root"         , PositiveDomain<value_type>()       , 1, 0), \
         ( sign , __Sign__, details::sign  ,"sign"                , UnboundedDomain<value_type>()      , 1, 0), \
         ( chi  , __Chi__ ,                ,"chi"                 , UnboundedDomain<value_type>()      , 0, 1) \
      ) \
   ) \
   /**/
#
#
# /* Generates code for all binary operators and integral type pairs. */
# define VF_UNARY_FUNCTIONS(_, O) \
      VF_UNARY_FUNCTIONS_CODE O \
   /**/

#if defined( HAVE_QD_H ) && defined(HAVE_MPFR)
# define VF_CHECK_ARITHMETIC_TYPE()                                        \
   BOOST_STATIC_ASSERT( (::boost::is_arithmetic<value_1_type>::value ||    \
                         ::boost::is_same<value_1_type, std::complex<float> >::value || \
                         ::boost::is_same<value_1_type, std::complex<double> >::value || \
                         ::boost::is_same<value_1_type,mp_type>::value ||  \
                         ::boost::is_same<value_1_type,dd_real>::value ||  \
                         ::boost::is_same<value_1_type,qd_real>::value) ); \
   /**/
#elif defined( HAVE_QD_H )
# define VF_CHECK_ARITHMETIC_TYPE()                                        \
   BOOST_STATIC_ASSERT( (::boost::is_arithmetic<value_1_type>::value ||    \
                         ::boost::is_same<value_1_type, std::complex<float> >::value || \
                         ::boost::is_same<value_1_type, std::complex<double> >::value || \
                         ::boost::is_same<value_1_type,dd_real>::value ||  \
                         ::boost::is_same<value_1_type,qd_real>::value) ); \
   /**/
#elif defined( HAVE_MPFR )
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

#define VF_UNARY_FUNCTIONS_CODE(O)                                      \
    template < typename ExprT1 >                                        \
class VF_FUNC_NAME( O ) : public UnaryFunctor<typename ExprT1::value_type>              \
    {                                                                   \
    public:                                                             \
                                                                        \
        static const size_type context = ExprT1::context;               \
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
            typedef typename detail::ExtractGm<Geo_t>::gmc_ptrtype gmc_ptrtype; \
            typedef typename detail::ExtractGm<Geo_t>::gmc_type gmc_type; \
            typedef typename tensor2_expr_type::shape shape;            \
                                                                        \
            struct is_zero { static const bool value = tensor2_expr_type::is_zero::value; }; \
                                                                        \
            tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) \
                :                                                       \
                _M_expr( expr.expression(), geom ),                     \
                _M_gmc( detail::ExtractGm<Geo_t>::get( geom ) )         \
                    {                                                   \
                        update( geom );                                 \
                    }                                                   \
            tensor( this_type const& expr,Geo_t const& geom, Basis_i_t const& /*fev*/ ) \
                :                                                       \
                _M_expr( expr.expression(), geom ),                     \
                _M_gmc( detail::ExtractGm<Geo_t>::get( geom ) )         \
                    {                                                   \
                        update( geom );                                 \
                    }                                                   \
            tensor( this_type const& expr, Geo_t const& geom )             \
                :                                                       \
                _M_expr( expr.expression(), geom ),                     \
                _M_gmc( detail::ExtractGm<Geo_t>::get( geom ) )         \
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
            template<typename IndexI, typename IndexJ>                  \
                value_type                                              \
                evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
            template<typename IndexI, typename IndexJ, int PatternContext> \
                value_type                                              \
                evalijq( IndexI const& /*i*/, IndexJ const& /*j*/, uint16_type c1, uint16_type c2, uint16_type q, \
                         mpl::int_<PatternContext> ) const              \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
            template<typename IndexI>                                   \
                value_type                                              \
                evaliq( IndexI const& /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
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
                                              mpl::identity<ExprT1> >::type::type > > \
    VF_FUNC_SYMBOL( O )( ExprT1 const& __e1 )                           \
    {                                                                   \
        typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,         \
            mpl::identity<Cst<ExprT1> >,                                \
            mpl::identity<ExprT1> >::type::type t1;                     \
        typedef VF_FUNC_NAME(O)<t1> expr_t;                             \
        return Expr< expr_t >(  expr_t( t1( __e1 ) ) );                 \
    }                                                                   \
    /**/
#

namespace Life
{
namespace vf
{

BOOST_PP_LIST_FOR_EACH_PRODUCT(VF_UNARY_FUNCTIONS, 1, (VF_APPLICATIVE_UNARY_FUNCS))
}
}
/// \endcond
#endif /* STD_MATH_UNARY_FUNCTORS_HPP */
