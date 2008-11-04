/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-01-17

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)

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
   \file operators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-01-17
 */
#if !defined( __LIFE_VF_OPERATORS_HPP )
#define __LIFE_VF_OPERATORS_HPP 1

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
# include <boost/preprocessor/stringize.hpp>
namespace Life
{
namespace vf
{
/// \cond detail
# /* Accessors for the operator datatype. */
# define VF_OPERATOR_NAME(O)           BOOST_PP_TUPLE_ELEM(9, 0, O)
# define VF_OPERATOR_SYMBOL(O)         BOOST_PP_TUPLE_ELEM(9, 1, O)
# define VF_OPERATOR_TERM(O)           BOOST_PP_TUPLE_ELEM(9, 2, O)
# // 0 means any rank, 1 or more is the rank for which the operator is defined
# define VF_OPERATOR_TYPE_RANK_DEF(O)  BOOST_PP_TUPLE_ELEM(9, 3, O)
# define VF_OPERATOR_TYPE_IS_COMP(O)   BOOST_PP_TUPLE_ELEM(9, 4, O)
# define VF_OPERATOR_TYPE_COMP(O)      BOOST_PP_TUPLE_ELEM(9, 5, O)
# define VF_OPERATOR_CONTEXT(O)        BOOST_PP_TUPLE_ELEM(9, 6, O)
# define VF_OPERATOR_RANK(O)           BOOST_PP_TUPLE_ELEM(9, 7, O)
# define VF_OPERATOR_TRANSPOSE(O)      BOOST_PP_TUPLE_ELEM(9, 8, O)

const size_type jk = vm::JACOBIAN|vm::KB;
const size_type jkd = jkd|vm::FIRST_DERIVATIVE;
const size_type jkn = jkd|vm::FIRST_DERIVATIVE_NORMAL;
const size_type jkg = jkd|vm::GRAD;
const size_type jkh = jkd|vm::HESSIAN;


# /* List of applicative operators. */
# define VF_OPERATORS \
   BOOST_PP_TUPLE_TO_LIST( \
      12, \
      (                                                                 \
       ( OpId   , id   , id   , 0, 0, 0, vm::JACOBIAN          , RankSame,false ), \
       ( OpDx   , dx   , dx   , 0, 1, 0, vm::JACOBIAN|vm::KB|vm::GRAD , RankSame,false ), \
       ( OpDy   , dy   , dy   , 0, 1, 1, vm::JACOBIAN|vm::KB|vm::GRAD , RankSame,false ), \
       ( OpDz   , dz   , dz   , 0, 1, 2, vm::JACOBIAN|vm::KB|vm::GRAD , RankSame,false ), \
       ( OpDn   , dn   , dn   , 0, 0, 0, vm::JACOBIAN|vm::KB|vm::NORMAL|vm::FIRST_DERIVATIVE|vm::FIRST_DERIVATIVE_NORMAL , RankSame,false ), \
       ( OpGrad , grad , grad , 0, 0, 0, vm::JACOBIAN|vm::KB|vm::GRAD , RankUp,true ), \
       ( OpDiv  , div  , div  , 1, 0, 0, vm::DIV|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false ), \
       ( OpCurl , curl , curl , 1, 0, 0, vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE , RankSame,false ), \
       ( OpCurlX, curlx, curlx, 1, 1, 0, vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false ), \
       ( OpCurlY, curly, curly, 1, 1, 1, vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false ), \
       ( OpCurlZ, curlz, curlz, 1, 1, 2, vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false ), \
       ( OpHess , hess , hess,  0, 0, 0, vm::JACOBIAN|vm::KB|vm::HESSIAN|vm::FIRST_DERIVATIVE , RankUp2,false ) \
      ) \
   ) \
   /**/
#
enum OperatorType { __TEST, __TRIAL, __VALUE };
# define VF_OP_TYPE_TYPE(T)         BOOST_PP_TUPLE_ELEM(6, 0, T)
# define VF_OP_TYPE_IS_GENERIC(T)   BOOST_PP_TUPLE_ELEM(6, 1, T)
# define VF_OP_TYPE_IS_TEST(T)      BOOST_PP_TUPLE_ELEM(6, 2, T)
# define VF_OP_TYPE_IS_TRIAL(T)     BOOST_PP_TUPLE_ELEM(6, 3, T)
# define VF_OP_TYPE_IS_VALUE(T)     BOOST_PP_TUPLE_ELEM(6, 4, T)
# define VF_OP_TYPE_SUFFIX(T)       BOOST_PP_TUPLE_ELEM(6, 5, T)
# define VF_OPERATORS_TYPE \
   BOOST_PP_TUPLE_TO_LIST( \
      4, \
      ( \
         ( OperatorType, 1, 0, 0, 0, g ),\
         ( __TEST      , 0, 1, 0, 0,   ),\
         ( __TRIAL     , 0, 0, 1, 0, t ),\
         ( __VALUE     , 0, 0, 0, 1, v ) \
      ) \
    ) \
   /**/
#

# define VF_TRIAL_NAME(T)         BOOST_PP_TUPLE_ELEM(3, 0, T)
# define VF_TRIAL_IS_GENERIC(T)   BOOST_PP_TUPLE_ELEM(3, 1, T)
# define VF_TRIAL_VALUE(T)        BOOST_PP_TUPLE_ELEM(3, 2, T)
#

# define VF_TRIAL \
   BOOST_PP_TUPLE_TO_LIST( \
      3, \
      ( \
         ( bool,  1, 0 ),\
         ( false, 0, 0 ),\
         ( true , 0, 1 ) \
      ) \
    ) \
   /**/
#
# define VF_OP_SPECIALIZATION_IF_NOT_GENERIC( E, T )        \
   BOOST_PP_IF ( BOOST_PP_NOT( VF_OP_TYPE_IS_GENERIC( T ) ),\
                 BOOST_PP_IDENTITY( <E ),                   \
                 BOOST_PP_EMPTY )()                         \
   BOOST_PP_IF ( BOOST_PP_NOT( VF_OP_TYPE_IS_GENERIC( T ) ),\
                 BOOST_PP_COMMA,                            \
                 BOOST_PP_EMPTY )()                         \
   BOOST_PP_IF ( BOOST_PP_NOT( VF_OP_TYPE_IS_GENERIC( T ) ),\
                 BOOST_PP_IDENTITY( VF_OP_TYPE_TYPE( T )> ),\
                 BOOST_PP_EMPTY )()                         \
   /**/
#
# define VF_OP_SWITCH( T, A, B )                \
  BOOST_PP_IF ( T,                              \
                BOOST_PP_IDENTITY( A ),         \
                BOOST_PP_IDENTITY( B )          \
              )()                               \
  /**/
# define VF_OP_SWITCH_ELSE_EMPTY( T, A )        \
  BOOST_PP_IF ( T,                              \
                BOOST_PP_IDENTITY( A ),         \
                BOOST_PP_EMPTY                  \
              )()                               \
  /**/
#
# define VF_OP_ADD_COMP( O )                                            \
  BOOST_PP_IF ( VF_OPERATOR_TYPE_IS_COMP( O ),                          \
                BOOST_PP_IDENTITY( ( VF_OPERATOR_TYPE_COMP( O ) ) ),    \
                BOOST_PP_EMPTY                                          \
              )()                                                       \
  /**/
#
#
#define VF_OP_TYPE_OBJECT(T)                                    \
   BOOST_PP_IF ( VF_OP_TYPE_IS_GENERIC( T ),                    \
                 BOOST_PP_IDENTITY( sw ),                       \
                 BOOST_PP_IDENTITY( VF_OP_TYPE_TYPE( T ) ) )()  \
  /**/
#
#
# /* Generates code for all binary operators and integral type pairs. */
# define VF_ARRAY_OPERATOR(_, OT) \
      VF_ARRAY_OPERATOR_CODE OT   \
   /**/
#define VF_FUSION1(T) BOOST_PP_IF( BOOST_PP_NOT(VF_OP_TYPE_IS_VALUE( T )), \
                                BOOST_PP_IDENTITY(typedef typename fusion::result_of::value_at_key<map_basis_context_type), \
                                BOOST_PP_EMPTY )()    \
                                   BOOST_PP_IF( BOOST_PP_NOT(VF_OP_TYPE_IS_VALUE( T )), \
                                BOOST_PP_COMMA,                         \
                                BOOST_PP_EMPTY )()                      \
                                   BOOST_PP_IF( BOOST_PP_NOT(VF_OP_TYPE_IS_VALUE( T )), \
                                 BOOST_PP_IDENTITY(key_type>::type::element_type basis_context_type), \
                                BOOST_PP_IDENTITY( typedef boost::none_t basis_context_type ))()                      \
                   /**/

#define VF_FUSION2(T) BOOST_PP_IF( BOOST_PP_NOT(VF_OP_TYPE_IS_VALUE( T )), \
                                BOOST_PP_IDENTITY(typedef typename fusion::result_of::value_at_key<map_basis_context_type), \
                                BOOST_PP_EMPTY )()    \
                                   BOOST_PP_IF( BOOST_PP_NOT(VF_OP_TYPE_IS_VALUE( T )), \
                                BOOST_PP_COMMA,                         \
                                BOOST_PP_EMPTY )()                      \
                                   BOOST_PP_IF( BOOST_PP_NOT(VF_OP_TYPE_IS_VALUE( T )), \
                                 BOOST_PP_IDENTITY(key_type>::type::pointer basis_context_ptrtype), \
                                   BOOST_PP_IDENTITY( typedef boost::none_t basis_context_ptrtype ) )() \
                   /**/

#define VF_ARRAY_OPERATOR_CODE(O,T)                                     \
    template <class Element                                             \
              BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ), BOOST_PP_COMMA, BOOST_PP_EMPTY )() \
        BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ), BOOST_PP_IDENTITY( VF_OP_TYPE_TYPE( T ) sw ), BOOST_PP_EMPTY )() > \
    class VF_OPERATOR_NAME( O ) VF_OP_SPECIALIZATION_IF_NOT_GENERIC( Element, T ) \
        {                                                               \
        public:                                                         \
                                                                        \
            static const size_type context = VF_OPERATOR_CONTEXT( O );  \
                                                                        \
            typedef Element element_type;                               \
            typedef boost::shared_ptr<element_type> element_ptrtype;    \
            typedef VF_OPERATOR_NAME( O )<element_type, VF_OP_TYPE_OBJECT(T)> this_type; \
            typedef this_type self_type;                                \
                                                                        \
            typedef typename element_type::functionspace_type functionspace_type; \
            typedef typename functionspace_type::reference_element_type* fe_ptrtype; \
            typedef typename functionspace_type::reference_element_type fe_type; \
            typedef typename functionspace_type::value_type value_type; \
            static const uint16_type rank = fe_type::rank;              \
            static const uint16_type nComponents1 = fe_type::nComponents1; \
            static const uint16_type nComponents2 = fe_type::nComponents2; \
                                                                        \
            template<typename Func>                                     \
                struct HasTestFunction                                  \
            {                                                           \
                static const bool result = VF_OP_SWITCH( BOOST_PP_OR( VF_OP_TYPE_IS_TRIAL( T ), VF_OP_TYPE_IS_VALUE( T ) ), false , (boost::is_same<Func,fe_type>::value) ); \
            };                                                          \
                                                                        \
            template<typename Func>                                     \
                struct HasTrialFunction                                 \
            {                                                           \
                static const bool result = VF_OP_SWITCH( VF_OP_TYPE_IS_TRIAL( T ), (boost::is_same<Func,fe_type>::value), false ); \
            };                                                          \
                                                                        \
                                                                        \
            VF_OPERATOR_NAME( O ) ( element_type const& v )             \
                : _M_v ( v )                                            \
            {                                                           \
                if ( VF_OP_TYPE_IS_VALUE( T ) )                         \
                    v.updateGlobalValues();                             \
                Debug( 5051 ) << "[" BOOST_PP_STRINGIZE(VF_OPERATOR_NAME( O )) "] default constructor\n"; \
            }                                                           \
            VF_OPERATOR_NAME( O )( VF_OPERATOR_NAME( O ) const& op )    \
                : _M_v ( op._M_v )                                      \
            {                                                           \
                Debug( 5051 ) << "[" BOOST_PP_STRINGIZE(VF_OPERATOR_NAME( O )) "] copy constructor\n"; \
            }                                                           \
                                                                        \
            element_type const& e() const { return _M_v; }              \
            template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t> \
                struct tensor                                           \
            {                                                           \
                typedef this_type expression_type;                      \
                typedef BOOST_PP_CAT( Basis_,BOOST_PP_CAT(VF_OP_SWITCH( BOOST_PP_NOT( VF_OP_TYPE_IS_TRIAL( T ) ), i ,j ), _t)) map_basis_context_type; \
                typedef typename mpl::if_<mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>, \
                    typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >, \
                    mpl::identity<detail::gmc<0> >,                     \
                    mpl::identity<detail::gmc<1> > >::type,             \
                    typename mpl::if_<fusion::result_of::has_key<map_basis_context_type, detail::gmc<0> >, \
                    mpl::identity<detail::gmc<0> >,                     \
                    mpl::identity<detail::gmc<1> > >::type>::type::type key_type; \
                typedef typename mpl::if_<fusion::result_of::has_key<map_basis_context_type, detail::gmc<0> >, \
                    mpl::identity<detail::gmc<0> >,                     \
                    mpl::identity<detail::gmc<1> > >::type::type basis_context_key_type;  \
                typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type; \
                BOOST_MPL_ASSERT_MSG( ( fusion::result_of::has_key<map_basis_context_type, basis_context_key_type >::value ), INVALID_BASISMAP_OP, (map_basis_context_type, key_type, basis_context_key_type, mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>, mpl::bool_<VF_OP_TYPE_IS_TRIAL( T )> )); \
                typedef typename mpl::if_<mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>, \
                    mpl::identity<mpl::int_<0> >,                       \
                    mpl::identity<typename fusion::result_of::value_at_key<map_basis_context_type,basis_context_key_type>::type::element_type > >::type::type basis_context_type; \
                typedef typename mpl::if_<mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>, \
                    mpl::identity<mpl::int_<0> >,                       \
                    mpl::identity<typename fusion::result_of::value_at_key<map_basis_context_type,basis_context_key_type>::type::pointer > >::type::type basis_context_ptrtype; \
                                                                        \
                typedef typename element_type::value_type value_type;   \
                typedef boost::multi_array<value_type,3> array_type;    \
                                                                        \
                typedef value_type result_type;                         \
                typedef typename element_type::polyset_type function_rank_type; \
                typedef typename VF_OPERATOR_RANK( O )<function_rank_type>::type return_value_type; \
                                                                        \
                typedef typename mpl::if_<mpl::equal_to<mpl::int_<return_value_type::rank>, \
                    mpl::int_<0> >,                                     \
                mpl::identity<Shape<gmc_type::NDim, Scalar, false> >,   \
                typename mpl::if_<mpl::equal_to<mpl::int_<return_value_type::rank>, \
                mpl::int_<1> >,                                         \
                    mpl::identity<Shape<gmc_type::NDim, Vectorial, VF_OPERATOR_TRANSPOSE(O)> >, \
                                  mpl::identity<Shape<gmc_type::NDim, Tensor2, false> > >::type>::type::type shape; \
                typedef typename fe_type::PreCompute pc_type;           \
                                                                        \
                                                                        \
                template<typename E>\
                struct ttt {                                                            \
                    typedef typename mpl::if_<boost::is_same<E,mpl::int_<0> >, \
                                              mpl::identity<functionspace_type>, \
                                              mpl::identity<basis_context_type> >::type::type type; \
                };                                                      \
                static const bool dim_ok  = (VF_OPERATOR_TYPE_COMP(O) < gmc_type::NDim); \
                static const bool fe_ok  = mpl::if_<mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>, \
                    mpl::bool_<true>,                                   \
                    boost::is_same<typename ttt<basis_context_type>::type::reference_element_type, fe_type> >::type::value; \
                struct is_zero {                                        \
                    /*static const bool value = !(dim_ok && fe_ok);*/   \
                    static const bool value = false;                    \
                };                                                      \
                                                                        \
                static const uint16_type rank = return_value_type::rank+1; \
                static const uint16_type nComponents = return_value_type::nComponents; \
                                                                        \
                tensor( this_type const& expr,                          \
                        Geo_t const& geom,                              \
                        Basis_i_t const& VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TEST( T ), fev ), \
                        Basis_j_t const& VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TRIAL( T ), feu ) ) \
                    :                                                   \
                    _M_expr( expr ),                                    \
                    _M_fec( VF_OP_SWITCH( VF_OP_TYPE_IS_TEST( T ),      \
                                          fusion::at_key<basis_context_key_type>( fev ).get() , \
                                          VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TRIAL( T ), \
                                                                   fusion::at_key<basis_context_key_type>( feu ).get() ) ) ), \
                    _M_np( fusion::at_key<key_type>( geom )->nPoints() ), \
                    _M_pc( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ), \
                    _M_loc(VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), expr.e().BOOST_PP_CAT(VF_OPERATOR_TERM( O ),Extents)(*fusion::at_key<key_type>( geom )) ) ) \
                        { update( geom ); }                             \
                tensor( this_type const& expr,                          \
                        Geo_t const& geom,                              \
                        Basis_i_t const& VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TEST( T ), fev  ) ) \
                    :                                                   \
                    _M_expr( expr ),                                    \
                    _M_fec( VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TEST( T ), \
                                                     fusion::at_key<basis_context_key_type>( fev ).get() ) ), \
                    _M_np( fusion::at_key<key_type>( geom )->nPoints() ), \
                    _M_pc( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ), \
                    _M_loc(VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), expr.e().BOOST_PP_CAT(VF_OPERATOR_TERM( O ),Extents)(*fusion::at_key<key_type>( geom )) ) ) \
                        { update( geom ); }                             \
                tensor( this_type const& expr,                          \
                        Geo_t const& geom )                             \
                    :                                                   \
                    _M_expr( expr ),                                    \
                    _M_np( fusion::at_key<key_type>( geom )->nPoints() ), \
                    _M_pc( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ), \
                    _M_loc(VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), expr.e().BOOST_PP_CAT(VF_OPERATOR_TERM( O ),Extents)(*fusion::at_key<key_type>( geom )) ) ) \
                        {                                               \
                            update( geom );                             \
                            BOOST_MPL_ASSERT_MSG( VF_OP_TYPE_IS_VALUE( T ), INVALID_CALL_TO_CONSTRUCTOR, ()); \
                        }                                               \
                                                                        \
                void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) \
                {                                                       \
                    update( geom, mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                       \
                void update( Geo_t const& geom, Basis_i_t const& /*fev*/ ) \
                {                                                       \
                    update( geom, mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                       \
                void update( Geo_t const& geom )                        \
                {                                                       \
                    update( geom, mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                       \
                void update( Geo_t const& geom, mpl::bool_<true> )      \
                {                                                       \
                    if ( dim_ok )                                       \
                        {                                               \
                            if ( fusion::at_key<key_type>( geom )->elementIsAFace() ) \
                                _M_pc.update( fusion::at_key<key_type>( geom )->xRefs() ); \
                            std::fill( _M_loc.data(), _M_loc.data()+_M_loc.num_elements(), value_type( 0 ) ); \
                            _M_expr.e().VF_OPERATOR_SYMBOL( O )( *fusion::at_key<key_type>( geom ), _M_pc, _M_loc ); \
                        }                                               \
                }                                                       \
                void update( Geo_t const& geom, mpl::bool_<false> )     \
                {                                                       \
                    Life::detail::ignore_unused_variable_warning(geom); \
                }                                                       \
                template<typename IndexI, typename IndexJ>              \
                    result_type                                         \
                    evalijq( IndexI const& i,                           \
                             IndexJ const& VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TRIAL( T ), j ), \
                             uint16_type c1, uint16_type c2, uint16_type q  ) const \
                {                                                       \
                    Life::detail::ignore_unused_variable_warning(i);    \
                    return evaliq( VF_OP_SWITCH( BOOST_PP_NOT( VF_OP_TYPE_IS_TRIAL( T ) ),i,j), c1, c2, q ); \
                }                                                       \
                                                                        \
                template<typename IndexI, typename IndexJ, int PatternContext> \
                    result_type                                         \
                    evalijq( IndexI const& i,                           \
                             IndexJ const& j,                           \
                             uint16_type c1, uint16_type c2, uint16_type q, \
                             mpl::int_<PatternContext> ) const          \
                {                                                       \
                    return evalijq( i, j, c1, c2, q );                  \
                }                                                       \
                template<typename IndexI>                               \
                    result_type                                         \
                    evaliq( IndexI const& i, uint16_type c1, uint16_type c2, uint16_type q  ) const \
                {                                                       \
                    return evaliq_( i, c1, c2, q, mpl::bool_<dim_ok && fe_ok>() ); \
                }                                                       \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q ) const \
                {                                                       \
                    BOOST_MPL_ASSERT_MSG( VF_OP_TYPE_IS_VALUE( T ), INVALID_CALL_TO_EVALQ, ()); \
                    return evalq( c1, c2, q, mpl::int_<shape::rank>() ); \
                }                                                       \
            private:                                                    \
                template<typename IndexI>                               \
                    result_type                                         \
                    evaliq_( IndexI const& /*i*/,                       \
                             uint16_type /*c1*/, uint16_type /*c2*/,    \
                             int /*q*/,                                 \
                             mpl::bool_<false> ) const                  \
                {                                                       \
                    return 0;                                           \
                }                                                       \
                template<typename IndexI>                               \
                    result_type                                         \
                    evaliq_( IndexI const& i, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<true> ) const \
                {                                                       \
                    return evaliq__( i, c1, c2, q, mpl::bool_<true>(), mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                       \
                template<typename IndexI>                               \
                    result_type                                         \
                    evaliq__( IndexI const& /*i*/, uint16_type c1, uint16_type c2, uint16_type q, \
                              mpl::bool_<true>, mpl::bool_<true> ) const \
                {                                                       \
                    return evalq( c1, c2, q, mpl::int_<shape::rank>() ); \
                }                                                       \
                template<typename IndexI>                               \
                    result_type                                         \
                    evaliq__( IndexI const& i, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<true>, mpl::bool_<false> ) const \
                {                                                       \
                    return  _M_fec->VF_OPERATOR_TERM( O )( i, c1, c2, q ); \
                }                                                       \
                                                                        \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0> ) const \
                {                                                       \
                    Life::detail::ignore_unused_variable_warning(c1);   \
                    Life::detail::ignore_unused_variable_warning(c2);   \
                    return _M_loc[0][0][q];                           \
                }                                                       \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1> ) const \
                {                                                       \
                    return evalq( c1, c2, q, mpl::int_<1>(), mpl::bool_<shape::is_transposed>() ); \
                }                                                       \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1>, mpl::bool_<false> ) const \
                {                                                       \
                    Life::detail::ignore_unused_variable_warning(c1);   \
                    Life::detail::ignore_unused_variable_warning(c2);   \
                    return _M_loc[c1][0][q];                            \
                }                                                       \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1>, mpl::bool_<true> ) const \
                {                                                       \
                    Life::detail::ignore_unused_variable_warning(c1);   \
                    Life::detail::ignore_unused_variable_warning(c2);   \
                    return _M_loc[0][c2][q];                            \
                }                                                       \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> ) const \
                {                                                       \
                    return _M_loc[c1][c2][q];                           \
                }                                                       \
                this_type const& _M_expr;                               \
                basis_context_ptrtype _M_fec;                           \
                const uint16_type _M_np;                                \
                pc_type _M_pc;                                          \
                array_type _M_loc;                                      \
                /*typename element_type::BOOST_PP_CAT( VF_OPERATOR_TERM( O ), _type) _M_loc;*/ \
            };                                                          \
                                                                        \
        protected:                                                      \
            VF_OPERATOR_NAME( O ) () {}                                 \
            element_type const& _M_v;                                   \
        };                                                              \
    template <class ELEM                                                \
              BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ),  BOOST_PP_COMMA, BOOST_PP_EMPTY )() \
        BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ),  BOOST_PP_IDENTITY( VF_OP_TYPE_TYPE( T ) sw ), BOOST_PP_EMPTY )() > \
    inline Expr< VF_OPERATOR_NAME( O )< ELEM, VF_OP_TYPE_OBJECT(T)> >   \
    BOOST_PP_CAT( VF_OPERATOR_SYMBOL(O), VF_OP_TYPE_SUFFIX(T) )( ELEM const& expr ) \
        {                                                               \
            typedef VF_OPERATOR_NAME( O )< ELEM, VF_OP_TYPE_OBJECT(T)> expr_t; \
            return Expr< expr_t >(  expr_t(expr) );                     \
        }                                                               \
    /**/
#
//
// Generate the code
//
BOOST_PP_LIST_FOR_EACH_PRODUCT(VF_ARRAY_OPERATOR, 2, (VF_OPERATORS, VF_OPERATORS_TYPE))
/// \endcond
}
}

#endif /* __LIFE_VF_OPERATORS_HPP */
