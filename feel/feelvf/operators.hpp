
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-17

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
   \file operators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-17
 */
#if !defined( __FEELPP_VF_OPERATORS_HPP )
#define __FEELPP_VF_OPERATORS_HPP 1

#include <feel/feelpoly/quadmapped.hpp>

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


namespace Feel
{
namespace vf
{
/// \cond detail
# /* Accessors for the operator datatype. */
# define VF_OPERATOR_NAME(O)           BOOST_PP_TUPLE_ELEM(11, 0, O)
# define VF_OPERATOR_SYMBOL(O)         BOOST_PP_TUPLE_ELEM(11, 1, O)
# define VF_OPERATOR_TERM(O)           BOOST_PP_TUPLE_ELEM(11, 2, O)
# // 0 means any rank, 1 or more is the rank for which the operator is defined
# define VF_OPERATOR_TYPE_RANK_DEF(O)  BOOST_PP_TUPLE_ELEM(11, 3, O)
# define VF_OPERATOR_TYPE_IS_COMP(O)   BOOST_PP_TUPLE_ELEM(11, 4, O)
# define VF_OPERATOR_TYPE_COMP(O)      BOOST_PP_TUPLE_ELEM(11, 5, O)
# define VF_OPERATOR_CONTEXT(O)        BOOST_PP_TUPLE_ELEM(11, 6, O)
# define VF_OPERATOR_RANK(O)           BOOST_PP_TUPLE_ELEM(11, 7, O)
# define VF_OPERATOR_TRANSPOSE(O)      BOOST_PP_TUPLE_ELEM(11, 8, O)
# define VF_OPERATOR_DIFFORDERIM(O)    BOOST_PP_TUPLE_ELEM(11, 9, O)
# define VF_OPERATOR_TERMINAL(O)       BOOST_PP_TUPLE_ELEM(11,10, O)



# /* List of applicative operators. */
# define VF_OPERATORS \
   BOOST_PP_TUPLE_TO_LIST( \
      13, \
      (                                                                 \
          ( OpId   , id   , id   , 0, 0, 0, vm::JACOBIAN          , RankSame,false, 0, 1 ), \
          ( OpDx   , dx   , dx   , 0, 1, 0, vm::JACOBIAN|vm::KB|vm::GRAD , RankSame,false,-1,1 ), \
          ( OpDy   , dy   , dy   , 0, 1, 1, vm::JACOBIAN|vm::KB|vm::GRAD , RankSame,false,-1,1 ), \
          ( OpDz   , dz   , dz   , 0, 1, 2, vm::JACOBIAN|vm::KB|vm::GRAD , RankSame,false,-1,1 ), \
          ( OpDn   , dn   , dn   , 0, 0, 0, vm::JACOBIAN|vm::KB|vm::NORMAL|vm::FIRST_DERIVATIVE|vm::FIRST_DERIVATIVE_NORMAL , RankSame,false,-1,1 ), \
          ( OpGrad , grad , grad , 0, 0, 0, vm::JACOBIAN|vm::KB|vm::GRAD , RankUp,true,-1,1 ), \
          ( OpDiv  , div  , div  , 1, 0, 0, vm::DIV|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false,-1,1 ), \
          ( OpCurl , curl , curl , 1, 0, 0, vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE , RankSame,false,-1,1 ), \
          ( OpCurlX, curlx, curlx, 1, 1, 0, vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false,-1,1 ), \
          ( OpCurlY, curly, curly, 1, 1, 1, vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false,-1,1 ), \
          ( OpCurlZ, curlz, curlz, 1, 1, 2, vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false,-1,1 ), \
          ( OpHess , hess , hess,  0, 0, 0, vm::JACOBIAN|vm::KB|vm::HESSIAN|vm::FIRST_DERIVATIVE , RankUp2,false,-2,1 ), \
          ( OpLap  , laplacian, laplacian,  0, 0, 0, vm::JACOBIAN|vm::KB|vm::LAPLACIAN|vm::FIRST_DERIVATIVE , RankUp2,false,-2,1 ) \
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
                                 BOOST_PP_IDENTITY(key_type>::type::element_type* basis_context_ptrtype), \
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
            typedef typename functionspace_type::mortar_fe_type mortar_fe_type; \
            typedef typename functionspace_type::geoelement_type geoelement_type; \
            typedef typename functionspace_type::gm_type gm_type; \
            typedef typename functionspace_type::value_type value_type; \
            typedef value_type evaluate_type;                           \
                                                                        \
            static const uint16_type rank = fe_type::rank;              \
            static const uint16_type nComponents1 = fe_type::nComponents1; \
            static const uint16_type nComponents2 = fe_type::nComponents2; \
            static const bool is_terminal = VF_OPERATOR_TERMINAL(O);    \
                                                                        \
            static const uint16_type imorder_test=element_type::functionspace_type::basis_type::nOrder + VF_OPERATOR_DIFFORDERIM(O); \
            static const uint16_type imorder = (imorder_test==invalid_uint16_type_value)?0:imorder_test; \
            static const bool imIsPoly = true;                          \
                                                                        \
            template<typename Func>                                     \
                struct HasTestFunction                                  \
            {                                                           \
                static const bool result = VF_OP_SWITCH( BOOST_PP_OR( VF_OP_TYPE_IS_TRIAL( T ), VF_OP_TYPE_IS_VALUE( T ) ), false , \
                                                         (boost::is_same<Func,fe_type>::value||(element_type::is_mortar&&boost::is_same<Func,mortar_fe_type>::value)) ); \
            };                                                          \
                                                                        \
            template<typename Func>                                     \
                struct HasTrialFunction                                 \
            {                                                           \
                static const bool result = VF_OP_SWITCH( VF_OP_TYPE_IS_TRIAL( T ), \
                                                         (boost::is_same<Func,fe_type>::value||(element_type::is_mortar&&boost::is_same<Func,mortar_fe_type>::value)), false ); \
            };                                                          \
                                                                        \
                                                                        \
            VF_OPERATOR_NAME( O ) ( element_type const& v, bool useInterpWithConfLoc=false ) \
              : M_v ( boost::cref(v) ),                                 \
              M_useInterpWithConfLoc( useInterpWithConfLoc )            \
            {                                                           \
                if ( VF_OP_TYPE_IS_VALUE( T ) )                         \
                    v.updateGlobalValues();                             \
                DVLOG(2) << "[" BOOST_PP_STRINGIZE(VF_OPERATOR_NAME( O )) "] default constructor\n"; \
            }                                                           \
            VF_OPERATOR_NAME( O )( VF_OPERATOR_NAME( O ) const& op )    \
              : M_v ( op.M_v ),                                         \
              M_useInterpWithConfLoc( op.M_useInterpWithConfLoc )       \
                                                                        \
            {                                                           \
                DVLOG(2) << "[" BOOST_PP_STRINGIZE(VF_OPERATOR_NAME( O )) "] copy constructor\n"; \
            }                                                           \
            template<typename... TheExpr>                               \
            struct Lambda                                               \
            {                                                           \
                typedef this_type type;                                 \
            };                                                          \
            template<typename... TheExpr>                               \
                typename Lambda<TheExpr...>::type                       \
                operator()( TheExpr... e) { return *this; }             \
                                                                        \
            element_type const& e() const { return M_v; }              \
            bool useInterpWithConfLoc() const { return M_useInterpWithConfLoc; } \
            template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t> \
                struct tensor                                           \
            {                                                           \
                typedef this_type expression_type;                      \
                typedef BOOST_PP_CAT( Basis_,BOOST_PP_CAT(VF_OP_SWITCH( BOOST_PP_NOT( VF_OP_TYPE_IS_TRIAL( T ) ), i ,j ), _t)) map_basis_context_type; \
                typedef typename mpl::if_<mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>, \
                    typename mpl::if_<fusion::result_of::has_key<Geo_t,vf::detail::gmc<0> >, \
                    mpl::identity<vf::detail::gmc<0> >,                     \
                    mpl::identity<vf::detail::gmc<1> > >::type,             \
                    typename mpl::if_<fusion::result_of::has_key<map_basis_context_type,vf::detail::gmc<0> >, \
                    mpl::identity<vf::detail::gmc<0> >,                     \
                    mpl::identity<vf::detail::gmc<1> > >::type>::type::type key_type; \
                typedef typename mpl::if_<fusion::result_of::has_key<map_basis_context_type,vf::detail::gmc<0> >, \
                    mpl::identity<vf::detail::gmc<0> >,                     \
                    mpl::identity<vf::detail::gmc<1> > >::type::type basis_context_key_type;  \
                typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type; \
                typedef boost::shared_ptr<gmc_type> gmc_ptrtype;        \
                typedef typename gmc_type::gm_type gm_type;             \
                BOOST_MPL_ASSERT_MSG( ( fusion::result_of::has_key<map_basis_context_type, basis_context_key_type >::value ), INVALID_BASISMAP_OP, (map_basis_context_type, key_type, basis_context_key_type, mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>, mpl::bool_<VF_OP_TYPE_IS_TRIAL( T )> )); \
                typedef typename mpl::if_<mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>, \
                    mpl::identity<mpl::int_<0> >,                       \
                    mpl::identity<typename fusion::result_of::value_at_key<map_basis_context_type,basis_context_key_type>::type::element_type > >::type::type basis_context_type; \
                typedef typename mpl::if_<mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>, \
                    mpl::identity<mpl::int_<0> >,                       \
                    mpl::identity<typename fusion::result_of::value_at_key<map_basis_context_type,basis_context_key_type>::type::element_type* > >::type::type basis_context_ptrtype; \
                typedef typename element_type::value_type value_type;   \
                typedef typename matrix_node<value_type>::type matrix_node_type; \
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
                typedef boost::shared_ptr<pc_type> pc_ptrtype;          \
                typedef typename fe_type::template Context<context, fe_type, gm_type,geoelement_type,gmc_type::context> ctx_type; \
                typedef boost::shared_ptr<ctx_type> ctx_ptrtype;        \
                typedef Eigen::Matrix<value_type,shape::M,shape::N> loc_type; \
                typedef Eigen::Matrix<value_type,shape::M,shape::N> ret_type; \
                typedef boost::multi_array<loc_type,1> array_type;    \
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
                                                    mpl::bool_<true>,   \
                                                    typename mpl::or_<boost::is_same<typename ttt<basis_context_type>::type::reference_element_type, fe_type>, \
                                                                      mpl::and_<mpl::bool_<element_type::is_mortar>,boost::is_same<typename ttt<basis_context_type>::type::reference_element_type, mortar_fe_type> > \
                                                                      >::type >::type::value; \
                struct is_zero {                                        \
                    /*static const bool value = !(dim_ok && fe_ok);*/   \
                    static const bool value = false;                    \
                };                                                      \
                                                                        \
                static const uint16_type rank = return_value_type::rank+1; \
                static const uint16_type nComponents = return_value_type::nComponents; \
                                                                        \
                static const bool isSameGeo = boost::is_same<typename gmc_type::element_type,geoelement_type>::value; \
                                                                        \
                tensor( tensor const& t )                               \
                    :                                                   \
                    M_expr( t.M_expr ),                                 \
                    M_geot( new gmc_type( *t.M_geot ) ),                \
                    M_fec( VF_OP_SWITCH( VF_OP_TYPE_IS_VALUE( T ), , t.M_fec/*new basis_context_type( *t.M_fec )*/ ) ), \
                    M_np( M_geot->nPoints() ),                          \
                    M_pc( new pc_type( M_expr.e().functionSpace()->fe(), M_geot->xRefs() )), \
                    /*M_pcf(),*/                                        \
                    M_ctx( VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), (new ctx_type( M_expr.e().functionSpace()->fe(), M_geot, (pc_ptrtype const&)M_pc ) ) ) ), \
                    M_loc(VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), M_expr.e().BOOST_PP_CAT(VF_OPERATOR_TERM( O ),Extents)(*M_geot) ) ), \
                    M_zero( ret_type::Zero() ),                         \
                    M_did_init( t.M_did_init ),                         \
                    M_hasRelationMesh( t.M_hasRelationMesh ),           \
                    M_same_mesh( t.M_same_mesh )                        \
                        {                                               \
                        }                                               \
                                                                        \
                tensor( this_type const& expr,                          \
                        Geo_t const& geom,                              \
                        Basis_i_t const& VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TEST( T ), fev ), \
                        Basis_j_t const& VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TRIAL( T ), feu ) ) \
                    :                                                   \
                    M_expr( expr ),                                     \
                    M_geot( fusion::at_key<key_type>( geom ) ),         \
                    M_fec( VF_OP_SWITCH( VF_OP_TYPE_IS_TEST( T ),       \
                                         fusion::at_key<basis_context_key_type>( fev ).get() , \
                                         VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TRIAL( T ), \
                                                                  fusion::at_key<basis_context_key_type>( feu ).get() ) ) ), \
                    M_np( fusion::at_key<key_type>( geom )->nPoints() ), \
                    M_pc( this->createPcIfSameGeom(expr,geom, mpl::bool_<isSameGeo>() ) ) /*new pc_type( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ))*/, \
                    /*M_pcf(),*/                                        \
                    M_ctx( VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), \
                                                ( this->createCtxIfSameGeom(expr,geom, mpl::bool_<isSameGeo>() )) ) ), \
                    M_loc(VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), expr.e().BOOST_PP_CAT(VF_OPERATOR_TERM( O ),Extents)(*fusion::at_key<key_type>( geom )) ) ), \
                    M_zero( ret_type::Zero() ),                         \
                    M_did_init( false ),                                \
                    M_hasRelationMesh( fusion::at_key<key_type>( geom )->element().mesh()->isRelatedTo( expr.e().functionSpace()->mesh()) ), \
                    M_same_mesh( M_hasRelationMesh && isSameGeo )         \
                        {                                               \
                            if(!M_same_mesh)                            \
                                    expr.e().functionSpace()->mesh()->tool_localization()->updateForUse(); \
                        }                                               \
                tensor( this_type const& expr,                          \
                        Geo_t const& geom,                              \
                        Basis_i_t const& VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TEST( T ), fev  ) ) \
                    :                                                   \
                    M_expr( expr ),                                    \
                    M_geot( fusion::at_key<key_type>( geom ) ),         \
                    M_fec( VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TEST( T ), \
                                                     fusion::at_key<basis_context_key_type>( fev ).get() ) ), \
                    M_np( fusion::at_key<key_type>( geom )->nPoints() ), \
                    M_pc( this->createPcIfSameGeom(expr,geom, mpl::bool_<isSameGeo>() ) ) /*new pc_type( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) )*/, \
                    /*M_pcf(),*/                                        \
                    M_ctx( VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), \
                                                ( this->createCtxIfSameGeom(expr,geom, mpl::bool_<isSameGeo>() )) ) ), \
                    M_loc(VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), expr.e().BOOST_PP_CAT(VF_OPERATOR_TERM( O ),Extents)(*fusion::at_key<key_type>( geom )) ) ), \
                    M_zero( ret_type::Zero() ),                         \
                    M_did_init( false ),                                \
                    M_hasRelationMesh( fusion::at_key<key_type>( geom )->element().mesh()->isRelatedTo( expr.e().functionSpace()->mesh()) ), \
                    M_same_mesh( M_hasRelationMesh && isSameGeo )         \
                        {                                               \
                            if(!M_same_mesh)                            \
                                expr.e().functionSpace()->mesh()->tool_localization()->updateForUse(); \
                            /*update( geom );*/                         \
                        }                                               \
                tensor( this_type const& expr,                          \
                        Geo_t const& geom )                             \
                    :                                                   \
                    M_expr( expr ),                                     \
                    M_geot( fusion::at_key<key_type>( geom ) ),         \
                    M_np( fusion::at_key<key_type>( geom )->nPoints() ), \
                    M_pc( this->createPcIfSameGeom(expr,geom, mpl::bool_<isSameGeo>() ) ) /* new pc_type( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) )*/, \
                    /*M_pcf(),*/                                        \
                    M_ctx( VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), \
                                                    ( this->createCtxIfSameGeom(expr,geom, mpl::bool_<isSameGeo>() )) ) ), \
                    M_loc(VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), expr.e().BOOST_PP_CAT(VF_OPERATOR_TERM( O ),Extents)(*fusion::at_key<key_type>( geom )) ) ), \
                    M_zero( ret_type::Zero() ),                         \
                    M_did_init( false ),                                \
                    M_hasRelationMesh( fusion::at_key<key_type>( geom )->element().mesh()->isRelatedTo( expr.e().functionSpace()->mesh()) ), \
                    M_same_mesh( M_hasRelationMesh && isSameGeo )         \
                        {                                               \
                            if(!M_same_mesh)                            \
                                expr.e().functionSpace()->mesh()->tool_localization()->updateForUse(); \
                            /*update( geom ); */                        \
                            BOOST_MPL_ASSERT_MSG( VF_OP_TYPE_IS_VALUE( T ), INVALID_CALL_TO_CONSTRUCTOR, ()); \
                        }                                               \
                template<typename IM>                                   \
                    void init( IM const& im )                           \
                {                                                       \
                    M_did_init = true;                                  \
                    /*                                                  \
                    QuadMapped<IM> qm;                                  \
                    typedef typename QuadMapped<IM>::permutation_type permutation_type; \
                    typename QuadMapped<IM>::permutation_points_type ppts( qm( im ) ); \
                                                                        \
                    M_pcf.resize( im.nFaces() );                       \
                    for ( uint16_type __f = 0; __f < im.nFaces(); ++__f ) \
                        {                                               \
                            for( permutation_type __p( permutation_type::IDENTITY ); \
                                 __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p ) \
                                {                                       \
                                    M_pcf[__f][__p.value()] = pc_ptrtype(  new pc_type( M_expr.e().functionSpace()->fe(), \
                                                                                         ppts[__f].find(__p)->second ) ); \
                                }                                       \
                        }                                               \
                    */                                                  \
                    /* expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ), */ \
                }                                                       \
                                                                        \
                void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) \
                {                                                       \
                    update( geom, fev, feu, mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                       \
                void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu , mpl::bool_<true> ) \
                {                                                       \
                    update( geom, mpl::bool_<true>() );                 \
                }                                                       \
                void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu , mpl::bool_<false> ) \
                {                                                       \
                    if (M_same_mesh)                                    \
                        updateInCaseOfInterpolate( geom, fev, feu, mpl::bool_<false>() ); \
                    else                                                \
                        updateInCaseOfInterpolate( geom, fev, feu, mpl::bool_<true>() ); \
                }                                                       \
                void updateInCaseOfInterpolate( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu , mpl::bool_<false> ) \
                {                                                       \
                    /*nothing : always same context*/                   \
                }                                                       \
                void updateInCaseOfInterpolate( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu , mpl::bool_<true> ) \
                {   /*with interp*/                                     \
                    VF_OP_SWITCH( VF_OP_TYPE_IS_TEST( T ),              \
                                  M_fec =fusion::at_key<basis_context_key_type>( fev ).get() , \
                                  VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TRIAL( T ), \
                                                           M_fec = fusion::at_key<basis_context_key_type>( feu ).get() ) ) ; \
                }                                                       \
                void update( Geo_t const& geom, Basis_i_t const& fev )  \
                {                                                       \
                    update( geom, fev, mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                       \
                void update( Geo_t const& geom, Basis_i_t const& fev, mpl::bool_<true> ) \
                {                                                       \
                    update( geom, mpl::bool_<true>() );                 \
                }                                                       \
                void update( Geo_t const& geom, Basis_i_t const& fev , mpl::bool_<false>) \
                {                                                       \
                    if (M_same_mesh)                                    \
                        updateInCaseOfInterpolate( geom, fev, mpl::bool_<false>() ); \
                    else                                                \
                        updateInCaseOfInterpolate( geom, fev, mpl::bool_<true>() ); \
                }                                                       \
                void updateInCaseOfInterpolate( Geo_t const& geom, Basis_i_t const& fev, mpl::bool_<false> ) \
                {                                                       \
                    /*no interp*/                                       \
                }                                                       \
                void updateInCaseOfInterpolate( Geo_t const& geom, Basis_i_t const& fev,  mpl::bool_<true> ) \
                {   /*with interp*/                                     \
                    VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TEST( T ),   \
                                             M_fec = fusion::at_key<basis_context_key_type>( fev ).get() ) ; \
                }                                                       \
                template <typename CTX>                                 \
                    void updateContext( CTX const& ctx )                \
                {                                                       \
                    std::fill( M_loc.data(), M_loc.data()+M_loc.num_elements(), loc_type::Zero() ); \
                    M_expr.e().VF_OPERATOR_SYMBOL( O )( *ctx, M_loc ); \
                }                                                       \
                void update( Geo_t const& geom )                        \
                {                                                       \
                    /*BOOST_STATIC_ASSERT( dim_ok );*/                  \
                    update( geom, mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                       \
                void update( Geo_t const& geom, uint16_type face )      \
                {                                                       \
                    /*BOOST_STATIC_ASSERT( dim_ok );*/                  \
                    update( geom, face, mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                       \
                void update( Geo_t const& geom, uint16_type face1, mpl::bool_<true> ) \
                {                                                       \
                    std::fill( M_loc.data(), M_loc.data()+M_loc.num_elements(), loc_type::Zero() ); \
                    this->updateCtxFaceIfSameGeom(geom,mpl::bool_<isSameGeo>() ); \
                    if (M_same_mesh)                                   \
                        M_expr.e().VF_OPERATOR_SYMBOL( O )( *M_ctx, M_loc ); \
                    else  {                                             \
                        matrix_node_type __ptsreal = M_expr.e().ptsInContext(*fusion::at_key<key_type>( geom ), mpl::int_<2>()); \
                        M_expr.e().BOOST_PP_CAT(VF_OPERATOR_SYMBOL( O ),Interpolate)( /**M_ctx,*/ __ptsreal, M_loc, M_expr.useInterpWithConfLoc()/*false*//*true*/, \
                                                                                      fusion::at_key<key_type>( geom )->element().face( fusion::at_key<key_type>( geom )->faceId() ).vertices() ); \
                    }                                                   \
                }                                                       \
                void update( Geo_t const& geom, mpl::bool_<true> )      \
                {                                                       \
                    std::fill( M_loc.data(), M_loc.data()+M_loc.num_elements(), loc_type::Zero() ); \
                    this->updateCtxIfSameGeom(geom,mpl::bool_<isSameGeo>() ); \
                    if (M_same_mesh) {                                  \
                        /*std::cout << "\n idv no interp \n";*/         \
                        M_expr.e().VF_OPERATOR_SYMBOL( O )( *M_ctx, M_loc ); \
                    }                                                   \
                    else {                                              \
                        /*std::cout << "\n idv with interp \n";*/       \
                        matrix_node_type __ptsreal = M_expr.e().ptsInContext(*fusion::at_key<key_type>( geom ), mpl::int_<1>()); \
                        auto setOfPts = ( fusion::at_key<key_type>( geom )->faceId() != invalid_uint16_type_value )? \
                          fusion::at_key<key_type>( geom )->element().face( fusion::at_key<key_type>( geom )->faceId() ).vertices() : \
                          fusion::at_key<key_type>( geom )->element().vertices(); \
                        M_expr.e().BOOST_PP_CAT(VF_OPERATOR_SYMBOL( O ),Interpolate)( /**M_ctx,*/ __ptsreal, M_loc, M_expr.useInterpWithConfLoc()/*false*//*true*/, setOfPts  ); \
                    }                                                   \
                }                                                       \
                void update( Geo_t const& geom, mpl::bool_<false> )     \
                {                                                       \
                    /*std::cout << "\n idv no interp \n";*/             \
                    Feel::detail::ignore_unused_variable_warning(geom); \
                }                                                       \
                                                                        \
                ret_type const&                                         \
                    evalijq( uint16_type i,                             \
                             uint16_type VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TRIAL( T ), j ), \
                             uint16_type q  ) const                     \
                {                                                       \
                    Feel::detail::ignore_unused_variable_warning(i);    \
                    return evaliq( VF_OP_SWITCH( BOOST_PP_NOT( VF_OP_TYPE_IS_TRIAL( T ) ),i,j), q ); \
                }                                                       \
                    result_type                                         \
                        evalijq( uint16_type i,                         \
                                 uint16_type VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TRIAL( T ), j ), \
                             uint16_type c1, uint16_type c2, uint16_type q  ) const \
                {                                                       \
                    Feel::detail::ignore_unused_variable_warning(i);    \
                    return evaliq( VF_OP_SWITCH( BOOST_PP_NOT( VF_OP_TYPE_IS_TRIAL( T ) ),i,j), c1, c2, q ); \
                }                                                       \
                                                                        \
                template<int PatternContext> \
                    result_type                                         \
                    evalijq( uint16_type i,                           \
                             uint16_type j,                           \
                             uint16_type c1, uint16_type c2, uint16_type q, \
                             mpl::int_<PatternContext> ) const          \
                {                                                       \
                    return evalijq( i, j, c1, c2, q );                  \
                }                                                       \
                                                                        \
                ret_type const&                                         \
                    evaliq( uint16_type i, uint16_type q  ) const       \
                {                                                       \
                    return evaliq_( i, q, mpl::bool_<dim_ok && fe_ok>() ); \
                }                                                       \
                result_type                                             \
                    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q  ) const \
                {                                                       \
                    return evaliq_( i, c1, c2, q, mpl::bool_<dim_ok && fe_ok>() ); \
                }                                                       \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q ) const \
                {                                                       \
                    BOOST_MPL_ASSERT_MSG( VF_OP_TYPE_IS_VALUE( T ), INVALID_CALL_TO_EVALQ, ()); \
                    return evalq( c1, c2, q, mpl::int_<shape::rank>() ); \
                }                                                       \
                ret_type const&                                         \
                    evalq( uint16_type q ) const                        \
                {                                                       \
                    BOOST_MPL_ASSERT_MSG( VF_OP_TYPE_IS_VALUE( T ), INVALID_CALL_TO_EVALQ, ()); \
                    return M_loc[q];                                    \
                }                                                       \
            private:                                                    \
                                                                        \
                    result_type                                         \
                    evaliq_( uint16_type /*i*/,                       \
                             uint16_type /*c1*/, uint16_type /*c2*/,    \
                             int /*q*/,                                 \
                             mpl::bool_<false> ) const                  \
                {                                                       \
                    return 0;                                           \
                }                                                       \
                    ret_type const&                                            \
                        evaliq_( uint16_type /*i*/,                     \
                             int /*q*/,                                 \
                             mpl::bool_<false> ) const                  \
                {                                                       \
                    return M_zero;                                      \
                }                                                       \
                                                                        \
                    ret_type const&                                     \
                        evaliq_( uint16_type i, uint16_type q, mpl::bool_<true> ) const \
                    {                                                   \
                        return evaliq__( i, q, mpl::bool_<true>(), mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                    }                                                   \
                    result_type                                         \
                    evaliq_( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<true> ) const \
                {                                                       \
                    return evaliq__( i, c1, c2, q, mpl::bool_<true>(), mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                       \
                                                                        \
                    ret_type  const&                                    \
                        evaliq__( uint16_type /*i*/,  uint16_type q,    \
                                  mpl::bool_<true>, mpl::bool_<true> ) const \
                    {                                                   \
                        return M_loc[q];                                \
                    }                                                   \
                    result_type                                         \
                    evaliq__( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q, \
                              mpl::bool_<true>, mpl::bool_<true> ) const \
                {                                                       \
                    return evalq( c1, c2, q, mpl::int_<shape::rank>() ); \
                }                                                       \
                                                                        \
                    ret_type const&                                         \
                        evaliq__( uint16_type i, uint16_type q, mpl::bool_<true>, mpl::bool_<false> ) const \
                    {                                                   \
                        return M_fec->VF_OPERATOR_TERM( O )( i, q ); \
                    }                                                   \
                    result_type                                         \
                    evaliq__( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<true>, mpl::bool_<false> ) const \
                    {                                                   \
                        return M_fec->VF_OPERATOR_TERM( O )( i, c1, c2, q ); \
                    }                                                   \
                                                                        \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0> ) const \
                {                                                       \
                    Feel::detail::ignore_unused_variable_warning(c1);   \
                    Feel::detail::ignore_unused_variable_warning(c2);   \
                    return M_loc[q](0,0);                           \
                }                                                       \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1> ) const \
                {                                                       \
                    return evalq( c1, c2, q, mpl::int_<1>(), mpl::bool_<shape::is_transposed>() ); \
                }                                                       \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1>, mpl::bool_<false> ) const \
                {                                                       \
                    Feel::detail::ignore_unused_variable_warning(c1);   \
                    Feel::detail::ignore_unused_variable_warning(c2);   \
                    return M_loc[q](c1,0);                                \
                }                                                       \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1>, mpl::bool_<true> ) const \
                {                                                       \
                    Feel::detail::ignore_unused_variable_warning(c1);   \
                    Feel::detail::ignore_unused_variable_warning(c2);   \
                    return M_loc[q](0,c2);                                \
                }                                                       \
                result_type                                             \
                    evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> ) const \
                {                                                       \
                    return M_loc[q](c1,c2);                               \
                }                                                      \
                pc_ptrtype createPcIfSameGeom(this_type const& expr, Geo_t const& geom,mpl::bool_<true>) \
                {                                                       \
                    return pc_ptrtype (new pc_type( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ); \
                }                                                       \
                pc_ptrtype createPcIfSameGeom(this_type const& expr, Geo_t const& geom,mpl::bool_<false>) \
                {                                                       \
                    return pc_ptrtype();                                \
                }                                                       \
                ctx_ptrtype createCtxIfSameGeom(this_type const& expr, Geo_t const& geom,mpl::bool_<true>) \
                {                                                       \
                    return ctx_ptrtype( new ctx_type( expr.e().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_ptrtype const&)M_pc ) ); \
                }                                                       \
                ctx_ptrtype createCtxIfSameGeom(this_type const& expr, Geo_t const& geom,mpl::bool_<false>) \
                {                                                       \
                    return ctx_ptrtype( /*new ctx_type( )*/ );          \
                }                                                       \
                void updateCtxIfSameGeom(Geo_t const& geom, mpl::bool_<true> )    \
                {   if (fusion::at_key<key_type>( geom )->faceId() != invalid_uint16_type_value ) /*face case*/ \
                        M_pc->update(fusion::at_key<key_type>( geom )->pc()->nodes() ); \
                    M_ctx->update( fusion::at_key<key_type>( geom ),  (pc_ptrtype const&) M_pc ); \
                }                                                       \
                void updateCtxIfSameGeom(Geo_t const& geom, mpl::bool_<false> ) \
                {                                                       \
                }                                                       \
                void updateCtxFaceIfSameGeom(Geo_t const& geom, mpl::bool_<true> )    \
                {                                                       \
                    /*uint16_type face = fusion::at_key<key_type>( geom )->faceId(); \
                    uint16_type perm = fusion::at_key<key_type>( geom )->permutation().value(); \
                    M_ctx->update( fusion::at_key<key_type>( geom ), (pc_ptrtype const&) M_pcf[face][perm] );*/ \
                    M_pc->update(fusion::at_key<key_type>( geom )->pc()->nodes() ); \
                    M_ctx->update( fusion::at_key<key_type>( geom ), (pc_ptrtype const&) M_pc ); \
                }                                                       \
                void updateCtxFaceIfSameGeom(Geo_t const& geom, mpl::bool_<false> ) \
                {                                                       \
                }                                                       \
                this_type const& M_expr;                               \
                gmc_ptrtype M_geot;                                    \
                basis_context_ptrtype M_fec;                           \
                const uint16_type M_np;                                \
                pc_ptrtype M_pc;                                       \
                /*std::vector<std::map<uint16_type, pc_ptrtype> > M_pcf;*/ \
                ctx_ptrtype M_ctx;                                      \
                array_type M_loc;                                      \
                ret_type M_zero;                                        \
                /*typename element_type::BOOST_PP_CAT( VF_OPERATOR_TERM( O ), _type) M_loc;*/ \
                bool M_did_init;                                        \
                const bool M_hasRelationMesh;                           \
                const bool M_same_mesh;                                 \
            };                                                          \
                                                                        \
        protected:                                                      \
            VF_OPERATOR_NAME( O ) () {}                                 \
            boost::reference_wrapper<const element_type>  M_v;         \
            bool M_useInterpWithConfLoc;                                     \
        };                                                              \
    template <class ELEM                                                \
              BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ),  BOOST_PP_COMMA, BOOST_PP_EMPTY )() \
        BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ),  BOOST_PP_IDENTITY( VF_OP_TYPE_TYPE( T ) sw ), BOOST_PP_EMPTY )() > \
    inline Expr< VF_OPERATOR_NAME( O )< ELEM, VF_OP_TYPE_OBJECT(T)> >   \
    BOOST_PP_CAT( VF_OPERATOR_SYMBOL(O), VF_OP_TYPE_SUFFIX(T) )( ELEM const& expr,bool useInterpWithConfLoc=false ) \
        {                                                               \
            typedef VF_OPERATOR_NAME( O )< ELEM, VF_OP_TYPE_OBJECT(T)> expr_t; \
            return Expr< expr_t >(  expr_t(expr,useInterpWithConfLoc) ); \
        }                                                               \
                                                                        \
    template <class ELEM                                                \
              BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ),  BOOST_PP_COMMA, BOOST_PP_EMPTY )() \
        BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ),  BOOST_PP_IDENTITY( VF_OP_TYPE_TYPE( T ) sw ), BOOST_PP_EMPTY )() > \
    inline Expr< VF_OPERATOR_NAME( O )< ELEM, VF_OP_TYPE_OBJECT(T)> >   \
    BOOST_PP_CAT( VF_OPERATOR_SYMBOL(O), VF_OP_TYPE_SUFFIX(T) )( boost::shared_ptr<ELEM> expr,bool useInterpWithConfLoc=false ) \
        {                                                               \
            typedef VF_OPERATOR_NAME( O )< ELEM, VF_OP_TYPE_OBJECT(T)> expr_t; \
            return Expr< expr_t >(  expr_t(*expr,useInterpWithConfLoc) ); \
        }                                                               \

/**/
#
//
// Generate the code
//
BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_ARRAY_OPERATOR, 2, ( VF_OPERATORS, VF_OPERATORS_TYPE ) )
/// \endcond

// try to add operators to Python library
/*
#define VF_DEF(_,OT) \
    VF_DEF2 OT;

#define VF_DEF2(O,T) \
    typedef VF_OPERATOR_NAME(O)<ELEM,VF_OP_TYPEOBJECT(T)> BOOST_PP_CAT(expr_t,BOOST_PP_CAT(VF_OPERATOR_SYMBOL(O) , VF_OP_TYPE_SUFFIX(T)))

BOOST_PP_LIST_FOR_EACH_PRODUCT(VF_DEF,2,(VF_OPERATORS,VF_OPERATORS_TYPE))
*/

}
}

#endif /* __FEELPP_VF_OPERATORS_HPP */
