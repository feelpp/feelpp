/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-17

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2016 Feel++ Consortium

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

#include <feel/feelvf/exproptionalconcat.hpp>
#include <feel/feelvf/one.hpp>

namespace Feel
{
struct ContextGeometricBase;

namespace vf
{

namespace detail
{
template <typename ValueType, int ShapeM,int ShapeN>
void convertEigenMatrixTensor( Eigen::Tensor<ValueType,2> const& input, Eigen::Matrix<ValueType,ShapeM,ShapeN> & output )
{
    output = Eigen::Map< const Eigen::Matrix<ValueType,ShapeM,ShapeN> >( input.data() );
}
template <typename ValueType, int ShapeM,int ShapeN>
void convertEigenMatrixTensor( Eigen::Tensor<ValueType,3> const& input, Eigen::Matrix<ValueType,ShapeM,ShapeN> & output )
{
    //if ( input.size() == ShapeM*ShapeN )
    if ( input.dimension(2)  == 1 )
        output = Eigen::Map< const Eigen::Matrix<ValueType,ShapeM,ShapeN> >( input.data() );
    else
        CHECK( false ) << "TODO dim "<< input.dimension(0) << ","<< input.dimension(1) << ","<<input.dimension(2) <<"\n";
}

template <typename ValueType, long ShapeM,long ShapeN>
void convertEigenMatrixTensor( Eigen::TensorFixedSize<ValueType,Eigen::Sizes<ShapeM,ShapeN>> const& input, Eigen::Matrix<ValueType,ShapeM,ShapeN> & output )
{
    output = Eigen::Map< const Eigen::Matrix<ValueType,ShapeM,ShapeN> >( input.data() );
}
template <typename ValueType, long ShapeM,long ShapeN>
void convertEigenMatrixTensor( Eigen::TensorFixedSize<ValueType,Eigen::Sizes<ShapeM,ShapeN,1>> const& input, Eigen::Matrix<ValueType,ShapeM,ShapeN> & output )
{
    //if ( input.size() == ShapeM*ShapeN )
    if ( input.dimension(2)  == 1 )
        output = Eigen::Map< const Eigen::Matrix<ValueType,ShapeM,ShapeN> >( input.data() );
    else
        CHECK( false ) << "TODO dim "<< input.dimension(0) << ","<< input.dimension(1) << ","<<input.dimension(2) <<"\n";
}


template <typename ValueType, long ShapeM,long ShapeN>
Eigen::Map<const Eigen::Matrix<ValueType,ShapeM,ShapeN>> convertEigenMatrixTensor( Eigen::TensorFixedSize<ValueType,Eigen::Sizes<ShapeM,ShapeN>> const& input )
{
    return Eigen::Map< const Eigen::Matrix<ValueType,ShapeM,ShapeN> >( input.data() );
}
template <typename ValueType, long ShapeM,long ShapeN>
Eigen::Map<const Eigen::Matrix<ValueType,ShapeM,ShapeN>>  convertEigenMatrixTensor( Eigen::TensorFixedSize<ValueType,Eigen::Sizes<ShapeM,ShapeN,1>> const& input )
{
    return Eigen::Map< const Eigen::Matrix<ValueType,ShapeM,ShapeN> >( input.data() );
}

template <typename ValueType, long ShapeM,long ShapeN,long ShapeP>
Eigen::Map<const Eigen::Matrix<ValueType,ShapeM*ShapeP,ShapeN>>  convertEigenMatrixTensor( Eigen::TensorFixedSize<ValueType,Eigen::Sizes<ShapeM,ShapeN,ShapeP>> const& input )
{
    return Eigen::Map< const Eigen::Matrix<ValueType,ShapeM*ShapeP,ShapeN> >( input.data() );
}

}

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
      16, \
      (                                                                 \
       ( OpId   , id   , id   , 0, 0, 0, 0          , RankSame,false, 0, 1 ), \
       ( OpN    , normal    , normalComponent    , 1, 0, 0, vm::NORMAL_COMPONENT|vm::KB|vm::NORMAL , RankDown,false, 0, 1 ), \
       ( OpDx   , dx   , dx   , 0, 1, 0, vm::KB|vm::GRAD , RankSame,false,-1,1 ), \
       ( OpDy   , dy   , dy   , 0, 1, 1, vm::KB|vm::GRAD , RankSame,false,-1,1 ), \
       ( OpDz   , dz   , dz   , 0, 1, 2, vm::KB|vm::GRAD , RankSame,false,-1,1 ), \
       ( OpDn   , dn   , dn   , 0, 0, 0, vm::KB|vm::NORMAL|vm::FIRST_DERIVATIVE|vm::FIRST_DERIVATIVE_NORMAL , RankSame,false,-1,1 ), \
       ( OpGrad , grad , grad , 0, 0, 0, vm::KB|vm::GRAD , RankUp,true,-1,1 ), \
       ( OpSymmGrad , symm_grad , symmetricGradient , 1, 0, 0, vm::KB|vm::GRAD|vm::SYMM , RankUp,true,-1,1 ), \
       ( OpDiv  , div  , div  , 1, 0, 0, vm::DIV|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false,-1,1 ), \
       ( OpCurl , curl , curl , 1, 0, 0, vm::CURL|vm::KB|vm::FIRST_DERIVATIVE , RankCurl,false,-1,1 ), \
       ( OpCurlX, curlx, curlx, 1, 1, 0, vm::CURL|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false,-1,1 ), \
       ( OpCurlY, curly, curly, 1, 1, 1, vm::CURL|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false,-1,1 ), \
       ( OpCurlZ, curlz, curlz, 1, 1, 2, vm::CURL|vm::KB|vm::FIRST_DERIVATIVE , RankDown,false,-1,1 ), \
       ( OpHess , hess , hess,  0, 0, 0, vm::KB|vm::HESSIAN|vm::SECOND_DERIVATIVE , RankUp2,false,-2,1 ), \
       ( OpLap  , laplacian, laplacian,  0, 0, 0, vm::KB|vm::LAPLACIAN|vm::SECOND_DERIVATIVE , RankSame,false,-2,1 ), \
       ( OpTrace  , trace, trace,  0, 0, 0, vm::TRACE , Rank0,false,0,1 ) \
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


# /* Generates code for all binary operators and integral type pairs. */
# define VF_ARRAY_OPERATOR_DECLARATION(_, OT) \
      VF_ARRAY_OPERATOR_CODE_DECLARATION OT   \
   /**/

#define VF_ARRAY_OPERATOR_CODE_DECLARATION(O,T)                                     \
    template <class Element                                             \
    BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ), BOOST_PP_COMMA, BOOST_PP_EMPTY )() \
              BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ), BOOST_PP_IDENTITY( VF_OP_TYPE_TYPE( T ) sw ), BOOST_PP_EMPTY )() > \
    class VF_OPERATOR_NAME( O ) VF_OP_SPECIALIZATION_IF_NOT_GENERIC( Element, T ); \
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
            typedef Element element_type;                               \
            typedef std::shared_ptr<element_type> element_ptrtype;    \
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
                                                                        \
            static const uint16_type rank = fe_type::rank;              \
            static const uint16_type nComponents1 = fe_type::nComponents1; \
            static const uint16_type nComponents2 = fe_type::nComponents2; \
            static const bool is_terminal = VF_OPERATOR_TERMINAL(O);    \
            static const bool is_hdiv_conforming = Feel::is_hdiv_conforming_v<fe_type>; \
            static const bool is_hcurl_conforming = Feel::is_hcurl_conforming_v<fe_type>; \
            inline static const size_type context = (is_hdiv_conforming_v<fe_type>?(VF_OPERATOR_CONTEXT( O )|vm::JACOBIAN|vm::KB)\
                                                     :(is_hcurl_conforming_v<fe_type>?(VF_OPERATOR_CONTEXT( O )|vm::KB):VF_OPERATOR_CONTEXT( O )))|(VF_OP_TYPE_IS_VALUE( T )?vm::INTERPOLANT:vm::BASIS_FUNCTION); \
                                                                        \
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
            template<typename Func>                                     \
                static const bool has_test_basis = VF_OP_SWITCH( BOOST_PP_OR( VF_OP_TYPE_IS_TRIAL( T ), VF_OP_TYPE_IS_VALUE( T ) ), false , \
                                                                 (boost::is_same<Func,fe_type>::value||(element_type::is_mortar&&boost::is_same<Func,mortar_fe_type>::value)) ); \
            template<typename Func>                                     \
                static const bool has_trial_basis = VF_OP_SWITCH( VF_OP_TYPE_IS_TRIAL( T ), \
                                                                  (boost::is_same<Func,fe_type>::value||(element_type::is_mortar&&boost::is_same<Func,mortar_fe_type>::value)), false ); \
            using basis_t = typename mpl::if_<mpl::bool_<element_type::is_mortar>,mortar_fe_type,fe_type>::type; \
            using test_basis = VF_OP_SWITCH( BOOST_PP_OR( VF_OP_TYPE_IS_TRIAL( T ), VF_OP_TYPE_IS_VALUE( T ) ), std::nullptr_t, basis_t); \
            using trial_basis = VF_OP_SWITCH( VF_OP_TYPE_IS_TRIAL( T ), basis_t, std::nullptr_t ); \
                                                                        \
            template <uint16_type TheDim>                               \
            struct EvaluateShape                                        \
            {                                                           \
                typedef typename element_type::polyset_type function_rank_type; \
                typedef typename VF_OPERATOR_RANK( O )<function_rank_type>::type return_value_type; \
                                                                        \
                typedef typename mpl::if_<mpl::equal_to<mpl::int_<return_value_type::rank>, \
                                                        mpl::int_<0> >, \
                                          mpl::identity<Shape<TheDim, Scalar, false> >, \
                                          typename mpl::if_<mpl::equal_to<mpl::int_<return_value_type::rank>, \
                                                                          mpl::int_<1> >, \
                                                            mpl::identity<Shape<TheDim, Vectorial, VF_OPERATOR_TRANSPOSE(O)> >, \
                                                            mpl::identity<Shape<TheDim, Tensor2, false> > >::type>::type::type shape_type; \
                using type = shape_type;                                \
            };                                                          \
            using evaluate_type = Eigen::Matrix<value_type,             \
                                                EvaluateShape<functionspace_type::nRealDim>::type::M, \
                                                EvaluateShape<functionspace_type::nRealDim>::type::N>; \
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
            fe_ptrtype  fe() const { return M_v.fe(); }                 \
            uint16_type polynomialOrder() const {                       \
                int imorder_test = element_type::functionspace_type::basis_type::nOrder + VF_OPERATOR_DIFFORDERIM(O); \
                return (imorder_test<0)?0:imorder_test; \
            }                                                           \
                                                                        \
            bool isPolynomial() const { return true; }                  \
                                                                        \
            element_type const& e() const { return M_v; }              \
            bool useInterpWithConfLoc() const { return M_useInterpWithConfLoc; } \
                                                                        \
            template <typename TheFeType = fe_type>                     \
            evaluate_type evaluate(bool p,  worldcomm_ptr_t const& worldcomm, \
                                   typename std::enable_if_t< isP0Continuous<TheFeType>::result && \
                                   std::is_same_v< this_type, OpId<element_type, VF_OP_TYPE_OBJECT(T)> > && \
                                   VF_OP_TYPE_IS_VALUE( T ) >* = nullptr ) const \
            {                                                           \
                evaluate_type res = evaluate_type::Constant( 0. );      \
                if ( this->e().functionSpace()->nLocalDofWithGhost() ) \
                    for ( uint16_type c1=0 ; c1<nComponents1 ; ++c1 )   \
                        for ( uint16_type c2=0 ; c2<nComponents2 ; ++c2 ) \
                            res( c1,c2 ) = (this->e()(c2+nComponents2*c1)); \
                                                                        \
                /* TODO : not apply broadcast if all proc have the dof */ \
                /* TODO : check compatibility with arg worldcomm and the one used by the functionspace */ \
                if ( p )                                                \
                    mpi::broadcast( *worldcomm, res, this->e().map().procOnGlobalCluster( 0 ) ); \
                return res;                                             \
            }                                                           \
                                                                        \
            template <typename SymbolsExprType>                         \
                self_type applySymbolsExpr( SymbolsExprType const& se ) const \
            {                                                           \
                return *this;                                           \
            }                                                           \
                                                                        \
            template <typename TheSymbolExprType>                       \
            bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const \
            {                                                           \
                if constexpr ( isP0Continuous<fe_type>::result || functionspace_type::nDim == 0 ) \
                                 return false;                          \
                else if constexpr ( functionspace_type::nDim == 1 )     \
                                      return symb == "x";               \
                else if constexpr ( functionspace_type::nDim == 2 )     \
                                      return symb == "x" || symb == "y"; \
                else                                                    \
                    return symb == "x" || symb == "y" || symb == "z";   \
            }                                                           \
                                                                        \
            template <int diffOrder, typename TheSymbolExprType>        \
                auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr, \
                           TheSymbolExprType const& se ) const          \
            {                                                           \
                if constexpr ( std::is_same_v< this_type, OpId<element_type, VF_OP_TYPE_OBJECT(T)> > ) \
                    {                                                   \
                        std::map<std::string,int> compNameToIndex = { {"x",0},{"y",1},{"z",2} }; \
                        if constexpr ( fe_type::nComponents == 1 )      \
                        {                                               \
                            using diff_expr_type = std::decay_t<decltype( Feel::vf::expr( OpGrad<element_type, VF_OP_TYPE_OBJECT(T)>( this->e(), this->useInterpWithConfLoc() ) )(0,0) )>; \
                            auto res = exprOptionalConcat<diff_expr_type>(); \
                            if ( diffVariable == "x" || diffVariable == "y" || diffVariable == "z" ) \
                                res.expression().add( Feel::vf::expr( OpGrad<element_type, VF_OP_TYPE_OBJECT(T)>( this->e(), this->useInterpWithConfLoc() ) )(0,compNameToIndex[diffVariable]) ); \
                            return res;                                 \
                        }                                               \
                        else if constexpr ( fe_type::is_vectorial )     \
                            {                                           \
                                using diff_expr_type = std::decay_t<decltype( Feel::vf::expr( OpGrad<element_type, VF_OP_TYPE_OBJECT(T)>( this->e(), this->useInterpWithConfLoc() ) )*Feel::vf::one<functionspace_type::nRealDim>(0) )>; \
                                auto res = exprOptionalConcat<diff_expr_type>(); \
                                if ( diffVariable == "x" || diffVariable == "y" || diffVariable == "z" ) \
                                    res.expression().add( Feel::vf::expr( OpGrad<element_type, VF_OP_TYPE_OBJECT(T)>( this->e(), this->useInterpWithConfLoc() ) )*Feel::vf::one<functionspace_type::nRealDim>(compNameToIndex[diffVariable]) ); \
                                return res;                             \
                            }                                           \
                        else                                            \
                        {                                               \
                            CHECK( false ) << "TODO tensorial case";    \
                            return *this;                               \
                        }                                               \
                    }                                                   \
                else                                                    \
                {                                                       \
                    CHECK( false ) << "TODO";                           \
                    return *this;                                       \
                }                                                       \
            }                                                           \
                                                                        \
            template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t> \
                struct tensor                                           \
            {                                                           \
                typedef this_type expression_type;                      \
                static constexpr size_type context = expression_type::context; \
                                                                        \
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
                typedef std::shared_ptr<gmc_type> gmc_ptrtype;        \
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
                using shape = typename EvaluateShape<gmc_type::NDim>::type; \
                typedef typename fe_type::PreCompute pc_type;           \
                typedef std::shared_ptr<pc_type> pc_ptrtype;          \
                typedef typename fe_type::template Context<context, fe_type, gm_type,geoelement_type,/*gmc_type::*/context, gmc_type::subEntityCoDim> ctx_type; \
                typedef std::shared_ptr<ctx_type> ctx_ptrtype;        \
                /*typedef Eigen::Matrix<value_type,shape::M,shape::N> loc_type;*/ \
                using loc_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape::M,shape::N>>; \
                using eigen_matrix_mn_type = eigen_matrix_type<shape::M,shape::N,value_type>; \
                using ret_type = Eigen::Map<const eigen_matrix_mn_type>; \
                typedef boost::multi_array<loc_type,1> array_type;    \
                                                                        \
                                                                        \
                template<typename E>                                    \
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
                static const bool isSameGeo = std::is_same_v<typename gmc_type::element_type,geoelement_type>; \
                                                                        \
                tensor( tensor const& t )                               \
                    :                                                   \
                    M_expr( t.M_expr ),                                 \
                    M_geot( t.M_geot/*new gmc_type( *t.M_geot )*/ ),    \
                    M_fec( VF_OP_SWITCH( VF_OP_TYPE_IS_VALUE( T ), , t.M_fec/*new basis_context_type( *t.M_fec )*/ ) ), \
                    M_np( M_geot->nPoints() ),                          \
                    /*M_pc( new pc_type( M_expr.e().functionSpace()->fe(), M_geot->xRefs() )),*/ \
                    M_pc( t.M_pc ), \
                    /*M_pcf(),*/                                        \
                    /*M_ctx( VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), (new ctx_type( M_expr.e().functionSpace()->fe(), M_geot, (pc_ptrtype const&)M_pc ) ) ) ),*/ \
                    M_ctx( VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), t.M_ctx ) ), \
                    M_loc(VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_VALUE( T ), M_expr.e().BOOST_PP_CAT(VF_OPERATOR_TERM( O ),Extents)(*M_geot) ) ), \
                    M_mzero(),                                          \
                    M_zero( eigen_matrix_mn_type::Zero() ),                         \
                    M_returnEigenMatrix( eigen_matrix_mn_type::Zero() ),            \
                    M_did_init( t.M_did_init ),                         \
                    M_hasRelationMesh( t.M_hasRelationMesh ),           \
                    M_same_mesh( t.M_same_mesh )                        \
                        {                                               \
                            M_mzero.setZero();                          \
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
                    M_mzero(),                                          \
                    M_zero( eigen_matrix_mn_type::Zero() ),                         \
                    M_returnEigenMatrix( eigen_matrix_mn_type::Zero() ),            \
                    M_did_init( false ),                                \
                    M_hasRelationMesh( fusion::at_key<key_type>( geom )->element().mesh()->isRelatedTo( expr.e().functionSpace()->mesh()) ), \
                    M_same_mesh( M_hasRelationMesh && isSameGeo )         \
                        {                                               \
                            M_mzero.setZero();                          \
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
                    M_mzero(),                                          \
                    M_zero( eigen_matrix_mn_type::Zero() ),                         \
                    M_returnEigenMatrix( eigen_matrix_mn_type::Zero() ),            \
                    M_did_init( false ),                                \
                    M_hasRelationMesh( fusion::at_key<key_type>( geom )->element().mesh()->isRelatedTo( expr.e().functionSpace()->mesh()) ), \
                    M_same_mesh( M_hasRelationMesh && isSameGeo )         \
                        {                                               \
                            M_mzero.setZero();                          \
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
                    M_mzero(),                                          \
                    M_zero( eigen_matrix_mn_type::Zero() ),                         \
                    M_returnEigenMatrix( eigen_matrix_mn_type::Zero() ),            \
                    M_did_init( false ),                                \
                    M_hasRelationMesh( ( expr.e().functionSpace() && expr.e().functionSpace()->mesh() )? fusion::at_key<key_type>( geom )->element().mesh()->isRelatedTo( expr.e().functionSpace()->mesh()) : false ), \
                    M_same_mesh( M_hasRelationMesh && isSameGeo )         \
                        {                                               \
                            M_mzero.setZero();                          \
                            if(!M_same_mesh && expr.e().functionSpace() && expr.e().functionSpace()->mesh() ) \
                                expr.e().functionSpace()->mesh()->tool_localization()->updateForUse(); \
                            /*update( geom ); */                        \
                            BOOST_MPL_ASSERT_MSG( VF_OP_TYPE_IS_VALUE( T ), INVALID_CALL_TO_CONSTRUCTOR, ()); \
                        }                                               \
                template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType> \
                    tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, \
                            this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs ) \
                    :                                                   \
                    tensor( expr, geom, theInitArgs... )                \
                {}                                                      \
                                                                        \
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
                    if constexpr ( VF_OP_TYPE_IS_VALUE( T ) )           \
                    {                                                       \
                        update( geom );                                     \
                    }                                                       \
                    else                                                    \
                    {                                                       \
                        if (!M_same_mesh)                                   \
                        {   /*with interp*/                                 \
                            VF_OP_SWITCH( VF_OP_TYPE_IS_TEST( T ),          \
                                  M_fec =fusion::at_key<basis_context_key_type>( fev ).get() , \
                                  VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TRIAL( T ), \
                                                           M_fec = fusion::at_key<basis_context_key_type>( feu ).get() ) ) ; \
                        }                                                   \
                    }                                                       \
                }                                                           \
                void update( Geo_t const& geom, Basis_i_t const& fev )  \
                {                                                       \
                    if constexpr ( VF_OP_TYPE_IS_VALUE( T ) ) \
                    {                                                       \
                        update( geom );                                     \
                    }                                                       \
                    else                                                    \
                    {                                                       \
                        if (!M_same_mesh)                                   \
                        {   /*with interp*/                                 \
                            VF_OP_SWITCH_ELSE_EMPTY( VF_OP_TYPE_IS_TEST( T ),   \
                                    M_fec = fusion::at_key<basis_context_key_type>( fev ).get() ) ; \
                        }                                                   \
                    }                                                       \
                }                                                           \
                template <typename ... CTX>                             \
                    void updateContext( CTX const& ... ctx )            \
                {                                                       \
                    typedef typename boost::remove_reference<typename boost::remove_const< decltype(*M_expr.e().selectContext( ctx...) ) >::type >::type ctxspace_type; \
                    static const bool ctxspace_is_geometricspace = boost::is_base_of<ContextGeometricBase/*GeometricSpaceBase*/,ctxspace_type>::type::value; \
                                                                        \
                    std::fill( M_loc.data(), M_loc.data()+M_loc.num_elements(), M_mzero.constant(0.) ); \
                    /*M_expr.e().VF_OPERATOR_SYMBOL( O )( *ctx, M_loc );*/ \
                    /*M_expr.e().VF_OPERATOR_SYMBOL( O )( *M_expr.e().selectContext( ctx...), M_loc );*/ \
                    this->updateContext( M_expr.e().selectContext( ctx...), mpl::bool_<ctxspace_is_geometricspace>() ); \
                }                                                       \
                template <typename CTX>                                 \
                    void updateContext( CTX const& ctx, mpl::true_ /**/ ) \
                {                                                       \
                    Geo_t geom( fusion::make_pair<vf::detail::gmc<0> >( ctx->gmContext() ) ); \
                    M_pc->update( fusion::at_key<key_type>( geom )->xRefs() ); \
                    M_ctx->update( fusion::at_key<key_type>( geom ),  (pc_ptrtype const&) M_pc ); \
                    M_expr.e().VF_OPERATOR_SYMBOL( O )( *M_ctx, M_loc ); \
                }                                                       \
                template <typename CTX>                                 \
                    void updateContext( CTX const& ctx, mpl::false_ /**/ ) \
                {                                                       \
                    CHECK( ctx ) << "no ctx";CHECK( ctx->gmContext() ) << "no gmctx";M_expr.e().VF_OPERATOR_SYMBOL( O )( *ctx, M_loc ); \
                }                                                       \
                                                                        \
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
                    std::fill( M_loc.data(), M_loc.data()+M_loc.num_elements(), M_mzero.constant(0.) ); \
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
                    std::fill( M_loc.data(), M_loc.data()+M_loc.num_elements(), M_mzero.constant(0.) ); \
                    this->updateCtxIfSameGeom(geom,mpl::bool_<isSameGeo>() ); \
                    if (M_same_mesh) {                                  \
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
                template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType> \
                    void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, \
                                 Geo_t const& geom, const TheArgsType&... theUpdateArgs ) \
                {                                                       \
                    this->update( geom,theUpdateArgs... );              \
                }                                                       \
                ret_type                                                \
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
                ret_type                                                \
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
                ret_type                                                \
                    evalq( uint16_type q ) const                        \
                {                                                       \
                    BOOST_MPL_ASSERT_MSG( VF_OP_TYPE_IS_VALUE( T ), INVALID_CALL_TO_EVALQ, ()); \
                    return Feel::vf::detail::convertEigenMatrixTensor( M_loc[q] ); \
                    /*return M_returnEigenMatrix;*/                     \
                    /*return M_loc[q];*/                                \
                }                                                       \
            private:                                                    \
                                                                        \
                result_type                                             \
                    evaliq_( uint16_type /*i*/,                         \
                    uint16_type /*c1*/, uint16_type /*c2*/,             \
                    int /*q*/,                                          \
                    mpl::bool_<false> ) const                           \
                {                                                       \
                    return 0;                                           \
                }                                                       \
                ret_type                                                \
                    evaliq_( uint16_type /*i*/,                         \
                    int /*q*/,                                          \
                    mpl::bool_<false> ) const                           \
                {                                                       \
                    return ret_type(M_zero.data());                     \
                }                                                       \
                                                                        \
                ret_type                                                \
                    evaliq_( uint16_type i, uint16_type q, mpl::bool_<true> ) const \
                {                                                       \
                    return evaliq__( i, q, mpl::bool_<true>(), mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                           \
                result_type                                             \
                    evaliq_( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<true> ) const \
                {                                                       \
                    return evaliq__( i, c1, c2, q, mpl::bool_<true>(), mpl::bool_<VF_OP_TYPE_IS_VALUE( T )>() ); \
                }                                                       \
                                                                        \
                ret_type                                                \
                    evaliq__( uint16_type /*i*/,  uint16_type q,        \
                    mpl::bool_<true>, mpl::bool_<true> ) const          \
                    {                                                   \
                        return Feel::vf::detail::convertEigenMatrixTensor( M_loc[q]); \
                        /*return M_returnEigenMatrix;*/                 \
                        /*return M_loc[q];*/                            \
                    }                                                   \
                    result_type                                         \
                    evaliq__( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q, \
                              mpl::bool_<true>, mpl::bool_<true> ) const \
                {                                                       \
                    return evalq( c1, c2, q, mpl::int_<shape::rank>() ); \
                }                                                       \
                                                                        \
                    ret_type                                            \
                        evaliq__( uint16_type i, uint16_type q, mpl::bool_<true>, mpl::bool_<false> ) const \
                    {                                                   \
                        return Feel::vf::detail::convertEigenMatrixTensor( M_fec->VF_OPERATOR_TERM( O )( i, q ) ); \
                        /*return M_returnEigenMatrix;*/                 \
                        /*return M_fec->VF_OPERATOR_TERM( O )( i, q );*/ \
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
                    if ( expr.e().functionSpace() && expr.e().functionSpace()->fe() ) \
                        return pc_ptrtype (new pc_type( expr.e().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ); \
                    else                                                \
                        return pc_ptrtype();                            \
                }                                                       \
                pc_ptrtype createPcIfSameGeom(this_type const& expr, Geo_t const& geom,mpl::bool_<false>) \
                {                                                       \
                    return pc_ptrtype();                                \
                }                                                       \
                ctx_ptrtype createCtxIfSameGeom(this_type const& expr, Geo_t const& geom,mpl::bool_<true>) \
                {                                                       \
                    if ( expr.e().functionSpace()  && expr.e().functionSpace()->fe() ) \
                        return ctx_ptrtype( new ctx_type( expr.e().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_ptrtype const&)M_pc ) ); \
                    else                                                \
                        return ctx_ptrtype();                           \
                }                                                       \
                ctx_ptrtype createCtxIfSameGeom(this_type const& expr, Geo_t const& geom,mpl::bool_<false>) \
                {                                                       \
                    return ctx_ptrtype( /*new ctx_type( )*/ );          \
                }                                                       \
                void updateCtxIfSameGeom(Geo_t const& geom, mpl::bool_<true> )    \
                {                                                       \
                    if constexpr ( gmc_type::subEntityCoDim > 0 )       \
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
                loc_type M_mzero;                                       \
                eigen_matrix_mn_type M_zero;                            \
                mutable eigen_matrix_mn_type M_returnEigenMatrix;       \
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
    BOOST_PP_CAT( VF_OPERATOR_SYMBOL(O), VF_OP_TYPE_SUFFIX(T) )( ELEM const& expr,bool useInterpWithConfLoc=false, std::enable_if_t<std::is_base_of<FunctionSpaceBase::ElementBase,ELEM>::value>* = nullptr ) \
        {                                                               \
            typedef VF_OPERATOR_NAME( O )< ELEM, VF_OP_TYPE_OBJECT(T)> expr_t; \
            return Expr< expr_t >(  expr_t(expr,useInterpWithConfLoc) ); \
        }                                                               \
                                                                        \
    template <class ELEM                                                \
              BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ),  BOOST_PP_COMMA, BOOST_PP_EMPTY )() \
        BOOST_PP_IF( VF_OP_TYPE_IS_GENERIC( T ),  BOOST_PP_IDENTITY( VF_OP_TYPE_TYPE( T ) sw ), BOOST_PP_EMPTY )() > \
    inline Expr< VF_OPERATOR_NAME( O )< ELEM, VF_OP_TYPE_OBJECT(T)> >   \
    BOOST_PP_CAT( VF_OPERATOR_SYMBOL(O), VF_OP_TYPE_SUFFIX(T) )( std::shared_ptr<ELEM> expr,bool useInterpWithConfLoc=false, std::enable_if_t<std::is_base_of<FunctionSpaceBase::ElementBase,ELEM>::value>* = nullptr ) \
        {                                                               \
            typedef VF_OPERATOR_NAME( O )< ELEM, VF_OP_TYPE_OBJECT(T)> expr_t; \
            return Expr< expr_t >(  expr_t(*expr,useInterpWithConfLoc) ); \
        }                                                               \

/**/
#
//
// Generate the code
//
BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_ARRAY_OPERATOR_DECLARATION, 2, ( VF_OPERATORS, VF_OPERATORS_TYPE ) )
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
