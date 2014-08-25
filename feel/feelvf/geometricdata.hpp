/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-03-14

  Copyright (C) 2006 EPFL
  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2014 Feel++ Consortium

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
   \file geometricdata.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-03-14
 */
#ifndef __GeometricData_H
#define __GeometricData_H 1

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


#include <feel/feelvf/expr.hpp>


namespace Feel { namespace vf {

/// \cond detail
# /* Information about C operators */
#
# /* Accessors for the operator datatype. */
# define VF_GD_SYMBOL(O)      BOOST_PP_TUPLE_ELEM(7, 0, O)
# define VF_GD_NAME(O)        BOOST_PP_TUPLE_ELEM(7, 1, O)
# define VF_GD_DIM(O)         BOOST_PP_TUPLE_ELEM(7, 2, O)
# define VF_GD_CONTEXT(O)     BOOST_PP_TUPLE_ELEM(7, 3, O)
# define VF_GD_RETURN(O)      BOOST_PP_TUPLE_ELEM(7, 4, O)
# define VF_GD_VALUE(O)       BOOST_PP_TUPLE_ELEM(7, 5, O)
# define VF_GD_IMORDER(O)     BOOST_PP_TUPLE_ELEM(7, 6, O)

#

const size_type jn = vm::JACOBIAN|vm::NORMAL;
const size_type jkbn = vm::JACOBIAN|vm::KB|vm::NORMAL;
const size_type jt = vm::JACOBIAN|vm::NORMAL|vm::TANGENT;
const size_type jp = vm::JACOBIAN|vm::POINT;
const size_type jkp = vm::KB|vm::JACOBIAN|vm::POINT;

# /* List of applicative unary operators. */
#if 1
# define VF_GD                                                          \
   BOOST_PP_TUPLE_TO_LIST(                                              \
       22,                                                              \
      (                                                                 \
       ( N       , GDN       , 0, jkbn, Vectorial, M_gmc->unitNormal( q )[ c1 ] , 0), \
       ( Nx      , GDNx      , 0, jkbn, Scalar   , M_gmc->unitNormal( q )[ 0 ]  , 0), \
       ( Ny      , GDNy      , 1, jkbn, Scalar   , M_gmc->unitNormal( q )[ 1 ]  , 0), \
       ( Nz      , GDNz      , 2, jkbn, Scalar   , M_gmc->unitNormal( q )[ 2 ]  , 0), \
       ( Nref    , GDNref    , 0, 0, Vectorial, M_gmc->refNormal( q )[ c1 ] , 0), \
       ( normalNorm, GDnormalNorm, 0, 0, Scalar, M_gmc->normalNorm( q ) , 0), \
       ( pT      , GDpT      , 0, jt, Tensor2, M_gmc->projectorTangent( c1, c2, q), 0), \
       ( T       , GDT       , 0, jt, Vectorial, M_gmc->unitTangent( q )[ c1 ], 0), \
       ( Tx      , GDTx      , 0, jt, Scalar   , M_gmc->unitTangent( q )[ 0 ] , 0), \
       ( Ty      , GDTy      , 1, jt, Scalar   , M_gmc->unitTangent( q )[ 1 ] , 0), \
       ( Tz      , GDTz      , 2, jt, Scalar   , M_gmc->unitTangent( q )[ 2 ] , 0), \
       ( detJ    , GDDetJ    , 0, jp, Scalar   , M_gmc->J( q )                , 0), \
       ( J       , GDJ       , 0, jkp,Tensor2  , M_gmc->K( c1, c2, q )        , 0), \
       ( JinvT   , GDJinv    , 0, jkp,Tensor2  , M_gmc->B( c1, c2, q )        , 0), \
       ( P       , GDP       , 0, jp, Vectorial, M_gmc->xReal( q )[ c1 ]      , 1), \
       ( Px      , GDPx      , 0, jp, Scalar   , M_gmc->xReal( q )[ 0 ]       , 1), \
       ( Py      , GDPy      , 1, jp, Scalar   , M_gmc->xReal( q )[ 1 ]       , 1), \
       ( Pz      , GDPz      , 2, jp, Scalar   , M_gmc->xReal( q )[ 2 ]       , 1), \
       ( C       , GDC       , 0, jp, Vectorial, M_gmc->barycenterReal()[c1]  , 0), \
       ( Cx      , GDCx      , 0, jp, Scalar   , M_gmc->barycenterReal()[0]   , 0), \
       ( Cy      , GDCy      , 1, jp, Scalar   , M_gmc->barycenterReal()[1]   , 0), \
       ( Cz      , GDCz      , 2, jp, Scalar   , M_gmc->barycenterReal()[2]   , 0) \
          )                                                             \
       )
/**/
# define VF_GD2                                                         \
    BOOST_PP_TUPLE_TO_LIST(                                             \
        10,                                                              \
        (                                                               \
            ( h       , GDH       , 0, 0 , Scalar   , M_gmc->h()                   , 0), \
            ( hMin    , GDHMin    , 0, 0 , Scalar   , M_gmc->hMin()                , 0), \
            ( hFace   , GDHFace   , 0, 0 , Scalar   , M_gmc->hFace()               , 0), \
            ( meas    , GDMeas    , 0, 0 , Scalar   , M_gmc->meas()                , 0), \
            ( measPEN , GDMeasPEN , 0, 0 , Scalar   , M_gmc->measurePointElementNeighbors(), 0), \
            ( nPEN    , GDNPEN    , 0, 0 , Scalar   , M_gmc->element().numberOfPointElementNeighbors(), 0), \
            ( measFace, GDHMeasFace,0, 0 , Scalar   , M_gmc->measFace()            , 0), \
            ( eid     , GDEid     , 0, 0 , Scalar   , M_gmc->id()                  , 0), \
            ( emarker , GDEmarker , 0, 0 , Scalar   , M_gmc->marker().value()      , 0), \
            ( emarker2, GDEmarker2, 0, 0 , Scalar   , M_gmc->marker2().value()     , 0) \
            )                                                           \
        )                                                               \
/**/
#else
# define VF_GD                                                          \
   BOOST_PP_TUPLE_TO_LIST(                                              \
      20,                                                               \
      (                                                                 \
       ( N       , GDN       , 0, jn, Vectorial, M_gmc->unitNormal( q )[ c1 ] , 0), \
       ( Nx      , GDNx      , 0, jn, Scalar   , M_gmc->unitNormal( q )[ 0 ]  , 0), \
       ( Ny      , GDNy      , 1, jn, Scalar   , M_gmc->unitNormal( q )[ 1 ]  , 0), \
       ( Nz      , GDNz      , 2, jn, Scalar   , M_gmc->unitNormal( q )[ 2 ]  , 0), \
       ( T       , GDT       , 0, jn, Vectorial, M_gmc->unitTangent( q )[ c1 ], 0), \
       ( Tx      , GDTx      , 0, jn, Scalar   , M_gmc->unitTangent( q )[ 0 ] , 0), \
       ( Ty      , GDTy      , 1, jn, Scalar   , M_gmc->unitTangent( q )[ 1 ] , 0), \
       ( Tz      , GDTz      , 2, jn, Scalar   , M_gmc->unitTangent( q )[ 2 ] , 0), \
       ( P       , GDP       , 0, jp, Vectorial, M_gmc->xReal( q )[ c1 ]      , 1), \
       ( Px      , GDPx      , 0, jp, Scalar   , M_gmc->xReal( q )[ 0 ]       , 1), \
       ( Py      , GDPy      , 1, jp, Scalar   , M_gmc->xReal( q )[ 1 ]       , 1), \
       ( Pz      , GDPz      , 2, jp, Scalar   , M_gmc->xReal( q )[ 2 ]       , 1), \
       ( C       , GDC       , 0, jp, Vectorial, M_gmc->barycenterReal()[c1]  , 0), \
       ( Cx      , GDCx      , 0, jp, Scalar   , M_gmc->barycenterReal()[0]   , 0), \
       ( Cy      , GDCy      , 1, jp, Scalar   , M_gmc->barycenterReal()[1]   , 0), \
       ( Cz      , GDCz      , 2, jp, Scalar   , M_gmc->barycenterReal()[2]   , 0), \
       ( h       , GDH       , 0, 0 , Scalar   , M_gmc->h()                   , 0), \
       ( hFace   , GDHFace   , 0, 0 , Scalar   , M_gmc->hFace()               , 0), \
       ( eid     , GDEid     , 0, 0 , Scalar   , M_gmc->id()                  , 0), \
       ( emarker , GDEmarker , 0, 0 , Scalar   , M_gmc->marker().value()      , 0) \
      )                                                                          \
   )                                                                             \
/**/

#endif
#

# /* Generates code for all binary operators and integral type pairs. */
# define VF_ARRAY_GD(_, O) \
      VF_ARRAY_GD_CODE O   \
   /**/

#define VF_ARRAY_GD_CODE(O)                                             \
    class VF_GD_NAME( O )                                               \
    {                                                                   \
    public:                                                             \
                                                                        \
        static const size_type context = VF_GD_CONTEXT( O );            \
        static const bool is_terminal = false;                          \
                                                                        \
        static const uint16_type imorder = VF_GD_IMORDER(O);            \
        static const bool imIsPoly = true;                              \
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
        typedef VF_GD_NAME(O) this_type;                                \
        typedef double value_type;                                      \
        typedef value_type evaluate_type;                               \
                                                                        \
        VF_GD_NAME(O) ()                                                \
        {                                                               \
        }                                                               \
        VF_GD_NAME(O) ( VF_GD_NAME(O) const& /*__vf*/ )                 \
        {                                                               \
        }                                                               \
        template<typename... TheExpr>                                      \
        struct Lambda                                                   \
        {                                                               \
            typedef this_type type;                                     \
        };                                                              \
        template<typename... TheExpr>                                      \
        typename Lambda<TheExpr...>::type                                  \
        operator()( TheExpr... e  ) { return this_type(); }         \
                                                                        \
                                                                        \
        template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t> \
            struct tensor                                               \
        {                                                               \
            typedef this_type expression_type;                          \
            typedef typename mpl::if_<fusion::result_of::has_key<Geo_t,vf::detail::gmc<0> >, \
                mpl::identity<vf::detail::gmc<0> >,                         \
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;  \
            typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype; \
            typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type; \
            typedef typename gmc_type::value_type value_type;           \
            typedef  value_type evaluate_type;                          \
            typedef VF_GD_RETURN(O)<gmc_type::NDim> return_value_type;  \
            typedef Shape<gmc_type::NDim, VF_GD_RETURN(O), false> shape; \
                                                                        \
            struct is_zero { static const bool value = false; };        \
                                                                        \
            tensor( this_type const& /*expr*/,                          \
                    Geo_t const& geom,                                  \
                    Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) \
                :                                                       \
                M_gmc( fusion::at_key<key_type>( geom ).get() )        \
                {}                                                      \
            tensor( this_type const& /*expr*/,                          \
                    Geo_t const& geom,                                  \
                    Basis_i_t const& /*fev*/ )                          \
                :                                                       \
                M_gmc( fusion::at_key<key_type>( geom ).get() )        \
                {}                                                      \
            tensor( this_type const& /*expr*/,                          \
                    Geo_t const& geom )                                 \
                :                                                       \
                M_gmc( fusion::at_key<key_type>( geom ).get() )        \
                {}                                                      \
            template<typename IM>                                       \
                void init( IM const& im )                               \
            {                                                           \
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
                M_gmc = fusion::at_key<key_type>( geom ).get();        \
            }                                                           \
            void update( Geo_t const& geom, uint16_type face )          \
            {                                                           \
                /*BOOST_STATIC_ASSERT( dim_ok );*/                      \
                update( geom );                                         \
            }                                                           \
            template<typename CTX>                                      \
                void updateContext( CTX const& ctx )                    \
            {                                                           \
                M_gmc = ctx->gmContext().get();                        \
            }                                                           \
                                                                        \
                value_type                                              \
                evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return evalq_( c1, c2, q, mpl::bool_<(VF_GD_DIM(O) < gmc_type::NDim)>(), mpl::int_<100>() ); \
            }                                                           \
                template<int PatternContext>                            \
                value_type                                              \
                evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<PatternContext> ) const \
            {                                                           \
                return evalq_( c1, c2, q, mpl::bool_<(VF_GD_DIM(O) < gmc_type::NDim)>(), mpl::int_<100>() ); \
            }                                                           \
                                                                        \
                value_type                                              \
                evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return evalq_( c1, c2, q, mpl::bool_<(VF_GD_DIM(O) < gmc_type::NDim)>(), mpl::int_<100>() ); \
            }                                                           \
            value_type                                                  \
                evalq( uint16_type c1, uint16_type c2, uint16_type q ) const   \
            {                                                           \
                return evalq_( c1, c2, q, mpl::bool_<(VF_GD_DIM(O) < gmc_type::NDim)>(), mpl::int_<100>() ); \
            }                                                           \
        private:                                                        \
            value_type                                                  \
                evalq_( uint16_type, uint16_type, uint16_type /*q*/, mpl::bool_<false>, mpl::int_<100> ) const \
            {                                                           \
                return 0;                                               \
            }                                                           \
            value_type                                                  \
                evalq_( uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<true>, mpl::int_<100> ) const \
            {                                                           \
                Feel::detail::ignore_unused_variable_warning(q);        \
                Feel::detail::ignore_unused_variable_warning(c1);       \
                Feel::detail::ignore_unused_variable_warning(c2);       \
                return VF_GD_VALUE( O );                                \
            }                                                           \
            gmc_ptrtype M_gmc;                                         \
        };                                                              \
    };                                                                  \
    inline                                                              \
    Expr< VF_GD_NAME( O )>                                              \
    VF_GD_SYMBOL(O)()                                                   \
    {                                                                   \
        typedef VF_GD_NAME( O ) expr_t;                                 \
        return Expr< expr_t >(  expr_t() );                             \
    }                                                                   \
    /**/

//
// Generate the code
//
BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_ARRAY_GD, 1, ( VF_GD ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_ARRAY_GD, 1, ( VF_GD2 ) )
/// \endcond
} // vf
} // Feel
#endif /* __GeometricData_H */
