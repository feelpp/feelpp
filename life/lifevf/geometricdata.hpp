/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-03-14

  Copyright (C) 2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file geometricdata.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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

namespace Life
{
namespace vf
{
/// \cond detail
# /* Information about C operators */
#
# /* Accessors for the operator datatype. */
# define VF_GD_SYMBOL(O)      BOOST_PP_TUPLE_ELEM(6, 0, O)
# define VF_GD_NAME(O)        BOOST_PP_TUPLE_ELEM(6, 1, O)
# define VF_GD_DIM(O)         BOOST_PP_TUPLE_ELEM(6, 2, O)
# define VF_GD_CONTEXT(O)     BOOST_PP_TUPLE_ELEM(6, 3, O)
# define VF_GD_RETURN(O)      BOOST_PP_TUPLE_ELEM(6, 4, O)
# define VF_GD_VALUE(O)       BOOST_PP_TUPLE_ELEM(6, 5, O)

#

const size_type jn = vm::JACOBIAN|vm::NORMAL;
const size_type jt = vm::JACOBIAN|vm::NORMAL|vm::TANGENT;
const size_type jp = vm::JACOBIAN|vm::POINT;

# /* List of applicative unary operators. */
# define VF_GD                                                          \
   BOOST_PP_TUPLE_TO_LIST(                                              \
      20,                                                               \
      (                                                                 \
       ( N       , GDN       , 0, jn, Vectorial, _M_gmc->unitNormal( q )[ c1 ] ), \
       ( Nx      , GDNx      , 0, jn, Scalar   , _M_gmc->unitNormal( q )[ 0 ] ), \
       ( Ny      , GDNy      , 1, jn, Scalar   , _M_gmc->unitNormal( q )[ 1 ] ), \
       ( Nz      , GDNz      , 2, jn, Scalar   , _M_gmc->unitNormal( q )[ 2 ] ), \
       ( T       , GDT       , 0, jn, Vectorial, _M_gmc->unitTangent( q )[ c1 ] ), \
       ( Tx      , GDTx      , 0, jn, Scalar   , _M_gmc->unitTangent( q )[ 0 ] ), \
       ( Ty      , GDTy      , 1, jn, Scalar   , _M_gmc->unitTangent( q )[ 1 ] ), \
       ( Tz      , GDTz      , 2, jn, Scalar   , _M_gmc->unitTangent( q )[ 2 ] ), \
       ( P       , GDP       , 0, jp, Vectorial, _M_gmc->xReal( q )[ c1 ]      ), \
       ( Px      , GDPx      , 0, jp, Scalar   , _M_gmc->xReal( q )[ 0 ]      ), \
       ( Py      , GDPy      , 1, jp, Scalar   , _M_gmc->xReal( q )[ 1 ]      ), \
       ( Pz      , GDPz      , 2, jp, Scalar   , _M_gmc->xReal( q )[ 2 ]      ), \
       ( C       , GDC       , 0, jp, Vectorial, _M_gmc->barycenterReal()[c1] ), \
       ( Cx      , GDCx      , 0, jp, Scalar   , _M_gmc->barycenterReal()[0]  ), \
       ( Cy      , GDCy      , 1, jp, Scalar   , _M_gmc->barycenterReal()[1]  ), \
       ( Cz      , GDCz      , 2, jp, Scalar   , _M_gmc->barycenterReal()[2]  ), \
       ( h       , GDH       , 0, 0 , Scalar   , _M_gmc->h()                  ), \
       ( hFace   , GDHFace   , 0, 0 , Scalar   , _M_gmc->hFace()              ), \
       ( eid     , GDEid     , 0, 0 , Scalar   , _M_gmc->id()                 ), \
       ( emarker , GDEmarker , 0, 0 , Scalar   , _M_gmc->marker().value()     ) \
      )                                                                          \
   )                                                                             \
/**/
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
                                                                        \
        VF_GD_NAME(O) ()                                                \
        {                                                               \
        }                                                               \
        VF_GD_NAME(O) ( VF_GD_NAME(O) const& /*__vf*/ )                 \
        {                                                               \
        }                                                               \
                                                                        \
        template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t> \
            struct tensor                                               \
        {                                                               \
            typedef this_type expression_type;                          \
            typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >, \
                mpl::identity<detail::gmc<0> >,                         \
                mpl::identity<detail::gmc<1> > >::type::type key_type;  \
            typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype; \
            typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type; \
            typedef typename gmc_type::value_type value_type;           \
            typedef VF_GD_RETURN(O)<gmc_type::NDim> return_value_type;  \
            typedef Shape<gmc_type::NDim, VF_GD_RETURN(O), false> shape; \
                                                                        \
            struct is_zero { static const bool value = false; };        \
                                                                        \
            tensor( this_type const& /*expr*/,                          \
                    Geo_t const& geom,                                  \
                    Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) \
                :                                                       \
                _M_gmc( fusion::at_key<key_type>( geom ).get() )        \
                {}                                                      \
            tensor( this_type const& /*expr*/,                          \
                    Geo_t const& geom,                                  \
                    Basis_i_t const& /*fev*/ )                          \
                :                                                       \
                _M_gmc( fusion::at_key<key_type>( geom ).get() )        \
                {}                                                      \
            tensor( this_type const& /*expr*/,                          \
                    Geo_t const& geom )                                 \
                :                                                       \
                _M_gmc( fusion::at_key<key_type>( geom ).get() )        \
                {}                                                      \
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
                _M_gmc = fusion::at_key<key_type>( geom ).get();        \
            }                                                           \
            template<typename AnyIndexI,typename AnyIndexJ>             \
                value_type                                              \
                evalijq( AnyIndexI const& /*i*/, AnyIndexJ const& /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return evalq_( c1, c2, q, mpl::bool_<(VF_GD_DIM(O) < gmc_type::NDim)>(), mpl::int_<100>() ); \
            }                                                           \
            template<typename AnyIndexI,typename AnyIndexJ, int PatternContext> \
                value_type                                              \
                evalijq( AnyIndexI const& /*i*/, AnyIndexJ const& /*j*/, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<PatternContext> ) const \
            {                                                           \
                return evalq_( c1, c2, q, mpl::bool_<(VF_GD_DIM(O) < gmc_type::NDim)>(), mpl::int_<100>() ); \
            }                                                           \
            template<typename AnyIndexI>                                \
                value_type                                              \
                evaliq( AnyIndexI const& /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
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
                Life::detail::ignore_unused_variable_warning(q);        \
                Life::detail::ignore_unused_variable_warning(c1);       \
                Life::detail::ignore_unused_variable_warning(c2);       \
                return VF_GD_VALUE( O );                                \
            }                                                           \
            gmc_ptrtype _M_gmc;                                         \
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
BOOST_PP_LIST_FOR_EACH_PRODUCT(VF_ARRAY_GD, 1, (VF_GD))
/// \endcond
} // vf
} // Life
#endif /* __GeometricData_H */
