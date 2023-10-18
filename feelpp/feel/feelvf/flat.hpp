//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 22 Aug 2017
//! @copyright 2017 Feel++ Consortium
//!

#ifndef FEELPP_FLAT_HPP
#define FEELPP_FLAT_HPP 1

namespace Feel {


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
#include <feel/feelcore/traits.hpp>
#include <feel/feelvf/unaryfunctor.hpp>


#if defined( FEELPP_HAS_QD_H ) && defined(FEELPP_HAS_MPFR)
# define VF_CHECK_ARITHMETIC_TYPE(VALUE_TYPE)                           \
   BOOST_STATIC_ASSERT( (::boost::is_arithmetic<VALUE_TYPE>::value ||    \
                         ::boost::is_same<VALUE_TYPE, std::complex<float> >::value || \
                         ::boost::is_same<VALUE_TYPE, std::complex<double> >::value || \
                         ::boost::is_same<VALUE_TYPE,mp_type>::value ||  \
                         ::boost::is_same<VALUE_TYPE,dd_real>::value ||  \
                         ::boost::is_same<VALUE_TYPE,qd_real>::value) ); \
   /**/
#elif defined( FEELPP_HAS_QD_H )
# define VF_CHECK_ARITHMETIC_TYPE(VALUE_TYPE)                           \
   BOOST_STATIC_ASSERT( (::boost::is_arithmetic<VALUE_TYPE>::value ||    \
                         ::boost::is_same<VALUE_TYPE, std::complex<float> >::value || \
                         ::boost::is_same<VALUE_TYPE, std::complex<double> >::value || \
                         ::boost::is_same<VALUE_TYPE,dd_real>::value ||  \
                         ::boost::is_same<VALUE_TYPE,qd_real>::value) ); \
   /**/
#elif defined( FEELPP_HAS_MPFR )
# define VF_CHECK_ARITHMETIC_TYPE(VALUE_TYPE)                           \
   BOOST_STATIC_ASSERT( (::boost::is_arithmetic<VALUE_TYPE>::value ||    \
                         ::boost::is_same<VALUE_TYPE, std::complex<float> >::value || \
                         ::boost::is_same<VALUE_TYPE, std::complex<double> >::value || \
                         ::boost::is_same<VALUE_TYPE,mp_type>::value) ); \
   /**/
#else
# define VF_CHECK_ARITHMETIC_TYPE(VALUE_TYPE)                           \
    BOOST_STATIC_ASSERT( ( ::boost::is_arithmetic<VALUE_TYPE>::value || \
                           ::boost::is_same<VALUE_TYPE, std::complex<float> >::value || \
                           ::boost::is_same<VALUE_TYPE, std::complex<double> >::value ) \
                         );                                             \
   /**/
#endif

namespace Feel
{
namespace vf
{

template < typename ExprT1 >
class FlatExpr
{
public:

    static const size_type context = ExprT1::context;
    static const bool is_terminal = ExprT1::is_terminal;

    static inline const uint16_type imorder = ExprT1::imorder;
    static const bool imIsPoly = ExprT1::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT1::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT1::template HasTrialFunction<Func>::result;
    };
    template<typename Func>
    static const bool has_test_basis = ExprT1::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = ExprT1::template has_trial_basis<Func>;
    using test_basis = typename ExprT1::test_basis;
    using trial_basis = typename ExprT1::trial_basis;

    typedef UnaryFunctor<typename ExprT1::value_type> super;
    typedef typename super::functordomain_type functordomain_type;
    typedef typename super::functordomain_ptrtype functordomain_ptrtype;
    typedef ExprT1 expression_1_type;
    typedef Val<ExprT1> this_type;
    typedef typename expression_1_type::value_type value_1_type;
    typedef value_1_type value_type;
    typedef value_type evaluate_type;

    VF_CHECK_ARITHMETIC_TYPE(value_1_type)

    explicit FlatExpr( expression_1_type const& __expr1  )
    :
    M_expr_1( __expr1 )
    {
        DVLOG(2) << "Val::Val default constructorn";
    }

    FlatExpr( FlatExpr const& __vfp  )
        :
        M_expr_1( __vfp.M_expr_1 )
    {
        DVLOG(2) << "Val::Val copy constructorn";
    }

    bool isSymetric() const
    {
        return false;
    }

    expression_1_type const& expression() const
    {
        return M_expr_1;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t,Basis_j_t> tensor2_expr_type;
        typedef typename tensor2_expr_type::value_type value_type;
        using key_type = key_t<Geo_t>;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef typename tensor2_expr_type::shape expr_shape;
        using shape = Shape<expr_shape::M*expr_shape::N,Vectorial,false>;
        struct is_zero
        {
            static const bool value = tensor2_expr_type::is_zero::value;
        };

        template<typename ExprT>
        tensor( ExprT const& expr, Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
            :
            M_expr( expr.expression(), geom ),
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
            update( geom );
        }
        template<typename ExprT>
        tensor( ExprT const& expr,Geo_t const& geom, Basis_i_t const& /*fev*/ )
            :
            M_expr( expr.expression(), geom ),
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
            update( geom );
        }
        template<typename ExprT>
        tensor( ExprT const& expr, Geo_t const& geom )
            :
            M_expr( expr.expression(), geom ),
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
            update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_expr.update( geom );

        }

        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return evalq( c1, c2, q );
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q, mpl::int_<shape::rank>() );
        }
    private:
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return M_expr.evalq( c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1> ) const
        {
            return M_expr.evalq( c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> ) const
        {
            return M_expr.evalq( c1, c2, q );
        }
    private:
        tensor2_expr_type M_expr;
        gmc_ptrtype M_gmc;
    };

protected:
    FlatExpr() {}

    expression_1_type M_expr_1;
};
/// \endcond

/**
 * \brief precompute expression tensor
 *
 * This allows for more efficient  bi/linear form assembly
 */
template<typename ExprT1>
inline
Expr<FlatExpr<ExprT1>>
flat( ExprT1 const& __e1 )
{
    typedef FlatExpr<ExprT1> expr_t;
    return Expr< expr_t >(  expr_t( t1( __e1 ) ) );
}
} // vf
} // Feel
#endif /* FEELPP_VF_VAL_HPP */


}
#endif
