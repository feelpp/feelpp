/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 17 Mar 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_VF_CAST_HPP
#define FEELPP_VF_CAST_HPP 1

#include <feel/feelvf/expr.hpp>

namespace Feel
{
namespace vf
{
class CastBase {};

template < typename T, class ExprT>
class Cast : public CastBase
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

    using value_type = T;
    using expression_type = ExprT;
    using evaluate_type = Eigen::Matrix<value_type,expression_type::evaluate_type::RowsAtCompileTime, expression_type::evaluate_type::ColsAtCompileTime>;
    using this_type = Cast<T,ExprT>;

    Cast() = delete;

    constexpr explicit Cast( ExprT&& value )
        :
        M_expr( std::move( value ) )
    {
    }
    constexpr explicit Cast( ExprT const& value )
        :
        M_expr( value )
    {
    }
    Cast( Cast const& c ) = default;
    Cast& operator=( Cast const& c ) = default;

    //! polynomial order
    uint16_type polynomialOrder() const { return 0; }

    //! expression is polynomial?
    constexpr bool isPolynomial() const { return true; }

    expression_type const& expression() const { return M_expr; }
    expression_type & expression()  { return M_expr; }

    constexpr value_type evaluate( bool p, worldcomm_ptr_t const& worldcomm) const
    {
        auto eval = M_expr.evaluate(p,worldcomm);
        evaluate_type res(eval.rows(),eval.cols());
        for ( uint16_type i=0;i< eval.rows();++i )
            for ( uint16_type j=0;j< eval.cols();++j )
                res(i,j) = (value_type) eval(i,j);
        return res;
    }

    template<typename... TheExpr>
    struct Lambda
    {
        typedef expression_type type;
    };
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  ) { return typename Lambda<TheExpr...>::type(M_expr); }

    template<typename Geo_t, typename Basis_i_t=mpl::void_, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename ExprT::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        using value_type = typename Cast<T,ExprT>::value_type;

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

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_expr( expr.expression(), geom, fev, feu )
        {
        }
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr.expression(), geom, fev )
        {
        }
        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr.expression(), geom )
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
        evalij( uint16_type i, uint16_type j ) const
        {
            return (value_type)M_expr.evalij( i, j );
        }


        constexpr value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q  ) const
        {
            return (value_type)M_expr.evalijq( i, j, c1, c2, q );
        }
        template<int PatternContext>
        constexpr value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return (value_type)M_expr.evalijq( i, j, c1, c2, q );
        }

        constexpr value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q  ) const
        {
            return (value_type)M_expr.evaliq( i, c1, c2, q );
        }
        constexpr value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return (value_type)M_expr.evalq( c1, c2, q );
        }
        tensor_expr_type M_expr;
    };
protected:
    ExprT M_expr;
};

template<typename ExprT>
inline
Expr<Cast<int,ExprT>>
to_int( ExprT&& v )
{
    typedef Cast<int,ExprT> cast_t;
    return Expr< cast_t >(  cast_t( v ) );
}

template<typename T, typename ExprT>
inline
Expr<Cast<T,ExprT>>
cast( ExprT&& v )
{
    typedef Cast<T,ExprT> cast_t;
    return Expr< cast_t >(  cast_t( v ) );
}


} // vf
} // Feel

#endif /* FEELPP_VF_CAST_HPP */
