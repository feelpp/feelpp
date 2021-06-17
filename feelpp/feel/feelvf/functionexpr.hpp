/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Thibaut Metivet <thibaut.metivet@inria.fr>
       Date: 2021-06-03

  Copyright (C) 2021 Feel++ Consortium

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
#ifndef _FUNCTION_EXPR_HPP
#define _FUNCTION_EXPR_HPP 1

#include <vector>
#include <type_traits>

#include <feel/feelvf/expr.hpp>

namespace Feel {

template< typename FunctionType, int Order, typename FunctionArgsTupleType >
class FunctionExpr: ExprDynamicBase
{
public:
    using self_type = FunctionExpr< FunctionType, Order, FunctionArgsTupleType >;
    using super_type = ExprDynamicBase;
    using function_type = FunctionType;
    using args_tuple_type = FunctionArgsTupleType;

    static constexpr auto make_type_tuple = []( auto&& ... xs ) {
        return hana::tuple<hana::type<std::decay_t<decltype(xs)>>...>{};
    };
    static constexpr auto args_types = decltype( hana::unpack( std::declval<args_tuple_type>(), make_type_tuple ) ){};

    typedef double value_type;

    // context
    template<typename ... ExprType>
    struct ExprContextProd
    {
        static constexpr size_type context = ( ... | ExprType::context );
    };
    static constexpr size_type context = decltype( hana::unpack( args_types, hana::template_<ExprContextProd> ) )::type::context;

    static constexpr bool is_terminal = false;

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
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    //--------------------------------------------------------------------//
    // Constructors
    FunctionExpr( function_type const& f, args_tuple_type const& args ) : 
        super_type( hana::fold( args, 0, []( size_type c, auto const& e ){ return c | Feel::vf::dynamicContext(e); } ) ),
        M_function( f ),
        M_argsTuple( args ),
        M_imOrder( -1 )
    {}
    FunctionExpr( function_type && f, args_tuple_type && args ) : 
        super_type( hana::fold( args, 0, []( size_type c, auto const& e ){ return c | Feel::vf::dynamicContext(e); } ) ),
        M_function( f ),
        M_argsTuple( args ),
        M_imOrder( -1 )
    {}
    FunctionExpr( FunctionExpr const& ) = default;
    FunctionExpr( FunctionExpr && ) = default;

    //--------------------------------------------------------------------//
    // Polynomial order
    uint16_type polynomialOrder() const 
    {
        if( M_imOrder >= 0 )
            return M_imOrder;
        else
            return Order;
    }
    void setPolynomialOrder( int imorder ) { M_imOrder = imorder; }
    bool isPolynomial() const { return false; }

    //--------------------------------------------------------------------//
    // Accessors
    function_type const& function() const { return M_function; }
    args_tuple_type const& argsTuple() const { return M_argsTuple; }

    //--------------------------------------------------------------------//
    // Expr tensor
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef self_type expr_type;
        using value_type = typename expr_type::value_type;
        // geomap context
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gmc_type::gm_type gm_type;
        // args tensors
        static constexpr auto arg_expr_tensor_type = []( auto argExpr ) { 
            typedef std::decay_t<decltype(argExpr)> arg_expr_type;
            typedef typename arg_expr_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> arg_tensor_type;
            return std::make_shared<arg_tensor_type>( argExpr, std::declval<Geo_t>() );
        };
        using args_tensors_tuple_type = decltype( 
                hana::transform( std::declval<args_tuple_type>(), arg_expr_tensor_type )
                );
        // args types
        static constexpr auto arg_expr_eval_type = []( auto const& argTensor ) {
            typedef Feel::decay_type<decltype(argTensor)> arg_tensor_type;
            if constexpr( arg_tensor_type::shape::is_scalar )
            {
                return argTensor->evalq( 0, 0, 0 );
            }
            else
            {
                typedef Feel::eigen_matrix_type<arg_tensor_type::shape::M, arg_tensor_type::shape::N, typename arg_tensor_type::value_type> tensor_eval_type;
                tensor_eval_type evalq;
                return evalq;
            }
        };
        using args_eval_tuple_type = decltype(
                hana::transform( std::declval<args_tensors_tuple_type>(), arg_expr_eval_type )
                );
        using function_invoke_result_type = std::invoke_result_t< decltype(hana::unpack), args_eval_tuple_type, function_type >;
        // shape
        template<typename T, bool isArithmetic>
        struct shapeType;
        template<typename T>
        struct shapeType<T, true> { 
            /* the function returns a scalar */
            using type = Shape<gmc_type::NDim, Scalar, false, false>; 
        };
        template<typename T>
        struct shapeType<T, false> { 
            /* the function returns a matrix (Eigen) */
            using type = ShapeGeneric<gmc_type::NDim, T::RowsAtCompileTime, T::ColsAtCompileTime>;
        };
        typedef typename shapeType<function_invoke_result_type, std::is_arithmetic_v<function_invoke_result_type>>::type shape;
        typedef Eigen::Matrix<value_type,shape::M,shape::N> matrix_shape_type;
        // is zero
        struct is_zero { static const bool value = false; };

        // constructors
        tensor( expr_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) :
            tensor( expr, geom )
        {}

        tensor( expr_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev ) :
            tensor( expr, geom )
        {}

        tensor( expr_type const& expr, Geo_t const& geom ) :
            M_expr( expr ),
            M_argsTensorsTuple( hana::transform( expr.argsTuple(), 
                        [&geom]( auto const& argExpr ) {
                            typedef std::decay_t<decltype(argExpr)> arg_expr_type;
                            typedef typename arg_expr_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> arg_tensor_type;
                            return std::make_shared<arg_tensor_type>( argExpr, geom );
                        } )
                    ),
            M_localMatrices( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {}

        tensor( tensor const& ) = default;

        template<typename IM>
        void init( IM const& im )
        {}

        // update
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
            hana::for_each( M_argsTensorsTuple, [&geom]( auto const& argTensor ) 
                    {
                        argTensor->update( geom );
                    } );
            this->updateImpl();
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            hana::for_each( M_argsTensorsTuple, [&geom,&face]( auto const& argTensor ) 
                    {
                        argTensor->update( geom, face );
                    } );
            this->updateImpl();
        }

        // eval
        matrix_shape_type const&
        evalijq( uint16_type, uint16_type, uint16_type q ) const
        {
            // Evaluation only
            return this->evalq( q );
        }
        value_type
        evalijq( uint16_type, uint16_type, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            // Evaluation only
            return this->evalq( c1, c2, q );
        }
        value_type
        evaliq( uint16_type, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            // Evaluation only
            return this->evalq( c1, c2, q );
        }
        matrix_shape_type const&
        evaliq( uint16_type, uint16_type q ) const
        {
            // Evaluation only
            return this->evalq( q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_localMatrices[q](c1,c2);
        }
        matrix_shape_type const&
        evalq( uint16_type q ) const
        {
            return M_localMatrices[q];
        }

    private:
        void updateImpl()
        {
            for( uint16_type q = 0; q < M_localMatrices.size(); ++q )
            {
                auto argsValuesTuple = hana::transform( M_argsTensorsTuple, [q]( auto const& argTensor ) {
                            typedef Feel::decay_type<decltype(argTensor)> arg_tensor_type;
                            if constexpr( arg_tensor_type::shape::is_scalar )
                            {
                                return argTensor->evalq( 0, 0, q );
                            }
                            else
                            {
                                typedef Feel::eigen_matrix_type<arg_tensor_type::shape::M, arg_tensor_type::shape::N, typename arg_tensor_type::value_type> tensor_eval_type;
                                tensor_eval_type evalq;
                                for( int c1 = 0; c1 < arg_tensor_type::shape::M; ++c1 )
                                    for( int c2 = 0; c2 < arg_tensor_type::shape::N; ++c2 )
                                        evalq(c1,c2) = argTensor->evalq( c1, c2, q );
                                return evalq; 
                                // TODO: improve with lazy-evaluation support (Eigen NullaryExpr)
                            }
                        } );
                if constexpr( std::is_arithmetic_v< function_invoke_result_type > )
                    M_localMatrices[q](0,0) = hana::unpack( argsValuesTuple, M_expr.function() );
                else
                    M_localMatrices[q] = hana::unpack( argsValuesTuple, M_expr.function() );
            }
        }

    private:
        expr_type const& M_expr;
        args_tensors_tuple_type M_argsTensorsTuple;
        std::vector<matrix_shape_type> M_localMatrices;
    };

private:
    function_type M_function;
    args_tuple_type M_argsTuple;
    int M_imOrder;

};

template< typename FunctionType, int Order, typename FunctionArgsTupleType >
auto functionExprWithArgsTuple( FunctionType const& f, FunctionArgsTupleType const& args )
{
    return Expr( FunctionExpr<FunctionType, Order, FunctionArgsTupleType>( f, args ) );
}

template< typename FunctionType, int Order, typename FunctionArgsTupleType >
auto functionExprWithArgsTuple( FunctionType && f, FunctionArgsTupleType && args )
{
    return Expr( FunctionExpr<FunctionType, Order, FunctionArgsTupleType>( f, args ) );
}

template<typename FunctionType, typename... FunctionArgsExprsType>
auto functionExpr( FunctionType f, const FunctionArgsExprsType & ...  argsExprs )
{
    // default Order = 2
    return functionExprWithArgsTuple<FunctionType, 2>( f, hana::make_tuple( argsExprs... ) );
}

template<int O, typename FunctionType, typename... FunctionArgsExprsType>
auto functionExpr( FunctionType f, const FunctionArgsExprsType & ...  argsExprs )
{
    return functionExprWithArgsTuple<FunctionType, O>( f, hana::make_tuple( argsExprs... ) );
}
    
} // namespace Feel

#endif //_FUNCTION_EXPR_HPP
