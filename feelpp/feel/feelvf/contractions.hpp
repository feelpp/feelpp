/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme
       Date: 2026-03-26

  Copyright (C) 2026 Feel++ Consortium

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
   \file contractions.hpp
   \author Christophe Prud'homme
   \date 2026-03-26
 */
#ifndef FEELPP_VF_CONTRACTIONS_HPP
#define FEELPP_VF_CONTRACTIONS_HPP 1

#include <feel/feelvf/voigt.hpp>

namespace Feel
{
namespace vf
{

namespace detail
{

template <typename T>
using normalized_scalar_expression_t =
    std::remove_cvref_t<decltype( storageExpr( std::declval<T const&>() ) )>;

template <int Dim, std::size_t... StorageIndex>
consteval auto
symmetricStoragePairsImpl( std::index_sequence<StorageIndex...> )
{
    return std::array{ symmetricPairFromStorageIndex<Dim, static_cast<int>( StorageIndex )>()... };
}

template <int Dim>
inline constexpr auto symmetric_storage_pairs_v =
    symmetricStoragePairsImpl<Dim>( std::make_index_sequence<symmetric_storage_size_v<Dim>>{} );

template <int Dim, SymmetricTensorNotation Notation, typename LambdaExprT, typename MuExprT>
class IsotropicStiffness
{
public:
    using lambda_expression_type = normalized_scalar_expression_t<LambdaExprT>;
    using mu_expression_type = normalized_scalar_expression_t<MuExprT>;

    static constexpr int storage_size = symmetric_storage_size_v<Dim>;
    static constexpr size_type context = lambda_expression_type::context | mu_expression_type::context;
    static inline const bool is_terminal = false;

    template <typename Func>
    struct HasTestFunction
    {
        static inline const bool result =
            lambda_expression_type::template HasTestFunction<Func>::result ||
            mu_expression_type::template HasTestFunction<Func>::result;
    };

    template <typename Func>
    struct HasTrialFunction
    {
        static inline const bool result =
            lambda_expression_type::template HasTrialFunction<Func>::result ||
            mu_expression_type::template HasTrialFunction<Func>::result;
    };

    template <typename Func>
    static inline const bool has_test_basis =
        lambda_expression_type::template has_test_basis<Func> ||
        mu_expression_type::template has_test_basis<Func>;

    template <typename Func>
    static inline const bool has_trial_basis =
        lambda_expression_type::template has_trial_basis<Func> ||
        mu_expression_type::template has_trial_basis<Func>;

    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;
    using value_type = typename lambda_expression_type::value_type;
    using evaluate_type = Eigen::Matrix<value_type, storage_size, storage_size>;
    using this_type = IsotropicStiffness<Dim, Notation, LambdaExprT, MuExprT>;

    explicit IsotropicStiffness( LambdaExprT const& lambda, MuExprT const& mu )
        :
        M_lambda( storageExpr( lambda ) ),
        M_mu( storageExpr( mu ) )
    {}

    constexpr uint16_type polynomialOrder() const
    {
        return std::max( M_lambda.polynomialOrder(), M_mu.polynomialOrder() );
    }

    constexpr bool isPolynomial() const
    {
        return M_lambda.isPolynomial() && M_mu.isPolynomial();
    }

    evaluate_type evaluate( bool p ) const
    {
        auto const lambdaValue = M_lambda.evaluate( p )( 0, 0 );
        auto const muValue = M_mu.evaluate( p )( 0, 0 );
        evaluate_type result;

        for ( int row = 0; row < storage_size; ++row )
            for ( int col = 0; col < storage_size; ++col )
                result( row, col ) = entryValue( row, col, lambdaValue, muValue );

        return result;
    }

    void setParameterValues( std::map<std::string,double> const& mp )
    {
        M_lambda.setParameterValues( mp );
        M_mu.setParameterValues( mp );
    }

    void updateParameterValues( std::map<std::string,double>& pv ) const
    {
        M_lambda.updateParameterValues( pv );
        M_mu.updateParameterValues( pv );
    }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto lambdaExpr = M_lambda.applySymbolsExpr( se );
        auto muExpr = M_mu.applySymbolsExpr( se );
        using result_type = IsotropicStiffness<Dim,
                                               Notation,
                                               std::decay_t<decltype( lambdaExpr )>,
                                               std::decay_t<decltype( muExpr )>>;
        return result_type( lambdaExpr, muExpr );
    }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        return M_lambda.hasSymbolDependency( symb, se ) ||
               M_mu.hasSymbolDependency( symb, se );
    }

    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>>& res,
                           TheSymbolExprType const& se ) const
    {
        M_lambda.dependentSymbols( symb, res, se );
        M_mu.dependentSymbols( symb, res, se );
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        auto lambdaDiff = M_lambda.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
        auto muDiff = M_mu.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
        using result_type = IsotropicStiffness<Dim,
                                               Notation,
                                               std::decay_t<decltype( lambdaDiff )>,
                                               std::decay_t<decltype( muDiff )>>;
        return result_type( lambdaDiff, muDiff );
    }

    lambda_expression_type const& lambdaExpr() const { return M_lambda; }
    mu_expression_type const& muExpr() const { return M_mu; }

    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using expression_type = this_type;
        using lambda_tensor_type = typename lambda_expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using mu_tensor_type = typename mu_expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using value_type = typename expression_type::value_type;
        using shape = ShapeGeneric<gmc_t<Geo_t>::nDim, storage_size, storage_size>;

        template <class Args> struct sig
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = false;
        };

        tensor( expression_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_lambda( expr.lambdaExpr(), geom, fev, feu ),
            M_mu( expr.muExpr(), geom, fev, feu )
        {}

        tensor( expression_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            M_lambda( expr.lambdaExpr(), geom, fev ),
            M_mu( expr.muExpr(), geom, fev )
        {}

        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_lambda( expr.lambdaExpr(), geom ),
            M_mu( expr.muExpr(), geom )
        {}

        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType& ttse,
                expression_type const& expr, Geo_t const& geom, TheArgsType const&... theInitArgs )
            :
            M_lambda( std::true_type{}, exprExpanded.expression().lambdaExpr(), ttse, expr.lambdaExpr(), geom, theInitArgs... ),
            M_mu( std::true_type{}, exprExpanded.expression().muExpr(), ttse, expr.muExpr(), geom, theInitArgs... )
        {}

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_lambda.update( geom, fev, feu );
            M_mu.update( geom, fev, feu );
        }

        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_lambda.update( geom, fev );
            M_mu.update( geom, fev );
        }

        void update( Geo_t const& geom )
        {
            M_lambda.update( geom );
            M_mu.update( geom );
        }

        template <typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            M_lambda.updateContext( ctx... );
            M_mu.updateContext( ctx... );
        }

        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType& ttse,
                     Geo_t const& geom, TheArgsType const&... theUpdateArgs )
        {
            M_lambda.update( std::true_type{}, exprExpanded.expression().lambdaExpr(), ttse, geom, theUpdateArgs... );
            M_mu.update( std::true_type{}, exprExpanded.expression().muExpr(), ttse, geom, theUpdateArgs... );
        }

        value_type evalij( uint16_type i, uint16_type j ) const
        {
            return entryValue( i, j, M_lambda.evalij( i, j ), M_mu.evalij( i, j ) );
        }

        value_type evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return entryValue( c1, c2, M_lambda.evalijq( i, j, 0, 0, q ), M_mu.evalijq( i, j, 0, 0, q ) );
        }

        template <int PatternContext>
        value_type evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                            mpl::int_<PatternContext> ) const
        {
            return entryValue( c1, c2,
                               M_lambda.evalijq( i, j, 0, 0, q, mpl::int_<PatternContext>{} ),
                               M_mu.evalijq( i, j, 0, 0, q, mpl::int_<PatternContext>{} ) );
        }

        value_type evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return entryValue( c1, c2, M_lambda.evaliq( i, 0, 0, q ), M_mu.evaliq( i, 0, 0, q ) );
        }

        value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return entryValue( c1, c2, M_lambda.evalq( 0, 0, q ), M_mu.evalq( 0, 0, q ) );
        }

        lambda_tensor_type M_lambda;
        mu_tensor_type M_mu;
    };

private:
    static value_type entryValue( int rowStorageIndex, int colStorageIndex, value_type lambda, value_type mu )
    {
        auto const& rowPair = symmetric_storage_pairs_v<Dim>[rowStorageIndex];
        auto const& colPair = symmetric_storage_pairs_v<Dim>[colStorageIndex];
        int const i = rowPair[0];
        int const j = rowPair[1];
        int const k = colPair[0];
        int const l = colPair[1];

        if ( i == j && k == l )
            return ( i == k ) ? ( lambda + 2.0*mu ) : lambda;

        if ( rowStorageIndex == colStorageIndex && i != j )
            return 2.0*mu;

        return value_type( 0 );
    }

    lambda_expression_type M_lambda;
    mu_expression_type M_mu;
};

template <typename T>
struct is_isotropic_stiffness_expr : std::false_type {};

template <int Dim, SymmetricTensorNotation Notation, typename LambdaExprT, typename MuExprT>
struct is_isotropic_stiffness_expr<Expr<IsotropicStiffness<Dim, Notation, LambdaExprT, MuExprT>>> : std::true_type {};

template <typename T>
inline constexpr bool is_isotropic_stiffness_expr_v = is_isotropic_stiffness_expr<std::remove_cvref_t<T>>::value;

template <typename T>
concept IsotropicStiffnessExpression = is_isotropic_stiffness_expr_v<T>;

template <typename ConstitutiveExprT, typename TensorExprT>
consteval bool
isValidSymmetricConstitutiveAction()
{
    constexpr int constitutiveRows = expression_rows_v<ConstitutiveExprT>;
    constexpr int constitutiveCols = expression_cols_v<ConstitutiveExprT>;
    constexpr int tensorDim = expression_rows_v<TensorExprT>;
    constexpr int storageSize = symmetric_storage_size_v<tensorDim>;

    return constitutiveRows == constitutiveCols &&
           constitutiveRows == storageSize;
}

template <typename ConstitutiveExprT, typename TensorExprT>
concept StaticSymmetricConstitutiveAction =
    StaticSquareMatrixExpression<ConstitutiveExprT> &&
    StaticSquareMatrixExpression<TensorExprT> &&
    ( expression_rows_v<TensorExprT> > 1 ) &&
    isValidSymmetricConstitutiveAction<ConstitutiveExprT, TensorExprT>();

template <typename ConstitutiveExprT, typename StorageExprT>
consteval bool
isValidSymmetricStorageConstitutiveAction()
{
    constexpr int constitutiveRows = expression_rows_v<ConstitutiveExprT>;
    constexpr int constitutiveCols = expression_cols_v<ConstitutiveExprT>;
    constexpr int storageSize = vector_length_v<StorageExprT>;

    return constitutiveRows == constitutiveCols &&
           constitutiveRows == storageSize;
}

template <typename ConstitutiveExprT, typename StorageExprT>
concept StaticSymmetricStorageConstitutiveAction =
    StaticSquareMatrixExpression<ConstitutiveExprT> &&
    StaticSymmetricVectorExpression<StorageExprT> &&
    isValidSymmetricStorageConstitutiveAction<ConstitutiveExprT, StorageExprT>();

template <int Dim, int RowStorageIndex, int ColStorageIndex, typename LambdaExprT, typename MuExprT>
[[nodiscard]] inline auto
isotropicStiffnessEntry( LambdaExprT const& lambda, MuExprT const& mu )
{
    constexpr int i = symmetric_storage_row_v<Dim, RowStorageIndex>;
    constexpr int j = symmetric_storage_col_v<Dim, RowStorageIndex>;
    constexpr int k = symmetric_storage_row_v<Dim, ColStorageIndex>;
    constexpr int l = symmetric_storage_col_v<Dim, ColStorageIndex>;

    if constexpr ( i == j && k == l )
    {
        if constexpr ( i == k )
            return lambda + 2.0*mu;
        else
            return lambda;
    }
    else if constexpr ( RowStorageIndex == ColStorageIndex && i != j )
    {
        return 2.0*mu;
    }
    else
    {
        return 0.0;
    }
}

template <int Dim, typename LambdaExprT, typename MuExprT, std::size_t... FlatIndices>
[[nodiscard]] inline auto
isotropicStiffnessImpl( LambdaExprT const& lambda, MuExprT const& mu,
                        std::index_sequence<FlatIndices...> )
{
    constexpr int storageSize = symmetric_storage_size_v<Dim>;
    return mat<storageSize, storageSize>(
        isotropicStiffnessEntry<Dim,
                                static_cast<int>( FlatIndices/storageSize ),
                                static_cast<int>( FlatIndices%storageSize )>( lambda, mu )... );
}

template <int Dim, typename ExprT, std::size_t... FlatIndices>
[[nodiscard]] inline auto
scaledIdentityImpl( ExprT const& expr, std::index_sequence<FlatIndices...> )
{
    return mat<Dim, Dim>(
        ([]( auto const& value )
         {
             if constexpr ( static_cast<int>( FlatIndices/Dim ) == static_cast<int>( FlatIndices%Dim ) )
                 return value;
             else
                 return cst( 0.0 );
         }( expr ))... );
}

template <int Dim, typename ExprT>
[[nodiscard]] inline auto
scaledIdentity( ExprT const& expr )
{
    return scaledIdentityImpl<Dim>( expr, std::make_index_sequence<Dim*Dim>{} );
}

template <int Dim, typename StorageExprT, std::size_t... StorageIndex>
[[nodiscard]] inline auto
symmetricStorageTraceImpl( StorageExprT const& storage, std::index_sequence<StorageIndex...> )
{
    return ( ([]( auto const& value )
               {
                   if constexpr ( symmetric_storage_row_v<Dim, static_cast<int>( StorageIndex )> ==
                                  symmetric_storage_col_v<Dim, static_cast<int>( StorageIndex )> )
                       return value;
                   else
                       return cst( 0.0 );
               }( vectorComponent<static_cast<int>( StorageIndex )>( storage ) )) + ... );
}

template <int Dim, typename StorageExprT>
[[nodiscard]] inline auto
symmetricStorageTrace( StorageExprT const& storage )
{
    return symmetricStorageTraceImpl<Dim>( storage, std::make_index_sequence<symmetric_storage_size_v<Dim>>{} );
}

template <SymmetricTensorNotation Notation, int Dim, int StorageIndex, typename TensorExprT>
[[nodiscard]] inline auto
symmetricStorageTerm( TensorExprT const& expr )
{
    if constexpr ( Notation == SymmetricTensorNotation::Mandel )
        return mandelTerm<Dim, StorageIndex>( expr );
    else
        return voigtTerm<Dim, StorageIndex>( expr );
}

template <SymmetricTensorNotation Notation, int Dim, int StorageIndex, typename ExprT>
[[nodiscard]] inline auto
applySymmetricReconstructionScale( ExprT const& expr )
{
    constexpr int i = symmetric_storage_row_v<Dim, StorageIndex>;
    constexpr int j = symmetric_storage_col_v<Dim, StorageIndex>;

    if constexpr ( Notation == SymmetricTensorNotation::Mandel && i != j )
        return cst( 1.0/std::numbers::sqrt2_v<double> )*expr;
    else
        return expr;
}

template <SymmetricTensorNotation Notation, int Dim, int StorageIndex, typename TensorExprT>
[[nodiscard]] inline auto
bilinearRightTerm( TensorExprT const& expr )
{
    auto term = symmetricStorageTerm<Notation, Dim, StorageIndex>( expr );

    if constexpr ( Notation == SymmetricTensorNotation::Voigt &&
                   symmetric_storage_row_v<Dim, StorageIndex> != symmetric_storage_col_v<Dim, StorageIndex> )
        return cst( 2.0 )*term;
    else
        return term;
}

template <SymmetricTensorNotation Notation, int Dim, int RowStorageIndex, typename ConstitutiveExprT, typename TensorExprT, std::size_t... StorageIndices>
[[nodiscard]] inline auto
constitutiveStorageRowSumImpl( ConstitutiveExprT const& C, TensorExprT const& eps, std::index_sequence<StorageIndices...> )
{
    return ( ( component<RowStorageIndex, static_cast<int>( StorageIndices )>( C )*
               symmetricStorageTerm<Notation, Dim, static_cast<int>( StorageIndices )>( eps ) ) + ... );
}

template <SymmetricTensorNotation Notation, int Dim, int Row, int Col, typename ConstitutiveExprT, typename TensorExprT>
[[nodiscard]] inline auto
constitutiveTensorEntry( ConstitutiveExprT const& C, TensorExprT const& eps )
{
    constexpr int storageIndex = symmetricStorageIndexFromComponent<Dim, Row, Col>();
    auto rowSum = constitutiveStorageRowSumImpl<Notation, Dim, storageIndex>( C, eps,
                                                                              std::make_index_sequence<symmetric_storage_size_v<Dim>>{} );
    return applySymmetricReconstructionScale<Notation, Dim, storageIndex>( rowSum );
}

template <SymmetricTensorNotation Notation, int Dim, typename ConstitutiveExprT, typename TensorExprT, std::size_t... FlatIndices>
[[nodiscard]] inline auto
constitutiveTensorActionImpl( ConstitutiveExprT const& C, TensorExprT const& eps, std::index_sequence<FlatIndices...> )
{
    return mat<Dim, Dim>(
        constitutiveTensorEntry<Notation,
                                Dim,
                                static_cast<int>( FlatIndices/Dim ),
                                static_cast<int>( FlatIndices%Dim )>( C, eps )... );
}

template <SymmetricTensorNotation Notation, int Dim, int RowStorageIndex, typename ConstitutiveExprT, typename LeftTensorExprT, typename RightTensorExprT>
[[nodiscard]] inline auto
constitutiveTensorBilinearRow( ConstitutiveExprT const& C, LeftTensorExprT const& left, RightTensorExprT const& right )
{
    return constitutiveStorageRowSumImpl<Notation, Dim, RowStorageIndex>( C, left,
                                                                          std::make_index_sequence<symmetric_storage_size_v<Dim>>{} )*
           bilinearRightTerm<Notation, Dim, RowStorageIndex>( right );
}

template <SymmetricTensorNotation Notation, int Dim, typename ConstitutiveExprT, typename LeftTensorExprT, typename RightTensorExprT, std::size_t... RowStorageIndices>
[[nodiscard]] inline auto
constitutiveTensorBilinearFormImpl( ConstitutiveExprT const& C,
                                    LeftTensorExprT const& left,
                                    RightTensorExprT const& right,
                                    std::index_sequence<RowStorageIndices...> )
{
    return ( constitutiveTensorBilinearRow<Notation, Dim, static_cast<int>( RowStorageIndices )>( C, left, right ) + ... );
}

template <SymmetricTensorNotation Notation, typename ConstitutiveExprT, typename TensorExprT>
[[nodiscard]] inline auto
constitutiveAction( ConstitutiveExprT const& C, TensorExprT const& eps )
{
    if constexpr ( Notation == SymmetricTensorNotation::Mandel )
    {
        constexpr int dim = expression_rows_v<TensorExprT>;
        return constitutiveTensorActionImpl<Notation, dim>( C, eps, std::make_index_sequence<dim*dim>{} );
    }
    else
    {
        return unvoigt( C*voigt( eps ) );
    }
}

template <SymmetricTensorNotation Notation, typename ConstitutiveExprT, typename LeftTensorExprT, typename RightTensorExprT>
[[nodiscard]] inline auto
constitutiveBilinearForm( ConstitutiveExprT const& C, LeftTensorExprT const& left, RightTensorExprT const& right )
{
    if constexpr ( Notation == SymmetricTensorNotation::Mandel )
    {
        constexpr int dim = expression_rows_v<LeftTensorExprT>;
        return constitutiveTensorBilinearFormImpl<Notation, dim>( C, left, right,
                                                                  std::make_index_sequence<symmetric_storage_size_v<dim>>{} );
    }
    else
    {
        return voigt_inner( C*voigt( left ), voigt( right ) );
    }
}

template <SymmetricTensorNotation Notation, typename ConstitutiveExprT, typename StorageExprT>
[[nodiscard]] inline auto
constitutiveStorageAction( ConstitutiveExprT const& C, StorageExprT const& eps )
{
    return C*eps;
}

template <SymmetricTensorNotation Notation, typename ConstitutiveExprT, typename LeftStorageExprT, typename RightStorageExprT>
[[nodiscard]] inline auto
constitutiveStorageBilinearForm( ConstitutiveExprT const& C, LeftStorageExprT const& left, RightStorageExprT const& right )
{
    if constexpr ( Notation == SymmetricTensorNotation::Mandel )
        return inner( C*left, right );
    else
        return voigt_inner( C*left, right );
}

template <int Dim, typename LambdaExprT, typename MuExprT, typename TensorExprT>
[[nodiscard]] inline auto
isotropicTensorAction( LambdaExprT const& lambda, MuExprT const& mu, TensorExprT const& eps )
{
    return scaledIdentity<Dim>( lambda*trace( eps ) ) + cst( 2.0 )*mu*eps;
}

template <typename LambdaExprT, typename MuExprT, typename LeftTensorExprT, typename RightTensorExprT>
[[nodiscard]] inline auto
isotropicTensorBilinearForm( LambdaExprT const& lambda, MuExprT const& mu,
                             LeftTensorExprT const& left, RightTensorExprT const& right )
{
    return lambda*trace( left )*trace( right ) + cst( 2.0 )*mu*inner( left, right );
}

template <int Dim, SymmetricTensorNotation Notation, typename LambdaExprT, typename MuExprT, typename StorageExprT>
[[nodiscard]] inline auto
isotropicStorageAction( LambdaExprT const& lambda, MuExprT const& mu, StorageExprT const& eps )
{
    auto traceExpr = symmetricStorageTrace<Dim>( eps );

    return [&]<std::size_t... StorageIndex>( std::index_sequence<StorageIndex...> )
    {
        return vec(
            ([]( auto const& traceValue,
                 auto const& lambdaValue,
                 auto const& muValue,
                 auto const& storage )
             {
                 if constexpr ( symmetric_storage_row_v<Dim, static_cast<int>( StorageIndex )> ==
                                symmetric_storage_col_v<Dim, static_cast<int>( StorageIndex )> )
                     return lambdaValue*traceValue + cst( 2.0 )*muValue*vectorComponent<static_cast<int>( StorageIndex )>( storage );
                 else
                     return cst( 2.0 )*muValue*vectorComponent<static_cast<int>( StorageIndex )>( storage );
             }( traceExpr, lambda, mu, eps ))... );
    }( std::make_index_sequence<symmetric_storage_size_v<Dim>>{} );
}

template <int Dim, SymmetricTensorNotation Notation, typename LambdaExprT, typename MuExprT, typename LeftStorageExprT, typename RightStorageExprT>
[[nodiscard]] inline auto
isotropicStorageBilinearForm( LambdaExprT const& lambda, MuExprT const& mu,
                              LeftStorageExprT const& left, RightStorageExprT const& right )
{
    auto traceLeft = symmetricStorageTrace<Dim>( left );
    auto traceRight = symmetricStorageTrace<Dim>( right );

    if constexpr ( Notation == SymmetricTensorNotation::Mandel )
        return lambda*traceLeft*traceRight + cst( 2.0 )*mu*inner( left, right );
    else
        return lambda*traceLeft*traceRight + cst( 2.0 )*mu*voigt_inner( left, right );
}

} // namespace detail

template <int Dim, SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel, typename LambdaExprT, typename MuExprT>
[[nodiscard]] inline auto
isotropic_stiffness( LambdaExprT const& lambda, MuExprT const& mu )
{
    static_assert( Dim > 0, "isotropic_stiffness() expects a strictly positive dimension" );
    using expr_t = detail::IsotropicStiffness<Dim, Notation, LambdaExprT, MuExprT>;
    return Expr<expr_t>( expr_t( lambda, mu ) );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSquareMatrixExpression TensorExprT>
[[nodiscard]] inline auto
contract( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
          TensorExprT const& eps )
{
    static_assert( Notation == ConstitutiveNotation,
                   "contract(C, eps) expects the notation to match isotropic_stiffness<Dim, Notation>" );
    static_assert( detail::expression_rows_v<TensorExprT> == Dim,
                   "contract(C, eps) expects a tensor dimension compatible with isotropic_stiffness<Dim, Notation>" );

    return detail::isotropicTensorAction<Dim>( C.expression().lambdaExpr(), C.expression().muExpr(), eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSymmetricVectorExpression StorageExprT>
[[nodiscard]] inline auto
contract( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
          StorageExprT const& eps )
{
    static_assert( Notation == ConstitutiveNotation,
                   "contract(C, eps) expects the notation to match isotropic_stiffness<Dim, Notation>" );
    static_assert( detail::vector_length_v<StorageExprT> == detail::symmetric_storage_size_v<Dim>,
                   "contract(C, eps) expects a storage vector compatible with isotropic_stiffness<Dim, Notation>" );

    return detail::isotropicStorageAction<Dim, Notation>( C.expression().lambdaExpr(), C.expression().muExpr(), eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSquareMatrixExpression LeftTensorExprT,
          detail::StaticSquareMatrixExpression RightTensorExprT>
[[nodiscard]] inline auto
contract( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
          LeftTensorExprT const& left,
          RightTensorExprT const& right )
{
    static_assert( Notation == ConstitutiveNotation,
                   "contract(C, left, right) expects the notation to match isotropic_stiffness<Dim, Notation>" );
    static_assert( detail::expression_rows_v<LeftTensorExprT> == Dim &&
                   detail::expression_rows_v<RightTensorExprT> == Dim,
                   "contract(C, left, right) expects tensors compatible with isotropic_stiffness<Dim, Notation>" );

    return detail::isotropicTensorBilinearForm( C.expression().lambdaExpr(), C.expression().muExpr(), left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSymmetricVectorExpression LeftStorageExprT,
          detail::StaticSymmetricVectorExpression RightStorageExprT>
[[nodiscard]] inline auto
contract( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
          LeftStorageExprT const& left,
          RightStorageExprT const& right )
{
    static_assert( Notation == ConstitutiveNotation,
                   "contract(C, left, right) expects the notation to match isotropic_stiffness<Dim, Notation>" );
    static_assert( detail::vector_length_v<LeftStorageExprT> == detail::symmetric_storage_size_v<Dim> &&
                   detail::vector_length_v<RightStorageExprT> == detail::symmetric_storage_size_v<Dim>,
                   "contract(C, left, right) expects storage vectors compatible with isotropic_stiffness<Dim, Notation>" );

    return detail::isotropicStorageBilinearForm<Dim, Notation>( C.expression().lambdaExpr(), C.expression().muExpr(), left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          detail::StaticSquareMatrixExpression ConstitutiveExprT,
          detail::StaticSquareMatrixExpression TensorExprT>
[[nodiscard]] inline auto
contract( ConstitutiveExprT const& C, TensorExprT const& eps )
{
    static_assert( detail::isValidSymmetricConstitutiveAction<ConstitutiveExprT, TensorExprT>(),
                   "contract(C, eps) expects a square constitutive matrix compatible with the symmetric storage of eps" );

    return detail::constitutiveAction<Notation>( C, eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          detail::StaticSquareMatrixExpression ConstitutiveExprT,
          detail::StaticSymmetricVectorExpression StorageExprT>
[[nodiscard]] inline auto
contract( ConstitutiveExprT const& C, StorageExprT const& eps )
{
    static_assert( detail::isValidSymmetricStorageConstitutiveAction<ConstitutiveExprT, StorageExprT>(),
                   "contract(C, eps) expects a square constitutive matrix compatible with the symmetric storage vector eps" );

    return detail::constitutiveStorageAction<Notation>( C, eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          detail::StaticSquareMatrixExpression ConstitutiveExprT,
          detail::StaticSquareMatrixExpression LeftTensorExprT,
          detail::StaticSquareMatrixExpression RightTensorExprT>
[[nodiscard]] inline auto
contract( ConstitutiveExprT const& C, LeftTensorExprT const& left, RightTensorExprT const& right )
{
    static_assert( detail::expression_rows_v<LeftTensorExprT> == detail::expression_rows_v<RightTensorExprT>,
                   "contract(C, left, right) expects tensors with matching dimensions" );
    static_assert( detail::isValidSymmetricConstitutiveAction<ConstitutiveExprT, LeftTensorExprT>(),
                   "contract(C, left, right) expects a square constitutive matrix compatible with the symmetric storage of left" );

    return detail::constitutiveBilinearForm<Notation>( C, left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          detail::StaticSquareMatrixExpression ConstitutiveExprT,
          detail::StaticSymmetricVectorExpression LeftStorageExprT,
          detail::StaticSymmetricVectorExpression RightStorageExprT>
[[nodiscard]] inline auto
contract( ConstitutiveExprT const& C, LeftStorageExprT const& left, RightStorageExprT const& right )
{
    static_assert( detail::vector_length_v<LeftStorageExprT> == detail::vector_length_v<RightStorageExprT>,
                   "contract(C, left, right) expects storage vectors with matching symmetric storage size" );
    static_assert( detail::isValidSymmetricStorageConstitutiveAction<ConstitutiveExprT, LeftStorageExprT>(),
                   "contract(C, left, right) expects a square constitutive matrix compatible with the storage vector left" );

    return detail::constitutiveStorageBilinearForm<Notation>( C, left, right );
}

template <typename ConstitutiveExprT, typename TensorExprT>
[[nodiscard]] inline auto
voigt_contract( ConstitutiveExprT const& C, TensorExprT const& eps )
{
    return contract<SymmetricTensorNotation::Voigt>( C, eps );
}

template <typename ConstitutiveExprT, typename LeftTensorExprT, typename RightTensorExprT>
[[nodiscard]] inline auto
voigt_contract( ConstitutiveExprT const& C, LeftTensorExprT const& left, RightTensorExprT const& right )
{
    return contract<SymmetricTensorNotation::Voigt>( C, left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSquareMatrixExpression TensorExprT>
[[nodiscard]] inline auto
ddot( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
      TensorExprT const& eps )
{
    return contract<Notation>( C, eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSymmetricVectorExpression StorageExprT>
[[nodiscard]] inline auto
ddot( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
      StorageExprT const& eps )
{
    return contract<Notation>( C, eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSquareMatrixExpression LeftTensorExprT,
          detail::StaticSquareMatrixExpression RightTensorExprT>
[[nodiscard]] inline auto
ddot( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
      LeftTensorExprT const& left,
      RightTensorExprT const& right )
{
    return contract<Notation>( C, left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSymmetricVectorExpression LeftStorageExprT,
          detail::StaticSymmetricVectorExpression RightStorageExprT>
[[nodiscard]] inline auto
ddot( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
      LeftStorageExprT const& left,
      RightStorageExprT const& right )
{
    return contract<Notation>( C, left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          typename ConstitutiveExprT,
          typename TensorExprT>
requires detail::StaticSymmetricConstitutiveAction<ConstitutiveExprT, TensorExprT>
[[nodiscard]] inline auto
ddot( ConstitutiveExprT const& C, TensorExprT const& eps )
{
    return contract<Notation>( C, eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          typename ConstitutiveExprT,
          typename StorageExprT>
requires detail::StaticSymmetricStorageConstitutiveAction<ConstitutiveExprT, StorageExprT>
[[nodiscard]] inline auto
ddot( ConstitutiveExprT const& C, StorageExprT const& eps )
{
    return contract<Notation>( C, eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          typename ConstitutiveExprT,
          typename LeftTensorExprT,
          typename RightTensorExprT>
requires detail::StaticSymmetricConstitutiveAction<ConstitutiveExprT, LeftTensorExprT> &&
         detail::StaticSquareMatrixExpression<RightTensorExprT> &&
         ( detail::expression_rows_v<LeftTensorExprT> == detail::expression_rows_v<RightTensorExprT> )
[[nodiscard]] inline auto
ddot( ConstitutiveExprT const& C, LeftTensorExprT const& left, RightTensorExprT const& right )
{
    return contract<Notation>( C, left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          typename ConstitutiveExprT,
          typename LeftStorageExprT,
          typename RightStorageExprT>
requires detail::StaticSymmetricStorageConstitutiveAction<ConstitutiveExprT, LeftStorageExprT> &&
         detail::StaticSymmetricVectorExpression<RightStorageExprT> &&
         ( detail::vector_length_v<LeftStorageExprT> == detail::vector_length_v<RightStorageExprT> )
[[nodiscard]] inline auto
ddot( ConstitutiveExprT const& C, LeftStorageExprT const& left, RightStorageExprT const& right )
{
    return contract<Notation>( C, left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          typename ConstitutiveExprT,
          typename TensorExprT>
requires detail::StaticSymmetricConstitutiveAction<ConstitutiveExprT, TensorExprT>
[[nodiscard]] inline auto
double_contract( ConstitutiveExprT const& C, TensorExprT const& eps )
{
    return ddot<Notation>( C, eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          typename ConstitutiveExprT,
          typename StorageExprT>
requires detail::StaticSymmetricStorageConstitutiveAction<ConstitutiveExprT, StorageExprT>
[[nodiscard]] inline auto
double_contract( ConstitutiveExprT const& C, StorageExprT const& eps )
{
    return ddot<Notation>( C, eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          typename ConstitutiveExprT,
          typename LeftTensorExprT,
          typename RightTensorExprT>
requires detail::StaticSymmetricConstitutiveAction<ConstitutiveExprT, LeftTensorExprT> &&
         detail::StaticSquareMatrixExpression<RightTensorExprT> &&
         ( detail::expression_rows_v<LeftTensorExprT> == detail::expression_rows_v<RightTensorExprT> )
[[nodiscard]] inline auto
double_contract( ConstitutiveExprT const& C, LeftTensorExprT const& left, RightTensorExprT const& right )
{
    return ddot<Notation>( C, left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          typename ConstitutiveExprT,
          typename LeftStorageExprT,
          typename RightStorageExprT>
requires detail::StaticSymmetricStorageConstitutiveAction<ConstitutiveExprT, LeftStorageExprT> &&
         detail::StaticSymmetricVectorExpression<RightStorageExprT> &&
         ( detail::vector_length_v<LeftStorageExprT> == detail::vector_length_v<RightStorageExprT> )
[[nodiscard]] inline auto
double_contract( ConstitutiveExprT const& C, LeftStorageExprT const& left, RightStorageExprT const& right )
{
    return ddot<Notation>( C, left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSquareMatrixExpression TensorExprT>
[[nodiscard]] inline auto
double_contract( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
                 TensorExprT const& eps )
{
    return ddot<Notation>( C, eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSymmetricVectorExpression StorageExprT>
[[nodiscard]] inline auto
double_contract( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
                 StorageExprT const& eps )
{
    return ddot<Notation>( C, eps );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSquareMatrixExpression LeftTensorExprT,
          detail::StaticSquareMatrixExpression RightTensorExprT>
[[nodiscard]] inline auto
double_contract( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
                 LeftTensorExprT const& left,
                 RightTensorExprT const& right )
{
    return ddot<Notation>( C, left, right );
}

template <SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel,
          int Dim,
          SymmetricTensorNotation ConstitutiveNotation,
          typename LambdaExprT,
          typename MuExprT,
          detail::StaticSymmetricVectorExpression LeftStorageExprT,
          detail::StaticSymmetricVectorExpression RightStorageExprT>
[[nodiscard]] inline auto
double_contract( Expr<detail::IsotropicStiffness<Dim, ConstitutiveNotation, LambdaExprT, MuExprT>> const& C,
                 LeftStorageExprT const& left,
                 RightStorageExprT const& right )
{
    return ddot<Notation>( C, left, right );
}

} // namespace vf
} // namespace Feel

#endif /* FEELPP_VF_CONTRACTIONS_HPP */
