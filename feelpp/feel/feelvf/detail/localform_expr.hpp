/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Copyright (C) 2026 Feel++ Consortium

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
   \file localform_expr.hpp
   \brief Lowered local-form expression wrappers.
 */
#ifndef FEELPP_VF_DETAIL_LOCALFORM_EXPR_HPP
#define FEELPP_VF_DETAIL_LOCALFORM_EXPR_HPP 1

#include <boost/fusion/include/at_key.hpp>

#include <feel/feelvf/detail/localform_base.hpp>
#include <feel/feelvf/detail/localform_coeff_ir.hpp>
#include <feel/feelvf/operations.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{

template<typename TestSpaceType, typename TrialSpaceType, typename CoefficientExprType>
class LoweredScalarMassExpr
{
public:
    using coefficient_expr_type = remove_cvref_t<CoefficientExprType>;
    static_assert( is_localform_coeff_ir_v<coefficient_expr_type>,
                   "LoweredScalarMassExpr expects a normalized localform coefficient IR type" );
    static const size_type context = coefficient_expr_type::context | vm::JACOBIAN;
    static inline const bool is_terminal = false;

    using this_type = LoweredScalarMassExpr<TestSpaceType, TrialSpaceType, CoefficientExprType>;
    using expression_type = this_type;
    using value_type = typename coefficient_expr_type::value_type;
    using evaluate_type = typename coefficient_expr_type::evaluate_type;
    using test_basis = typename TestSpaceType::reference_element_type;
    using trial_basis = typename TrialSpaceType::reference_element_type;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = std::is_same_v<Func, test_basis>;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = std::is_same_v<Func, trial_basis>;
    };
    template<typename Func>
    static inline const bool has_test_basis = HasTestFunction<Func>::result;
    template<typename Func>
    static inline const bool has_trial_basis = HasTrialFunction<Func>::result;

    LoweredScalarMassExpr( coefficient_expr_type coefficientExpr, uint16_type polynomialOrder )
        :
        M_coefficientExpr( std::move( coefficientExpr ) ),
        M_polynomialOrder( polynomialOrder )
    {}

    size_type dynamicContext() const { return Feel::vf::dynamicContext( M_coefficientExpr ); }
    uint16_type polynomialOrder() const { return M_polynomialOrder; }
    bool isPolynomial() const { return M_coefficientExpr.isPolynomial(); }
    constexpr bool isSymmetric() const { return std::is_same_v<TestSpaceType, TrialSpaceType>; }
    coefficient_expr_type const& coefficientExpr() const { return M_coefficientExpr; }
    evaluate_type evaluate( bool p ) const { return M_coefficientExpr.evaluate( p ); }

    template<typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto newExpr = M_coefficientExpr.applySymbolsExpr( se );
        using new_expr_type = remove_cvref_t<decltype( newExpr )>;
        return LoweredScalarMassExpr<TestSpaceType, TrialSpaceType, new_expr_type>( newExpr, M_polynomialOrder );
    }

    template<typename Geo_t, typename Basis_i_t = mpl::void_, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using expression_type = this_type;
        using coefficient_tensor_type = typename coefficient_expr_type::template tensor<Geo_t>;
        using value_type = typename expression_type::value_type;
        using gmc_type = gmc_t<Geo_t>;
        using shape = Shape<gmc_type::nDim, Scalar, false, false>;

        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_coefficientExpr( expr.coefficientExpr(), geom ),
            M_testFec( fev ),
            M_trialFec( feu )
        {}
        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_coefficientExpr( expr.coefficientExpr(), geom ),
            M_testFec( fev )
        {}
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_coefficientExpr( expr.coefficientExpr(), geom )
        {}

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_coefficientExpr.update( geom );
            M_testFec = fev;
            M_trialFec = feu;
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_coefficientExpr.update( geom );
            M_testFec = fev;
        }
        void update( Geo_t const& geom ) { M_coefficientExpr.update( geom ); }
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx ) { M_coefficientExpr.updateContext( ctx... ); }

        value_type evalij( uint16_type /*i*/, uint16_type /*j*/ ) const
        {
            return value_type( 0 );
        }
        value_type evalijq( uint16_type i, uint16_type j,
                            uint16_type /*c1*/, uint16_type /*c2*/,
                            uint16_type q ) const
        {
            auto const& testFec = fusion::at_key<vf::detail::gmc<0> >( M_testFec );
            auto const& trialFec = fusion::at_key<vf::detail::gmc<0> >( M_trialFec );
            return M_coefficientExpr.evalq( 0, 0, q ) * testFec->id( i, 0, 0, q ) * trialFec->id( j, 0, 0, q );
        }
        template<int PatternContext>
        value_type evalijq( uint16_type i, uint16_type j,
                            uint16_type c1, uint16_type c2, uint16_type q,
                            mpl::int_<PatternContext> ) const
        {
            return evalijq( i, j, c1, c2, q );
        }
        value_type evaliq( uint16_type i,
                           uint16_type /*c1*/, uint16_type /*c2*/,
                           uint16_type q ) const
        {
            auto const& testFec = fusion::at_key<vf::detail::gmc<0> >( M_testFec );
            return M_coefficientExpr.evalq( 0, 0, q ) * testFec->id( i, 0, 0, q );
        }
        value_type evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_coefficientExpr.evalq( 0, 0, q );
        }

        coefficient_tensor_type M_coefficientExpr;
        Basis_i_t M_testFec;
        Basis_j_t M_trialFec;
    };

private:
    coefficient_expr_type M_coefficientExpr;
    uint16_type M_polynomialOrder;
};

template<typename SpaceType, typename CoefficientExprType>
class LoweredScalarSourceExpr
{
public:
    using coefficient_expr_type = remove_cvref_t<CoefficientExprType>;
    static_assert( is_localform_coeff_ir_v<coefficient_expr_type>,
                   "LoweredScalarSourceExpr expects a normalized localform coefficient IR type" );
    static const size_type context = coefficient_expr_type::context | vm::JACOBIAN;
    static inline const bool is_terminal = false;

    using this_type = LoweredScalarSourceExpr<SpaceType, CoefficientExprType>;
    using expression_type = this_type;
    using value_type = typename coefficient_expr_type::value_type;
    using evaluate_type = typename coefficient_expr_type::evaluate_type;
    using test_basis = typename SpaceType::reference_element_type;
    using trial_basis = std::nullptr_t;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = std::is_same_v<Func, test_basis>;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    static inline const bool has_test_basis = HasTestFunction<Func>::result;
    template<typename Func>
    static inline const bool has_trial_basis = false;

    LoweredScalarSourceExpr( coefficient_expr_type coefficientExpr, uint16_type polynomialOrder )
        :
        M_coefficientExpr( std::move( coefficientExpr ) ),
        M_polynomialOrder( polynomialOrder )
    {}

    size_type dynamicContext() const { return Feel::vf::dynamicContext( M_coefficientExpr ); }
    uint16_type polynomialOrder() const { return M_polynomialOrder; }
    bool isPolynomial() const { return M_coefficientExpr.isPolynomial(); }
    coefficient_expr_type const& coefficientExpr() const { return M_coefficientExpr; }
    evaluate_type evaluate( bool p ) const { return M_coefficientExpr.evaluate( p ); }

    template<typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto newExpr = M_coefficientExpr.applySymbolsExpr( se );
        using new_expr_type = remove_cvref_t<decltype( newExpr )>;
        return LoweredScalarSourceExpr<SpaceType, new_expr_type>( newExpr, M_polynomialOrder );
    }

    template<typename Geo_t, typename Basis_i_t = mpl::void_, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using expression_type = this_type;
        using coefficient_tensor_type = typename coefficient_expr_type::template tensor<Geo_t>;
        using value_type = typename expression_type::value_type;
        using gmc_type = gmc_t<Geo_t>;
        using shape = Shape<gmc_type::nDim, Scalar, false, false>;

        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_coefficientExpr( expr.coefficientExpr(), geom ),
            M_testFec( fev )
        {}
        tensor( expression_type const& expr, Geo_t const& geom )
            :
            M_coefficientExpr( expr.coefficientExpr(), geom )
        {}

        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_coefficientExpr.update( geom );
            M_testFec = fev;
        }
        void update( Geo_t const& geom ) { M_coefficientExpr.update( geom ); }
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx ) { M_coefficientExpr.updateContext( ctx... ); }

        value_type evalij( uint16_type /*i*/, uint16_type /*j*/ ) const
        {
            return value_type( 0 );
        }
        value_type evalijq( uint16_type i, uint16_type /*j*/,
                            uint16_type /*c1*/, uint16_type /*c2*/,
                            uint16_type q ) const
        {
            return evaliq( i, 0, 0, q );
        }
        template<int PatternContext>
        value_type evalijq( uint16_type i, uint16_type j,
                            uint16_type c1, uint16_type c2, uint16_type q,
                            mpl::int_<PatternContext> ) const
        {
            return evalijq( i, j, c1, c2, q );
        }
        value_type evaliq( uint16_type i,
                           uint16_type /*c1*/, uint16_type /*c2*/,
                           uint16_type q ) const
        {
            auto const& testFec = fusion::at_key<vf::detail::gmc<0> >( M_testFec );
            return M_coefficientExpr.evalq( 0, 0, q ) * testFec->id( i, 0, 0, q );
        }
        value_type evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_coefficientExpr.evalq( 0, 0, q );
        }

        coefficient_tensor_type M_coefficientExpr;
        Basis_i_t M_testFec;
    };

private:
    coefficient_expr_type M_coefficientExpr;
    uint16_type M_polynomialOrder;
};

} // namespace detail
} // namespace vf
} // namespace Feel

#endif
