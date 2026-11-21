/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <feel/feelvf/detail/localform_tensor.hpp>
#include <feel/feelvf/operations.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{

/** @brief Whether a contraction consumes a tensor-valued coefficient. */
template<typename ContractionType>
inline constexpr bool localform_accepts_tensor_coefficient_v = []
{
    if constexpr ( requires { ContractionType::accepts_tensor_coefficient; } )
        return ContractionType::accepts_tensor_coefficient;
    else
        return false;
}();

/**
 * @brief Lowered scalar bilinear expression with a contraction strategy.
 *
 * Mass and H1 stiffness share coefficient evaluation, FE-context binding,
 * quadrature batching, workspace lifetime, and assembler integration. The
 * contraction strategy supplies only the mathematical basis operator.
 *
 * @tparam TestSpaceType scalar test function space
 * @tparam TrialSpaceType scalar trial function space
 * @tparam CoefficientExprType normalized scalar coefficient IR
 * @tparam ContractionType value-value or gradient-gradient strategy
 */
template<typename TestSpaceType, typename TrialSpaceType,
         typename CoefficientExprType, typename ContractionType>
class LoweredScalarBilinearExpr
{
public:
    using coefficient_expr_type = remove_cvref_t<CoefficientExprType>;
    using contraction_type = ContractionType;
    using test_basis = typename TestSpaceType::reference_element_type;
    using trial_basis = typename TrialSpaceType::reference_element_type;
    static constexpr bool accepts_tensor_coefficient =
        localform_accepts_tensor_coefficient_v<contraction_type>;
    static_assert( is_localform_coeff_ir_v<coefficient_expr_type> ||
                   accepts_tensor_coefficient,
                   "LoweredScalarBilinearExpr expects a normalized scalar coefficient IR or a tensor-aware contraction" );
    static constexpr size_type basis_mapping_context =
        ( is_hdiv_conforming_v<test_basis> ||
          is_hdiv_conforming_v<trial_basis> )
            ? vm::JACOBIAN | vm::KB
            : ( is_hcurl_conforming_v<test_basis> ||
                is_hcurl_conforming_v<trial_basis> )
                  ? vm::KB
                  : 0;
    static const size_type context = coefficient_expr_type::context |
                                     contraction_type::context |
                                     basis_mapping_context;
    static inline const bool is_terminal = false;

    using this_type = LoweredScalarBilinearExpr<TestSpaceType, TrialSpaceType,
                                                CoefficientExprType,
                                                ContractionType>;
    using expression_type = this_type;
    using value_type = typename coefficient_expr_type::value_type;
    using evaluate_type = typename coefficient_expr_type::evaluate_type;
    static constexpr int test_workspace_rows =
        test_basis::is_order_dynamic
            ? Eigen::Dynamic
            : finiteElementFixedLocalDofCountOrZero<test_basis>( 1 );
    static constexpr int trial_workspace_rows =
        trial_basis::is_order_dynamic
            ? Eigen::Dynamic
            : finiteElementFixedLocalDofCountOrZero<trial_basis>( 1 );
    using test_workspace_type =
        Eigen::Matrix<value_type, test_workspace_rows, Eigen::Dynamic,
                      Eigen::RowMajor>;
    using trial_workspace_type =
        Eigen::Matrix<value_type, trial_workspace_rows, Eigen::Dynamic,
                      Eigen::RowMajor>;

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

    /** @brief Construct from a normalized coefficient and quadrature order. */
    LoweredScalarBilinearExpr( coefficient_expr_type coefficientExpr,
                               uint16_type polynomialOrder )
        :
        M_coefficientExpr( std::move( coefficientExpr ) ),
        M_polynomialOrder( polynomialOrder )
    {}

    size_type dynamicContext() const { return Feel::vf::dynamicContext( M_coefficientExpr ); }
    uint16_type polynomialOrder() const { return M_polynomialOrder; }
    bool isPolynomial() const { return M_coefficientExpr.isPolynomial(); }
    /**
     * @return whether symmetry follows from the lowered representation
     *
     * A tensor coefficient is conservatively reported as nonsymmetric because
     * its pointwise symmetry cannot be established from the expression type.
     */
    constexpr bool isSymmetric() const
    {
        return std::is_same_v<TestSpaceType, TrialSpaceType> &&
               !accepts_tensor_coefficient;
    }
    coefficient_expr_type const& coefficientExpr() const { return M_coefficientExpr; }
    evaluate_type evaluate( bool p ) const { return M_coefficientExpr.evaluate( p ); }

    template<typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto newExpr = M_coefficientExpr.applySymbolsExpr( se );
        using new_expr_type = remove_cvref_t<decltype( newExpr )>;
        return LoweredScalarBilinearExpr<TestSpaceType, TrialSpaceType,
                                         new_expr_type, ContractionType>(
            newExpr, M_polynomialOrder );
    }

    template<typename Geo_t, typename Basis_i_t = mpl::void_, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using expression_type = this_type;
        using coefficient_tensor_type = typename coefficient_expr_type::template tensor<Geo_t>;
        using value_type = typename expression_type::value_type;
        using gmc_type = gmc_t<Geo_t>;
        using shape = Shape<gmc_type::nDim, Scalar, false, false>;
        using tensor_coefficient_workspace_type =
            Eigen::Matrix<value_type, gmc_type::nDim*gmc_type::nDim,
                          Eigen::Dynamic, Eigen::RowMajor>;

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
            M_coefficientBatchValid = false;
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_coefficientExpr.update( geom );
            M_testFec = fev;
            M_coefficientBatchValid = false;
        }
        void update( Geo_t const& geom )
        {
            M_coefficientExpr.update( geom );
            M_coefficientBatchValid = false;
        }
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
        {
            M_coefficientExpr.updateContext( ctx... );
            M_coefficientBatchValid = false;
        }

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
            if constexpr ( accepts_tensor_coefficient )
                return contraction_type::point(
                    *testFec, *trialFec, M_coefficientExpr, i, j, q );
            else
                return M_coefficientExpr.evalq( 0, 0, q ) *
                       contraction_type::point( *testFec, *trialFec, i, j, q );
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

        /**
         * @brief Evaluate one lowered scalar tensor entry over all quadrature points.
         *
         * Coefficients are evaluated once per updated element and reused for
         * every local `(i,j)` pair. Basis values come from the lazy
         * component-major FE tables and the product is evaluated by Eigen.
         *
         * @param request value, test, or bilinear local indices
         * @param values caller-owned point-block output
         */
        void evaluateQuadratureBatch(
            quadrature::evaluation_request const& request,
            Eigen::Ref<Eigen::Matrix<value_type, Eigen::Dynamic, 1>> values ) const
            requires ( !std::is_same_v<Basis_i_t, mpl::void_> )
        {
            auto const& testFec = fusion::at_key<vf::detail::gmc<0> >( M_testFec );
            CHECK_EQ( values.size(), testFec->nPoints() );
            this->ensureCoefficientBatch( values.size() );

            switch ( request.kind )
            {
            case quadrature::evaluation_kind::value:
                values = M_coefficientBatch;
                break;
            case quadrature::evaluation_kind::test:
                values.array() = M_coefficientBatch.array() *
                                 testFec->idBatch( request.test ).transpose().array();
                break;
            case quadrature::evaluation_kind::bilinear:
            {
                auto const& trialFec =
                    fusion::at_key<vf::detail::gmc<0> >( M_trialFec );
                CHECK_EQ( values.size(), trialFec->nPoints() );
                contraction_type::pointBlock(
                    *testFec, *trialFec, request.test, request.trial,
                    M_coefficientBatch, values );
                break;
            }
            }
        }

        /**
         * @brief Form the complete local bilinear tensor with Eigen.
         *
         * This is the backend-neutral assembler hook. Physical quadrature
         * weights are combined with the normalized coefficient once and the
         * selected mathematical strategy performs the local contraction.
         *
         * @tparam WeightVector Eigen-compatible quadrature-weight vector
         * @tparam LocalMatrix assembler-owned fixed or dynamic Eigen matrix
         * @param weights physical quadrature weights for the current element
         * @param output complete local matrix before orientation transforms
         */
        template<typename WeightVector, typename LocalMatrix>
        void evaluateLocalBilinearTensor( WeightVector const& weights,
                                          LocalMatrix& output ) const
            requires ( !std::is_same_v<Basis_i_t, mpl::void_> &&
                       !std::is_same_v<Basis_j_t, mpl::void_> )
        {
            auto const& testFec =
                fusion::at_key<vf::detail::gmc<0> >( M_testFec );
            auto const& trialFec =
                fusion::at_key<vf::detail::gmc<0> >( M_trialFec );
            CHECK_EQ( weights.size(), testFec->nPoints() );
            CHECK_EQ( weights.size(), trialFec->nPoints() );
            if constexpr ( accepts_tensor_coefficient )
            {
                contraction_type::contractTensor(
                    *testFec, *trialFec, M_coefficientExpr, weights, output,
                    M_testComponent, M_weightedTrial,
                    M_tensorCoefficientBatch );
            }
            else
            {
                this->ensureCoefficientBatch( weights.size() );
                resizeEigenVectorWorkspace(
                    M_weightedCoefficient, weights.size() );
                M_weightedCoefficient.array() =
                    M_coefficientBatch.array()*weights.array();
                contraction_type::contract(
                    *testFec, *trialFec, M_weightedCoefficient, output,
                    M_testComponent, M_weightedTrial );
            }
        }

    private:
        /** @brief Build or reuse the coefficient point block for this element. */
        void ensureCoefficientBatch( Eigen::Index pointCount ) const
        {
            if ( M_coefficientBatchValid &&
                 M_coefficientBatch.size() == pointCount )
                return;
            resizeEigenVectorWorkspace( M_coefficientBatch, pointCount );
            quadrature::evaluateBatchOrFallback(
                M_coefficientExpr,
                quadrature::evaluation_request::value( 0, 0 ),
                M_coefficientBatch );
            M_coefficientBatchValid = true;
        }

    public:

        coefficient_tensor_type M_coefficientExpr;
        Basis_i_t M_testFec;
        Basis_j_t M_trialFec;
        mutable Eigen::Matrix<value_type, Eigen::Dynamic, 1> M_coefficientBatch;
        mutable Eigen::Matrix<value_type, Eigen::Dynamic, 1> M_weightedCoefficient;
        mutable test_workspace_type M_testComponent;
        mutable trial_workspace_type M_weightedTrial;
        mutable tensor_coefficient_workspace_type M_tensorCoefficientBatch;
        mutable bool M_coefficientBatchValid = false;
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
            M_coefficientBatchValid = false;
        }
        void update( Geo_t const& geom )
        {
            M_coefficientExpr.update( geom );
            M_coefficientBatchValid = false;
        }
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
        {
            M_coefficientExpr.updateContext( ctx... );
            M_coefficientBatchValid = false;
        }

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

        /**
         * @brief Evaluate the lowered scalar source tensor over a point block.
         * @param request value or test-function local indices
         * @param values caller-owned point-block output
         */
        void evaluateQuadratureBatch(
            quadrature::evaluation_request const& request,
            Eigen::Ref<Eigen::Matrix<value_type, Eigen::Dynamic, 1>> values ) const
            requires ( !std::is_same_v<Basis_i_t, mpl::void_> )
        {
            auto const& testFec = fusion::at_key<vf::detail::gmc<0> >( M_testFec );
            CHECK_EQ( values.size(), testFec->nPoints() );
            this->ensureCoefficientBatch( values.size() );
            if ( request.kind == quadrature::evaluation_kind::value )
                values = M_coefficientBatch;
            else
                values.array() = M_coefficientBatch.array() *
                                 testFec->idBatch( request.test ).transpose().array();
        }

        /**
         * @brief Form the complete scalar source vector with Eigen.
         * @tparam WeightVector Eigen-compatible quadrature weights
         * @tparam LocalVector assembler-owned fixed or dynamic Eigen vector
         * @param weights physical quadrature weights for the current element
         * @param output complete local source vector
         */
        template<typename WeightVector, typename LocalVector>
        void evaluateLocalLinearTensor( WeightVector const& weights,
                                        LocalVector& output ) const
            requires ( !std::is_same_v<Basis_i_t, mpl::void_> )
        {
            auto const& testFec =
                fusion::at_key<vf::detail::gmc<0> >( M_testFec );
            CHECK_EQ( weights.size(), testFec->nPoints() );
            CHECK_EQ( output.size(), testFec->nDofs() );
            this->ensureCoefficientBatch( weights.size() );
            resizeEigenVectorWorkspace(
                M_weightedCoefficient, weights.size() );
            M_weightedCoefficient.array() =
                M_coefficientBatch.array()*weights.array();
            output.noalias() =
                testFec->valueBatchTable()*M_weightedCoefficient;
        }

    private:
        /** @brief Build or reuse the coefficient point block for this element. */
        void ensureCoefficientBatch( Eigen::Index pointCount ) const
        {
            if ( M_coefficientBatchValid &&
                 M_coefficientBatch.size() == pointCount )
                return;
            resizeEigenVectorWorkspace( M_coefficientBatch, pointCount );
            quadrature::evaluateBatchOrFallback(
                M_coefficientExpr,
                quadrature::evaluation_request::value( 0, 0 ),
                M_coefficientBatch );
            M_coefficientBatchValid = true;
        }

    public:

        coefficient_tensor_type M_coefficientExpr;
        Basis_i_t M_testFec;
        mutable Eigen::Matrix<value_type, Eigen::Dynamic, 1> M_coefficientBatch;
        mutable Eigen::Matrix<value_type, Eigen::Dynamic, 1> M_weightedCoefficient;
        mutable bool M_coefficientBatchValid = false;
    };

private:
    coefficient_expr_type M_coefficientExpr;
    uint16_type M_polynomialOrder;
};

} // namespace detail
} // namespace vf
} // namespace Feel

#endif
