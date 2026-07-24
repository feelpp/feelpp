/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#ifndef FEELPP_VF_DETAIL_LOCALFORM_TENSOR_HPP
#define FEELPP_VF_DETAIL_LOCALFORM_TENSOR_HPP 1

#include <concepts>

#include <Eigen/Core>

#include <feel/feelcore/context.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{

/**
 * @brief Tensor capable of forming a complete local bilinear-form matrix.
 *
 * The tensor consumes the quadrature rule's current physical integration
 * weights and writes into the assembler-owned local matrix. Implementations
 * may reuse private Eigen workspaces but must not retain either argument.
 *
 * @tparam Tensor updated variational-expression tensor
 * @tparam WeightVector Eigen column vector containing quadrature weights times
 * the geometry integration Jacobian
 * @tparam LocalMatrix assembler-owned Eigen local matrix
 */
template<typename Tensor, typename WeightVector, typename LocalMatrix>
concept LocalBilinearTensorEvaluable =
    requires( Tensor const& tensor,
              WeightVector const& weights,
              LocalMatrix& matrix )
    {
        { tensor.evaluateLocalBilinearTensor( weights, matrix ) } ->
            std::same_as<void>;
    };

/**
 * @brief Tensor capable of forming a complete local linear-form vector.
 * @tparam Tensor updated variational-expression tensor
 * @tparam WeightVector physical quadrature-weight vector
 * @tparam LocalVector assembler-owned Eigen local vector
 */
template<typename Tensor, typename WeightVector, typename LocalVector>
concept LocalLinearTensorEvaluable =
    requires( Tensor const& tensor,
              WeightVector const& weights,
              LocalVector& vector )
    {
        { tensor.evaluateLocalLinearTensor( weights, vector ) } ->
            std::same_as<void>;
    };

/** @brief Resize a reusable Eigen matrix only when its shape changes. */
template<typename Matrix>
void resizeEigenWorkspace( Matrix& matrix, Eigen::Index rows,
                           Eigen::Index columns )
{
    if ( matrix.rows() != rows || matrix.cols() != columns )
        matrix.resize( rows, columns );
}

/** @brief Resize a reusable Eigen vector only when its length changes. */
template<typename Vector>
void resizeEigenVectorWorkspace( Vector& vector, Eigen::Index size )
{
    if ( vector.size() != size )
        vector.resize( size );
}

/**
 * @brief Accumulate a component-wise local matrix contraction.
 *
 * The accessors return non-owning DoF-by-quadrature component tables. This
 * helper concatenates those views in reusable Eigen workspaces and issues one
 * matrix product for the complete value, gradient, curl, divergence, or mixed
 * contraction, without a per-DoF row traversal.
 */
template<typename TestContext, typename TrialContext, typename Vector,
         typename OutputMatrix, typename ScratchMatrix,
         typename TestTableAccessor, typename TrialTableAccessor>
void contractComponentTables( TestContext const& test,
                              TrialContext const& trial,
                              uint16_type componentCount,
                              Vector const& weightedCoefficient,
                              OutputMatrix& output,
                              ScratchMatrix& testComponents,
                              ScratchMatrix& weightedTrial,
                              TestTableAccessor testTable,
                              TrialTableAccessor trialTable )
{
    auto const testDofs = static_cast<Eigen::Index>( test.nDofs() );
    auto const trialDofs = static_cast<Eigen::Index>( trial.nDofs() );
    auto const points = weightedCoefficient.size();
    CHECK_EQ( output.rows(), testDofs );
    CHECK_EQ( output.cols(), trialDofs );

    auto const blockPoints = static_cast<Eigen::Index>( componentCount )*points;
    resizeEigenWorkspace( testComponents, testDofs, blockPoints );
    resizeEigenWorkspace( weightedTrial, trialDofs, blockPoints );
    for ( uint16_type component = 0; component < componentCount; ++component )
    {
        auto const testValues = testTable( component );
        auto const trialValues = trialTable( component );
        CHECK_EQ( testValues.rows(), testDofs );
        CHECK_EQ( trialValues.rows(), trialDofs );
        CHECK_EQ( testValues.cols(), points );
        CHECK_EQ( trialValues.cols(), points );
        auto const firstPoint = static_cast<Eigen::Index>( component )*points;
        testComponents.middleCols( firstPoint, points ) = testValues;
        weightedTrial.middleCols( firstPoint, points ).array() =
            trialValues.array().rowwise()*
            weightedCoefficient.transpose().array();
    }
    output.noalias() = testComponents*weightedTrial.transpose();
}

/** @brief Strategy tag for scalar value-value mass contractions. */
struct scalar_mass_contraction
{
    static constexpr size_type context = vm::JACOBIAN;

    /** @return one scalar mass integrand value. */
    template<typename TestContext, typename TrialContext>
    static auto point( TestContext const& test, TrialContext const& trial,
                       uint16_type i, uint16_type j, uint32_type q )
    {
        return test.id( i, 0, 0, q )*trial.id( j, 0, 0, q );
    }

    /** @brief Evaluate one local mass entry over a quadrature block. */
    template<typename TestContext, typename TrialContext, typename Vector>
    static void pointBlock( TestContext const& test, TrialContext const& trial,
                            uint16_type i, uint16_type j,
                            Vector const& coefficient, Vector& values )
    {
        values.array() = coefficient.array() *
                         test.idBatch( i ).transpose().array() *
                         trial.idBatch( j ).transpose().array();
    }

    /** @brief Form the complete scalar mass matrix with one Eigen product. */
    template<typename TestContext, typename TrialContext, typename Vector,
             typename OutputMatrix, typename ScratchMatrix>
    static void contract( TestContext const& test, TrialContext const& trial,
                          Vector const& weightedCoefficient,
                          OutputMatrix& output,
                          ScratchMatrix& /*testComponent*/,
                          ScratchMatrix& weightedTrial )
    {
        auto const& testValues = test.valueBatchTable();
        auto const& trialValues = trial.valueBatchTable();
        CHECK_EQ( output.rows(), testValues.rows() );
        CHECK_EQ( output.cols(), trialValues.rows() );
        CHECK_EQ( weightedCoefficient.size(), testValues.cols() );
        CHECK_EQ( weightedCoefficient.size(), trialValues.cols() );

        resizeEigenWorkspace(
            weightedTrial, trialValues.rows(), trialValues.cols() );
        weightedTrial.array() =
            trialValues.array().rowwise()*weightedCoefficient.transpose().array();
        output.noalias() = testValues*weightedTrial.transpose();
    }
};

/** @brief Strategy tag for vector value inner-product mass contractions. */
struct vector_mass_contraction
{
    static constexpr size_type context = vm::JACOBIAN;

    /** @return one vector basis inner product. */
    template<typename TestContext, typename TrialContext>
    static auto point( TestContext const& test, TrialContext const& trial,
                       uint16_type i, uint16_type j, uint32_type q )
    {
        static_assert( TestContext::nComponents1 == TrialContext::nComponents1 );
        using value_type = typename TestContext::value_type;
        value_type result = 0;
        for ( uint16_type component = 0;
              component < TestContext::nComponents1; ++component )
            result += test.id( i, component, 0, q )*
                      trial.id( j, component, 0, q );
        return result;
    }

    /** @brief Evaluate one vector mass entry over a quadrature block. */
    template<typename TestContext, typename TrialContext, typename Vector>
    static void pointBlock( TestContext const& test, TrialContext const& trial,
                            uint16_type i, uint16_type j,
                            Vector const& coefficient, Vector& values )
    {
        static_assert( TestContext::nComponents1 == TrialContext::nComponents1 );
        values.setZero();
        for ( uint16_type component = 0;
              component < TestContext::nComponents1; ++component )
            values.array() += coefficient.array()*
                              test.idBatch( i, component ).transpose().array()*
                              trial.idBatch( j, component ).transpose().array();
    }

    /** @brief Form the complete vector mass matrix. */
    template<typename TestContext, typename TrialContext, typename Vector,
             typename OutputMatrix, typename ScratchMatrix>
    static void contract( TestContext const& test, TrialContext const& trial,
                          Vector const& weightedCoefficient,
                          OutputMatrix& output,
                          ScratchMatrix& testComponent,
                          ScratchMatrix& weightedTrial )
    {
        static_assert( TestContext::nComponents1 == TrialContext::nComponents1 );
        contractComponentTables(
            test, trial, TestContext::nComponents1, weightedCoefficient,
            output, testComponent, weightedTrial,
            [&]( uint16_type c )
            { return test.valueComponentBatchTable( c ); },
            [&]( uint16_type c )
            { return trial.valueComponentBatchTable( c ); } );
    }
};

/** @brief Strategy tag for isotropic scalar H1 stiffness contractions. */
struct scalar_h1_stiffness_contraction
{
    static constexpr size_type context =
        vm::JACOBIAN | vm::KB | vm::GRAD;

    /** @return one scalar gradient inner product. */
    template<typename TestContext, typename TrialContext>
    static auto point( TestContext const& test, TrialContext const& trial,
                       uint16_type i, uint16_type j, uint32_type q )
    {
        using value_type = typename TestContext::value_type;
        value_type result = 0;
        for ( uint16_type direction = 0; direction < TestContext::NDim;
              ++direction )
            result += test.grad( i, 0, direction, q )*
                      trial.grad( j, 0, direction, q );
        return result;
    }

    /** @brief Evaluate one stiffness entry over a quadrature block. */
    template<typename TestContext, typename TrialContext, typename Vector>
    static void pointBlock( TestContext const& test, TrialContext const& trial,
                            uint16_type i, uint16_type j,
                            Vector const& coefficient, Vector& values )
    {
        static_assert( TestContext::NDim == TrialContext::NDim );
        values.setZero();
        for ( uint16_type direction = 0; direction < TestContext::NDim;
              ++direction )
            values.array() += coefficient.array() *
                              test.gradientBatch( i, 0, direction ).transpose().array() *
                              trial.gradientBatch( j, 0, direction ).transpose().array();
    }

    /** @brief Form the complete isotropic scalar H1 stiffness matrix. */
    template<typename TestContext, typename TrialContext, typename Vector,
             typename OutputMatrix, typename ScratchMatrix>
    static void contract( TestContext const& test, TrialContext const& trial,
                          Vector const& weightedCoefficient,
                          OutputMatrix& output,
                          ScratchMatrix& testComponent,
                          ScratchMatrix& weightedTrial )
    {
        static_assert( TestContext::NDim == TrialContext::NDim );
        contractComponentTables(
            test, trial, TestContext::NDim, weightedCoefficient, output,
            testComponent, weightedTrial,
            [&]( uint16_type direction )
            { return test.gradientComponentBatchTable( 0, direction ); },
            [&]( uint16_type direction )
            { return trial.gradientComponentBatchTable( 0, direction ); } );
    }
};

/**
 * @brief Strategy for anisotropic scalar H1 diffusion.
 *
 * The coefficient is the physical tensor @f$K_{\alpha\beta}(x_q)@f$ and
 * the complete local matrix is
 * @f[
 * A_{ij}=\sum_q w_qJ_q\sum_{\alpha,\beta}
 * \partial_\alpha v_i(x_q)K_{\alpha\beta}(x_q)
 * \partial_\beta u_j(x_q).
 * @f]
 * Coefficient components are evaluated once per element and the transformed
 * trial-gradient blocks are contracted with one Eigen matrix product.
 */
struct anisotropic_scalar_h1_stiffness_contraction
{
    static constexpr size_type context =
        vm::JACOBIAN | vm::KB | vm::GRAD;
    static constexpr bool accepts_tensor_coefficient = true;

    /** @return one anisotropic gradient contraction at one point. */
    template<typename TestContext, typename TrialContext,
             typename CoefficientTensor>
    static auto point( TestContext const& test, TrialContext const& trial,
                       CoefficientTensor const& coefficient,
                       uint16_type i, uint16_type j, uint32_type q )
    {
        static_assert( TestContext::NDim == TrialContext::NDim );
        using value_type = typename TestContext::value_type;
        value_type result = 0;
        for ( uint16_type alpha = 0; alpha < TestContext::NDim; ++alpha )
            for ( uint16_type beta = 0; beta < TestContext::NDim; ++beta )
                result += test.grad( i, 0, alpha, q )*
                          coefficient.evalq( alpha, beta, q )*
                          trial.grad( j, 0, beta, q );
        return result;
    }

    /**
     * @brief Form the complete anisotropic scalar H1 matrix.
     *
     * @param coefficient updated tensor-valued DSEL coefficient
     * @param weights physical quadrature weights @f$w_qJ_q@f$
     * @param coefficientComponents reusable `(alpha,beta)`-by-point storage
     */
    template<typename TestContext, typename TrialContext,
             typename CoefficientTensor, typename WeightVector,
             typename OutputMatrix, typename ScratchMatrix,
             typename CoefficientMatrix>
    static void contractTensor( TestContext const& test,
                                TrialContext const& trial,
                                CoefficientTensor const& coefficient,
                                WeightVector const& weights,
                                OutputMatrix& output,
                                ScratchMatrix& testComponents,
                                ScratchMatrix& weightedTrial,
                                CoefficientMatrix& coefficientComponents )
    {
        static_assert( TestContext::NDim == TrialContext::NDim );
        constexpr uint16_type dimension = TestContext::NDim;
        auto const testDofs = static_cast<Eigen::Index>( test.nDofs() );
        auto const trialDofs = static_cast<Eigen::Index>( trial.nDofs() );
        auto const points = weights.size();
        CHECK_EQ( output.rows(), testDofs );
        CHECK_EQ( output.cols(), trialDofs );

        resizeEigenWorkspace(
            coefficientComponents, dimension*dimension, points );
        for ( uint16_type alpha = 0; alpha < dimension; ++alpha )
            for ( uint16_type beta = 0; beta < dimension; ++beta )
            {
                auto const row = static_cast<Eigen::Index>( alpha*dimension + beta );
                for ( Eigen::Index q = 0; q < points; ++q )
                    coefficientComponents( row, q ) =
                        coefficient.evalq( alpha, beta,
                                           static_cast<uint16_type>( q ) );
            }

        resizeEigenWorkspace(
            testComponents, testDofs, dimension*points );
        resizeEigenWorkspace(
            weightedTrial, trialDofs, dimension*points );
        for ( uint16_type alpha = 0; alpha < dimension; ++alpha )
        {
            auto const firstPoint = static_cast<Eigen::Index>( alpha )*points;
            testComponents.middleCols( firstPoint, points ) =
                test.gradientComponentBatchTable( 0, alpha );
            auto transformed = weightedTrial.middleCols( firstPoint, points );
            transformed.setZero();
            for ( uint16_type beta = 0; beta < dimension; ++beta )
            {
                auto const coefficientRow =
                    static_cast<Eigen::Index>( alpha*dimension + beta );
                transformed.array() +=
                    trial.gradientComponentBatchTable( 0, beta ).array().rowwise()*
                    ( weights.array()*coefficientComponents.row( coefficientRow ).transpose().array() ).transpose();
            }
        }
        output.noalias() = testComponents*weightedTrial.transpose();
    }
};

/** @brief Strategy tag for H(curl) curl-curl contractions. */
struct curl_curl_contraction
{
    static constexpr size_type context =
        vm::JACOBIAN | vm::KB | vm::CURL | vm::FIRST_DERIVATIVE;

    /** @return one physical curl inner product. */
    template<typename TestContext, typename TrialContext>
    static auto point( TestContext const& test, TrialContext const& trial,
                       uint16_type i, uint16_type j, uint32_type q )
    {
        constexpr uint16_type components =
            TestContext::nComponents1 == 3 ? 3 : 1;
        using value_type = typename TestContext::value_type;
        value_type result = 0;
        for ( uint16_type component = 0; component < components; ++component )
            result += test.curl( i, component, 0, q )*
                      trial.curl( j, component, 0, q );
        return result;
    }

    /** @brief Evaluate one curl-curl entry over a quadrature block. */
    template<typename TestContext, typename TrialContext, typename Vector>
    static void pointBlock( TestContext const& test, TrialContext const& trial,
                            uint16_type i, uint16_type j,
                            Vector const& coefficient, Vector& values )
    {
        constexpr uint16_type components =
            TestContext::nComponents1 == 3 ? 3 : 1;
        values.setZero();
        for ( uint16_type component = 0; component < components; ++component )
            values.array() += coefficient.array()*
                              test.curlBatch( i, component ).transpose().array()*
                              trial.curlBatch( j, component ).transpose().array();
    }

    /** @brief Form the complete curl-curl matrix. */
    template<typename TestContext, typename TrialContext, typename Vector,
             typename OutputMatrix, typename ScratchMatrix>
    static void contract( TestContext const& test, TrialContext const& trial,
                          Vector const& weightedCoefficient,
                          OutputMatrix& output,
                          ScratchMatrix& testComponent,
                          ScratchMatrix& weightedTrial )
    {
        constexpr uint16_type components =
            TestContext::nComponents1 == 3 ? 3 : 1;
        contractComponentTables(
            test, trial, components, weightedCoefficient, output,
            testComponent, weightedTrial,
            [&]( uint16_type c )
            { return test.curlComponentBatchTable( c ); },
            [&]( uint16_type c )
            { return trial.curlComponentBatchTable( c ); } );
    }
};

/** @brief Strategy tag for H(div) divergence-divergence contractions. */
struct divergence_divergence_contraction
{
    static constexpr size_type context =
        vm::JACOBIAN | vm::KB | vm::DIV | vm::FIRST_DERIVATIVE;

    /** @return one physical divergence product. */
    template<typename TestContext, typename TrialContext>
    static auto point( TestContext const& test, TrialContext const& trial,
                       uint16_type i, uint16_type j, uint32_type q )
    {
        return test.div( i, 0, 0, q )*trial.div( j, 0, 0, q );
    }

    /** @brief Evaluate one divergence-divergence quadrature block. */
    template<typename TestContext, typename TrialContext, typename Vector>
    static void pointBlock( TestContext const& test, TrialContext const& trial,
                            uint16_type i, uint16_type j,
                            Vector const& coefficient, Vector& values )
    {
        values.array() = coefficient.array()*
                         test.divergenceBatch( i ).transpose().array()*
                         trial.divergenceBatch( j ).transpose().array();
    }

    /** @brief Form the complete divergence-divergence matrix. */
    template<typename TestContext, typename TrialContext, typename Vector,
             typename OutputMatrix, typename ScratchMatrix>
    static void contract( TestContext const& test, TrialContext const& trial,
                          Vector const& weightedCoefficient,
                          OutputMatrix& output,
                          ScratchMatrix& testComponent,
                          ScratchMatrix& weightedTrial )
    {
        contractComponentTables(
            test, trial, 1, weightedCoefficient, output,
            testComponent, weightedTrial,
            [&]( uint16_type )
            { return test.divergenceComponentBatchTable(); },
            [&]( uint16_type )
            { return trial.divergenceComponentBatchTable(); } );
    }
};

/** @brief Strategy tag for scalar-test/H(div)-trial mixed divergence forms. */
struct value_divergence_contraction
{
    static constexpr size_type context =
        vm::JACOBIAN | vm::KB | vm::DIV | vm::FIRST_DERIVATIVE;

    /** @return test value times trial divergence. */
    template<typename TestContext, typename TrialContext>
    static auto point( TestContext const& test, TrialContext const& trial,
                       uint16_type i, uint16_type j, uint32_type q )
    {
        return test.id( i, 0, 0, q )*trial.div( j, 0, 0, q );
    }

    /** @brief Evaluate one mixed divergence quadrature block. */
    template<typename TestContext, typename TrialContext, typename Vector>
    static void pointBlock( TestContext const& test, TrialContext const& trial,
                            uint16_type i, uint16_type j,
                            Vector const& coefficient, Vector& values )
    {
        values.array() = coefficient.array()*
                         test.idBatch( i ).transpose().array()*
                         trial.divergenceBatch( j ).transpose().array();
    }

    /** @brief Form the complete mixed divergence matrix. */
    template<typename TestContext, typename TrialContext, typename Vector,
             typename OutputMatrix, typename ScratchMatrix>
    static void contract( TestContext const& test, TrialContext const& trial,
                          Vector const& weightedCoefficient,
                          OutputMatrix& output,
                          ScratchMatrix& testComponent,
                          ScratchMatrix& weightedTrial )
    {
        contractComponentTables(
            test, trial, 1, weightedCoefficient, output,
            testComponent, weightedTrial,
            [&]( uint16_type )
            { return test.valueComponentBatchTable(); },
            [&]( uint16_type )
            { return trial.divergenceComponentBatchTable(); } );
    }
};

/** @brief Strategy tag for H(div)-test/scalar-trial mixed divergence forms. */
struct divergence_value_contraction
{
    static constexpr size_type context =
        vm::JACOBIAN | vm::KB | vm::DIV | vm::FIRST_DERIVATIVE;

    /** @return test divergence times trial value. */
    template<typename TestContext, typename TrialContext>
    static auto point( TestContext const& test, TrialContext const& trial,
                       uint16_type i, uint16_type j, uint32_type q )
    {
        return test.div( i, 0, 0, q )*trial.id( j, 0, 0, q );
    }

    /** @brief Evaluate one adjoint mixed-divergence quadrature block. */
    template<typename TestContext, typename TrialContext, typename Vector>
    static void pointBlock( TestContext const& test, TrialContext const& trial,
                            uint16_type i, uint16_type j,
                            Vector const& coefficient, Vector& values )
    {
        values.array() = coefficient.array()*
                         test.divergenceBatch( i ).transpose().array()*
                         trial.idBatch( j ).transpose().array();
    }

    /** @brief Form the complete adjoint mixed-divergence matrix. */
    template<typename TestContext, typename TrialContext, typename Vector,
             typename OutputMatrix, typename ScratchMatrix>
    static void contract( TestContext const& test, TrialContext const& trial,
                          Vector const& weightedCoefficient,
                          OutputMatrix& output,
                          ScratchMatrix& testComponent,
                          ScratchMatrix& weightedTrial )
    {
        contractComponentTables(
            test, trial, 1, weightedCoefficient, output,
            testComponent, weightedTrial,
            [&]( uint16_type )
            { return test.divergenceComponentBatchTable(); },
            [&]( uint16_type )
            { return trial.valueComponentBatchTable(); } );
    }
};

} // namespace detail
} // namespace vf
} // namespace Feel

#endif
