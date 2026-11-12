/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#pragma once

#include <concepts>
#include <cstdint>

#include <Eigen/Core>

namespace Feel
{
namespace quadrature
{

/**
 * @brief Kind of scalar tensor evaluation requested over a quadrature block.
 *
 * The kind distinguishes an expression value, a test-function contribution,
 * and a bilinear test/trial contribution without multiplying the public tensor
 * API by three nearly identical batch concepts.
 */
enum class evaluation_kind
{
    value,
    test,
    bilinear
};

/**
 * @brief Indices associated with one quadrature-block tensor evaluation.
 *
 * A request is constructed once for a complete quadrature block. It carries
 * the same mathematical indices as `evalq`, `evaliq`, or `evalijq`; only the
 * quadrature-point index is replaced by the output vector position.
 */
struct evaluation_request
{
    evaluation_kind kind = evaluation_kind::value;
    std::uint32_t test = 0;
    std::uint32_t trial = 0;
    std::uint16_t component1 = 0;
    std::uint16_t component2 = 0;

    /** @brief Construct a value-only quadrature request. */
    static constexpr evaluation_request value( std::uint16_t c1, std::uint16_t c2 ) noexcept
    {
        return { evaluation_kind::value, 0, 0, c1, c2 };
    }

    /** @brief Construct a test-function quadrature request. */
    static constexpr evaluation_request testFunction( std::uint32_t i,
                                                       std::uint16_t c1,
                                                       std::uint16_t c2 ) noexcept
    {
        return { evaluation_kind::test, i, 0, c1, c2 };
    }

    /** @brief Construct a bilinear test/trial quadrature request. */
    static constexpr evaluation_request bilinearForm( std::uint32_t i,
                                                       std::uint32_t j,
                                                       std::uint16_t c1,
                                                       std::uint16_t c2 ) noexcept
    {
        return { evaluation_kind::bilinear, i, j, c1, c2 };
    }
};

} // namespace quadrature

/**
 * @brief Tensor supporting evaluation of one complete quadrature-point block.
 *
 * This is the single capability concept used by the quadrature layer. A tensor
 * satisfying it writes values into caller-owned Eigen storage; the operation
 * must not resize or retain the output view.
 *
 * @tparam Tensor variational-expression tensor type
 * @tparam Value scalar value type stored by the quadrature rule
 */
template<typename Tensor, typename Value>
concept QuadratureBatchEvaluable =
    requires( Tensor const& tensor,
              quadrature::evaluation_request const& request,
              Eigen::Matrix<Value, Eigen::Dynamic, 1>& values )
    {
        { tensor.evaluateQuadratureBatch( request, values ) } -> std::same_as<void>;
    };

namespace quadrature
{

/**
 * @brief Evaluate a quadrature block through the batch capability or fallback.
 *
 * The capability decision is compile-time. Unsupported tensors retain the
 * established scalar `evalq`/`evaliq`/`evalijq` contract, while supported
 * tensors can use Eigen operations over the caller-owned output vector.
 *
 * @tparam Tensor expression tensor type
 * @tparam Derived Eigen column-vector expression type
 * @param tensor updated expression tensor
 * @param request mathematical test/trial/component indices
 * @param values output vector, one entry per quadrature point
 */
template<typename Tensor, typename Derived>
void evaluateBatchOrFallback(
    Tensor const& tensor,
    evaluation_request const& request,
    Eigen::MatrixBase<Derived>& values )
{
    using value_type = typename Derived::Scalar;
    static_assert( Derived::IsVectorAtCompileTime &&
                   Derived::ColsAtCompileTime == 1,
                   "quadrature batch output must be a column vector" );
    if constexpr ( QuadratureBatchEvaluable<Tensor, value_type> )
    {
        tensor.evaluateQuadratureBatch( request, values.derived() );
    }
    else
    {
        for ( Eigen::Index q = 0; q < values.size(); ++q )
        {
            switch ( request.kind )
            {
            case evaluation_kind::value:
                values[q] = tensor.evalq( request.component1, request.component2,
                                          static_cast<std::uint32_t>( q ) );
                break;
            case evaluation_kind::test:
                values[q] = tensor.evaliq( request.test, request.component1,
                                           request.component2,
                                           static_cast<std::uint32_t>( q ) );
                break;
            case evaluation_kind::bilinear:
                values[q] = tensor.evalijq( request.test, request.trial,
                                            request.component1, request.component2,
                                            static_cast<std::uint32_t>( q ) );
                break;
            }
        }
    }
}

} // namespace quadrature
} // namespace Feel
