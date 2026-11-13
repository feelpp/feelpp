/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELVF_FACET_HPP
#define FEELPP_FEELVF_FACET_HPP 1

#include <concepts>
#include <type_traits>
#include <utility>

#include <feel/feelmesh/filters.hpp>
#include <feel/feelvf/basis.hpp>
#include <feel/feelvf/twovalued.hpp>

namespace Feel::vf::facet
{

namespace detail
{
/**
 * @brief Basis-role aliases carried by a DSEL expression.
 *
 * Non-expression objects, including value-like objects before promotion, are
 * treated as role-free.
 */
template <typename ExprType>
struct FacetExpressionBasis
{
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;
};

/**
 * @brief Basis-role aliases extracted from an expression exposing them.
 */
template <typename ExprType>
    requires requires {
        typename std::remove_cvref_t<ExprType>::test_basis;
        typename std::remove_cvref_t<ExprType>::trial_basis;
    }
struct FacetExpressionBasis<ExprType>
{
    using expression_type = std::remove_cvref_t<ExprType>;
    using test_basis = typename expression_type::test_basis;
    using trial_basis = typename expression_type::trial_basis;
};

template <typename ExprType>
using facet_test_basis_t = typename FacetExpressionBasis<ExprType>::test_basis;

template <typename ExprType>
using facet_trial_basis_t = typename FacetExpressionBasis<ExprType>::trial_basis;

/**
 * @brief True when an expression reports a test function for a basis type.
 */
template <typename ExprType, typename BasisType>
concept FacetExpressionHasTestFunction =
    !std::same_as<BasisType, std::nullptr_t> &&
    requires {
        { std::remove_cvref_t<ExprType>::template HasTestFunction<BasisType>::result }
            -> std::convertible_to<bool>;
    } &&
    std::remove_cvref_t<ExprType>::template HasTestFunction<BasisType>::result;

/**
 * @brief True when an expression reports a trial function for a basis type.
 */
template <typename ExprType, typename BasisType>
concept FacetExpressionHasTrialFunction =
    !std::same_as<BasisType, std::nullptr_t> &&
    requires {
        { std::remove_cvref_t<ExprType>::template HasTrialFunction<BasisType>::result }
            -> std::convertible_to<bool>;
    } &&
    std::remove_cvref_t<ExprType>::template HasTrialFunction<BasisType>::result;

/**
 * @brief Expression carries or reports a test role.
 */
template <typename ExprType>
concept FacetExpressionWithTestRole =
    !std::same_as<facet_test_basis_t<ExprType>, std::nullptr_t> ||
    FacetExpressionHasTestFunction<ExprType, facet_test_basis_t<ExprType>> ||
    FacetExpressionHasTestFunction<ExprType, facet_trial_basis_t<ExprType>>;

/**
 * @brief Expression carries or reports a trial role.
 */
template <typename ExprType>
concept FacetExpressionWithTrialRole =
    !std::same_as<facet_trial_basis_t<ExprType>, std::nullptr_t> ||
    FacetExpressionHasTrialFunction<ExprType, facet_test_basis_t<ExprType>> ||
    FacetExpressionHasTrialFunction<ExprType, facet_trial_basis_t<ExprType>>;

/**
 * @brief Expression carries only a test role.
 */
template <typename ExprType>
concept FacetExpressionWithOnlyTestRole =
    FacetExpressionWithTestRole<ExprType> &&
    !FacetExpressionWithTrialRole<ExprType>;

/**
 * @brief Expression carries only a trial role.
 */
template <typename ExprType>
concept FacetExpressionWithOnlyTrialRole =
    FacetExpressionWithTrialRole<ExprType> &&
    !FacetExpressionWithTestRole<ExprType>;

/**
 * @brief Expression carries both test and trial roles and needs explicit side choice.
 */
template <typename ExprType>
concept FacetExpressionWithAmbiguousRole =
    FacetExpressionWithTestRole<ExprType> &&
    FacetExpressionWithTrialRole<ExprType>;

/**
 * @brief Expression can be routed by inference or by value/test-side fallback.
 */
template <typename ExprType>
concept FacetExpressionWithInferredRole =
    !FacetExpressionWithAmbiguousRole<ExprType>;

} // namespace detail

/**
 * @brief Interior facets of a mesh.
 *
 * This is a naming helper for skeleton formulations; it forwards to the
 * existing `internalfaces(mesh)` range and does not change range semantics.
 */
template <typename MeshType>
decltype(auto)
interior( MeshType&& mesh )
{
    return internalfaces( std::forward<MeshType>( mesh ) );
}

/**
 * @brief Boundary facets of a mesh.
 *
 * This forwards to the existing `boundaryfaces(mesh)` range.
 */
template <typename MeshType>
decltype(auto)
boundary( MeshType&& mesh )
{
    return boundaryfaces( std::forward<MeshType>( mesh ) );
}

/**
 * @brief Full mesh skeleton, including interior and boundary facets.
 *
 * This forwards to the existing `faces(mesh)` range.
 */
template <typename MeshType>
decltype(auto)
skeleton( MeshType&& mesh )
{
    return faces( std::forward<MeshType>( mesh ) );
}

/**
 * @brief Left trace of a test-side or value expression on an interior facet.
 */
template <typename ExprType>
auto
leftTest( ExprType const& expr )
{
    return leftface( expr );
}

/**
 * @brief Right trace of a test-side or value expression on an interior facet.
 */
template <typename ExprType>
auto
rightTest( ExprType const& expr )
{
    return rightface( expr );
}

/**
 * @brief Left trace of a trial-side expression on an interior facet.
 */
template <typename ExprType>
auto
leftTrial( ExprType const& expr )
{
    return leftfacet( expr );
}

/**
 * @brief Right trace of a trial-side expression on an interior facet.
 */
template <typename ExprType>
auto
rightTrial( ExprType const& expr )
{
    return rightfacet( expr );
}

/**
 * @brief Left trace with side inferred from the DSEL expression role.
 *
 * Test expressions use `leftface`, trial expressions use `leftfacet`, and
 * value-like expressions use the value/test-side convention.
 */
template <typename ExprType>
    requires detail::FacetExpressionWithInferredRole<ExprType>
auto
left( ExprType const& expr )
{
    if constexpr ( detail::FacetExpressionWithOnlyTrialRole<ExprType> )
        return leftTrial( expr );
    else
        return leftTest( expr );
}

/**
 * @brief Reject role-ambiguous left traces.
 */
template <typename ExprType>
    requires detail::FacetExpressionWithAmbiguousRole<ExprType>
auto
left( ExprType const& expr ) = delete;

/**
 * @brief Right trace with side inferred from the DSEL expression role.
 *
 * Test expressions use `rightface`, trial expressions use `rightfacet`, and
 * value-like expressions use the value/test-side convention.
 */
template <typename ExprType>
    requires detail::FacetExpressionWithInferredRole<ExprType>
auto
right( ExprType const& expr )
{
    if constexpr ( detail::FacetExpressionWithOnlyTrialRole<ExprType> )
        return rightTrial( expr );
    else
        return rightTest( expr );
}

/**
 * @brief Reject role-ambiguous right traces.
 */
template <typename ExprType>
    requires detail::FacetExpressionWithAmbiguousRole<ExprType>
auto
right( ExprType const& expr ) = delete;

/**
 * @brief Sum of left and right test-side traces on an interior facet.
 */
template <typename ExprType>
auto
testSum( ExprType const& expr )
{
    return leftTest( expr ) + rightTest( expr );
}

/**
 * @brief Sum of left and right trial-side traces on an interior facet.
 */
template <typename ExprType>
auto
trialSum( ExprType const& expr )
{
    return leftTrial( expr ) + rightTrial( expr );
}

/**
 * @brief Sum of left and right traces with side inferred from the expression.
 */
template <typename ExprType>
    requires detail::FacetExpressionWithInferredRole<ExprType>
auto
sum( ExprType const& expr )
{
    return left( expr ) + right( expr );
}

/**
 * @brief Reject role-ambiguous trace sums.
 */
template <typename ExprType>
    requires detail::FacetExpressionWithAmbiguousRole<ExprType>
auto
sum( ExprType const& expr ) = delete;

/**
 * @brief Mean of left and right test-side traces on an interior facet.
 */
template <typename ExprType>
auto
testMean( ExprType const& expr )
{
    return 0.5 * testSum( expr );
}

/**
 * @brief Mean of left and right trial-side traces on an interior facet.
 */
template <typename ExprType>
auto
trialMean( ExprType const& expr )
{
    return 0.5 * trialSum( expr );
}

/**
 * @brief Mean of left and right traces with side inferred from the expression.
 */
template <typename ExprType>
    requires detail::FacetExpressionWithInferredRole<ExprType>
auto
mean( ExprType const& expr )
{
    return 0.5 * sum( expr );
}

/**
 * @brief Reject role-ambiguous trace means.
 */
template <typename ExprType>
    requires detail::FacetExpressionWithAmbiguousRole<ExprType>
auto
mean( ExprType const& expr ) = delete;

/**
 * @brief Side-matched product of a test trace and a trial trace.
 *
 * This is the explicit two-sided pairing
 * `leftface(v)*leftfacet(u) + rightface(v)*rightfacet(u)`.
 */
template <typename TestExprType, typename TrialExprType>
auto
testTrialSum( TestExprType const& testExpr, TrialExprType const& trialExpr )
{
    return leftTest( testExpr ) * leftTrial( trialExpr ) +
           rightTest( testExpr ) * rightTrial( trialExpr );
}

/**
 * @brief Transparent trace notation for facet-space expressions.
 *
 * This helper keeps the mathematical skeleton syntax visible without adding a
 * projection or changing the expression.
 */
template <typename ExprType>
auto
trace( ExprType const& expr )
{
    return expr;
}

/**
 * @brief Trace test expression living on a facet space.
 *
 * The helper is intentionally transparent: it returns the expression unchanged
 * and lets the formulation still show where the trace test appears.
 */
template <typename ExprType>
auto
traceTest( ExprType const& expr )
{
    return trace( expr );
}

/**
 * @brief Trace trial expression living on a facet space.
 *
 * The helper is intentionally transparent: it returns the expression unchanged
 * and does not introduce a projection or boundary condition.
 */
template <typename ExprType>
auto
traceTrial( ExprType const& expr )
{
    return trace( expr );
}

/**
 * @brief Boundary normal trace with role inferred from the expression.
 */
template <typename ExprType>
    requires ( detail::FacetExpressionWithInferredRole<ExprType> &&
               requires ( ExprType const& expr ) { Feel::vf::normal( expr ); } )
auto
normal( ExprType const& expr )
{
    return Feel::vf::normal( expr );
}

/**
 * @brief Reject role-ambiguous normal traces.
 */
template <typename ExprType>
    requires detail::FacetExpressionWithAmbiguousRole<ExprType>
auto
normal( ExprType const& expr ) = delete;

/**
 * @brief Boundary normal trace of a test basis expression.
 */
template <typename TestBasisType>
auto
normalTest( TestBasisType const& basis )
{
    return Feel::vf::facet::normal( basis );
}

/**
 * @brief Boundary normal trace of a trial basis expression.
 */
template <typename TrialBasisType>
auto
normalTrial( TrialBasisType const& basis )
{
    return Feel::vf::facet::normal( basis );
}

/**
 * @brief Sum of left and right normal traces with side inferred by role.
 */
template <typename ExprType>
    requires ( detail::FacetExpressionWithInferredRole<ExprType> &&
               requires ( ExprType const& expr ) { Feel::vf::facet::normal( expr ); } )
auto
normalSum( ExprType const& expr )
{
    return sum( Feel::vf::facet::normal( expr ) );
}

/**
 * @brief Reject role-ambiguous normal trace sums.
 */
template <typename ExprType>
    requires detail::FacetExpressionWithAmbiguousRole<ExprType>
auto
normalSum( ExprType const& expr ) = delete;

/**
 * @brief Sum of left and right normal test traces on an interior facet.
 */
template <typename TestBasisType>
auto
normalTestSum( TestBasisType const& basis )
{
    return testSum( normalTest( basis ) );
}

/**
 * @brief Sum of left and right normal trial traces on an interior facet.
 */
template <typename TrialBasisType>
auto
normalTrialSum( TrialBasisType const& basis )
{
    return trialSum( normalTrial( basis ) );
}

} // namespace Feel::vf::facet

#endif // FEELPP_FEELVF_FACET_HPP
