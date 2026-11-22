/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/** @file precomputetabulation.hpp
 * @brief Persistent reference-table ownership for existing FE precompute adapters.
 */
#ifndef FEELPP_FEELPOLY_PRECOMPUTETABULATION_HPP
#define FEELPP_FEELPOLY_PRECOMPUTETABULATION_HPP 1
#include <feel/feelpoly/tabulation.hpp>
#include <optional>
#include <variant>

namespace Feel::detail
{

/** @brief Canonical column-major owners consumed by geometry and compatibility layouts. */
template <typename Source>
using PrecomputeReferenceTables = ReferenceTabulation<typename Source::value_type,
                                                      Source::nDim, orderSpecification<Source>(), Dynamic, Dynamic, Evaluation::All,
                                                      Source::nComponents, Eigen::ColMajor>;

/** @brief Reuse one selected derivative-order workspace and canonical reference tables.
 * @tparam Source Existing FE polynomial representation; no alternative evaluator is introduced.
 * @tparam Native Whether this source supplies a native modal tabulation kernel.
 */
template <typename Source, bool Native = NativeTabulationSource<Source>>
class PrecomputeTabulationState
{

    /** Compile-time requests spanning each supported maximum derivative order. */
    using Workspaces = std::variant<TabulationWorkspace<Source, Evaluation::Values>,
                                    TabulationWorkspace<Source, Evaluation::Values | Evaluation::Gradients>,
                                    TabulationWorkspace<Source, Evaluation::All>>;

  public:
    /** Native table owner; only requested groups acquire coefficient storage. */
    using Tables = PrecomputeReferenceTables<Source>;

    /** @brief Initialize degree metadata without allocating point tables. */
    explicit PrecomputeTabulationState( RuntimeOrder degree = RuntimeOrder{
                                            orderSpecification<Source>() == Dynamic ? 0 : orderSpecification<Source>() } )
        : M_tables( degree ) {}

    /** @brief Fill selected tables once, retaining all buffers for subsequent updates. */
    template <typename Points>
    void compute( Source const& source, Points const& points, Evaluation request )
    {
        // Validate the selection before changing persistent table storage.
        if ( !validEvaluation( request ) )
        {
            throw std::invalid_argument( "invalid precompute request" );
        }

        auto const degree = polynomialOrder( source );

        // Rebind degree-dependent scratch only when the effective degree changes.
        if ( M_tables.order() != degree )
        {
            M_tables = Tables{ RuntimeOrder{ degree } };
            M_workspace.reset();
        }

        // Dispatch once and reuse the workspace for the highest requested derivative.
        auto evaluateSelected = [&]( auto selection )
        {
            constexpr auto selected = decltype( selection )::value;
            constexpr auto derivative = derivativeOrder( selected );

            if ( !M_workspace || M_workspace->index() != derivative )
            {
                M_workspace.emplace( std::in_place_index<derivative>, RuntimeOrder{ degree } );
            }

            auto output = M_tables.template select<selected>();
            tabulateInto( source, points, output, std::get<derivative>( *M_workspace ) );
        };

        withEvaluationRequest( request, evaluateSelected );

        // Release groups excluded by the new request after evaluation succeeds.
        auto clearUnrequested = [&]<unsigned Derivative>( std::size_t slot, int, int )
        {
            constexpr auto group = static_cast<Evaluation>( 1u << Derivative );

            if ( !includes( request, group ) )
            {
                M_tables.table( slot ).resizeForOverwrite( 0, 0 );
            }
        };

        forEachDerivative<Source::nDim, Evaluation::All>( clearUnrequested );
        M_request = request;
    }

    /** @brief Borrow canonical owners while the state remains unchanged. */
    [[nodiscard]] Tables const& tables() const& noexcept { return M_tables; }

    /** @brief Current request controlling native and compatibility access. */
    [[nodiscard]] Evaluation request() const noexcept { return M_request; }

    /** @brief Reject access to an unrequested derivative group. */
    void require( Evaluation group ) const
    {
        if ( !includes( M_request, group ) ) throw std::logic_error( "precompute table was not requested" );
    }

  private:
    /** Sole canonical reference-table owners. */
    Tables M_tables;

    /** Only the required maximum derivative order is constructed. */
    std::optional<Workspaces> M_workspace;

    /** Last successfully computed table selection. */
    Evaluation M_request = Evaluation::All;
};

/** @brief Adapt an unmigrated basis or zero-dimensional reference through its existing evaluators. */
template <typename Source>
class PrecomputeTabulationState<Source, false>
{
  public:
    /** Canonical owners shared with native precompute consumers. */
    using Tables = PrecomputeReferenceTables<Source>;

    /** @brief Retain the source degree while leaving table storage empty. */
    explicit PrecomputeTabulationState( RuntimeOrder degree = RuntimeOrder{
                                            orderSpecification<Source>() == Dynamic ? 0 : orderSpecification<Source>() } )
        : M_tables( degree ) {}

    /** @brief Fill requested owners from existing family evaluators without new numerical kernels. */
    template <typename Points>
    void compute( Source const& source, Points const& points, Evaluation request )
    {
        // Preserve the same request and degree checks for unmigrated basis families.
        if ( !validEvaluation( request ) ) throw std::invalid_argument( "invalid precompute request" );

        auto const degree = polynomialOrder( source );
        if ( M_tables.order() != degree ) M_tables = Tables{ RuntimeOrder{ degree } };

        // Adapt existing value evaluation directly into the canonical owner.
        if ( includes( request, Evaluation::Values ) )
            M_tables.values().assign( source.values( points ) );
        else
            M_tables.values().resizeForOverwrite( 0, 0 );

        // Zero-dimensional references have no coordinate derivatives to evaluate.
        if constexpr ( Source::nDim > 0 )
        {
            // Keep each legacy numerical evaluator behind the shared table layout.
            if ( includes( request, Evaluation::Gradients ) )
            {
                auto gradients = source.gradients( points );
                for ( int direction = 0; direction < Source::nDim; ++direction )
                    M_tables.gradient( direction ).assign( gradients[direction] );
            }

            if ( includes( request, Evaluation::Hessians ) )
            {
                auto hessians = source.hessians( points );
                for ( int i = 0; i < Source::nDim; ++i )
                    for ( int j = 0; j < Source::nDim; ++j )
                        M_tables.hessian( i, j ).assign( hessians[i * Source::nDim + j] );
            }
        }

        // Publish the new request only after all selected evaluations succeed.
        M_request = request;
    }

    /** @brief Borrow the canonical value owner. */
    [[nodiscard]] Tables const& tables() const& noexcept { return M_tables; }

    /** @brief Current table request. */
    [[nodiscard]] Evaluation request() const noexcept { return M_request; }

    /** @brief Reject access to an unrequested group. */
    void require( Evaluation group ) const
    {
        if ( !includes( M_request, group ) ) throw std::logic_error( "precompute table was not requested" );
    }

  private:
    /** Canonical owners filled through legacy family evaluators. */
    Tables M_tables;

    /** Last computed request. */
    Evaluation M_request = Evaluation::All;
};
} // namespace Feel::detail
#endif // FEELPP_FEELPOLY_PRECOMPUTETABULATION_HPP
