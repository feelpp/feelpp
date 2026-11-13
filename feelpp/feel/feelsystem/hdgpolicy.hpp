/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_HDGPOLICY_HPP
#define FEELPP_FEELSYSTEM_HDGPOLICY_HPP 1

#include <string>
#include <type_traits>
#include <utility>

#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/productspaceconcepts.hpp>
#include <feel/options.hpp>

namespace Feel::systems
{

/**
 * @brief Static-condensation policy carried by an HDG problem descriptor.
 *
 * The condenser is a template parameter because the existing block-form solve
 * path dispatches on condenser tag type. The policy owns no matrix or backend;
 * it only records the default strategy and named-argument values that an
 * `HdgSystem` should pass to `blockform2::solve()`.
 *
 * @tparam CondenserType Feel++ condenser tag such as `condenser_poisson`.
 */
template <CondenserTag CondenserType = condenser_poisson>
class HdgCondensationPolicy
{
public:
    using condenser_type = std::remove_cvref_t<CondenserType>;

    /**
     * @brief Construct a disabled monolithic policy.
     */
    HdgCondensationPolicy() = default;

    /**
     * @brief Construct a policy with an explicit condenser tag.
     */
    explicit HdgCondensationPolicy( condenser_type condenser )
        : M_condenser( std::move( condenser ) )
    {}

    /**
     * @brief Return true when solve should pass `_condense = true`.
     */
    bool enabled() const noexcept { return M_enabled; }

    /**
     * @brief Enable or disable solve-time static condensation.
     */
    HdgCondensationPolicy& enabled( bool value ) noexcept
    {
        M_enabled = value;
        M_strategy = value ? solve::strategy::static_condensation : solve::strategy::monolithic;
        return *this;
    }

    /**
     * @brief Alias for `.enabled(value)` matching block-form terminology.
     */
    HdgCondensationPolicy& condense( bool value ) noexcept
    {
        return this->enabled( value );
    }

    /**
     * @brief Block assembly strategy implied by the policy.
     */
    solve::strategy strategy() const noexcept { return M_strategy; }

    /**
     * @brief Override the block assembly strategy explicitly.
     */
    HdgCondensationPolicy& strategy( solve::strategy value ) noexcept
    {
        M_strategy = value;
        M_enabled = ( value == solve::strategy::static_condensation );
        return *this;
    }

    /**
     * @brief Condenser tag passed to `blockform2::solve()`.
     */
    condenser_type const& condenser() const noexcept { return M_condenser; }

    /**
     * @brief Return a copy of this policy with another condenser tag.
     */
    template <CondenserTag NewCondenserType>
    [[nodiscard]] auto condenser( NewCondenserType condenser ) const
    {
        HdgCondensationPolicy<std::remove_cvref_t<NewCondenserType>> policy( std::move( condenser ) );
        policy.enabled( M_enabled );
        policy.strategy( M_strategy );
        if ( M_hasBackendName )
            policy.backendName( M_backendName );
        return policy;
    }

    /**
     * @brief Return true when a condensed backend name was specified.
     */
    bool hasBackendName() const noexcept { return M_hasBackendName; }

    /**
     * @brief Set the final condensed backend name, e.g. `prefixvm(name,"sc")`.
     */
    HdgCondensationPolicy& backendName( std::string name )
    {
        M_backendName = std::move( name );
        M_hasBackendName = true;
        return *this;
    }

    /**
     * @brief Final backend name used by the condensed trace solve.
     */
    std::string backendName( std::string const& problemName ) const
    {
        if ( M_hasBackendName )
            return M_backendName;
        return prefixvm( problemName, "sc" );
    }

    /**
     * @brief Base name to pass to blockform solves for the condensed backend.
     *
     * Current block forms create the trace/global backend as
     * `prefixvm(_name,"sc")`. If the policy stores a final backend name with
     * that suffix, this helper returns the corresponding base name.
     */
    std::string solveBackendBaseName( std::string const& problemName ) const
    {
        if ( !M_hasBackendName )
            return problemName;

        std::string const marker = "__feelpp_hdg_backend";
        std::string const markerWithSuffix = prefixvm( marker, "sc" );
        if ( !markerWithSuffix.starts_with( marker ) )
            return problemName;

        std::string const suffix = markerWithSuffix.substr( marker.size() );
        if ( suffix.empty() || !M_backendName.ends_with( suffix ) )
            return problemName;

        return M_backendName.substr( 0, M_backendName.size() - suffix.size() );
    }

    /**
     * @brief Default backend name for local postprocessing solves.
     */
    std::string postBackendName( std::string const& problemName ) const
    {
        return prefixvm( problemName, "post" );
    }

private:
    /** @brief True when solve-time static condensation should be requested. */
    bool M_enabled = false;
    /** @brief Block-form strategy implied or selected by the policy. */
    solve::strategy M_strategy = solve::strategy::monolithic;
    /** @brief Condenser tag forwarded to `blockform2::solve()`. */
    condenser_type M_condenser{};
    /** @brief True when `M_backendName` was explicitly set. */
    bool M_hasBackendName = false;
    /** @brief Final condensed backend name. */
    std::string M_backendName;
};

} // namespace Feel::systems

namespace Feel::hdg
{

/**
 * @brief Create the default HDG/static-condensation policy.
 */
[[nodiscard]] inline auto
condensation()
{
    return systems::HdgCondensationPolicy<condenser_poisson>{};
}

} // namespace Feel::hdg

#endif // FEELPP_FEELSYSTEM_HDGPOLICY_HPP
