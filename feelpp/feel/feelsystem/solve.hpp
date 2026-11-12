/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_SOLVE_HPP
#define FEELPP_FEELSYSTEM_SOLVE_HPP 1

/**
 * @file solve.hpp
 * @brief Customization-point facade that assembles when needed and solves.
 */

#include <utility>

#include <feel/feelsystem/assembly.hpp>
#include <feel/feelsystem/concepts.hpp>
#include <feel/feelsystem/namedarguments.hpp>

namespace Feel::systems
{

namespace detail
{
/**
 * @brief Assemble a scalar problem while forwarding only assembly arguments.
 */
template <ProblemLike ProblemType, typename... Args>
auto
assembleForSolve( ProblemType const& problem, Args&&... args )
{
    if constexpr ( hasNamedArgument<na::backend, Args...> && hasNamedArgument<na::name, Args...> )
    {
        return systems::assemble( problem,
                                  _backend = namedArgumentValue<na::backend>( std::forward<Args>( args )... ),
                                  _name = namedArgumentValue<na::name>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::backend, Args...> )
    {
        return systems::assemble( problem,
                                  _backend = namedArgumentValue<na::backend>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::name, Args...> )
    {
        return systems::assemble( problem,
                                  _name = namedArgumentValue<na::name>( std::forward<Args>( args )... ) );
    }
    else
    {
        return systems::assemble( problem );
    }
}

/**
 * @brief Assemble a block problem while forwarding backend, name, and strategy.
 */
template <BlockProblemLike ProblemType, typename... Args>
    requires ( !HdgProblemLike<ProblemType> )
auto
assembleForSolve( ProblemType const& problem, Args&&... args )
{
    if constexpr ( hasNamedArgument<na::backend, Args...> &&
                   hasNamedArgument<na::name, Args...> &&
                   hasNamedArgument<na::strategy, Args...> )
    {
        return systems::assemble( problem,
                                  _backend = namedArgumentValue<na::backend>( std::forward<Args>( args )... ),
                                  _name = namedArgumentValue<na::name>( std::forward<Args>( args )... ),
                                  _strategy = namedArgumentValue<na::strategy>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::backend, Args...> &&
                        hasNamedArgument<na::name, Args...> )
    {
        return systems::assemble( problem,
                                  _backend = namedArgumentValue<na::backend>( std::forward<Args>( args )... ),
                                  _name = namedArgumentValue<na::name>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::backend, Args...> &&
                        hasNamedArgument<na::strategy, Args...> )
    {
        return systems::assemble( problem,
                                  _backend = namedArgumentValue<na::backend>( std::forward<Args>( args )... ),
                                  _strategy = namedArgumentValue<na::strategy>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::name, Args...> &&
                        hasNamedArgument<na::strategy, Args...> )
    {
        return systems::assemble( problem,
                                  _name = namedArgumentValue<na::name>( std::forward<Args>( args )... ),
                                  _strategy = namedArgumentValue<na::strategy>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::backend, Args...> )
    {
        return systems::assemble( problem,
                                  _backend = namedArgumentValue<na::backend>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::name, Args...> )
    {
        return systems::assemble( problem,
                                  _name = namedArgumentValue<na::name>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::strategy, Args...> )
    {
        return systems::assemble( problem,
                                  _strategy = namedArgumentValue<na::strategy>( std::forward<Args>( args )... ) );
    }
    else
    {
        return systems::assemble( problem );
    }
}

/**
 * @brief Assemble an HDG problem while forwarding backend, name, and strategy.
 */
template <HdgProblemLike ProblemType, typename... Args>
auto
assembleForSolve( ProblemType const& problem, Args&&... args )
{
    if constexpr ( hasNamedArgument<na::backend, Args...> &&
                   hasNamedArgument<na::name, Args...> &&
                   hasNamedArgument<na::strategy, Args...> )
    {
        return systems::assemble( problem,
                                  _backend = namedArgumentValue<na::backend>( std::forward<Args>( args )... ),
                                  _name = namedArgumentValue<na::name>( std::forward<Args>( args )... ),
                                  _strategy = namedArgumentValue<na::strategy>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::backend, Args...> &&
                        hasNamedArgument<na::name, Args...> )
    {
        return systems::assemble( problem,
                                  _backend = namedArgumentValue<na::backend>( std::forward<Args>( args )... ),
                                  _name = namedArgumentValue<na::name>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::backend, Args...> &&
                        hasNamedArgument<na::strategy, Args...> )
    {
        return systems::assemble( problem,
                                  _backend = namedArgumentValue<na::backend>( std::forward<Args>( args )... ),
                                  _strategy = namedArgumentValue<na::strategy>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::name, Args...> &&
                        hasNamedArgument<na::strategy, Args...> )
    {
        return systems::assemble( problem,
                                  _name = namedArgumentValue<na::name>( std::forward<Args>( args )... ),
                                  _strategy = namedArgumentValue<na::strategy>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::backend, Args...> )
    {
        return systems::assemble( problem,
                                  _backend = namedArgumentValue<na::backend>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::name, Args...> )
    {
        return systems::assemble( problem,
                                  _name = namedArgumentValue<na::name>( std::forward<Args>( args )... ) );
    }
    else if constexpr ( hasNamedArgument<na::strategy, Args...> )
    {
        return systems::assemble( problem,
                                  _strategy = namedArgumentValue<na::strategy>( std::forward<Args>( args )... ) );
    }
    else
    {
        return systems::assemble( problem );
    }
}

} // namespace detail

/**
 * @brief Function object solving systems or problem descriptors.
 *
 * For already assembled systems, `solve_t` forwards directly to
 * `system.solve(...)`. For problem descriptors, it first assembles an
 * appropriate system while preserving assembly-specific named arguments, then
 * forwards the original argument pack to the system solve.
 */
struct solve_t
{
    /**
     * @brief Solve an already assembled system object.
     */
    template <SystemLike SystemType, typename... Args>
    auto operator()( SystemType&& system, Args&&... args ) const
    {
        return std::forward<SystemType>( system ).solve( std::forward<Args>( args )... );
    }

    /**
     * @brief Assemble then solve a scalar problem descriptor.
     */
    template <ProblemLike ProblemType, typename... Args>
    auto operator()( ProblemType const& problem, Args&&... args ) const
    {
        auto system = detail::assembleForSolve( problem, args... );
        return ( *this )( system, std::forward<Args>( args )... );
    }

    /**
     * @brief Assemble then solve a block problem descriptor.
     */
    template <BlockProblemLike ProblemType, typename... Args>
        requires ( !HdgProblemLike<ProblemType> )
    auto operator()( ProblemType const& problem, Args&&... args ) const
    {
        auto system = detail::assembleForSolve( problem, args... );
        return ( *this )( system, std::forward<Args>( args )... );
    }

    /**
     * @brief Assemble then solve an HDG problem descriptor.
     */
    template <HdgProblemLike ProblemType, typename... Args>
    auto operator()( ProblemType const& problem, Args&&... args ) const
    {
        auto system = detail::assembleForSolve( problem, args... );
        return ( *this )( system, std::forward<Args>( args )... );
    }
};

/**
 * @brief Global FeelSystem solve customization-point object.
 */
inline constexpr solve_t solve{};

/**
 * @brief Member implementation of `Problem::solve()`.
 */
template <typename SpaceType, typename... Tags>
template <typename... Args>
auto
Problem<SpaceType, Tags...>::solve( Args&&... args ) const
{
    return systems::solve( *this, std::forward<Args>( args )... );
}

/**
 * @brief Member implementation of `BlockProblem::solve()`.
 */
template <typename ProductSpaceType, typename... Tags>
    requires StaticProductSpacesType<ProductSpaceType>
template <typename... Args>
auto
BlockProblem<ProductSpaceType, Tags...>::solve( Args&&... args ) const
{
    return systems::solve( *this, std::forward<Args>( args )... );
}

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_SOLVE_HPP
