/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_ASSEMBLY_HPP
#define FEELPP_FEELSYSTEM_ASSEMBLY_HPP 1

/**
 * @file assembly.hpp
 * @brief Customization-point facade that assembles problem descriptors.
 */

#include <string>
#include <utility>

#include <feel/feelcore/parameter.hpp>
#include <feel/feelsystem/eigen.hpp>
#include <feel/feelsystem/linear.hpp>
#include <feel/feelsystem/block.hpp>
#include <feel/feelsystem/hdg.hpp>
#include <feel/feelsystem/nonlinear.hpp>
#include <feel/feelsystem/picard.hpp>
#include <feel/feelsystem/problem.hpp>

namespace Feel::systems
{

/**
 * @brief Function object assembling FeelSystem problem descriptors.
 *
 * `assemble_t` selects the concrete system wrapper from the problem concept:
 * `LinearSystem`, `PicardSystem`, `NonlinearSystem`, `EigenSystem`,
 * `BlockSystem`, or `HdgSystem`. Optional named arguments include `_backend`,
 * `_name`, and for block-like systems `_strategy`.
 */
struct assemble_t
{
    /**
     * @brief Assemble a scalar stationary linear problem.
     */
    template <LinearProblemLike ProblemType, typename... Args>
        requires ( !GeneralizedEigenvalueProblemLike<ProblemType> )
    auto operator()( ProblemType const& problem, Args&&... args ) const
    {
        using system_type = LinearSystem<ProblemType>;
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto backend = namedArgs.get_else( _backend, typename system_type::backend_ptrtype{} );
        std::string backendName = namedArgs.get_else( _name, std::string{} );
        auto system = system_type( problem, backend, std::move( backendName ) );
        system.assemble();
        return system;
    }

    /**
     * @brief Assemble a Picard fixed-point problem.
     */
    template <PicardProblemLike ProblemType, typename... Args>
    auto operator()( ProblemType const& problem, Args&&... args ) const
    {
        using system_type = PicardSystem<ProblemType>;
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto backend = namedArgs.get_else( _backend, typename system_type::backend_ptrtype{} );
        std::string backendName = namedArgs.get_else( _name, std::string{} );
        auto system = system_type( problem, backend, std::move( backendName ) );
        system.assemble();
        return system;
    }

    /**
     * @brief Assemble a residual/Jacobian nonlinear problem.
     */
    template <NonlinearProblemLike ProblemType, typename... Args>
    auto operator()( ProblemType const& problem, Args&&... args ) const
    {
        using system_type = NonlinearSystem<ProblemType>;
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto backend = namedArgs.get_else( _backend, typename system_type::backend_ptrtype{} );
        std::string backendName = namedArgs.get_else( _name, std::string{} );
        auto system = system_type( problem, backend, std::move( backendName ) );
        system.assemble();
        return system;
    }

    /**
     * @brief Assemble a standard or generalized eigenvalue problem.
     */
    template <GeneralizedEigenvalueProblemLike ProblemType, typename... Args>
    auto operator()( ProblemType const& problem, Args&&... args ) const
    {
        using system_type = EigenSystem<ProblemType>;
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto backend = namedArgs.get_else( _backend, typename system_type::backend_ptrtype{} );
        std::string backendName = namedArgs.get_else( _name, std::string{} );
        auto system = system_type( problem, backend, std::move( backendName ) );
        system.assemble();
        return system;
    }

    /**
     * @brief Assemble a product-space block problem.
     */
    template <BlockProblemLike ProblemType, typename... Args>
        requires ( !HdgProblemLike<ProblemType> )
    auto operator()( ProblemType const& problem, Args&&... args ) const
    {
        using system_type = BlockSystem<ProblemType>;
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto backend = namedArgs.get_else( _backend, typename system_type::backend_ptrtype{} );
        std::string backendName = namedArgs.get_else( _name, std::string{} );
        auto strategy = namedArgs.get_else( _strategy, problem.strategy() );
        auto system = system_type( problem, backend, std::move( backendName ), strategy );
        system.assemble();
        return system;
    }

    /**
     * @brief Assemble an HDG problem through its policy-aware member facade.
     */
    template <HdgProblemLike ProblemType, typename... Args>
    auto operator()( ProblemType const& problem, Args&&... args ) const
    {
        return problem.assemble( std::forward<Args>( args )... );
    }
};

/**
 * @brief Global FeelSystem assembly customization-point object.
 */
inline constexpr assemble_t assemble{};

/**
 * @brief Member implementation of `Problem::assemble()`.
 */
template <typename SpaceType, typename... Tags>
template <typename... Args>
auto
Problem<SpaceType, Tags...>::assemble( Args&&... args ) const
{
    return systems::assemble( *this, std::forward<Args>( args )... );
}

/**
 * @brief Member implementation of `BlockProblem::assemble()`.
 */
template <typename ProductSpaceType, typename... Tags>
    requires StaticProductSpacesType<ProductSpaceType>
template <typename... Args>
auto
BlockProblem<ProductSpaceType, Tags...>::assemble( Args&&... args ) const
{
    return systems::assemble( *this, std::forward<Args>( args )... );
}

} // namespace Feel::systems

namespace Feel
{
inline constexpr auto const& assemble = systems::assemble;
}

#endif // FEELPP_FEELSYSTEM_ASSEMBLY_HPP
