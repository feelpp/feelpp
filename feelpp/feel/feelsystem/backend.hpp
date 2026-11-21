/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_BACKEND_HPP
#define FEELPP_FEELSYSTEM_BACKEND_HPP 1

/**
 * @file backend.hpp
 * @brief Backend creation helpers shared by FeelSystem system wrappers.
 */

#include <complex>
#include <string>
#include <type_traits>

#include <feel/feelalg/backend.hpp>
#include <feel/options.hpp>

namespace Feel::systems::detail
{

/**
 * @brief Dependent `false` used to reject unsupported scalar types in templates.
 */
template <typename>
inline constexpr bool dependent_false_v = false;

/** @brief Ensure that all options for a backend prefix have defaults. */
inline void
ensureBackendOptions( std::string const& name, po::variables_map& vm )
{
    auto const backendOptionName = prefixvm( name, "backend" );
    if ( vm.count( backendOptionName ) )
        return;

    auto options = Feel::backend_options( name );
    po::variables_map defaultVm;
    char const* defaultArgv[] = { "feelsystem-backend-defaults" };
    po::store( po::parse_command_line( 1, const_cast<char**>( defaultArgv ), options ), defaultVm );
    po::notify( defaultVm );
    for ( auto const& option : defaultVm )
    {
        if ( !vm.count( option.first ) )
            vm.insert( option );
    }
}

/**
 * @brief Select a backend prefix from an explicit request, problem name, or fallback.
 *
 * @param requestedName User-provided backend name.
 * @param problemName Name carried by the problem descriptor.
 * @param fallbackName Hard-coded wrapper fallback.
 * @return First non-empty backend name in that order.
 */
inline std::string
defaultBackendName( std::string const& requestedName,
                    std::string const& problemName,
                    std::string const& fallbackName )
{
    if ( !requestedName.empty() )
        return requestedName;
    if ( !problemName.empty() )
        return problemName;
    return fallbackName;
}

/**
 * @brief Create or reuse a Feel++ backend for a FeelSystem wrapper.
 *
 * The helper registers the backend-manager deleter once per scalar type, adds
 * default backend options when a prefix has not been seen by the global
 * variables map, and reuses the backend manager unless
 * `<name>.backend.rebuild` requests a rebuild.
 *
 * @tparam ValueType Supported scalar type, currently `double` or
 * `std::complex<double>`.
 * @param name Backend option prefix.
 * @param worldComm Communicator associated with the problem space.
 * @return Shared backend instance keyed by backend kind, prefix, and world size.
 */
template <typename ValueType>
std::shared_ptr<Backend<ValueType>>
makeBackend( std::string const& name, worldcomm_ptr_t const& worldComm )
{
    if constexpr ( !( std::is_same_v<ValueType, double> || std::is_same_v<ValueType, std::complex<double>> ) )
    {
        static_assert( dependent_false_v<ValueType>, "Feel::systems default backend creation supports double and std::complex<double>" );
    }
    else
    {
        static bool observed = false;
        if ( !observed )
        {
            Environment::addDeleteObserver( Feel::detail::BackendManagerDeleter<ValueType>::instance() );
            observed = true;
        }

        auto vm = Environment::vm();
        ensureBackendOptions( name, vm );
        auto const backendOptionName = prefixvm( name, "backend" );
        std::string const kind = vm[backendOptionName].template as<std::string>();
        bool const rebuild = vm[prefixvm( name, "backend.rebuild" )].template as<bool>();
        auto key = boost::make_tuple( kind, name, worldComm->globalSize() );
        auto& manager = Feel::detail::BackendManager<ValueType>::instance();
        auto it = manager.find( key );
        if ( it != manager.end() && !rebuild )
            return it->second;

        if ( it != manager.end() && rebuild )
            it->second->clear();

        auto backend = Backend<ValueType>::build( kind, name, worldComm, vm );
        manager[key] = backend;
        return backend;
    }
}

/**
 * @brief Create a sequential backend with the same implementation as another backend.
 */
template <typename ValueType>
std::shared_ptr<Backend<ValueType>>
makeBackendLike( std::shared_ptr<Backend<ValueType>> const& reference,
                 std::string const& name,
                 worldcomm_ptr_t const& worldComm )
{
    CHECK( reference ) << "a reference backend is required";
    auto vm = Environment::vm();
    ensureBackendOptions( name, vm );
    auto const kind = Backend<ValueType>::enumToKind( reference->type() );
    return Backend<ValueType>::build( kind, name, worldComm, vm );
}

} // namespace Feel::systems::detail

#endif // FEELPP_FEELSYSTEM_BACKEND_HPP
