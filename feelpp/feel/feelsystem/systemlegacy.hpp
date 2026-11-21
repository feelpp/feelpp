/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_SYSTEMLEGACY_HPP
#define FEELPP_FEELSYSTEM_SYSTEMLEGACY_HPP 1

/**
 * @file systemlegacy.hpp
 * @brief Historical `Feel::System` base class kept for compatibility.
 */

#include <memory>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feeltypes.hpp>

namespace Feel
{
/**
 * Legacy base class kept for source compatibility with the historical
 * `SystemExplicit` and `SystemImplicit*` headers.
 *
 * New problem/system work should use the concept-constrained API in
 * `Feel::systems` through `<feel/feelsystem/system.hpp>`.
 */
template<typename SpaceType>
class System
{
public:
    /** @brief Topological dimension of the function space. */
    static inline const uint16_type Dim = SpaceType::nDim;
    /** @brief Self type. */
    using system_type = System<SpaceType>;
    /** @brief Scalar value type. */
    using value_type = typename SpaceType::value_type;
    /** @brief Function-space type. */
    using functionspace_type = SpaceType;
    /** @brief Shared function-space pointer type. */
    using functionspace_ptrtype = std::shared_ptr<SpaceType>;
    /** @brief Finite-element field type. */
    using element_type = typename functionspace_type::element_type;

    /**
     * @brief Construct from a function space and variables map.
     */
    System( functionspace_ptrtype const& Xh, po::variables_map const& vm )
        : M_vm( vm ), M_Xh( Xh )
    {}
    /** @brief Copy constructor. */
    System( System const& s ) = default;
    /** @brief Virtual destructor for derived legacy systems. */
    virtual ~System() = default;

    /** @brief Copy assignment. */
    System& operator=( System const& s ) = default;

    /** @brief Variables map used by the legacy system. */
    po::variables_map const& vm() const { return M_vm; }
    /** @brief Function space used by the legacy system. */
    functionspace_ptrtype const& functionSpace() const { return M_Xh; }

    /** @brief Replace the variables map. */
    void setVm( po::variables_map const& vm ) { M_vm = vm; }
    /** @brief Replace the function space. */
    void setFunctionSpace( functionspace_ptrtype const& Xh ) { M_Xh = Xh; }

    /** @brief Assemble the legacy system. */
    virtual void assemble() = 0;
    /** @brief Solve into finite-element field `u`. */
    virtual void solve( element_type& u ) = 0;

private:
    /** @brief Variables map copied into the system. */
    po::variables_map M_vm;
    /** @brief Function space used by the system. */
    functionspace_ptrtype M_Xh;
};

} // namespace Feel

#endif // FEELPP_FEELSYSTEM_SYSTEMLEGACY_HPP
