/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#ifndef FEELPP_FEELPOLY_BROKEN_HPP
#define FEELPP_FEELPOLY_BROKEN_HPP 1

#include <feel/feeldiscr/doflayout.hpp>
#include <feel/feelpoly/discontinuous.hpp>

namespace Feel
{
namespace fem
{
/**
 * @brief Decorate a reference finite element with cell-local assembly topology.
 *
 * The wrapped Ciarlet element, tabulation, interpolation, and physical mapping
 * are preserved. Only continuity and DoF ownership are changed: every local
 * DoF is attached to the cell and no entity orientation is shared globally.
 *
 * @tparam FE wrapped concrete reference finite-element type
 */
template<class FE>
class BrokenElement : public FE
{
    using super = FE;

public:
    using wrapped_fe_type = FE; /**< Wrapped conforming or nonconforming FE. */
    using continuity_type = Discontinuous; /**< Totally discontinuous continuity policy. */

    static constexpr bool isContinuous = false; /**< Broken functions are not globally continuous. */
    static constexpr bool is_broken = true; /**< Identifies the assembly-topology modifier. */

    using super::super;
    BrokenElement() = default;

    /** @return wrapped FE instance. */
    [[nodiscard]] super const& wrapped() const noexcept { return *this; }
    /** @return wrapped FE instance. */
    [[nodiscard]] super& wrapped() noexcept { return *this; }

    /** @return cell attachment for every local DoF. */
    [[nodiscard]] typename super::DofAttachment dofAttachment( uint16_type localDofId ) const override
    {
        return typename super::DofAttachment{
            .entityDim = static_cast<int8_type>( super::nDim ),
            .entityId = 0,
            .ordinal = this->dofParent( localDofId ),
            .kind = this->dofType( localDofId ) };
    }

    /** @return cell-only entity cardinality. */
    [[nodiscard]] uint16_type localDofCountOnEntity( uint16_type topologicalDim,
                                                      uint16_type localEntity,
                                                      bool perComponent = false ) const override
    {
        return topologicalDim == super::nDim && localEntity == 0
                   ? this->localDofCount( perComponent )
                   : 0;
    }

    /** @return zero because broken DoFs do not belong to a shared facet. */
    [[nodiscard]] uint16_type localDofCountOnFacet( uint16_type = 0,
                                                    bool = false ) const override
    {
        return 0;
    }

    /** @return zero broken DoFs per vertex. */
    [[nodiscard]] constexpr uint16_type dofPerVertex() const noexcept { return 0; }
    /** @return zero broken DoFs per edge. */
    [[nodiscard]] constexpr uint16_type dofPerEdge() const noexcept { return 0; }
    /** @return zero broken DoFs per face. */
    [[nodiscard]] constexpr uint16_type dofPerFace() const noexcept { return 0; }
    /** @return all per-component DoFs on the cell. */
    [[nodiscard]] uint16_type dofPerVolume() const noexcept { return this->localDofCount( true ); }

    /** @return compatibility alias for dofPerVertex(). */
    [[nodiscard]] constexpr uint16_type runtimeDofPerVertex() const noexcept { return 0; }
    /** @return compatibility alias for dofPerEdge(). */
    [[nodiscard]] constexpr uint16_type runtimeDofPerEdge() const noexcept { return 0; }
    /** @return compatibility alias for dofPerFace(). */
    [[nodiscard]] constexpr uint16_type runtimeDofPerFace() const noexcept { return 0; }
    /** @return compatibility alias for dofPerVolume(). */
    [[nodiscard]] uint16_type runtimeDofPerVolume() const noexcept { return dofPerVolume(); }

    /** @return identity because cell-local DoFs are never shared. */
    template<class Element>
    [[nodiscard]] DofTransform dofTransform( Element const&, uint16_type ) const
    {
        return DofTransform::identity();
    }
};
} // namespace fem

/**
 * @brief Descriptor modifier producing a broken version of any FE family.
 * @tparam Family wrapped public finite-element family descriptor
 */
template<class Family>
class Broken : public Family
{
public:
    using wrapped_family_type = Family; /**< Wrapped family descriptor. */
    using component_basis_type = Broken<typename Family::component_basis_type>; /**< Broken scalar component family. */
    static constexpr bool is_broken = true; /**< Identifies the descriptor modifier. */

    /** @brief Bind the wrapped descriptor, then decorate its concrete FE. */
    template<uint16_type N, uint16_type R = N, typename T = double, typename Convex = Simplex<N>>
    struct apply
    {
        using wrapped_type = typename Family::template apply<N, R, T, Convex>::type;
        using type = fem::BrokenElement<wrapped_type>;
        using result_type = type;
    };

    /** @brief Preserve brokenness while changing a legacy family tag. */
    template<uint16_type NewTag>
    struct ChangeTag
    {
        using type = Broken<typename Family::template ChangeTag<NewTag>::type>;
    };
};

/** @brief True when a descriptor or concrete FE carries broken assembly topology. */
template<class T>
inline constexpr bool is_broken_v = []
{
    if constexpr ( requires { T::is_broken; } )
        return static_cast<bool>( T::is_broken );
    else
        return false;
}();

} // namespace Feel
#endif
