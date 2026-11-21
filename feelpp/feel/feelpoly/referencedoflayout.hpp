// SPDX-FileCopyrightText: 2005-2026 University of Strasbourg
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <feel/feelpoly/order.hpp>

namespace Feel::detail
{
/**
 * @brief Shared finite-element layout protocol with a cell-owned modal policy.
 * @tparam Derived Reference element exposing dimension, shape and cardinality.
 * @tparam Modal Whether modal functionals belong entirely to the cell.
 *
 * Families override virtual ownership hooks. Component numbering, local layout,
 * entity counts and facet closure use one implementation for nodal and modal
 * elements on simplices and hypercubes.
 */
template<class Derived, bool Modal = false>
class ReferenceDofLayout
{
public:
    /** @brief Destroy a layout through its interface. */
    virtual ~ReferenceDofLayout() = default;
    /** @brief Reference-entity ownership of one local functional. */
    struct DofAttachment
    {
        /** Sentinel for an unknown entity or ordinal. */
        static constexpr uint16_type invalid_id = std::numeric_limits<uint16_type>::max();
        /** Entity dimension; negative means unknown. */
        int8_t entityDim = -1;
        /** Reference entity number. */
        uint16_type entityId = invalid_id;
        /** Functional ordinal within the entity. */
        uint16_type ordinal = invalid_id;
        /** Compatibility functional kind. */
        uint16_type kind = 0;

        /** @brief Whether the functional has a known entity attachment. */
        [[nodiscard]] bool isValid() const { return entityDim >= 0; }
    };

    /** @brief Component and entity coordinates of a local functional. */
    struct LocalDofLayout
    {
        /** Flattened local functional index. */
        uint16_type localDofId = 0;
        /** Index within one component. */
        uint16_type parentLocalDofId = 0;
        /** Field component. */
        uint16_type component = 0;
        /** Entity attachment metadata. */
        DofAttachment attachment;
    };

    /**
     * \return number of local dofs per component.
     *
     * This is FE-owned so consumers do not assume local dof layout.
     */
    virtual uint16_type localDofPerComponent() const
    {
        if constexpr ( Modal ) return derived().localDof();
        else return Derived::nLocalDof;
    }

    /**
     * \return number of local dofs owned by the FE layout.
     *
     * When \p perComponent is false this is the flattened local cardinality
     * exposed to local-to-global tables. When it is true this is the parent
     * local cardinality for one component.
     */
    virtual uint16_type localDofCount( bool perComponent = false ) const
    {
        const uint16_type localDofPerComp = this->localDofPerComponent();
        if ( perComponent )
            return localDofPerComp;

        if constexpr ( Derived::is_product )
            return ::Feel::detail::checkedLocalIndexValue(
                static_cast<size_type>( Derived::nComponents ) * localDofPerComp );
        else
            return localDofPerComp;
    }

    /** \return parent local cardinality for one field component. */
    [[nodiscard]] uint16_type localDof() const
    {
        return this->localDofCount( true );
    }

    /** \return per-component cardinality on one reference vertex. */
    [[nodiscard]] uint16_type dofPerVertex() const
    {
        if constexpr ( Modal ) return this->localDofCountOnEntity( 0, 0, true );
        else return Derived::nDofPerVertex;
    }

    /** \return per-component cardinality on one reference edge. */
    [[nodiscard]] uint16_type dofPerEdge() const
    {
        if constexpr ( Modal ) return this->localDofCountOnEntity( 1, 0, true );
        else return Derived::nDofPerEdge;
    }

    /** \return per-component cardinality on one reference face. */
    [[nodiscard]] uint16_type dofPerFace() const
    {
        if constexpr ( Modal ) return this->localDofCountOnEntity( 2, 0, true );
        else return Derived::nDofPerFace;
    }

    /** \return per-component cardinality on one reference volume. */
    [[nodiscard]] uint16_type dofPerVolume() const
    {
        if constexpr ( Modal ) return this->localDofCountOnEntity( 3, 0, true );
        else return Derived::nDofPerVolume;
    }

    /** \return per-component cardinality on one reference entity. */
    [[nodiscard]] uint16_type dofPerEntity( uint16_type topologicalDim,
                                            uint16_type localEntity = 0 ) const
    {
        return this->localDofCountOnEntity( topologicalDim, localEntity, true );
    }

    /**
     * \return local dof id from (parent local dof id, component)
     *
     * Layout ownership stays in FE implementation.
     */
    virtual uint16_type localDofId( uint16_type parentLocalDofId, uint16_type component = 0 ) const
    {
        if ( parentLocalDofId >= this->localDofPerComponent() || component >= Derived::nComponents )
            throw std::out_of_range("local degree-of-freedom index or component is outside the element");
        if constexpr ( Derived::is_product )
        {
            const uint16_type localDofPerComp = this->localDofPerComponent();
            const auto localDof = checkedLocalIndexValue(
                static_cast<size_type>(component) * localDofPerComp + parentLocalDofId );
            if ( this->dofParent( localDof ) != parentLocalDofId )
                throw std::logic_error("inconsistent component/local degree-of-freedom layout");
            return localDof;
        }
        else
        {
            if ( component != 0 )
                throw std::out_of_range("non-product finite element supports component zero only");
            return parentLocalDofId;
        }
    }

    /**
     * \return entity attachment metadata for a local dof.
     *
     * Default implementation returns unknown entity attachment and keeps
     * dof type as functional kind.
     */
    virtual DofAttachment dofAttachment( uint16_type localDofId ) const
    {
        if constexpr ( Modal )
            return DofAttachment{ .entityDim = Derived::nDim, .entityId = 0,
                                  .ordinal = this->dofParent( localDofId ), .kind = 1 };
        return DofAttachment{ .entityDim = -1,
                              .entityId = DofAttachment::invalid_id,
                              .ordinal = DofAttachment::invalid_id,
                              .kind = this->dofType( localDofId ) };
    }

    /**
     * \return full FE-owned layout descriptor for a local dof.
     */
    LocalDofLayout localDofLayout( uint16_type localDofId ) const
    {
        return LocalDofLayout{ .localDofId = localDofId,
                               .parentLocalDofId = this->dofParent( localDofId ),
                               .component = this->component( localDofId ),
                               .attachment = this->dofAttachment( localDofId ) };
    }

    /**
     * \return number of local dofs attached to a reference entity.
     *
     * The default implementation derives the count from FE-owned local dof
     * layouts. FE families may override when a cheaper/runtime-aware path is
     * available.
     */
    virtual uint16_type localDofCountOnEntity( uint16_type topologicalDim,
                                               uint16_type localEntity,
                                               bool perComponent = false ) const
    {
        if constexpr ( Modal )
            return topologicalDim == Derived::nDim && localEntity == 0
                ? this->localDofCount( perComponent ) : 0;
        uint16_type count = 0;
        const uint16_type nLocalDof = this->localDofCount();
        for ( uint16_type localDof = 0; localDof < nLocalDof; ++localDof )
        {
            auto const layout = this->localDofLayout( localDof );
            auto const& attachment = layout.attachment;
            if ( !attachment.isValid() )
                continue;
            if ( perComponent && layout.component != 0 )
                continue;
            if ( attachment.entityDim == static_cast<int8_t>( topologicalDim ) &&
                 attachment.entityId == localEntity )
                ++count;
        }
        return count;
    }

    /**
     * \return number of local dofs attached to a reference codim-1 facet.
     */
    virtual uint16_type localDofCountOnFacet( uint16_type localFacet = 0,
                                              bool perComponent = false ) const
    {
        if constexpr ( Modal ) return 0;
        if constexpr ( Derived::nDim == 1 )
            return this->localDofCountOnEntity( 0, localFacet, perComponent );
        else
        {
            using convex_type = typename Derived::convex_type;
            using face_type = typename convex_type::topological_face_type;

            auto isInFacetClosure = [localFacet]( DofAttachment const& attachment ) -> bool
            {
                if ( !attachment.isValid() )
                    return false;

                if ( attachment.entityDim == 0 )
                {
                    for ( uint16_type localVertex = 0; localVertex < face_type::numVertices; ++localVertex )
                        if ( static_cast<uint16_type>( convex_type::f2p( localFacet, localVertex ) ) == attachment.entityId )
                            return true;
                    return false;
                }

                if ( attachment.entityDim == 1 )
                {
                    if constexpr ( Derived::nDim == 2 )
                        return attachment.entityId == localFacet;
                    else
                    {
                        for ( uint16_type localEdge = 0; localEdge < face_type::numEdges; ++localEdge )
                            if ( static_cast<uint16_type>( convex_type::f2e( localFacet, localEdge ) ) == attachment.entityId )
                                return true;
                        return false;
                    }
                }

                if ( attachment.entityDim == 2 )
                {
                    if constexpr ( Derived::nDim == 3 )
                        return attachment.entityId == localFacet;
                    else
                        return false;
                }

                return false;
            };

            uint16_type count = 0;
            const uint16_type nLocalDof = this->localDofCount();
            for ( uint16_type localDof = 0; localDof < nLocalDof; ++localDof )
            {
                auto const layout = this->localDofLayout( localDof );
                if ( perComponent && layout.component != 0 )
                    continue;
                if ( isInFacetClosure( layout.attachment ) )
                    ++count;
            }
            return count;
        }
    }

    /**
     * \return true when a local dof has a representative interpolation point.
     *
     * Point-evaluation finite elements keep the historical behavior: the
     * representative point is the FE point indexed by dofParent(localDofId).
     * Moment/modal finite elements should override this method when no such
     * geometric representative exists.
     */
    virtual bool dofHasRepresentativePoint( uint16_type localDofId ) const
    {
        (void)localDofId;
        return !Modal;
    }

    /**
     * \return point index in points() used to represent localDofId.
     *
     * This method is only meaningful when dofHasRepresentativePoint(localDofId)
     * returns true.
     */
    virtual uint16_type dofRepresentativePointIndex( uint16_type localDofId ) const
    {
        return this->dofParent( localDofId );
    }

    /**
     * \return finite-element owned functional kind for a local dof.
     *
     * The default preserves legacy DofAttachment::kind values. FE families with
     * moment functionals can override this without changing their public
     * DofAttachment compatibility value.
     */
    virtual uint16_type dofFunctionalKind( uint16_type localDofId ) const
    {
        return this->dofType( localDofId );
    }

    //! \return the component of a local dof
    virtual uint16_type component( uint16_type localDofId ) const
    {
        validateLocalDof( localDofId );
        return Derived::is_product ? localDofId / this->localDofPerComponent() : 0;
    }

    //! \return a parent local dof id for each component (for example, the first component)
    virtual uint16_type dofParent( uint16_type localDofId ) const
    {
        validateLocalDof( localDofId );
        return Derived::is_product ? localDofId % this->localDofPerComponent() : localDofId;
    }

    //! \return the type of a local dof
    virtual uint16_type dofType( uint16_type localDofId ) const
    {
        validateLocalDof( localDofId );
        return 1;
    }

    //! give an unsymmetric dof index i, provide the symmetric one
    virtual uint16_type unsymmToSymm( uint16_type i ) const { return i; }

private:
    /** @brief Access the reference element's effective cardinality. */
    Derived const& derived() const { return static_cast<Derived const&>(*this); }
    /** @brief Reject an invalid flattened functional index. */
    void validateLocalDof( uint16_type id ) const
    {
        if ( id >= this->localDofCount() )
            throw std::out_of_range("local degree-of-freedom index is outside the element");
    }
};
}
