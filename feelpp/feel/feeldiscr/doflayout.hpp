/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2026-05-09

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef FEELPP_FEELDISCR_DOFLAYOUT_HPP
#define FEELPP_FEELDISCR_DOFLAYOUT_HPP 1

#include <concepts>
#include <limits>
#include <span>
#include <tuple>
#include <vector>

#include <feel/feelcore/feeltypes.hpp>

namespace Feel
{

enum class DofEntityKind : uint8_type
{
    Vertex = 0,
    Edge = 1,
    Facet = 2,
    Cell = 3,
    Unknown = std::numeric_limits<uint8_type>::max()
};

enum class DofFunctionalKind : uint16_type
{
    PointValue = 0,
    ComponentPointValue = 1,
    NormalMoment = 2,
    TangentialMoment = 3,
    InteriorMoment = 4,
    Modal = 5,
    Other = std::numeric_limits<uint16_type>::max()
};

enum class DofTransformKind : uint8_type
{
    Identity = 0,
    Sign = 1,
    Permutation = 2,
    Matrix = 3
};

struct DofEntityAttachment
{
    static constexpr uint16_type invalid_id = std::numeric_limits<uint16_type>::max();

    DofEntityKind kind = DofEntityKind::Unknown;
    int8_type topologicalDim = -1;
    uint16_type localEntity = invalid_id;
    uint16_type ordinal = invalid_id;
    bool shared = false;

    [[nodiscard]] bool isValid() const noexcept
    {
        return kind != DofEntityKind::Unknown &&
               topologicalDim >= 0 &&
               localEntity != invalid_id &&
               ordinal != invalid_id;
    }
};

struct LocalDofDescriptor
{
    uint16_type localDof = 0;
    uint16_type parentLocalDof = 0;
    uint16_type component = 0;
    DofEntityAttachment attachment;
    DofFunctionalKind functional = DofFunctionalKind::Other;
};

template<typename SizeT = uint32_type>
struct DofKey
{
    using size_type = SizeT;
    static constexpr size_type invalid_id = std::numeric_limits<size_type>::max();

    uint8_type topologicalDim = std::numeric_limits<uint8_type>::max();
    rank_type canonicalPartition = 0;
    size_type canonicalEntityId = invalid_id;
    uint16_type ordinal = DofEntityAttachment::invalid_id;
    uint16_type component = 0;
    uint16_type familyTag = 0;
    uint16_type variant = 0;

    [[nodiscard]] bool isValid() const noexcept
    {
        return canonicalEntityId != invalid_id &&
               ordinal != DofEntityAttachment::invalid_id &&
               topologicalDim != std::numeric_limits<uint8_type>::max();
    }

    friend bool operator==( DofKey const& a, DofKey const& b ) = default;

    friend bool operator<( DofKey const& a, DofKey const& b ) noexcept
    {
        return std::tie( a.topologicalDim, a.canonicalPartition, a.canonicalEntityId, a.ordinal,
                         a.component, a.familyTag, a.variant ) <
               std::tie( b.topologicalDim, b.canonicalPartition, b.canonicalEntityId, b.ordinal,
                         b.component, b.familyTag, b.variant );
    }
};

struct DofTransform
{
    DofTransformKind kind = DofTransformKind::Identity;
    int16_type sign = 1;
    std::vector<uint16_type> permutation = {};

    [[nodiscard]] bool isIdentity() const noexcept
    {
        return kind == DofTransformKind::Identity ||
               ( kind == DofTransformKind::Sign && sign == 1 );
    }

    [[nodiscard]] std::span<uint16_type const> permutationView() const noexcept
    {
        return permutation;
    }

    [[nodiscard]] static DofTransform identity() noexcept
    {
        return {};
    }

    [[nodiscard]] static DofTransform signedOrientation( int16_type orientation ) noexcept
    {
        return DofTransform{ .kind = DofTransformKind::Sign,
                             .sign = orientation };
    }
};

template<class Attachment>
concept LegacyDofAttachment =
    requires( Attachment const& attachment )
{
    { attachment.entityDim } -> std::convertible_to<int>;
    { attachment.entityId } -> std::convertible_to<uint16_type>;
    { attachment.ordinal } -> std::convertible_to<uint16_type>;
    { attachment.kind } -> std::convertible_to<uint16_type>;
    { attachment.isValid() } -> std::convertible_to<bool>;
};

template<class FE>
concept FiniteElementDofCardinalityProvider =
    requires( FE const& fe, uint16_type entityDim, uint16_type entityId, uint16_type facetId, bool perComponent )
{
    { fe.localDofCount( perComponent ) } -> std::convertible_to<uint16_type>;
    { fe.localDofCountOnEntity( entityDim, entityId, perComponent ) } -> std::convertible_to<uint16_type>;
    { fe.localDofCountOnFacet( facetId, perComponent ) } -> std::convertible_to<uint16_type>;
};

template<class FE>
concept FiniteElementDofLayoutProvider =
    FiniteElementDofCardinalityProvider<FE> &&
    requires( FE const& fe, uint16_type localDof, uint16_type component )
{
    { fe.localDofPerComponent() } -> std::convertible_to<uint16_type>;
    { fe.localDofId( localDof, component ) } -> std::convertible_to<uint16_type>;
    { fe.localDofLayout( localDof ).localDofId } -> std::convertible_to<uint16_type>;
    { fe.localDofLayout( localDof ).parentLocalDofId } -> std::convertible_to<uint16_type>;
    { fe.localDofLayout( localDof ).component } -> std::convertible_to<uint16_type>;
    { fe.localDofLayout( localDof ).attachment.entityDim } -> std::convertible_to<int>;
    { fe.localDofLayout( localDof ).attachment.entityId } -> std::convertible_to<uint16_type>;
    { fe.localDofLayout( localDof ).attachment.ordinal } -> std::convertible_to<uint16_type>;
    { fe.localDofLayout( localDof ).attachment.kind } -> std::convertible_to<uint16_type>;
    { fe.localDofLayout( localDof ).attachment.isValid() } -> std::convertible_to<bool>;
};

template<class FE>
concept FiniteElementDofRepresentativePointProvider =
    requires( FE const& fe, uint16_type localDof )
{
    { fe.dofHasRepresentativePoint( localDof ) } -> std::convertible_to<bool>;
    { fe.dofRepresentativePointIndex( localDof ) } -> std::convertible_to<uint16_type>;
};

template<class FE>
concept FiniteElementDofFunctionalProvider =
    requires( FE const& fe, uint16_type localDof )
{
    { fe.dofFunctionalKind( localDof ) } -> std::convertible_to<uint16_type>;
};

template<class FE, class Element>
concept FiniteElementDofTransformProvider =
    requires( FE const& fe, Element const& element, uint16_type localDof )
{
    { fe.dofTransform( element, localDof ) } -> std::same_as<DofTransform>;
};

template<class FE, class Element>
concept FiniteElementDofKeyProvider =
    requires( FE const& fe, Element const& element, uint16_type localDof )
{
    { fe.dofKey( element, localDof ) };
    requires FiniteElementDofTransformProvider<FE, Element>;
};

template<class FE>
concept FiniteElementReferenceTopologyProvider =
    requires
{
    typename FE::reference_convex_type;
    { FE::nDim } -> std::convertible_to<uint16_type>;
    { FE::reference_convex_type::numVertices } -> std::convertible_to<uint16_type>;
};

template<class FE>
concept FiniteElementRuntimeSizedDofProvider =
    FiniteElementReferenceTopologyProvider<FE> &&
    FiniteElementDofCardinalityProvider<FE> &&
    FiniteElementDofLayoutProvider<FE> &&
    FiniteElementDofRepresentativePointProvider<FE> &&
    FiniteElementDofFunctionalProvider<FE> &&
    requires( FE const& fe )
{
    { fe.order() } -> std::convertible_to<uint16_type>;
    { fe.runtimeOrder() } -> std::convertible_to<uint16_type>;
};

template<FiniteElementDofLayoutProvider FE, class Element>
[[nodiscard]] DofTransform
finiteElementEntityOrientationTransform( FE const& fe, Element const& element, uint16_type localDof )
{
    auto const layout = fe.localDofLayout( localDof );
    auto const& attachment = layout.attachment;
    if ( !attachment.isValid() )
        return DofTransform::identity();

    switch ( attachment.entityDim )
    {
    case 1:
    {
        if constexpr ( Element::nDim > 1 )
        {
            using edge_permutation_type = typename Element::edge_permutation_type;
            auto const edgePermutation = element.edgePermutation( attachment.entityId );
            if ( edgePermutation.value() == edge_permutation_type::REVERSE_PERMUTATION )
                return DofTransform::signedOrientation( -1 );
        }
        break;
    }
    case 2:
    {
        if constexpr ( Element::nDim == 3 )
        {
            auto const& face = element.face( attachment.entityId );
            if ( face.ad_first() != element.id() )
                return DofTransform::signedOrientation( -1 );
        }
        break;
    }
    default:
        break;
    }
    return DofTransform::identity();
}

[[nodiscard]] inline DofEntityKind
dofEntityKindFromTopologicalDim( int8_type entityDim, uint16_type referenceDim ) noexcept
{
    if ( entityDim < 0 )
        return DofEntityKind::Unknown;
    if ( entityDim == 0 )
        return DofEntityKind::Vertex;
    if ( static_cast<uint16_type>( entityDim ) == referenceDim )
        return DofEntityKind::Cell;
    if ( static_cast<uint16_type>( entityDim + 1 ) == referenceDim )
        return DofEntityKind::Facet;
    if ( entityDim == 1 )
        return DofEntityKind::Edge;
    return DofEntityKind::Unknown;
}

template<class FE, LegacyDofAttachment Attachment>
[[nodiscard]] DofEntityAttachment
makeDofEntityAttachment( FE const&, Attachment const& attachment )
{
    if ( !attachment.isValid() )
        return {};

    auto const entityDim = static_cast<int8_type>( attachment.entityDim );
    return DofEntityAttachment{
        .kind = dofEntityKindFromTopologicalDim( entityDim, FE::nDim ),
        .topologicalDim = entityDim,
        .localEntity = static_cast<uint16_type>( attachment.entityId ),
        .ordinal = static_cast<uint16_type>( attachment.ordinal ),
        .shared = static_cast<uint16_type>( entityDim ) < FE::nDim
    };
}

template<FiniteElementDofLayoutProvider FE>
[[nodiscard]] uint16_type
finiteElementLocalDofCount( FE const& fe )
{
    return fe.localDofCount( false );
}

template<FiniteElementDofLayoutProvider FE>
[[nodiscard]] DofFunctionalKind
finiteElementDofFunctionalKind( FE const& fe, uint16_type localDof )
{
    if constexpr ( FiniteElementDofFunctionalProvider<FE> )
        return static_cast<DofFunctionalKind>( fe.dofFunctionalKind( localDof ) );
    else
    {
        auto const layout = fe.localDofLayout( localDof );
        return static_cast<DofFunctionalKind>( layout.attachment.kind );
    }
}

template<FiniteElementDofLayoutProvider FE>
[[nodiscard]] LocalDofDescriptor
makeLocalDofDescriptor( FE const& fe, uint16_type localDof )
{
    auto const layout = fe.localDofLayout( localDof );
    return LocalDofDescriptor{
        .localDof = static_cast<uint16_type>( layout.localDofId ),
        .parentLocalDof = static_cast<uint16_type>( layout.parentLocalDofId ),
        .component = static_cast<uint16_type>( layout.component ),
        .attachment = makeDofEntityAttachment( fe, layout.attachment ),
        .functional = finiteElementDofFunctionalKind( fe, localDof )
    };
}

template<FiniteElementDofLayoutProvider FE>
[[nodiscard]] std::vector<LocalDofDescriptor>
makeLocalDofDescriptors( FE const& fe )
{
    std::vector<LocalDofDescriptor> descriptors;
    auto const nLocalDof = finiteElementLocalDofCount( fe );
    descriptors.reserve( nLocalDof );
    for ( uint16_type localDof = 0; localDof < nLocalDof; ++localDof )
        descriptors.push_back( makeLocalDofDescriptor( fe, localDof ) );
    return descriptors;
}

enum class DofLayoutValidationFailure : uint8_type
{
    None = 0,
    EmptyLocalDofCount,
    InvalidLocalDofPerComponent,
    InvalidLocalDofCount,
    InvalidLocalDofId,
    DuplicateLocalDofId,
    MissingLocalDofId,
    InvalidParentLocalDofId,
    InvalidComponent,
    InvalidLocalDofIdMapping,
    InvalidAttachment,
    InvalidEntityDim,
    InvalidEntityId,
    InvalidOrdinal
};

struct DofLayoutValidationResult
{
    bool valid = true;
    DofLayoutValidationFailure failure = DofLayoutValidationFailure::None;
    uint16_type localDof = 0;
    int8_type entityDim = -1;
    uint16_type localEntity = DofEntityAttachment::invalid_id;
    uint16_type ordinal = DofEntityAttachment::invalid_id;

    [[nodiscard]] explicit operator bool() const noexcept { return valid; }
};

template<FiniteElementReferenceTopologyProvider FE>
[[nodiscard]] constexpr uint16_type
finiteElementReferenceEntityCount( uint16_type topologicalDim ) noexcept
{
    if ( topologicalDim == 0 )
        return FE::reference_convex_type::numVertices;
    if ( topologicalDim == FE::nDim )
        return 1;
    if ( topologicalDim == 1 )
    {
        if constexpr ( requires { FE::reference_convex_type::numEdges; } )
            return FE::reference_convex_type::numEdges;
        else
            return 0;
    }
    if ( topologicalDim == 2 )
    {
        if constexpr ( requires { FE::reference_convex_type::numFaces; } )
            return FE::reference_convex_type::numFaces;
        else
            return 0;
    }
    if ( topologicalDim == 3 )
    {
        if constexpr ( FE::nDim == 3 )
            return 1;
        else
            return 0;
    }
    return 0;
}

template<FiniteElementDofLayoutProvider FE>
    requires FiniteElementReferenceTopologyProvider<FE>
[[nodiscard]] DofLayoutValidationResult
validateFiniteElementDofLayout( FE const& fe )
{
    auto makeFailure = []( DofLayoutValidationFailure failure,
                           uint16_type localDof = 0,
                           int8_type entityDim = -1,
                           uint16_type localEntity = DofEntityAttachment::invalid_id,
                           uint16_type ordinal = DofEntityAttachment::invalid_id )
    {
        return DofLayoutValidationResult{ .valid = false,
                                          .failure = failure,
                                          .localDof = localDof,
                                          .entityDim = entityDim,
                                          .localEntity = localEntity,
                                          .ordinal = ordinal };
    };

    auto const nLocalDof = finiteElementLocalDofCount( fe );
    auto const nParentLocalDof = fe.localDofPerComponent();

    if ( nLocalDof == 0 )
        return makeFailure( DofLayoutValidationFailure::EmptyLocalDofCount );
    if ( nParentLocalDof == 0 )
        return makeFailure( DofLayoutValidationFailure::InvalidLocalDofPerComponent );
    if ( ( nLocalDof % nParentLocalDof ) != 0 )
        return makeFailure( DofLayoutValidationFailure::InvalidLocalDofCount );

    auto const nComponents = static_cast<uint16_type>( nLocalDof / nParentLocalDof );
    std::vector<bool> seen( nLocalDof, false );

    for ( uint16_type localDof = 0; localDof < nLocalDof; ++localDof )
    {
        auto const layout = fe.localDofLayout( localDof );
        auto const& attachment = layout.attachment;
        if ( layout.localDofId >= nLocalDof )
            return makeFailure( DofLayoutValidationFailure::InvalidLocalDofId, localDof );
        if ( seen[layout.localDofId] )
            return makeFailure( DofLayoutValidationFailure::DuplicateLocalDofId, layout.localDofId );
        seen[layout.localDofId] = true;

        if ( layout.parentLocalDofId >= nParentLocalDof )
            return makeFailure( DofLayoutValidationFailure::InvalidParentLocalDofId, localDof );
        if ( layout.component >= nComponents )
            return makeFailure( DofLayoutValidationFailure::InvalidComponent, localDof );
        if ( fe.localDofId( layout.parentLocalDofId, layout.component ) != layout.localDofId )
            return makeFailure( DofLayoutValidationFailure::InvalidLocalDofIdMapping, localDof );

        if ( !attachment.isValid() )
            return makeFailure( DofLayoutValidationFailure::InvalidAttachment, localDof );

        auto const entityDim = static_cast<int8_type>( attachment.entityDim );
        auto const localEntity = static_cast<uint16_type>( attachment.entityId );
        auto const ordinal = static_cast<uint16_type>( attachment.ordinal );
        if ( entityDim < 0 || static_cast<uint16_type>( entityDim ) > FE::nDim )
            return makeFailure( DofLayoutValidationFailure::InvalidEntityDim, localDof, entityDim, localEntity, ordinal );

        auto const nEntities = finiteElementReferenceEntityCount<FE>( static_cast<uint16_type>( entityDim ) );
        if ( nEntities == 0 || localEntity >= nEntities )
            return makeFailure( DofLayoutValidationFailure::InvalidEntityId, localDof, entityDim, localEntity, ordinal );

        if constexpr ( FiniteElementDofCardinalityProvider<FE> )
        {
            auto const nEntityDof = fe.localDofCountOnEntity( static_cast<uint16_type>( entityDim ), localEntity, true );
            if ( nEntityDof == 0 || ordinal >= nEntityDof )
                return makeFailure( DofLayoutValidationFailure::InvalidOrdinal, localDof, entityDim, localEntity, ordinal );
        }
    }

    for ( uint16_type localDof = 0; localDof < nLocalDof; ++localDof )
        if ( !seen[localDof] )
            return makeFailure( DofLayoutValidationFailure::MissingLocalDofId, localDof );

    return {};
}

template<FiniteElementDofLayoutProvider FE>
    requires FiniteElementReferenceTopologyProvider<FE>
[[nodiscard]] bool
finiteElementDofLayoutIsComplete( FE const& fe )
{
    return static_cast<bool>( validateFiniteElementDofLayout( fe ) );
}

} // namespace Feel

#endif // FEELPP_FEELDISCR_DOFLAYOUT_HPP
