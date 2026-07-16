/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#ifndef FEELPP_FEELDISCR_DOFLAYOUT_HPP
#define FEELPP_FEELDISCR_DOFLAYOUT_HPP 1

#include <concepts>
#include <limits>
#include <optional>
#include <memory>
#include <span>
#include <tuple>
#include <vector>

#include <feel/feelcore/feeltypes.hpp>
#include <feel/feelpoly/hcurlpolynomialset.hpp>
#include <feel/feelpoly/hdivpolynomialset.hpp>

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
    SignedPermutation = 3,
    Matrix = 4
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
    /** @brief Shared row-major storage for a dense entity-block operator. */
    struct MatrixStorage
    {
        uint16_type rows = 0;
        uint16_type cols = 0;
        std::vector<double> values;

        [[nodiscard]] double operator()( uint16_type row, uint16_type col ) const noexcept
        {
            return values[static_cast<std::size_t>( row ) * cols + col];
        }
    };

    DofTransformKind kind = DofTransformKind::Identity;
    int16_type sign = 1;
    std::vector<uint16_type> permutation = {};
    std::shared_ptr<MatrixStorage const> matrix = {};
    uint16_type matrixRow = 0;

    [[nodiscard]] bool isIdentity() const noexcept
    {
        return kind == DofTransformKind::Identity ||
               ( kind == DofTransformKind::Sign && sign == 1 );
    }

    /** @return canonical entity ordinal corresponding to a local ordinal. */
    [[nodiscard]] uint16_type mapOrdinal( uint16_type ordinal ) const noexcept
    {
        return ordinal < permutation.size() ? permutation[ordinal] : ordinal;
    }

    [[nodiscard]] std::span<uint16_type const> permutationView() const noexcept
    {
        return permutation;
    }

    /** @return one row of the dense entity operator, or an empty span. */
    [[nodiscard]] std::span<double const> matrixRowView() const noexcept
    {
        if ( !matrix || matrixRow >= matrix->rows )
            return {};
        return { matrix->values.data() + static_cast<std::size_t>( matrixRow ) * matrix->cols,
                 matrix->cols };
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

    /** @brief Build a possibly signed entity permutation operator. */
    [[nodiscard]] static DofTransform signedPermutation( int16_type orientation,
                                                         std::vector<uint16_type> permutation )
    {
        const auto kind = orientation == 1 ? DofTransformKind::Permutation
                                           : DofTransformKind::SignedPermutation;
        return DofTransform{ .kind = kind,
                             .sign = orientation,
                             .permutation = std::move( permutation ) };
    }

    /** @brief Build one row descriptor of a shared dense entity operator. */
    [[nodiscard]] static DofTransform matrixOperator( std::shared_ptr<MatrixStorage const> storage,
                                                      uint16_type row )
    {
        return DofTransform{ .kind = DofTransformKind::Matrix,
                             .matrix = std::move( storage ),
                             .matrixRow = row };
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

/**
 * @brief Minimal mathematical finite-element family contract used by Feel++.
 *
 * DoF topology requirements are composed from the canonical provider concepts
 * above; they are deliberately not repeated here.
 */
template<class FE>
concept FiniteElementFamily =
    FiniteElementReferenceTopologyProvider<FE> &&
    FiniteElementDofLayoutProvider<FE> &&
    requires( FE const& fe )
{
    typename FE::value_type;
    typename FE::primal_space_type;
    typename FE::dual_space_type;
    typename FE::polyset_type;
    typename FE::continuity_type;
    { FE::nComponents } -> std::convertible_to<uint16_type>;
    { FE::is_scalar } -> std::convertible_to<bool>;
    { FE::is_vectorial } -> std::convertible_to<bool>;
    { fe.order() } -> std::convertible_to<uint16_type>;
    { fe.polynomialDegree() } -> std::convertible_to<uint16_type>;
    { fe.localDof() } -> std::convertible_to<uint16_type>;
    { fe.dofPerEntity( uint16_type{}, uint16_type{} ) } -> std::convertible_to<uint16_type>;
};

/** @brief Ciarlet FE exposing its primal space and ordered dual set. */
template<class FE>
concept CiarletFiniteElement =
    FiniteElementFamily<FE> &&
    requires( FE const& fe )
{
    { fe.primal() } -> std::same_as<typename FE::primal_space_type const&>;
    { fe.dual() } -> std::same_as<typename FE::dual_space_type const&>;
    fe.dual()( fe.primal() );
};

/**
 * @brief Local reference finite element consumable by FunctionSpace.
 *
 * This composed contract intentionally reuses the canonical Ciarlet and DoF
 * provider concepts.  The additional requirements are the compatibility
 * facade currently consumed by FunctionSpace and basis contexts; they can be
 * retired individually as those consumers move to dedicated customization
 * points.
 */
template<class FE>
concept ReferenceFiniteElement =
    CiarletFiniteElement<FE> &&
    FiniteElementRuntimeSizedDofProvider<FE> &&
    requires( FE const& fe )
{
    typename FE::SSpace;
    typename FE::PreCompute;
    typename FE::points_type;
    { FE::rank } -> std::convertible_to<uint16_type>;
    { FE::nComponents1 } -> std::convertible_to<uint16_type>;
    { FE::nComponents2 } -> std::convertible_to<uint16_type>;
    { FE::isContinuous } -> std::convertible_to<bool>;
    { FE::is_modal } -> std::convertible_to<bool>;
    { FE::is_product } -> std::convertible_to<bool>;
    { fe.points() } -> std::same_as<typename FE::points_type const&>;
    fe.localInterpolant();
};

template<class FE>
inline constexpr bool finiteElementIsNonconforming =
    []
    {
        if constexpr ( requires { FE::is_nonconforming; } )
            return static_cast<bool>( FE::is_nonconforming );
        else
            return false;
    }();

/** @brief True when FE assembly topology is explicitly cell-local. */
template<class FE>
inline constexpr bool finiteElementIsBroken =
    []
    {
        if constexpr ( requires { FE::is_broken; } )
            return static_cast<bool>( FE::is_broken );
        else
            return false;
    }();

/** @brief Discontinuous L2 finite element. */
template<class FE>
concept L2FiniteElement =
    CiarletFiniteElement<FE> &&
    FE::continuity_type::is_discontinuous_totally;

/** @brief Conforming H1 finite element. */
template<class FE>
concept H1FiniteElement =
    CiarletFiniteElement<FE> &&
    FE::isContinuous &&
    !finiteElementIsNonconforming<FE> &&
    !is_hdiv_conforming_v<FE> &&
    !is_hcurl_conforming_v<FE>;

/** @brief Broken/nonconforming H1 finite element (for example CR). */
template<class FE>
concept NonconformingH1FiniteElement =
    CiarletFiniteElement<FE> && finiteElementIsNonconforming<FE>;

/** @brief H(curl)-conforming finite element. */
template<class FE>
concept HCurlFiniteElement =
    CiarletFiniteElement<FE> && FE::is_vectorial && is_hcurl_conforming_v<FE> &&
    !finiteElementIsBroken<FE>;

/** @brief H(div)-conforming finite element. */
template<class FE>
concept HDivFiniteElement =
    CiarletFiniteElement<FE> && FE::is_vectorial && is_hdiv_conforming_v<FE> &&
    !finiteElementIsBroken<FE>;

namespace detail
{
[[nodiscard]] inline std::vector<uint16_type>
composeDofPermutations( std::vector<uint16_type> const& lhs,
                        std::vector<uint16_type> const& rhs )
{
    if ( lhs.size() != rhs.size() )
        return {};
    std::vector<uint16_type> result( lhs.size() );
    for ( uint16_type i = 0; i < lhs.size(); ++i )
        result[i] = rhs[lhs[i]];
    return result;
}

template<class FacePermutation>
[[nodiscard]] std::optional<std::vector<uint16_type>>
simplexTriangularDofPermutation( FacePermutation permutation, uint16_type nDof )
{
    uint16_type nSide = 0;
    while ( nSide * ( nSide + 1 ) / 2 < nDof )
        ++nSide;
    if ( nSide * ( nSide + 1 ) / 2 != nDof )
        return std::nullopt;

    std::vector<uint16_type> reverseHypotenuse( nDof );
    std::vector<uint16_type> reverseBase( nDof );
    for ( uint16_type i = 0; i < nSide; ++i )
    {
        for ( uint16_type j = 0; j <= i; ++j )
        {
            const uint16_type first = i + j * ( j - 1 ) / 2 + j * ( nSide - j );
            const uint16_type last = i + ( i - j ) * ( i - j - 1 ) / 2 + ( i - j ) * ( nSide - ( i - j ) );
            reverseHypotenuse[first] = last;
        }
        const uint16_type begin = i * nSide - i * ( i - 1 ) / 2;
        const uint16_type end = ( i + 1 ) * nSide - ( i + 1 ) * i / 2 - 1;
        for ( uint16_type j = 0; j <= nSide - i - 1; ++j )
            reverseBase[begin + j] = end - j;
    }

    if ( permutation == FacePermutation( FacePermutation::IDENTITY ) )
    {
        std::vector<uint16_type> identity( nDof );
        for ( uint16_type i = 0; i < nDof; ++i )
            identity[i] = i;
        return identity;
    }
    if ( permutation == FacePermutation( FacePermutation::REVERSE_HYPOTENUSE ) )
        return reverseHypotenuse;
    if ( permutation == FacePermutation( FacePermutation::REVERSE_BASE ) )
        return reverseBase;

    auto const rotationAntiClock = composeDofPermutations( reverseBase, reverseHypotenuse );
    auto const rotationClockwise = composeDofPermutations( reverseHypotenuse, reverseBase );
    if ( permutation == FacePermutation( FacePermutation::ROTATION_ANTICLOCK ) )
        return rotationAntiClock;
    if ( permutation == FacePermutation( FacePermutation::ROTATION_CLOCKWISE ) )
        return rotationClockwise;
    if ( permutation == FacePermutation( FacePermutation::REVERSE_HEIGHT ) )
        return composeDofPermutations( reverseBase, rotationClockwise );
    return std::nullopt;
}
} // namespace detail

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
            {
                const auto nEntityDof = fe.localDofCountOnEntity( 1, attachment.entityId, true );
                std::vector<uint16_type> permutation( nEntityDof );
                for ( uint16_type ordinal = 0; ordinal < nEntityDof; ++ordinal )
                    permutation[ordinal] = static_cast<uint16_type>( nEntityDof - 1 - ordinal );
                return DofTransform::signedPermutation( -1, std::move( permutation ) );
            }
        }
        break;
    }
    case 2:
    {
        if constexpr ( Element::nDim == 3 )
        {
            auto const& face = element.face( attachment.entityId );
            const int16_type sign = face.ad_first() == element.id() ? 1 : -1;
            const auto nEntityDof = fe.localDofCountOnEntity( 2, attachment.entityId, true );
            if ( nEntityDof > 1 )
            {
                auto permutation = Feel::detail::simplexTriangularDofPermutation(
                    element.facePermutation( attachment.entityId ), nEntityDof );
                if ( permutation )
                    return DofTransform::signedPermutation( sign, std::move( *permutation ) );
            }
            if ( sign < 0 )
                return DofTransform::signedOrientation( sign );
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
