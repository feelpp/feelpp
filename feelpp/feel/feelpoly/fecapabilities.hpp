/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file fecapabilities.hpp
 * @brief Published finite-element and geometry operation capability matrix.
 */
#if !defined( FEELPP_FEELPOLY_FECAPABILITIES_HPP )
#define FEELPP_FEELPOLY_FECAPABILITIES_HPP 1

#include <feel/feelpoly/geometrycoefficients.hpp>
#include <feel/feelpoly/order.hpp>

#include <cstdint>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace Feel::fe
{
/** Public finite-element family represented by one capability row. */
enum class Family : std::uint8_t
{
    Lagrange,
    RaviartThomas,
    BrezziDouglasMarini,
    NedelecFirstKind,
    NedelecSecondKind
};

/** Operation whose support envelope is being queried. */
enum class Operation : std::uint8_t
{
    SpaceConstruction,
    VolumeAssembly,
    FacetAssembly,
    LinearSolve,
    UnrelatedMeshTransfer,
    CheckpointRestart,
    Export,
    Submesh
};

/** Strength of the published support statement. */
enum class CapabilityLevel : std::uint8_t
{
    Unsupported,
    Available,
    Validated
};

/** Execution context that can change operation support. */
enum class Execution : std::uint8_t
{
    Serial,
    Distributed
};

/** Operand applied to the domain field during unrelated-mesh transfer. */
enum class TransferOperand : std::uint8_t
{
    Identity,
    Differential
};

/** One complete family/cell/order/geometry/operation query. */
struct CapabilityQuery
{
    Family family;
    geometry::GeometryCellKind cell;
    uint16_type dimension;
    uint16_type order;
    uint16_type geometry_order;
    Operation operation;
    Execution execution = Execution::Serial;
    /** Transfer-only metadata; `family` describes the domain FE. */
    TransferOperand transfer_operand = TransferOperand::Identity;
    /** An omitted image family requests transfer within the domain family. */
    std::optional<Family> transfer_image_family = std::nullopt;
};

/** Result and stable diagnostic associated with a capability query. */
struct Capability
{
    CapabilityLevel level;
    std::string_view diagnostic;

    [[nodiscard]] constexpr bool available() const noexcept
    {
        return level != CapabilityLevel::Unsupported;
    }

    [[nodiscard]] constexpr bool validated() const noexcept
    {
        return level == CapabilityLevel::Validated;
    }
};

[[nodiscard]] constexpr std::string_view familyName( Family family ) noexcept
{
    switch ( family )
    {
    case Family::Lagrange: return "Lagrange";
    case Family::RaviartThomas: return "Raviart-Thomas";
    case Family::BrezziDouglasMarini: return "Brezzi-Douglas-Marini";
    case Family::NedelecFirstKind: return "Nedelec first kind";
    case Family::NedelecSecondKind: return "Nedelec second kind";
    }
    return "unknown";
}

[[nodiscard]] constexpr std::string_view operationName(
    Operation operation ) noexcept
{
    switch ( operation )
    {
    case Operation::SpaceConstruction: return "space construction";
    case Operation::VolumeAssembly: return "volume assembly";
    case Operation::FacetAssembly: return "facet assembly";
    case Operation::LinearSolve: return "linear solve";
    case Operation::UnrelatedMeshTransfer: return "unrelated-mesh transfer";
    case Operation::CheckpointRestart: return "checkpoint/restart";
    case Operation::Export: return "export";
    case Operation::Submesh: return "submesh";
    }
    return "unknown";
}

/**
 * @brief Query the centralized validated/available finite-element envelope.
 *
 * `Available` deliberately differs from `Validated`: it permits existing
 * mathematically implemented higher orders without advertising end-to-end PDE
 * validation beyond the tested G1--G3/order-0--3 envelope.
 */
[[nodiscard]] constexpr Capability
capability( CapabilityQuery const& query ) noexcept
{
    bool const lagrange = query.family == Family::Lagrange;
    bool const compatible = !lagrange;
    bool const supportedDimension =
        lagrange ? ( query.dimension >= 1 && query.dimension <= 3 )
                 : ( query.dimension == 2 || query.dimension == 3 );
    if ( !supportedDimension )
        return { CapabilityLevel::Unsupported,
                 "family is unavailable in the requested dimension" };
    if ( compatible &&
         query.cell != geometry::GeometryCellKind::Simplex )
        return { CapabilityLevel::Unsupported,
                 "H(div)/H(curl) tensor-product cells are not implemented" };
    if ( query.geometry_order < 1 )
        return { CapabilityLevel::Unsupported,
                 "geometry order must be positive" };

    bool const validatedPdeEnvelope =
        query.order <= 3 && query.geometry_order <= 3;
    switch ( query.operation )
    {
    case Operation::SpaceConstruction:
    case Operation::VolumeAssembly:
    case Operation::FacetAssembly:
        return validatedPdeEnvelope
                   ? Capability{ CapabilityLevel::Validated, {} }
                   : Capability{ CapabilityLevel::Available,
                                 "available outside the end-to-end validated envelope" };

    case Operation::LinearSolve:
    {
        bool const validatedSolve =
            ( lagrange && query.order == 2 && query.geometry_order <= 3 ) ||
            ( query.family == Family::RaviartThomas &&
              query.order == 1 && query.geometry_order == 2 ) ||
            ( query.family == Family::NedelecFirstKind &&
              query.order == 1 &&
              ( query.geometry_order == 1 || query.geometry_order == 3 ) );
        return validatedSolve
                   ? Capability{ CapabilityLevel::Validated, {} }
                   : Capability{ CapabilityLevel::Available,
                                 "solve is available but lacks end-to-end validation for this combination" };
    }

    case Operation::Submesh:
        return lagrange && query.order == 2 && query.geometry_order == 2
                   ? Capability{ CapabilityLevel::Validated, {} }
                   : Capability{ CapabilityLevel::Available,
                                 "submesh use is available but lacks end-to-end PDE validation for this combination" };

    case Operation::UnrelatedMeshTransfer:
        if ( query.execution == Execution::Distributed )
        {
            if ( compatible ||
                 query.transfer_operand != TransferOperand::Identity )
                return { CapabilityLevel::Unsupported,
                         "distributed unrelated-mesh transfer does not support mapped H(div)/H(curl) or differential operands" };
            if ( query.transfer_image_family.value_or( query.family ) !=
                 Family::Lagrange )
                return { CapabilityLevel::Unsupported,
                         "distributed unrelated-mesh transfer requires a point-functional image FE; moment sample transport is not implemented" };
            return { CapabilityLevel::Available,
                     "distributed nodal identity transfer is available but lacks end-to-end validation for this combination" };
        }
        if ( query.transfer_operand != TransferOperand::Identity ||
             query.transfer_image_family.value_or( query.family ) != query.family )
            return { CapabilityLevel::Available,
                     "transfer operand or family pair is outside the fully validated moment envelope" };
        if ( lagrange && query.order <= 2 && query.geometry_order <= 2 )
            return { CapabilityLevel::Validated, {} };
        if ( ( query.family == Family::RaviartThomas ||
               query.family == Family::NedelecFirstKind ) &&
             query.order <= 1 && query.geometry_order <= 2 )
            return { CapabilityLevel::Validated, {} };
        return { CapabilityLevel::Available,
                 "transfer is available outside the fully validated moment envelope" };

    case Operation::CheckpointRestart:
        if ( lagrange && query.order == 2 && query.geometry_order == 2 )
            return { CapabilityLevel::Validated, {} };
        return { CapabilityLevel::Available,
                 "checkpoint is available but lacks full family/order PDE coverage" };

    case Operation::Export:
        if ( query.geometry_order == 1 && query.order <= 3 )
            return { CapabilityLevel::Validated, {} };
        if ( query.geometry_order <= 3 )
            return { CapabilityLevel::Available,
                     "curved export support depends on the selected format" };
        return { CapabilityLevel::Available,
                 "high-order export is format-dependent above G3" };
    }
    return { CapabilityLevel::Unsupported, "unknown operation" };
}

/** Reject an unsupported query at a public factory or operation boundary. */
inline void requireCapability( CapabilityQuery const& query )
{
    auto const result = capability( query );
    if ( result.available() )
        return;
    std::ostringstream message;
    message << familyName( query.family ) << " "
            << operationName( query.operation ) << " is unsupported on "
            << ( query.cell == geometry::GeometryCellKind::Simplex
                     ? "simplex"
                     : "hypercube" )
            << " dimension " << query.dimension
            << ", FE order " << query.order
            << ", geometry order " << query.geometry_order
            << ": " << result.diagnostic;
    throw std::invalid_argument( message.str() );
}

/** Build and require a query directly from a static- or runtime-order mesh. */
template <Family FamilyValue, typename MeshType>
inline void requireCapability(
    MeshType const& mesh, RuntimeOrder order,
    Operation operation = Operation::SpaceConstruction )
{
    requireCapability( CapabilityQuery{
        FamilyValue,
        MeshType::shape_type::is_simplex
            ? geometry::GeometryCellKind::Simplex
            : geometry::GeometryCellKind::Hypercube,
        MeshType::nDim, order.value, mesh.order(), operation,
        mesh.worldComm().globalSize() > 1
            ? Execution::Distributed
            : Execution::Serial } );
}

} // namespace Feel::fe

#endif // FEELPP_FEELPOLY_FECAPABILITIES_HPP
