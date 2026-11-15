/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/** @file pointsetpolicy.hpp
 * @brief Describe the actual mathematics, support and provenance of node policies.
 */
#pragma once
#include <string_view>
#include <type_traits>

namespace Feel
{
template <class, int, class>
class PointSetEquiSpaced;
template <class, int, class>
class PointSetFekete;
template <class, int, class>
class PointSetWarpBlend;
template <class, int, class>
class PointSetElectrostatic;
template <class, int, class>
class PointSetGaussLobatto;
template <class, int, class>
class PointSetRecursiveLGL;

/** @brief Actual node construction, independent of historical policy aliases. */
enum class PointSetConstruction
{
    Equispaced,         /**< Rational simplex/product lattice. */
    TensorLGL,          /**< One-dimensional LGL factors on a product cell. */
    TriangularTable,    /**< Legacy triangular Fekete interior tables. */
    WarpBlend,          /**< Warburton-style explicit warp-and-blend construction. */
    ElectrostaticTable, /**< Legacy triangular electrostatic orbit tables. */
    RecursiveLGL,       /**< Isaac recursion from one-dimensional LGL factors. */
    Custom              /**< A user-supplied policy owns its support contract. */
};

/** @brief Public support/provenance description without renaming legacy policies. */
struct PointSetPolicyInfo
{
    /** Requested policy name, distinct from its actual backend. */
    std::string_view name;
    /** Construction selected for this reference cell. */
    PointSetConstruction construction;
    /** Human-readable dimension and degree envelope. */
    std::string_view coverage;
    /** Source of the algorithm or existing coordinate artifacts. */
    std::string_view provenance;
    /** Whether the family accepts a runtime degree without selecting another family. */
    bool runtimeDegree;
    /** Whether full boundary nodes are provided; interior-only selections remain DG data. */
    bool boundaryNodes;
};

/** @brief Inspect a point policy's actual construction for one reference-cell type.
 * @note Metadata does not certify external data or silently select a fallback.
 * Legacy coordinates/names are preserved, including electrostatic's old "fekete" name.
 */
template <template <class, int, class> class Policy, typename Convex>
[[nodiscard]] constexpr PointSetPolicyInfo pointSetPolicy()
{
    using Selected = Policy<Convex, 1, double>;
    if constexpr ( std::is_same_v<Selected, PointSetEquiSpaced<Convex, 1, double>> )
        return { "equispaced", PointSetConstruction::Equispaced, "simplex/product dimensions 1-3; checked nodal cardinality", "Feel++ rational lattice", true, true };
    else if constexpr ( std::is_same_v<Selected, PointSetRecursiveLGL<Convex, 1, double>> )
        return { "recursive-lgl", PointSetConstruction::RecursiveLGL, "simplex dimensions 1-3; degrees 0-32", "Isaac (2020), DOI 10.1137/20M1321802", true, true };
    else if constexpr ( std::is_same_v<Selected, PointSetFekete<Convex, 1, double>> )
    {
        if constexpr ( Convex::is_hypercube || Convex::nDim == 1 )
            return { "fekete", PointSetConstruction::TensorLGL, "product cells/segments; runtime degrees 0-15", "Feel++ Gauss-Lobatto backend", true, true };
        else if constexpr ( Convex::nDim == 3 )
            return { "fekete", PointSetConstruction::WarpBlend, "tetrahedra; runtime degrees 0-15; static tuned degrees 0-15 then existing asymptotic parameter", "Feel++ warp-and-blend backend", true, true };
        else
            return { "fekete", PointSetConstruction::TriangularTable, "triangles; static degrees 0-19, runtime 0-15; base nodes at 0-2", "Legacy Feel++ tables; boundary from warp-and-blend", true, true };
    }
    else if constexpr ( std::is_same_v<Selected, PointSetWarpBlend<Convex, 1, double>> )
        return { "warpblend", PointSetConstruction::WarpBlend, "simplex dimensions 1-3; static tuned degrees 0-15 then existing asymptotic parameter", "Warburton construction in Feel++", false, true };
    else if constexpr ( std::is_same_v<Selected, PointSetElectrostatic<Convex, 1, double>> )
        return { "electrostatic", PointSetConstruction::ElectrostaticTable, "simplex base nodes at degrees 0-2; triangular static interior tables at 3-7", "Legacy Feel++ electrostatic orbit coordinates", false, true };
    else if constexpr ( std::is_same_v<Selected, PointSetGaussLobatto<Convex, 1, double>> )
        return { "gausslobatto", PointSetConstruction::TensorLGL, "product cells/segments; static degrees subject to nodal cardinality", "Existing Feel++ Gauss-Lobatto-Jacobi root construction", false, true };
    else
        return { "custom", PointSetConstruction::Custom, "defined by the selected policy", "defined by the selected policy", false, false };
}
} // namespace Feel
