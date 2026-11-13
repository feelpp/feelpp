/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_TERMS_HPP
#define FEELPP_FEELSYSTEM_TERMS_HPP 1

/**
 * @file terms.hpp
 * @brief Term roles, metadata, diagnostics, and assembly-planning summaries.
 */

#include <algorithm>
#include <cstddef>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace Feel::systems
{

/**
 * @brief Supported time-scheme families for evolution assembly.
 */
enum class time_scheme_kind
{
    none,    /**< No time scheme is attached. */
    bdf,     /**< First-order BDF time scheme. */
    newmark  /**< Second-order Newmark time scheme. */
};

/**
 * @brief Lightweight metadata used internally by time-step assembly.
 *
 * The context intentionally exposes scalar time information and scheme
 * coefficients only. History fields and restart/save ownership remain in the
 * concrete `feelts` time-scheme object.
 */
struct time_step_context
{
    /** @brief Active time-scheme family. */
    time_scheme_kind scheme = time_scheme_kind::none;
    /** @brief Current target time. */
    double time = 0.0;
    /** @brief Current time-step size. */
    double timeStep = 0.0;
    /** @brief Current time iteration index. */
    int iteration = 0;
    /** @brief First-order derivative coefficient used by BDF-like schemes. */
    double derivativeCoefficient = 0.0;
    /** @brief First derivative coefficient used by second-order schemes. */
    double firstDerivativeCoefficient = 0.0;
    /** @brief Second derivative coefficient used by second-order schemes. */
    double secondDerivativeCoefficient = 0.0;

    /** @brief Return true when the context describes a BDF step. */
    bool isBdf() const noexcept { return scheme == time_scheme_kind::bdf; }
    /** @brief Return true when the context describes a Newmark step. */
    bool isNewmark() const noexcept { return scheme == time_scheme_kind::newmark; }
    /** @brief Return true for first-order evolution contexts. */
    bool isFirstOrder() const noexcept { return this->isBdf(); }
    /** @brief Return true for second-order evolution contexts. */
    bool isSecondOrder() const noexcept { return this->isNewmark(); }
};

/**
 * @brief Storage role used by problem descriptors to route terms.
 */
enum class term_role
{
    lhs,         /**< Stationary lhs or generic matrix contribution. */
    rhs,         /**< Stationary rhs or generic vector contribution. */
    mass,        /**< Mass/time-derivative matrix contribution. */
    damping,     /**< Damping matrix contribution. */
    stiffness,   /**< Stiffness matrix contribution. */
    residual,    /**< Nonlinear residual contribution. */
    jacobian,    /**< Nonlinear Jacobian contribution. */
    iteration,   /**< Picard iteration contribution. */
    operator_a,  /**< Eigenproblem operator A contribution. */
    operator_b,  /**< Eigenproblem operator B contribution. */
    constraint,  /**< Deferred constraint or boundary-condition contribution. */
    time_scheme  /**< Attached time-scheme descriptor. */
};

/** @brief Tag value type selecting `term_role::lhs`. */
struct lhs_t {};
/** @brief Tag value type selecting `term_role::rhs`. */
struct rhs_t {};
/** @brief Tag value type selecting `term_role::mass`. */
struct mass_t {};
/** @brief Tag value type selecting `term_role::damping`. */
struct damping_t {};
/** @brief Tag value type selecting `term_role::stiffness`. */
struct stiffness_t {};
/** @brief Tag value type selecting `term_role::residual`. */
struct residual_t {};
/** @brief Tag value type selecting `term_role::jacobian`. */
struct jacobian_t {};
/** @brief Tag value type selecting `term_role::iteration`. */
struct iteration_t {};
/** @brief Tag value type selecting `term_role::operator_a`. */
struct operator_a_t {};
/** @brief Tag value type selecting `term_role::operator_b`. */
struct operator_b_t {};
/** @brief Tag value type selecting `term_role::constraint`. */
struct constraint_t {};
/** @brief Tag value type selecting `term_role::time_scheme`. */
struct time_scheme_t {};

/** @brief Value object selecting `term_role::lhs`. */
inline constexpr lhs_t lhs{};
/** @brief Value object selecting `term_role::rhs`. */
inline constexpr rhs_t rhs{};
/** @brief Value object selecting `term_role::mass`. */
inline constexpr mass_t mass{};
/** @brief Value object selecting `term_role::damping`. */
inline constexpr damping_t damping{};
/** @brief Value object selecting `term_role::stiffness`. */
inline constexpr stiffness_t stiffness{};
/** @brief Value object selecting `term_role::residual`. */
inline constexpr residual_t residual{};
/** @brief Value object selecting `term_role::jacobian`. */
inline constexpr jacobian_t jacobian{};
/** @brief Value object selecting `term_role::iteration`. */
inline constexpr iteration_t iteration{};
/** @brief Value object selecting `term_role::operator_a`. */
inline constexpr operator_a_t operator_a{};
/** @brief Value object selecting `term_role::operator_b`. */
inline constexpr operator_b_t operator_b{};
/** @brief Value object selecting `term_role::constraint`. */
inline constexpr constraint_t constraint{};
/** @brief Value object selecting `term_role::time_scheme`. */
inline constexpr time_scheme_t time_scheme{};

/** @brief Map `lhs_t` to its storage role. */
constexpr term_role role( lhs_t ) noexcept
{
    return term_role::lhs;
}
/** @brief Map `rhs_t` to its storage role. */
constexpr term_role role( rhs_t ) noexcept
{
    return term_role::rhs;
}
/** @brief Map `mass_t` to its storage role. */
constexpr term_role role( mass_t ) noexcept
{
    return term_role::mass;
}
/** @brief Map `damping_t` to its storage role. */
constexpr term_role role( damping_t ) noexcept
{
    return term_role::damping;
}
/** @brief Map `stiffness_t` to its storage role. */
constexpr term_role role( stiffness_t ) noexcept
{
    return term_role::stiffness;
}
/** @brief Map `residual_t` to its storage role. */
constexpr term_role role( residual_t ) noexcept
{
    return term_role::residual;
}
/** @brief Map `jacobian_t` to its storage role. */
constexpr term_role role( jacobian_t ) noexcept
{
    return term_role::jacobian;
}
/** @brief Map `iteration_t` to its storage role. */
constexpr term_role role( iteration_t ) noexcept
{
    return term_role::iteration;
}
/** @brief Map `operator_a_t` to its storage role. */
constexpr term_role role( operator_a_t ) noexcept
{
    return term_role::operator_a;
}
/** @brief Map `operator_b_t` to its storage role. */
constexpr term_role role( operator_b_t ) noexcept
{
    return term_role::operator_b;
}
/** @brief Map `constraint_t` to its storage role. */
constexpr term_role role( constraint_t ) noexcept
{
    return term_role::constraint;
}
/** @brief Map `time_scheme_t` to its storage role. */
constexpr term_role role( time_scheme_t ) noexcept
{
    return term_role::time_scheme;
}

/**
 * @brief Human-readable name for a term storage role.
 */
constexpr std::string_view
roleName( term_role r ) noexcept
{
    switch ( r )
    {
    case term_role::lhs:
        return "lhs";
    case term_role::rhs:
        return "rhs";
    case term_role::mass:
        return "mass";
    case term_role::damping:
        return "damping";
    case term_role::stiffness:
        return "stiffness";
    case term_role::residual:
        return "residual";
    case term_role::jacobian:
        return "jacobian";
    case term_role::iteration:
        return "iteration";
    case term_role::operator_a:
        return "operator_a";
    case term_role::operator_b:
        return "operator_b";
    case term_role::constraint:
        return "constraint";
    case term_role::time_scheme:
        return "time_scheme";
    }
    return "unknown";
}

/**
 * @brief Assembly stage currently applying a term or constraint.
 */
enum class AssemblyMode
{
    Linear,        /**< Full stationary linear assembly. */
    Affine,        /**< Affine matrix/vector assembly reused by other systems. */
    MatrixTerm,    /**< Assembly of one matrix-like role. */
    VectorTerm,    /**< Assembly of one vector-like role. */
    Residual,      /**< Nonlinear residual assembly. */
    Jacobian,      /**< Nonlinear Jacobian assembly. */
    Picard,        /**< Picard iteration assembly. */
    EvolutionStep, /**< Time-step assembly. */
    Mixed,         /**< Product-space mixed/block assembly. */
    Constraint,    /**< Constraint application. */
    Eigen,         /**< Eigen operator assembly. */
    Unknown        /**< Unspecified assembly mode. */
};

/**
 * @brief Human-readable name for an assembly mode.
 */
constexpr std::string_view
assemblyModeName( AssemblyMode mode ) noexcept
{
    switch ( mode )
    {
    case AssemblyMode::Linear:
        return "linear";
    case AssemblyMode::Affine:
        return "affine";
    case AssemblyMode::MatrixTerm:
        return "matrix_term";
    case AssemblyMode::VectorTerm:
        return "vector_term";
    case AssemblyMode::Residual:
        return "residual";
    case AssemblyMode::Jacobian:
        return "jacobian";
    case AssemblyMode::Picard:
        return "picard";
    case AssemblyMode::EvolutionStep:
        return "evolution_step";
    case AssemblyMode::Mixed:
        return "mixed";
    case AssemblyMode::Constraint:
        return "constraint";
    case AssemblyMode::Eigen:
        return "eigen";
    case AssemblyMode::Unknown:
        return "unknown";
    }
    return "unknown";
}

/**
 * @brief Type-erased assembly context passed to internal term applicators.
 *
 * The context is intentionally lightweight and non-owning. Concrete builders
 * set the references that are meaningful for the current assembly path and
 * callbacks recover them with typed accessors such as `matrixAs<T>()` or
 * `bilinearFormAs<T>()`.
 */
struct AssemblyContext
{
    /** @brief Current assembly mode. */
    AssemblyMode mode = AssemblyMode::Unknown;
    /** @brief Current term storage role, if known. */
    std::optional<term_role> role;
    /** @brief Optional product-space row index. */
    std::optional<int> row;
    /** @brief Optional product-space column index. */
    std::optional<int> column;

    /** @brief Non-owning pointer to the current problem descriptor. */
    void const* problem = nullptr;
    /** @brief Non-owning pointer to the current scalar or product space. */
    void const* space = nullptr;
    /** @brief Non-owning pointer to `time_step_context`, when available. */
    void const* timeStep = nullptr;
    /** @brief Non-owning pointer to the current nonlinear/Picard state. */
    void const* state = nullptr;
    /** @brief Non-owning pointer to the backend, when available. */
    void const* backend = nullptr;
    /** @brief Non-owning pointer to the concrete time scheme, when available. */
    void const* timeScheme = nullptr;

    /** @brief Non-owning pointer to a matrix pointer object. */
    void const* matrix = nullptr;
    /** @brief Non-owning pointer to a vector pointer object. */
    void const* vector = nullptr;
    /** @brief Non-owning pointer to a live scalar bilinear form. */
    void* bilinearForm = nullptr;
    /** @brief Non-owning pointer to a live scalar linear form. */
    void* linearForm = nullptr;
    /** @brief Non-owning pointer to a live block bilinear form. */
    void* blockBilinearForm = nullptr;
    /** @brief Non-owning pointer to a live block linear form. */
    void* blockLinearForm = nullptr;

    /**
     * @brief Return true when a term role is available.
     */
    bool hasRole() const noexcept { return static_cast<bool>( role ); }

    /**
     * @brief Return true when a block row index is available.
     */
    bool hasRow() const noexcept { return static_cast<bool>( row ); }

    /**
     * @brief Return true when a block column index is available.
     */
    bool hasColumn() const noexcept { return static_cast<bool>( column ); }

    /**
     * @brief Return true when time-step metadata is available.
     */
    bool hasTimeStep() const noexcept { return timeStep != nullptr; }

    /**
     * @brief Return true when a nonlinear or Picard state is available.
     */
    bool hasState() const noexcept { return state != nullptr; }

    /** @brief Cast and dereference `problem` as `T const&`. */
    template <typename T>
    T const& problemAs() const { return *static_cast<T const*>( problem ); }

    /** @brief Cast and dereference `space` as `T const&`. */
    template <typename T>
    T const& spaceAs() const { return *static_cast<T const*>( space ); }

    /** @brief Cast and dereference `timeStep` as `T const&`. */
    template <typename T>
    T const& timeStepAs() const { return *static_cast<T const*>( timeStep ); }

    /** @brief Cast and dereference `state` as `T const&`. */
    template <typename T>
    T const& stateAs() const { return *static_cast<T const*>( state ); }

    /** @brief Cast and dereference `backend` as `T const&`. */
    template <typename T>
    T const& backendAs() const { return *static_cast<T const*>( backend ); }

    /** @brief Cast and dereference `timeScheme` as `T const&`. */
    template <typename T>
    T const& timeSchemeAs() const { return *static_cast<T const*>( timeScheme ); }

    /** @brief Cast and dereference `matrix` as `T const&`. */
    template <typename T>
    T const& matrixAs() const { return *static_cast<T const*>( matrix ); }

    /** @brief Cast and dereference `vector` as `T const&`. */
    template <typename T>
    T const& vectorAs() const { return *static_cast<T const*>( vector ); }

    /** @brief Cast and dereference `bilinearForm` as `T&`. */
    template <typename T>
    T& bilinearFormAs() const { return *static_cast<T*>( bilinearForm ); }

    /** @brief Cast and dereference `linearForm` as `T&`. */
    template <typename T>
    T& linearFormAs() const { return *static_cast<T*>( linearForm ); }

    /** @brief Cast and dereference `blockBilinearForm` as `T&`. */
    template <typename T>
    T& blockBilinearFormAs() const { return *static_cast<T*>( blockBilinearForm ); }

    /** @brief Cast and dereference `blockLinearForm` as `T&`. */
    template <typename T>
    T& blockLinearFormAs() const { return *static_cast<T*>( blockLinearForm ); }
};

namespace detail
{
/**
 * @brief Fill missing active-system handles with owned algebraic storage.
 */
template <typename ActiveSystemType, typename MatrixPtr, typename VectorPtr>
ActiveSystemType
completeActiveSystem( ActiveSystemType activeSystem,
                      MatrixPtr const& fallbackMatrix,
                      VectorPtr const& fallbackVector )
{
    if ( !activeSystem.first )
        activeSystem.first = fallbackMatrix;
    if ( !activeSystem.second )
        activeSystem.second = fallbackVector;
    return activeSystem;
}
} // namespace detail

/**
 * @brief Mathematical category associated with a stored term.
 *
 * The kind is intentionally conservative: it describes how FeelSystem should
 * reason about a term, not a fully introspected symbolic expression.
 */
enum class TermKind
{
    Lhs,              /**< Generic lhs matrix term. */
    Rhs,              /**< Generic rhs vector term. */
    Mass,             /**< Mass/time-derivative term. */
    Damping,          /**< Damping term. */
    Stiffness,        /**< Stiffness/diffusion term. */
    Reaction,         /**< Reaction term. */
    Source,           /**< Source term. */
    Residual,         /**< Nonlinear residual term. */
    Jacobian,         /**< Nonlinear Jacobian term. */
    PicardIteration,  /**< Picard fixed-point iteration term. */
    OperatorA,        /**< Eigen operator A term. */
    OperatorB,        /**< Eigen operator B term. */
    Block,            /**< Product-space matrix block term. */
    Row,              /**< Product-space vector row term. */
    Constraint,       /**< Algebraic constraint term. */
    BoundaryCondition, /**< Boundary-condition term. */
    TimeScheme,       /**< Time-scheme descriptor term. */
    HdgVolume,        /**< HDG volume term. */
    HdgSkeleton,      /**< HDG skeleton term. */
    HdgStabilization, /**< HDG stabilization term. */
    HdgFlux,          /**< HDG flux term. */
    HdgConstraint,    /**< HDG constraint term. */
    Unknown           /**< Unknown or user-defined term kind. */
};

/**
 * @brief Geometric or algebraic support associated with a term.
 */
enum class TermSupport
{
    Volume,        /**< Volume mesh support. */
    BoundaryFacet, /**< Boundary-facet support. */
    InteriorFacet, /**< Interior-facet support. */
    Skeleton,      /**< Skeleton/trace support. */
    Algebraic,     /**< Algebraic support without direct mesh range. */
    Time,          /**< Time-scheme support. */
    Unknown        /**< Unknown or user-defined support. */
};

/**
 * @brief Linearity classification used by assembly planning.
 */
enum class TermLinearity
{
    Linear,    /**< Linear vector contribution. */
    Bilinear,  /**< Bilinear matrix contribution. */
    Affine,    /**< Affine contribution. */
    Nonlinear, /**< State-dependent nonlinear contribution. */
    Unknown    /**< Unknown or user-defined linearity. */
};

/**
 * @brief Algebraic role played by a stored term.
 */
enum class TermRole
{
    Matrix,         /**< Matrix contribution. */
    Vector,         /**< Vector contribution. */
    Residual,       /**< Residual contribution. */
    Jacobian,       /**< Jacobian contribution. */
    Constraint,     /**< Constraint contribution. */
    TimeDerivative, /**< Time-derivative contribution. */
    LocalCell,      /**< Local cell block contribution. */
    Trace,          /**< Trace/skeleton contribution. */
    Coupling,       /**< Coupling contribution. */
    Stabilization,  /**< Stabilization contribution. */
    Flux,           /**< Flux contribution. */
    Unknown         /**< Unknown or user-defined algebraic role. */
};

/**
 * @brief Side-car metadata carried by each erased FeelSystem term.
 *
 * This metadata is deliberately descriptive. It is not a replacement for the
 * Feel++ DSEL expression and it does not try to infer complete symbolic
 * structure. Unknown or empty values are valid when a builder cannot classify a
 * term safely.
 */
struct TermMetadata
{
    /**
     * @brief Return true when a block row index is available.
     */
    bool hasRow() const noexcept { return static_cast<bool>( row ); }

    /**
     * @brief Return true when a block column index is available.
     */
    bool hasColumn() const noexcept { return static_cast<bool>( column ); }

    /** @brief Stable short name for diagnostics and planning. */
    std::string name;
    /** @brief Optional equation or model component label. */
    std::string equation;
    /** @brief Optional human-readable description. */
    std::string description;

    /** @brief Mathematical category of the term. */
    TermKind kind = TermKind::Unknown;
    /** @brief Mesh or algebraic support touched by the term. */
    TermSupport support = TermSupport::Unknown;
    /** @brief Linearity classification. */
    TermLinearity linearity = TermLinearity::Unknown;
    /** @brief Algebraic role used by diagnostics and planning. */
    TermRole role = TermRole::Unknown;

    /** @brief Optional product-space row index. */
    std::optional<int> row;
    /** @brief Optional product-space column index. */
    std::optional<int> column;

    /** @brief True when the term is affine in the unknown. */
    bool isAffine = false;
    /** @brief True when the term depends on time. */
    bool isTimeDependent = false;
    /** @brief True when the term depends on symbolic parameters. */
    bool isParameterDependent = false;
    /** @brief True when the term depends on the nonlinear/Picard state. */
    bool isStateDependent = false;

    /** @brief True when the sparse CPU assembly path is a candidate. */
    bool isSparseAssemblyCandidate = true;
    /** @brief True when device execution may be legal. */
    bool isDeviceCandidate = false;
    /** @brief True when matrix-free evaluation may be legal. */
    bool isMatrixFreeCandidate = false;
    /** @brief True when the term is shaped like an HDG local operation. */
    bool isHdgCandidate = false;
    /** @brief True when the term may participate in a trace Schur complement. */
    bool isTraceSchurCandidate = false;
    /** @brief True when automatic differentiation may generate related terms. */
    bool isAdCandidate = false;
    /** @brief True when reduced-order-model projection may be legal. */
    bool isRomCandidate = false;
    /** @brief True when ensemble reuse may be legal. */
    bool isEnsembleCandidate = false;

    /** @brief Field names read or written by the term. */
    std::vector<std::string> fields;
    /** @brief Symbolic parameter names read by the term. */
    std::vector<std::string> parameters;
    /** @brief Marker names touched by the term. */
    std::vector<std::string> markers;
};

/**
 * @brief Human-readable name for a term kind.
 */
constexpr std::string_view
termKindName( TermKind kind ) noexcept
{
    switch ( kind )
    {
    case TermKind::Lhs:
        return "lhs";
    case TermKind::Rhs:
        return "rhs";
    case TermKind::Mass:
        return "mass";
    case TermKind::Damping:
        return "damping";
    case TermKind::Stiffness:
        return "stiffness";
    case TermKind::Reaction:
        return "reaction";
    case TermKind::Source:
        return "source";
    case TermKind::Residual:
        return "residual";
    case TermKind::Jacobian:
        return "jacobian";
    case TermKind::PicardIteration:
        return "picard_iteration";
    case TermKind::OperatorA:
        return "operator_a";
    case TermKind::OperatorB:
        return "operator_b";
    case TermKind::Block:
        return "block";
    case TermKind::Row:
        return "row";
    case TermKind::Constraint:
        return "constraint";
    case TermKind::BoundaryCondition:
        return "boundary_condition";
    case TermKind::TimeScheme:
        return "time_scheme";
    case TermKind::HdgVolume:
        return "hdg_volume";
    case TermKind::HdgSkeleton:
        return "hdg_skeleton";
    case TermKind::HdgStabilization:
        return "hdg_stabilization";
    case TermKind::HdgFlux:
        return "hdg_flux";
    case TermKind::HdgConstraint:
        return "hdg_constraint";
    case TermKind::Unknown:
        return "unknown";
    }
    return "unknown";
}

/**
 * @brief Human-readable name for a term support.
 */
constexpr std::string_view
termSupportName( TermSupport support ) noexcept
{
    switch ( support )
    {
    case TermSupport::Volume:
        return "volume";
    case TermSupport::BoundaryFacet:
        return "boundary_facet";
    case TermSupport::InteriorFacet:
        return "interior_facet";
    case TermSupport::Skeleton:
        return "skeleton";
    case TermSupport::Algebraic:
        return "algebraic";
    case TermSupport::Time:
        return "time";
    case TermSupport::Unknown:
        return "unknown";
    }
    return "unknown";
}

/**
 * @brief Human-readable name for term linearity.
 */
constexpr std::string_view
termLinearityName( TermLinearity linearity ) noexcept
{
    switch ( linearity )
    {
    case TermLinearity::Linear:
        return "linear";
    case TermLinearity::Bilinear:
        return "bilinear";
    case TermLinearity::Affine:
        return "affine";
    case TermLinearity::Nonlinear:
        return "nonlinear";
    case TermLinearity::Unknown:
        return "unknown";
    }
    return "unknown";
}

/**
 * @brief Human-readable name for an algebraic term role.
 */
constexpr std::string_view
termRoleName( TermRole role ) noexcept
{
    switch ( role )
    {
    case TermRole::Matrix:
        return "matrix";
    case TermRole::Vector:
        return "vector";
    case TermRole::Residual:
        return "residual";
    case TermRole::Jacobian:
        return "jacobian";
    case TermRole::Constraint:
        return "constraint";
    case TermRole::TimeDerivative:
        return "time_derivative";
    case TermRole::LocalCell:
        return "local_cell";
    case TermRole::Trace:
        return "trace";
    case TermRole::Coupling:
        return "coupling";
    case TermRole::Stabilization:
        return "stabilization";
    case TermRole::Flux:
        return "flux";
    case TermRole::Unknown:
        return "unknown";
    }
    return "unknown";
}

/**
 * @brief Structured, side-effect-free summary of one stored term.
 *
 * The summary copies the metadata fields that are useful for diagnostics and
 * planning. It never owns or invokes the erased assembly callback.
 */
struct TermSummary
{
    /** @brief Insertion-order term index. */
    std::size_t index = 0;
    /** @brief Storage role used by the problem descriptor. */
    term_role storageRole = term_role::lhs;
    /** @brief Legacy descriptor label derived from the storage role. */
    std::string label;
    /** @brief Metadata name, when provided. */
    std::string name;
    /** @brief Mathematical category. */
    TermKind kind = TermKind::Unknown;
    /** @brief Mesh or algebraic support. */
    TermSupport support = TermSupport::Unknown;
    /** @brief Linearity classification. */
    TermLinearity linearity = TermLinearity::Unknown;
    /** @brief Algebraic role. */
    TermRole role = TermRole::Unknown;
    /** @brief Optional product-space row index. */
    std::optional<int> row;
    /** @brief Optional product-space column index. */
    std::optional<int> column;
    /** @brief True when the term is affine in the unknown. */
    bool isAffine = false;
    /** @brief True when the term depends on time. */
    bool isTimeDependent = false;
    /** @brief True when the term depends on symbolic parameters. */
    bool isParameterDependent = false;
    /** @brief True when the term depends on the nonlinear/Picard state. */
    bool isStateDependent = false;
    /** @brief True when sparse CPU assembly is a candidate. */
    bool isSparseAssemblyCandidate = true;
    /** @brief True when device execution may be legal. */
    bool isDeviceCandidate = false;
    /** @brief True when matrix-free evaluation may be legal. */
    bool isMatrixFreeCandidate = false;
    /** @brief True when the term is shaped like an HDG local operation. */
    bool isHdgCandidate = false;
    /** @brief True when the term may participate in a trace Schur complement. */
    bool isTraceSchurCandidate = false;
    /** @brief True when automatic differentiation may generate related terms. */
    bool isAdCandidate = false;
    /** @brief True when reduced-order-model projection may be legal. */
    bool isRomCandidate = false;
    /** @brief True when ensemble reuse may be legal. */
    bool isEnsembleCandidate = false;
    /** @brief Field names read or written by the term. */
    std::vector<std::string> fields;
    /** @brief Symbolic parameter names read by the term. */
    std::vector<std::string> parameters;
    /** @brief Marker names touched by the term. */
    std::vector<std::string> markers;

    /**
     * @brief Return true when a block row index is available.
     */
    bool hasRow() const noexcept { return static_cast<bool>( row ); }

    /**
     * @brief Return true when a block column index is available.
     */
    bool hasColumn() const noexcept { return static_cast<bool>( column ); }
};

/**
 * @brief One conservative diagnostic issue found in stored term metadata.
 */
struct TermConsistencyIssue
{
    /** @brief Insertion-order term index that raised the diagnostic. */
    std::size_t index = 0;
    /** @brief Diagnostic severity such as `"warning"`. */
    std::string severity = "warning";
    /** @brief Human-readable diagnostic message. */
    std::string message;
};

/**
 * @brief Strategy candidate flags grouped in an assembly plan batch.
 *
 * These flags describe strategy eligibility for diagnostics and task-graph
 * construction, based on a batch's term metadata.
 */
struct AssemblyStrategyFlags
{
    /** @brief True when a batch can use the current sparse CPU assembly path. */
    bool sparseCpuCandidate = true;
    /** @brief True when a batch can be evaluated as an operator action. */
    bool matrixFreeCandidate = false;
    /** @brief True when a batch may be suitable for device execution. */
    bool deviceCandidate = false;
    /** @brief True when a batch is shaped like an HDG local dense operation. */
    bool hdgLocalDenseCandidate = false;
    /** @brief True when a batch can contribute to a trace Schur complement. */
    bool traceSchurCandidate = false;
    /** @brief True when a batch may be projected into a reduced basis. */
    bool romCandidate = false;
    /** @brief True when a batch can be shared across ensemble realizations. */
    bool ensembleCandidate = false;
    /** @brief True when a Jacobian batch can be generated from AD metadata. */
    bool adGeneratedJacobianCandidate = false;

    friend bool operator==( AssemblyStrategyFlags const&, AssemblyStrategyFlags const& ) = default;
};

/**
 * @brief One metadata batch in an assembly plan.
 */
struct AssemblyBatch
{
    /** @brief Geometric or algebraic support shared by this batch. */
    TermSupport support = TermSupport::Unknown;
    /** @brief Algebraic role shared by this batch. */
    TermRole role = TermRole::Unknown;
    /** @brief Optional block row targeted by this batch. */
    std::optional<int> row;
    /** @brief Optional block column targeted by this batch. */
    std::optional<int> column;
    /** @brief Strategy candidates shared by this batch. */
    AssemblyStrategyFlags strategy;
    /** @brief Insertion-order term ids belonging to this metadata batch. */
    std::vector<std::size_t> termIds;
    /** @brief Field names read by the terms in this batch. */
    std::vector<std::string> fields;
    /** @brief Parameter names read by the terms in this batch. */
    std::vector<std::string> parameters;
    /** @brief Marker names touched by the terms in this batch. */
    std::vector<std::string> markers;

    /**
     * @brief Return true when a block row index is available.
     */
    bool hasRow() const noexcept { return static_cast<bool>( row ); }

    /**
     * @brief Return true when a block column index is available.
     */
    bool hasColumn() const noexcept { return static_cast<bool>( column ); }
};

/**
 * @brief Debug-oriented assembly plan built from term metadata.
 */
struct AssemblyPlan
{
    /** @brief Metadata batches sorted by support, role, block indices, and strategy flags. */
    std::vector<AssemblyBatch> batches;
    /** @brief Conservative metadata diagnostics collected while building the plan. */
    std::vector<TermConsistencyIssue> issues;

    /**
     * @brief Return true when the plan has no batches.
     */
    bool empty() const noexcept { return batches.empty(); }

    /**
     * @brief Number of term ids referenced by all batches.
     */
    std::size_t nTerms() const noexcept
    {
        std::size_t count = 0;
        for ( auto const& batch : batches )
            count += batch.termIds.size();
        return count;
    }
};

/**
 * @brief Options controlling how a debug assembly plan is derived.
 */
struct AssemblyPlanOptions
{
    /** @brief Treat generic block metadata as an HDG assembly context. */
    bool hdgContext = false;
};

struct term_descriptor;

namespace detail
{
/**
 * @brief Join a list of strings for compact diagnostics.
 */
inline std::string
joinDiagnosticStrings( std::vector<std::string> const& values )
{
    std::ostringstream os;
    for ( std::size_t i = 0; i != values.size(); ++i )
    {
        if ( i != 0 )
            os << ",";
        os << values[i];
    }
    return os.str();
}

/**
 * @brief Append strings that are not already present.
 */
inline void
appendUniqueDiagnosticStrings( std::vector<std::string>& target,
                               std::vector<std::string> const& source )
{
    for ( auto const& value : source )
    {
        if ( std::find( target.begin(), target.end(), value ) == target.end() )
            target.push_back( value );
    }
}

/**
 * @brief Join term ids for compact diagnostics.
 */
inline std::string
joinDiagnosticIds( std::vector<std::size_t> const& values )
{
    std::ostringstream os;
    for ( std::size_t i = 0; i != values.size(); ++i )
    {
        if ( i != 0 )
            os << ",";
        os << values[i];
    }
    return os.str();
}

/**
 * @brief Return true when an algebraic role is matrix-like enough for lhs storage.
 */
constexpr bool
isMatrixLikeRole( TermRole role ) noexcept
{
    return role == TermRole::Matrix ||
           role == TermRole::LocalCell ||
           role == TermRole::Trace ||
           role == TermRole::Coupling ||
           role == TermRole::Stabilization ||
           role == TermRole::Flux;
}

/**
 * @brief Return true when metadata role is compatible with a storage role.
 */
constexpr bool
isCompatibleStorageRole( term_role storageRole, TermRole role ) noexcept
{
    switch ( storageRole )
    {
    case term_role::lhs:
    case term_role::damping:
    case term_role::stiffness:
    case term_role::operator_a:
    case term_role::operator_b:
        return isMatrixLikeRole( role );
    case term_role::rhs:
        return role == TermRole::Vector;
    case term_role::mass:
    case term_role::time_scheme:
        return role == TermRole::TimeDerivative || isMatrixLikeRole( role );
    case term_role::residual:
        return role == TermRole::Residual;
    case term_role::jacobian:
        return role == TermRole::Jacobian;
    case term_role::iteration:
        return isMatrixLikeRole( role );
    case term_role::constraint:
        return role == TermRole::Constraint;
    }
    return true;
}

/**
 * @brief Return true when a kind belongs to the HDG metadata family.
 */
constexpr bool
isHdgTermKind( TermKind kind ) noexcept
{
    return kind == TermKind::HdgVolume ||
           kind == TermKind::HdgSkeleton ||
           kind == TermKind::HdgStabilization ||
           kind == TermKind::HdgFlux ||
           kind == TermKind::HdgConstraint;
}

/**
 * @brief Return true when a role is naturally local/trace-like.
 */
constexpr bool
isLocalOrTraceRole( TermRole role ) noexcept
{
    return role == TermRole::LocalCell ||
           role == TermRole::Trace ||
           role == TermRole::Coupling ||
           role == TermRole::Stabilization ||
           role == TermRole::Flux;
}
} // namespace detail

/**
 * @brief Build a structured summary for one term descriptor.
 */
inline TermSummary
termSummary( std::size_t index, term_descriptor const& descriptor );

/**
 * @brief Check stored term metadata for conservative consistency issues.
 */
template <typename TermContainer>
std::vector<TermConsistencyIssue>
checkTermConsistency( TermContainer const& terms );

/**
 * @brief Build structured summaries for all term descriptors in insertion order.
 */
template <typename TermContainer>
std::vector<TermSummary>
termSummaries( TermContainer const& terms )
{
    std::vector<TermSummary> summaries;
    summaries.reserve( terms.size() );
    for ( std::size_t i = 0; i != terms.size(); ++i )
        summaries.push_back( termSummary( i, terms[i] ) );
    return summaries;
}

/**
 * @brief Derive strategy candidate flags from one term summary.
 */
inline AssemblyStrategyFlags
assemblyStrategyFlags( TermSummary const& summary,
                       AssemblyPlanOptions const& options = {} )
{
    AssemblyStrategyFlags flags;
    flags.sparseCpuCandidate = summary.isSparseAssemblyCandidate;
    flags.matrixFreeCandidate = summary.isMatrixFreeCandidate;
    flags.deviceCandidate = summary.isDeviceCandidate;
    flags.hdgLocalDenseCandidate =
        options.hdgContext ||
        summary.isHdgCandidate ||
        detail::isHdgTermKind( summary.kind ) ||
        summary.support == TermSupport::Skeleton ||
        detail::isLocalOrTraceRole( summary.role );
    flags.traceSchurCandidate =
        options.hdgContext ||
        summary.isTraceSchurCandidate ||
        summary.kind == TermKind::HdgSkeleton ||
        summary.kind == TermKind::HdgConstraint ||
        summary.role == TermRole::Trace;
    flags.romCandidate = summary.isRomCandidate;
    flags.ensembleCandidate = summary.isEnsembleCandidate;
    flags.adGeneratedJacobianCandidate = summary.isAdCandidate;
    return flags;
}

/**
 * @brief Return true when a summary belongs to an existing assembly batch.
 */
inline bool
sameAssemblyBatchKey( AssemblyBatch const& batch,
                      TermSummary const& summary,
                      AssemblyStrategyFlags const& flags )
{
    return batch.support == summary.support &&
           batch.role == summary.role &&
           batch.row == summary.row &&
           batch.column == summary.column &&
           batch.strategy == flags;
}

namespace detail
{
/**
 * @brief Return a sortable key for strategy flags.
 */
inline auto
assemblyStrategySortKey( AssemblyStrategyFlags const& strategy )
{
    return std::tuple{
        strategy.sparseCpuCandidate,
        strategy.matrixFreeCandidate,
        strategy.deviceCandidate,
        strategy.hdgLocalDenseCandidate,
        strategy.traceSchurCandidate,
        strategy.romCandidate,
        strategy.ensembleCandidate,
        strategy.adGeneratedJacobianCandidate };
}

/**
 * @brief Deterministic ordering for assembly-plan batches.
 */
inline bool
assemblyBatchLess( AssemblyBatch const& lhs, AssemblyBatch const& rhs )
{
    return std::tuple{
               lhs.support,
               lhs.role,
               lhs.row,
               lhs.column,
               assemblyStrategySortKey( lhs.strategy ) } <
           std::tuple{
               rhs.support,
               rhs.role,
               rhs.row,
               rhs.column,
               assemblyStrategySortKey( rhs.strategy ) };
}
} // namespace detail

/**
 * @brief Build a debug assembly plan from term descriptors.
 *
 * The plan is built only from metadata. It never creates a backend, never
 * assembles a form, and never evaluates an expression.
 */
template <typename TermContainer>
AssemblyPlan
buildAssemblyPlan( TermContainer const& terms,
                   AssemblyPlanOptions const& options = {} )
{
    AssemblyPlan plan;
    plan.issues = checkTermConsistency( terms );

    for ( auto const& summary : termSummaries( terms ) )
    {
        auto flags = assemblyStrategyFlags( summary, options );
        bool inserted = false;
        for ( auto& batch : plan.batches )
        {
            if ( sameAssemblyBatchKey( batch, summary, flags ) )
            {
                batch.termIds.push_back( summary.index );
                detail::appendUniqueDiagnosticStrings( batch.fields, summary.fields );
                detail::appendUniqueDiagnosticStrings( batch.parameters, summary.parameters );
                detail::appendUniqueDiagnosticStrings( batch.markers, summary.markers );
                inserted = true;
                break;
            }
        }

        if ( !inserted )
        {
            AssemblyBatch batch;
            batch.support = summary.support;
            batch.role = summary.role;
            batch.row = summary.row;
            batch.column = summary.column;
            batch.strategy = flags;
            batch.termIds.push_back( summary.index );
            batch.fields = summary.fields;
            batch.parameters = summary.parameters;
            batch.markers = summary.markers;
            plan.batches.push_back( std::move( batch ) );
        }
    }

    std::stable_sort( plan.batches.begin(), plan.batches.end(), detail::assemblyBatchLess );
    return plan;
}

/**
 * @brief Check stored term metadata for conservative consistency issues.
 *
 * The checker is intentionally permissive: unknown metadata is allowed and
 * custom callbacks are not rejected. The returned issues are diagnostics, not
 * assembly or solve policy.
 */
template <typename TermContainer>
std::vector<TermConsistencyIssue>
checkTermConsistency( TermContainer const& terms )
{
    std::vector<TermConsistencyIssue> issues;
    for ( std::size_t i = 0; i != terms.size(); ++i )
    {
        auto const& descriptor = terms[i];
        auto const& metadata = descriptor.metadata();

        if ( metadata.role != TermRole::Unknown &&
             !detail::isCompatibleStorageRole( descriptor.role, metadata.role ) )
        {
            issues.push_back( TermConsistencyIssue{
                i,
                "warning",
                "metadata role is unusual for the term storage role" } );
        }

        if ( metadata.kind == TermKind::Block )
        {
            if ( !metadata.hasRow() || !metadata.hasColumn() )
                issues.push_back( TermConsistencyIssue{
                    i,
                    "warning",
                    "block term should carry row and column metadata" } );
        }

        if ( metadata.kind == TermKind::Row && !metadata.hasRow() )
        {
            issues.push_back( TermConsistencyIssue{
                i,
                "warning",
                "row term should carry row metadata" } );
        }

        if ( metadata.hasColumn() && !metadata.hasRow() )
        {
            issues.push_back( TermConsistencyIssue{
                i,
                "warning",
                "column metadata without row metadata is ambiguous" } );
        }
    }
    return issues;
}

/**
 * @brief Write one compact term summary line.
 */
inline std::ostream&
describeTermSummary( std::ostream& os, TermSummary const& summary )
{
    os << "  [" << summary.index << "]"
       << " name=" << ( summary.name.empty() ? summary.label : summary.name )
       << " term=" << roleName( summary.storageRole )
       << " kind=" << termKindName( summary.kind )
       << " support=" << termSupportName( summary.support )
       << " role=" << termRoleName( summary.role )
       << " linearity=" << termLinearityName( summary.linearity );

    if ( summary.hasRow() )
        os << " row=" << *summary.row;
    if ( summary.hasColumn() )
        os << " col=" << *summary.column;
    if ( summary.isAffine )
        os << " affine=true";
    if ( summary.isTimeDependent )
        os << " timeDependent=true";
    if ( summary.isParameterDependent )
        os << " parameterDependent=true";
    if ( summary.isStateDependent )
        os << " stateDependent=true";
    if ( summary.isDeviceCandidate )
        os << " deviceCandidate=true";
    if ( summary.isMatrixFreeCandidate )
        os << " matrixFreeCandidate=true";
    if ( summary.isHdgCandidate )
        os << " hdgCandidate=true";
    if ( summary.isAdCandidate )
        os << " adCandidate=true";
    if ( summary.isRomCandidate )
        os << " romCandidate=true";
    if ( summary.isEnsembleCandidate )
        os << " ensembleCandidate=true";
    if ( !summary.fields.empty() )
        os << " fields=" << detail::joinDiagnosticStrings( summary.fields );
    if ( !summary.parameters.empty() )
        os << " parameters=" << detail::joinDiagnosticStrings( summary.parameters );
    if ( !summary.markers.empty() )
        os << " markers=" << detail::joinDiagnosticStrings( summary.markers );
    return os;
}

/**
 * @brief Write a side-effect-free diagnostic description for a term list.
 */
template <typename TermContainer>
std::ostream&
describeTerms( std::ostream& os, std::string_view problemName, TermContainer const& terms )
{
    os << "Problem: " << ( problemName.empty() ? std::string_view{ "<unnamed>" } : problemName ) << '\n';
    os << "Terms:\n";
    for ( auto const& summary : termSummaries( terms ) )
        describeTermSummary( os, summary ) << '\n';

    auto issues = checkTermConsistency( terms );
    if ( !issues.empty() )
    {
        os << "Issues:\n";
        for ( auto const& issue : issues )
            os << "  [" << issue.index << "] " << issue.severity << ": " << issue.message << '\n';
    }
    return os;
}

/**
 * @brief Return a side-effect-free diagnostic description as a string.
 */
template <typename TermContainer>
std::string
describeTerms( std::string_view problemName, TermContainer const& terms )
{
    std::ostringstream os;
    describeTerms( os, problemName, terms );
    return os.str();
}

/**
 * @brief Write one compact assembly-plan batch line.
 */
inline std::ostream&
describeAssemblyBatch( std::ostream& os, std::size_t index, AssemblyBatch const& batch )
{
    os << "  [" << index << "]"
       << " support=" << termSupportName( batch.support )
       << " role=" << termRoleName( batch.role );
    if ( batch.hasRow() )
        os << " row=" << *batch.row;
    if ( batch.hasColumn() )
        os << " col=" << *batch.column;
    os << " terms=" << detail::joinDiagnosticIds( batch.termIds );
    if ( batch.strategy.sparseCpuCandidate )
        os << " sparseCpu=true";
    if ( batch.strategy.matrixFreeCandidate )
        os << " matrixFree=true";
    if ( batch.strategy.deviceCandidate )
        os << " device=true";
    if ( batch.strategy.hdgLocalDenseCandidate )
        os << " hdgLocalDense=true";
    if ( batch.strategy.traceSchurCandidate )
        os << " traceSchur=true";
    if ( batch.strategy.romCandidate )
        os << " rom=true";
    if ( batch.strategy.ensembleCandidate )
        os << " ensemble=true";
    if ( batch.strategy.adGeneratedJacobianCandidate )
        os << " adGeneratedJacobian=true";
    if ( !batch.fields.empty() )
        os << " fields=" << detail::joinDiagnosticStrings( batch.fields );
    if ( !batch.parameters.empty() )
        os << " parameters=" << detail::joinDiagnosticStrings( batch.parameters );
    if ( !batch.markers.empty() )
        os << " markers=" << detail::joinDiagnosticStrings( batch.markers );
    return os;
}

/**
 * @brief Write a side-effect-free assembly plan description.
 */
inline std::ostream&
describeAssemblyPlan( std::ostream& os, std::string_view problemName, AssemblyPlan const& plan )
{
    os << "AssemblyPlan: " << ( problemName.empty() ? std::string_view{ "<unnamed>" } : problemName ) << '\n';
    os << "Batches:\n";
    for ( std::size_t i = 0; i != plan.batches.size(); ++i )
        describeAssemblyBatch( os, i, plan.batches[i] ) << '\n';

    if ( !plan.issues.empty() )
    {
        os << "Issues:\n";
        for ( auto const& issue : plan.issues )
            os << "  [" << issue.index << "] " << issue.severity << ": " << issue.message << '\n';
    }
    return os;
}

/**
 * @brief Return a side-effect-free assembly plan description as a string.
 */
inline std::string
describeAssemblyPlan( std::string_view problemName, AssemblyPlan const& plan )
{
    std::ostringstream os;
    describeAssemblyPlan( os, problemName, plan );
    return os.str();
}

/**
 * @brief Return conservative default metadata for an existing term role.
 */
inline TermMetadata
defaultTermMetadata( term_role role )
{
    TermMetadata metadata;
    metadata.name = std::string( roleName( role ) );

    switch ( role )
    {
    case term_role::lhs:
        metadata.kind = TermKind::Lhs;
        metadata.role = TermRole::Matrix;
        metadata.linearity = TermLinearity::Bilinear;
        metadata.isAffine = true;
        break;
    case term_role::rhs:
        metadata.kind = TermKind::Rhs;
        metadata.role = TermRole::Vector;
        metadata.linearity = TermLinearity::Linear;
        metadata.isAffine = true;
        break;
    case term_role::mass:
        metadata.kind = TermKind::Mass;
        metadata.role = TermRole::TimeDerivative;
        metadata.linearity = TermLinearity::Bilinear;
        metadata.isAffine = true;
        metadata.isTimeDependent = true;
        break;
    case term_role::damping:
        metadata.kind = TermKind::Damping;
        metadata.role = TermRole::Matrix;
        metadata.linearity = TermLinearity::Bilinear;
        metadata.isAffine = true;
        metadata.isTimeDependent = true;
        break;
    case term_role::stiffness:
        metadata.kind = TermKind::Stiffness;
        metadata.role = TermRole::Matrix;
        metadata.linearity = TermLinearity::Bilinear;
        metadata.isAffine = true;
        break;
    case term_role::residual:
        metadata.kind = TermKind::Residual;
        metadata.role = TermRole::Residual;
        metadata.linearity = TermLinearity::Nonlinear;
        metadata.isStateDependent = true;
        break;
    case term_role::jacobian:
        metadata.kind = TermKind::Jacobian;
        metadata.role = TermRole::Jacobian;
        metadata.linearity = TermLinearity::Nonlinear;
        metadata.isStateDependent = true;
        break;
    case term_role::iteration:
        metadata.kind = TermKind::PicardIteration;
        metadata.role = TermRole::Matrix;
        metadata.linearity = TermLinearity::Nonlinear;
        metadata.isStateDependent = true;
        break;
    case term_role::operator_a:
        metadata.kind = TermKind::OperatorA;
        metadata.role = TermRole::Matrix;
        metadata.linearity = TermLinearity::Bilinear;
        metadata.isAffine = true;
        break;
    case term_role::operator_b:
        metadata.kind = TermKind::OperatorB;
        metadata.role = TermRole::Matrix;
        metadata.linearity = TermLinearity::Bilinear;
        metadata.isAffine = true;
        break;
    case term_role::constraint:
        metadata.kind = TermKind::Constraint;
        metadata.support = TermSupport::Algebraic;
        metadata.role = TermRole::Constraint;
        metadata.linearity = TermLinearity::Unknown;
        break;
    case term_role::time_scheme:
        metadata.kind = TermKind::TimeScheme;
        metadata.support = TermSupport::Time;
        metadata.role = TermRole::TimeDerivative;
        metadata.linearity = TermLinearity::Unknown;
        metadata.isTimeDependent = true;
        break;
    }

    return metadata;
}

/**
 * @brief Descriptor associated with an erased term callback.
 *
 * `role` and `label` remain public for compatibility with existing tests and
 * counting code. `metadata()` exposes the richer side-car used by newer
 * planning APIs.
 */
struct term_descriptor
{
    /** @brief Construct a default lhs descriptor. */
    term_descriptor() = default;

    /**
     * @brief Construct a descriptor with default metadata for `roleValue`.
     */
    term_descriptor( term_role roleValue, std::string_view labelValue )
        : term_descriptor( roleValue, labelValue, defaultTermMetadata( roleValue ) )
    {
    }

    /**
     * @brief Construct a descriptor with explicit metadata.
     */
    term_descriptor( term_role roleValue, std::string_view labelValue, TermMetadata metadataValue )
        : role( roleValue ),
          label( labelValue ),
          M_metadata( std::move( metadataValue ) )
    {
        if ( M_metadata.name.empty() )
            M_metadata.name = std::string( label );
    }

    /**
     * @brief Metadata carried by this descriptor.
     */
    TermMetadata const& metadata() const noexcept { return M_metadata; }

    /**
     * @brief Mutable metadata for internal builders.
     */
    TermMetadata& metadata() noexcept { return M_metadata; }

    /** @brief Storage role associated with the erased term callback. */
    term_role role = term_role::lhs;
    /** @brief Legacy short label associated with the storage role. */
    std::string_view label = {};

  private:
    /** @brief Rich metadata used by diagnostics and planning APIs. */
    TermMetadata M_metadata = defaultTermMetadata( term_role::lhs );
};

inline TermSummary
termSummary( std::size_t index, term_descriptor const& descriptor )
{
    auto const& metadata = descriptor.metadata();
    TermSummary summary;
    summary.index = index;
    summary.storageRole = descriptor.role;
    summary.label = std::string( descriptor.label );
    summary.name = metadata.name;
    summary.kind = metadata.kind;
    summary.support = metadata.support;
    summary.linearity = metadata.linearity;
    summary.role = metadata.role;
    summary.row = metadata.row;
    summary.column = metadata.column;
    summary.isAffine = metadata.isAffine;
    summary.isTimeDependent = metadata.isTimeDependent;
    summary.isParameterDependent = metadata.isParameterDependent;
    summary.isStateDependent = metadata.isStateDependent;
    summary.isSparseAssemblyCandidate = metadata.isSparseAssemblyCandidate;
    summary.isDeviceCandidate = metadata.isDeviceCandidate;
    summary.isMatrixFreeCandidate = metadata.isMatrixFreeCandidate;
    summary.isHdgCandidate = metadata.isHdgCandidate;
    summary.isTraceSchurCandidate = metadata.isTraceSchurCandidate;
    summary.isAdCandidate = metadata.isAdCandidate;
    summary.isRomCandidate = metadata.isRomCandidate;
    summary.isEnsembleCandidate = metadata.isEnsembleCandidate;
    summary.fields = metadata.fields;
    summary.parameters = metadata.parameters;
    summary.markers = metadata.markers;
    return summary;
}

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_TERMS_HPP
