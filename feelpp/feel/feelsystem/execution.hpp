/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_EXECUTION_HPP
#define FEELPP_FEELSYSTEM_EXECUTION_HPP 1

/**
 * @file execution.hpp
 * @brief Scheduler-neutral assembly task metadata and diagnostic task graphs.
 */

#include <cstddef>
#include <memory>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <feel/feelsystem/terms.hpp>

namespace Feel::systems
{

/**
 * @brief Resource categories used by scheduler-ready assembly task metadata.
 */
enum class AssemblyResourceKind
{
    Matrix,
    Vector,
    Field,
    Coefficient,
    Mesh,
    GeometryCache,
    TimeScheme,
    Backend,
    Scratch,
    Unknown
};

/**
 * @brief Human-readable name for an assembly resource kind.
 */
constexpr std::string_view
assemblyResourceKindName( AssemblyResourceKind kind ) noexcept
{
    switch ( kind )
    {
    case AssemblyResourceKind::Matrix:
        return "matrix";
    case AssemblyResourceKind::Vector:
        return "vector";
    case AssemblyResourceKind::Field:
        return "field";
    case AssemblyResourceKind::Coefficient:
        return "coefficient";
    case AssemblyResourceKind::Mesh:
        return "mesh";
    case AssemblyResourceKind::GeometryCache:
        return "geometry_cache";
    case AssemblyResourceKind::TimeScheme:
        return "time_scheme";
    case AssemblyResourceKind::Backend:
        return "backend";
    case AssemblyResourceKind::Scratch:
        return "scratch";
    case AssemblyResourceKind::Unknown:
        return "unknown";
    }
    return "unknown";
}

/**
 * @brief Scheduler access mode for one assembly resource.
 */
enum class AssemblyAccessMode
{
    Read,
    Write,
    ReadWrite,
    CommutativeAdd,
    AtomicAdd
};

/**
 * @brief Human-readable name for an assembly access mode.
 */
constexpr std::string_view
assemblyAccessModeName( AssemblyAccessMode mode ) noexcept
{
    switch ( mode )
    {
    case AssemblyAccessMode::Read:
        return "read";
    case AssemblyAccessMode::Write:
        return "write";
    case AssemblyAccessMode::ReadWrite:
        return "read_write";
    case AssemblyAccessMode::CommutativeAdd:
        return "commutative_add";
    case AssemblyAccessMode::AtomicAdd:
        return "atomic_add";
    }
    return "unknown";
}

/**
 * @brief Scheduler identifier stored in assembly execution options.
 */
enum class AssemblySchedulerKind
{
    Serial,
    Tbb,
    Specx,
    Starpu
};

/**
 * @brief Human-readable name for a scheduler kind.
 */
constexpr std::string_view
assemblySchedulerKindName( AssemblySchedulerKind kind ) noexcept
{
    switch ( kind )
    {
    case AssemblySchedulerKind::Serial:
        return "serial";
    case AssemblySchedulerKind::Tbb:
        return "tbb";
    case AssemblySchedulerKind::Specx:
        return "specx";
    case AssemblySchedulerKind::Starpu:
        return "starpu";
    }
    return "unknown";
}

/**
 * @brief Conservative range identity for scheduler-ready assembly metadata.
 *
 * Records the support, marker names, and optional partition label associated
 * with assembly task metadata.
 */
struct AssemblyRangeId
{
    /** @brief Geometric or algebraic support represented by the range. */
    TermSupport support = TermSupport::Unknown;
    /** @brief Optional marker names associated with the range. */
    std::vector<std::string> markers;
    /** @brief Optional stable label identifying a range partition. */
    std::string name;

    /**
     * @brief Return true when marker metadata is available.
     */
    bool hasMarkers() const noexcept { return !markers.empty(); }

    /**
     * @brief Return true when a range label is available.
     */
    bool hasName() const noexcept { return !name.empty(); }
};

/**
 * @brief Conservative degree-of-freedom target metadata.
 */
struct AssemblyDofRange
{
    /** @brief Optional block index associated with this dof range. */
    std::optional<int> block;
    /** @brief Algebraic role associated with this dof range. */
    TermRole role = TermRole::Unknown;
    /** @brief Optional stable label identifying a degree-of-freedom partition. */
    std::string name;

    /**
     * @brief Return true when a block index is available.
     */
    bool hasBlock() const noexcept { return static_cast<bool>( block ); }

    /**
     * @brief Return true when a dof-range label is available.
     */
    bool hasName() const noexcept { return !name.empty(); }
};

/**
 * @brief One resource read or written by a scheduler-ready assembly task.
 */
struct AssemblyResource
{
    /** @brief Resource category. */
    AssemblyResourceKind kind = AssemblyResourceKind::Unknown;
    /** @brief Algebraic role targeted by this resource. */
    TermRole role = TermRole::Unknown;
    /** @brief Optional block row targeted by this resource. */
    std::optional<int> row;
    /** @brief Optional block column targeted by this resource. */
    std::optional<int> column;
    /** @brief Range metadata associated with this resource access. */
    AssemblyRangeId range;
    /** @brief Row dof metadata associated with matrix-like resources. */
    AssemblyDofRange rowDofs;
    /** @brief Column dof metadata associated with matrix-like resources. */
    AssemblyDofRange columnDofs;
    /** @brief Scheduler access mode requested for this resource. */
    AssemblyAccessMode access = AssemblyAccessMode::Read;
    /** @brief Optional stable resource name. */
    std::string name;

    /**
     * @brief Return true when a block row index is available.
     */
    bool hasRow() const noexcept { return static_cast<bool>( row ); }

    /**
     * @brief Return true when a block column index is available.
     */
    bool hasColumn() const noexcept { return static_cast<bool>( column ); }

    /**
     * @brief Return true when a resource name is available.
     */
    bool hasName() const noexcept { return !name.empty(); }
};

/**
 * @brief Scheduler-facing metadata for one assembly task.
 */
struct AssemblyTaskInfo
{
    /** @brief Assembly batch id used to derive the task. */
    std::size_t batchId = 0;
    /** @brief Term ids carried by the originating assembly batch. */
    std::vector<std::size_t> termIds;
    /** @brief Geometric or algebraic support assembled by the task. */
    TermSupport support = TermSupport::Unknown;
    /** @brief Algebraic role assembled by the task. */
    TermRole role = TermRole::Unknown;
    /** @brief Optional block row targeted by the task. */
    std::optional<int> row;
    /** @brief Optional block column targeted by the task. */
    std::optional<int> column;
    /** @brief Strategy candidates inherited from the assembly batch. */
    AssemblyStrategyFlags strategy;
    /** @brief Conservative range identity associated with the task. */
    AssemblyRangeId range;
    /** @brief Conservative floating-point operation estimate, zero when unknown. */
    std::size_t estimatedFlops = 0;
    /** @brief Conservative byte-traffic estimate, zero when unknown. */
    std::size_t estimatedBytes = 0;
    /** @brief True when reductions should preserve deterministic merge policy. */
    bool deterministicReductionRequired = true;

    /**
     * @brief Return true when a block row index is available.
     */
    bool hasRow() const noexcept { return static_cast<bool>( row ); }

    /**
     * @brief Return true when a block column index is available.
     */
    bool hasColumn() const noexcept { return static_cast<bool>( column ); }

    /**
     * @brief Return true when either cost estimate is available.
     */
    bool hasCostEstimate() const noexcept { return estimatedFlops != 0 || estimatedBytes != 0; }
};

/**
 * @brief Scheduler-ready task descriptor derived from an assembly batch.
 */
struct AssemblyTask
{
    /** @brief Task metadata. */
    AssemblyTaskInfo info;
    /** @brief Resources read or written by the task. */
    std::vector<AssemblyResource> resources;

    /**
     * @brief Return true when the task has no resource descriptors.
     */
    bool empty() const noexcept { return resources.empty(); }
};

/**
 * @brief Side-effect-free task graph derived from an assembly plan.
 */
struct AssemblyTaskGraph
{
    /** @brief Scheduler-ready task descriptors. */
    std::vector<AssemblyTask> tasks;
    /** @brief Metadata diagnostics inherited from the assembly plan. */
    std::vector<TermConsistencyIssue> issues;

    /**
     * @brief Return true when the graph has no tasks.
     */
    bool empty() const noexcept { return tasks.empty(); }

    /**
     * @brief Number of task descriptors in the graph.
     */
    std::size_t nTasks() const noexcept { return tasks.size(); }

    /**
     * @brief Number of term ids referenced by all tasks.
     */
    std::size_t nTerms() const noexcept
    {
        std::size_t count = 0;
        for ( auto const& task : tasks )
            count += task.info.termIds.size();
        return count;
    }
};

/**
 * @brief Scheduler and resource preferences for assembly task descriptors.
 */
struct AssemblyExecutionOptions
{
    /** @brief Requested scheduler identifier. */
    AssemblySchedulerKind scheduler = AssemblySchedulerKind::Serial;
    /** @brief Optional batch grain size, zero when not specified. */
    std::size_t grainSize = 0;
    /** @brief True when deterministic reduction/merge behavior is required. */
    bool deterministic = true;
    /** @brief True when the execution policy permits device-capable batches. */
    bool allowDevice = false;

    /**
     * @brief Return true when a grain size has been requested.
     */
    bool hasGrainSize() const noexcept { return grainSize != 0; }
};

/**
 * @brief Scheduler-neutral interface for assembly task submission.
 *
 * Executors accept task descriptors through `submit()` and complete submitted
 * work through `wait()`. The recording implementation retains descriptors for
 * inspection; problem assembly uses the serial form replay paths.
 */
class AssemblyExecutor
{
  public:
    virtual ~AssemblyExecutor() = default;

    /**
     * @brief Submit one scheduler-ready task descriptor.
     */
    virtual void submit( AssemblyTask const& task ) = 0;

    /**
     * @brief Wait for all submitted descriptors to finish.
     */
    virtual void wait() = 0;
};

/**
 * @brief Recording executor used to validate scheduler-neutral task handoff.
 */
class RecordingAssemblyExecutor : public AssemblyExecutor
{
  public:
    /**
     * @brief Store a task descriptor without executing assembly.
     */
    void submit( AssemblyTask const& task ) override { M_tasks.push_back( task ); }

    /**
     * @brief Complete the no-op recording executor.
     */
    void wait() override {}

    /**
     * @brief Submitted task descriptors.
     */
    std::vector<AssemblyTask> const& tasks() const noexcept { return M_tasks; }

    /**
     * @brief Remove all recorded task descriptors.
     */
    void clear() { M_tasks.clear(); }

  private:
    /** @brief Task descriptors submitted to the recording executor. */
    std::vector<AssemblyTask> M_tasks;
};

namespace detail
{
/**
 * @brief Return true when support needs mesh and geometry resources.
 */
constexpr bool
assemblySupportNeedsGeometry( TermSupport support ) noexcept
{
    return support == TermSupport::Volume ||
           support == TermSupport::BoundaryFacet ||
           support == TermSupport::InteriorFacet ||
           support == TermSupport::Skeleton;
}

/**
 * @brief Target resource kind for an algebraic term role.
 */
constexpr AssemblyResourceKind
assemblyTargetResourceKind( TermRole role ) noexcept
{
    switch ( role )
    {
    case TermRole::Matrix:
    case TermRole::Jacobian:
    case TermRole::TimeDerivative:
    case TermRole::LocalCell:
    case TermRole::Trace:
    case TermRole::Coupling:
    case TermRole::Stabilization:
    case TermRole::Flux:
        return AssemblyResourceKind::Matrix;
    case TermRole::Vector:
    case TermRole::Residual:
        return AssemblyResourceKind::Vector;
    case TermRole::Constraint:
        return AssemblyResourceKind::Matrix;
    case TermRole::Unknown:
        return AssemblyResourceKind::Unknown;
    }
    return AssemblyResourceKind::Unknown;
}

/**
 * @brief Target access mode for an algebraic term role.
 */
constexpr AssemblyAccessMode
assemblyTargetAccessMode( TermRole role ) noexcept
{
    if ( role == TermRole::Constraint )
        return AssemblyAccessMode::ReadWrite;
    if ( role == TermRole::Unknown )
        return AssemblyAccessMode::Read;
    return AssemblyAccessMode::CommutativeAdd;
}

/**
 * @brief Return true when a resource kind is an algebraic assembly target.
 */
constexpr bool
isAssemblyTargetResource( AssemblyResourceKind kind ) noexcept
{
    return kind == AssemblyResourceKind::Matrix ||
           kind == AssemblyResourceKind::Vector;
}

/**
 * @brief Build conservative range metadata from one batch.
 */
inline AssemblyRangeId
assemblyRangeId( AssemblyBatch const& batch )
{
    AssemblyRangeId range;
    range.support = batch.support;
    range.markers = batch.markers;
    return range;
}

/**
 * @brief Build conservative row dof metadata from one batch.
 */
inline AssemblyDofRange
assemblyRowDofRange( AssemblyBatch const& batch )
{
    AssemblyDofRange dofs;
    dofs.block = batch.row;
    dofs.role = batch.role;
    return dofs;
}

/**
 * @brief Build conservative column dof metadata from one batch.
 */
inline AssemblyDofRange
assemblyColumnDofRange( AssemblyBatch const& batch )
{
    AssemblyDofRange dofs;
    dofs.block = batch.column;
    dofs.role = batch.role;
    return dofs;
}
} // namespace detail

/**
 * @brief Build a target resource descriptor for one assembly batch.
 */
inline AssemblyResource
assemblyTargetResource( AssemblyBatch const& batch )
{
    AssemblyResource resource;
    resource.kind = detail::assemblyTargetResourceKind( batch.role );
    resource.role = batch.role;
    resource.row = batch.row;
    resource.column = batch.column;
    resource.range = detail::assemblyRangeId( batch );
    resource.rowDofs = detail::assemblyRowDofRange( batch );
    resource.columnDofs = detail::assemblyColumnDofRange( batch );
    resource.access = detail::assemblyTargetAccessMode( batch.role );
    return resource;
}

/**
 * @brief Build scheduler-ready resource descriptors for one assembly batch.
 */
inline std::vector<AssemblyResource>
assemblyResources( AssemblyBatch const& batch )
{
    std::vector<AssemblyResource> resources;

    auto target = assemblyTargetResource( batch );
    if ( target.kind != AssemblyResourceKind::Unknown )
        resources.push_back( std::move( target ) );

    if ( batch.role == TermRole::Constraint )
    {
        AssemblyResource vectorTarget;
        vectorTarget.kind = AssemblyResourceKind::Vector;
        vectorTarget.role = batch.role;
        vectorTarget.row = batch.row;
        vectorTarget.column = batch.column;
        vectorTarget.range = detail::assemblyRangeId( batch );
        vectorTarget.rowDofs = detail::assemblyRowDofRange( batch );
        vectorTarget.access = AssemblyAccessMode::ReadWrite;
        resources.push_back( std::move( vectorTarget ) );
    }

    if ( detail::assemblySupportNeedsGeometry( batch.support ) )
    {
        AssemblyResource mesh;
        mesh.kind = AssemblyResourceKind::Mesh;
        mesh.range = detail::assemblyRangeId( batch );
        mesh.access = AssemblyAccessMode::Read;
        resources.push_back( std::move( mesh ) );

        AssemblyResource geometry;
        geometry.kind = AssemblyResourceKind::GeometryCache;
        geometry.range = detail::assemblyRangeId( batch );
        geometry.access = AssemblyAccessMode::Read;
        resources.push_back( std::move( geometry ) );
    }

    if ( batch.support == TermSupport::Time )
    {
        AssemblyResource timeScheme;
        timeScheme.kind = AssemblyResourceKind::TimeScheme;
        timeScheme.range = detail::assemblyRangeId( batch );
        timeScheme.access = AssemblyAccessMode::Read;
        resources.push_back( std::move( timeScheme ) );
    }

    for ( auto const& field : batch.fields )
    {
        AssemblyResource resource;
        resource.kind = AssemblyResourceKind::Field;
        resource.role = batch.role;
        resource.row = batch.row;
        resource.column = batch.column;
        resource.range = detail::assemblyRangeId( batch );
        resource.access = AssemblyAccessMode::Read;
        resource.name = field;
        resources.push_back( std::move( resource ) );
    }

    for ( auto const& parameter : batch.parameters )
    {
        AssemblyResource resource;
        resource.kind = AssemblyResourceKind::Coefficient;
        resource.role = batch.role;
        resource.row = batch.row;
        resource.column = batch.column;
        resource.range = detail::assemblyRangeId( batch );
        resource.access = AssemblyAccessMode::Read;
        resource.name = parameter;
        resources.push_back( std::move( resource ) );
    }

    return resources;
}

/**
 * @brief Build scheduler-facing task metadata for one assembly batch.
 */
inline AssemblyTaskInfo
assemblyTaskInfo( std::size_t batchId,
                  AssemblyBatch const& batch,
                  AssemblyExecutionOptions const& options = {} )
{
    AssemblyTaskInfo info;
    info.batchId = batchId;
    info.termIds = batch.termIds;
    info.support = batch.support;
    info.role = batch.role;
    info.row = batch.row;
    info.column = batch.column;
    info.strategy = batch.strategy;
    info.range = detail::assemblyRangeId( batch );
    info.deterministicReductionRequired =
        options.deterministic &&
        detail::isAssemblyTargetResource( detail::assemblyTargetResourceKind( batch.role ) );
    return info;
}

/**
 * @brief Build one scheduler-ready task from one assembly batch.
 */
inline AssemblyTask
assemblyTask( std::size_t batchId,
              AssemblyBatch const& batch,
              AssemblyExecutionOptions const& options = {} )
{
    AssemblyTask task;
    task.info = assemblyTaskInfo( batchId, batch, options );
    task.resources = assemblyResources( batch );
    return task;
}

/**
 * @brief Build a scheduler-ready task graph from an assembly plan.
 */
inline AssemblyTaskGraph
buildAssemblyTaskGraph( AssemblyPlan const& plan,
                        AssemblyExecutionOptions const& options = {} )
{
    AssemblyTaskGraph graph;
    graph.issues = plan.issues;
    graph.tasks.reserve( plan.batches.size() );
    for ( std::size_t i = 0; i != plan.batches.size(); ++i )
        graph.tasks.push_back( assemblyTask( i, plan.batches[i], options ) );
    return graph;
}

/**
 * @brief Build a scheduler-ready task graph directly from term descriptors.
 */
template <typename TermContainer>
AssemblyTaskGraph
buildAssemblyTaskGraph( TermContainer const& terms,
                        AssemblyPlanOptions const& planOptions = {},
                        AssemblyExecutionOptions const& executionOptions = {} )
{
    return buildAssemblyTaskGraph( buildAssemblyPlan( terms, planOptions ), executionOptions );
}

/**
 * @brief Create a recording executor for inspecting assembly task submission.
 */
inline std::unique_ptr<AssemblyExecutor>
makeAssemblyExecutor( AssemblyExecutionOptions const& )
{
    return std::make_unique<RecordingAssemblyExecutor>();
}

/**
 * @brief Write one compact range id description.
 */
inline std::ostream&
describeAssemblyRangeId( std::ostream& os, AssemblyRangeId const& range )
{
    os << "support=" << termSupportName( range.support );
    if ( range.hasName() )
        os << " name=" << range.name;
    if ( range.hasMarkers() )
        os << " markers=" << detail::joinDiagnosticStrings( range.markers );
    return os;
}

/**
 * @brief Write one compact resource description.
 */
inline std::ostream&
describeAssemblyResource( std::ostream& os, AssemblyResource const& resource )
{
    os << "resource=" << assemblyResourceKindName( resource.kind )
       << " access=" << assemblyAccessModeName( resource.access );
    if ( resource.role != TermRole::Unknown )
        os << " role=" << termRoleName( resource.role );
    if ( resource.hasRow() )
        os << " row=" << *resource.row;
    if ( resource.hasColumn() )
        os << " col=" << *resource.column;
    if ( resource.hasName() )
        os << " name=" << resource.name;
    os << " ";
    describeAssemblyRangeId( os, resource.range );
    return os;
}

/**
 * @brief Write one compact task description.
 */
inline std::ostream&
describeAssemblyTask( std::ostream& os, std::size_t index, AssemblyTask const& task )
{
    os << "  [" << index << "]"
       << " batch=" << task.info.batchId
       << " support=" << termSupportName( task.info.support )
       << " role=" << termRoleName( task.info.role )
       << " terms=" << detail::joinDiagnosticIds( task.info.termIds );
    if ( task.info.hasRow() )
        os << " row=" << *task.info.row;
    if ( task.info.hasColumn() )
        os << " col=" << *task.info.column;
    if ( task.info.deterministicReductionRequired )
        os << " deterministicReduction=true";
    if ( task.info.hasCostEstimate() )
        os << " estimatedFlops=" << task.info.estimatedFlops
           << " estimatedBytes=" << task.info.estimatedBytes;
    os << '\n';

    for ( auto const& resource : task.resources )
    {
        os << "    ";
        describeAssemblyResource( os, resource ) << '\n';
    }
    return os;
}

/**
 * @brief Write a side-effect-free scheduler-ready task graph description.
 */
inline std::ostream&
describeAssemblyTaskGraph( std::ostream& os,
                           std::string_view problemName,
                           AssemblyTaskGraph const& graph )
{
    os << "AssemblyTaskGraph: " << ( problemName.empty() ? std::string_view{ "<unnamed>" } : problemName ) << '\n';
    os << "Tasks:\n";
    for ( std::size_t i = 0; i != graph.tasks.size(); ++i )
        describeAssemblyTask( os, i, graph.tasks[i] );

    if ( !graph.issues.empty() )
    {
        os << "Issues:\n";
        for ( auto const& issue : graph.issues )
            os << "  [" << issue.index << "] " << issue.severity << ": " << issue.message << '\n';
    }
    return os;
}

/**
 * @brief Return a side-effect-free task graph description as a string.
 */
inline std::string
describeAssemblyTaskGraph( std::string_view problemName,
                           AssemblyTaskGraph const& graph )
{
    std::ostringstream os;
    describeAssemblyTaskGraph( os, problemName, graph );
    return os.str();
}

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_EXECUTION_HPP
