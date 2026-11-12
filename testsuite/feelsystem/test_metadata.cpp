/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file
 * @brief FeelSystem regression test.
 */

#define BOOST_TEST_MODULE feelsystem_metadata
#include <feel/feelcore/testsuite.hpp>

#include <boost/hana/integral_constant.hpp>

#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
bool contains( std::string const& text, std::string_view needle )
{
    return text.find( needle ) != std::string::npos;
}

template <typename Predicate>
bool hasAssemblyBatch( Feel::systems::AssemblyPlan const& plan, Predicate predicate )
{
    for ( auto const& batch : plan.batches )
    {
        if ( predicate( batch ) )
            return true;
    }
    return false;
}

template <typename Predicate>
bool hasAssemblyTask( Feel::systems::AssemblyTaskGraph const& graph, Predicate predicate )
{
    for ( auto const& task : graph.tasks )
    {
        if ( predicate( task ) )
            return true;
    }
    return false;
}

template <typename Predicate>
bool hasAssemblyResource( Feel::systems::AssemblyTask const& task, Predicate predicate )
{
    for ( auto const& resource : task.resources )
    {
        if ( predicate( resource ) )
            return true;
    }
    return false;
}
} // namespace

BOOST_AUTO_TEST_SUITE( feelsystem_metadata_suite )

BOOST_AUTO_TEST_CASE( assembly_context_exposes_mode_role_indices_and_typed_storage )
{
    Feel::systems::AssemblyContext context;
    int problem = 7;
    int space = 11;
    int timeStep = 13;
    int state = 17;
    int backend = 19;
    int timeScheme = 23;
    int matrix = 29;
    int vector = 31;
    int bilinearForm = 37;
    int linearForm = 41;
    int blockBilinearForm = 43;
    int blockLinearForm = 47;

    context.mode = Feel::systems::AssemblyMode::Mixed;
    context.role = Feel::systems::term_role::lhs;
    context.row = 2;
    context.column = 3;
    context.problem = &problem;
    context.space = &space;
    context.timeStep = &timeStep;
    context.state = &state;
    context.backend = &backend;
    context.timeScheme = &timeScheme;
    context.matrix = &matrix;
    context.vector = &vector;
    context.bilinearForm = &bilinearForm;
    context.linearForm = &linearForm;
    context.blockBilinearForm = &blockBilinearForm;
    context.blockLinearForm = &blockLinearForm;

    BOOST_CHECK_EQUAL( std::string( Feel::systems::assemblyModeName( context.mode ) ), "mixed" );
    BOOST_CHECK( context.hasRole() );
    BOOST_CHECK( context.hasRow() );
    BOOST_CHECK( context.hasColumn() );
    BOOST_CHECK( context.hasTimeStep() );
    BOOST_CHECK( context.hasState() );
    BOOST_CHECK( *context.role == Feel::systems::term_role::lhs );
    BOOST_CHECK_EQUAL( *context.row, 2 );
    BOOST_CHECK_EQUAL( *context.column, 3 );
    BOOST_CHECK_EQUAL( context.problemAs<int>(), 7 );
    BOOST_CHECK_EQUAL( context.spaceAs<int>(), 11 );
    BOOST_CHECK_EQUAL( context.timeStepAs<int>(), 13 );
    BOOST_CHECK_EQUAL( context.stateAs<int>(), 17 );
    BOOST_CHECK_EQUAL( context.backendAs<int>(), 19 );
    BOOST_CHECK_EQUAL( context.timeSchemeAs<int>(), 23 );
    BOOST_CHECK_EQUAL( context.matrixAs<int>(), 29 );
    BOOST_CHECK_EQUAL( context.vectorAs<int>(), 31 );
    BOOST_CHECK_EQUAL( context.bilinearFormAs<int>(), 37 );
    BOOST_CHECK_EQUAL( context.linearFormAs<int>(), 41 );
    BOOST_CHECK_EQUAL( context.blockBilinearFormAs<int>(), 43 );
    BOOST_CHECK_EQUAL( context.blockLinearFormAs<int>(), 47 );
}

BOOST_AUTO_TEST_CASE( scalar_linear_terms_expose_conservative_metadata )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto problem =
        linearProblem( _space = Xh, _name = "metadata-linear" )
            .lhs( integrate( _range = elements( mesh ),
                             _expr = inner( grad( u ), grad( v ) ) ) )
            .rhs( integrate( _range = elements( mesh ),
                             _expr = v ) );

    BOOST_CHECK_EQUAL( problem.lhsTerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.rhsTerms().size(), 1 );

    auto const& lhsMetadata = problem.lhsTerms().metadata( 0 );
    BOOST_CHECK( lhsMetadata.kind == systems::TermKind::Lhs );
    BOOST_CHECK( lhsMetadata.role == systems::TermRole::Matrix );
    BOOST_CHECK( lhsMetadata.linearity == systems::TermLinearity::Bilinear );
    BOOST_CHECK( lhsMetadata.isAffine );
    BOOST_CHECK_EQUAL( lhsMetadata.name, "lhs" );

    auto const& rhsMetadata = problem.rhsTerms()[0].metadata();
    BOOST_CHECK( rhsMetadata.kind == systems::TermKind::Rhs );
    BOOST_CHECK( rhsMetadata.role == systems::TermRole::Vector );
    BOOST_CHECK( rhsMetadata.linearity == systems::TermLinearity::Linear );
    BOOST_CHECK( rhsMetadata.isAffine );
    BOOST_CHECK_EQUAL( rhsMetadata.name, "rhs" );
}

BOOST_AUTO_TEST_CASE( linear_problem_describe_exposes_stable_term_fields_without_assembly )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto problem =
        linearProblem( _space = Xh, _name = "metadata-linear-diagnostics" )
            .lhs( integrate( _range = elements( mesh ),
                             _expr = inner( grad( u ), grad( v ) ) ) )
            .rhs( integrate( _range = elements( mesh ),
                             _expr = v ) );

    auto summaries = problem.termSummaries();
    BOOST_REQUIRE_EQUAL( summaries.size(), 2 );
    BOOST_CHECK_EQUAL( summaries[0].index, 0 );
    BOOST_CHECK( summaries[0].kind == systems::TermKind::Lhs );
    BOOST_CHECK( summaries[0].role == systems::TermRole::Matrix );
    BOOST_CHECK( summaries[1].kind == systems::TermKind::Rhs );
    BOOST_CHECK( summaries[1].role == systems::TermRole::Vector );
    BOOST_CHECK( problem.checkTerms().empty() );

    std::ostringstream os;
    problem.describe( os );
    auto const text = os.str();
    BOOST_CHECK( contains( text, "Problem: metadata-linear-diagnostics" ) );
    BOOST_CHECK( contains( text, "Terms:" ) );
    BOOST_CHECK( contains( text, "kind=lhs" ) );
    BOOST_CHECK( contains( text, "role=matrix" ) );
    BOOST_CHECK( contains( text, "linearity=bilinear" ) );
    BOOST_CHECK( contains( text, "kind=rhs" ) );
    BOOST_CHECK( contains( text, "role=vector" ) );
    BOOST_CHECK( !contains( text, "Issues:" ) );
    BOOST_CHECK( contains( problem.describe(), "metadata-linear-diagnostics" ) );
}

BOOST_AUTO_TEST_CASE( assembly_plan_strategy_flags_and_order_are_metadata_driven )
{
    using namespace Feel::systems;

    auto matrixMetadata = defaultTermMetadata( term_role::lhs );
    matrixMetadata.name = "matrix-strategy";
    matrixMetadata.support = TermSupport::Volume;
    matrixMetadata.isDeviceCandidate = true;
    matrixMetadata.isMatrixFreeCandidate = true;
    matrixMetadata.isRomCandidate = true;

    auto vectorMetadata = defaultTermMetadata( term_role::rhs );
    vectorMetadata.name = "vector-strategy";
    vectorMetadata.support = TermSupport::Volume;
    vectorMetadata.isSparseAssemblyCandidate = false;
    vectorMetadata.isEnsembleCandidate = true;

    std::vector<term_descriptor> firstTerms;
    firstTerms.emplace_back( term_role::rhs, "rhs", vectorMetadata );
    firstTerms.emplace_back( term_role::lhs, "lhs", matrixMetadata );

    std::vector<term_descriptor> secondTerms;
    secondTerms.emplace_back( term_role::lhs, "lhs", matrixMetadata );
    secondTerms.emplace_back( term_role::rhs, "rhs", vectorMetadata );

    auto const firstPlan = buildAssemblyPlan( firstTerms );
    auto const secondPlan = buildAssemblyPlan( secondTerms );

    BOOST_REQUIRE_EQUAL( firstPlan.batches.size(), 2 );
    BOOST_REQUIRE_EQUAL( secondPlan.batches.size(), 2 );
    for ( std::size_t i = 0; i != firstPlan.batches.size(); ++i )
    {
        BOOST_CHECK( firstPlan.batches[i].support == secondPlan.batches[i].support );
        BOOST_CHECK( firstPlan.batches[i].role == secondPlan.batches[i].role );
    }

    BOOST_CHECK( hasAssemblyBatch( firstPlan,
                                   []( auto const& batch )
                                   {
                                       return batch.role == TermRole::Matrix &&
                                              batch.strategy.sparseCpuCandidate &&
                                              batch.strategy.deviceCandidate &&
                                              batch.strategy.matrixFreeCandidate &&
                                              batch.strategy.romCandidate;
                                   } ) );
    BOOST_CHECK( hasAssemblyBatch( firstPlan,
                                   []( auto const& batch )
                                   {
                                       return batch.role == TermRole::Vector &&
                                              !batch.strategy.sparseCpuCandidate &&
                                              batch.strategy.ensembleCandidate;
                                   } ) );
}

BOOST_AUTO_TEST_CASE( assembly_task_graph_resources_are_scheduler_ready_metadata )
{
    using namespace Feel::systems;

    AssemblyExecutionOptions defaultOptions;
    BOOST_CHECK( defaultOptions.scheduler == AssemblySchedulerKind::Serial );
    BOOST_CHECK( defaultOptions.deterministic );
    BOOST_CHECK( !defaultOptions.hasGrainSize() );
    BOOST_CHECK_EQUAL( std::string( assemblySchedulerKindName( defaultOptions.scheduler ) ), "serial" );

    auto matrixMetadata = defaultTermMetadata( term_role::lhs );
    matrixMetadata.name = "matrix-resource";
    matrixMetadata.support = TermSupport::BoundaryFacet;
    matrixMetadata.fields = { "u" };
    matrixMetadata.parameters = { "mu" };
    matrixMetadata.markers = { "Gamma" };

    auto constraintMetadata = defaultTermMetadata( term_role::constraint );
    constraintMetadata.name = "constraint-resource";

    std::vector<term_descriptor> terms;
    terms.emplace_back( term_role::lhs, "lhs", matrixMetadata );
    terms.emplace_back( term_role::constraint, "constraint", constraintMetadata );

    auto const graph = buildAssemblyTaskGraph( terms );
    BOOST_REQUIRE_EQUAL( graph.nTasks(), 2 );
    BOOST_REQUIRE_EQUAL( graph.nTerms(), 2 );
    BOOST_CHECK( graph.issues.empty() );

    BOOST_CHECK( hasAssemblyTask( graph,
                                  []( auto const& task )
                                  {
                                      return task.info.role == TermRole::Matrix &&
                                             task.info.deterministicReductionRequired &&
                                             task.info.range.hasMarkers() &&
                                             hasAssemblyResource( task,
                                                                  []( auto const& resource )
                                                                  {
                                                                      return resource.kind == AssemblyResourceKind::Matrix &&
                                                                             resource.access == AssemblyAccessMode::CommutativeAdd;
                                                                  } ) &&
                                             hasAssemblyResource( task,
                                                                  []( auto const& resource )
                                                                  {
                                                                      return resource.kind == AssemblyResourceKind::Mesh &&
                                                                             resource.access == AssemblyAccessMode::Read;
                                                                  } ) &&
                                             hasAssemblyResource( task,
                                                                  []( auto const& resource )
                                                                  {
                                                                      return resource.kind == AssemblyResourceKind::GeometryCache &&
                                                                             resource.access == AssemblyAccessMode::Read;
                                                                  } ) &&
                                             hasAssemblyResource( task,
                                                                  []( auto const& resource )
                                                                  {
                                                                      return resource.kind == AssemblyResourceKind::Field &&
                                                                             resource.name == "u" &&
                                                                             resource.access == AssemblyAccessMode::Read;
                                                                  } ) &&
                                             hasAssemblyResource( task,
                                                                  []( auto const& resource )
                                                                  {
                                                                      return resource.kind == AssemblyResourceKind::Coefficient &&
                                                                             resource.name == "mu" &&
                                                                             resource.access == AssemblyAccessMode::Read;
                                                                  } );
                                  } ) );

    BOOST_CHECK( hasAssemblyTask( graph,
                                  []( auto const& task )
                                  {
                                      return task.info.role == TermRole::Constraint &&
                                             hasAssemblyResource( task,
                                                                  []( auto const& resource )
                                                                  {
                                                                      return resource.kind == AssemblyResourceKind::Matrix &&
                                                                             resource.access == AssemblyAccessMode::ReadWrite;
                                                                  } ) &&
                                             hasAssemblyResource( task,
                                                                  []( auto const& resource )
                                                                  {
                                                                      return resource.kind == AssemblyResourceKind::Vector &&
                                                                             resource.access == AssemblyAccessMode::ReadWrite;
                                                                  } );
                                  } ) );

    AssemblyExecutionOptions nondeterministicOptions;
    nondeterministicOptions.deterministic = false;
    auto const nondeterministicGraph = buildAssemblyTaskGraph( buildAssemblyPlan( terms ), nondeterministicOptions );
    BOOST_CHECK( hasAssemblyTask( nondeterministicGraph,
                                  []( auto const& task )
                                  {
                                      return task.info.role == TermRole::Matrix &&
                                             !task.info.deterministicReductionRequired;
                                  } ) );

    RecordingAssemblyExecutor recorder;
    for ( auto const& task : graph.tasks )
        recorder.submit( task );
    recorder.wait();
    BOOST_CHECK_EQUAL( recorder.tasks().size(), graph.tasks.size() );

    auto executor = makeAssemblyExecutor( defaultOptions );
    auto* recordingExecutor = dynamic_cast<RecordingAssemblyExecutor*>( executor.get() );
    BOOST_REQUIRE( recordingExecutor != nullptr );
    for ( auto const& task : graph.tasks )
        executor->submit( task );
    executor->wait();
    BOOST_CHECK_EQUAL( recordingExecutor->tasks().size(), graph.tasks.size() );

    auto const text = describeAssemblyTaskGraph( "metadata-task-graph", graph );
    BOOST_CHECK( contains( text, "AssemblyTaskGraph: metadata-task-graph" ) );
    BOOST_CHECK( contains( text, "resource=matrix" ) );
    BOOST_CHECK( contains( text, "access=commutative_add" ) );
    BOOST_CHECK( contains( text, "markers=Gamma" ) );
}

BOOST_AUTO_TEST_CASE( linear_problem_assembly_plan_is_side_effect_free )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto problem =
        linearProblem( _space = Xh, _name = "metadata-linear-plan" )
            .lhs( integrate( _range = elements( mesh ),
                             _expr = inner( grad( u ), grad( v ) ) ) )
            .rhs( integrate( _range = elements( mesh ),
                             _expr = v ) );

    auto const plan = problem.assemblyPlan();
    BOOST_REQUIRE_EQUAL( plan.nTerms(), 2 );
    BOOST_REQUIRE_EQUAL( plan.batches.size(), 2 );
    BOOST_CHECK( plan.issues.empty() );
    BOOST_CHECK( hasAssemblyBatch( plan,
                                   []( auto const& batch )
                                   {
                                       return batch.role == systems::TermRole::Matrix &&
                                              batch.strategy.sparseCpuCandidate;
                                   } ) );
    BOOST_CHECK( hasAssemblyBatch( plan,
                                   []( auto const& batch )
                                   {
                                       return batch.role == systems::TermRole::Vector &&
                                              batch.strategy.sparseCpuCandidate;
                                   } ) );

    auto const text = problem.describeAssemblyPlan();
    BOOST_CHECK( contains( text, "AssemblyPlan: metadata-linear-plan" ) );
    BOOST_CHECK( contains( text, "Batches:" ) );
    BOOST_CHECK( contains( text, "sparseCpu=true" ) );

    auto const graph = problem.assemblyTaskGraph();
    BOOST_REQUIRE_EQUAL( graph.nTerms(), 2 );
    BOOST_REQUIRE_EQUAL( graph.nTasks(), 2 );
    BOOST_CHECK( graph.issues.empty() );
    BOOST_CHECK( hasAssemblyTask( graph,
                                  []( auto const& task )
                                  {
                                      return task.info.role == systems::TermRole::Matrix &&
                                             hasAssemblyResource( task,
                                                                  []( auto const& resource )
                                                                  {
                                                                      return resource.kind == systems::AssemblyResourceKind::Matrix &&
                                                                             resource.access == systems::AssemblyAccessMode::CommutativeAdd;
                                                                  } );
                                  } ) );
    BOOST_CHECK( hasAssemblyTask( graph,
                                  []( auto const& task )
                                  {
                                      return task.info.role == systems::TermRole::Vector &&
                                             hasAssemblyResource( task,
                                                                  []( auto const& resource )
                                                                  {
                                                                      return resource.kind == systems::AssemblyResourceKind::Vector &&
                                                                             resource.access == systems::AssemblyAccessMode::CommutativeAdd;
                                                                  } );
                                  } ) );
    BOOST_CHECK( contains( problem.describeAssemblyTaskGraph(), "AssemblyTaskGraph: metadata-linear-plan" ) );
}

BOOST_AUTO_TEST_CASE( nonlinear_affine_and_state_terms_are_distinguishable )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto problem =
        nonlinearProblem( _space = Xh, _name = "metadata-nonlinear" )
            .reaction( integrate( _range = elements( mesh ),
                                  _expr = u * v ) )
            .source( integrate( _range = elements( mesh ),
                                _expr = v ) )
            .residual( []( auto const&, auto& ) {} )
            .jacobian( []( auto const&, auto& ) {} );

    auto const& reactionMetadata = problem.lhsTerms().metadata( 0 );
    BOOST_CHECK( reactionMetadata.kind == systems::TermKind::Reaction );
    BOOST_CHECK( reactionMetadata.role == systems::TermRole::Matrix );
    BOOST_CHECK( reactionMetadata.isAffine );
    BOOST_CHECK_EQUAL( reactionMetadata.name, "reaction" );

    auto const& sourceMetadata = problem.rhsTerms().metadata( 0 );
    BOOST_CHECK( sourceMetadata.kind == systems::TermKind::Source );
    BOOST_CHECK( sourceMetadata.role == systems::TermRole::Vector );
    BOOST_CHECK( sourceMetadata.isAffine );
    BOOST_CHECK_EQUAL( sourceMetadata.name, "source" );

    auto const& residualMetadata = problem.residualTerms().metadata( 0 );
    BOOST_CHECK( residualMetadata.kind == systems::TermKind::Residual );
    BOOST_CHECK( residualMetadata.role == systems::TermRole::Residual );
    BOOST_CHECK( residualMetadata.linearity == systems::TermLinearity::Nonlinear );
    BOOST_CHECK( residualMetadata.isStateDependent );

    auto const& jacobianMetadata = problem.jacobianTerms().metadata( 0 );
    BOOST_CHECK( jacobianMetadata.kind == systems::TermKind::Jacobian );
    BOOST_CHECK( jacobianMetadata.role == systems::TermRole::Jacobian );
    BOOST_CHECK( jacobianMetadata.linearity == systems::TermLinearity::Nonlinear );
    BOOST_CHECK( jacobianMetadata.isStateDependent );
}

BOOST_AUTO_TEST_CASE( nonlinear_diagnostics_keep_affine_and_state_terms_visible )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto problem =
        nonlinearProblem( _space = Xh, _name = "metadata-nonlinear-diagnostics" )
            .reaction( integrate( _range = elements( mesh ),
                                  _expr = u * v ) )
            .source( integrate( _range = elements( mesh ),
                                _expr = v ) )
            .residual( []( auto const&, auto& ) {} )
            .jacobian( []( auto const&, auto& ) {} );

    BOOST_CHECK( problem.checkTerms().empty() );
    auto const text = problem.describe();
    BOOST_CHECK( contains( text, "kind=reaction" ) );
    BOOST_CHECK( contains( text, "kind=source" ) );
    BOOST_CHECK( contains( text, "kind=residual" ) );
    BOOST_CHECK( contains( text, "kind=jacobian" ) );
    BOOST_CHECK( contains( text, "stateDependent=true" ) );
}

BOOST_AUTO_TEST_CASE( nonlinear_problem_assembly_plan_separates_residual_and_jacobian )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto problem =
        nonlinearProblem( _space = Xh, _name = "metadata-nonlinear-plan" )
            .reaction( integrate( _range = elements( mesh ),
                                  _expr = u * v ) )
            .source( integrate( _range = elements( mesh ),
                                _expr = v ) )
            .residual( []( auto const&, auto& ) {} )
            .jacobian( []( auto const&, auto& ) {} );

    auto const plan = problem.assemblyPlan();
    BOOST_REQUIRE_EQUAL( plan.nTerms(), 4 );
    BOOST_CHECK( plan.issues.empty() );
    BOOST_CHECK( hasAssemblyBatch( plan,
                                   []( auto const& batch )
                                   { return batch.role == systems::TermRole::Residual; } ) );
    BOOST_CHECK( hasAssemblyBatch( plan,
                                   []( auto const& batch )
                                   { return batch.role == systems::TermRole::Jacobian; } ) );
    BOOST_CHECK( contains( problem.describeAssemblyPlan(), "role=jacobian" ) );
}

BOOST_AUTO_TEST_CASE( picard_evolution_and_constraint_terms_are_classified )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto picardProblem =
        nonlinearProblem( _space = Xh, _name = "metadata-picard" )
            .iteration( []( auto const&, auto&, auto& ) {} );

    auto const& iterationMetadata = picardProblem.iterationTerms().metadata( 0 );
    BOOST_CHECK( iterationMetadata.kind == systems::TermKind::PicardIteration );
    BOOST_CHECK( iterationMetadata.role == systems::TermRole::Matrix );
    BOOST_CHECK( iterationMetadata.isStateDependent );

    auto evolutionProblem =
        parabolicProblem( _space = Xh, _name = "metadata-evolution" )
            .mass( integrate( _range = elements( mesh ),
                              _expr = u * v ) )
            .stiffness( integrate( _range = elements( mesh ),
                                   _expr = inner( grad( u ), grad( v ) ) ) );

    auto const& massMetadata = evolutionProblem.massTerms().metadata( 0 );
    BOOST_CHECK( massMetadata.kind == systems::TermKind::Mass );
    BOOST_CHECK( massMetadata.role == systems::TermRole::TimeDerivative );
    BOOST_CHECK( massMetadata.isAffine );
    BOOST_CHECK( massMetadata.isTimeDependent );

    auto const& stiffnessMetadata = evolutionProblem.stiffnessTerms().metadata( 0 );
    BOOST_CHECK( stiffnessMetadata.kind == systems::TermKind::Stiffness );
    BOOST_CHECK( stiffnessMetadata.role == systems::TermRole::Matrix );

    auto constrainedProblem =
        linearProblem( _space = Xh, _name = "metadata-constraint" )
            .constraint( []( auto&, auto& ) {} );

    auto const& constraintMetadata = constrainedProblem.constraintTerms().metadata( 0 );
    BOOST_CHECK( constraintMetadata.kind == systems::TermKind::Constraint );
    BOOST_CHECK( constraintMetadata.role == systems::TermRole::Constraint );
    BOOST_CHECK( constraintMetadata.support == systems::TermSupport::Algebraic );
}

BOOST_AUTO_TEST_CASE( evolution_and_picard_diagnostics_are_side_effect_free )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto picardProblem =
        nonlinearProblem( _space = Xh, _name = "metadata-picard-diagnostics" )
            .iteration( []( auto const&, auto&, auto& ) {} );
    BOOST_CHECK( picardProblem.checkTerms().empty() );
    BOOST_CHECK( contains( picardProblem.describe(), "kind=picard_iteration" ) );
    BOOST_CHECK( contains( picardProblem.describe(), "stateDependent=true" ) );

    auto evolutionProblem =
        parabolicProblem( _space = Xh, _name = "metadata-evolution-diagnostics" )
            .mass( integrate( _range = elements( mesh ),
                              _expr = u * v ) )
            .stiffness( integrate( _range = elements( mesh ),
                                   _expr = inner( grad( u ), grad( v ) ) ) );
    BOOST_CHECK( evolutionProblem.checkTerms().empty() );
    auto const text = evolutionProblem.describe();
    BOOST_CHECK( contains( text, "kind=mass" ) );
    BOOST_CHECK( contains( text, "role=time_derivative" ) );
    BOOST_CHECK( contains( text, "timeDependent=true" ) );
    BOOST_CHECK( contains( text, "kind=stiffness" ) );
}

BOOST_AUTO_TEST_CASE( evolution_problem_assembly_plan_keeps_time_derivative_role )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto evolutionProblem =
        parabolicProblem( _space = Xh, _name = "metadata-evolution-plan" )
            .mass( integrate( _range = elements( mesh ),
                              _expr = u * v ) )
            .stiffness( integrate( _range = elements( mesh ),
                                   _expr = inner( grad( u ), grad( v ) ) ) );

    auto const plan = evolutionProblem.assemblyPlan();
    BOOST_REQUIRE_EQUAL( plan.nTerms(), 2 );
    BOOST_CHECK( hasAssemblyBatch( plan,
                                   []( auto const& batch )
                                   { return batch.role == systems::TermRole::TimeDerivative; } ) );
    BOOST_CHECK( contains( evolutionProblem.describeAssemblyPlan(), "role=time_derivative" ) );
}

BOOST_AUTO_TEST_CASE( natural_boundary_terms_keep_boundary_support )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto problem =
        linearProblem( _space = Xh, _name = "metadata-boundary" )
            .neumann( _range = boundaryfaces( mesh ),
                      _test = v,
                      _expr = cst( 2.0 ) )
            .robin( _range = boundaryfaces( mesh ),
                    _trial = u,
                    _test = v,
                    _alpha = cst( 3.0 ),
                    _expr = cst( 4.0 ) )
            .fourier( _range = boundaryfaces( mesh ),
                      _trial = u,
                      _test = v,
                      _alpha = cst( 5.0 ),
                      _expr = cst( 6.0 ) );

    BOOST_CHECK_EQUAL( problem.lhsTerms().size(), 2 );
    BOOST_CHECK_EQUAL( problem.rhsTerms().size(), 3 );

    BOOST_CHECK( problem.rhsTerms().metadata( 0 ).kind == systems::TermKind::BoundaryCondition );
    BOOST_CHECK( problem.rhsTerms().metadata( 0 ).support == systems::TermSupport::BoundaryFacet );
    BOOST_CHECK_EQUAL( problem.rhsTerms().metadata( 0 ).name, "neumann" );
    BOOST_CHECK_EQUAL( problem.lhsTerms().metadata( 0 ).name, "robin" );
    BOOST_CHECK_EQUAL( problem.lhsTerms().metadata( 1 ).name, "fourier" );
}

BOOST_AUTO_TEST_CASE( block_terms_expose_row_and_column_metadata )
{
    using namespace Feel;
    using namespace Feel::vf;
    using namespace boost::hana::literals;

    auto mesh = unitSquare();
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = product( Vh, Qh );

    auto u = trial( Vh );
    auto v = test( Vh );
    auto p = trial( Qh );
    auto q = test( Qh );

    auto problem =
        mixedProblem( _space = Xh, _name = "metadata-block" )
            .block( 0_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = inner( grad( u ), grad( v ) ) ) )
            .block( 0_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = -p * div( v ) ) )
            .row( 1_c,
                  integrate( _range = elements( mesh ),
                             _expr = q ) )
            .constraint( []( auto&, auto& ) {} );

    BOOST_CHECK_EQUAL( problem.blockTerms().size(), 2 );
    BOOST_CHECK_EQUAL( problem.rowTerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.constraintTerms().size(), 1 );

    auto const& firstBlock = problem.blockTerms().metadata( 0 );
    BOOST_CHECK( firstBlock.kind == systems::TermKind::Block );
    BOOST_REQUIRE( firstBlock.hasRow() );
    BOOST_REQUIRE( firstBlock.hasColumn() );
    BOOST_CHECK_EQUAL( *firstBlock.row, 0 );
    BOOST_CHECK_EQUAL( *firstBlock.column, 0 );

    auto const& secondBlock = problem.blockTerms().metadata( 1 );
    BOOST_REQUIRE( secondBlock.hasRow() );
    BOOST_REQUIRE( secondBlock.hasColumn() );
    BOOST_CHECK_EQUAL( *secondBlock.row, 0 );
    BOOST_CHECK_EQUAL( *secondBlock.column, 1 );

    auto const& rowMetadata = problem.rowTerms().metadata( 0 );
    BOOST_CHECK( rowMetadata.kind == systems::TermKind::Row );
    BOOST_REQUIRE( rowMetadata.hasRow() );
    BOOST_CHECK_EQUAL( *rowMetadata.row, 1 );
    BOOST_CHECK( !rowMetadata.hasColumn() );

    auto const& constraintMetadata = problem.constraintTerms().metadata( 0 );
    BOOST_CHECK( constraintMetadata.kind == systems::TermKind::Constraint );
    BOOST_CHECK( constraintMetadata.role == systems::TermRole::Constraint );
}

BOOST_AUTO_TEST_CASE( block_and_hdg_diagnostics_expose_row_column_metadata )
{
    using namespace Feel;
    using namespace Feel::vf;
    using namespace boost::hana::literals;
    namespace role = Feel::systems::roles;

    auto mesh = unitSquare();
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = product( Vh, Qh );

    auto u = trial( Vh );
    auto v = test( Vh );
    auto p = trial( Qh );
    auto q = test( Qh );

    auto blockProblem =
        mixedProblem( _space = Xh, _name = "metadata-block-diagnostics" )
            .block( 0_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = inner( grad( u ), grad( v ) ) ) )
            .block( 0_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = -p * div( v ) ) )
            .row( 1_c,
                  integrate( _range = elements( mesh ),
                             _expr = q ) )
            .constraint( []( auto&, auto& ) {} );

    BOOST_CHECK( blockProblem.checkTerms().empty() );
    auto const blockText = blockProblem.describe();
    BOOST_CHECK( contains( blockText, "Problem: metadata-block-diagnostics" ) );
    BOOST_CHECK( contains( blockText, "kind=block" ) );
    BOOST_CHECK( contains( blockText, "row=0" ) );
    BOOST_CHECK( contains( blockText, "col=1" ) );
    BOOST_CHECK( contains( blockText, "kind=row" ) );
    BOOST_CHECK( contains( blockText, "row=1" ) );
    BOOST_CHECK( contains( blockText, "kind=constraint" ) );

    auto H = hdg::poissonSpaces<1>( mesh );
    auto uh = H.flux().trial();
    auto vh = H.flux().test();
    auto phat = H.trace().trial();
    auto mu = H.trace().test();

    auto hdgDiagnosticProblem =
        hdgProblem( _spaces = H, _name = "metadata-hdg-diagnostics" )
            .block( role::flux, role::flux,
                    integrate( _range = elements( mesh ),
                               _expr = inner( uh, vh ) ) )
            .block( role::trace, role::trace,
                    integrate( _range = boundaryfaces( mesh ),
                               _expr = phat * mu ) );

    BOOST_CHECK( hdgDiagnosticProblem.checkTerms().empty() );
    auto const hdgText = hdgDiagnosticProblem.describe();
    BOOST_CHECK( contains( hdgText, "Problem: metadata-hdg-diagnostics" ) );
    BOOST_CHECK( contains( hdgText, "kind=block" ) );
    BOOST_CHECK( contains( hdgText, "row=0" ) );
    BOOST_CHECK( contains( hdgText, "col=0" ) );
    BOOST_CHECK( contains( hdgText, "row=2" ) );
    BOOST_CHECK( contains( hdgText, "col=2" ) );
}

BOOST_AUTO_TEST_CASE( block_and_hdg_assembly_plans_keep_block_and_hdg_strategy_metadata )
{
    using namespace Feel;
    using namespace Feel::vf;
    using namespace boost::hana::literals;
    namespace role = Feel::systems::roles;

    auto mesh = unitSquare();
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = product( Vh, Qh );

    auto u = trial( Vh );
    auto v = test( Vh );
    auto p = trial( Qh );
    auto q = test( Qh );

    auto blockProblem =
        mixedProblem( _space = Xh, _name = "metadata-block-plan" )
            .block( 0_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = inner( grad( u ), grad( v ) ) ) )
            .block( 0_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = -p * div( v ) ) )
            .row( 1_c,
                  integrate( _range = elements( mesh ),
                             _expr = q ) )
            .constraint( []( auto&, auto& ) {} );

    auto const blockPlan = blockProblem.assemblyPlan();
    BOOST_REQUIRE_EQUAL( blockPlan.nTerms(), 4 );
    BOOST_CHECK( hasAssemblyBatch( blockPlan,
                                   []( auto const& batch )
                                   {
                                       return batch.role == systems::TermRole::Matrix &&
                                              batch.hasRow() &&
                                              batch.hasColumn() &&
                                              *batch.row == 0 &&
                                              *batch.column == 1;
                                   } ) );
    BOOST_CHECK( hasAssemblyBatch( blockPlan,
                                   []( auto const& batch )
                                   {
                                       return batch.role == systems::TermRole::Vector &&
                                              batch.hasRow() &&
                                              !batch.hasColumn() &&
                                              *batch.row == 1;
                                   } ) );
    BOOST_CHECK( contains( blockProblem.describeAssemblyPlan(), "row=1" ) );

    auto const blockGraph = blockProblem.assemblyTaskGraph();
    BOOST_REQUIRE_EQUAL( blockGraph.nTerms(), 4 );
    BOOST_CHECK( hasAssemblyTask( blockGraph,
                                  []( auto const& task )
                                  {
                                      return task.info.role == systems::TermRole::Matrix &&
                                             task.info.hasRow() &&
                                             task.info.hasColumn() &&
                                             *task.info.row == 0 &&
                                             *task.info.column == 1;
                                  } ) );
    BOOST_CHECK( contains( blockProblem.describeAssemblyTaskGraph(), "AssemblyTaskGraph: metadata-block-plan" ) );

    auto H = hdg::poissonSpaces<1>( mesh );
    auto uh = H.flux().trial();
    auto vh = H.flux().test();
    auto phat = H.trace().trial();
    auto mu = H.trace().test();

    auto hdgDiagnosticProblem =
        hdgProblem( _spaces = H, _name = "metadata-hdg-plan" )
            .block( role::flux, role::flux,
                    integrate( _range = elements( mesh ),
                               _expr = inner( uh, vh ) ) )
            .block( role::trace, role::trace,
                    integrate( _range = boundaryfaces( mesh ),
                               _expr = phat * mu ) );

    auto const hdgPlan = hdgDiagnosticProblem.assemblyPlan();
    BOOST_REQUIRE_EQUAL( hdgPlan.nTerms(), 2 );
    BOOST_CHECK( !hdgPlan.batches.empty() );
    for ( auto const& batch : hdgPlan.batches )
    {
        BOOST_CHECK( batch.strategy.hdgLocalDenseCandidate );
        BOOST_CHECK( batch.strategy.traceSchurCandidate );
    }
    auto const hdgPlanText = hdgDiagnosticProblem.describeAssemblyPlan();
    BOOST_CHECK( contains( hdgPlanText, "hdgLocalDense=true" ) );
    BOOST_CHECK( contains( hdgPlanText, "traceSchur=true" ) );

    auto const hdgGraph = hdgDiagnosticProblem.assemblyTaskGraph();
    BOOST_REQUIRE_EQUAL( hdgGraph.nTerms(), 2 );
    for ( auto const& task : hdgGraph.tasks )
    {
        BOOST_CHECK( task.info.strategy.hdgLocalDenseCandidate );
        BOOST_CHECK( task.info.strategy.traceSchurCandidate );
    }
    BOOST_CHECK( contains( hdgDiagnosticProblem.describeAssemblyTaskGraph(), "resource=matrix" ) );
}

BOOST_AUTO_TEST_SUITE_END()
