/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

   SPDX-FileCopyrightText: 2026 University of Strasbourg

   SPDX-License-Identifier: LGPL-3.0-or-later
*/

#define BOOST_TEST_MODULE test_functionspace_manager
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/functionspacebuildinstrumentation.hpp>
#include <feel/feeldiscr/functionspacemanager.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelmesh/meshmover.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
FunctionSpaceManagerConfig managerConfig( bool enabled = false,
                                                std::size_t maxEntries = 64,
                                                std::size_t maxEntriesPerMesh = 16 )
{
    return FunctionSpaceManagerConfig{ enabled, maxEntries, maxEntriesPerMesh, true };
}

void resetManager( FunctionSpaceManagerConfig config = managerConfig() )
{
    auto& manager = FunctionSpaceManager::instance();
    manager.clear();
    manager.configure( config );
}

class InstrumentationScope
{
  public:
    InstrumentationScope()
        : M_wasEnabled( FunctionSpaceBuildInstrumentation::enabled() )
    {
        FunctionSpaceBuildInstrumentation::reset();
        FunctionSpaceBuildInstrumentation::setEnabled( true );
    }

    ~InstrumentationScope()
    {
        FunctionSpaceBuildInstrumentation::reset();
        FunctionSpaceBuildInstrumentation::setEnabled( M_wasEnabled );
    }

  private:
    bool M_wasEnabled;
};

/**
 * @brief Request the same managed space twice and verify that it is reused.
 * @tparam Factory Callable returning a managed function-space pointer.
 * @param factory Explicit-reuse factory request.
 * @return The first managed space pointer.
 */
template <typename Factory>
auto checkExplicitReuse( Factory&& factory )
{
    auto first = factory();
    auto second = factory();
    BOOST_TEST( first.get() == second.get() );
    return first;
}

/**
 * @brief Verify the pointer-reuse result of an automatic-policy factory.
 * @tparam Factory Callable returning an automatic-policy function-space pointer.
 * @param factory Automatic-policy factory request.
 * @param shouldReuse Whether the manager configuration enables reuse.
 * @return The first function-space pointer.
 */
template <typename Factory>
auto checkAutomaticPolicy( Factory&& factory, bool shouldReuse )
{
    auto first = factory();
    auto second = factory();
    BOOST_TEST( ( first.get() == second.get() ) == shouldReuse );
    return first;
}

/**
 * @brief Verify that two managed factories retain distinct concrete spaces.
 * @tparam LeftFactory Callable returning the first managed function-space pointer.
 * @tparam RightFactory Callable returning the second managed function-space pointer.
 * @param leftFactory Explicit-reuse request for the first concrete space type.
 * @param rightFactory Explicit-reuse request for the second concrete space type.
 */
template <typename LeftFactory, typename RightFactory>
void checkDistinctManagedTypes( LeftFactory&& leftFactory, RightFactory&& rightFactory )
{
    auto left = leftFactory();
    auto right = rightFactory();
    BOOST_TEST( static_cast<void*>( left.get() ) != static_cast<void*>( right.get() ) );
}
} // namespace

BOOST_AUTO_TEST_CASE( explicit_reuse_normalizes_options_and_reports_stats )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );
    InstrumentationScope instrumentation;

    auto XhDefault = checkExplicitReuse( [&mesh]()
    {
        return Pch<2>( _mesh=mesh,
                       _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );
    auto XhVertices = Pch<2>(
        _mesh=mesh,
        _extended_doftable=DofTableExtendedType::VERTICES,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto XhNone = Pch<2>(
        _mesh=mesh,
        _extended_doftable=DofTableExtendedType::NONE,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto Dh = checkExplicitReuse( [&mesh]()
    {
        return Pdh<2>( _mesh=mesh,
                       _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );

    BOOST_TEST( XhDefault.get() == XhVertices.get() );
    BOOST_TEST( XhDefault.get() != XhNone.get() );

    auto const counts = FunctionSpaceBuildInstrumentation::counts();
    BOOST_TEST( counts.functionSpaceConstructions == 3 );
    BOOST_TEST( counts.dofTableBuilds == 3 );

    auto const stats = FunctionSpaceManager::instance().stats( mesh );
    BOOST_TEST( stats.lookups == 6 );
    BOOST_TEST( stats.hits == 3 );
    BOOST_TEST( stats.misses == 3 );
    BOOST_TEST( stats.builds == 3 );
    BOOST_TEST( stats.retainedEntries == 3 );
    BOOST_TEST( stats.retainedDofs ==
                XhDefault->nDof() + XhNone->nDof() + Dh->nDof() );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( concrete_function_space_types_are_isolated )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );

    checkDistinctManagedTypes(
        [&mesh]()
        {
            return Pch<2>( _mesh=mesh,
                           _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
        },
        [&mesh]()
        {
            return Pdh<2>( _mesh=mesh,
                           _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
        } );

    auto const stats = FunctionSpaceManager::instance().stats( mesh );
    BOOST_TEST( stats.retainedEntries == 2 );
    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( vector_factories_reuse_and_preserve_range_semantics )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );

    auto Vh = checkExplicitReuse( [&mesh]()
    {
        return Pchv<1>( _mesh=mesh,
                         _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );
    auto VhVertices = Pchv<1>(
        _mesh=mesh,
        _extended_doftable=DofTableExtendedType::VERTICES,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto VhNone = Pchv<1>(
        _mesh=mesh,
        _extended_doftable=DofTableExtendedType::NONE,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto Dvh = checkExplicitReuse( [&mesh]()
    {
        return Pdhv<1>( _mesh=mesh,
                         _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );
    auto Xh = Pch<1>( _mesh=mesh,
                      _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto Dh = Pdh<1>( _mesh=mesh,
                      _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto VhOrder = Pchv<2>( _mesh=mesh,
                            _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto VhTag = Pchv<1,PointSetFekete,double,1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto VhPointSet = Pchv<2,PointSetEquiSpaced>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto VhFloat = checkExplicitReuse( [&mesh]()
    {
        return Pchv<1,PointSetFekete,float>(
            _mesh=mesh,
            _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );
    auto DvhOrder = Pdhv<2>( _mesh=mesh,
                             _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto DvhTag = Pdhv<1,PointSetFekete,1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto DvhPointSet = Pdhv<2,PointSetEquiSpaced>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );

    BOOST_TEST( Vh.get() == VhVertices.get() );
    BOOST_TEST( Vh.get() != VhNone.get() );
    BOOST_TEST( static_cast<void*>( Vh.get() ) != static_cast<void*>( Dvh.get() ) );
    BOOST_TEST( static_cast<void*>( Vh.get() ) != static_cast<void*>( Xh.get() ) );
    BOOST_TEST( static_cast<void*>( Vh.get() ) != static_cast<void*>( Dh.get() ) );
    BOOST_TEST( static_cast<void*>( Dvh.get() ) != static_cast<void*>( Xh.get() ) );
    BOOST_TEST( static_cast<void*>( Dvh.get() ) != static_cast<void*>( Dh.get() ) );
    BOOST_TEST( static_cast<void*>( Xh.get() ) != static_cast<void*>( Dh.get() ) );
    BOOST_TEST( static_cast<void*>( Vh.get() ) != static_cast<void*>( VhOrder.get() ) );
    BOOST_TEST( static_cast<void*>( Vh.get() ) != static_cast<void*>( VhTag.get() ) );
    BOOST_TEST( static_cast<void*>( Vh.get() ) != static_cast<void*>( VhPointSet.get() ) );
    BOOST_TEST( static_cast<void*>( Vh.get() ) != static_cast<void*>( VhFloat.get() ) );
    BOOST_TEST( static_cast<void*>( Dvh.get() ) != static_cast<void*>( DvhOrder.get() ) );
    BOOST_TEST( static_cast<void*>( Dvh.get() ) != static_cast<void*>( DvhTag.get() ) );
    BOOST_TEST( static_cast<void*>( Dvh.get() ) != static_cast<void*>( DvhPointSet.get() ) );

    auto continuousRange1 = Pchv<1>( mesh, elements( mesh ) );
    auto continuousRange2 = Pchv<1>( mesh, elements( mesh ) );
    auto discontinuousRange1 = Pdhv<1>( mesh, elements( mesh ) );
    auto discontinuousRange2 = Pdhv<1>( mesh, elements( mesh ) );
    BOOST_TEST( continuousRange1.get() != continuousRange2.get() );
    BOOST_TEST( discontinuousRange1.get() != discontinuousRange2.get() );

    auto const stats = FunctionSpaceManager::instance().stats( mesh );
    BOOST_TEST( stats.retainedEntries == 12 );
    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( vector_automatic_policy_follows_manager_configuration )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );

    checkAutomaticPolicy( [&mesh]() { return Pchv<1>( mesh ); }, false );
    checkAutomaticPolicy( [&mesh]() { return Pdhv<1>( mesh ); }, false );

    resetManager( managerConfig( true ) );
    checkAutomaticPolicy( [&mesh]() { return Pchv<1>( mesh ); }, true );
    checkAutomaticPolicy( [&mesh]() { return Pdhv<1>( mesh ); }, true );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( vector_named_policies_replace_and_bypass_matching_entries )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );

    auto Vh = Pchv<1>( _mesh=mesh,
                       _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto rebuiltVh = Pchv<1>( _mesh=mesh,
                              _fspace_reuse_policy=FunctionSpaceReusePolicy::rebuild );
    auto reusedVh = Pchv<1>( _mesh=mesh,
                             _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto bypassedVh = Pchv<1>( _mesh=mesh,
                               _fspace_reuse_policy=FunctionSpaceReusePolicy::bypass );
    BOOST_TEST( rebuiltVh.get() != Vh.get() );
    BOOST_TEST( reusedVh.get() == rebuiltVh.get() );
    BOOST_TEST( bypassedVh.get() != rebuiltVh.get() );

    auto Dvh = Pdhv<1>( _mesh=mesh,
                        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto rebuiltDvh = Pdhv<1>( _mesh=mesh,
                               _fspace_reuse_policy=FunctionSpaceReusePolicy::rebuild );
    auto reusedDvh = Pdhv<1>( _mesh=mesh,
                              _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto bypassedDvh = Pdhv<1>( _mesh=mesh,
                                _fspace_reuse_policy=FunctionSpaceReusePolicy::bypass );
    BOOST_TEST( rebuiltDvh.get() != Dvh.get() );
    BOOST_TEST( reusedDvh.get() == rebuiltDvh.get() );
    BOOST_TEST( bypassedDvh.get() != rebuiltDvh.get() );
    BOOST_TEST( FunctionSpaceManager::instance().stats( mesh ).retainedEntries == 2 );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( tensor_factories_reuse_and_keep_concrete_types_isolated )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );

    auto Mh = checkExplicitReuse( [&mesh]()
    {
        return Pchm<1>( _mesh=mesh,
                         _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );
    auto Mhs = checkExplicitReuse( [&mesh]()
    {
        return Pchms<1>( _mesh=mesh,
                          _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );
    auto Dmh = checkExplicitReuse( [&mesh]()
    {
        return Pdhm<1>( _mesh=mesh,
                         _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );
    auto Dmhs = checkExplicitReuse( [&mesh]()
    {
        return Pdhms<1>( _mesh=mesh,
                          _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );
    auto MhVertices = Pchm<1>(
        _mesh=mesh,
        _extended_doftable=DofTableExtendedType::VERTICES,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto MhNone = Pchm<1>(
        _mesh=mesh,
        _extended_doftable=DofTableExtendedType::NONE,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto MhTag = Pchm<1,double,PointSetFekete,1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto MhPointSet = Pchm<2,double,PointSetEquiSpaced>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto MhFloat = checkExplicitReuse( [&mesh]()
    {
        return Pchm<1,float>( _mesh=mesh,
                               _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );
    auto DmhTag = Pdhm<1,double,PointSetFekete,1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );

    BOOST_TEST( Mh.get() == MhVertices.get() );
    BOOST_TEST( static_cast<void*>( Mh.get() ) != static_cast<void*>( Mhs.get() ) );
    BOOST_TEST( static_cast<void*>( Mh.get() ) != static_cast<void*>( Dmh.get() ) );
    BOOST_TEST( static_cast<void*>( Mh.get() ) != static_cast<void*>( Dmhs.get() ) );
    BOOST_TEST( static_cast<void*>( Mhs.get() ) != static_cast<void*>( Dmh.get() ) );
    BOOST_TEST( static_cast<void*>( Mhs.get() ) != static_cast<void*>( Dmhs.get() ) );
    BOOST_TEST( static_cast<void*>( Dmh.get() ) != static_cast<void*>( Dmhs.get() ) );
    BOOST_TEST( static_cast<void*>( Mh.get() ) != static_cast<void*>( MhNone.get() ) );
    BOOST_TEST( static_cast<void*>( Mh.get() ) != static_cast<void*>( MhTag.get() ) );
    BOOST_TEST( static_cast<void*>( Mh.get() ) != static_cast<void*>( MhPointSet.get() ) );
    BOOST_TEST( static_cast<void*>( Mh.get() ) != static_cast<void*>( MhFloat.get() ) );
    BOOST_TEST( static_cast<void*>( Dmh.get() ) != static_cast<void*>( DmhTag.get() ) );
    BOOST_TEST( FunctionSpaceManager::instance().stats( mesh ).retainedEntries == 9 );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( tensor_automatic_policy_reuses_by_default )
{
    resetManager( managerConfig( true ) );
    auto mesh = unitSquare( 0.25 );

    checkAutomaticPolicy( [&mesh]() { return Pchm<1>( mesh ); }, true );
    checkAutomaticPolicy( [&mesh]() { return Pchms<1>( mesh ); }, true );
    checkAutomaticPolicy( [&mesh]() { return Pdhm<1>( mesh ); }, true );
    checkAutomaticPolicy( [&mesh]() { return Pdhms<1>( mesh ); }, true );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( tensor_rebuild_and_bypass_affect_only_the_matching_entry )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );

    auto Mh = Pchm<1>( _mesh=mesh,
                       _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto Dmh = Pdhm<1>( _mesh=mesh,
                        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto rebuiltMh = Pchm<1>( _mesh=mesh,
                              _fspace_reuse_policy=FunctionSpaceReusePolicy::rebuild );
    auto reusedDmh = Pdhm<1>( _mesh=mesh,
                              _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( rebuiltMh.get() != Mh.get() );
    BOOST_TEST( reusedDmh.get() == Dmh.get() );
    BOOST_TEST( FunctionSpaceManager::instance().stats( mesh ).retainedEntries == 2 );

    auto const statsBeforeBypass = FunctionSpaceManager::instance().stats( mesh );
    auto bypassedMh = Pchm<1>( _mesh=mesh,
                               _fspace_reuse_policy=FunctionSpaceReusePolicy::bypass );
    auto const statsAfterBypass = FunctionSpaceManager::instance().stats( mesh );
    BOOST_TEST( bypassedMh.get() != rebuiltMh.get() );
    BOOST_TEST( statsAfterBypass.lookups == statsBeforeBypass.lookups );
    BOOST_TEST( statsAfterBypass.hits == statsBeforeBypass.hits );
    BOOST_TEST( statsAfterBypass.misses == statsBeforeBypass.misses );
    BOOST_TEST( statsAfterBypass.builds == statsBeforeBypass.builds );
    BOOST_TEST( statsAfterBypass.rebuilds == statsBeforeBypass.rebuilds );
    BOOST_TEST( statsAfterBypass.retainedEntries == statsBeforeBypass.retainedEntries );
    BOOST_TEST( statsAfterBypass.retainedDofs == statsBeforeBypass.retainedDofs );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( policies_preserve_default_and_direct_new_semantics )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );
    InstrumentationScope instrumentation;

    auto automatic1 = Pch<1>( mesh );
    auto automatic2 = Pch<1>( mesh );
    BOOST_TEST( automatic1.get() != automatic2.get() );

    auto reused1 = checkExplicitReuse( [&mesh]()
    {
        return Pch<1>( _mesh=mesh,
                       _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    } );

    auto rebuilt = Pch<1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::rebuild );
    BOOST_TEST( rebuilt.get() != reused1.get() );
    auto reusedAfterRebuild = Pch<1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( reusedAfterRebuild.get() == rebuilt.get() );

    auto bypassed = Pch<1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::bypass );
    BOOST_TEST( bypassed.get() != rebuilt.get() );

    using space_type = Pch_type<typename decltype( mesh )::element_type,1>;
    auto direct = space_type::New( _mesh=mesh );
    BOOST_TEST( direct.get() != rebuilt.get() );

    auto const stats = FunctionSpaceManager::instance().stats( mesh );
    BOOST_TEST( stats.lookups == 4 );
    BOOST_TEST( stats.hits == 2 );
    BOOST_TEST( stats.misses == 2 );
    BOOST_TEST( stats.builds == 2 );
    BOOST_TEST( stats.rebuilds == 1 );
    BOOST_TEST( stats.retainedEntries == 1 );

    resetManager( managerConfig( true ) );
    auto enabled1 = Pch<1>( mesh );
    auto enabled2 = Pch<1>( mesh );
    BOOST_TEST( enabled1.get() == enabled2.get() );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( bounded_lru_and_per_mesh_clear )
{
    resetManager( managerConfig( false, 2, 1 ) );
    auto mesh1 = unitSquare( 0.25 );
    auto mesh2 = unitSquare( 0.20 );

    auto Xh1 = Pch<1>(
        _mesh=mesh1,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto Xh2 = Pch<2>(
        _mesh=mesh1,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( FunctionSpaceManager::instance().stats( mesh1 ).retainedEntries == 1 );

    auto Xh1Again = Pch<1>(
        _mesh=mesh1,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( Xh1Again.get() != Xh1.get() );

    auto Yh = Pch<1>(
        _mesh=mesh2,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( FunctionSpaceManager::instance().stats().retainedEntries == 2 );

    FunctionSpaceManager::instance().clear( mesh1 );
    BOOST_TEST( FunctionSpaceManager::instance().stats( mesh1 ).retainedEntries == 0 );
    BOOST_TEST( FunctionSpaceManager::instance().stats( mesh2 ).retainedEntries == 1 );

    auto const mesh1Stats = FunctionSpaceManager::instance().stats( mesh1 );
    BOOST_TEST( mesh1Stats.evictions >= 2 );
    BOOST_TEST( Xh1->nDof() > 0 );
    BOOST_TEST( Xh2->nDof() > 0 );
    BOOST_TEST( Yh->nDof() > 0 );

    FunctionSpaceManager::instance().clear();

    resetManager( managerConfig( false, 2, 2 ) );
    auto global1 = Pch<1>(
        _mesh=mesh1,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto global2 = Pch<1>(
        _mesh=mesh2,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto global1Hit = Pch<1>(
        _mesh=mesh1,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( global1Hit.get() == global1.get() );
    auto global3 = Pch<2>(
        _mesh=mesh1,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto global2AfterEviction = Pch<1>(
        _mesh=mesh2,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( global2AfterEviction.get() != global2.get() );
    BOOST_TEST( FunctionSpaceManager::instance().stats().retainedEntries == 2 );
    BOOST_TEST( global3->nDof() > 0 );

    FunctionSpaceManager::instance().clear();
}

BOOST_AUTO_TEST_CASE( structural_changes_invalidate_but_ranges_stay_uncached )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );
    InstrumentationScope instrumentation;

    auto Xh1 = Pch<1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto const revisionBefore = mesh->functionSpaceStructuralRevision();
    mesh->meshModified();
    BOOST_TEST( mesh->functionSpaceStructuralRevision() > revisionBefore );

    auto Xh2 = Pch<1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( Xh2.get() != Xh1.get() );
    auto const stats = FunctionSpaceManager::instance().stats( mesh );
    BOOST_TEST( stats.invalidations == 1 );

    auto range1 = Pch<1>( mesh, elements( mesh ) );
    auto range2 = Pch<1>( mesh, elements( mesh ) );
    BOOST_TEST( range1.get() != range2.get() );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( failed_factory_does_not_install_an_entry )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );
    using mesh_type = typename decltype( mesh )::element_type;
    using space_type = Pch_type<mesh_type,1>;
    FunctionSpaceManagerOptions const options{
        DofTableExtendedType::VERTICES, MESH_RENUMBER | MESH_CHECK };

    BOOST_CHECK_THROW(
        getOrCreateFunctionSpace<space_type>(
            mesh, options, FunctionSpaceReusePolicy::reuse,
            []() -> std::shared_ptr<space_type>
            { throw std::runtime_error( "expected test failure" ); } ),
        std::runtime_error );
    BOOST_TEST( FunctionSpaceManager::instance().stats( mesh ).retainedEntries == 0 );

    auto Xh = Pch<1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( Xh->nDof() > 0 );
    BOOST_TEST( FunctionSpaceManager::instance().stats( mesh ).retainedEntries == 1 );
    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( geometry_changes_refresh_live_reused_spaces )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );
    auto Xh = Pch<1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    auto const structuralRevision = mesh->functionSpaceStructuralRevision();
    auto const geometryRevision = mesh->geometryRevision();
    auto const pointBefore = Xh->dof()->dofPoint( 0 ).template get<0>();

    auto Vh = Pchv<1>( mesh );
    auto displacement = Vh->element();
    displacement.setConstant( 1.0 );
    meshMove( mesh, displacement );

    auto const pointAfter = Xh->dof()->dofPoint( 0 ).template get<0>();
    BOOST_TEST( mesh->geometryRevision() == geometryRevision + 1 );
    BOOST_TEST( mesh->functionSpaceStructuralRevision() == structuralRevision );
    BOOST_TEST( pointAfter[0] == pointBefore[0] + 1.0, boost::test_tools::tolerance( 1e-12 ) );
    BOOST_TEST( pointAfter[1] == pointBefore[1] + 1.0, boost::test_tools::tolerance( 1e-12 ) );

    auto XhAgain = Pch<1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( XhAgain.get() == Xh.get() );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( mpi_local_miss_forces_collective_rebuild )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );
    InstrumentationScope instrumentation;

    auto Xh1 = Pch<2>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    if ( mesh->worldComm().globalRank() == 0 )
        FunctionSpaceManager::instance().clear( mesh );
    mesh->worldComm().globalComm().barrier();

    auto Xh2 = Pch<2>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( Xh2.get() != Xh1.get() );

    auto const counts = FunctionSpaceBuildInstrumentation::counts();
    BOOST_TEST( counts.functionSpaceConstructions == 2 );
    BOOST_TEST( counts.dofTableBuilds == 2 );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( mpi_vector_local_miss_forces_collective_rebuild )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );
    InstrumentationScope instrumentation;

    auto Xh1 = Pchv<2>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    if ( mesh->worldComm().globalRank() == 0 )
        FunctionSpaceManager::instance().clear( mesh );
    mesh->worldComm().globalComm().barrier();

    auto Xh2 = Pchv<2>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( Xh2.get() != Xh1.get() );

    auto const counts = FunctionSpaceBuildInstrumentation::counts();
    BOOST_TEST( counts.functionSpaceConstructions == 2 );
    BOOST_TEST( counts.dofTableBuilds == 2 );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( mpi_tensor_local_miss_forces_collective_rebuild )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );
    InstrumentationScope instrumentation;

    auto Xh1 = Pchm<1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    if ( mesh->worldComm().globalRank() == 0 )
        FunctionSpaceManager::instance().clear( mesh );
    mesh->worldComm().globalComm().barrier();

    auto Xh2 = Pchm<1>(
        _mesh=mesh,
        _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
    BOOST_TEST( Xh2.get() != Xh1.get() );

    auto const counts = FunctionSpaceBuildInstrumentation::counts();
    BOOST_TEST( counts.functionSpaceConstructions == 2 );
    BOOST_TEST( counts.dofTableBuilds == 2 );

    FunctionSpaceManager::instance().clear( mesh );
}

BOOST_AUTO_TEST_CASE( mpi_request_signature_mismatch_is_rejected )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );
    if ( mesh->worldComm().globalSize() == 1 )
        return;

    auto const dte = mesh->worldComm().globalRank() % 2 == 0
                         ? DofTableExtendedType::VERTICES
                         : DofTableExtendedType::NONE;
    BOOST_CHECK_THROW(
        Pch<1>(
            _mesh=mesh,
            _extended_doftable=dte,
            _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse ),
        std::runtime_error );
    BOOST_TEST( FunctionSpaceManager::instance().stats( mesh ).retainedEntries == 0 );
}

BOOST_AUTO_TEST_CASE( mpi_vector_request_signature_mismatch_is_rejected )
{
    resetManager();
    auto mesh = unitSquare( 0.25 );
    if ( mesh->worldComm().globalSize() == 1 )
        return;

    auto const dte = mesh->worldComm().globalRank() % 2 == 0
                         ? DofTableExtendedType::VERTICES
                         : DofTableExtendedType::NONE;
    BOOST_CHECK_THROW(
        Pchv<1>(
            _mesh=mesh,
            _extended_doftable=dte,
            _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse ),
        std::runtime_error );
    BOOST_TEST( FunctionSpaceManager::instance().stats( mesh ).retainedEntries == 0 );
}
