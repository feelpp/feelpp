// SPDX-License-Identifier: LGPL-3.0-or-later
/** @file test_exporter_transient.cpp
 * @brief EnSight temporal metadata, immutable snapshots and lifecycle regression tests.
 *
 * Run serially and with MPI, with unmerged output and merge/pack settings.
 * Numeric reader integration is exercised by the exporter benchmark readback tool.
 */
#define BOOST_TEST_MODULE test_exporter_transient
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/exportfieldset.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelmesh/meshmover.hpp>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** @brief Detect accidental temporal metadata in the shared field-storage model. */
template <typename T>
concept HasTime = requires( T value ) { value.time(); };

//! \brief Detect temporal or publication policy leaking into field storage.
template <typename T>
concept HasStepPolicy =
    requires( T value ) { value.index(); } || requires( T value ) { value.isIgnored(); } ||
    requires( T value ) { value.isNew(); } || requires( T value ) { value.isOnDisk(); } ||
    requires( T value ) { value.state(); } || requires( T value ) { value.setState( 0 ); } ||
    requires( T value ) { value.load(); };

/** @brief Read a small metadata file after its writer has closed it. */
std::string readMetadata( fs::path const &path )
{
    std::ifstream input( path.string() );
    BOOST_REQUIRE_MESSAGE( input, "cannot read " << path );
    return { std::istreambuf_iterator<char>( input ), std::istreambuf_iterator<char>() };
}

//! \brief Restore process-wide manager configuration after a test, including failures.
class ManagerConfigScope
{
  public:
    //! \brief Select automatic reuse without changing the caller's permanent settings.
    explicit ManagerConfigScope( bool enabled )
        : M_previous( FunctionSpaceManager::instance().config() )
    {
        auto config = M_previous;
        config.enabled = enabled;
        config.maxEntries = 64;
        config.maxEntriesPerMesh = 16;
        FunctionSpaceManager::instance().configure( config );
    }

    //! \brief Restore the settings used by the surrounding regression suite.
    ~ManagerConfigScope() { FunctionSpaceManager::instance().configure( M_previous ); }

  private:
    FunctionSpaceManagerConfig M_previous; //!< Configuration to restore on scope exit.
};

//! \brief Check warm/cold managed export spaces without sharing mutable field buffers.
//! \tparam MeshType Linear or quadratic geometric mesh type.
//! \param mesh Mesh shared by application spaces, exporters and transient steps.
template <typename MeshType>
void checkManagedConversionSpaces( std::shared_ptr<MeshType> const &mesh )
{
    auto &manager = FunctionSpaceManager::instance();
    for ( bool enabled : { true, false } )
    {
        ManagerConfigScope config( enabled );
        manager.clear( mesh );
        auto nodalSpace = Pch<MeshType::nOrder>( mesh );
        auto elementSpace = Pdh<0>( mesh );
        auto const before = manager.stats( mesh );
        auto e = exporter( _mesh = mesh, _name = "managed_conversion", _geo = "static" );
        auto other = exporter( _mesh = mesh, _name = "other_conversion", _geo = "static" );

        // Expressions force conversion-space lookup instead of adopting an input FE space.
        e->add( "reference", cst( 3. ) );
        e->add( "cell_reference", cst( 7. ), "element" );
        other->add( "reference", cst( 4. ) );
        other->add( "cell_reference", cst( 8. ), "element" );
        auto first = e->step( 0. );
        first->add( "signal", cst( 5. ) );
        first->add( "cell_signal", cst( 9. ), elements( mesh ), "element" );
        auto second = e->step( 1. );
        second->add( "signal", cst( 6. ) );
        second->add( "cell_signal", cst( 10. ), elements( mesh ), "element" );

        for ( auto const &fields : { e->staticFields(), other->staticFields(),
                                     first->fieldSet(), second->fieldSet() } )
        {
            auto const &nodal = fields->beginNodal()->second.second[0][0];
            auto const &element = fields->beginElement()->second.second[0][0];
            BOOST_CHECK_EQUAL( nodal->functionSpace() == nodalSpace, enabled );
            BOOST_CHECK_EQUAL( element->functionSpace() == elementSpace, enabled );
        }
        BOOST_CHECK( first->fieldSet()->spaceCache() == second->fieldSet()->spaceCache() );
        BOOST_CHECK( first->fieldSet()->spaceCache()->M_p1 ==
                     second->fieldSet()->beginNodal()->second.second[0][0]->functionSpace() );
        BOOST_CHECK( first->fieldSet()->spaceCache()->M_p0 ==
                     second->fieldSet()->beginElement()->second.second[0][0]->functionSpace() );
        BOOST_CHECK( e->staticFields()->beginNodal()->second.second[0][0] !=
                     first->fieldSet()->beginNodal()->second.second[0][0] );
        BOOST_CHECK_CLOSE( e->staticFields()->beginNodal()->second.second[0][0]->min(), 3., 1e-8 );
        BOOST_CHECK_CLOSE( first->fieldSet()->beginNodal()->second.second[0][0]->min(), 5., 1e-8 );
        BOOST_CHECK_CLOSE( second->fieldSet()->beginNodal()->second.second[0][0]->min(), 6., 1e-8 );
        BOOST_CHECK_CLOSE( e->staticFields()->beginElement()->second.second[0][0]->min(),
                           7., 1e-8 );
        auto const after = manager.stats( mesh );
        BOOST_CHECK_EQUAL( after.builds, before.builds );
        BOOST_CHECK_EQUAL( after.hits - before.hits, enabled ? 6 : 0 );

        // A cold exporter registers spaces which subsequent application factories reuse.
        manager.clear( mesh );
        auto cold = exporter( _mesh = mesh, _name = "cold_conversion", _geo = "static" );
        cold->add( "reference", cst( 1. ) );
        cold->add( "cell_reference", cst( 2. ), "element" );
        BOOST_CHECK_EQUAL( cold->staticFields()->spaceCache()->M_p1 ==
                               Pch<MeshType::nOrder>( mesh ), enabled );
        BOOST_CHECK_EQUAL( cold->staticFields()->spaceCache()->M_p0 == Pdh<0>( mesh ), enabled );
        manager.clear( mesh );
    }
}
} // namespace

//! \brief Reuse P1/P2 nodal and P0 element spaces across static and transient fields.
//! Both enabled and disabled automatic policies run on every MPI rank.
BOOST_AUTO_TEST_CASE( managed_conversion_spaces )
{
    checkManagedConversionSpaces( unitSquare( 0.25 ) );
#if FEELPP_MESH_MAX_ORDER >= 2
    auto quadratic = createGMSHMesh(
        _mesh = new Mesh<Simplex<2, 2>>,
        _desc = domain( _name = "managed_quadratic", _shape = "hypercube", _dim = 2,
                        _order = 2, _h = 0.25 ) );
    checkManagedConversionSpaces( quadratic );
#endif
}

//! \brief Separate field residency from scheduling, publication and legacy archives.
//! Ignored samples must bypass every registration path before requiring a mesh.
BOOST_AUTO_TEST_CASE( fieldset_lifecycle )
{
    auto mesh = unitSquare( 0.25 );
    using mesh_type = typename decltype( mesh )::element_type;
    using field_set_type = Feel::detail::ExportFieldSet<mesh_type>;
    using time_set_type = Feel::detail::TimeSet<mesh_type>;
    static_assert( !HasTime<field_set_type> && !HasStepPolicy<field_set_type> );
    static_assert( HasTime<typename time_set_type::Step> );
    static_assert( HasStepPolicy<typename time_set_type::Step> );

    field_set_type fields;
    BOOST_CHECK( !fields.hasData() && !fields.isInMemory() && !fields.hasMesh() );
    BOOST_CHECK_EQUAL( fields.revision(), 0 );
    fields.add( "constant", 2.5 );
    auto revision = fields.revision();
    BOOST_CHECK( fields.hasData() && fields.isInMemory() );
    fields.cleanup();
    BOOST_CHECK( fields.hasData() && !fields.isInMemory() );
    BOOST_CHECK_EQUAL( fields.revision(), revision );
    BOOST_CHECK_EQUAL( fields.scalar( "constant" ), 2.5 );

    time_set_type sequence( "fieldset_lifecycle" );
    sequence.setMesh( mesh );
    auto active = sequence.step( 0., 2 );
    auto ignored = sequence.step( 1., 2 );
    auto last = sequence.step( 2., 2 );
    auto scalar = Pch<1>( mesh )->element( cst( 4. ) );
    auto const &cache = ignored->fieldSet()->spaceCache();
    BOOST_CHECK( active->fieldSet() != ignored->fieldSet() );
    BOOST_CHECK( active->fieldSet()->spaceCache() == cache );
    BOOST_CHECK( !cache->M_p0 && !cache->M_p1 );

    ignored->add( "scalar", 1. );
    ignored->addComplex( "complex", complex_type( 1., 2. ) );
    ignored->add( "field", scalar );
    ignored->add( { "named_field" }, scalar );
    ignored->add( "separate_name", "separate_file", scalar );
    ignored->add( "both", scalar, std::set<std::string>{ "nodal", "element" } );
    ignored->add( "expression", Px() );
    ignored->add( "range", cst( 2. ), elements( mesh ), "element" );
    ignored->addRegions();
    ignored->addRegions( "region", "region_file" );
    BOOST_CHECK( ignored->isIgnored() && !ignored->hasMesh() );
    BOOST_CHECK( !ignored->hasData() && !ignored->isInMemory() && !ignored->isOnDisk() );
    BOOST_CHECK( ignored->fieldSet()->names().empty() );
    BOOST_CHECK_EQUAL( ignored->fieldSet()->revision(), 0 );
    BOOST_CHECK( !cache->M_p0 && !cache->M_p1 );

    active->add( "field", scalar );
    BOOST_CHECK( cache->M_p1 );
    active->setState( STEP_ON_DISK );
    BOOST_CHECK( active->isOnDisk() );
    revision = active->fieldSet()->revision();
    active->cleanup();
    BOOST_CHECK( active->isOnDisk() && !active->isInMemory() && active->hasData() );
    BOOST_CHECK_EQUAL( active->fieldSet()->revision(), revision );
    BOOST_CHECK_EQUAL( active->fieldName( "field", true ), "field" );
    BOOST_CHECK( active->nodal( "field" ).second.empty() );
    active->load(); // Legacy metadata flag only: no payload is restored.
    BOOST_CHECK( active->isInMemory() && active->isOnDisk() );
    BOOST_CHECK( active->nodal( "field" ).second.empty() );
    active->fieldSet()->add( "direct", 7. );
    BOOST_CHECK( !active->isOnDisk() );
    BOOST_CHECK_EQUAL( active->state() & STEP_ON_DISK, 0 );
    active->setState( STEP_ON_DISK );
    active->setMesh( mesh );
    BOOST_CHECK( !active->isOnDisk() );
    active->setState( STEP_ON_DISK );
    active->fieldSet()->materialize( fields, true );
    BOOST_CHECK( !active->isOnDisk() );
    BOOST_CHECK_EQUAL( active->scalar( "constant" ), 2.5 );
    active->fieldSet()->materialize( fields, false );
    BOOST_CHECK( !active->fieldSet()->names().count( "constant" ) );
    active->setState( STEP_ON_DISK );
    active->cleanup();
    BOOST_CHECK_EQUAL( sequence.stepsToWriteOnDisk().size(), 1 );
    BOOST_CHECK( *sequence.stepsToWriteOnDisk().begin() == last );

    // The wire format still contains exactly the legacy state bitset.
    std::stringstream metadata;
    {
        boost::archive::text_oarchive archive( metadata );
        archive << sequence;
    }
    time_set_type restored( "restored" );
    {
        boost::archive::text_iarchive archive( metadata );
        archive >> restored;
    }
    BOOST_CHECK_EQUAL( restored.numberOfSteps(), 3 );
    BOOST_CHECK_EQUAL( restored.numberOfActiveSteps(), 2 );
    for ( double time : { 0., 1., 2. } )
    {
        auto expected = sequence.step( time, 2 );
        auto actual = restored.step( time, 2 );
        BOOST_CHECK_EQUAL( actual->state(), expected->state() );
        BOOST_CHECK_EQUAL( actual->index(), expected->index() );
        BOOST_CHECK_EQUAL( actual->activeIndex(), expected->activeIndex() );
    }
    BOOST_CHECK( restored.step( 0. )->isOnDisk() );
    BOOST_CHECK( !restored.step( 0. )->isInMemory() );
    BOOST_CHECK( restored.step( 1. )->isIgnored() );
}

//! \brief Keep the minimal I/O correctness policy without exposing tuning controls.
//! All ranks reject invalid/divergent policies before I/O; successful output
//! freezes the policy. Global mesh numbering must not silently wrap.
BOOST_AUTO_TEST_CASE( io_policy_contract )
{
    auto mesh = unitSquare( 0.25 );
    using mesh_type = typename decltype( mesh )::element_type;
    using numbering_type = Feel::detail::MeshContiguousNumberingMapping<mesh_type, float>;
    auto const largest = std::numeric_limits<typename mesh_type::index_type>::max();
    BOOST_CHECK_EQUAL( numbering_type::checkedNumberingSum( largest - 1, 1 ), largest );
    BOOST_CHECK_THROW( numbering_type::checkedNumberingSum( largest, 1 ), std::overflow_error );

    auto invalid = exporter( _mesh = mesh, _name = "invalid_io_policy", _geo = "static" );
    invalid->setIOPolicy( static_cast<ExporterIOPolicy>( 99 ) );
    BOOST_CHECK_THROW( invalid->save(), std::invalid_argument );

    if ( Environment::worldComm().size() > 1 )
    {
        auto divergent = exporter( _mesh = mesh, _name = "divergent_io_policy", _geo = "static" );
        divergent->setIOPolicy( Environment::isMasterRank() ? ExporterIOPolicy::Automatic
                                                          : ExporterIOPolicy::Root );
        BOOST_CHECK_THROW( divergent->save(), std::invalid_argument );
    }

    auto pressure = Pch<1>( mesh )->element( cst( 1. ) );
    auto e = exporter( _mesh = mesh, _name = "io_policy_contract", _geo = "static" );
    e->setIOPolicy( ExporterIOPolicy::Root );
    e->step( 0. )->add( "pressure", pressure );
    e->save();
    BOOST_CHECK_THROW( e->setIOPolicy( ExporterIOPolicy::Automatic ), std::logic_error );
}

/** @brief Static snapshots survive source mutation and dynamic-step cleanup. */
BOOST_AUTO_TEST_CASE( immutable_fields )
{
    auto mesh = unitSquare(0.25);
    auto Qh = Pdh<0>(mesh);
    auto Vh = Pchv<1>(mesh);
    auto pid = Qh->element(cst(7.));
    auto velocity = Vh->element(vec(Px(),Py()));
    // Force an internal time-set ID other than one. Case-file IDs are local
    // to the standalone case, and static variables have an implicit ID of one
    // in the supported legacy reader.
    auto sibling = exporter(_mesh=mesh, _name="unused_sibling", _geo="static");
    auto e = exporter(_mesh=mesh, _name="immutable_fields", _geo="static");
    BOOST_REQUIRE_EQUAL(e->type(), "ensightgold");
    bool const native=e->supportsNativeStaticFields();
    e->add("pid",pid);
    e->add("reference_velocity",velocity);
    BOOST_CHECK_EQUAL(e->defaultTimeSet()->numberOfSteps(),0);
    BOOST_CHECK_THROW(e->add("pid",pid), std::invalid_argument);
    pid.on(_range=elements(mesh),_expr=cst(99.));
    auto snapshot = e->staticFields();
    BOOST_CHECK_CLOSE(snapshot->beginElement()->second.second[0][0]->min(), 7., 1e-8);
    e->save(); // static-only output, before any temporal Step exists
    if (Environment::isMasterRank())
    {
        auto text = readMetadata(fs::path(e->path()) / "immutable_fields.case");
        std::ofstream fixture((fs::path(e->path()) / "immutable_fields_static_only.case").string());
        fixture << text;
    }
    BOOST_CHECK_THROW(e->restart(0), std::logic_error);
    auto first = e->step(0);
    BOOST_CHECK_THROW(first->add("pid",pid), std::invalid_argument);
    BOOST_CHECK_THROW(e->add("late",pid), std::invalid_argument);
    for (int k=0; k<3; ++k)
    {
        velocity.on(_range=elements(mesh),_expr=(1.+k)*vec(Px(),Py()));
        e->step(k)->add("velocity",velocity);
        e->step(k)->add("viscosity",1.25,true);
        e->step(k)->add("energy",double(k));
        e->save();
        BOOST_CHECK(e->staticFieldsWritten());
        BOOST_CHECK_EQUAL(snapshot->beginElement()->second.second.empty(),native);
    }
    e->save(); // no new fields: static payload must not be rewritten
    Environment::worldComm().barrier();
    if (Environment::isMasterRank())
    {
        auto base = fs::path(e->path());
        auto text = readMetadata(base / "immutable_fields.case");
        BOOST_CHECK(text.find("VARIABLE:\n") != std::string::npos);
        BOOST_CHECK(text.find("time set:        1\n") != std::string::npos);
        if (native)
        {
            BOOST_CHECK(text.find("scalar per element: pid immutable_fields.pid.scl\n") != std::string::npos);
            BOOST_CHECK(text.find("vector per node: reference_velocity immutable_fields.reference_velocity.vec\n") != std::string::npos);
        }
        else BOOST_CHECK(text.find("scalar per element: 1 1 pid ") != std::string::npos);
        BOOST_CHECK(text.find("constant per case: viscosity 1.25\n") != std::string::npos);
        BOOST_CHECK_EQUAL(readMetadata(base / "immutable_fields.energy.cst"), "0\n1\n2\n");
        BOOST_CHECK(fs::exists(base / "immutable_fields.geo"));
        BOOST_CHECK(!fs::exists(base / "immutable_fields.1.geo"));
        BOOST_CHECK(fs::exists(base / "immutable_fields.static-fields"));
        BOOST_CHECK(!fs::exists(base / "immutable_fields.timeset.static-fields"));
        BOOST_CHECK(!fs::exists(base / "immutable_fields.pid.scl.0001"));
        BOOST_CHECK(text.find("APPENDED_CASEFILES") == std::string::npos);
    }
    auto restarted = exporter(_mesh=mesh,_name="immutable_fields",_geo="static");
    BOOST_CHECK_THROW(restarted->restart(1),std::logic_error);
    auto displacement = Vh->element(vec(cst(.01),cst(0.)));
    meshMove(mesh,displacement);
    BOOST_CHECK_THROW(e->save(),std::logic_error);
}

/** @brief MPI schema disagreement is rejected collectively, before interpolation. */
BOOST_AUTO_TEST_CASE( inconsistent_static_schema )
{
    auto mesh = unitSquare(0.25);
    auto field = Pdh<0>(mesh)->element(cst(1.));
    auto e = exporter(_mesh=mesh,_name="schema",_geo="static");
    if (Environment::worldComm().globalSize() > 1)
        BOOST_CHECK_THROW(e->add(Environment::isMasterRank() ? "a" : "b",field),std::invalid_argument);
    BOOST_CHECK_THROW(exporter(_mesh=mesh,_geo="misspelled"),std::invalid_argument);
    auto moving = exporter(_mesh=mesh,_geo="change");
    moving->add("pid",field); // Geometry policy does not set field lifetime.
    BOOST_CHECK_EQUAL(moving->defaultTimeSet()->numberOfSteps(),0);
}

/** @brief Different time sets keep distinct filenames, schemas and active counts. */
BOOST_AUTO_TEST_CASE( timesets_and_frequency )
{
    auto mesh = unitSquare(0.25);
    auto field = Pdh<0>(mesh)->element(cst(1.));
    auto e = exporter(_mesh=mesh,_name="first",_geo="static");
    e->setFreq(2);
    auto second = e->addTimeSet("second");
    e->timeSet(second)->setMesh(mesh);
    for (int k=0; k<6; ++k)
    {
        field.on(_range=elements(mesh),_expr=cst(double(k)));
        e->step(k,0)->add("value",field);
        e->step(k,0)->add("energy",double(k));
        if (k < 4) e->step(k,second)->add("other",field);
        e->save(); // the second time set's cleaned schema is used again at k=4
    }
    Environment::worldComm().barrier();
    if (Environment::isMasterRank())
    {
        auto base = fs::path(e->path());
        auto first = readMetadata(base / "first.case");
        auto other = readMetadata(base / "second.case");
        BOOST_CHECK(first.find("number of steps: 3\n") != std::string::npos);
        BOOST_CHECK(other.find("number of steps: 2\n") != std::string::npos);
        BOOST_CHECK(fs::exists(base / "first.timeset"));
        BOOST_CHECK(fs::exists(base / "second.timeset"));
        BOOST_CHECK_EQUAL(readMetadata(base / "first.energy.cst"),"0\n2\n4\n");
        auto sos = readMetadata(base / ("second-paraview-" + std::to_string(Environment::worldComm().globalSize()) + ".sos"));
        BOOST_CHECK(sos.find("casefile: second.case") != std::string::npos);
        BOOST_CHECK(sos.find("data_path: " + fs::absolute(base).string()) != std::string::npos);
    }
}

/** @brief All exporter add overloads share a time-neutral dataset and preserve their values. */
BOOST_AUTO_TEST_CASE( dataset_overloads )
{
    auto mesh=unitSquare(.25);
    using mesh_type=typename decltype(mesh)::element_type;
    using field_set_type=Feel::detail::ExportFieldSet<mesh_type>;
    static_assert(!std::is_same_v<field_set_type,typename Feel::detail::TimeSet<mesh_type>::step_type>);
    static_assert(!HasTime<field_set_type>);
    auto scalar=Pch<1>(mesh)->element(Px()+2*Py());
    auto tensor=Pchm<1>(mesh)->element(mat<2,2>(cst(1.),cst(2.),cst(3.),cst(4.)));
    auto symmetric=Pchms<1>(mesh)->element(mat<2,2>(cst(5.),cst(6.),cst(6.),cst(7.)));
    using mixed_space_type=FunctionSpace<mesh_type,bases<Lagrange<1,Scalar>,Lagrange<1,Vectorial>>>;
    auto mixedSpace=mixed_space_type::New(_mesh=mesh);
    auto mixed=mixedSpace->element();
    mixed.template element<0>().on(_range=elements(mesh),_expr=cst(4.));
    mixed.template element<1>().on(_range=elements(mesh),_expr=vec(cst(2.),cst(3.)));
    auto e=exporter(_mesh=mesh,_name="dataset_overloads",_geo="static");
    e->add("scalar",scalar);
    e->add("pointer",std::make_shared<decltype(scalar)>(scalar));
    e->add("tensor",tensor);
    e->add("symmetric",symmetric);
    e->add("mixed",mixed);
    e->add("collision_1",scalar);
    BOOST_CHECK_THROW(e->add("collision",mixed),std::invalid_argument);
    BOOST_CHECK(!e->staticFields()->names().count("collision_0"));
    e->add("expression",Px()+2*Py());
    e->add("vector_expression",vec(Px(),Py()));
    e->add("tensor_expression",mat<2,2>(cst(1.),cst(2.),cst(3.),cst(4.)));
    e->add("element_expression",cst(8.),"element");
    e->add("range_expression",cst(9.),elements(mesh),"element");
    e->add("boundary_expression",cst(10.),boundaryfaces(mesh));
    e->add("both",scalar,std::set<std::string>{"nodal","element"});
    e->add("u ;)[ ",scalar);
    e->add("viscosity",1.25);
    e->add("explicit_constant",2.5,true);
    e->addRegions();
    BOOST_CHECK_EQUAL(e->defaultTimeSet()->numberOfSteps(),0);
    BOOST_CHECK(e->staticFields()->names().count("u"));
    BOOST_CHECK(e->staticFields()->names().count("both_n"));
    BOOST_CHECK(e->staticFields()->names().count("both_e"));
    BOOST_CHECK_THROW(e->add("pid",scalar),std::invalid_argument);
    BOOST_CHECK_THROW(e->add("both_n",scalar),std::invalid_argument);
    BOOST_CHECK_THROW(e->add("not_constant",1.,false),std::invalid_argument);
    scalar.on(_range=elements(mesh),_expr=cst(-123.));
    tensor.on(_range=elements(mesh),_expr=mat<2,2>(cst(-1.),cst(-1.),cst(-1.),cst(-1.)));
    symmetric.on(_range=elements(mesh),_expr=mat<2,2>(cst(-1.),cst(-1.),cst(-1.),cst(-1.)));
    e->save();
    BOOST_CHECK_EQUAL(e->defaultTimeSet()->numberOfSteps(),e->supportsNativeStaticFields()?0:1);
    if (Environment::isMasterRank())
    {
        auto text=readMetadata(fs::path(e->path())/"dataset_overloads.case");
        BOOST_CHECK(text.find("constant per case: viscosity 1.25\n")!=std::string::npos);
        BOOST_CHECK(text.find("constant per case: explicit_constant 2.5\n")!=std::string::npos);
    }
}

/** @brief Invalid registration, layout changes and direct time-set collisions fail coherently. */
BOOST_AUTO_TEST_CASE( dataset_contract_errors )
{
    auto mesh=unitSquare(.25);
    auto other=unitSquare(.3);
    auto field=Pdh<0>(mesh)->element(cst(1.));
    auto foreign=Pdh<0>(other)->element(cst(1.));
    auto e=exporter(_mesh=mesh,_name="contract_errors",_geo="static");
    BOOST_CHECK_THROW(e->add("foreign",foreign),std::invalid_argument);
    std::shared_ptr<decltype(field)> nullField;
    BOOST_CHECK_THROW(e->add("null",nullField),std::invalid_argument);
    BOOST_CHECK_THROW(e->add("",field),std::invalid_argument);
    BOOST_CHECK_THROW(e->add("bad_rep",field,"invalid"),std::invalid_argument);
    BOOST_CHECK_THROW(e->add("bad_range",cst(1.),elements(other)),std::invalid_argument);
    BOOST_CHECK_THROW(e->add("nan",std::numeric_limits<double>::quiet_NaN()),std::invalid_argument);
    if (Environment::worldComm().globalSize()>1)
    {
        BOOST_CHECK_THROW(e->add("different_rep",field,Environment::isMasterRank()?"nodal":"element"),std::invalid_argument);
        BOOST_CHECK_THROW(e->add("different_constant",double(Environment::worldComm().globalRank())),std::invalid_argument);
        BOOST_CHECK_THROW(e->add("different_validity",field,Environment::isMasterRank()?"invalid":"element"),std::invalid_argument);
    }
    e->add("fixed",field);
    using mesh_type=typename decltype(mesh)::element_type;
    auto gold=std::dynamic_pointer_cast<ExporterEnsightGold<mesh_type,1>>(e);
    BOOST_CHECK_THROW((std::make_shared<ExporterEnsightGold<mesh_type,1>>(*gold)),std::logic_error);
    BOOST_CHECK_THROW(e->setMesh(other),std::logic_error);
    BOOST_CHECK_THROW(e->setMesh(mesh,EXPORTER_GEOMETRY_CHANGE),std::logic_error);
    BOOST_CHECK_THROW(e->setPrefix("renamed"),std::logic_error);
    BOOST_CHECK_THROW(e->setType("gmsh"),std::logic_error);
    BOOST_CHECK_THROW(e->setMeshFragmentation(e->meshFragmentation()),std::logic_error);
    auto direct=e->defaultTimeSet()->step(0);
    direct->add("fixed",field); // Direct access bypasses exporter.step's eager name reservation.
    BOOST_CHECK_THROW(e->save(),std::invalid_argument); // Still rejected before I/O.

    auto late=exporter(_mesh=mesh,_name="late_fields",_geo="static");
    late->step(0)->add("value",field);
    BOOST_CHECK_THROW(late->add("late",field),std::invalid_argument);
    auto extra=exporter(_mesh=mesh,_name="extra_timeset",_geo="static");
    extra->add("fixed",field);
    extra->timeSet(extra->addTimeSet("other"))->setMesh(mesh);
    BOOST_CHECK_THROW(extra->save(),std::invalid_argument);
    auto rebound=exporter(_mesh=mesh,_name="rebound",_geo="static");
    rebound->add("fixed",field);
    rebound->defaultTimeSet()->setMesh(other);
    BOOST_CHECK_THROW(rebound->save(),std::invalid_argument);
    auto renamed=exporter(_mesh=mesh,_name="renamed_sequence",_geo="static");
    renamed->add("fixed",field);
    renamed->defaultTimeSet()->setName("different_sequence");
    BOOST_CHECK_THROW(renamed->save(),std::invalid_argument);
}

//! \brief Partial FE snapshots remain static, own their values and zero untouched DOFs.
BOOST_AUTO_TEST_CASE( dataset_partial_support )
{
    auto mesh = unitSquare( .25 );
    using mesh_type = typename decltype( mesh )::element_type;
    auto range = elements( mesh, Px() < cst( .5 ),
                           _selector = select_elements_from_expression::with_value, _value = 1 );
    auto scalar = Pch<1>( mesh, range )->element( cst( 2. ) );
    auto quadratic = Pch<2>( mesh, range )->element( cst( 3. ) );
    auto cell = Pdh<0>( mesh, range )->element( cst( 7. ) );
    auto vector = Pchv<1>( mesh, range )->element( vec( cst( 4. ), cst( 5. ) ) );
    auto tensorSpace = Pchm_type<mesh_type, 1>::New( _mesh = mesh, _range = range );
    auto tensor = tensorSpace->element(
        mat<2, 2>( cst( 1. ), cst( 2. ), cst( 3. ), cst( 4. ) ) );
    using mixed_space_type = FunctionSpace<mesh_type, bases<Lagrange<1, Scalar>, Lagrange<1, Vectorial>>>;
    auto partialSupport = std::make_shared<MeshSupport<mesh_type>>( mesh, range );
    auto fullSupport = std::make_shared<MeshSupport<mesh_type>>( mesh );
    auto mixedSpace = mixed_space_type::New(
        _mesh = mesh, _range = boost::fusion::make_vector( partialSupport, fullSupport ) );
    auto mixed = mixedSpace->element();
    mixed.template element<0>().on( _range = range, _expr = cst( 6. ), _close = true );
    mixed.template element<1>().on( _range = elements( mesh ), _expr = vec( cst( 8. ), cst( 9. ) ), _close = true );

    auto e = exporter( _mesh = mesh, _name = "partial_support", _geo = "static" );
    e->add( "scalar", scalar );
    e->add( "quadratic", quadratic );
    e->add( "cell", cell );
    e->add( "vector", vector );
    e->add( "tensor", tensor );
    e->add( "mixed", mixed );
    BOOST_CHECK_EQUAL( e->defaultTimeSet()->numberOfSteps(), 0 );

    // Changing the partial sources must leave the registered snapshots unchanged.
    scalar.on( _range = range, _expr = cst( 20. ) );
    quadratic.on( _range = range, _expr = cst( 30. ) );
    cell.on( _range = range, _expr = cst( 70. ) );
    vector.on( _range = range, _expr = vec( cst( 40. ), cst( 50. ) ) );
    tensor.on( _range = range, _expr = mat<2, 2>( cst( 10. ), cst( 20. ), cst( 30. ), cst( 40. ) ) );
    mixed.template element<0>().on( _range = range, _expr = cst( 60. ) );

    auto expected = Pch<1>( mesh )->element();
    auto checkNodal = [&]( auto const& field, double value )
    {
        BOOST_CHECK( support( field->functionSpace() )->isFullSupport() );
        expected.zero();
        // Synchronize supported DOFs across full-mesh partition boundaries.
        expected.on( _range = range, _expr = cst( value ), _close = true );
        BOOST_CHECK_SMALL( normL2( _range = range, _expr = idv( *field ) - cst( value ) ), 1e-12 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ), _expr = idv( *field ) - idv( expected ) ), 1e-12 );
    };
    auto fields = e->staticFields();
    checkNodal( fields->nodal( "scalar" ).second[0][0], 2. );
    checkNodal( fields->nodal( "quadratic" ).second[0][0], 3. );
    checkNodal( fields->nodal( "vector" ).second[0][0], 4. );
    checkNodal( fields->nodal( "vector" ).second[1][0], 5. );
    for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j )
            checkNodal( fields->nodal( "tensor" ).second[i][j], 1. + 2 * i + j );
    checkNodal( fields->nodal( "mixed_0" ).second[0][0], 6. );
    BOOST_CHECK_CLOSE( fields->nodal( "mixed_1" ).second[0][0]->min(), 8., 1e-8 );
    BOOST_CHECK_CLOSE( fields->nodal( "mixed_1" ).second[1][0]->min(), 9., 1e-8 );
    auto expectedCell = Pdh<0>( mesh )->element();
    expectedCell.on( _range = range, _expr = cst( 7. ), _close = true );
    auto cellSnapshot = fields->element( "cell" ).second[0][0];
    BOOST_CHECK( support( cellSnapshot->functionSpace() )->isFullSupport() );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ), _expr = idv( *cellSnapshot ) - idv( expectedCell ) ), 1e-12 );
    e->save();
}

/** @brief Stationary add/save works for all geometry hints and native/per-step storage. */
BOOST_AUTO_TEST_CASE( stationary_geometry_policies )
{
    auto mesh=unitSquare(.25);
    auto field=Pdh<0>(mesh)->element(cst(11.));
    for (auto const& geometry : {"static","change","change_coords_only"})
    {
        std::string name=std::string("stationary_")+geometry;
        auto e=exporter(_mesh=mesh,_name=name,_geo=geometry);
        e->restart(0); // Missing metadata is a no-op, not an actual restart.
        e->add("value",field);
        e->add("constant",12.5);
        BOOST_CHECK_EQUAL(e->defaultTimeSet()->numberOfSteps(),0);
        e->save();
        e->save();
        BOOST_CHECK_THROW(e->setPath(e->path()+"_moved"),std::logic_error);
        BOOST_CHECK_THROW(e->add("late",field),std::invalid_argument);
        if (Environment::isMasterRank())
        {
            auto text=readMetadata(fs::path(e->path())/(name+".case"));
            BOOST_CHECK(text.find(" change_coords_only")==std::string::npos);
            BOOST_CHECK(text.find("constant per case: constant 12.5\n")!=std::string::npos);
            if (e->supportsNativeStaticFields()) BOOST_CHECK(text.find("TIME:")==std::string::npos);
        }
    }
}

/** @brief Dataset fields remain available with ignored samples and pack boundaries. */
BOOST_AUTO_TEST_CASE( dataset_frequency )
{
    auto mesh=unitSquare(.25);
    auto field=Pdh<0>(mesh)->element(cst(3.));
    auto e=exporter(_mesh=mesh,_name="dataset_frequency",_geo="static");
    e->setFreq(2);
    e->add("fixed",field);
    for (int k=0;k<6;++k)
    {
        field.on(_range=elements(mesh),_expr=cst(double(k)));
        e->step(k)->add("signal",field);
        e->save();
    }
    if (Environment::isMasterRank())
    {
        auto text=readMetadata(fs::path(e->path())/"dataset_frequency.case");
        BOOST_CHECK(text.find("number of steps: 3\n")!=std::string::npos);
        BOOST_CHECK(text.find("fixed")!=std::string::npos);
    }
}

/** @brief FILE_INDEX roundtrip covers batched offsets and reuse on empty files. */
BOOST_AUTO_TEST_CASE( transient_index )
{
    auto path = (fs::path(Environment::exportsRepository()) / "index-roundtrip.bin").string();
    MPI_File file;
    BOOST_REQUIRE_EQUAL(MPI_File_open(Environment::worldComm().comm(),path.c_str(),
                                     MPI_MODE_RDWR|MPI_MODE_CREATE,MPI_INFO_NULL,&file),MPI_SUCCESS);
    BOOST_REQUIRE_EQUAL(MPI_File_set_size(file,0),MPI_SUCCESS);
    Feel::detail::FileIndex written(Environment::worldCommPtr());
    for (int k=0; k<1000; ++k) written.add(80+8*k);
    MPI_Offset offset = 9000;
    written.write(file,offset);
    BOOST_REQUIRE_EQUAL(MPI_File_close(&file),MPI_SUCCESS);
    BOOST_REQUIRE_EQUAL(MPI_File_open(Environment::worldComm().comm(),path.c_str(),
                                     MPI_MODE_RDWR,MPI_INFO_NULL,&file),MPI_SUCCESS);
    Feel::detail::FileIndex loaded(Environment::worldCommPtr());
    loaded.read(file);
    BOOST_CHECK_EQUAL(loaded.numberOfBlock(),1000);
    BOOST_CHECK_EQUAL(loaded.nextFreePosFile(),9000);
    BOOST_CHECK(loaded.fileBlocks() == written.fileBlocks());
    BOOST_REQUIRE_EQUAL(MPI_File_set_size(file,0),MPI_SUCCESS);
    loaded.read(file);
    BOOST_CHECK(!loaded.defined());
    BOOST_CHECK_EQUAL(loaded.nextFreePosFile(),-1);
    BOOST_REQUIRE_EQUAL(MPI_File_close(&file),MPI_SUCCESS);
}
