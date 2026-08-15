#define BOOST_TEST_MODULE test_hpddm_two_level

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/petschpddm.hpp>
#include <feel/feelalg/preconditionerpetsc.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{
po::options_description makeOptions()
{
    po::options_description options( "test_hpddm_two_level options" );
    options.add( backend_options( "test-hpddm-one-level" ) );
    options.add( backend_options( "test-hpddm-algebraic-two-level" ) );
    return options;
}

AboutData makeAbout()
{
    return AboutData( "test_hpddm_two_level", "test_hpddm_two_level", "0.1",
                      "test one-level and algebraic two-level HPDDM",
                      Feel::AboutData::License_GPL,
                      "Copyright (c) 2026 Feel++ Consortium" );
}

bool hpddmAvailable()
{
    if ( petscHasHpddmSymbol( "PCCreate_HPDDM" ) )
        return true;

    BOOST_TEST_MESSAGE( "Skipping HPDDM runtime test: PCHPDDM is unavailable" );
    return false;
}

po::variables_map hpddmOptions( std::string const& prefix,
                                bool algebraicTwoLevel )
{
    po::variables_map options;
    if ( algebraicTwoLevel )
        options.insert( {
            prefixvm( prefix, "pc-hpddm-levels-1-eps-nev" ),
            po::variable_value( boost::any( 2 ), false )
        } );
    if ( algebraicTwoLevel )
        options.insert( {
            prefixvm( prefix, "pc-hpddm-levels-1-st-pc-type" ),
            po::variable_value( boost::any( std::string( "mat" ) ), false )
        } );
    return options;
}

void runHpddmSetup( std::string const& prefix,
                    bool expectTwoLevels )
{
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto u = Xh->element();
    auto v = Xh->element();
    auto b = backend();
    auto A = b->newMatrix( _trial=Xh, _test=Xh );
    form2( _trial=Xh, _test=Xh, _matrix=A ) =
        integrate( _range=elements( mesh ),
                   _expr=inner( gradt( u ), grad( v ) ) + idt( u )*id( v ) );
    A->close();

    auto p = Feel::preconditioner( _prefix=prefix,
                                   _matrix=A,
                                   _pc=HPDDM_PRECOND,
                                   _backend=b );
    auto pp = toPETSc( p );
    BOOST_REQUIRE( pp );

    BOOST_CHECK( !pp->hasHpddmAuxiliaryMatrix() );
    BOOST_CHECK( !pp->hasHpddmAuxiliaryIS() );

    p->init();
    BOOST_REQUIRE( pp->pc() );

    const char* pcType = nullptr;
    int ierr = PCGetType( pp->pc(), &pcType );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    BOOST_REQUIRE( pcType );
    BOOST_CHECK_EQUAL( std::string( pcType ), std::string( PCHPDDM ) );

    ierr = PCSetUp( pp->pc() );
    BOOST_REQUIRE_MESSAGE( ierr == 0, prefix << " PCHPDDM setup failed" );

    if ( Environment::worldComm().globalSize() > 1 )
    {
        using get_complexities_type = PetscErrorCode (*)( PC, PetscReal*, PetscReal* );
        auto getComplexities = petscHpddmSymbol<get_complexities_type>(
            "PCHPDDMGetComplexities" );
        BOOST_REQUIRE_MESSAGE( getComplexities,
                               "PCHPDDMGetComplexities is unavailable" );
        PetscReal gridComplexity = 0, operatorComplexity = 0;
        ierr = getComplexities( pp->pc(), &gridComplexity, &operatorComplexity );
        BOOST_REQUIRE_EQUAL( ierr, 0 );
        if ( expectTwoLevels )
        {
            BOOST_CHECK_GT( gridComplexity, 1.0 );
            BOOST_CHECK_GT( operatorComplexity, 1.0 );
        }
        else
        {
            BOOST_CHECK_EQUAL( gridComplexity, 1.0 );
            BOOST_CHECK_EQUAL( operatorComplexity, 1.0 );
        }
    }

    // Neither one-level nor algebraic two-level HPDDM requires an auxiliary
    // local Neumann operator.
    BOOST_CHECK( !pp->hasHpddmAuxiliaryMatrix() );
    BOOST_CHECK( !pp->hasHpddmAuxiliaryIS() );
}
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_CASE( test_one_level_hpddm_needs_no_auxiliary_matrix )
{
    auto const options = hpddmOptions( "test-hpddm-one-level", false );
    BOOST_CHECK( options.empty() );

    if ( !hpddmAvailable() )
        return;

    runHpddmSetup( "test-hpddm-one-level", false );
}

BOOST_AUTO_TEST_CASE( test_algebraic_two_level_hpddm_needs_no_auxiliary_matrix )
{
    auto const options = hpddmOptions(
        "test-hpddm-algebraic-two-level", true );
    BOOST_REQUIRE( options.count(
        prefixvm( "test-hpddm-algebraic-two-level", "pc-hpddm-levels-1-eps-nev" ) ) );
    auto const stPcType = prefixvm(
        "test-hpddm-algebraic-two-level", "pc-hpddm-levels-1-st-pc-type" );
    BOOST_REQUIRE( options.count( stPcType ) );
    BOOST_CHECK_EQUAL( options[stPcType].as<std::string>(), "mat" );
    if ( !hpddmAvailable() )
        return;

    runHpddmSetup( "test-hpddm-algebraic-two-level", true );
}
