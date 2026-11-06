/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#define BOOST_TEST_MODULE ginac jit mpi testsuite

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelvf/vf.hpp>

#include <boost/filesystem.hpp>

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
namespace fs = boost::filesystem;

// Use the expression repository, which is shared by the ranks, rather than /tmp.
struct JitDirectory
{
    fs::path path;

    JitDirectory()
    {
        auto const& world = Feel::Environment::worldComm();
        std::string name;
        if ( world.isMasterRank() )
            name = ( fs::path( Feel::Environment::exprRepository() ) / fs::unique_path( "jit-mpi-%%%%-%%%%" ) ).string();
        boost::mpi::broadcast( world.globalComm(), name, world.masterRank() );
        path = name;
    }

    ~JitDirectory()
    {
        auto const& world = Feel::Environment::worldComm();
        world.globalComm().barrier();
        if ( world.isMasterRank() )
            fs::remove_all( path );
    }
};

int sourcePrints = 0;

void countSourcePrints( GiNaC::ex const& value, GiNaC::print_context const& context )
{
    ++sourcePrints;
    context.s << "(";
    value.print( context );
    context.s << ")";
}

GiNaC::ex trackedExpression( GiNaC::ex const& value )
{
    static unsigned serial = GiNaC::function::register_new(
        GiNaC::function_options( "jit_source_writer", 1 ).print_func<GiNaC::print_csrc_double>( countSourcePrints ) );
    return GiNaC::function( serial, value );
}
}

BOOST_AUTO_TEST_SUITE( ginac_jit_mpi )

BOOST_AUTO_TEST_CASE( persistent_expression_is_built_on_ensemble_worldcomm )
{
    JitDirectory root;
    auto& world = Feel::Environment::worldComm();
    int nsplit = ( world.globalSize() > 1 && world.globalSize() % 2 == 0 ) ? 2 : 1;
    auto [color, w, wglob] = Feel::Environment::worldCommPtr()->split( nsplit );

    std::vector<GiNaC::symbol> symbols{ GiNaC::symbol( "x" ), GiNaC::symbol( "y" ), GiNaC::symbol( "z" ) };
    std::string expression = "x+" + std::to_string( color + 1 );

    auto e = Feel::vf::expr<2>( expression, symbols, "", *w, root.path.string() );
    e.setParameterValues( { { "x", 3.0 }, { "y", 0.0 }, { "z", 0.0 } } );

    BOOST_CHECK_CLOSE( e.evaluate()( 0, 0 ), 4.0 + color, 1e-12 );
}

BOOST_AUTO_TEST_CASE( only_group_master_generates_source_for_each_overload )
{
    JitDirectory root;
    int nsplit = Feel::Environment::worldComm().globalSize() > 1 ? 2 : 1;
    auto [color, w, wglob] = Feel::Environment::worldCommPtr()->split( nsplit );
    auto const& comm = w->localComm();
    std::string stem = ( root.path / ( "group-" + std::to_string( color ) ) ).string();
    GiNaC::symbol x( "x" ), y( "y" );
    sourcePrints = 0;

    GiNaC::FUNCP_1P one = nullptr;
    GiNaC::compile_ex( trackedExpression( x + color + 1 ), x, one, stem + "-one", comm );
    BOOST_CHECK_CLOSE( one( 3.0 ), 4.0 + color, 1e-12 );

    GiNaC::FUNCP_2P two = nullptr;
    GiNaC::compile_ex( trackedExpression( x + y + color ), x, y, two, stem + "-two", comm );
    BOOST_CHECK_CLOSE( two( 3.0, 4.0 ), 7.0 + color, 1e-12 );

    GiNaC::FUNCP_CUBA cuba = nullptr;
    GiNaC::compile_ex( GiNaC::lst{ trackedExpression( x + y + color ) }, GiNaC::lst{ x, y }, cuba, stem + "-cuba", true, comm );
    int n = 2, nf = 1;
    double args[] = { 3.0, 4.0 }, result = 0;
    cuba( &n, args, &nf, &result );
    BOOST_CHECK_CLOSE( result, 7.0 + color, 1e-12 );

    BOOST_CHECK_EQUAL( sourcePrints, comm.rank() == 0 ? 3 : 0 );
}

BOOST_AUTO_TEST_CASE( ensemble_groups_build_different_numbers_of_expressions )
{
    JitDirectory root;
    int nsplit = Feel::Environment::worldComm().globalSize() > 1 ? 2 : 1;
    auto [color, w, wglob] = Feel::Environment::worldCommPtr()->split( nsplit );
    std::vector<GiNaC::symbol> symbols{ GiNaC::symbol( "x" ) };

    // A collective on MPI_COMM_WORLD would hang when group 1 builds its second expression.
    for ( int i = 0; i <= color; ++i )
    {
        int offset = 100 + 10 * color + i;
        auto e = Feel::vf::expr<2>( "x+" + std::to_string( offset ), symbols, "", *w, root.path.string() );
        e.setParameterValues( { { "x", 3.0 } } );
        BOOST_CHECK_CLOSE( e.evaluate()( 0, 0 ), 3.0 + offset, 1e-12 );
    }
}

BOOST_AUTO_TEST_CASE( temporary_expression_is_built_on_ensemble_worldcomm )
{
    int nsplit = Feel::Environment::worldComm().globalSize() > 1 ? 2 : 1;
    auto [color, w, wglob] = Feel::Environment::worldCommPtr()->split( nsplit );
    GiNaC::symbol x( "x" );
    GiNaC::FUNCP_CUBA fp = nullptr;
    GiNaC::compile_ex( GiNaC::lst{ x + color + 1 }, GiNaC::lst{ x }, fp, "", false, w->localComm() );
    int n = 1;
    double xval = 3, result = 0;
    fp( &n, &xval, &n, &result );
    BOOST_CHECK_CLOSE( result, 4.0 + color, 1e-12 );
}

BOOST_AUTO_TEST_CASE( compiler_and_filesystem_errors_reach_every_group_rank )
{
    JitDirectory root;
    int nsplit = Feel::Environment::worldComm().globalSize() > 1 ? 2 : 1;
    auto [color, w, wglob] = Feel::Environment::worldCommPtr()->split( nsplit );
    auto const& comm = w->localComm();
    std::string stem = ( root.path / ( "failure-" + std::to_string( color ) ) ).string();
    GiNaC::symbol x( "x" ), invalid( "this is invalid C++" );
    GiNaC::FUNCP_CUBA fp = nullptr;

    BOOST_CHECK_EXCEPTION(
        GiNaC::compile_ex( GiNaC::lst{ invalid }, GiNaC::lst{ x }, fp, stem, true, comm ),
        std::runtime_error,
        []( std::runtime_error const& e ) { return std::string( e.what() ).find( "compile_src_file" ) != std::string::npos; } );

    if ( comm.rank() == 0 )
        fs::create_directory( stem + "-locked.lock" );
    comm.barrier();
    BOOST_CHECK_EXCEPTION(
        GiNaC::compile_ex( GiNaC::lst{ x }, GiNaC::lst{ x }, fp, stem + "-locked", true, comm ),
        std::runtime_error,
        []( std::runtime_error const& e ) { return std::string( e.what() ).find( "lock file" ) != std::string::npos; } );

    if ( comm.rank() == 0 )
    {
        std::ofstream blocked( stem + "-parent" );
        blocked << "a file cannot be used as the source directory";
    }
    comm.barrier();
    BOOST_CHECK_THROW(
        GiNaC::compile_ex( GiNaC::lst{ x }, GiNaC::lst{ x }, fp, stem + "-parent/expression", true, comm ),
        std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END()
