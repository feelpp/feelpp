/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#define BOOST_TEST_MODULE ginac jit cache testsuite

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelvf/vf.hpp>

#include <boost/filesystem.hpp>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace fs = boost::filesystem;

namespace
{
double
evalExpression( std::string const& expression, std::string const& filename, double x )
{
    std::vector<GiNaC::symbol> symbols{ GiNaC::symbol( "x" ) };
    auto e = Feel::vf::expr<2>( expression, symbols, filename );
    e.setParameterValues( { { "x", x } } );
    return e.evaluate()( 0, 0 );
}

std::string
requiredEnv( char const* name )
{
    char const* value = std::getenv( name );
    if ( !value || value[0] == '\0' )
        throw std::runtime_error( std::string( "missing environment variable " ) + name );
    return value;
}

void
waitForStartFile( fs::path const& startFile )
{
    for ( int i = 0; i < 500; ++i )
    {
        if ( fs::exists( startFile ) )
            return;
        std::this_thread::sleep_for( std::chrono::milliseconds( 10 ) );
    }
    throw std::runtime_error( "timed out waiting for concurrent JIT start file" );
}

int
runWorkerProcess( std::string const& executable, std::string const& filename, std::string const& startFile,
                  std::string const& expression, std::string const& expected )
{
    pid_t pid = fork();
    if ( pid == -1 )
        return -1;

    if ( pid == 0 )
    {
        setenv( "FEELPP_GINAC_JIT_WORKER", "1", 1 );
        setenv( "FEELPP_GINAC_JIT_FILENAME", filename.c_str(), 1 );
        setenv( "FEELPP_GINAC_JIT_START_FILE", startFile.c_str(), 1 );
        setenv( "FEELPP_GINAC_JIT_EXPR", expression.c_str(), 1 );
        setenv( "FEELPP_GINAC_JIT_EXPECTED", expected.c_str(), 1 );

        execl( executable.c_str(), executable.c_str(),
               "--run_test=ginac_jit_cache/worker",
               "--log_level=all",
               "--",
               "--directory=testsuite/test_ginac_jit_cache_worker",
               static_cast<char*>( nullptr ) );
        _exit( 127 );
    }

    return pid;
}
}

BOOST_AUTO_TEST_SUITE( ginac_jit_cache )

BOOST_AUTO_TEST_CASE( sequential_same_filename_rebuilds_changed_expression )
{
    fs::path root = fs::temp_directory_path() / ( "feelpp-ginac-jit-sequential-" + std::to_string( getpid() ) );
    fs::create_directories( root );
    fs::path filename = root / "shared_expr";

    double valueA = evalExpression( "x+1", filename.string(), 3.0 );
    BOOST_CHECK_CLOSE( valueA, 4.0, 1e-12 );

    double valueB = evalExpression( "2*x+5", filename.string(), 3.0 );
    BOOST_CHECK_CLOSE( valueB, 11.0, 1e-12 );

    fs::remove_all( root );
}

BOOST_AUTO_TEST_CASE( worker )
{
    if ( !std::getenv( "FEELPP_GINAC_JIT_WORKER" ) )
        return;

    std::string filename = requiredEnv( "FEELPP_GINAC_JIT_FILENAME" );
    std::string startFile = requiredEnv( "FEELPP_GINAC_JIT_START_FILE" );
    std::string expression = requiredEnv( "FEELPP_GINAC_JIT_EXPR" );
    double expected = std::stod( requiredEnv( "FEELPP_GINAC_JIT_EXPECTED" ) );

    waitForStartFile( startFile );

    for ( int i = 0; i < 8; ++i )
    {
        double value = evalExpression( expression, filename, 3.0 );
        BOOST_CHECK_CLOSE( value, expected, 1e-12 );
    }
}

BOOST_AUTO_TEST_CASE( concurrent_same_filename_different_expressions )
{
    std::string executable = boost::unit_test::framework::master_test_suite().argv[0];
    fs::path root = fs::temp_directory_path() / ( "feelpp-ginac-jit-concurrent-" + std::to_string( getpid() ) );
    fs::create_directories( root );

    fs::path filename = root / "shared_expr";
    fs::path startFile = root / "start";

    std::vector<int> children;
    for ( int i = 0; i < 8; ++i )
    {
        bool useA = ( i % 2 == 0 );
        int pid = runWorkerProcess( executable, filename.string(), startFile.string(),
                                    useA ? "x+1" : "2*x+5",
                                    useA ? "4" : "11" );
        BOOST_REQUIRE( pid > 0 );
        children.push_back( pid );
    }

    {
        std::ofstream start( startFile.string(), std::ios::out | std::ios::trunc );
        start << "start\n";
    }

    for ( int pid : children )
    {
        int status = 0;
        BOOST_REQUIRE_EQUAL( waitpid( pid, &status, 0 ), pid );
        BOOST_REQUIRE_MESSAGE( WIFEXITED( status ), "worker " << pid << " did not exit normally" );
        BOOST_CHECK_EQUAL( WEXITSTATUS( status ), 0 );
    }

    int compiledModules = 0;
    for ( fs::directory_iterator it( root ), end; it != end; ++it )
    {
        std::string name = it->path().filename().string();
        if ( name.find( "shared_expr." ) == 0 && it->path().extension() == ".so" )
            ++compiledModules;
    }
    BOOST_CHECK_GE( compiledModules, 2 );

    fs::remove_all( root );
}

BOOST_AUTO_TEST_SUITE_END()
