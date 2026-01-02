/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2025-01-09

  Copyright (C) 2025 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_repository.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2025-01-09
 */
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE repository testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/repository.hpp>
#include <fstream>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description repositoryoptions( "Test Repository options" );
    repositoryoptions.add_options()
        ( "test.custom.basedir", Feel::po::value<std::string>()->default_value( "/tmp/feelpp-test-custom" ), "base directory for custom tests" )
        ( "test.custom.subdir", Feel::po::value<std::string>()->default_value( "mysubdir" ), "subdirectory for custom tests" )
    ;
    return repositoryoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_repository" ,
                           "test_repository" ,
                           "0.1",
                           "Repository class tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2025 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( repository )

/**
 * Test basic Repository construction and default values
 */
BOOST_AUTO_TEST_CASE( test_repository_default )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_default" );

    Repository repo;
    // Need to configure before accessing root()
    repo.configure();
    
    BOOST_CHECK( repo.root().empty() == false );
    BOOST_CHECK( fs::exists( repo.root() ) );
    BOOST_CHECK( fs::is_directory( repo.root() ) );
    
    BOOST_TEST_MESSAGE( "test_repository_default done" );
}

/**
 * Test global repository location
 */
BOOST_AUTO_TEST_CASE( test_repository_global )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_global" );

    auto config = globalRepository( "test-global-repo" );
    Repository repo( config );
    repo.configure();

    BOOST_CHECK( repo.isGlobal() );
    BOOST_CHECK( !repo.isRelative() );
    BOOST_CHECK( !repo.isAbsolute() );
    BOOST_CHECK( !repo.isGit() );
    BOOST_CHECK( !repo.isCustom() );
    
    // Global repository should be under home directory or configured global_root
    fs::path root = repo.root();
    BOOST_CHECK( fs::exists( root ) );
    BOOST_CHECK( fs::is_directory( root ) );
    
    // Check subdirectories are created
    BOOST_CHECK( fs::exists( repo.geo() ) );
    BOOST_CHECK( fs::exists( repo.exprs() ) );
    BOOST_CHECK( fs::exists( repo.logs() ) );

    BOOST_TEST_MESSAGE( fmt::format("Global repository root: {}", root.string() ) );
    BOOST_TEST_MESSAGE( "test_repository_global done" );
}

/**
 * Test relative repository location
 */
BOOST_AUTO_TEST_CASE( test_repository_relative )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_relative" );

    auto config = localRepository( "test-relative-repo" );
    Repository repo( config );
    repo.configure();

    BOOST_CHECK( !repo.isGlobal() );
    BOOST_CHECK( repo.isRelative() );
    BOOST_CHECK( repo.isLocal() ); // isLocal() should match isRelative()
    BOOST_CHECK( !repo.isAbsolute() );
    BOOST_CHECK( !repo.isGit() );
    BOOST_CHECK( !repo.isCustom() );
    
    // Relative repository should be under current working directory
    fs::path root = repo.root();
    BOOST_CHECK( fs::exists( root ) );
    BOOST_CHECK( fs::is_directory( root ) );
    
    // Root should contain "feelppdb" (the default)
    BOOST_CHECK( root.string().find("feelppdb") != std::string::npos );

    BOOST_TEST_MESSAGE( fmt::format("Relative repository root: {}", root.string() ) );
    BOOST_TEST_MESSAGE( "test_repository_relative done" );
}

/**
 * Test absolute repository location
 */
BOOST_AUTO_TEST_CASE( test_repository_absolute )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_absolute" );

    fs::path abs_path = "/tmp/feelpp-test-absolute-repo";
    
    // Clean up if exists
    if ( fs::exists( abs_path ) )
        fs::remove_all( abs_path );

    Repository::Config config( abs_path, Location::absolute );
    Repository repo( config );
    repo.configure();

    BOOST_CHECK( !repo.isGlobal() );
    BOOST_CHECK( !repo.isRelative() );
    BOOST_CHECK( repo.isAbsolute() );
    BOOST_CHECK( !repo.isGit() );
    BOOST_CHECK( !repo.isCustom() );
    
    // Absolute repository should match exactly what we specified
    BOOST_CHECK_EQUAL( repo.root(), abs_path );
    BOOST_CHECK( fs::exists( repo.root() ) );
    BOOST_CHECK( fs::is_directory( repo.root() ) );

    // Clean up
    if ( Environment::isMasterRank() )
        fs::remove_all( abs_path );

    BOOST_TEST_MESSAGE( "test_repository_absolute done" );
}

/**
 * Test git repository location
 */
BOOST_AUTO_TEST_CASE( test_repository_git )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_git" );

    // Create a temporary git repository for testing
    fs::path test_git_dir = "/tmp/feelpp-test-git-repo";
    
    if ( Environment::isMasterRank() )
    {
        // Clean up if exists
        if ( fs::exists( test_git_dir ) )
            fs::remove_all( test_git_dir );

        fs::create_directories( test_git_dir );
        fs::create_directories( test_git_dir / ".git" );
        
        // Create a subdirectory where we'll run the test
        fs::create_directories( test_git_dir / "subdir" / "deeper" );
    }
    Environment::worldComm().barrier();

    // Save current directory
    fs::path original_dir = fs::current_path();
    
    try
    {
        // Change to subdirectory within git repo
        fs::current_path( test_git_dir / "subdir" / "deeper" );
        
        Repository::Config config( "test-git-repo", Location::git );
        Repository repo( config );
        repo.configure();

        BOOST_CHECK( !repo.isGlobal() );
        BOOST_CHECK( !repo.isRelative() );
        BOOST_CHECK( !repo.isAbsolute() );
        BOOST_CHECK( repo.isGit() );
        BOOST_CHECK( !repo.isCustom() );
        
        // Git repository should be at the git root
        fs::path root = repo.root();
        BOOST_CHECK( root.string().find( test_git_dir.string() ) != std::string::npos );
        BOOST_CHECK( fs::exists( root ) );

        BOOST_TEST_MESSAGE( fmt::format("Git repository root: {}", root.string() ) );
    }
    catch ( const std::exception& e )
    {
        BOOST_TEST_MESSAGE( fmt::format("Exception in git test: {}", e.what() ) );
        fs::current_path( original_dir );
        throw;
    }
    
    // Restore original directory
    fs::current_path( original_dir );
    
    // Clean up
    if ( Environment::isMasterRank() && fs::exists( test_git_dir ) )
        fs::remove_all( test_git_dir );

    BOOST_TEST_MESSAGE( "test_repository_git done" );
}

/**
 * Test git repository error when not in git directory
 */
// BOOST_AUTO_TEST_CASE( test_repository_git_error )
// {
//     using namespace Feel;
//     BOOST_TEST_MESSAGE( "test_repository_git_error" );
// 
//     fs::path test_non_git_dir = "/tmp/feelpp-test-non-git";
//     fs::path original_dir = fs::current_path();
//     
//     if ( Environment::isMasterRank() )
//     {
//         if ( fs::exists( test_non_git_dir ) )
//             fs::remove_all( test_non_git_dir );
//         fs::create_directories( test_non_git_dir );
//     }
//     Environment::worldComm().barrier();
// 
//     bool exception_caught = false;
//     
//     try
//     {
//         fs::current_path( test_non_git_dir );
//         
//         Repository::Config config( "test-should-fail", Location::git );
//         Repository repo( config );
//         
//         // This should throw because we're not in a git repository
//         try 
//         {
//             repo.configure();
//             BOOST_TEST_MESSAGE( "ERROR: Expected exception was not thrown!" );
//         }
//         catch ( const std::invalid_argument& e )
//         {
//             BOOST_TEST_MESSAGE( fmt::format("Expected exception caught: {}", e.what() ) );
//             exception_caught = true;
//         }
//     }
//     catch ( ... )
//     {
//         BOOST_TEST_MESSAGE( "Unexpected exception during test setup" );
//         fs::current_path( original_dir );
//         throw;
//     }
//     
//     // Restore directory before any checks
//     fs::current_path( original_dir );
//     
//     // Now do the actual boost check
//     BOOST_CHECK( exception_caught );
//     
//     // Clean up
//     Environment::worldComm().barrier();
//     if ( Environment::isMasterRank() && fs::exists( test_non_git_dir ) )
//         fs::remove_all( test_non_git_dir );
// 
//     BOOST_TEST_MESSAGE( "test_repository_git_error done" );
// }

/**
 * Test custom repository location with lambda
 */
BOOST_AUTO_TEST_CASE( test_repository_custom )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_custom" );

    // Test custom repository with a lambda that uses options
    auto config = customRepository( "fallback-dir", []() {
        // Access Feel++ options within the lambda
        std::string basedir = Feel::soption(_name="test.custom.basedir");
        std::string subdir = Feel::soption(_name="test.custom.subdir");
        fs::path custom_path = fs::path(basedir) / subdir;
        
        VLOG(2) << fmt::format("Custom lambda computed path: {}", custom_path.string());
        return custom_path;
    });

    Repository repo( config );
    repo.configure();

    BOOST_CHECK( !repo.isGlobal() );
    BOOST_CHECK( !repo.isRelative() );
    BOOST_CHECK( !repo.isAbsolute() );
    BOOST_CHECK( !repo.isGit() );
    BOOST_CHECK( repo.isCustom() );
    
    // Custom repository should use the computed path
    fs::path root = repo.root();
    std::string expected_basedir = Feel::soption(_name="test.custom.basedir");
    std::string expected_subdir = Feel::soption(_name="test.custom.subdir");
    
    BOOST_CHECK( root.string().find( expected_basedir ) != std::string::npos );
    BOOST_CHECK( root.string().find( expected_subdir ) != std::string::npos );
    BOOST_CHECK( fs::exists( root ) );
    BOOST_CHECK( fs::is_directory( root ) );

    BOOST_TEST_MESSAGE( fmt::format("Custom repository root: {}", root.string() ) );
    
    // Clean up
    if ( Environment::isMasterRank() )
    {
        fs::path cleanup_path = fs::path(expected_basedir);
        if ( fs::exists( cleanup_path ) )
            fs::remove_all( cleanup_path );
    }

    BOOST_TEST_MESSAGE( "test_repository_custom done" );
}

/**
 * Test custom repository with environment-based location
 */
BOOST_AUTO_TEST_CASE( test_repository_custom_env )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_custom_env" );

    // Set an environment variable for this test
    setenv("FEELPP_TEST_CUSTOM_ROOT", "/tmp/feelpp-custom-env-test", 1);

    auto config = customRepository( "fallback", []() {
        const char* env_path = getenv("FEELPP_TEST_CUSTOM_ROOT");
        if ( env_path )
            return fs::path(env_path) / "custom-from-env";
        else
            return fs::path("/tmp/fallback-custom");
    });

    Repository repo( config );
    repo.configure();

    BOOST_CHECK( repo.isCustom() );
    
    fs::path root = repo.root();
    BOOST_CHECK( root.string().find("feelpp-custom-env-test") != std::string::npos );
    BOOST_CHECK( root.string().find("custom-from-env") != std::string::npos );
    BOOST_CHECK( fs::exists( root ) );

    BOOST_TEST_MESSAGE( fmt::format("Custom env repository root: {}", root.string() ) );
    
    // Clean up
    if ( Environment::isMasterRank() && fs::exists( root ) )
        fs::remove_all( root.parent_path() );
    
    unsetenv("FEELPP_TEST_CUSTOM_ROOT");

    BOOST_TEST_MESSAGE( "test_repository_custom_env done" );
}

/**
 * Test that custom repository callback takes precedence over directory parameter
 */
BOOST_AUTO_TEST_CASE( test_repository_custom_ignores_directory )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_custom_ignores_directory" );

    std::string callback_result = "/tmp/feelpp-callback-path";
    
    auto config = customRepository( "this-should-be-ignored", [callback_result]() {
        return fs::path(callback_result);
    });

    Repository repo( config );
    
    // Configure with a different directory - should be ignored for custom location
    repo.configure( "different-directory", Location::custom );

    BOOST_CHECK( repo.isCustom() );
    
    fs::path root = repo.root();
    
    // Should use callback result, not the directory parameter
    BOOST_CHECK( root == fs::path(callback_result) );
    BOOST_CHECK( root.string().find("callback-path") != std::string::npos );
    BOOST_CHECK( root.string().find("this-should-be-ignored") == std::string::npos );
    BOOST_CHECK( root.string().find("different-directory") == std::string::npos );
    
    BOOST_TEST_MESSAGE( fmt::format("Callback correctly took precedence: {}", root.string() ) );
    
    // Clean up
    if ( Environment::isMasterRank() && fs::exists( root ) )
        fs::remove_all( root );

    BOOST_TEST_MESSAGE( "test_repository_custom_ignores_directory done" );
}

/**
 * Test custom repository error handling
 */
BOOST_AUTO_TEST_CASE( test_repository_custom_error )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_custom_error" );

    // Test with a lambda that throws
    auto config = customRepository( "fallback-dir", []() -> fs::path {
        throw std::runtime_error("Intentional error in custom callback");
    });

    Repository repo( config );
    
    // Should not throw, but use fallback
    BOOST_CHECK_NO_THROW( repo.configure() );
    
    // Should have used fallback
    fs::path root = repo.root();
    BOOST_CHECK( root.string().find("fallback-dir") != std::string::npos );
    BOOST_CHECK( fs::exists( root ) );

    BOOST_TEST_MESSAGE( fmt::format("Custom error fallback root: {}", root.string() ) );

    BOOST_TEST_MESSAGE( "test_repository_custom_error done" );
}

/**
 * Test repository directory structure and appenders
 */
BOOST_AUTO_TEST_CASE( test_repository_structure )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_structure" );

    auto config = localRepository( "test-structure" );
    Repository repo( config );
    repo.configure();

    // Check that directory structure is correct
    fs::path root = repo.root();
    fs::path dir = repo.directory();
    fs::path dir_no_appenders = repo.directoryWithoutAppenders();
    fs::path rel_dir = repo.relativeDirectory();

    BOOST_CHECK( fs::exists( root ) );
    BOOST_CHECK( fs::exists( dir ) );
    
    // With append_np enabled by default, directory should have np_X
    std::string dir_str = dir.string();
    if ( repo.config().append_np )
    {
        BOOST_CHECK( dir_str.find("np_") != std::string::npos );
    }
    
    // directoryWithoutAppenders should not have np_X
    std::string dir_no_app_str = dir_no_appenders.string();
    
    BOOST_TEST_MESSAGE( fmt::format("Directory: {}", dir.string() ) );
    BOOST_TEST_MESSAGE( fmt::format("Directory without appenders: {}", dir_no_appenders.string() ) );
    BOOST_TEST_MESSAGE( fmt::format("Relative directory: {}", rel_dir.string() ) );

    BOOST_TEST_MESSAGE( "test_repository_structure done" );
}

/**
 * Test repository configuration modification
 */
BOOST_AUTO_TEST_CASE( test_repository_reconfig )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_reconfig" );

    Repository repo;
    
    // Initially configure as relative
    repo.configure( "test-reconfig-1", Location::relative );
    fs::path root1 = repo.root();
    BOOST_CHECK( repo.isRelative() );
    BOOST_CHECK( fs::exists( root1 ) );
    
    // Reconfigure as absolute
    fs::path abs_path = "/tmp/feelpp-reconfig-test";
    if ( Environment::isMasterRank() && fs::exists( abs_path ) )
        fs::remove_all( abs_path );
    Environment::worldComm().barrier();
    
    repo.configure( abs_path, Location::absolute );
    fs::path root2 = repo.root();
    BOOST_CHECK( repo.isAbsolute() );
    BOOST_CHECK_EQUAL( root2, abs_path );
    BOOST_CHECK( fs::exists( root2 ) );
    
    // Roots should be different
    BOOST_CHECK_NE( root1, root2 );
    
    // Clean up
    if ( Environment::isMasterRank() && fs::exists( abs_path ) )
        fs::remove_all( abs_path );

    BOOST_TEST_MESSAGE( "test_repository_reconfig done" );
}

/**
 * Test repository user information
 */
BOOST_AUTO_TEST_CASE( test_repository_user_info )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_user_info" );

    Repository repo;
    repo.configure();

    std::string user_name = repo.userName();
    std::string user_email = repo.userEmail();
    
    // User name should be set (from system or config)
    BOOST_CHECK( !user_name.empty() );
    
    BOOST_TEST_MESSAGE( fmt::format("User name: {}", user_name ) );
    BOOST_TEST_MESSAGE( fmt::format("User email: {}", user_email ) );

    BOOST_TEST_MESSAGE( "test_repository_user_info done" );
}

/**
 * Test repository cd() functionality
 */
BOOST_AUTO_TEST_CASE( test_repository_cd )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_repository_cd" );

    fs::path original_dir = fs::current_path();

    auto config = localRepository( "test-cd" );
    Repository repo( config );
    repo.configure();
    
    fs::path repo_dir = repo.directory();
    
    // Change to repository directory
    repo.cd();
    
    fs::path current_dir = fs::current_path();
    BOOST_CHECK_EQUAL( current_dir, repo_dir );
    
    // Restore original directory
    fs::current_path( original_dir );

    BOOST_TEST_MESSAGE( "test_repository_cd done" );
}

BOOST_AUTO_TEST_SUITE_END()
