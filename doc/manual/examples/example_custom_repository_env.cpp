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
   \file example_custom_repository_env.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2025-01-09
   
   Example: Custom repository with different location strategies
   Demonstrates:
   - Default global repository
   - Custom location from environment variable  
   - Git repository location
   - Absolute directory location
   
   Note: Only ONE example can be run at a time since Feel++ Environment
   can only be initialized once per MPI execution.
   
   Usage:
     # Default repository
     ./example_custom_repository_env
     
     # Custom repository from environment variable
     export FEELPP_CUSTOM_ROOT=/path/to/custom
     ./example_custom_repository_env --test=custom
     
     # Git repository
     ./example_custom_repository_env --test=git
     
     # Absolute directory
     ./example_custom_repository_env --test=absolute
 */

#include <feel/feelcore/environment.hpp>
#include <cstdlib>

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "custom_repo_env",
                           "custom_repo_env",
                           "0.1",
                           "Custom repository locations demonstration",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2025 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", 
                     "christophe.prudhomme@feelpp.org", "" );
    return about;
}

int main( int argc, char** argv )
{
    using namespace Feel;

    try {
        // Determine which test to run based on command-line argument
        std::string test_type = "default";
        for (int i = 1; i < argc; ++i) {
            std::string arg(argv[i]);
            if (arg.find("--test=") == 0) {
                test_type = arg.substr(7);
            }
        }

        std::cout << "\n=== Custom Repository Example ===" << std::endl;
        std::cout << "Test type: " << test_type << "\n" << std::endl;

        if (test_type == "custom") {
            //
            // Test: Custom repository from environment variable
            //
            std::cout << "Using custom repository with environment variable callback\n" << std::endl;
            
            auto config = customRepository( "fallback-dir", []() -> fs::path {
                const char* custom_root = std::getenv("FEELPP_CUSTOM_ROOT");
                if ( custom_root )
                {
                    Feel::cout << "  Found FEELPP_CUSTOM_ROOT: " << custom_root << std::endl;
                    return fs::path(custom_root);
                }
                else
                {
                    Feel::cout << "  FEELPP_CUSTOM_ROOT not set, using fallback" << std::endl;
                    return fs::path("/tmp/feelpp-fallback");
                }
            });
            
            Environment env( _argc=argc, _argv=argv, _about=makeAbout(), _config=config );
            
            Feel::cout << "\nRepository paths:" << std::endl;
            Feel::cout << "  Root:  " << Environment::rootRepository() << std::endl;
            Feel::cout << "  App:   " << Environment::appRepository() << std::endl;
            Feel::cout << "  Logs:  " << Environment::logsRepository() << std::endl;
            Feel::cout << "  Exprs: " << Environment::exprRepository() << std::endl;
        }
        else if (test_type == "git") {
            //
            // Test: Git repository location
            // Feel++ will detect .git directory and use <git-root>/feelppdb
            //
            std::cout << "Using git repository location\n" << std::endl;
            
            auto cfg = gitRepository();
            Environment env( _argc=argc, _argv=argv, _about=makeAbout(), _config=cfg );
            
            Feel::cout << "\nRepository paths:" << std::endl;
            Feel::cout << "  Root:  " << Environment::rootRepository() << std::endl;
            Feel::cout << "  App:   " << Environment::appRepository() << std::endl;
            Feel::cout << "  Logs:  " << Environment::logsRepository() << std::endl;
            Feel::cout << "  Exprs: " << Environment::exprRepository() << std::endl;
        }
        else if (test_type == "absolute") {
            //
            // Test: Absolute directory location
            // Useful for temporary or per-run storage
            //
            std::cout << "Using absolute directory location\n" << std::endl;
            
            auto cfg = absoluteRepository("/tmp/feelpp-absolute");
            Environment env( _argc=argc, _argv=argv, _about=makeAbout(), _config=cfg );
            
            Feel::cout << "\nRepository paths:" << std::endl;
            Feel::cout << "  Root:  " << Environment::rootRepository() << std::endl;
            Feel::cout << "  App:   " << Environment::appRepository() << std::endl;
            Feel::cout << "  Logs:  " << Environment::logsRepository() << std::endl;
            Feel::cout << "  Exprs: " << Environment::exprRepository() << std::endl;
        }
        else {
            //
            // Test: Default repository (global location)
            //
            std::cout << "Using default global repository\n" << std::endl;
            
            Environment env( _argc=argc, _argv=argv, _about=makeAbout() );
            
            Feel::cout << "\nRepository paths:" << std::endl;
            Feel::cout << "  Root:  " << Environment::rootRepository() << std::endl;
            Feel::cout << "  App:   " << Environment::appRepository() << std::endl;
            Feel::cout << "  Logs:  " << Environment::logsRepository() << std::endl;
            Feel::cout << "  Exprs: " << Environment::exprRepository() << std::endl;
        }
        
        Feel::cout << "\n=== Test completed successfully ===" << std::endl;
    }
    catch( std::exception const& e )
    {
        Feel::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

/*
 * Additional Usage Examples:
 *
 * 1. Default repository (global location):
 *    ./example_custom_repository_env
 *
 * 2. Custom location via environment variable:
 *    export FEELPP_CUSTOM_ROOT=/mnt/data/feelpp
 *    ./example_custom_repository_env --test=custom
 *
 * 3. Git repository (finds .git and uses git-root/feelppdb):
 *    ./example_custom_repository_env --test=git
 *
 * 4. Absolute directory:
 *    ./example_custom_repository_env --test=absolute
 *
 * 5. In a batch script (SLURM example):
 *    #!/bin/bash
 *    #SBATCH --job-name=feelpp_job
 *    export FEELPP_CUSTOM_ROOT=$SCRATCH/feelpp
 *    srun ./example_custom_repository_env --test=custom
 *
 * 6. With MPI:
 *    mpirun -np 4 ./example_custom_repository_env --test=git
 */
