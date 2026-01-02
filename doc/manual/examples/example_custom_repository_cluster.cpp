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
   \file example_custom_repository_cluster.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2025-01-09
   
   Example: Conditional repository selection for cluster vs local execution
 */

#include <feel/feelcore/environment.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "Cluster Repository Example options" );
    options.add_options()
        ( "storage.project", Feel::po::value<std::string>()->default_value( "myproject" ), 
          "project name for storage organization" )
    ;
    return options.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "custom_repo_cluster",
                           "custom_repo_cluster",
                           "0.1",
                           "Cluster-aware custom repository",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2025 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", 
                     "christophe.prudhomme@feelpp.org", "" );
    return about;
}

int main( int argc, char** argv )
{
    using namespace Feel;

    //
    // Conditional repository selection based on execution environment
    // Automatically detects cluster and uses appropriate storage
    //
    auto config = customRepository( "fallback-dir", []() {
        // Check if we're on a cluster
        const char* cluster = getenv("CLUSTER_NAME");
        
        fs::path result;
        
        if ( cluster )
        {
            Feel::cout << tc::cyan 
                       << "Cluster detected: " << cluster
                       << tc::reset << std::endl;
            
            // On cluster: use scratch space (fast parallel filesystem)
            const char* scratch = getenv("SCRATCH");
            if ( scratch )
            {
                result = fs::path(scratch) / "feelpp" / soption(_name="storage.project");
                Feel::cout << tc::green 
                           << "Using SCRATCH storage: " << result.string()
                           << tc::reset << std::endl;
                return result;
            }
            
            // Fallback for cluster without SCRATCH
            const char* work = getenv("WORK");
            if ( work )
            {
                result = fs::path(work) / "feelpp" / soption(_name="storage.project");
                Feel::cout << tc::green 
                           << "Using WORK storage: " << result.string()
                           << tc::reset << std::endl;
                return result;
            }
        }
        
        // Local machine: use home directory
        result = findHome() / "feelpp-results" / soption(_name="storage.project");
        Feel::cout << tc::green 
                   << "Using local storage: " << result.string()
                   << tc::reset << std::endl;
        
        return result;
    });

    //
    // Initialize Feel++ environment
    //
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout(),
                     _config=config );

    //
    // Display repository information
    //
    Feel::cout << "\nRepository configuration:" << std::endl;
    Feel::cout << "  Root: " << Environment::rootRepository() << std::endl;
    Feel::cout << "  App:  " << Environment::appRepository() << std::endl;
    Feel::cout << "  Logs: " << Environment::logsRepository() << std::endl;
    
    Feel::cout << "\nProject: " << soption(_name="storage.project") << std::endl;

    return 0;
}

/*
 * Usage examples:
 *
 * 1. Local execution (uses home directory):
 *    ./example_custom_repository_cluster --storage.project=simulation1
 *
 * 2. On a cluster with SCRATCH (e.g., SLURM):
 *    export CLUSTER_NAME=mycluster
 *    export SCRATCH=/scratch/username
 *    ./example_custom_repository_cluster --storage.project=large-simulation
 *
 * 3. On a cluster with WORK storage:
 *    export CLUSTER_NAME=mycluster
 *    export WORK=/work/username
 *    ./example_custom_repository_cluster --storage.project=my-project
 *
 * 4. In a SLURM batch script:
 *    #!/bin/bash
 *    #SBATCH --job-name=feelpp_sim
 *    export CLUSTER_NAME=$SLURM_CLUSTER_NAME
 *    # SCRATCH is usually already set by the system
 *    srun ./example_custom_repository_cluster --storage.project=$SLURM_JOB_NAME
 *
 * The application automatically adapts to the execution environment!
 */
