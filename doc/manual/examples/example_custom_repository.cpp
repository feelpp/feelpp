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
   \file example_custom_repository.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2025-01-09
   
   Example demonstrating custom repository location using callbacks
 */

#include <feel/feelcore/environment.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "Custom Repository Example options" );
    options.add_options()
        ( "storage.base", Feel::po::value<std::string>()->default_value( "/tmp/feelpp-storage" ), 
          "base directory for storage" )
        ( "storage.project", Feel::po::value<std::string>()->default_value( "myproject" ), 
          "project name for storage organization" )
        ( "storage.user", Feel::po::value<std::string>()->default_value( "" ), 
          "user name (empty = auto-detect)" )
    ;
    return options.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "custom_repo_example",
                           "custom_repo_example",
                           "0.1",
                           "Example of custom repository location",
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
    // Custom repository using parsed options
    // This is the most powerful approach - the callback can access all Feel++ options
    //
    auto config = customRepository( "fallback-dir", []() {
        // Access parsed options within the callback
        std::string base = soption(_name="storage.base");
        std::string project = soption(_name="storage.project");
        std::string user = soption(_name="storage.user");
        
        // Auto-detect user if not specified
        if ( user.empty() )
            user = findUser();
        
        // Construct path: /storage/base/user/project
        fs::path result = fs::path(base) / user / project;
        
        Feel::cout << tc::green 
                   << "Custom repository callback computed: " << result.string() 
                   << tc::reset << std::endl;
        
        return result;
    });

    //
    // Initialize Feel++ environment with custom repository
    // The callback will be invoked after options are parsed
    //
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout(),
                     _config=config );

    //
    // Now the repository is configured and we can use it
    //
    Feel::cout << "Repository root: " << Environment::rootRepository() << std::endl;
    Feel::cout << "App repository: " << Environment::appRepository() << std::endl;
    Feel::cout << "Logs directory: " << Environment::logsRepository() << std::endl;

    //
    // Example: Demonstrate that options affect the repository location
    //
    Feel::cout << "\nConfiguration used:" << std::endl;
    Feel::cout << "  storage.base = " << soption(_name="storage.base") << std::endl;
    Feel::cout << "  storage.project = " << soption(_name="storage.project") << std::endl;
    Feel::cout << "  storage.user = " << soption(_name="storage.user") << std::endl;

    //
    // Example: Change repository at runtime (if needed)
    //
    if ( Environment::vm().count("storage.project") )
    {
        std::string new_project = soption(_name="storage.project") + "-v2";
        Feel::cout << "\nChanging to variant project: " << new_project << std::endl;
        
        // We can reconfigure the repository
        Environment::changeRepository(
            _directory=boost::format("variant/%1%") % new_project
        );
        
        Feel::cout << "New app repository: " << Environment::appRepository() << std::endl;
    }

    return 0;
}

/*
 * Usage examples:
 *
 * 1. Use default options:
 *    ./example_custom_repository
 *
 * 2. Override storage location:
 *    ./example_custom_repository --storage.base=/scratch/myproject --storage.user=john
 *
 * 3. Specify project and base directory:
 *    ./example_custom_repository --storage.base=/work --storage.project=simulation1
 *
 * 4. Combine with standard Feel++ options:
 *    ./example_custom_repository --storage.base=/work --storage.project=sim1 --v=2
 *
 * 5. With config file:
 *    ./example_custom_repository --config-file=mysim.cfg
 *    (where mysim.cfg contains: storage.base=/data storage.project=test1)
 *
 * See also:
 *  - example_custom_repository_env.cpp: Environment variable based selection
 *  - example_custom_repository_cluster.cpp: Automatic cluster detection
 */
