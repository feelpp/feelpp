/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2006-02-01

  Copyright (C) 2006 EPFL

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2006-02-01
 */
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <iostream>

#include <boost/foreach.hpp>
#include <boost/plugin.hpp>

#include <life/lifecore/application.hpp>
#include <life/lifepoly/lagrange.hpp>
#include <life/lifepoly/polynomialset.hpp>



#include <polyvis.hpp>


Life::AboutData
makeAbout()
{
    Life::AboutData about( "life_polyvis" ,
                            "life_polyvis" ,
                            "0.1",
                            "Generate ensight files to visualize polynomials",
                            Life::AboutData::License_LGPL,
                            "Copyright (c) 2006 EPFL");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@epfl.ch", "");
    return about;

}

class PolyvisApp : public Life::Application
{
public:
    typedef Life::Application super;

    PolyvisApp( int argc,  char** argv, Life::AboutData const& ad, Life::po::options_description const& od )
        :
        super( argc, argv, ad, od )
    {}

    void run() ;
};
void
PolyvisApp::run()
{
    try {

        std::string poly = vm()["poly"].as<std::string>();
        std::string lib = vm()["lib"].as<std::string>();
#if 1
        /* get the handle of the library */
        boost::plugin::dll d ( lib );

        boost::plugin::plugin_factory <Life::Polyvis> pf (d);


        std::cout << "*** Creating an instance of plugin class " << poly << " out of lib " << lib << "\n";
        std::auto_ptr <Life::Polyvis> p (pf.create ( poly,
                                                     "Lagrange",
                                                     vm()["dim"].as<int>(),
                                                     vm()["order"].as<int>()
                                                     ));
#else

#endif
        std::cout << "*** Calling method of the created instance\n";
    }
    catch ( std::logic_error const & e )
        {
            /* report error, and skip the library */
            std::cerr << "Could not load polynomial family: " << e.what () << std::endl;
        }

}

int
main( int argc, char** argv )
{
    Life::po::options_description desc("Specific options");
    desc.add_options()
        ("lib", Life::po::value<std::string>()->default_value("stdpoly.so"), "library of polynomials")
        //("list", "list components in pluging")
        ("poly,p", Life::po::value<std::string>()->default_value("Lagrange"), "polynomials to display")
        ("dim,d", Life::po::value<int>()->default_value(2), "dimension")
        ("order,o", Life::po::value<int>()->default_value(1), "polynomials order")
        ("convex,c", Life::po::value<int>()->default_value(1), "convex")
        ;

    PolyvisApp app( argc, argv, makeAbout(), desc );

    app.run();



}


