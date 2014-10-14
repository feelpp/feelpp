/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-17

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file polyvis.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-17
 */
/** include predefined feel command line options */
#include <feel/options.hpp>

#include <feel/feelcore/application.hpp>

#include <polyvisbase.hpp>

/** use Feel namespace */
using namespace Feel;


/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description polyvisoptions( "Polyvis options" );
    polyvisoptions.add_options()
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ( "dim", po::value<uint16_type>()->default_value( 2 ), "dimension" )
        ( "poly", po::value<std::string>()->default_value( "lagrange" ), "polynomial family" )
        ( "order", po::value<uint16_type>()->default_value( 2 ), "polynomial order" )
        ( "convex", po::value<std::string>()->default_value( "Simplex" ), "Convex type (Simplex, Hypercube" )

        ( "meshes-2d", po::value< std::string >(), "mesh name" )
        ( "meshes-3d", po::value< std::string >(), "mesh name" )


        ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ;
    return polyvisoptions;
}


/**
 * @brief Generates data to plot polynomial sets over a simplex or hypercube
 *
 */
class PolyvisApp
    :
    public Application
{
    typedef Application super;

public:

    typedef PolyvisBase polyvisbase_type;
    typedef boost::shared_ptr<PolyvisBase> polyvisbase_ptrtype;
    PolyvisApp()
        :
        super(),
        polyvis( PolyvisBase::New( Environment::vm() ) )
        {
            std::cout << "[PolyvisApp] init\n";
        }
    void run();

private:
    polyvisbase_ptrtype polyvis;
}; // PolyvisApp




void
PolyvisApp::run()
{
    std::cout << "[PolyvisApp::run] changeRepo\n";
    /**
     * we change to the directory where the results and logs will be
     * stored
     */
    /** \code */
    this->changeRepository( boost::format( "%1%/%2%/h_%3%/" )
                            % this->about().appName()
                            % polyvis->name()
                            % Environment::vm()["hsize"].as<double>()
                          );
    std::cout << "[PolyvisApp::run] run\n";
    polyvis->run();

    /** \endcode */
} // Polyvis::run

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about( _name="polyvis",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org") );
    PolyvisApp polyvis;
    polyvis.run();
}
