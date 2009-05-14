/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-04-17

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file polyvis.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-04-17
 */
/** include predefined life command line options */
#include <life/options.hpp>

#include <life/lifecore/application.hpp>

#include <polyvisbase.hpp>

/** use Life namespace */
using namespace Life;


/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Life::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description polyvisoptions("Polyvis options");
    polyvisoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.1 ), "mesh size")
        ("dim", po::value<uint16_type>()->default_value( 2 ), "dimension")
        ("poly", po::value<std::string>()->default_value( "lagrange" ), "polynomial family")
        ("order", po::value<uint16_type>()->default_value( 2 ), "polynomial order")
        ("convex", po::value<std::string>()->default_value( "Simplex" ), "Convex type (Simplex, SimplexProduct")
        ;
    return polyvisoptions.add( Life::life_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Life::Application subclass.
 *
 * \return some data about the application.
 */
inline
AboutData
makeAbout()
{
    AboutData about( "polyvis" ,
                     "polyvis" ,
                     "0.2",
                     "nD(n=1,2,3) Polynomial on simplices or simplex products",
                     Life::AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


/**
 * \class Polyvis
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

    /**
     * Constructor
     */
    PolyvisApp( int argc, char** argv,
                AboutData const& ad,
                po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        polyvis( PolyvisBase::New( this->vm() ) )
    {
        std::cout << "[PolyvisApp] init\n";
    }

    /**
     * run the application
     */
    void run();

private:
    polyvisbase_ptrtype polyvis;
}; // PolyvisApp




void
PolyvisApp::run()
{
    /**
     * print help if --help is passed to the command line
     */
    /** \code */
    if ( this->vm().count( "help" ) )
        {
            std::cout << "[PolyvisApp::run] help\n";
            std::cout << this->optionsDescription() << "\n";
            return;
        }
    /** \endcode */

    std::cout << "[PolyvisApp::run] changeRepo\n";
    /**
     * we change to the directory where the results and logs will be
     * stored
     */
    /** \code */
    this->changeRepository( boost::format( "%1%/%2%/h_%3%/" )
                            % this->about().appName()
                            % polyvis->name()
                            % this->vm()["hsize"].as<double>()
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
    //try {
    /**
     * intantiate a Polyvis<Dim> class with Dim=2 (e.g. geometric dimension is 2)
     */
    /** \code */
    PolyvisApp polyvis( argc, argv, makeAbout(), makeOptions() );
    /** \encode */

    /**
     * run the application
     */
    /** \code */
    polyvis.run();
    /** \endcode */
#if 0
    }
    catch( std::exception const& e )
        {
            std::cout << "Caught exception " << e.what() << std::endl;
        }
    catch( ... )
        {
            std::cout << "Caught an unknown exception " << std::endl;
        }
#endif
}






