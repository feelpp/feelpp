/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2011-11-18

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file mcs_example.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2011-11-18
 */

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/backend.hpp>
using namespace Feel;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
//# marker2 #
inline
po::options_description
makeOptions()
{
    po::options_description myappoptions( "MyMCSApp options" );
    myappoptions.add_options()
    ( "dt", po::value<double>()->default_value( 1 ), "time step value" )
    ;

    // return the options myappoptions and the feel_options defined
    // internally by Feel
    return myappoptions.add( feel_options() ).add( backend_options( "myapp" ) );
}
//# endmarker2 #


/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
//# marker3 #
inline
AboutData
makeAbout()
{
    AboutData about( "myapp" ,
                     "myapp" ,
                     "0.1",
                     "my first Feel application",
                     AboutData::License_GPL,
                     "Copyright (c) 2008 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme",
                     "developer",
                     "christophe.prudhomme@ujf-grenoble.fr", "" );
    return about;
}
//# endmarker3 #

/**
 * \class MyMCSApp
 *
 * This is a demo class to illustrate what is done (at the very least)
 * in subclasses of Feel::Application
 *
 */
//# marker4 #
class MyMCSApp: public Application
{
public:

    /**
     * constructor only about data and no options description
     */
    MyMCSApp( int argc, char** argv, AboutData const& );

    /**
     * constructor about data and options description
     */
    MyMCSApp( int argc, char** argv,
              AboutData const&,
              po::options_description const&  );

    /**
     * This function is responsible for the actual work done by MyMCSApp.
     */
    void run();
};
//# endmarker4 #

//# marker5 #
MyMCSApp::MyMCSApp( int argc, char** argv,
                    AboutData const& ad )
    :
    Application( argc, argv, ad )
{}
MyMCSApp::MyMCSApp( int argc, char** argv,
                    AboutData const& ad,
                    po::options_description const& od )
    :
    Application( argc, argv, ad, od )
{}
//# endmarker5 #

//# marker6 #
void MyMCSApp::run()
{
    /**
     * print the help if --help is passed as an argument
     */
    /** \code */
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    /** \endcode */
    //# endmarker6 #
    /**
     * store all subsequent data files in a HOME/feel/doc/tutorial/myapp/
     */
    /** \code */
    //# marker8 #
    this->changeRepository( boost::format( "doc/tutorial/%1%/" )
                            % this->about().appName() );
    //# endmarker8 #
    /** \endcode */

    /**
     * print some information that will be written in the log file in
     * HOME/feel/doc/tutorial/myapp/myapp-1.0
     */
    /** \code */
    Log() << "the value of dt is " << this->vm()["dt"].as<double>() << "\n";
    Log() << "the value of myapp-solver-type is " << this->vm()["myapp.ksp-type"].as<std::string>() << "\n";
    Log() << "the value of myapp-pc-type is " << this->vm()["myapp.pc-type"].as<std::string>() << "\n";
    /** \endcode */
}


/**
 * main function: entry point of the program
 */
//# marker7 #
int main( int argc, char** argv )
{
    Feel::Environment env( argc, argv );

    /**
     * intantiate a MyMCSApp class
     */
    /** \code */
    MyMCSApp app( argc, argv, makeAbout(), makeOptions() );
    /** \endcode */

    /**
     * run the application
     */
    /** \code */
    app.run();
    /** \endcode */
}
//# endmarker7 #

