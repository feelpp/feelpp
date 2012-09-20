/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-04

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file myapp.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-04
 */
//# marker1 #
#include <feel/feel.hpp>
//# endmarker1 #

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
    po::options_description myappoptions( "MyApp options" );
    myappoptions.add_options()
    ( "dt", po::value<double>()->default_value( 1 ), "time step value" )
    ;

    // return the options myappoptions and the feel_options defined
    // internally by Feel
    return myappoptions.add( backend_options( "myapp" ) );
}
//# endmarker2 #


/**
 * \class MyApp
 *
 * This is a demo class to illustrate what is done (at the very least)
 * in subclasses of Feel::Application
 *
 */
//# marker4 #
class MyApp: public Application
{
public:
    /**
     * This function is responsible for the actual work done by MyApp.
     */
    void run();
};
//# endmarker4 #

//# marker6 #
void MyApp::run()
{
    /**
     * store all subsequent data files in a HOME/feel/doc/tutorial/myapp/
     */
    /** \code */
    //# marker8 #
    Environment::changeRepository( boost::format( "doc/manual/tutorial/%1%/" )
                                   % this->about().appName() );
    //# endmarker8 #
    /** \endcode */

    /**
     * print some information that will be written in the log file in
     * HOME/feel/doc/tutorial/myapp/myapp-1.0
     */
    /** \code */
    LOG(INFO) << "the value of dt is " << Environment::vm()["dt"].as<double>() << "\n";
    LOG(INFO) << "the value of myapp-solver-type is " << Environment::vm()["myapp.ksp-type"].as<std::string>() << "\n";
    LOG(INFO) << "the value of myapp-pc-type is " << Environment::vm()["myapp.pc-type"].as<std::string>() << "\n";
    /** \endcode */
}


/**
 * main function: entry point of the program
 */
//# marker7 #
int main( int argc, char** argv )
{
    /**
     * Initialize Feel++ Environment
     */
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="myapp",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    /**
     * intantiate a MyApp class
     */
    /** \code */
    MyApp app;
    /** \endcode */

    /**
     * run the application
     */
    /** \code */
    app.run();
    /** \endcode */
}
//# endmarker7 #
