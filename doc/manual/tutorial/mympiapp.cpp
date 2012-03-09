/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-11

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
   \file mympiapp.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-11
 */
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/applicationmpi.hpp>

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description mympiappoptions("MyMpiApp options");
    mympiappoptions.add_options()
        ("dt", po::value<double>()->default_value( 1 ), "time step value")
        ;

    // return the options mympiappoptions and the feel_options defined
    // internally by Feel
    return mympiappoptions.add( feel_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "mympiapp" ,
                     "mympiapp" ,
                     "0.1",
                     "my first Feel application",
                     AboutData::License_GPL,
                     "Copyright (c) 2008-2012 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme",
                    "developer",
                    "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;
}


class MyMpiApp: public ApplicationMpi
{
public:

    /**
     * constructor only about data and no options description
     */
    MyMpiApp( int argc, char** argv, AboutData const& );

    /**
     * constructor about data and options description
     */
    MyMpiApp( int argc, char** argv,
           AboutData const&,
           po::options_description const&  );

    /**
     * must be redefined by ApplicationMpi subclass
     */
    void run();
};

MyMpiApp::MyMpiApp(int argc, char** argv,
             AboutData const& ad )
    :
    ApplicationMpi( argc, argv, ad )
{}
MyMpiApp::MyMpiApp(int argc, char** argv,
             AboutData const& ad,
             po::options_description const& od )
    :
    ApplicationMpi( argc, argv, ad, od )
{}
void MyMpiApp::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    this->changeRepository( boost::format( "%1%/" )
                            % this->about().appName() );

    Log() << "the value of dt is " << this->vm()["dt"].as<double>()
          << "\n";

    Log() << "we are on processor " << this->processorName() << "\n";
    Log() << "this is process number " << this->processId()
          << " out of " << this->nProcess() << "\n";
}

//
// main function: entry point of the program
//
int main( int argc, char** argv )
{
    MyMpiApp app( argc, argv, makeAbout(), makeOptions() );

    app.run();
}


