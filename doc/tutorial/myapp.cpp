/* -*- mode: c++ coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \file application.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-04
 */
#include <life/options.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>

using namespace Life;

inline
po::options_description
makeOptions()
{
    po::options_description myappoptions("MyApp options");
    myappoptions.add_options()
        ("dt", po::value<double>()->default_value( 1 ), "time step value")
        ;

    // return the options myappoptions and the life_options defined
    // internally by Life
    return myappoptions.add( life_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "myapp" ,
                     "myapp" ,
                     "0.1",
                     "my first Life application",
                     AboutData::License_GPL,
                     "Copyright (c) 2008 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme",
                    "developer",
                    "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;
}


class MyApp: public Application
{
public:

    /**
     * constructor only about data and no options description
     */
    MyApp( int argc, char** argv, AboutData const& );

    /**
     * constructor about data and options description
     */
    MyApp( int argc, char** argv,
           AboutData const&,
           po::options_description const&  );

    /**
     * must be redefined by Application subclass
     */
    void run();
};

MyApp::MyApp(int argc, char** argv,
             AboutData const& ad )
    :
    Application( argc, argv, ad )
{}
MyApp::MyApp(int argc, char** argv,
             AboutData const& ad,
             po::options_description const& od )
    :
    Application( argc, argv, ad, od )
{}
void MyApp::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }
    this->changeRepository( boost::format( "doc/tutorial/%1%/" )
                            % this->about().appName() );

    Log() << "the value of dt is " << this->vm()["dt"].as<double>() << "\n";
}

//
// main function: entry point of the program
//
int main( int argc, char** argv )
{
    MyApp app( argc, argv, makeAbout(), makeOptions() );

    app.run();
}

