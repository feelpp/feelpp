/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-06-12

  Copyright (C) 2008 Universite Joseph Fourier (Grenoble I)
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
   \file test_bdf.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-06-12
 */
#include <feel/feelcore/feel.hpp>
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelts/tsbase.hpp>

Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_bdf" ,
                           "test_bdf" ,
                           "0.1",
                           "Bdf test",
                           Feel::AboutData::License_LGPL,
                           "Copyright (c) 2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}
inline
Feel::po::options_description
makeOptions()
{
    return Feel::feel_options().add( Feel::bdf_options( "test_bdf" ) );
}

using namespace Feel;
class MyApp: public Application
{
    typedef Application super;
public:
    MyApp()
        :
        super(),
        bdf( this->vm(), "bdf","test_bdf", Environment::worldComm() )
    {}
    void run()
    {
        for ( bdf.start(); bdf.isFinished() == false; bdf.next() )
        {
            bdf.shiftRight();
            TSBaseMetadata meta( bdf );
            meta.save();
            meta.load();
        }

    }
private:
    TSBase bdf;
};
int main( int argc, char** argv )
{
    Feel::Environment env( _argc=argc,_argv=argv,_desc=makeOptions(),_about=makeAbout() );
    MyApp myapp;
    myapp.run();



}
