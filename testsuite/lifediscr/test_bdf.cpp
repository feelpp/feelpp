/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-06-12

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-06-12
 */
#include <life/lifecore/life.hpp>
#include <life/options.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifediscr/bdf2.hpp>

Life::AboutData
makeAbout()
{
    Life::AboutData about( "test_bdf" ,
                           "test_bdf" ,
                           "0.1",
                           "Bdf test",
                           Life::AboutData::License_LGPL,
                           "Copyright (c) 2008 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}
inline
Life::po::options_description
makeOptions()
{
    return Life::life_options();
}

using namespace Life;
class MyApp: public Application
{
    typedef Application super;
public:
    MyApp( int argc, char** argv, Life::AboutData const& about, po::options_description const& od )
        :
        super( argc, argv, about, od ),
        bdf( this->vm(), "bdf" )
    {}
    void run()
    {
        for( bdf.start(); bdf.isFinished() == false; bdf.next() )
            {
                bdf.shiftRight();
                BdfBaseMetadata meta( bdf );
                meta.save();
                meta.load();
            }

    }
private:
    BdfBase bdf;
};
int main( int argc, char** argv )
{
    MyApp myapp( argc, argv, makeAbout(), makeOptions() );
    myapp.run();



}
