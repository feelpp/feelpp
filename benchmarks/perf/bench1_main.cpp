/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-06-14

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file bench1_main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-06-14
 */
#include <bench1.hpp>

namespace Life
{
AboutData
makeAbout()
{
    AboutData about( "bench1" ,
                     "bench1" ,
                     "0.2",
                     "assembly performance",
                     AboutData::License_LGPL,
                     "Copyright (c) 2005,2006 EPFL"
                     "Copyright (c) 2006-2010 Université Joseph Fourier (Grenoble 1)"
        );

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}

po::options_description
makeOptions()
{
    po::options_description desc("Specific options");
    desc.add_options()
        ("dim", po::value<int>()->default_value( 1 ), "dimension (1,2,3)")
        ("hsize", po::value<double>()->default_value( 0.1 ), "element size")
        ("shape", po::value<std::string>()->default_value( "simplex" ), "type of domain shape: simplex, hypercube ellipsoid")
        ;
    return desc.add( life_options() );
}


}
int
main( int argc, char** argv )
{
    Life::Bench1 bench1(argc, argv, Life::makeAbout(), Life::makeOptions() );

    bench1.run();
}

