/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2011-06-04

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
   \file thermalblock.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2011-06-04
 */
#include <thermalblock.hpp>


namespace Feel
{

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
po::options_description
makeThermalBlockOptions()
{
    po::options_description thermalblockoptions("ThermalBlock options");
    thermalblockoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.05 ), "mesh size")
        ("nx", po::value<int>()->default_value( 3 ), "number of blocks in the x direction")
        ("ny", po::value<int>()->default_value( 3 ), "number of blocks in the y direction")
        ("gamma_dir", Feel::po::value<double>()->default_value( 10 ),
         "penalisation parameter for the weak boundary Dirichlet formulation")
        ;
    return thermalblockoptions.add( Feel::feel_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
AboutData
makeThermalBlockAbout()
{
    AboutData about( "thermalblock" ,
                     "thermalblock" ,
                     "0.1",
                     "2D Heterogeneous Thermal Block Problem",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier");

    about.addAuthor("Abdoulaye Samake", "main developer", "abdoulaye.samake@ujf-grenoble.fr", "");
    about.addAuthor("Christophe Prud'homme", "contributor", "christophe.prudhomme@ujf-grenoble.fr", "");
    about.addAuthor("Stephane Veys", "contributor", "stephane.veys@imag.fr", "");
    return about;

}


}
