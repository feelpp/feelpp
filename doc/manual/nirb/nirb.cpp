/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-02-01

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file nirb.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-02-01
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

//#include <feel/feelalg/solvereigen.hpp>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <vector>
#include <algorithm>

#include <boost/range/algorithm/max_element.hpp>
#include "nirb.hpp"

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline po::options_description makeOptions()
{
    po::options_description NIRBoptions( "Nirb options" );
    NIRBoptions.add_options()
    // meshes parameters
    ( "hfinsize", po::value<double>()->default_value( 0.05 ), "fine mesh size" )
    ( "hcoarsesize", po::value<double>()->default_value( 0.1 ), "coarse mesh size" )
    ( "ReadingMeshes",po::value<int>()->default_value( 0 ),"Reading meshes in file if set to 1 -- default value = 0" )
    // Reduced basis parameters
    ( "NbSnapshot",po::value<int>()->default_value( 100 ),"numbers of snapshot computed" )
    ( "sizeRB",po::value<int>()->default_value( 15 ),"size of reduced basis" )

    ( "muMin", po::value<double>()->default_value( 0 ), "angle in [0,pi/2]" )
    ( "muMax", po::value<double>()->default_value( M_PI/2. ),"angle in [0,pi/2]" )

    ( "mu",po::value<double>()->default_value( 1. ),"angle in [0,pi/2]" )

    ( "Sampling",po::value<int>()->default_value( 1 ),"Does not compute sampling if set to 0, (set to 1 by default)" )
    ( "SamplingCoarse",po::value<int>()->default_value( 0 ),"Compute Coarse sampling if set to 1, (set to 0 by default)" )

    ( "Offline",po::value<int>()->default_value( 1 ),"integer equal to 0 if the offline has not  to be done" )

    ( "ComputeError",po::value<int>()->default_value( 1 ),"integer equal to 0 if the error computation has not to be done" )

    ( "polynomialOrder",po::value<int>()->default_value( 3 ),"polynomial order" );
    ;
    return NIRBoptions.add( Feel::feel_options() );
}

//-----------------------------------------
//-----------------------------------------
/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline AboutData makeAbout()
{
    AboutData about( "nirb-test" ,
                     "nirb-test" ,
                     "0.2",
                     "Non intrusive reduced basis",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


//-----------------------------------------
//-----------------------------------------
/**
 * main function: entry point of the program
 */
int main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    Application app;

    /**
     * register the simgets
     */

    const int polynomialOrder =  app.vm()["polynomialOrder"].as<int>();


    if ( polynomialOrder == 1 )
    {
        app.add( new NIRBTEST<1>() );
    }

    else if ( polynomialOrder == 2 )
    {
        app.add( new NIRBTEST<2>() );
    }

    else if ( polynomialOrder == 3 )
    {
        app.add( new NIRBTEST<3>() );
    }

    else
    {
        throw std::logic_error( "Error with polynomialOrder variable, this application allows only P1 P2 and P3" );
    }

    /**
     * run the application
     */
    app.run();
}






