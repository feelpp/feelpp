/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-11

  Copyright (C) 2006 EPFL

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file levelset2d.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-11
 */
#include "levelset.hpp"

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description levelsetoptions("LevelSet options");
    levelsetoptions.add_options()
        ("dt", Life::po::value<double>()->default_value( 0.02 ),
         "time step value")
        ("ft", Life::po::value<double>()->default_value( 1 ),
         "Final time value")
        ("hsize", Life::po::value<double>()->default_value( 0.1 ),
         "first h value to start convergence")
        ("export", Life::po::value<int>()->default_value( 0 ),
         "stride for result export (0=no export)")
        ("stabcoeff", Life::po::value<double>()->default_value( 0.1 ),
         "interior penalty stabilization coefficient")
        ("theta", Life::po::value<double>()->default_value( 0.5 ),
         "coefficient of theta-scheme in time (0=FE, 0.5=CN, 1=BE)")
        ("reinitevery", Life::po::value<int>()->default_value( 0 ),
         "number of timesteps between reinitializations (0=adaptive)")
        ("mingrad", Life::po::value<double>()->default_value( 0.2 ),
         "minimal gradient for adaptive reinitialization")
        ("bdf2", "use BDF2 time stepping")
        ;

    Life::po::options_description solveroptions("algebraic solver options");
    solveroptions.add_options()
        ("solver", Life::po::value<std::string>()->default_value( "gmres" ),
         "solver type (gmres, bicgstab)")
        ("tolerance", Life::po::value<double>()->default_value( 2.e-10 ),
         "solver tolerance")
        ("verbose", Life::po::value<int>()->default_value( 0 ),
         "(=0,1,2) print solver iterations")
        ("maxiter", Life::po::value<int>()->default_value( 1000 ),
         "set maximum number of iterations")
        ("fillin", Life::po::value<int>()->default_value( 2 ),
         "fill-in for incomplete factorizations")
        ("threshold", Life::po::value<double>()->default_value( 1.e-3 ),
         "threshold for incomplete factorizations")
        ;
    return levelsetoptions.add( solveroptions );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "levelset" ,
                            "levelset" ,
                            "0.1",
                            "2D and 3D Level Set Test Problem",
                            Life::AboutData::License_GPL,
                            "Copyright (c) 2006 EPFL");

    about.addAuthor("Christoph Winkelmann", "developer",
                    "christoph.winkelmann@epfl.ch", "");
    return about;

}



int
main( int argc, char** argv )
{
    Life::LevelSet levelset( argc, argv, makeAbout(), makeOptions());
    levelset.run();
}
