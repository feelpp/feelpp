/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-11

  Copyright (C) 2006 EPFL

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file levelset2d.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-11
 */
#include "levelset.hpp"

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description levelsetoptions( "LevelSet options" );
    levelsetoptions.add_options()
    ( "dt", Feel::po::value<double>()->default_value( 0.02 ),
      "time step value" )
    ( "ft", Feel::po::value<double>()->default_value( 1 ),
      "Final time value" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ),
      "first h value to start convergence" )
    ( "export", Feel::po::value<int>()->default_value( 0 ),
      "stride for result export (0=no export)" )
    ( "stabcoeff", Feel::po::value<double>()->default_value( 0.1 ),
      "interior penalty stabilization coefficient" )
    ( "theta", Feel::po::value<double>()->default_value( 0.5 ),
      "coefficient of theta-scheme in time (0=FE, 0.5=CN, 1=BE)" )
    ( "reinitevery", Feel::po::value<int>()->default_value( 0 ),
      "number of timesteps between reinitializations (0=adaptive)" )
    ( "mingrad", Feel::po::value<double>()->default_value( 0.2 ),
      "minimal gradient for adaptive reinitialization" )
    ( "bdf2", "use BDF2 time stepping" )
    ;

    Feel::po::options_description solveroptions( "algebraic solver options" );
    solveroptions.add_options()
    ( "solver", Feel::po::value<std::string>()->default_value( "gmres" ),
      "solver type (gmres, bicgstab)" )
    ( "tolerance", Feel::po::value<double>()->default_value( 2.e-10 ),
      "solver tolerance" )
    ( "verbose", Feel::po::value<int>()->default_value( 0 ),
      "(=0,1,2) print solver iterations" )
    ( "maxiter", Feel::po::value<int>()->default_value( 1000 ),
      "set maximum number of iterations" )
    ( "fillin", Feel::po::value<int>()->default_value( 2 ),
      "fill-in for incomplete factorizations" )
    ( "threshold", Feel::po::value<double>()->default_value( 1.e-3 ),
      "threshold for incomplete factorizations" )
    ;
    return levelsetoptions.add( solveroptions );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "levelset" ,
                           "levelset" ,
                           "0.1",
                           "2D and 3D Level Set Test Problem",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2006 EPFL" );

    about.addAuthor( "Christoph Winkelmann", "developer",
                     "christoph.winkelmann@epfl.ch", "" );
    return about;

}



int
main( int argc, char** argv )
{
    Feel::LevelSet levelset( argc, argv, makeAbout(), makeOptions() );
    levelset.run();
}
