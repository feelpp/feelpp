/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-21

  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)

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
   \file options.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-21
 */
#include <feel/options.hpp>
#include <feel/feelalg/enums.hpp>
#include <feel/feelfilters/gmshenums.hpp>

namespace Feel
{

std::string
prefixvm( std::string const& prefix,
          std::string const& opt,
          std::string const& sep );


po::options_description
file_options( std::string const& appname )
{
    po::options_description file( "File options" );
    file.add_options()
        ( "config-file", po::value<std::string>()->default_value(appname+".cfg"), "specify .cfg file" )
        ( "result-file", po::value<std::string>()->default_value(appname+".res"), "specify .res file" )
        ( "response-file", po::value<std::string>()->default_value(appname), "can be specified with '@name', too" )
        ;
    return file;
}

po::options_description
generic_options()
{
    po::options_description generic( "Generic options" );
    generic.add_options()
        ( "authors", "prints the authors list" )
        ( "copyright", "prints the copyright statement" )
        ( "help", "prints this help message" )
        ( "license", "prints the license text" )
        ( "version", "prints the version" )
        ( "v", po::value<int>()->default_value(0), "verbosity level" )
        ( "feelinfo", "prints feel libraries information" )
        ( "nochdir", "Don't change repository directory even though it is called" )
        ( "directory", po::value<std::string>(), "change directory to specified one" )
        ;
    return generic;
}
po::options_description
functions_options( std::string const& prefix )
{
    po::options_description _options( "Functions " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"functions.f" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "f" )
        ( prefixvm( prefix,"functions.g" ).c_str(), Feel::po::value<std::string>()->default_value( "0" ), "g" )
        ( prefixvm( prefix,"functions.alpha" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "alpha" )
        ( prefixvm( prefix,"functions.beta" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "beta" )
        ( prefixvm( prefix,"functions.beta_x" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "beta x" )
        ( prefixvm( prefix,"functions.beta_y" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "beta y" )
        ( prefixvm( prefix,"functions.beta_z" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "beta z" )
        ( prefixvm( prefix,"functions.epsilon" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "epsilon" )
        ( prefixvm( prefix,"functions.gamma" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "gamma" )
        ( prefixvm( prefix,"functions.delta" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "delta" )
        ( prefixvm( prefix,"functions.nu" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "nu" )
        ( prefixvm( prefix,"functions.mu" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "mu" )
        ( prefixvm( prefix,"functions.rho" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), "rho" )
        ;
    return _options;
}
po::options_description
parameters_options( std::string const& prefix )
{
    po::options_description _options( "Parameters " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"parameters.f" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "f" )
        ( prefixvm( prefix,"parameters.g" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "g" )
        ( prefixvm( prefix,"parameters.alpha" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "alpha" )
        ( prefixvm( prefix,"parameters.beta" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "beta" )
        ( prefixvm( prefix,"parameters.beta_x" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "beta x" )
        ( prefixvm( prefix,"parameters.beta_y" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "beta y" )
        ( prefixvm( prefix,"parameters.beta_z" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "beta z" )
        ( prefixvm( prefix,"parameters.epsilon" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "epsilon" )
        ( prefixvm( prefix,"parameters.gamma" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "gamma" )
        ( prefixvm( prefix,"parameters.delta" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "delta" )
        ( prefixvm( prefix,"parameters.nu" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "nu" )
        ( prefixvm( prefix,"parameters.mu" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "mu" )
        ( prefixvm( prefix,"parameters.rho" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "rho" )
        ;
    return _options;
}

po::options_description
gmsh_options( std::string const& prefix )
{
    po::options_description _options( "Gmsh " + prefix + " options" );
    _options.add_options()
    // solver options
        ( prefixvm( prefix,"gmsh.filename" ).c_str(), Feel::po::value<std::string>()->default_value( "untitled.geo" ), "Gmsh filename" )
        ( prefixvm( prefix,"gmsh.depends" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "list of files separated by , or ; that are dependencies of a loaded Gmsh geometry" )
        ( prefixvm( prefix,"gmsh.hsize" ).c_str(), Feel::po::value<double>()->default_value( 0.1 ), "default characteristic mesh size" )
        ( prefixvm( prefix,"gmsh.geo-variables-list" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "modify a list of geo variables (ex : alpha=1:beta=2)" )
        ( prefixvm( prefix,"gmsh.refine" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "refinement by splitting level" )
        ( prefixvm( prefix,"gmsh.straighten" ).c_str(), Feel::po::value<bool>()->default_value( true ), "straighten high order mesh" )
        ( prefixvm( prefix,"gmsh.structured" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "generated a structured mesh" )
        ( prefixvm( prefix,"gmsh.rebuild" ).c_str(), Feel::po::value<bool>()->default_value( true ), "force rebuild msh file from geo file" )
        ( prefixvm( prefix,"gmsh.physical_are_elementary_regions" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Physical regions are defined by elementary regions, useful for medit format" )
        ( prefixvm( prefix,"gmsh.partition" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Partition Gmsh mesh once generated or loaded" )
#if defined(HAVE_METIS)
        ( prefixvm( prefix,"gmsh.partitioner" ).c_str(), Feel::po::value<int>()->default_value( GMSH_PARTITIONER_DEFAULT ), "Gmsh partitioner (1=CHACO, 2=METIS)" )
#else
        ( prefixvm( prefix,"gmsh.partitioner" ).c_str(), Feel::po::value<int>()->default_value( GMSH_PARTITIONER_DEFAULT ), "Gmsh partitioner (1=CHACO)" )
#endif
        ( prefixvm( prefix,"gmsh.format" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "Gmsh file format (0=ASCII, 1=BINARY)" )
        ( prefixvm( prefix,"gmsh.substructuring" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "Generate substructuration (0=false, 1=true)" )
        // domain
        ( prefixvm( prefix,"gmsh.domain.dim" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "Gmsh domain dimension" )
        ( prefixvm( prefix,"gmsh.domain.xmin" ).c_str(), Feel::po::value<double>()->default_value( 0. ), "Gmsh domain xmin" )
        ( prefixvm( prefix,"gmsh.domain.xmax" ).c_str(), Feel::po::value<double>()->default_value( 1. ), "Gmsh domain xmax" )
        ( prefixvm( prefix,"gmsh.domain.ymin" ).c_str(), Feel::po::value<double>()->default_value( 0. ), "Gmsh domain ymin" )
        ( prefixvm( prefix,"gmsh.domain.ymax" ).c_str(), Feel::po::value<double>()->default_value( 1. ), "Gmsh domain ymax" )
        ( prefixvm( prefix,"gmsh.domain.zmin" ).c_str(), Feel::po::value<double>()->default_value( 0. ), "Gmsh domain zmin" )
        ( prefixvm( prefix,"gmsh.domain.zmax" ).c_str(), Feel::po::value<double>()->default_value( 1. ), "Gmsh domain zmax" )

        ( prefixvm( prefix,"gmsh.domain.usenames" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "Gmsh domain use string name for physicals" )

        ( prefixvm( prefix,"gmsh.domain.substructuring" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "Generate substructuration (0=false, 1=true)" );


    return _options;

}
po::options_description
gmsh_domain_options( std::string const& prefix )
{
    po::options_description _options( "Gmsh Domain " + prefix + " options" );
    _options.add_options()
    // solver options
        ( prefixvm( prefix,"gmsh.domain.shape" ).c_str(), Feel::po::value<std::string>()->default_value( "hypercube" ), "Domain shape" )
        ( prefixvm( prefix,"gmsh.domain.convex" ).c_str(), Feel::po::value<std::string>()->default_value( "Simplex" ), "Convex type for Domain mesh (Simplex or Hypercube)" )
        ( prefixvm( prefix,"gmsh.domain.shear" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "shear value for hypercube domain" )
        ( prefixvm( prefix,"gmsh.domain.recombine" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "recombine elements to generate hypercube" )
        ( prefixvm( prefix,"gmsh.domain.recombine" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "recombine elements to generate hypercube" )

        ( prefixvm( prefix,"gmsh.domain.substructuring" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "generate substructuring markers for hypercube domain" )
        ( prefixvm( prefix,"gmsh.domain.usenames" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "mark boundaries with names" )
        ( prefixvm( prefix,"gmsh.domain.addmidpoint" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "add mid point on geometrical edges" )


        ( prefixvm( prefix,"gmsh.domain.xmin" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "minimum value in x-direction" )
        ( prefixvm( prefix,"gmsh.domain.ymin" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "minimum value in y-direction" )
        ( prefixvm( prefix,"gmsh.domain.zmin" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "minimum value in z-direction" )
        ( prefixvm( prefix,"gmsh.domain.xmax" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "maximum value in x-direction" )
        ( prefixvm( prefix,"gmsh.domain.ymax" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "maximum value in y-direction" )
        ( prefixvm( prefix,"gmsh.domain.zmax" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "maximum value in z-direction" )
        ( prefixvm( prefix,"gmsh.domain.nx" ).c_str(), Feel::po::value<double>()->default_value( 10 ), "number of subdivison in in x-direction" )
        ( prefixvm( prefix,"gmsh.domain.ny" ).c_str(), Feel::po::value<double>()->default_value( 10 ), "number of subdivison in in y-direction" )
        ;


    return _options;

}

po::options_description
on_options( std::string const& prefix )
{
    po::options_description _options( "Dirichlet treatment options " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"on.type" ).c_str(), Feel::po::value<int>()->default_value( ON_ELIMINATION ), "Strong Dirichlet conditions treatment type" )
        ( prefixvm( prefix,"on.verbose" ).c_str(), Feel::po::value<bool>()->default_value( false ), "print in logfiles information about Dirichlet conditions treatment" )
        ;
    return _options;
}

po::options_description
parallel_options( std::string const& prefix )
{
    po::options_description _options( "Parallel " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"parallel.cpu.enable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Enable the use of additional cores for parallelization" )
        ( prefixvm( prefix,"parallel.cpu.impl" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "Specify the implementation for multithreading" )
        ( prefixvm( prefix,"parallel.cpu.restrict" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "Restrict the multithreading to N additional cores per MPI process (0: guess the maximum number of usable cores)" )
        ( prefixvm( prefix,"parallel.gpu.enable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Enable the use of GPU for parallelization" )
        ;
    return _options;
}

po::options_description
ginac_options( std::string const& prefix )
{
    po::options_description _options( "GiNaC " + prefix + " options" );
    _options.add_options()
    // solver options
        ( prefixvm( prefix,"ginac.strict-parser" ).c_str(), Feel::po::value<bool>()->default_value( false ), "enable strict parsing of GiNaC expressions, no extra variables/symbols can be defined if set to true" )
        ;
    return _options;
}

po::options_description
error_options( std::string const& prefix )
{
    po::options_description _options( "Error options (" + prefix + ")" );
    _options.add_options()
    // error options
        ( prefixvm( prefix, "error.exact" ).c_str(), Feel::po::value<std::string>()->default_value(""), "exact solution" )
        ( prefixvm( prefix, "error.params" ).c_str(), Feel::po::value<std::string>()->default_value(""), "exact solution parameters" )
        ( prefixvm( prefix, "error.rhs" ).c_str(), Feel::po::value<std::string>()->default_value(""), "rhs" )
        ( prefixvm( prefix, "error.rhs.computed" ).c_str(), Feel::po::value<bool>()->default_value( false ), "rhs computed" )
        ( prefixvm( prefix, "error.convergence" ).c_str(), Feel::po::value<bool>()->default_value( false ), "convergence" )
        ( prefixvm( prefix, "error.convergence.steps" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "number of convergence steps" )
    ;
    return _options;
}

po::options_description
feel_options( std::string const& prefix  )
{
    auto opt = benchmark_options( prefix )
        .add( mesh_options( 1, prefix ) )
        .add( mesh_options( 2, prefix ) )
        .add( mesh_options( 3, prefix ) )
        /* alg options */
        .add( backend_options() )
#if defined(FEELPP_HAS_PETSC_H)
        .add( backendpetsc_options( prefix ) )
#endif
        .add( solvereigen_options( prefix ) )
#if defined( FEELPP_HAS_TRILINOS_EPETRA )
        .add( backendtrilinos_options( prefix ) )
#endif
        /* nonlinear solver options */
        .add( nlsolver_options() )

        /* discr options */
        .add( bdf_options( prefix ) )

        /* exporter options */
        .add( exporter_options( prefix ) )

        /* gmsh options */
        .add( gmsh_options( prefix ) )

        /* gmsh domain options */
        .add( gmsh_domain_options( prefix ) )
        #
#if defined(FEELPP_HAS_HARTS)
        .add( parallel_options( prefix ) )
#endif

        /* ginac options */
        .add( ginac_options( prefix ) )

        /* material options */
        .add( material_options( prefix ) )

        /* error options */
        .add( error_options( prefix ) )

        /* functions options */
        .add( functions_options( prefix ) )

        /* parameters options */
        .add( parameters_options( prefix ) )

        /* functions options */
        .add( on_options( prefix ) )

        ;

    if ( prefix.empty() )
        opt.add( generic_options() );
    return opt;

}
}
