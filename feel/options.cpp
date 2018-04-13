/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
#include <boost/algorithm/string.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/enums.hpp>
#include <feel/feelcrb/crbenums.hpp>
#include <feel/feelfilters/gmshenums.hpp>

#if defined(FEELPP_HAS_HARTS)
#include "hartsconfig.h"
#include "HARTS.h"
#endif

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
        ( "config-files", po::value<std::vector<std::string> >()->multitoken(), "specify a list of .cfg file" )
        ( "bc-file", po::value<std::string>()->default_value(appname+".bc"), "specify boundary condition (.bc) file" )
        ( "mod-file", po::value<std::string>()->default_value(appname+".mod"), "specify model (.mod) file" )
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
        ( "help", "prints this help message associated with the Feel++ application options" )
        ( "help-lib", "prints the help message associated with the Feel++ library options" )
        ( "license", "prints the license text" )
        ( "version", "prints the version" )
        ( "v", po::value<int>()->default_value(0), "Show all VLOG(m) messages for m <= this."
          " Overridable by --vmodule." )
        ( "vmodule", po::value<std::string>()->default_value(""), "per-module verbose level."
          " Argument is a comma-separated list of <module name>=<log level>."
          " <module name> is a glob pattern, matched against the filename base"
          " (that is, name ignoring .cc/.h./-inl.h)."
          " <log level> overrides any value given by --v." )
        ( "feelinfo", "prints feel libraries information" )
        ( "nochdir", "Don't change repository directory even though it is called" )
        ( "rmlogs", "remove logs after execution" )
        ( "rm", "remove application repository after execution" )
        ( "directory", po::value<std::string>(), "change directory to specified one" )
        ( "npdir", po::value<bool>()->default_value(true), "enable/disable sub-directory np_<number of processors>")
        ( "fail-on-unknown-option", po::value<bool>()->default_value(false), "exit feel++ application if unknown option found" )
        ( "show-preconditioner-options", "show on the fly the preconditioner options used" )
        ( "serialization-library", po::value<std::string>()->default_value("boost"), "Library used for serialization" )
        ( "display-stats", po::value<bool>()->default_value(false), "display statistics (timers, iterations counts...)" )
        ( "subdir.expr", po::value<std::string>()->default_value("exprs"), "subdirectory for expressions" )
        ;
    return generic;
}

po::options_description
functionspace_options( std::string const& prefix )
{
    po::options_description _options( "Function Space options" );
    _options.add_options()
        ( prefixvm( prefix, "connect").c_str(), Feel::po::value<bool>()->default_value(false), "Update dof when MESH_CHANGE_COORD ?" )
        ;
    return _options;
}

po::options_description
onelab_options( std::string const& prefix )
{
    po::options_description onelab( "Onelab options" );
    onelab.add_options()
        ( prefixvm( prefix, "onelab.enable" ).c_str(), Feel::po::value<int>()->default_value(0), "Generate OneLab files for interaction with Gmsh" )
        ( prefixvm( prefix, "onelab.remote" ).c_str(), Feel::po::value<std::string>()->default_value("localhost"), "Remote host for Onelab interface" )
        ( prefixvm( prefix, "onelab.chroot" ).c_str(), Feel::po::value<std::string>()->default_value(""), "Chroot to use on remote host" )
        ( prefixvm( prefix, "onelab.np" ).c_str(), Feel::po::value<int>()->default_value(1), "Number of MPI processes to use" )
        ( prefixvm( prefix, "onelab.sync.script" ).c_str(), Feel::po::value<std::string>()->default_value(""), "Script used for syncing data" )
        ;
    return onelab;
}

po::options_description
functions_options( std::string const& prefix )
{
    po::options_description _options( "Functions " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"x" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "x coordinate value " )
        ( prefixvm( prefix,"y" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "y coordinate value " )
        ( prefixvm( prefix,"z" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "z coordinate value " );
    std::vector<std::string> alphabet { "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega" };
    for(std::string& name : alphabet) {
        _options.add_options() ( prefixvm( prefix,"functions." + name ).c_str(), Feel::po::value<std::string>()->default_value( name.size() == 1 ? "0" : "1" ), name.c_str() );
        if(name == "beta") {
            _options.add_options()
            ( prefixvm( prefix,"functions." + name + "_x" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), std::string(name + " x").c_str() )
            ( prefixvm( prefix,"functions." + name + "_y" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), std::string(name + " y").c_str() )
            ( prefixvm( prefix,"functions." + name + "_z" ).c_str(), Feel::po::value<std::string>()->default_value( "1" ), std::string(name + " z").c_str() );
        }
    }
    return _options;
}
po::options_description
parameters_options( std::string const& prefix )
{
    po::options_description _options( "Parameters " + prefix + " options" );
    std::vector<std::string> alphabet { "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega" };
    for(std::string& name : alphabet) {
        _options.add_options() ( prefixvm( prefix,"parameters." + name ).c_str(), Feel::po::value<double>()->default_value( 0 ), name.c_str() );
        if(name == "beta") {
            _options.add_options()
            ( prefixvm( prefix,"parameters." + name + "_x" ).c_str(), Feel::po::value<double>()->default_value( 0 ), std::string(name + " x").c_str() )
            ( prefixvm( prefix,"parameters." + name + "_y" ).c_str(), Feel::po::value<double>()->default_value( 0 ), std::string(name + " y").c_str() )
            ( prefixvm( prefix,"parameters." + name + "_z" ).c_str(), Feel::po::value<double>()->default_value( 0 ), std::string(name + " z").c_str() );
        }
    }
    return _options;
}

po::options_description
nlopt_options( std::string const& prefix )
{
    po::options_description _options( "NLopt " + prefix + " options" );
    _options.add_options()
    // solver options
        ( prefixvm( prefix,"nlopt.ftol_rel" ).c_str(), Feel::po::value<double>()->default_value( 1e-4 ), "NLopt objective function relative tolerance" )
        ( prefixvm( prefix,"nlopt.ftol_abs" ).c_str(), Feel::po::value<double>()->default_value( 1e-10 ), "NLopt objective function  absolute tolerance" )
        ( prefixvm( prefix,"nlopt.xtol_rel" ).c_str(), Feel::po::value<double>()->default_value( 1e-4 ), "NLopt variables  relative tolerance" )
        ( prefixvm( prefix,"nlopt.xtol_abs" ).c_str(), Feel::po::value<double>()->default_value( 1e-10 ), "NLopt variables  absolute tolerance" )
        ( prefixvm( prefix,"nlopt.maxeval" ).c_str(), Feel::po::value<int>()->default_value( 30 ), "NLopt maximum number of evaluations" )
        ;
    return _options;
}

po::options_description
mesh_options( std::string const& prefix )
{
    po::options_description _options( "Mesh " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"mesh.filename").c_str(), po::value<std::string>()->default_value( "untitled.geo" ), "mesh filename" )
        ( prefixvm( prefix,"mesh.partition.enable").c_str(), po::value<bool>()->default_value(0), "partition the mesh using Feel++ partitioners" )
        ( prefixvm( prefix,"mesh.partition.size").c_str(), po::value<int>()->default_value(1), "number of partitions" )
        ( prefixvm( prefix,"mesh.partition.type").c_str(), po::value<std::string>()->default_value("metis"), "mesh partitioner: metis, (more to come)" )
        ( prefixvm( prefix,"mesh.save.enable" ).c_str(), Feel::po::value<bool>()->default_value( true ), "enable saving mesh to disk" )
        ( prefixvm( prefix,"mesh.save.formats" ).c_str(), Feel::po::value<std::vector<std::string>>()->default_value( {"json+h5","msh"} ), "format of the mesh: json+h5, msh" )
        ( prefixvm( prefix,"mesh.load.enable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "enable loading mesh from disk, overriding file name extension" );
        ( prefixvm( prefix,"mesh.load.format" ).c_str(), Feel::po::value<std::string>()->default_value( "json+h5" ), "file format to load: msh, json" );

    return _options;
}

po::options_description
gmsh_options( std::string const& prefix )
{
    po::options_description _options( "Gmsh " + prefix + " options" );

    _options.add_options()
    // gmsh options
        ( prefixvm( prefix,"gmsh.filename" ).c_str(), Feel::po::value<std::string>()->default_value( "untitled.geo" ), "Gmsh filename" )
        ( prefixvm( prefix,"gmsh.depends" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "list of files separated by , or ; that are dependencies of a loaded Gmsh geometry" )
        ( prefixvm( prefix,"gmsh.hsize" ).c_str(), Feel::po::value<double>()->default_value( 0.1 ), "default characteristic mesh size" )
        ( prefixvm( prefix,"gmsh.scale" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "scale the mesh after loading" )
        ( prefixvm( prefix,"gmsh.hsize2" ).c_str(), Feel::po::value<double>()->default_value( 0.1 ), "characteristic mesh size" )
        ( prefixvm( prefix,"gmsh.geo-variables-list" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "modify a list of geo variables (ex : alpha=1:beta=2)" )
        ( prefixvm( prefix,"gmsh.save" ).c_str(), Feel::po::value<bool>()->default_value( true ), "save msh file to disk once generated" )
        ( prefixvm( prefix,"gmsh.straighten" ).c_str(), Feel::po::value<bool>()->default_value( true ), "straighten high order mesh" )
        ( prefixvm( prefix,"gmsh.structured" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "generated a structured mesh" )
        ( prefixvm( prefix,"gmsh.rebuild" ).c_str(), Feel::po::value<bool>()->default_value( true ), "force rebuild msh file from geo file" )
        ( prefixvm( prefix,"gmsh.refine" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "refinement by splitting level" )
        ( prefixvm( prefix,"gmsh.physical_are_elementary_regions" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Physical regions are defined by elementary regions, useful for medit format" )
        ( prefixvm( prefix,"gmsh.savehdf5" ).c_str(), Feel::po::value<bool>()->default_value( true ), "save msh file to disk once generated in HDF5 format" )
        ( prefixvm( prefix,"gmsh.use-json" ).c_str(), Feel::po::value<bool>()->default_value( false ), "use json/hdf5 file if it exists, instead of the Gmsh files (geo or msh)" )
        ( prefixvm( prefix,"gmsh.partition" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Partition Gmsh mesh once generated or loaded" )
        ( prefixvm( prefix,"gmsh.respect_partition" ).c_str(), Feel::po::value<bool>()->default_value( false ), "true to respect paritioning when mesh is loaded, false to ensure that partition is within the number of processors" )
        ( prefixvm( prefix,"gmsh.npartitions" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "Number of partitions" )
        ( prefixvm( prefix,"gmsh.partitioner" ).c_str(), Feel::po::value<int>()->default_value( GMSH_PARTITIONER_DEFAULT ), "Gmsh partitioner (1=CHACO, 2=METIS)" )
        ( prefixvm( prefix,"gmsh.verbosity" ).c_str(), Feel::po::value<int>()->default_value( 5 ), "Gmsh verbosity level (0:silent except fatal errors, 1:+errors, 2:+warnings, 3:+direct, 4:+info except status bar, 5:normal, 99:debug)" )
        ( prefixvm( prefix,"gmsh.format" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "Gmsh file format (0=ASCII, 1=BINARY)" )

        ( prefixvm( prefix,"gmsh.in-memory" ).c_str(), Feel::po::value<bool>()->default_value( false ), "false to save on disk, true to read geometry directly from memory" )
        ( prefixvm( prefix,"gmsh.algo2d" ).c_str(), Feel::po::value<int>(), "used algorithm to mesh\n"
          "2D : \n"
          " MESHADAPT       1\n"
          " AUTO            2\n"
          " MESHADAPT_OLD   4\n"
          " DELAUNAY        5\n"
          " FRONTAL         6\n"
          " BAMG            7\n"
          " FRONTAL_QUAD    8\n"
          " PACK_PRLGRMS    9\n" )
        ( prefixvm( prefix,"gmsh.algo3d" ).c_str(), Feel::po::value<int>(), "used algorithm to mesh\n"
          "3D :\n"
          " DELAUNAY          1\n"
          " FRONTAL           4\n"
          " FRONTAL_DEL       5\n"
          " FRONTAL_HEX       6\n"
          " MMG3D             7\n"
          " RTREE             9" )
        ( prefixvm( prefix,"gmsh.randFactor" ).c_str(), Feel::po::value<double>()->default_value( -1 ), "Mesh.RandomFactor. -1 stand for default value" )

        ( prefixvm( prefix,"partition.linear" ).c_str(), Feel::po::value<bool>()->default_value( false ), "linear partitioning if true (false otherwise)" );

    return _options;


}
po::options_description
gmsh_domain_options( std::string const& prefix )
{
    po::options_description _options( "Gmsh Domain " + prefix + " options" );
#if defined( FEELPP_HAS_GMSH_H )
    _options.add_options()
    // solver options
        ( prefixvm( prefix,"gmsh.domain.shape" ).c_str(), Feel::po::value<std::string>()->default_value( "hypercube" ), "Domain shape" )
        ( prefixvm( prefix,"gmsh.domain.convex" ).c_str(), Feel::po::value<std::string>()->default_value( "Simplex" ), "Convex type for Domain mesh (Simplex or Hypercube)" )
        ( prefixvm( prefix,"gmsh.domain.shear" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "shear value for hypercube domain" )
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
        ( prefixvm( prefix,"gmsh.domain.nx" ).c_str(), Feel::po::value<double>()->default_value( 2 ), "number of subdivison in in x-direction" )
        ( prefixvm( prefix,"gmsh.domain.ny" ).c_str(), Feel::po::value<double>()->default_value( 2 ), "number of subdivison in in y-direction" )
        ( prefixvm( prefix,"gmsh.domain.nz" ).c_str(), Feel::po::value<double>()->default_value( 2 ), "number of subdivison in in z-direction" )
        ;

#endif // FEELPP_HAS_GMSH_H
    return _options;

}


po::options_description
arm_options( std::string const& prefix )
{
    po::options_description _options( "Acusim Raw Mesh (ARM) " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"arm.filename" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "ARM filename" )
        ;
    return _options;

}

po::options_description
on_options( std::string const& prefix )
{
    po::options_description _options( "Dirichlet treatment options " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"on.type" ).c_str(), Feel::po::value<std::string>()->default_value( "elimination" ), "Strong Dirichlet conditions treatment type: elimination, elimination_keep_diagonal, elimination_symmetric, elimination_symmetric_keep_diagonal, penalisation" )
        ( prefixvm( prefix,"on.verbose" ).c_str(), Feel::po::value<bool>()->default_value( false ), "print in logfiles information about Dirichlet conditions treatment" )
        ( prefixvm( prefix,"on.value_on_diagonal" ).c_str(), Feel::po::value<double>()->default_value( 1.0 ), "value on diagonal of operator when eliminating" )
        ;
    return _options;
}

po::options_description
parallel_options( std::string const& prefix )
{
    po::options_description _options( "Parallel " + prefix + " options" );
    _options.add_options()
#if defined(FEELPP_HAS_HARTS)
        ( prefixvm( prefix,"parallel.cpu.enable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Enable the use of additional cores for parallelization" )
        ( prefixvm( prefix,"parallel.cpu.impl" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "Specify the implementation for multithreading" )
        ( prefixvm( prefix,"parallel.cpu.restrict" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "Restrict the multithreading to N additional cores per MPI process (0: guess the maximum number of usable cores)" )
#if defined(HARTS_HAS_OPENCL)
        ( prefixvm( prefix,"parallel.opencl.enable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Enable the use of OpenCL for parallelization" )
        ( prefixvm( prefix,"parallel.opencl.device" ).c_str(), Feel::po::value<std::string>()->default_value("cpu"), "Specify the device to use for OpenCL (Valid entries: cpu or gpu" )
#endif
        ( prefixvm( prefix,"parallel.debug" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "Enable debugging for parallelization" )
#endif
        ;
    return _options;
}

po::options_description cnab2_options( std::string const& prefix )
{
    po::options_description _options( "CNAB2 (Crank-Nicolson Adam-Bashforth order 2 adaptive scheme) options (" + prefix + ")" );
    _options.add_options()
        ( prefixvm( prefix, "cnab2.time-initial" ).c_str(), po::value<double>()->default_value( 1e-8 ), "initial time step" )
        ( prefixvm( prefix, "cnab2.time-final" ).c_str(), po::value<double>()->default_value( 1 ), "target final time" )
        ( prefixvm( prefix, "cnab2.keps" ).c_str(), po::value<double>()->default_value( 1e-3 ), "time tolerance" )
        ( prefixvm( prefix, "cnab2.nstar" ).c_str(), po::value<int>()->default_value( 10 ), "averaging period" )
        ( prefixvm( prefix, "cnab2.steady" ).c_str(), po::value<bool>()->default_value( 0 ), "look for steady solution" )
        ( prefixvm( prefix, "cnab2.steady-tol" ).c_str(), po::value<double>()->default_value( 1e-6 ), "steady solution tolerance" )
        ( prefixvm( prefix, "cnab2.restart" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart " )
        ( prefixvm( prefix, "cnab2.restart.path" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "path where we reload old data" )
        ( prefixvm( prefix, "cnab2.restart.at-last-save" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart with ti the last save " )
        ( prefixvm( prefix, "cnab2.restart.step-before-last-save" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "do a restart with ti the ieme step before last save " )
        ( prefixvm( prefix, "cnab2.save" ).c_str(), Feel::po::value<bool>()->default_value( true ), "save elements in file " )
        ( prefixvm( prefix, "cnab2.save.freq" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "freq for save elements in file " );

    return _options;
}
/**
 * \return the command lines options for BDF
 */
po::options_description bdf_options( std::string const& prefix )
{
    po::options_description _options( "BDF (Backward Differences time discretization) options (" + prefix + ")" );
    _options.add_options()
    // solver options
    ( prefixvm( prefix, "bdf.time-initial" ).c_str(), Feel::po::value<double>()->default_value( 0.0 ), "initial time" )
    ( prefixvm( prefix, "bdf.time-final" ).c_str(), Feel::po::value<double>()->default_value( 1.0 ), "final time" )
    ( prefixvm( prefix, "bdf.time-step" ).c_str(), Feel::po::value<double>()->default_value( 1.0 ), "time step" )
    ( prefixvm( prefix, "bdf.strategy" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "strategy, 0=constant time steps, 1=adaptive time steps" )
    ( prefixvm( prefix, "bdf.steady" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "false: unsteady, true:steady" )
    ( prefixvm( prefix, "bdf.reverse" ).c_str(), Feel::po::value<bool>()->default_value( false ), "reverse time" )
    ( prefixvm( prefix, "bdf.restart" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart " )
    ( prefixvm( prefix, "bdf.restart.path" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "path where we reload old data" )
    ( prefixvm( prefix, "bdf.restart.at-last-save" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart with ti the last save " )
    ( prefixvm( prefix, "bdf.restart.step-before-last-save" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "do a restart with ti the ieme step before last save " )
    ( prefixvm( prefix, "bdf.iterations-between-order-change" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "iteration between time order change" )
    ( prefixvm( prefix, "bdf.save" ).c_str(), Feel::po::value<bool>()->default_value( true ), "save elements in file " )
    ( prefixvm( prefix, "bdf.save.freq" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "freq for save elements in file " )
    ( prefixvm( prefix, "bdf.rank-proc-in-files-name" ).c_str(), Feel::po::value<bool>()->default_value( false ), "the name of files generated has the rank of the processor automatically if true" )
    ( prefixvm( prefix, "bdf.file-format" ).c_str(), Feel::po::value<std::string>()->default_value( "binary" ), "save elements in file " )

    ;
    _options.add_options()
    ( prefixvm( prefix, "bdf.order" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "order in time" )
    ( prefixvm( prefix, "bdf.strategy-high-order-start" ).c_str(), Feel::po::value<int>()->default_value( 0 ), " 0 : fixe order, 1 : increase step by step order" )
    ;
    return _options;
}


po::options_description ts_options( std::string const& prefix )
{
    po::options_description _options( "BDF (Backward Differences time discretization) options (" + prefix + ")" );
    _options.add_options()
        // solver options
        ( prefixvm( prefix, "ts.time-initial" ).c_str(), Feel::po::value<double>()->default_value( 0.0 ), "initial time" )
        ( prefixvm( prefix, "ts.time-final" ).c_str(), Feel::po::value<double>()->default_value( 1.0 ), "final time" )
        ( prefixvm( prefix, "ts.time-step" ).c_str(), Feel::po::value<double>()->default_value( 1.0 ), "time step" )
        //( prefixvm( prefix, "ts.strategy" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "strategy, 0=constant time steps, 1=adaptive time steps" )
        ( prefixvm( prefix, "ts.steady" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "false: unsteady, true:steady" )
        ( prefixvm( prefix, "ts.reverse" ).c_str(), Feel::po::value<bool>()->default_value( false ), "reverse time" )
        ( prefixvm( prefix, "ts.restart" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart " )
        ( prefixvm( prefix, "ts.restart.path" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "path where we reload old data" )
        ( prefixvm( prefix, "ts.restart.at-last-save" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart with ti the last save " )
        ( prefixvm( prefix, "ts.restart.step-before-last-save" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "do a restart with ti the ieme step before last save " )
        ( prefixvm( prefix, "ts.save" ).c_str(), Feel::po::value<bool>()->default_value( true ), "save elements in file " )
        ( prefixvm( prefix, "ts.save.freq" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "freq for save elements in file " )
        ( prefixvm( prefix, "ts.rank-proc-in-files-name" ).c_str(), Feel::po::value<bool>()->default_value( false ), "the name of files generated has the rank of the processor automatically if true" )
        ( prefixvm( prefix, "ts.file-format" ).c_str(), Feel::po::value<std::string>()->default_value( "binary" ), "save elements in file " )
        ( prefixvm( prefix, "ts.display-stats").c_str(), po::value<bool>()->default_value(false), "display statistics (timers, iterations counts...) at each ts iteration" )
        ;
    return _options;
}

po::options_description stabilization_options( std::string const& prefix )
{
    po::options_description _options("Options for stabilization");
    _options.add_options()
        ( prefixvm( prefix, "stab.use" ).c_str(), po::value<bool>()->default_value( false ), "use or not stabilization method" )
        ( prefixvm( prefix, "stab.div" ).c_str(), po::value<bool>()->default_value( true ), "add stabilization on divergence" )
        ( prefixvm( prefix, "stab.compute-lambdaK" ).c_str(), po::value<bool>()->default_value( true ), "true : compute obtimized stabilization coefficient (highly recommended) / false use arbitrary constant" )
        ( prefixvm( prefix, "stab.type" ).c_str(), po::value<std::string>()->default_value( "douglas-wang"), "douglas-wang, gls or supg" )
        ( prefixvm( prefix, "stab.force" ).c_str(), po::value<bool>()->default_value( true ), "force the computation of stabilization terms at first iteration" )
        ( "stab.starting-Re", po::value<double>()->default_value( 0.9 ), "maximal value of local Reynolds number before the stabilization starts" )
        ( prefixvm( prefix, "stab.penal-lambdaK" ).c_str(), po::value<double>()->default_value( 1e-6 ), "value of the penalisation used in the computation of lambdaK" )

        ( prefixvm( prefix, "stab.export" ).c_str(), po::value<bool>()->default_value( true ), "export fields" )

        ( prefixvm( prefix, "stab.parameter" ).c_str(), po::value<int>()->default_value( 1 ), "different way to evaluate the stabilization parameter" )
        ;

    return _options.add( backend_options( prefixvm(prefix,"stab") ) );
}

po::options_description sc_options( std::string const& prefix )
{
    po::options_description _options("Options for static condensation");
    _options.add_options()
        ( prefixvm( prefix, "sc.condense" ).c_str(), po::value<bool>()->default_value( false ), "enable/disable static condensation" )
        ;

    return _options;
}


Feel::po::options_description
eimOptions( std::string const& prefix )
{
    Feel::po::options_description eimoptions( "EIM Options" );
    eimoptions.add_options()
        ( "eim.sampling-size"   , Feel::po::value<int>()->default_value( 30 ), "Offline  sampling size " )
        ( "eim.sampling-mode"   , Feel::po::value<std::string>()->default_value( "random" ), "EIM Offline : random, log-random, log-equidistribute, equidistribute " )
        ( "eim.error-max"   , Feel::po::value<double>()->default_value( 1e-6 ),       "Offline  tolerance" )
        ( "eim.online-tolerance"   , Feel::po::value<double>()->default_value( 1e-2 ),       "Online  tolerance" )
        ( "eim.dimension-max"   , Feel::po::value<int>()->default_value( 50 ),       "Offline  max WN size" )
        //( "eim.error-type"   , Feel::po::value<int>()->default_value( ( int )EIM_RESIDUAL ),       "EIM error type to be computed" )
        ( "eim.check.rb"   , Feel::po::value<int>()->default_value( 0 ),       "check reduced basis" )
        ( "eim.norm-used-for-residual" , Feel::po::value<std::string>()->default_value( "Linfty" ), "name of the norm used : Linfty - L2 - LinftyVec" )
        ( "eim.check.residual"   , Feel::po::value<int>()->default_value( 0 ),  "check residual" )
        ( "eim.reuse-prec"   , Feel::po::value<bool>()->default_value( 0 ),       "reuse or not the preconditioner" )
        ( "eim.rebuild-database" , Feel::po::value<bool>()->default_value( false ), "rebuild database (if it already exists)" )
        ( "eim.enrich-database" , Feel::po::value<bool>()->default_value( true ), "enrich the database (if it already exists)" )
        ( "eim.cvg-study" , Feel::po::value<bool>()->default_value( 0 ), "for convergence study" )
        ( "eim.compute-error-with-truth-expression" , Feel::po::value<bool>()->default_value( false ), "compute the error with the truth expression ( not its projection ) if true" )
        ( "eim.use-dimension-max-functions" , Feel::po::value<bool>()->default_value( 0 ), "force to use dimension-max basis functions" )
        ( "eim.computational-time-neval",Feel::po::value<int>()->default_value( 0 )," number of evaluation to perform to have the computational time of eim online step" )
        ( "eim.compute-expansion-of-expression",Feel::po::value<bool>()->default_value( false )," use true expression if true, else use the projection of the expression" )
        ( "eim.show-mu-selection",Feel::po::value<bool>()->default_value( false )," print list of parameters selected during offline step" )
        ( "eim.show-t-selection",Feel::po::value<bool>()->default_value( false )," print list of interpolation points selected during offline step" )
        ( "eim.show-offline-error",Feel::po::value<bool>()->default_value( false )," print list of error associated to mu selected during offline step" )

        ( "eim.elements.write", Feel::po::value<bool>()->default_value( true ), "Write evaluated nl solutions on disk"  )
        ( "eim.elements.directory", Feel::po::value<std::string>()->default_value( "nlsolutions" ), "directory were nl solutions are stored"  )
        ( "eim.elements.clean-directory", Feel::po::value<bool>()->default_value( false ), ""  )
        ;

    return eimoptions;
}

Feel::po::options_description
deimOptions( std::string const& prefix )
{
    Feel::po::options_description deimoptions( "DEIM Options" );
    deimoptions.add_options()
        ( prefixvm( prefix, "deim.dimension-max" ).c_str(), Feel::po::value<int>()->default_value( 20 ), "Offline  max WN size" )
        ( prefixvm( prefix, "deim.default-sampling-size" ).c_str(), Feel::po::value<int>()->default_value( 50 ), "Offline  sampling size"  )
        ( prefixvm( prefix, "deim.default-sampling-mode" ).c_str(), Feel::po::value<std::string>()->default_value( "equidistribute" ), "DEIM Offline : random, log-random, log-equidistribute, equidistribute "  )
        ( prefixvm( prefix, "deim.rebuild-database" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Rebuild the database from beginning if true"  )
        ( prefixvm( prefix, "deim.greedy.rtol" ).c_str(), Feel::po::value<double>()->default_value( 1e-8 ), "Asbolute Tolerance for greedy algorithm"  )
        ( prefixvm( prefix, "deim.greedy.atol" ).c_str(), Feel::po::value<double>()->default_value( 1e-16 ), "Relative Tolerance for greedy algorithm"  )
        ( prefixvm( prefix, "deim.store-vectors" ).c_str(), Feel::po::value<bool>()->default_value(true ), "Store Vectors for the parameters in the trainset in DEIM"  )
        ( prefixvm( prefix, "deim.store-matrices" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Store Matrices for the parameters in the trainset in MDEIM"  )

        ( prefixvm( prefix, "deim.elements.write" ).c_str(), Feel::po::value<bool>()->default_value( true ), "Write evaluated nl solutions on disk"  )
        ( prefixvm( prefix, "deim.elements.clean-directory" ).c_str(), Feel::po::value<bool>()->default_value( false ), ""  )
        ( prefixvm( prefix, "deim.elements.directory" ).c_str(), Feel::po::value<std::string>()->default_value( "nlsolutions" ), "directory were nl solutions are stored"  )

        ( prefixvm( prefix, "deim.optimized-online" ).c_str(), Feel::po::value<bool>()->default_value( true ), "Use optimized version for online assembly. DEBUG, this option has to be removed !"  )
        ;

    return deimoptions.add(backend_options(prefixvm(prefix,"deim-online")));
}

Feel::po::options_description
crbSEROptions( std::string const& prefix )
{
    Feel::po::options_description seroptions( "SER Options" );
    seroptions.add_options()
        ( "ser.rb-frequency", Feel::po::value<int>()->default_value( 0 ), "Number of RB basis built per step of SER process, 0 : no SER" )
        ( "ser.eim-frequency", Feel::po::value<int>()->default_value( 0 ), "Number of EIM basis built per step of SER process, 0 : no SER")
        ( "ser.use-rb-in-eim-mu-selection", Feel::po::value<bool>()->default_value( false ), "Use RB approx. to select parameters during EIM offline step")
        ( "ser.use-rb-in-eim-basis-build", Feel::po::value<bool>()->default_value( false ), "Use RB approx. to build EIM basis")
        ( "ser.rb-rebuild-freq", Feel::po::value<int>()->default_value( -1 ), "rebuild database with a given frequency" )
        ( "ser.error-estimation", Feel::po::value<bool>()->default_value( false ), "Use SER error estimation = Norm of residual (Riesz) as error indicator")
        ( "ser.use-greedy-in-rb", Feel::po::value<bool>()->default_value( false ), "Use SER error indicator to build RB basis (Greedy)")
        ( "ser.eim-subtrainset-method", Feel::po::value<int>()->default_value( 0 ), "Build subtrainset from error estimation : full trainset (=0), select mu for which residual < tol (=1), select mu for which residual > tol (=2)")
        ( "ser.print-rb-iterations_info", Feel::po::value<bool>()->default_value( false ), "Write online rb iterations (Picard) needed for Greedy")
        ( "ser.radapt-eim-rtol", Feel::po::value<double>()->default_value( 0.0 ), "Relative tolerance criterion for EIM r-adaptation - default(0.0) = no adaptation")
        ( "ser.radapt-rb-rtol", Feel::po::value<double>()->default_value( 0.0 ), "Relative tolerance criterion for RB r-adaptation - default(0.0) = no adaptation")
        ( "ser.eim-greedy-rtol", Feel::po::value<double>()->default_value( 0.0 ), "Relative tolerance criterion for the error indicator constraint in eim Greedy algorithm : choose mu from those which satisfies this criterion - default(0.0) = no error indicator constraint")
        ( "ser.corrected-rb-rtol", Feel::po::value<double>()->default_value( 0.0 ), "Relative tolerance criterion for RB correction (from error estimation) to be used - default(0.0) = no correction")
        ( "ser.nb-levels", Feel::po::value<int>()->default_value( 1 ), "number of SER levels = number of passages into SER algorithm (default 1)")
        ;
    return seroptions;
}

po::options_description
solvereigen_options( std::string const& prefix )
{
    std::string _prefix = prefix;
    boost::algorithm::trim( _prefix );

    if ( !_prefix.empty() && !boost::algorithm::ends_with( _prefix, "-" ) )
        _prefix += "-";


    //int nev,                  // number of requested eigenpairs
    //int ncv,                  // number of basis vectors
    //const double tol,         // solver tolerance
    //const unsigned int m_its) // maximum number of iterations
    po::options_description _options( "Solver EigenValue SLEPc -- " + prefix + " solver options" );
    _options.add_options()
    // solver options
        ( ( _prefix+"solvereigen.solver" ).c_str(), Feel::po::value<std::string>()->default_value( "krylovschur" ), "type of eigenvalue solver. Choice: power,lapack,subspace,arnoldi,krylovschur,arpack" )
        ( ( _prefix+"solvereigen.problem" ).c_str(), Feel::po::value<std::string>()->default_value( "ghep" ), "type of eigenvalue problem. Choice: nhep, hep, gnhep, ghep, pgnhep" )
        ( ( _prefix+"solvereigen.spectrum" ).c_str(), Feel::po::value<std::string>()->default_value( "largest_magnitude" ), "eigenvalue solver position in spectrum. Choice: largest_magnitude, smallest_magnitude, largest_real, smallest_real, largest_imaginary, smallest_imaginary" )
        ( ( _prefix+"solvereigen.transform" ).c_str(), Feel::po::value<std::string>()->default_value( "shift" ), "spectral transformation. Choice: shift, shift_invert, fold, cayley" )
        ( ( _prefix+"solvereigen.nev" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "number of requested eigenpairs" )
        ( ( _prefix+"solvereigen.ncv" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "dimension of the subspace" )
        ( ( _prefix+"solvereigen.mpd" ).c_str(), Feel::po::value<int>()->default_value( -2 ), "maximum projected dimension" )
        ( ( _prefix+"solvereigen.eps-target" ).c_str(), Feel::po::value<double>()->default_value( NAN ), "sets the value of the target" )
        ( ( _prefix+"solvereigen.interval-a" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "start of the interval in which all the eigenvalues are found" )
        ( ( _prefix+"solvereigen.interval-b" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "end of the interval in which all the eigenvalues are found" )
#if defined(PETSC_HAVE_MUMPS)
        ( ( _prefix+"solvereigen.st-pc-factor-mat-solver-package-type" ).c_str(), Feel::po::value<std::string>()->default_value( "mumps" ), "sets the software that is used to perform the factorization for the spectral transformation (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#else
        ( ( _prefix+"solvereigen.st-pc-factor-mat-solver-package-type" ).c_str(), Feel::po::value<std::string>()->default_value( "petsc" ), "sets the software that is used to perform the factorization for the spectral transformation (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#endif
        ( ( _prefix+"solvereigen.krylovschur-partitions" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "number of communicator to use for interval" )
        ( ( _prefix+"solvereigen.tolerance" ).c_str(), Feel::po::value<double>()->default_value( 1e-10 ), "solver tolerance" )
        ( ( _prefix+"solvereigen.maxiter" ).c_str(), Feel::po::value<int>()->default_value( 10000 ), "maximum number of iterations" )
        ( ( _prefix+"solvereigen.verbose" ).c_str(), Feel::po::value<bool>()->default_value( false ), "verbose eigen solver" )
        ( ( _prefix+"solvereigen.eps-monitor" ).c_str(), Feel::po::value<bool>()->default_value( false ), "monitor eigen problem solver" );

    return _options;
}

Feel::po::options_description
crbSCMOptions( std::string const& prefix = "")
{
    Feel::po::options_description crbscmoptions( "CRB SCM Options" );
    crbscmoptions.add_options()
    ( "crb.scm.sampling-size"   , Feel::po::value<int>()->default_value( 100 ),       "Offline SCM sampling size " )
    ( "crb.scm.tol"   , Feel::po::value<double>()->default_value( 1e-6 ),       "Offline SCM tolerance" )
    ( "crb.scm.iter-max"   , Feel::po::value<int>()->default_value( 10 ),       "Offline SCM max number of K iterations" )
    ( "crb.scm.Mplus" , Feel::po::value<int>()->default_value( 10 ),       "M+ value" )
    ( "crb.scm.Malpha" , Feel::po::value<int>()->default_value( 10 ),       "M_alpha value" )
    ( "crb.scm.level" , Feel::po::value<int>()->default_value( 1 ),       "level for recursion in lower bound computations" )
    ( "crb.scm.strategy" , Feel::po::value<int>()->default_value( 2 ),       "scm strategy (0=patera, 1=maday, 2=prudhomme" )
    ( "crb.scm.rebuild-database" , Feel::po::value<bool>()->default_value( 0 ), "rebuild database (if it already exists)" )
    ( "crb.scm.do-scm-for-mass-matrix" , Feel::po::value<bool>()->default_value( 0 ), "do scm for bilinear form m(.,.;mu) and not for a(.,.;mu) " )
    ( "crb.scm.print-matrix" , Feel::po::value<bool>()->default_value( 0 ), "print matrix " )
    ( "crb.scm.solvereigen-tol" ,  Feel::po::value<double>()->default_value( 1e-10 ), "solver eigen tolerance " )
    ( "crb.scm.solvereigen-maxiter" ,  Feel::po::value<int>()->default_value( 10000 ), "solver eigen maxiter " )
    ( "crb.scm.solvereigen-nev" ,  Feel::po::value<int>()->default_value( 1 ), "solver eigen nev " )
    ( "crb.scm.solvereigen-ncv" ,  Feel::po::value<int>()->default_value( 3 ), "solver eigen ncv " )
    ( "crb.scm.solvereigen-solver-type" ,  Feel::po::value<int>()->default_value( 5 ), "solver eigen type " )
    ( "crb.scm.cvg-study",Feel::po::value<bool>()->default_value( false ), "convergence study if true")
    ( "crb.scm.run-on-C",Feel::po::value<bool>()->default_value( false ), "use parameters selected in offline step if true ( in that case, Lb=Ub=FEM )")
    ( "crb.scm.use-logEquidistributed-C",Feel::po::value<int>()->default_value( 0 ), "parameters are log-equidistributed for the offline step (the value indicates the number of parameters)")
    ( "crb.scm.use-equidistributed-C",Feel::po::value<int>()->default_value( 0 ), "parameters are equidistributed for the offline step (the value indicates the number of parameters)")
    ( "crb.scm.use-predefined-C",Feel::po::value<bool>()->default_value( false ), "use a predefined sampling C ( parameters written on the file SamplingC")
    ( "crb.scm.use-scm",Feel::po::value<bool>()->default_value( false ), "use scm if true")
    ( "crb.scm.check-eigenvector",Feel::po::value<bool>()->default_value( true ), "check that eigenvector and eigenvalue are solution of the generalized eiganvalue problem if true")
    ( "crb.scm.check-eigenvector-tol",Feel::po::value<double>()->default_value( 1e-11 ), "tolerance of check-eigenvector")
    ;

    crbscmoptions
    .add( solvereigen_options( "crb.scm" ) );

    return crbscmoptions;
}

Feel::po::options_description
crbOptions( std::string const& prefix )
{
    Feel::po::options_description crboptions( "CRB Options" );
    crboptions.add_options()
        ( "crb.db.load"   , Feel::po::value<int>()->default_value( 2 ), "=0 use db.filename, =1 use last DB created =2 use last DB modified =3 use db.id" )
        ( "crb.db.update"   , Feel::po::value<int>()->default_value( 2 ), "=0 use db.filename, =1 use last DB created =2 use last DB modified =3 use db.id" )
        ( "crb.db.id"   , Feel::po::value<std::string>(), "DB unique id" )
        ( "crb.db.filename"   , Feel::po::value<std::string>(), "DB filename" )
        ( "crb.output-index"   , Feel::po::value<int>()->default_value( 0 ), "output index (default is right hand side = 0)" )
        ( "crb.sampling-size"   , Feel::po::value<int>()->default_value( 1000 ), "Offline  sampling size " )
        ( "crb.sampling-mode"   , Feel::po::value<std::string>()->default_value( "log-random" ), "Offline  sampling mode, equidistribute, log-equidistribute or log-random " )
        ( "crb.all-procs-have-same-sampling" , Feel::po::value<bool>()->default_value( false ), "all procs have the same sampling if true" )
        ( "crb.error-max"   , Feel::po::value<double>()->default_value( 1e-6 ),       "Offline  tolerance" )
        ( "crb.online-tolerance"   , Feel::po::value<double>()->default_value( 1e-2 ),       "Online  tolerance" )
        ( "crb.absolute-error" , Feel::po::value<bool>()->default_value( false ), "Impose to compute absolute error PFEM/CRB instead of relative" )
        ( "crb.dimension-max"   , Feel::po::value<int>()->default_value( 1 ),       "Offline max WN size, set to 1 by default to avoid to enrich the existing database if this option doesn't appear in onefeel interface or in the config file." )
        ( "crb.dimension"   , Feel::po::value<int>()->default_value( -1 ),       "Online  WN size" )
        ( "crb.error-type"   , Feel::po::value<int>()->default_value( ( int )CRB_RESIDUAL_SCM ),       "CRB error type to be computed" )
        ( "crb.compute-apee-for-each-time-step",Feel::po::value<bool>()->default_value( true ),"Compute error estimation for each time step (parabolic problems) is true, else compute only for the last one")
        ( "crb.factor"   , Feel::po::value<int>()->default_value( -1 ),  "factor useful to estimate error by empirical method" )
        ( "crb.Nm"   , Feel::po::value<int>()->default_value( 1 ),       "Offline  number of modes per mu (for the POD) " )
        ( "crb.apply-POD-to-WN"   , Feel::po::value<bool>()->default_value( false ), "apply a POD on approximation functions spaces (primal and dual) if true and if deal with a transient problem " )
        ( "crb.check.rb"   , Feel::po::value<int>()->default_value( 0 ),       "check reduced basis" )
        ( "crb.orthonormality-tol" , Feel::po::value<double>()->default_value( 1e-13 ),"tolerance of orthonormalisation : i.e. norm of matrix A(i,j)=scalarProduct( Wn[j], Wn[i] )" )
        ( "crb.orthonormality-max-iter" , Feel::po::value<int>()->default_value( 10 ),"while the tolerance is not reached, the orthonormalization step is done or until max-iter is reached" )

        ( "crb.check.residual"   , Feel::po::value<bool>()->default_value( false ),  "check residual" )
        ( "crb.reuse-prec"   , Feel::po::value<bool>()->default_value( 0 ),       "reuse or not the preconditioner" )
        ( "crb.use-primal-pc",Feel::po::value<bool>()->default_value( true ), "use specific preconditioner for the primal problem")
        ( "crb.orthonormalize-primal" , Feel::po::value<bool>()->default_value( 1 ), "orthonormalize or not " )
        ( "crb.orthonormalize-dual" , Feel::po::value<bool>()->default_value( 1 ), "orthonormalize or not " )
        ( "crb.solve-dual-problem" , Feel::po::value<bool>()->default_value( 1 ), "solve or not the dual problem (this bool will be ignored if error-type=CRB_RESIDUAL) " )
        ( "crb.visualize-basis" , Feel::po::value<bool>()->default_value( 0 ), "visualize elements of the reduced basis " )
        ( "crb.save-output-behavior" , Feel::po::value<bool>()->default_value( 0 ), "save output behavior in time" )
        ( "crb.seek-mu-in-complement" , Feel::po::value<bool>()->default_value( 1 ), "during the offline basis construction, see mu in M the complement of Wn" )
        ( "crb.rebuild-database" , Feel::po::value<bool>()->default_value( 0 ), "rebuild database (if it already exists)" )
        ( "crb.restart-from-N" , Feel::po::value<int>()->default_value( -1 ), "restart the database from specified N (note that when N=0 the complete approximation space is rebuilt and so it is equivalent to option crb.rebuild-database=true). In the case where N > Nmax we do nothing. By default it is set to a negative number in order to not interfer with option rebuild-database" )
        ( "crb.show-mu-selection" , Feel::po::value<bool>()->default_value( 0 ), " show mu selection during offline step to build RB space" )
        ( "crb.show-residual" , Feel::po::value<bool>()->default_value( 0 ), " show mu residuals values (used for the error estimation)" )
        ( "crb.print-error-during-rb-construction" , Feel::po::value<bool>()->default_value( 0 ), " print the max error (absolute) obtained during the offline step" )
        ( "crb.print-iterations-info", Feel::po::value<bool>()->default_value( 0 ), "Write offline Picard iterations summary" )
        ( "crb.compute-variance" , Feel::po::value<bool>()->default_value( 0 ), " if true the output is the variance and not l(v)" )
        ( "crb.save-information-for-variance",Feel::po::value<bool>()->default_value( 0 ), "if true will build variance matrix but it takes some times" )

        ( "crb.use-newton",Feel::po::value<bool>()->default_value( false ), "use newton algorithm (need to provide a jacobian and a residual)" )
        ( "crb.fixedpoint.aitken",Feel::po::value<bool>()->default_value( true ), "use Aitken relaxtion algorithm in nonlinear fixpoint solver" )

        ( "crb.fixedpoint.maxit",Feel::po::value<int>()->default_value( 20 ), "nb iteration max for the fixed point (online part)" )
        ( "crb.fixedpoint.increment-tol",Feel::po::value<double>()->default_value( 1e-10 ), "tolerance on solution for fixed point (online part)" )
        ( "crb.fixedpoint.output-tol",Feel::po::value<double>()->default_value( 1e-10 ), "tolerance on output for fixed point (online part)" )
        ( "crb.fixedpoint.verbose",Feel::po::value<bool>()->default_value( false ), "fixed point verbose if true" )
        ( "crb.fixedpoint.critical-value",Feel::po::value<double>()->default_value(1000 ), "will crash if increment error at the end of fixed point is greater than this critical value" )

        ( "crb.use-continuity",Feel::po::value<bool>()->default_value(true), "when apply Newton method, will use continuity method (like for natural convection problem for example)" )
        ( "crb.compute-error-on-reduced-residual-jacobian",Feel::po::value<bool>()->default_value( false ), "only for crb_trilinear")
        ( "crb.enable-convection-terms",Feel::po::value<bool>()->default_value( true ), "only for crb_trilinear")

        ( "crb.is-model-executed-in-steady-mode",Feel::po::value<bool>()->default_value( false ), "true if model is executed in steady mode, else turn it to false")
        ( "crb.use-ginac-for-beta-expressions",Feel::po::value<bool>()->default_value( false ), "use ginac to compute expression of beta coefficients if true")
        ( "crb.use-linear-model",Feel::po::value<bool>()->default_value( false ), "do not iterate in fixed point if true")

        ( "crb.use-predefined-WNmu",Feel::po::value<bool>()->default_value( false ), "read parameters to take for the offline step from a file named SamplingWNmu if true")
        ( "crb.use-predefined-test-sampling",Feel::po::value<bool>()->default_value( false ), "read parameters from file named SamplingForTest if true to run the test")
        ( "crb.use-logEquidistributed-WNmu",Feel::po::value<int>()->default_value( 0 ), "parameters are log-equidistributed for the offline step (the value indicates the number of parameters)")
        ( "crb.use-random-WNmu",Feel::po::value<int>()->default_value( 0 ), "parameters are log-equidistributed for the offline step (the value indicates the number of parameters)")
        ( "crb.use-equidistributed-WNmu",Feel::po::value<int>()->default_value( 0 ), "parameters are equidistributed for the offline step (the value indicates the number of parameters)")

        ( "crb.compute-stat",Feel::po::value<bool>()->default_value( true ), "compute statistics on the run if true")
        ( "crb.cvg-study",Feel::po::value<bool>()->default_value( false ), "convergence study if true")
        ( "crb.computational-time-neval",Feel::po::value<int>()->default_value( 0 )," number of evaluation to perform to have the computational time of crb online step" )

        ( "crb.run-on-WNmu",Feel::po::value<bool>()->default_value( false ), "use mu taken for build the reduced basis, so for steady problems we are very accurate")
        ( "crb.run-on-scm-parameters",Feel::po::value<bool>()->default_value( false ), "use mu taken during the SCM offline step ( for a(.,.;mu) ), so the coercivity constant is exact")
        ( "crb.script-mode",Feel::po::value<bool>()->default_value( false ), "disable error computation (need FEM computation) if true")
        ( "crb.db.format", Feel::po::value<std::string>()->default_value("boost"), "format in which the crb database is saved, either boost of hdf5")
        ( "crb.results-repo-name", Feel::po::value<std::string>()->default_value("default_repo"), "name for results repository, and also use for database storage")
        ( "crb.compute-fem-during-online",Feel::po::value<bool>()->default_value( true ), "compute fem during online step, necessary to compute the error between fem and crb")

        ( "crb.compute-matrix-information",Feel::po::value<bool>()->default_value( false ), "compute matrix information (i.e. conditioning, determinant) of reduced matrix if true")
        ( "crb.print-rb-matrix",Feel::po::value<bool>()->default_value( false ), "write rb matrix (octave format) in a file if true")

        ( "crb.use-symmetric-matrix",Feel::po::value<bool>()->default_value( true ), "don't transpose to have the matrix associated to the dual problem if true")
        ( "crb.stock-matrices",Feel::po::value<bool>()->default_value( true ), "assemble and stock all matrices/vectors if true, but it can takes a lot of memory")
        ( "crb.system-memory-evolution",Feel::po::value<bool>()->default_value( false ), "generate a file to plot memory evolution during offline step only on the master processor (file written : MemoryEvolution)")
        ( "crb.system-memory-evolution-on-all-procs",Feel::po::value<bool>()->default_value( false ), "same than system-memory-evolution but on all processors")

        ( "crb.use-accurate-apee",Feel::po::value<bool>()->default_value( false ), "use a posteriori error estimators from F.Casenave's paper if true, classic one else")
        ( "crb.optimize-offline-residual",Feel::po::value<bool>()->default_value( false ), "use optimize way for offline residual computation if true (temporary option)")
        ( "crb.offline-residual-version",Feel::po::value<int>()->default_value( 0 ), "offline residual version 0 : old no storage, 1 : new faster with storage")

        ( "crb.user-parameters",Feel::po::value<std::string>()->default_value( "" ), "values of parameters (used for one feel)")
        ( "crb.vary-only-parameter-components",Feel::po::value<std::string>()->default_value( "" ), "specify which parameter component vary (max : 2 components + time 't') and how many values we take in each direction. For example 0 10 1 20 means that component 0 will take 10 values and component 1 will take 20 values. We can write also t 1 10 and in this case the time will vary and also component 1")
        ( "crb.load-elements-database",Feel::po::value<bool>()->default_value( false ), "load database of elements if true, need to be true for visualization, need to be false to run CRB approximation on a different number of processors than this was used to build the reduced basis ")

        ( "crb.solve-fem-monolithic",Feel::po::value<bool>()->default_value( false ), "solve FEM problem without using EIM and without affine decomposition ")
        ( "crb.export-name-max-size",Feel::po::value<int>()->default_value( 30 ), "maximum size for variable names in export (truncature)")

        ("crb.minimization-func", Feel::po::value<std::string>(), "giving a functional f(output) - give the output which minimizes f(output)" )
        ("crb.minimization-param-name", Feel::po::value<std::string>()->default_value( "output" ), "name of the parameter to be replaced by the output in expression given by crb.minimization-func")

        ( "crb.use-fast-eim",Feel::po::value<bool>()->default_value( true ), "use fast eim algo (with rbspace context)")

        ;

    crboptions
        .add( crbSCMOptions() ).add(deimOptions(prefix));

    return crboptions;
}

    Feel::po::options_description
    crbSaddlePointOptions( std::string const& prefix )
{
    Feel::po::options_description crboptions( "CRB Options" );
    crboptions.add_options()
        ( "crb.saddlepoint.add-supremizer",Feel::po::value<bool>()->default_value( false ), "add the supremizer function to the first reduced basis")
        ( "crb.saddlepoint.orthonormalize0",Feel::po::value<bool>()->default_value( true ), "orthonormalize reduce basis for rbspace #0")
        ( "crb.saddlepoint.orthonormalize1",Feel::po::value<bool>()->default_value( true ), "orthonormalize reduce basis for rbspace #1")
        ( "crb.saddlepoint.version",Feel::po::value<int>()->default_value( 1 ), "test residual evaluation")
        ;

    crboptions.add( backend_options("backend-Xh0") );
    crboptions.add( backend_options("backend-Xh1") );

    return crboptions;
}


Feel::po::options_description
podOptions( std::string const& prefix )
{
    Feel::po::options_description podoptions( "POD Options" );
    podoptions.add_options()
    ( "pod.store-pod-matrix"   , Feel::po::value<bool>()->default_value( false ), "indicate if we store the pod matrix on a file" )
    ( "pod.store-pod-matrix-format-octave"   , Feel::po::value<bool>()->default_value( false ), "indicate if we store the pod matrix on a file with octave format" )
    ("pod.check-orthogonality",Feel::po::value<bool>()->default_value( true ), "check orthogonality of modes")
    ("pod.minimum-eigenvalue",Feel::po::value<double>()->default_value( 1e-11 ), "minimum acceptable value for eigenvalues")
    ("pod.check-tol",Feel::po::value<double>()->default_value( 1e-10 ), "when solving A w = lambda w, check that norm(A w) = norm( lambda w)")
    ;

    return podoptions;
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
    po::options_description _options( "Error " + prefix + " options" );
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

// preconditioner for Navier-Stokes
po::options_description
blockns_options( std::string const& prefix )
{
    po::options_description _options( "BLOCKNS options (" + prefix + ")" );
    _options.add_options()
        // error options
        ( prefixvm( prefix, "blockns" ).c_str(), Feel::po::value<bool>()->default_value(false), "enable BLOCKNS preconditioner" )
        ( prefixvm( prefix, "blockns.cd" ).c_str(), Feel::po::value<bool>()->default_value(false), "enable BLOCKNS/Velocity CD preconditioner" )
        ( prefixvm( prefix, "blockns.pcd" ).c_str(), Feel::po::value<bool>()->default_value(false), "enable BLOCKNS/Pressure CD preconditioner" )
        ( prefixvm( prefix, "blockns.pcd.inflow" ).c_str(), Feel::po::value<std::string>()->default_value("Robin"), "Type of boundary conditions at inflow: Robin or Dirichlet" )
        ( prefixvm( prefix, "blockns.pcd.outflow" ).c_str(), Feel::po::value<std::string>()->default_value("Dirichlet"), "Type of boundary conditions at inflow: Neumann or Dirichlet" )
        ( prefixvm( prefix, "blockns.pcd.order" ).c_str(), Feel::po::value<int>()->default_value(1), "order for pcd operator 1:Ap^-1 Fp Mp^-1 other: Mp^-1 Fp Ap^-1" )
        ( prefixvm( prefix, "blockns.pcd.diffusion" ).c_str(), Feel::po::value<std::string>()->default_value("Laplacian"), "Laplacian or BTBt" )
        ( prefixvm( prefix, "blockns.weakdir" ).c_str(), Feel::po::value<bool>()->default_value(0), "set to true for Weak dirichlet conditions for Fp and Ap, false otherwise" )
        ( prefixvm( prefix, "blockns.weakdir.penaldir" ).c_str(), Feel::po::value<double>()->default_value(10.), "Penalisation parameter for weak bc" )
        // options for pmm
        ( prefixvm( prefix, "blockns.pmm.diag" ).c_str(), Feel::po::value<bool>()->default_value(1), "set to true to use diagonal of the pressure mass matrix, false otherwise" )
        ;
    return _options;
}

// preconditioner for Magneto Static - regularized
po::options_description
ams_options( std::string const& prefix )
{
    po::options_description _options( "AMS options (" + prefix + ")" );
    _options.add_options()
        ( prefixvm( prefix, "useEdge" ).c_str(), Feel::po::value<bool>()->default_value(true), "true: SetConstantEdgeVector, false: SetCoordinates" )
        //( prefixvm( prefix, "threshold" ).c_str(), Feel::po::value<bool>()->default_value(false), "remove near null values in Grad" )
        ( prefixvm( prefix, "setAlphaBeta" ).c_str(), Feel::po::value<bool>()->default_value(false), "Use locally constructed A_alpha and A_beta" )
        //( prefixvm( prefix, "singular" ).c_str(), Feel::po::value<bool>()->default_value(false), "Force the relaxation parameter to be zero" )
        ( prefixvm( prefix, "relax" ).c_str(), Feel::po::value<double>()->default_value(1.), "Relaxation parameter" )
        ;
    return _options
        .add( backend_options( prefix.c_str() ))  // for AMS
        /* if ams.pc-type == AS */
        .add( backend_options( prefixvm(prefix, "1").c_str() )) // the (1,1).1 block
        .add( backend_options( prefixvm(prefix, "2").c_str() )) // the (1,1).2 block
        ;
}

// preconditioner for Magneto Static - saddle point
po::options_description
blockms_options( std::string const& prefix )
{
    po::options_description _options( "BLOCKMS options (" + prefix + ")" );
    _options.add_options()
        ( prefixvm( prefix, "penaldir" ).c_str(), Feel::po::value<double>()->default_value(10.), "Penalisation parameter for weak bc" );
    return _options
        .add(     ams_options( prefixvm(prefix, "11").c_str() ))  // the (1,1) block
        .add( backend_options( prefixvm(prefix, "22").c_str() )); // the (2,2) block
}

/**
 * \return the command lines options for the exporter
 */
po::options_description
exporter_options( std::string const& prefix )
{
    po::options_description _options( "Exporter " + prefix + " options" );
    _options.add_options()
        // do export
        ( prefixvm( prefix,"export" ).c_str(), Feel::po::value<bool>()->default_value( true ), "true if export, false otherwise" )
        // do export
        ( prefixvm( prefix,"exporter.export" ).c_str(), Feel::po::value<bool>()->default_value( true ), "true if export, false otherwise" )

        // exporter type
        ( prefixvm( prefix,"exporter.format" ).c_str(), Feel::po::value<std::string>()->default_value( "ensightgold" ), "type of exporter. Choices: ensight, ensightgold, gmsh" )



        //  geometry
        ( prefixvm( prefix,"exporter.geometry" ).c_str(), Feel::po::value<std::string>()->default_value( "change_coords_only" ), "Mesh change type, this option tells the exporter whether the mesh does not change(static), changes only the coordinates of the vertices (change_coords_only) or changes entirely (change). Choices: change_coords_only, change, static" )

        // prefix options
        ( prefixvm( prefix,"exporter.prefix" ).c_str(), Feel::po::value<std::string>()->default_value( prefix ), "prefix for exported files" )

        // directory options
        ( prefixvm( prefix,"exporter.directory" ).c_str(), Feel::po::value<std::string>()->default_value( "results" ), "directory for exported files" )

        // frequency options
        ( prefixvm( prefix,"exporter.freq" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "frequency at which results are exported" )

        // file type options
        ( prefixvm( prefix,"exporter.file-type" ).c_str(), Feel::po::value<std::string>()->default_value( "ascii" ), "file type in which the results are exported. Choices: ascii, binary" )

        // matlab options
        ( prefixvm( prefix,"exporter.matlab" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "export matrices and vectors to matlab files" )

        // spaces options (P0 to P1 smoothing)
        ( prefixvm( prefix,"exporter.element-spaces" ).c_str(), Feel::po::value<std::string>()->default_value( "P0" ), "spaces for P0 fields export. Choices : P0 (=only P0 visu), P1 (=smoothed P0, only P1 visu), P0+P1 (P0 + smoothed P0 visu)" )

        //
        // ensightgold
        //
        ( prefixvm( prefix,"exporter.ensightgold.use-sos" ).c_str(), Feel::po::value<bool>()->default_value( true ), "use sos  (true) or first case file (false) for multiple case files" )
        ( prefixvm( prefix,"exporter.ensightgold.save-face" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Save nodal functions on the face" )
        //  single
        ( prefixvm( prefix,"exporter.fileset" ).c_str(), Feel::po::value<bool>()->default_value( false ), "use fileset for transient simulations" )

        // merge timeteps or domains into single files
        ( prefixvm( prefix,"exporter.ensightgold.merge.markers" ).c_str(), Feel::po::value<bool>()->default_value( true ), "Merge exported data from different markers into a single file (reduces the number of output files )" )
        ( prefixvm( prefix,"exporter.ensightgold.merge.timesteps" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Merge exported timesteps into a single file (reduces the number of output files)" )
        ( prefixvm( prefix,"exporter.ensightgold.pack.timesteps" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "Allows to set the number of timesteps that will be stored in a single file, before switching to a new one. This option is meant to be used with --exporter.ensightgold.merge.timesteps. A value <= 0 means that all timesteps will go in the same file" )
        ( prefixvm( prefix,"exporter.gmsh.merge" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Merge exported data from different into a single file (reduces the number of output files)" )
#if defined(FEELPP_HAS_VTK) && defined(FEELPP_VTK_INSITU_ENABLED)
        ( prefixvm( prefix,"exporter.vtk.insitu.enable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Enable In-situ visualization with VTK exporter (Data won't be written to disk any longer, see --exporter.vtk.insitu.save to enable it)." )
        ( prefixvm( prefix,"exporter.vtk.insitu.pyscript" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "Specify a python user script for visualization." )
        ( prefixvm( prefix,"exporter.vtk.insitu.save" ).c_str(), Feel::po::value<bool>()->default_value( false ), "When this option is enabled, data will be written to disk when in using In-situ visualization." )
        ( prefixvm( prefix,"exporter.vtk.insitu.hostname" ).c_str(), Feel::po::value<std::string>()->default_value( "localhost" ), "Specify a hostname to which the simulation will connect for coprocessing." )
        ( prefixvm( prefix,"exporter.vtk.insitu.port" ).c_str(), Feel::po::value<int>()->default_value( 22222 ), "Specify the connection port used for coprocessing." )
#endif

        ;
    return _options;
}

po::options_description aitken_options( std::string const& prefix )
{
    po::options_description _options( "Aitken " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"aitken.type" ).c_str(), Feel::po::value<std::string>()->default_value( "method1" ), "standard,method1,fixed-relaxation" )
        ( prefixvm( prefix,"aitken.maxit" ).c_str(), Feel::po::value<int>()->default_value( 1000 ), "maximum number of iteration" )
        ( prefixvm( prefix,"aitken.initial_theta" ).c_str(), Feel::po::value<double>()->default_value( 1.0  ), "initial theta" )
        ( prefixvm( prefix,"aitken.min_theta" ).c_str(), Feel::po::value<double>()->default_value( 1e-4 ), "if theta computed < min_theta else theta=initial_theta" )
        ( prefixvm( prefix,"aitken.tol" ).c_str(), Feel::po::value<double>()->default_value( 1e-6 ), "fix-point tolerance" )
        ;
    return _options;
}

po::options_description
msi_options( std::string const& prefix )
{
    po::options_description _options( "Multiscale image " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"msi.level" ).c_str(), Feel::po::value<int>()->default_value( 0  ), "Coarsening level to pass from a fine grid to a coarse one" )
        ( prefixvm( prefix,"msi.pixelsize" ).c_str(), Feel::po::value<double>()->default_value( 8.9e-3  ), " Pixel size " )

               ;
    return _options;
}

// fit options
po::options_description
fit_options( std::string const& prefix )
{
    po::options_description _options( "Fit " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"fit.datafile" ).c_str(), Feel::po::value<std::string>(), "X - f(X) data file measures" )
        ( prefixvm( prefix,"fit.kind" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "Kind of interpolator : P0 (=0), P1 (=1), Spline (=2), Akima (=3)" )
        ( prefixvm( prefix,"fit.type" ).c_str(), Feel::po::value<std::string>()->default_value( "Akima" ), "type of interpolator : P0 , P1, Spline, Akima" )
        ( prefixvm( prefix,"fit.P0" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "left = 0, right = 1, center = 2" )
        ( prefixvm( prefix,"fit.P1_right" ).c_str(), Feel::po::value<int>()->default_value( 0  ), "zero = 0, constant = 1, extrapol = 2" )
        ( prefixvm( prefix,"fit.P1_left" ).c_str(), Feel::po::value<int>()->default_value( 1  ), "zero = 0, constant = 1, extrapol = 2" )
        ( prefixvm( prefix,"fit.Spline_right" ).c_str(), Feel::po::value<int>()->default_value( 0  ), "natural = 0, clamped = 1" )
        ( prefixvm( prefix,"fit.Spline_left" ).c_str(), Feel::po::value<int>()->default_value( 0  ), "natural = 0, clamped = 1" )

               ;
    return _options;
}

po::options_description
checker_options( std::string const& prefix )
{
    po::options_description _options( "Checker " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"checker.check" ).c_str(), Feel::po::value<bool>()->default_value(false), "run the check" )
        ( prefixvm( prefix,"checker.solution" ).c_str(), Feel::po::value<std::string>()->default_value("0"), "solution against which to check" )
        ( prefixvm( prefix,"checker.tolerance.exact" ).c_str(), Feel::po::value<double>()->default_value(1e-15), "tolerance for numerical exact solution check" )
        ( prefixvm( prefix,"checker.tolerance.order" ).c_str(), Feel::po::value<double>()->default_value(1e-1), "tolerance for order check" )
        ( prefixvm( prefix,"checker.name" ).c_str(), Feel::po::value<std::string>()->default_value("checker"), "name of the test" )
        ( prefixvm( prefix,"checker.filename" ).c_str(), Feel::po::value<std::string>()->default_value("checker.json"), "name of the test" )
        ;
    return _options;
}

po::options_description
journal_options( std::string const& prefix )
{
    po::options_description _options( "Journal " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"journal.enable" ).c_str(), Feel::po::value<bool>()->default_value(true), "enable simulation info" )
        ( prefixvm( prefix,"journal.filename" ).c_str(), Feel::po::value<std::string>()->default_value("simulation.info.json"), "name of the json file" )
        ( prefixvm( prefix,"journal.enable.mesh" ).c_str(), Feel::po::value<bool>()->default_value(true), "enable simulation info" )
        ( prefixvm( prefix,"journal.enable.functionspace" ).c_str(), Feel::po::value<bool>()->default_value(true), "enable simulation info" )
        ( prefixvm( prefix,"journal.enable.timer" ).c_str(), Feel::po::value<bool>()->default_value(true), "enable timer stats" )
        ( prefixvm( prefix,"journal.enable.hwsys" ).c_str(), Feel::po::value<bool>()->default_value(true), "enable hardware info" )
        ;
    return _options;
}

po::options_description
fmu_options( std::string const& prefix )
{
    po::options_description _options( "FMU " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"fmu.verbose" ).c_str(), Feel::po::value<bool>()->default_value(true), "Run with verbosity" )
        ( prefixvm( prefix,"fmu.display-variables-info" ).c_str(), Feel::po::value<bool>()->default_value(false), "Display the list of variables and description after loading FMU" )
        ( prefixvm( prefix,"fmu.filename" ).c_str(), Feel::po::value<std::string>()->default_value(""), "The fmu to load" )

        ( prefixvm( prefix,"fmu.solver.time-step" ).c_str(), Feel::po::value<double>()->default_value(0.1), "Time step for FMU Solver" )
        ( prefixvm( prefix,"fmu.solver.rtol" ).c_str(), Feel::po::value<double>()->default_value(1e-4), "Relative tolerance for FMU Solver" )

        ( prefixvm( prefix,"fmu.exported-variables" ).c_str(), po::value<std::vector<std::string>>()->multitoken(), "List of variables which have to be exported" )
        ( prefixvm( prefix,"fmu.export-directory" ).c_str(), po::value<std::string>()->default_value(""), "Location of the exported data. Default is the app directory." )

        ( prefixvm( prefix,"fmu.time-initial" ).c_str(), Feel::po::value<double>()->default_value(-1), "inital time for the simulation. Default is taken from the model" )
        ( prefixvm( prefix,"fmu.time-final" ).c_str(), Feel::po::value<double>()->default_value(-1), "Final time for the simulation. Default is taken from the model" )
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
        .add( backend_options("Ap") )
        .add( backend_options("Fp") )
        .add( backend_options("Mp") )
        .add( backend_options("Fu") )
        .add( backend_options("Bt") )
        .add( blockns_options( prefix ) )
        .add( sc_options( prefix ) )
        //.add( blockms_options( prefix ) )

        /* nonlinear solver options */
        .add( nlsolver_options() )

        /* discr options */
        .add( ts_options( prefix ) )
        .add( bdf_options( prefix ) )
        .add( cnab2_options( prefix ) )

        /* exporter options */
        .add( exporter_options( prefix ) )

        /* nlopt options */
#if defined(FEELPP_HAS_NLOPT)
        .add( nlopt_options( prefix ) )
#endif

        .add( mesh_options( prefix ) )
        /* arm options */
        .add( arm_options( prefix ) )
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

        /* onelab options */
        .add( onelab_options( prefix ) )

        /* function space options */
#if !defined( FEELPP_HAS_TRILINOS_EPETRA )
        .add( functionspace_options( prefix ) )
#endif
        .add( aitken_options( prefix ) )

        .add( msi_options(prefix) )
        .add( fit_options(prefix) )
        .add( checker_options(prefix) )
        .add( journal_options(prefix) )
        .add( fmu_options(prefix) )
        ;

    return opt;

}
}
