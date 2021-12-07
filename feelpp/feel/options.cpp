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
        ( "case", po::value<std::string>(), "specify a case directory" )
        ( "case.config-file", po::value<std::string>(), "specify the config-file in the case directory" )
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
        ( "dirs", "list standard feelpp directories" )
        ( "directory", po::value<std::string>(), "change directory to specified one" )
        ( "repository.prefix", po::value<std::string>(), "change directory to specified one" )
        ( "repository.case", po::value<std::string>(), "change directory to specified one relative to repository.prefix" )
        ( "repository.npdir", po::value<bool>()->default_value(true), "enable/disable sub-directory np_<number of processors>")
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
case_options( int default_dim, std::string const& default_discr, std::string const& prefix )
{
    po::options_description file( "Case options" );
    file.add_options()
        ( prefixvm( prefix, "case.dim").c_str(), po::value<int>()->default_value(default_dim), "case dimenstion" )
        ( prefixvm( prefix, "case.discretization").c_str(), po::value<std::string>()->default_value(default_discr), "case discretization" )
        ;
    return file;
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
    /*_options.add_options()
        ( prefixvm( prefix,"x" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "x coordinate value " )
        ( prefixvm( prefix,"y" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "y coordinate value " )
     ( prefixvm( prefix,"z" ).c_str(), Feel::po::value<double>()->default_value( 0 ), "z coordinate value " );*/
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
glpk_options( std::string const& prefix )
{
    po::options_description _options( "Glpk " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"glpk.verbosity").c_str(), po::value<int>()->default_value(2), "level of verbosity (between 1 and 4)" )
        ( prefixvm( prefix,"glpk.method").c_str(), po::value<int>()->default_value(1), "method (1=primal,2=dualprimal,3=dual)" )
        ( prefixvm( prefix,"glpk.tolerance-bounds").c_str(), po::value<double>()->default_value(1e-7), "Tolerance used to check if the basic solution is primal feasible" )
        ( prefixvm( prefix,"glpk.tolerance-dual").c_str(), po::value<double>()->default_value(1e-7), "Tolerance used to check if the basic solution is dual feasible" )
        ( prefixvm( prefix,"glpk.tolerance-pivot").c_str(), po::value<double>()->default_value(1e-10), "Tolerance used to choose eligble pivotal elements of the simplex table" )
        ( prefixvm( prefix,"glpk.iteration-limit").c_str(), po::value<int>()->default_value(INT_MAX), "Simplex iteration limit" )
        ( prefixvm( prefix,"glpk.time-limit").c_str(), po::value<int>()->default_value(INT_MAX), "Searching time limit, in milliseconds" )
        ( prefixvm( prefix,"glpk.presolve").c_str(), po::value<int>()->default_value(0), "enable presolver (0 or 1)" )
        ( prefixvm( prefix,"glpk.scaling").c_str(), po::value<int>()->default_value(128), "scaling options (1=geometric scaling, 16=equilibration scaling, 32=round scale factors to power of two, 64=skip scaling, 0x80=choose automatically)")
        ;
    return _options;
}

po::options_description
mesh_options( std::string const& prefix )
{
    po::options_description _options( "Mesh " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"mesh.filename").c_str(), po::value<std::string>(), "mesh filename" )
        ( prefixvm( prefix,"mesh.partition.enable").c_str(), po::value<bool>()->default_value(0), "partition the mesh using Feel++ partitioners" )
        ( prefixvm( prefix,"mesh.partition.size").c_str(), po::value<int>()->default_value(1), "number of partitions" )
        ( prefixvm( prefix,"mesh.partition.type").c_str(), po::value<std::string>()->default_value("metis"), "mesh partitioner: metis, (more to come)" )
        ( prefixvm( prefix,"mesh.save.enable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "enable saving mesh to disk" )
        ( prefixvm( prefix,"mesh.save.formats" ).c_str(), Feel::po::value<std::vector<std::string>>()->default_value( {"json+h5","msh"} ), "format of the mesh: json+h5, msh" )
        ( prefixvm( prefix,"mesh.load.enable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "enable loading mesh from disk, overriding file name extension" )
        ( prefixvm( prefix,"mesh.load.format" ).c_str(), Feel::po::value<std::string>()->default_value( "json+h5" ), "file format to load: msh, json" )
        ( prefixvm( prefix,"mesh.scale" ).c_str(), Feel::po::value<double>()->default_value( 1 ), "scale the mesh after loading" );

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
        ( prefixvm( prefix,"gmsh.hsize2" ).c_str(), Feel::po::value<double>()->default_value( 0.1 ), "characteristic mesh size" )
        ( prefixvm( prefix,"gmsh.geo-variables-list" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "modify a list of geo variables (ex : alpha=1:beta=2)" )
        ( prefixvm( prefix,"gmsh.save" ).c_str(), Feel::po::value<bool>()->default_value( true ), "save msh file to disk once generated" )
        ( prefixvm( prefix,"gmsh.straighten" ).c_str(), Feel::po::value<bool>()->default_value( true ), "straighten high order mesh" )
        ( prefixvm( prefix,"gmsh.structured" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "generated a structured mesh" )
        ( prefixvm( prefix,"gmsh.rebuild" ).c_str(), Feel::po::value<bool>()->default_value( true ), "force rebuild msh file from geo file" )
        ( prefixvm( prefix,"gmsh.refine" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "refinement by splitting level" )
        ( prefixvm( prefix,"gmsh.physical_are_elementary_regions" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Physical regions are defined by elementary regions, useful for medit format" )
        ( prefixvm( prefix,"gmsh.savehdf5" ).c_str(), Feel::po::value<bool>()->default_value( false ), "save msh file to disk once generated in HDF5 format" )
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

po::options_description remesh_options( std::string const& prefix = "" );

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
#if defined(FEELPP_HAS_HDF5)
    ( prefixvm( prefix, "bdf.file-format" ).c_str(), Feel::po::value<std::string>()->default_value( "hdf5" ), "save elements in file " )
#else
    ( prefixvm( prefix, "bdf.file-format" ).c_str(), Feel::po::value<std::string>()->default_value( "binary" ), "save elements in file " )
#endif
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
        ( prefixvm( prefix, "ts.adaptive" ).c_str(), Feel::po::value<bool>()->default_value( false ), "enable/disable adaptive time step if available" ) 
        ( prefixvm( prefix, "ts.tol" ).c_str(), Feel::po::value<double>()->default_value( 1e-3 ), "time step tolerance" )
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
#if defined(FEELPP_HAS_HDF5)
        ( prefixvm( prefix, "ts.file-format" ).c_str(), Feel::po::value<std::string>()->default_value( "hdf5" ), "save elements in file " )
#else
        ( prefixvm( prefix, "ts.file-format" ).c_str(), Feel::po::value<std::string>()->default_value( "binary" ), "save elements in file " )
#endif
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
        //( prefixvm( prefix, "sc.backend" ).c_str(), po::value<std::string>()->default_value( "petsc" ), "enable/disable static condensation" )
        ( prefixvm( prefix, "sc.condense.parallel" ).c_str(), po::value<bool>()->default_value( true ), "enable/disable parallel condense in static condensation" )
        ( prefixvm( prefix, "sc.condense.grain" ).c_str(), po::value<int>()->default_value( 100 ), "grain size for parallel local solve in static condensation" )
        ( prefixvm( prefix, "sc.condense.parallel.n" ).c_str(), po::value<int>()->default_value( 2 ), "number of tasks for parallel local solve in static condensation" )
        ( prefixvm( prefix, "sc.localsolve.parallel" ).c_str(), po::value<bool>()->default_value( true ), "enable/disable parallel local solve in static condensation" )
        ( prefixvm( prefix, "sc.localsolve.grain" ).c_str(), po::value<int>()->default_value( 100 ), "grain size for parallel local solve in static condensation" )
        ( prefixvm( prefix, "sc.localsolve.parallel.n" ).c_str(), po::value<int>()->default_value( 2 ), "number of tasks for parallel local solve in static condensation" )
        ( prefixvm( prefix, "sc.ibc_partitioning" ).c_str(), po::value<bool>()->default_value( false ), "enable/disable special partitioning with IBC: all faces with IBC imposed with the same pid" )
        ( prefixvm( prefix, "sc.ibc_partitioning.marker" ).c_str(), po::value<std::string>()->default_value( "Ibc" ), "marker of the IBC to apply the special partitioning" )
        ;

    return _options.add( backend_options( prefixvm(prefix,"sc") ) ).add( backend_options( prefixvm(prefix,"sc.post") ) );
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
        ( ( _prefix+"solvereigen.spectrum" ).c_str(), Feel::po::value<std::string>()->default_value( "largest_magnitude" ), "eigenvalue solver position in spectrum. Choice: largest_magnitude, smallest_magnitude, largest_real, smallest_real, largest_imaginary, smallest_imaginary, target_magnitude, target_real, target_imaginary" )
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

po::options_description
pcd_options( std::string const& prefix )
{
    po::options_description _options( "PCD options (" + prefix + ")" );
    _options.add_options()
        ( prefixvm( prefix, "pcd.bc-type-with-Dirichlet" ).c_str(), Feel::po::value<std::string>()->default_value("Robin"), "Type of boundary conditions with Dirichlet in NS : Robin or Dirichlet" )
        ( prefixvm( prefix, "pcd.bc-type-with-Neumann" ).c_str(), Feel::po::value<std::string>()->default_value("Dirichlet"), "Type of boundary conditions with Neumann in NS : Neumann or Dirichlet" )
        ( prefixvm( prefix, "pcd.order" ).c_str(), Feel::po::value<int>()->default_value(1), "order for pcd operator 1:Ap^-1 Fp Mp^-1 other: Mp^-1 Fp Ap^-1" )
        ( prefixvm( prefix, "pcd.diffusion" ).c_str(), Feel::po::value<std::string>()->default_value("Laplacian"), "Laplacian or BTBt" )
        ( prefixvm( prefix, "pcd.diffusion.weakdir" ).c_str(), Feel::po::value<bool>()->default_value(0), "set to true for Weak dirichlet conditions for Fp and Ap, false otherwise" )
        ( prefixvm( prefix, "pcd.diffusion.weakdir.penaldir" ).c_str(), Feel::po::value<double>()->default_value(10.), "Penalisation parameter for weak bc" )
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
        ( prefixvm( prefix,"checker.gradient" ).c_str(), Feel::po::value<std::string>()->default_value(""), "solution gradient against which to check" )
        ( prefixvm( prefix,"checker.tolerance.exact" ).c_str(), Feel::po::value<double>()->default_value(1e-15), "tolerance for numerical exact solution check" )
        ( prefixvm( prefix,"checker.tolerance.order" ).c_str(), Feel::po::value<double>()->default_value(1e-1), "tolerance for order check" )
        ( prefixvm( prefix,"checker.name" ).c_str(), Feel::po::value<std::string>()->default_value("checker"), "name of the test" )
        ( prefixvm( prefix,"checker.filename" ).c_str(), Feel::po::value<std::string>()->default_value("checker.json"), "name of the test" )
        ( prefixvm( prefix,"checker.script" ).c_str(), Feel::po::value<std::string>(), "script to run" )
        ( prefixvm( prefix,"checker.compute-pde-coefficients" ).c_str(), Feel::po::value<bool>()->default_value(true), "run script" )
        ( prefixvm( prefix,"checker.verbose" ).c_str(), Feel::po::value<bool>()->default_value(false), "checker prints more information" )
        ;
    return _options;
}

po::options_description
journal_options( std::string const& prefix )
{
    po::options_description options( "Journal " + prefix + " options" );
    options.add_options()
        ( prefixvm( prefix,"journal" ).c_str(), Feel::po::value<bool>()->default_value(true), "Enable simulation info" )
        ( prefixvm( prefix,"journal.filename" ).c_str(), Feel::po::value<std::string>(), "name of the json file" )
        ( prefixvm( prefix,"journal.auto" ).c_str(), Feel::po::value<bool>()->default_value(true), "Generate the journal automatically based on the code." )
        ( prefixvm( prefix,"journal.database" ).c_str(), Feel::po::value<bool>()->default_value(false), "Enable database usage" )
        ( prefixvm( prefix,"journal.database.name" ).c_str(), Feel::po::value<std::string>()->default_value("feelpp"), "Database name (MongoDB)" )
        ( prefixvm( prefix,"journal.database.host" ).c_str(), Feel::po::value<std::string>()->default_value("localhost"), "Database host or IP (MongoDB)" )
        ( prefixvm( prefix,"journal.database.port" ).c_str(), Feel::po::value<std::string>()->default_value("27017"), "Database port (MongoDB)" )
        ( prefixvm( prefix,"journal.database.user" ).c_str(), Feel::po::value<std::string>()->default_value(""), "Database username (MongoDB)" )
        ( prefixvm( prefix,"journal.database.password" ).c_str(), Feel::po::value<std::string>()->default_value(""), "Database password (MongoDB)" )
        ( prefixvm( prefix,"journal.database.authsrc" ).c_str(), Feel::po::value<std::string>()->default_value("admin"), "Database authenticate source collection [admin] (MongoDB)" )
        ( prefixvm( prefix,"journal.database.collection" ).c_str(), Feel::po::value<std::string>()->default_value("journal"), "Database journal collection [journal] (MongoDB)" )
        ;
    return options;
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
ptree_options( std::string const& prefix )
{
    po::options_description _options( "Ptree " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"json-editions" ).c_str(), po::value<std::vector<std::string> >()->multitoken(), "specify a list of entries to modified in json. format= key:value " )
        ;
    return _options;
}

po::options_description
eq_options( std::string const& prefix )
{
    po::options_description _options( "EQ " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"eq.verbosity" ).c_str(), po::value<int>()->default_value(0), "level of verbosity" )
        ( prefixvm( prefix,"eq.db.load").c_str(), po::value<bool>()->default_value(true), "load the database" )
        ( prefixvm( prefix,"eq.db.filename").c_str(), po::value<std::string>()->default_value("eq"), "directory in which save/load the database" )
        ( prefixvm( prefix,"eq.tolerance" ).c_str(), po::value<double>()->default_value(1e-10), "tolerance to reach for offline empirical quadrature" )
        ( prefixvm( prefix, "eq.sampling-size").c_str(), po::value<int>()->default_value(100), "sampling size for offline empirical quadrature" )
        ( prefixvm( prefix, "eq.max-order").c_str(), po::value<int>()->default_value(30), "maximum order of quadrature" )
        ( prefixvm( prefix, "eq.order").c_str(), po::value<int>()->default_value(-1), "force order of quadrature" )
        ( prefixvm( prefix, "eq.tolerance-zero").c_str(), po::value<double>()->default_value(1e-12), "tolerance to decide if weight is zero or not")
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
                   .add( backend_options( "Ap" ) )
                   .add( backend_options( "Fp" ) )
                   .add( backend_options( "Mp" ) )
                   .add( backend_options( "Fu" ) )
                   .add( backend_options( "Bt" ) )
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
        /* glpk options */
#if defined(FEELPP_HAS_GLPK_H)
                   .add( glpk_options( prefix ) )
#endif

                   .add( mesh_options( prefix ) )
                   /* arm options */
                   .add( arm_options( prefix ) )
                   /* gmsh options */
                   .add( gmsh_options( prefix ) )
                   /* remeshing options */
                   .add( remesh_options( prefix ) )
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

                   .add( msi_options( prefix ) )
                   .add( fit_options( prefix ) )
                   .add( checker_options( prefix ) )
                   .add( journal_options( prefix ) )
                   .add( fmu_options( prefix ) )
                   .add( ptree_options( prefix ) );

    return opt;

}
}
