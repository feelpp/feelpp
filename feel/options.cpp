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
        ( "directory", po::value<std::string>(), "change directory to specified one" )
        ( "npdir", po::value<bool>()->default_value(true), "enable/disable sub-directory np_<number of processors>")
        ( "fail-on-unknown-option", po::value<bool>()->default_value(false), "exit feel++ application if unknown option found" )
        ( "show-preconditioner-options", "show on the fly the preconditioner options used" )
        ( "serialization-library", po::value<std::string>()->default_value("boost"), "Library used for serialization" )
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
gmsh_options( std::string const& prefix )
{
    po::options_description _options( "Gmsh " + prefix + " options" );
    _options.add_options()
    // solver options
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
        ( prefixvm( prefix,"gmsh.partition" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Partition Gmsh mesh once generated or loaded" )
        ( prefixvm( prefix,"gmsh.respect_partition" ).c_str(), Feel::po::value<bool>()->default_value( false ), "true to respect paritioning when mesh is loaded, false to ensure that partition is within the number of processors" )
        ( prefixvm( prefix,"gmsh.npartitions" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "Number of partitions" )
#if defined(HAVE_METIS)
        ( prefixvm( prefix,"gmsh.partitioner" ).c_str(), Feel::po::value<int>()->default_value( GMSH_PARTITIONER_DEFAULT ), "Gmsh partitioner (1=CHACO, 2=METIS)" )
#else
        ( prefixvm( prefix,"gmsh.partitioner" ).c_str(), Feel::po::value<int>()->default_value( GMSH_PARTITIONER_DEFAULT ), "Gmsh partitioner (1=CHACO)" )
#endif
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


    return _options;

}

po::options_description
on_options( std::string const& prefix )
{
    po::options_description _options( "Dirichlet treatment options " + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix,"on.type" ).c_str(), Feel::po::value<std::string>()->default_value( "elimination" ), "Strong Dirichlet conditions treatment type: elimination, elimination_keep_diagonal, elimination_symmetric, elimination_symmetric_keep_diagonal, penalisation" )
        ( prefixvm( prefix,"on.verbose" ).c_str(), Feel::po::value<bool>()->default_value( false ), "print in logfiles information about Dirichlet conditions treatment" )
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
    ( prefixvm( prefix, "bdf.restart" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart " )
    ( prefixvm( prefix, "bdf.restart.path" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "path where we reload old data" )
    ( prefixvm( prefix, "bdf.restart.at-last-save" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart with ti the last save " )
    ( prefixvm( prefix, "bdf.restart.step-before-last-save" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "do a restart with ti the ieme step before last save " )
    ( prefixvm( prefix, "bdf.iterations-between-order-change" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "iteration between time order change" )
    ( prefixvm( prefix, "bdf.save" ).c_str(), Feel::po::value<bool>()->default_value( true ), "save elements in file " )
    ( prefixvm( prefix, "bdf.save.freq" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "freq for save elements in file " )
    ( prefixvm( prefix, "bdf.rank-proc-in-files-name" ).c_str(), Feel::po::value<bool>()->default_value( false ), "the name of files generated has the rank of the processor automatically if true" )
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
    ( prefixvm( prefix, "ts.restart" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart " )
    ( prefixvm( prefix, "ts.restart.path" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "path where we reload old data" )
    ( prefixvm( prefix, "ts.restart.at-last-save" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart with ti the last save " )
    ( prefixvm( prefix, "ts.restart.step-before-last-save" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "do a restart with ti the ieme step before last save " )
    ( prefixvm( prefix, "ts.save" ).c_str(), Feel::po::value<bool>()->default_value( true ), "save elements in file " )
    ( prefixvm( prefix, "ts.save.freq" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "freq for save elements in file " )
    ( prefixvm( prefix, "ts.rank-proc-in-files-name" ).c_str(), Feel::po::value<bool>()->default_value( false ), "the name of files generated has the rank of the processor automatically if true" )
    ;
    return _options;
}


Feel::po::options_description
eimOptions( std::string const& prefix )
{
    Feel::po::options_description eimoptions( "EIM Options" );
    eimoptions.add_options()
        ( "eim.sampling-size"   , Feel::po::value<int>()->default_value( 30 ), "Offline  sampling size " )
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
        ;

    return eimoptions;
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
        ( ( _prefix+"solvereigen.ncv" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "number of basis vectors" )
        ( ( _prefix+"solvereigen.tolerance" ).c_str(), Feel::po::value<double>()->default_value( 1e-10 ), "solver tolerance" )
        ( ( _prefix+"solvereigen.maxiter" ).c_str(), Feel::po::value<int>()->default_value( 10000 ), "maximum number of iterations" )
        ( ( _prefix+"solvereigen.verbose" ).c_str(), Feel::po::value<bool>()->default_value( false ), "verbose eigen solver" );

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
    ( "crb.output-index"   , Feel::po::value<int>()->default_value( 0 ), "output index (default is right hand side = 0)" )
    ( "crb.sampling-size"   , Feel::po::value<int>()->default_value( 1000 ), "Offline  sampling size " )
    ( "crb.sampling-mode"   , Feel::po::value<std::string>()->default_value( "log-random" ), "Offline  sampling mode, equidistribute, log-equidistribute or log-random " )
    ( "crb.all-procs-have-same-sampling" , Feel::po::value<bool>()->default_value( "false" ), "all procs have the same sampling if true" )
    ( "crb.error-max"   , Feel::po::value<double>()->default_value( 1e-6 ),       "Offline  tolerance" )
    ( "crb.online-tolerance"   , Feel::po::value<double>()->default_value( 1e-2 ),       "Online  tolerance" )
    ( "crb.dimension-max"   , Feel::po::value<int>()->default_value( 1 ),       "Offline max WN size, set to 1 by default to avoid to enrich the existing database if this option doesn't appear in onefeel interface or in the config file." )
    ( "crb.dimension"   , Feel::po::value<int>()->default_value( -1 ),       "Online  WN size" )
    ( "crb.error-type"   , Feel::po::value<int>()->default_value( ( int )CRB_RESIDUAL_SCM ),       "CRB error type to be computed" )
    ( "crb.compute-apee-for-each-time-step",Feel::po::value<bool>()->default_value( true ),"Compute error estimation for each time step (parabolic problems) is true, else compute only for the last one")
    ( "crb.factor"   , Feel::po::value<int>()->default_value( -1 ),  "factor useful to estimate error by empirical method" )
    ( "crb.Nm"   , Feel::po::value<int>()->default_value( 1 ),       "Offline  number of modes per mu (for the POD) " )
    ( "crb.apply-POD-to-WN"   , Feel::po::value<bool>()->default_value( false ), "apply a POD on approximation functions spaces (primal and dual) if true and if deal with a transient problem " )
    ( "crb.check.rb"   , Feel::po::value<int>()->default_value( 0 ),       "check reduced basis" )
        //( "crb.check.gs"   , Feel::po::value<int>()->default_value( 0 ),       "check Gram-Schmidt orthonormalisation" )
    ( "crb.orthonormality-tol" , Feel::po::value<double>()->default_value( 1e-13 ),"tolerance of orthonormalisation : i.e. norm of matrix A(i,j)=scalarProduct( Wn[j], Wn[i] )" )
    ( "crb.orthonormality-max-iter" , Feel::po::value<int>()->default_value( 10 ),"while the tolerance is not reached, the orthonormalization step is done or until max-iter is reached" )

    ( "crb.check.residual"   , Feel::po::value<bool>()->default_value( false ),  "check residual" )
    ( "crb.reuse-prec"   , Feel::po::value<bool>()->default_value( 0 ),       "reuse or not the preconditioner" )
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
    ( "crb.compute-variance" , Feel::po::value<bool>()->default_value( 0 ), " if true the output is the variance and not l(v)" )
    ( "crb.save-information-for-variance",Feel::po::value<bool>()->default_value( 0 ), "if true will build variance matrix but it takes some times" )

    ( "crb.use-newton",Feel::po::value<bool>()->default_value( false ), "use newton algorithm (need to provide a jacobian and a residual)" )
    ( "crb.max-fixedpoint-iterations",Feel::po::value<int>()->default_value( 10 ), "nb iteration max for the fixed point (online part)" )
    ( "crb.increment-fixedpoint-tol",Feel::po::value<double>()->default_value( 1e-10 ), "tolerance on solution for fixed point (online part)" )
    ( "crb.output-fixedpoint-tol",Feel::po::value<double>()->default_value( 1e-10 ), "tolerance on output for fixed point (online part)" )
    ( "crb.fixedpoint-verbose",Feel::po::value<bool>()->default_value( false ), "fixed point verbose if true" )
    ( "crb.fixedpoint-critical-value",Feel::po::value<double>()->default_value(1000 ), "will crash if increment error at the end of fixed point is greater than this critical value" )
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

    ( "crb.user-parameters",Feel::po::value<std::string>()->default_value( "" ), "values of parameters (used for one feel)")
    ( "crb.vary-only-parameter-components",Feel::po::value<std::string>()->default_value( "" ), "specify which parameter component vary (max : 2 components + time 't') and how many values we take in each direction. For example 0 10 1 20 means that component 0 will take 10 values and component 1 will take 20 values. We can write also t 1 10 and in this case the time will vary and also component 1")
    ( "crb.load-elements-database",Feel::po::value<bool>()->default_value( false ), "load database of elements if true, need to be true for visualization, need to be false to run CRB approximation on a different number of processors than this was used to build the reduced basis ")

    ( "crb.solve-fem-monolithic",Feel::po::value<bool>()->default_value( false ), "solve FEM problem without using EIM and without affine decomposition ")
    ( "crb.export-name-max-size",Feel::po::value<int>()->default_value( 30 ), "maximum size for variable names in export (truncature)")
    ;

    crboptions
        .add( crbSCMOptions() );

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
        // options for pmm
        ( prefixvm( prefix, "blockns.pmm.diag" ).c_str(), Feel::po::value<bool>()->default_value(1), "set to true to use diagonal of the pressure mass matrix, false otherwise" )
        ;
    return _options;
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
        ( prefixvm( prefix,"exporter.vtk.insitu.enable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Enable In-situ visualization with VTK exporter (Data won't be written to disk any longer)." )
        ( prefixvm( prefix,"exporter.vtk.insitu.pyscript" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "Specify a python user script for visualization." )
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

        /* nonlinear solver options */
        .add( nlsolver_options() )

        /* discr options */
        .add( ts_options( prefix ) )
        .add( bdf_options( prefix ) )

        /* exporter options */
        .add( exporter_options( prefix ) )

        /* nlopt options */
#if defined(FEELPP_HAS_NLOPT)
        .add( nlopt_options( prefix ) )
#endif

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
        ;

    return opt;

}
}
