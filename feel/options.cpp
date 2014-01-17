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
    ( prefixvm( prefix, "bdf.order" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "order in time" )
    ( prefixvm( prefix, "bdf.strategy" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "strategy, 0=constant time steps, 1=adaptive time steps" )
    ( prefixvm( prefix, "bdf.steady" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "false: unsteady, true:steady" )
    ( prefixvm( prefix, "bdf.restart" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart " )
    ( prefixvm( prefix, "bdf.restart.path" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "path where we reload old data" )
    ( prefixvm( prefix, "bdf.restart.at-last-save" ).c_str(), Feel::po::value<bool>()->default_value( false ), "do a restart with ti the last save " )
    ( prefixvm( prefix, "bdf.iterations-between-order-change" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "iteration between time order change" )
    ( prefixvm( prefix, "bdf.save" ).c_str(), Feel::po::value<bool>()->default_value( true ), "save elements in file " )
    ( prefixvm( prefix, "bdf.save.freq" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "freq for save elements in file " )
    ( prefixvm( prefix, "bdf.rank-proc-in-files-name" ).c_str(), Feel::po::value<bool>()->default_value( false ), "the name of files generated has the rank of the processor automatically if true" )
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
        ( "eim.rebuild-database" , Feel::po::value<bool>()->default_value( 0 ), "rebuild database (if it already exists)" )
        ( "eim.cvg-study" , Feel::po::value<bool>()->default_value( 0 ), "for convergence study" )
        ( "eim.compute-error-with-truth-expression" , Feel::po::value<bool>()->default_value( true ), "compute the error with the truth expression ( not its projection ) if true" )
        ( "eim.use-dimension-max-functions" , Feel::po::value<bool>()->default_value( 0 ), "force to use dimension-max basis functions" )
        ( "eim.computational-time-neval",Feel::po::value<int>()->default_value( 0 )," number of evaluation to perform to have the computational time of eim online step" )
        ( "eim.compute-expansion-of-expression",Feel::po::value<bool>()->default_value( false )," use true expression if true, else use the projection of the expression" )
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
    po::options_description _options( "Solver EigenValue Slepc -- " + prefix + " solver options" );
    _options.add_options()
    // solver options
    ( ( _prefix+"solvereigen.solver-type" ).c_str(), Feel::po::value<int>()->default_value( KRYLOVSCHUR ), "type of eigenvalue solver" )
    ( ( _prefix+"solvereigen.problem-type" ).c_str(), Feel::po::value<int>()->default_value( NHEP ), "type of eigenvalue problem" )
    ( ( _prefix+"solvereigen.position" ).c_str(), Feel::po::value<int>()->default_value( LARGEST_MAGNITUDE ), "eigenvalue solver position in spectrum: LARGEST_MAGNITUDE=0, SMALLEST_MAGNITUDE=1, LARGEST_REAL=2, SMALLEST_REAL=3, LARGEST_IMAGINARY=4, SMALLEST_IMAGINARY=5" )
    ( ( _prefix+"solvereigen.nev" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "number of requested eigenpairs" )
    ( ( _prefix+"solvereigen.ncv" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "number of basis vectors" )
    ( ( _prefix+"solvereigen.tol" ).c_str(), Feel::po::value<double>()->default_value( 1e-10 ), "solver tolerance" )
    ( ( _prefix+"solvereigen.maxiter" ).c_str(), Feel::po::value<int>()->default_value( 10000 ), "maximum number of iterations" );

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
    ( "crb.scm.use-scm",Feel::po::value<bool>()->default_value( true ), "use scm if true")
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
    ( "crb.dimension-max"   , Feel::po::value<int>()->default_value( 50 ),       "Offline  max WN size" )
    ( "crb.dimension"   , Feel::po::value<int>()->default_value( -1 ),       "Online  WN size" )
    ( "crb.error-type"   , Feel::po::value<int>()->default_value( ( int )CRB_RESIDUAL ),       "CRB error type to be computed" )
    ( "crb.factor"   , Feel::po::value<int>()->default_value( -1 ),  "factor useful to estimate error by empirical method" )
    ( "crb.Nm"   , Feel::po::value<int>()->default_value( 1 ),       "Offline  number of modes per mu (for the POD) " )
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

    ( "crb.is-model-executed-in-steady-mode",Feel::po::value<bool>()->default_value( true ), "true if model is executed in steady mode, else turn it to false")
    ( "crb.use-ginac-for-beta-expressions",Feel::po::value<bool>()->default_value( false ), "use ginac to compute expression of beta coefficients if true")
    ( "crb.use-linear-model",Feel::po::value<bool>()->default_value( false ), "do not iterate in fixed point if true")

    ( "crb.use-predefined-WNmu",Feel::po::value<bool>()->default_value( false ), "read parameters to take for the offline step from a file named SamplingWNmu if true")
    ( "crb.use-predefined-test-sampling",Feel::po::value<bool>()->default_value( false ), "read parameters from file named SamplingForTest if true to run the test")
    ( "crb.use-logEquidistributed-WNmu",Feel::po::value<int>()->default_value( 0 ), "parameters are log-equidistributed for the offline step (the value indicates the number of parameters)")
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
    ;

    crboptions
        .add( crbSCMOptions() );

    return crboptions;
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
