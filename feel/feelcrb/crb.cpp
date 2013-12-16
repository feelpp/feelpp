/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-24

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file opuscrb.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-24
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbscm.hpp>
#include <feel/feelcrb/eim.hpp>

namespace Feel
{
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
}

