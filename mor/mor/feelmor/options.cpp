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
#include <fmt/core.h>
#include <boost/algorithm/string.hpp>
#include <feel/feelcrb/crbenums.hpp>
#include <feel/options.hpp>

namespace Feel
{
std::string
prefixvm( std::string const& prefix,
          std::string const& opt,
          std::string const& sep = "." );

Feel::po::options_description
geim_options( std::string const& prefix )
{
    Feel::po::options_description geimoptions( "GEIM Options" );
    geimoptions.add_options()
        ( "geim.dimension-max", po::value<int>()->default_value(10), "maximum number of basis" )
        ( "geim.tolerance", po::value<double>()->default_value(1e-10), "tolerance" )
        ( "geim.rebuild-database", po::value<bool>()->default_value(false), "rebuild the database" )
        ( "geim.db.load", po::value<int>()->default_value(2), "=0 use db.filename, =1 use last DB created =2 use last DB modified =3 use db.id =4 create new db" )
        ( "geim.db.filename", po::value<std::string>()->default_value(""), "path to the db when db.load or db.update = 0" )
        ( "geim.db.id", po::value<std::string>()->default_value(""), "id of the db when db.load or db.update = 3" )
        ;
    return geimoptions.add( backend_options( prefixvm(prefix,"geim") ) );
}

Feel::po::options_description
eimOptions( std::string const& prefix )
{
    Feel::po::options_description eimoptions( "EIM Options" );
    eimoptions.add_options()
        ( prefixvm( prefix, "eim.sampling-size").c_str()   , Feel::po::value<int>()->default_value( 30 ), "Offline  sampling size " )
        ( prefixvm( prefix, "eim.sampling-mode").c_str()   , Feel::po::value<std::string>()->default_value( "random" ), "EIM Offline : random, log-random, log-equidistribute, equidistribute " )
        ( prefixvm( prefix, "eim.error-max").c_str()   , Feel::po::value<double>()->default_value( 1e-6 ),       "Offline  tolerance" )
        ( prefixvm( prefix, "eim.online-tolerance").c_str()   , Feel::po::value<double>()->default_value( 1e-2 ),       "Online  tolerance" )
        ( prefixvm( prefix, "eim.dimension-max").c_str()   , Feel::po::value<int>()->default_value( 50 ),       "Offline  max WN size" )
        //( prefixvm( prefix, "eim.error-type").c_str()   , Feel::po::value<int>()->default_value( ( int )EIM_RESIDUAL ),       "EIM error type to be computed" )
        ( prefixvm( prefix, "eim.check.rb").c_str()   , Feel::po::value<int>()->default_value( 0 ),       "check reduced basis" )
        ( prefixvm( prefix, "eim.norm-used-for-residual").c_str() , Feel::po::value<std::string>()->default_value( "Linfty" ), "name of the norm used : Linfty - L2 - LinftyVec" )
        ( prefixvm( prefix, "eim.check.residual").c_str()   , Feel::po::value<int>()->default_value( 0 ),  "check residual" )
        ( prefixvm( prefix, "eim.reuse-prec").c_str()   , Feel::po::value<bool>()->default_value( 0 ),       "reuse or not the preconditioner" )
        ( prefixvm( prefix, "eim.rebuild-database").c_str() , Feel::po::value<bool>()->default_value( false ), "rebuild database (if it already exists)" )
        ( prefixvm( prefix, "eim.enrich-database").c_str() , Feel::po::value<bool>()->default_value( true ), "enrich the database (if it already exists)" )
        ( prefixvm( prefix, "eim.cvg-study").c_str() , Feel::po::value<bool>()->default_value( 0 ), "for convergence study" )
        ( prefixvm( prefix, "eim.compute-error-with-truth-expression").c_str() , Feel::po::value<bool>()->default_value( false ), "compute the error with the truth expression ( not its projection ) if true" )
        ( prefixvm( prefix, "eim.use-dimension-max-functions").c_str() , Feel::po::value<bool>()->default_value( 0 ), "force to use dimension-max basis functions" )
        ( prefixvm( prefix, "eim.computational-time-neval").c_str(),Feel::po::value<int>()->default_value( 0 )," number of evaluation to perform to have the computational time of eim online step" )
        ( prefixvm( prefix, "eim.compute-expansion-of-expression").c_str(),Feel::po::value<bool>()->default_value( false )," use true expression if true, else use the projection of the expression" )
        ( prefixvm( prefix, "eim.show-mu-selection").c_str(),Feel::po::value<bool>()->default_value( false )," print list of parameters selected during offline step" )
        ( prefixvm( prefix, "eim.show-t-selection").c_str(),Feel::po::value<bool>()->default_value( false )," print list of interpolation points selected during offline step" )
        ( prefixvm( prefix, "eim.show-offline-error").c_str(),Feel::po::value<bool>()->default_value( false )," print list of error associated to mu selected during offline step" )

        ( prefixvm( prefix, "eim.elements.write").c_str(), Feel::po::value<bool>()->default_value( true ), "Write evaluated nl solutions on disk"  )
        ( prefixvm( prefix, "eim.elements.directory").c_str(), Feel::po::value<std::string>()->default_value( "nlsolutions" ), "directory were nl solutions are stored"  )
        ( prefixvm( prefix, "eim.elements.clean-directory").c_str(), Feel::po::value<bool>()->default_value( false ), ""  )
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
        ;

    return deimoptions.add(backend_options(prefixvm(prefix,"deim-online")));
}

Feel::po::options_description
pbdw_options( std::string const& prefix )
{
    Feel::po::options_description pbdwoptions( "PBDW Options" );
    pbdwoptions.add_options()
        ( "pbdw.rebuild-database", po::value<bool>()->default_value(false), "rebuild the database" )
        ( "pbdw.db.load", po::value<int>()->default_value(2), "=0 use db.filename, =1 use last DB created =2 use last DB modified =3 use db.id =4 create new db" )
        ( "pbdw.db.filename", po::value<std::string>()->default_value(""), "path to the db when db.load or db.update = 0" )
        ( "pbdw.db.id", po::value<std::string>()->default_value(""), "id of the db when db.load or db.update = 3" )
        ;
    return pbdwoptions.add( backend_options( prefixvm(prefix,"pbdw")) );
}

Feel::po::options_description
crbSEROptions( std::string const& prefix )
{
    Feel::po::options_description seroptions( "SER Options" );
    seroptions.add_options()
        ( prefixvm( prefix, "ser.rb-frequency").c_str(), Feel::po::value<int>()->default_value( 0 ), "Number of RB basis built per step of SER process, 0 : no SER" )
        ( prefixvm( prefix, "ser.eim-frequency").c_str(), Feel::po::value<int>()->default_value( 0 ), "Number of EIM basis built per step of SER process, 0 : no SER")
        ( prefixvm( prefix, "ser.use-rb-in-eim-mu-selection").c_str(), Feel::po::value<bool>()->default_value( false ), "Use RB approx. to select parameters during EIM offline step")
        ( prefixvm( prefix, "ser.use-rb-in-eim-basis-build").c_str(), Feel::po::value<bool>()->default_value( false ), "Use RB approx. to build EIM basis")
        ( prefixvm( prefix, "ser.rb-rebuild-freq").c_str(), Feel::po::value<int>()->default_value( -1 ), "rebuild database with a given frequency" )
        ( prefixvm( prefix, "ser.error-estimation").c_str(), Feel::po::value<bool>()->default_value( false ), "Use SER error estimation = Norm of residual (Riesz) as error indicator")
        ( prefixvm( prefix, "ser.use-greedy-in-rb").c_str(), Feel::po::value<bool>()->default_value( false ), "Use SER error indicator to build RB basis (Greedy)")
        ( prefixvm( prefix, "ser.eim-subtrainset-method").c_str(), Feel::po::value<int>()->default_value( 0 ), "Build subtrainset from error estimation : full trainset (=0), select mu for which residual < tol (=1), select mu for which residual > tol (=2)")
        ( prefixvm( prefix, "ser.print-rb-iterations_info").c_str(), Feel::po::value<bool>()->default_value( false ), "Write online rb iterations (Picard) needed for Greedy")
        ( prefixvm( prefix, "ser.radapt-eim-rtol").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "Relative tolerance criterion for EIM r-adaptation - default(0.0) = no adaptation")
        ( prefixvm( prefix, "ser.radapt-rb-rtol").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "Relative tolerance criterion for RB r-adaptation - default(0.0) = no adaptation")
        ( prefixvm( prefix, "ser.eim-greedy-rtol").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "Relative tolerance criterion for the error indicator constraint in eim Greedy algorithm : choose mu from those which satisfies this criterion - default(0.0) = no error indicator constraint")
        ( prefixvm( prefix, "ser.corrected-rb-rtol").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "Relative tolerance criterion for RB correction (from error estimation) to be used - default(0.0) = no correction")
        ( prefixvm( prefix, "ser.nb-levels").c_str(), Feel::po::value<int>()->default_value( 1 ), "number of SER levels = number of passages into SER algorithm (default 1)")
        ;
    return seroptions;
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
        ( prefixvm( prefix, "crb.db.load").c_str()   , Feel::po::value<int>()->default_value( 2 ), "=0 use db.filename, =1 use last DB created =2 use last DB modified =3 use db.id" )
        ( prefixvm( prefix, "crb.db.update").c_str()   , Feel::po::value<int>()->default_value( 2 ), "=0 use db.filename, =1 use last DB created =2 use last DB modified =3 use db.id" )
        ( prefixvm( prefix, "crb.db.id").c_str()   , Feel::po::value<std::string>()->default_value(""), "DB unique id" )
        ( prefixvm( prefix, "crb.db.filename").c_str()   , Feel::po::value<std::string>()->default_value(""), "DB filename" )
        ( prefixvm( prefix, "crb.output-index").c_str()   , Feel::po::value<int>()->default_value( 0 ), "output index (default is right hand side = 0)" )
        ( prefixvm( prefix, "crb.sampling-size").c_str()   , Feel::po::value<int>()->default_value( 1000 ), "Offline  sampling size " )
        ( prefixvm( prefix, "crb.randomize.use-log").c_str()   , Feel::po::value<bool>()->default_value(true ), "" )
        ( prefixvm( prefix, "crb.sampling-mode").c_str()   , Feel::po::value<std::string>()->default_value( "log-random" ), "Offline  sampling mode, equidistribute, log-equidistribute or log-random " )
        ( prefixvm( prefix, "crb.all-procs-have-same-sampling").c_str() , Feel::po::value<bool>()->default_value( false ), "all procs have the same sampling if true" )
        ( prefixvm( prefix, "crb.error-max").c_str()   , Feel::po::value<double>()->default_value( 1e-6 ),       "Offline  tolerance" )
        ( prefixvm( prefix, "crb.online-tolerance").c_str()   , Feel::po::value<double>()->default_value( 1e-2 ),       "Online  tolerance" )
        ( prefixvm( prefix, "crb.absolute-error").c_str() , Feel::po::value<bool>()->default_value( false ), "Impose to compute absolute error PFEM/CRB instead of relative" )
        ( prefixvm( prefix, "crb.dimension-max").c_str()   , Feel::po::value<int>()->default_value( 1 ),       "Offline max WN size, set to 1 by default to avoid to enrich the existing database if this option doesn't appear in onefeel interface or in the config file." )
        ( prefixvm( prefix, "crb.dimension").c_str()   , Feel::po::value<int>()->default_value( -1 ),       "Online  WN size" )
        ( prefixvm( prefix, "crb.error-type").c_str()   , Feel::po::value<int>()->default_value( ( int )CRB_RESIDUAL_SCM ),       "CRB error type to be computed" )
        ( prefixvm( prefix, "crb.compute-apee-for-each-time-step").c_str(),Feel::po::value<bool>()->default_value( true ),"Compute error estimation for each time step (parabolic problems) is true, else compute only for the last one")
        ( prefixvm( prefix, "crb.factor").c_str()   , Feel::po::value<int>()->default_value( -1 ),  "factor useful to estimate error by empirical method" )
        ( prefixvm( prefix, "crb.Nm").c_str()   , Feel::po::value<int>()->default_value( 1 ),       "Offline  number of modes per mu (for the POD) " )
        ( prefixvm( prefix, "crb.apply-POD-to-WN").c_str()   , Feel::po::value<bool>()->default_value( false ), "apply a POD on approximation functions spaces (primal and dual) if true and if deal with a transient problem " )
        ( prefixvm( prefix, "crb.check.rb").c_str()   , Feel::po::value<int>()->default_value( 0 ),       "check reduced basis" )
        ( prefixvm( prefix, "crb.orthonormality-tol").c_str() , Feel::po::value<double>()->default_value( 1e-13 ),"tolerance of orthonormalisation : i.e. norm of matrix A(i,j)=scalarProduct( Wn[j], Wn[i] )" )
        ( prefixvm( prefix, "crb.orthonormality-max-iter").c_str() , Feel::po::value<int>()->default_value( 10 ),"while the tolerance is not reached, the orthonormalization step is done or until max-iter is reached" )
        ( prefixvm( prefix, "crb.gram-schmidt.selection").c_str() , Feel::po::value<bool>()->default_value( false ),"Use selection in Gram-Schmidt orthonormalizatin to erase useless basis vectors" )
        ( prefixvm( prefix, "crb.gram-schmidt.selection.tol").c_str() , Feel::po::value<double>()->default_value( 1e-8 ),"Selective Gram-Schmidt alogrithm tolerance" )

        ( prefixvm( prefix, "crb.check.residual").c_str()   , Feel::po::value<bool>()->default_value( false ),  "check residual" )
        ( prefixvm( prefix, "crb.reuse-prec").c_str()   , Feel::po::value<bool>()->default_value( 0 ),       "reuse or not the preconditioner" )
        ( prefixvm( prefix, "crb.use-primal-pc").c_str(),Feel::po::value<bool>()->default_value( true ), "use specific preconditioner for the primal problem")
        ( prefixvm( prefix, "crb.orthonormalize-primal").c_str() , Feel::po::value<bool>()->default_value( 1 ), "orthonormalize or not " )
        ( prefixvm( prefix, "crb.orthonormalize-dual").c_str() , Feel::po::value<bool>()->default_value( 1 ), "orthonormalize or not " )
        ( prefixvm( prefix, "crb.solve-dual-problem").c_str() , Feel::po::value<bool>()->default_value( 1 ), "solve or not the dual problem (this bool will be ignored if error-type=CRB_RESIDUAL) " )
        ( prefixvm( prefix, "crb.visualize-basis").c_str() , Feel::po::value<bool>()->default_value( 0 ), "visualize elements of the reduced basis " )
        ( prefixvm( prefix, "crb.save-output-behavior").c_str() , Feel::po::value<bool>()->default_value( 0 ), "save output behavior in time" )
        ( prefixvm( prefix, "crb.seek-mu-in-complement").c_str() , Feel::po::value<bool>()->default_value( 1 ), "during the offline basis construction, see mu in M the complement of Wn" )
        ( prefixvm( prefix, "crb.rebuild-database").c_str() , Feel::po::value<bool>()->default_value( 0 ), "rebuild database (if it already exists)" )
        ( prefixvm( prefix, "crb.restart-from-N").c_str() , Feel::po::value<int>()->default_value( -1 ), "restart the database from specified N (note that when N=0 the complete approximation space is rebuilt and so it is equivalent to option crb.rebuild-database=true). In the case where N > Nmax we do nothing. By default it is set to a negative number in order to not interfer with option rebuild-database" )
        ( prefixvm( prefix, "crb.show-mu-selection").c_str() , Feel::po::value<bool>()->default_value( 0 ), " show mu selection during offline step to build RB space" )
        ( prefixvm( prefix, "crb.show-residual").c_str() , Feel::po::value<bool>()->default_value( 0 ), " show mu residuals values (used for the error estimation)" )
        ( prefixvm( prefix, "crb.print-error-during-rb-construction").c_str() , Feel::po::value<bool>()->default_value( 0 ), " print the max error (absolute) obtained during the offline step" )
        ( prefixvm( prefix, "crb.print-iterations-info").c_str(), Feel::po::value<bool>()->default_value( 0 ), "Write offline Picard iterations summary" )
        ( prefixvm( prefix, "crb.compute-variance").c_str() , Feel::po::value<bool>()->default_value( 0 ), " if true the output is the variance and not l(v)" )
        ( prefixvm( prefix, "crb.save-information-for-variance").c_str(),Feel::po::value<bool>()->default_value( 0 ), "if true will build variance matrix but it takes some times" )

        ( prefixvm( prefix, "crb.use-newton").c_str(),Feel::po::value<bool>()->default_value( false ), "use newton algorithm (need to provide a jacobian and a residual)" )
        ( prefixvm( prefix, "crb.fixedpoint.aitken").c_str(),Feel::po::value<bool>()->default_value( true ), "use Aitken relaxtion algorithm in nonlinear fixpoint solver" )

        ( prefixvm( prefix, "crb.fixedpoint.maxit").c_str(),Feel::po::value<int>()->default_value( 20 ), "nb iteration max for the fixed point (online part)" )
        ( prefixvm( prefix, "crb.fixedpoint.increment-tol").c_str(),Feel::po::value<double>()->default_value( 1e-10 ), "tolerance on solution for fixed point (online part)" )
        ( prefixvm( prefix, "crb.fixedpoint.output-tol").c_str(),Feel::po::value<double>()->default_value( 1e-10 ), "tolerance on output for fixed point (online part)" )
        ( prefixvm( prefix, "crb.fixedpoint.verbose").c_str(),Feel::po::value<bool>()->default_value( false ), "fixed point verbose if true" )
        ( prefixvm( prefix, "crb.fixedpoint.critical-value").c_str(),Feel::po::value<double>()->default_value(1000 ), "will crash if increment error at the end of fixed point is greater than this critical value" )

        ( prefixvm( prefix, "crb.use-continuity").c_str(),Feel::po::value<bool>()->default_value(true), "when apply Newton method, will use continuity method (like for natural convection problem for example)" )
        ( prefixvm( prefix, "crb.compute-error-on-reduced-residual-jacobian").c_str(),Feel::po::value<bool>()->default_value( false ), "only for crb_trilinear")
        ( prefixvm( prefix, "crb.enable-convection-terms").c_str(),Feel::po::value<bool>()->default_value( true ), "only for crb_trilinear")

        ( prefixvm( prefix, "crb.is-model-executed-in-steady-mode").c_str(),Feel::po::value<bool>()->default_value( false ), "true if model is executed in steady mode, else turn it to false")
        ( prefixvm( prefix, "crb.use-ginac-for-beta-expressions").c_str(),Feel::po::value<bool>()->default_value( false ), "use ginac to compute expression of beta coefficients if true")
        ( prefixvm( prefix, "crb.use-linear-model").c_str(),Feel::po::value<bool>()->default_value( false ), "do not iterate in fixed point if true")

        ( prefixvm( prefix, "crb.use-predefined-WNmu").c_str(),Feel::po::value<bool>()->default_value( false ), "read parameters to take for the offline step from a file named SamplingWNmu if true")
        ( prefixvm( prefix, "crb.use-predefined-test-sampling").c_str(),Feel::po::value<bool>()->default_value( false ), "read parameters from file named SamplingForTest if true to run the test")
        ( prefixvm( prefix, "crb.use-logEquidistributed-WNmu").c_str(),Feel::po::value<int>()->default_value( 0 ), "parameters are log-equidistributed for the offline step (the value indicates the number of parameters)")
        ( prefixvm( prefix, "crb.use-random-WNmu").c_str(),Feel::po::value<int>()->default_value( 0 ), "parameters are log-equidistributed for the offline step (the value indicates the number of parameters)")
        ( prefixvm( prefix, "crb.use-equidistributed-WNmu").c_str(),Feel::po::value<int>()->default_value( 0 ), "parameters are equidistributed for the offline step (the value indicates the number of parameters)")

        ( prefixvm( prefix, "crb.compute-stat").c_str(),Feel::po::value<bool>()->default_value( true ), "compute statistics on the run if true")
        ( prefixvm( prefix, "crb.cvg-study").c_str(),Feel::po::value<bool>()->default_value( false ), "convergence study if true")
        ( prefixvm( prefix, "crb.computational-time-neval").c_str(),Feel::po::value<int>()->default_value( 0 )," number of evaluation to perform to have the computational time of crb online step" )
        ( prefixvm( prefix, "crb.reload-last-sampling").c_str(),Feel::po::value<bool>()->default_value( false ), "")
        ( prefixvm( prefix, "crb.run-on-WNmu").c_str(),Feel::po::value<bool>()->default_value( false ), "use mu taken for build the reduced basis, so for steady problems we are very accurate")
        ( prefixvm( prefix, "crb.run-on-scm-parameters").c_str(),Feel::po::value<bool>()->default_value( false ), "use mu taken during the SCM offline step ( for a(.,.;mu) ), so the coercivity constant is exact")
        ( prefixvm( prefix, "crb.script-mode").c_str(),Feel::po::value<bool>()->default_value( false ), "disable error computation (need FEM computation) if true")
        ( "crb.db.format", Feel::po::value<std::string>()->default_value("boost"), "format in which the crb database is saved, either boost of hdf5")
        ( prefixvm( prefix, "crb.results-repo-name").c_str(), Feel::po::value<std::string>()->default_value("default_repo"), "name for results repository, and also use for database storage")
        ( prefixvm( prefix, "crb.compute-fem-during-online").c_str(),Feel::po::value<bool>()->default_value( true ), "compute fem during online step, necessary to compute the error between fem and crb")

        ( prefixvm( prefix, "crb.compute-matrix-information").c_str(),Feel::po::value<bool>()->default_value( false ), "compute matrix information (i.e. conditioning, determinant) of reduced matrix if true")
        ( prefixvm( prefix, "crb.print-rb-matrix").c_str(),Feel::po::value<bool>()->default_value( false ), "write rb matrix (octave format) in a file if true")

        ( prefixvm( prefix, "crb.use-symmetric-matrix").c_str(),Feel::po::value<bool>()->default_value( true ), "don't transpose to have the matrix associated to the dual problem if true")
        ( prefixvm( prefix, "crb.stock-matrices").c_str(),Feel::po::value<bool>()->default_value( true ), "assemble and stock all matrices/vectors if true, but it can takes a lot of memory")
        ( prefixvm( prefix, "crb.system-memory-evolution").c_str(),Feel::po::value<bool>()->default_value( false ), "generate a file to plot memory evolution during offline step only on the master processor (file written : MemoryEvolution)")
        ( prefixvm( prefix, "crb.system-memory-evolution-on-all-procs").c_str(),Feel::po::value<bool>()->default_value( false ), "same than system-memory-evolution but on all processors")

        ( prefixvm( prefix, "crb.use-accurate-apee").c_str(),Feel::po::value<bool>()->default_value( false ), "use a posteriori error estimators from F.Casenave's paper if true, classic one else")
        ( prefixvm( prefix, "crb.optimize-offline-residual").c_str(),Feel::po::value<bool>()->default_value( false ), "use optimize way for offline residual computation if true (temporary option)")
        ( prefixvm( prefix, "crb.offline-residual-version").c_str(),Feel::po::value<int>()->default_value( 0 ), "offline residual version 0 : old no storage, 1 : new faster with storage")

        ( prefixvm( prefix, "crb.user-parameters").c_str(),Feel::po::value<std::string>()->default_value( "" ), "values of parameters (used for one feel)")
        ( prefixvm( prefix, "crb.vary-only-parameter-components").c_str(),Feel::po::value<std::string>()->default_value( "" ), "specify which parameter component vary (max : 2 components + time 't') and how many values we take in each direction. For example 0 10 1 20 means that component 0 will take 10 values and component 1 will take 20 values. We can write also t 1 10 and in this case the time will vary and also component 1")
        ( prefixvm( prefix, "crb.load-elements-database").c_str(),Feel::po::value<bool>()->default_value( false ), "load database of elements if true, need to be true for visualization, need to be false to run CRB approximation on a different number of processors than this was used to build the reduced basis ")

        ( prefixvm( prefix, "crb.solve-fem-monolithic").c_str(),Feel::po::value<bool>()->default_value( false ), "solve FEM problem without using EIM and without affine decomposition ")
        ( prefixvm( prefix, "crb.export-name-max-size").c_str(),Feel::po::value<int>()->default_value( 30 ), "maximum size for variable names in export (truncature)")

        ("crb.minimization-func", Feel::po::value<std::string>()->default_value(""), "giving a functional f(output) - give the output which minimizes f(output)" )
        ("crb.minimization-param-name", Feel::po::value<std::string>()->default_value( "output" ), "name of the parameter to be replaced by the output in expression given by crb.minimization-func")

        ( prefixvm( prefix, "crb.use-fast-eim").c_str(),Feel::po::value<bool>()->default_value( true ), "use fast eim algo (with rbspace context)")



        ;

    crboptions
        .add( crbSCMOptions() ).add(deimOptions(prefix));

    return crboptions;
}

Feel::po::options_description
crbBlockOptions( int const& n_block )
{
    Feel::po::options_description crboptions( "CRB Block Options" );
    for ( int i=0; i<n_block; i++ )
    {
        crboptions.add_options()
            ( fmt::format("crb.block.orthonormalize{}", i ).c_str(), Feel::po::value<bool>()->default_value( true ), fmt::format("orthonormalize reduce basis for rbspace {}", i ).c_str() )
            ;
        crboptions.add( backend_options("backend-Xh"+std::to_string(i)) );
    }

    return crboptions;
}


Feel::po::options_description
crbSaddlePointOptions( std::string const& prefix, int const& n_block )
{
    Feel::po::options_description crboptions( "CRB Options" );
    crboptions.add_options()
        ( prefixvm( prefix, "crb.saddlepoint.add-supremizer").c_str(),Feel::po::value<bool>()->default_value( false ), "add the supremizer function to the first reduced basis")
        ( prefixvm( prefix, "crb.saddlepoint.orthonormalize0").c_str(),Feel::po::value<bool>()->default_value( true ), "orthonormalize reduce basis for rbspace #0")
        ( prefixvm( prefix, "crb.saddlepoint.orthonormalize1").c_str(),Feel::po::value<bool>()->default_value( true ), "orthonormalize reduce basis for rbspace #1")
        ( prefixvm( prefix, "crb.saddlepoint.version").c_str(),Feel::po::value<int>()->default_value( 1 ), "test residual evaluation")
        ;

    return crboptions.add( crbBlockOptions(n_block) );
}

Feel::po::options_description
crbAeroOptions( std::string const& prefix )
{
    Feel::po::options_description crboptions( "CRB Aero Options" );
    crboptions.add_options()
        ( prefixvm( prefix, "crb.aero.add-supremizer").c_str(),Feel::po::value<bool>()->default_value( false ), "add the supremizer function to the first reduced basis")
        ( prefixvm( prefix, "crb.aero.fix-mean-pressure").c_str(),Feel::po::value<bool>()->default_value( false ), "")
        ( prefixvm( prefix, "crb.aero.use-psit").c_str(),Feel::po::value<bool>()->default_value( false ), "")
        ( prefixvm( prefix, "crb.aero.psit.delta0").c_str(),Feel::po::value<double>()->default_value( 1. ), "")
        ( prefixvm( prefix, "crb.aero.linear-solve").c_str(),Feel::po::value<bool>()->default_value(false ), "")
        ( prefixvm( prefix, "crb.aero.use-newton").c_str(),Feel::po::value<bool>()->default_value(true ), "")
        ( prefixvm( prefix, "crb.aero.init-online").c_str(),Feel::po::value<bool>()->default_value(true ), "")
        ( prefixvm( prefix, "crb.aero.snes.rtol").c_str(),Feel::po::value<double>()->default_value( 1e-8 ), "")
        ( prefixvm( prefix, "crb.aero.online-continuation").c_str(),Feel::po::value<int>()->default_value( 10 ), "")
        ( prefixvm( prefix, "crb.aero.log-continuation").c_str(),Feel::po::value<bool>()->default_value( true ), "")
        ( prefixvm( prefix, "crb.aero.store-rb-sol").c_str(),Feel::po::value<bool>()->default_value( false ), "")
        ( prefixvm( prefix, "crb.aero.assemble-version").c_str(),Feel::po::value<int>()->default_value( 1 ), "")
        ;
    crboptions.add( crbSaddlePointOptions(prefix,3) );

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
}