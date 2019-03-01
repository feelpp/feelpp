/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2012-01-19

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
 \file modelalgebraicfactory.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>

namespace Feel
{
namespace FeelModels
{

    ModelAlgebraicFactory::ModelAlgebraicFactory(model_ptrtype const& model, backend_ptrtype const& __backend)
        :
        M_model(model),
        M_backend( __backend ),
        M_hasBuildLinearJacobian(false),
        M_hasBuildResidualCst(false),
        M_hasBuildLinearSystemCst(false),
        M_usePseudoTransientContinuation( boption(_prefix=model->prefix(),_name="pseudo-transient-continuation") ),
        M_pseudoTransientContinuationEvolutionMethod( soption(_prefix=model->prefix(),_name="pseudo-transient-continuation.evolution") ),
        M_pseudoTransientContinuationDelta0( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.delta0") ),
        M_pseudoTransientContinuationDeltaMax( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.delta-max") ),
        M_pseudoTransientContinuationSerVariant( soption(_prefix=model->prefix(),_name="pseudo-transient-continuation.ser-variant") ),
        M_pseudoTransientContinuationExpurThresholdHigh( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.expur.threshold-high") ),
        M_pseudoTransientContinuationExpurThresholdLow( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.expur.threshold-low") ),
        M_pseudoTransientContinuationExpurBetaHigh( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.expur.beta-high") ),
        M_pseudoTransientContinuationExpurBetaLow( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.expur.beta-low") )
    {
        model->timerTool("Constructor").start();
        auto graph = model->buildMatrixGraph();
        model->timerTool("Constructor").elapsed("graph");
        if ( graph )
        {
            model->timerTool("Constructor").restart();
            this->buildMatrixVector( graph,graph->mapRow().indexSplit() );
            model->timerTool("Constructor").elapsed("matrixVector");

            model->timerTool("Constructor").restart();
            this->buildOthers();
            model->timerTool("Constructor").elapsed("algebraicOthers");
        }
        model->timerTool("Constructor").stop();
    }

    //---------------------------------------------------------------------------------------------------------------//

    ModelAlgebraicFactory::ModelAlgebraicFactory(model_ptrtype const& model,
                           backend_ptrtype const& __backend,
                           graph_ptrtype const& graph,
                           indexsplit_ptrtype const& indexSplit )
        :
        M_model(model),
        M_backend(__backend ),
        M_hasBuildLinearJacobian(false),
        M_hasBuildResidualCst(false),
        M_hasBuildLinearSystemCst(false),
        M_usePseudoTransientContinuation( boption(_prefix=model->prefix(),_name="pseudo-transient-continuation") ),
        M_pseudoTransientContinuationEvolutionMethod( soption(_prefix=model->prefix(),_name="pseudo-transient-continuation.evolution") ),
        M_pseudoTransientContinuationDelta0( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.delta0") ),
        M_pseudoTransientContinuationDeltaMax( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.delta-max") ),
        M_pseudoTransientContinuationSerVariant( soption(_prefix=model->prefix(),_name="pseudo-transient-continuation.ser-variant") ),
        M_pseudoTransientContinuationExpurThresholdHigh( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.expur.threshold-high") ),
        M_pseudoTransientContinuationExpurThresholdLow( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.expur.threshold-low") ),
        M_pseudoTransientContinuationExpurBetaHigh( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.expur.beta-high") ),
        M_pseudoTransientContinuationExpurBetaLow( doption(_prefix=model->prefix(),_name="pseudo-transient-continuation.expur.beta-low") )
    {
        if (this->model()->verbose()) Feel::FeelModels::Log(this->model()->prefix()+".MethodNum","constructor1", "start",
                                                            this->model()->worldComm(),this->model()->verboseAllProc());

        this->init(graph,indexSplit);

        if (this->model()->verbose()) Feel::FeelModels::Log(this->model()->prefix()+".MethodNum","constructor", "finish",
                                                            this->model()->worldComm(),this->model()->verboseAllProc());
    }

    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::init(graph_ptrtype const& graph,
                     indexsplit_ptrtype const& indexSplit)
    {
        this->model()->timerTool("Constructor").start();
        this->buildMatrixVector(graph,indexSplit);
        this->model()->timerTool("Constructor").elapsed("matrixVector");

        this->model()->timerTool("Constructor").restart();
        this->buildOthers();
        this->model()->timerTool("Constructor").stop("algebraicOthers");
    }

    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::buildMatrixVector(graph_ptrtype const& graph,
                                             indexsplit_ptrtype const& indexSplit)
    {
        //vectors
        M_R = M_backend->newVector(graph->mapRowPtr()); //  size1,size1 ;
        if (this->model()->useCstVector())
            M_CstR = M_backend->newVector(graph->mapRowPtr());//size1,size1 );
        else M_CstR = M_R;
        if (this->model()->verbose()) Feel::FeelModels::Log(this->model()->prefix()+".MethodNum","buildMatrixVector","build all vectors",
                                                           this->model()->worldComm(),this->model()->verboseAllProc());

        size_type size1=0,size2=0;
        //matrix with standart pattern
        M_J = M_backend->newMatrix(size1,size2,size1,size2,graph,indexSplit);
        if (this->model()->useCstMatrix())
            M_CstJ = M_backend->newMatrix(size1,size2,size1,size2,graph,indexSplit);
        else M_CstJ = M_J;
        if (this->model()->buildMatrixPrecond())
            M_Prec = M_backend->newMatrix( size1,size2,size1,size2,graph,indexSplit );
        else M_Prec = M_J;
        //matrix with extended pattern
        if (this->model()->hasExtendedPattern() && this->model()->buildMatrixPrecond())
            M_Extended = M_backend->newMatrix( size1,size2,size1,size2,graph,indexSplit );
        else
            M_Extended = M_J;

        if (this->model()->verbose()) Feel::FeelModels::Log(this->model()->prefix()+".MethodNum","buildMatrixVector","build all matrix",
                                                           this->model()->worldComm(),this->model()->verboseAllProc());

        // M_J = M_backend->newMatrix( _trial=M_Xh, _test=M_Xh, _pattern_block=myblockpattern );
        // M_Extended = M_backend->newZeroMatrix( size1,size2,size1,size2 );

        if (this->model()->printGraph())
            graph->printPython( this->model()->printGraphFileName() );
    }

    void
    ModelAlgebraicFactory::buildOthers()
    {
        M_PrecondManage = preconditioner(_pc=(PreconditionerType) M_backend->pcEnumType() /*LU_PRECOND*/,
                                         _matrix=M_Prec,
                                         _backend= M_backend,
                                         _pcfactormatsolverpackage=(MatSolverPackageType) M_backend->matSolverPackageEnumType(),
                                         _worldcomm=M_backend->comm(),
                                         _prefix=M_backend->prefix() );

        // necessary to use the nonlinear solver
        M_backend->nlSolver()->jacobian = std::bind( &self_type::updateJacobian,
                                                     std::ref( *this ), std::placeholders::_1, std::placeholders::_2/*, M_R*/ );
        M_backend->nlSolver()->residual = std::bind( &self_type::updateResidual,
                                                     std::ref( *this ), std::placeholders::_1, std::placeholders::_2 );

        if ( M_usePseudoTransientContinuation && M_pseudoTransientContinuationSerVariant == "solution" )
        {
            M_pseudoTransientContinuationPreviousSolution = M_backend->newVector( M_J->mapColPtr() );
        }

    }
    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::rebuildMatrixVector(graph_ptrtype const& graph,
                                               indexsplit_ptrtype const& indexSplit)
    {
        M_hasBuildLinearJacobian=false;
        M_hasBuildResidualCst=false;
        M_hasBuildLinearSystemCst=false;
        this->buildMatrixVector(graph,indexSplit);

        M_PrecondManage->setMatrix(M_Prec);
    }

    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::reset(backend_ptrtype __backend,
                                 graph_ptrtype const& graph,
                                 indexsplit_ptrtype const& indexSplit)
    {

        M_backend=__backend;
        M_hasBuildLinearJacobian=false;
        M_hasBuildResidualCst=false;
        M_hasBuildLinearSystemCst=false;

        this->init(graph,indexSplit);
    }

    void
    ModelAlgebraicFactory::attachNullSpace( NullSpace<value_type> const& nullSpace )
    {
        CHECK( this->backend() ) << "backend not init\n";
        std::shared_ptr<NullSpace<value_type> > mynullspace( new NullSpace<value_type>(this->backend(),nullSpace) );
        this->backend()->attachNullSpace( mynullspace );
    }
    void
    ModelAlgebraicFactory::attachNearNullSpace( NullSpace<value_type> const& nearNullSpace )
    {
        CHECK( this->backend() ) << "backend not init\n";
        std::shared_ptr<NullSpace<value_type> > myNearNullSpace( new NullSpace<value_type>(this->backend(),nearNullSpace) );
        this->backend()->attachNearNullSpace( myNearNullSpace );
    }
    void
    ModelAlgebraicFactory::attachNearNullSpace( int k, NullSpace<value_type> const& nearNullSpace )
    {
        CHECK( this->backend() ) << "backend not init\n";
        CHECK( M_PrecondManage ) << "preconditioner not init\n";
        std::shared_ptr<NullSpace<value_type> > myNearNullSpace( new NullSpace<value_type>(this->backend(),nearNullSpace) );
        M_PrecondManage->attachNearNullSpace( k, myNearNullSpace );
    }

    void
    ModelAlgebraicFactory::addFunctionLinearAssembly( function_assembly_linear_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionLinearAssembly.size()).str() : key;
        M_addFunctionLinearAssembly[ keyUsed ] = func;
    }
    void
    ModelAlgebraicFactory::addFunctionLinearPostAssembly( function_assembly_linear_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionLinearPostAssembly.size()).str() : key;
        M_addFunctionLinearPostAssembly[ keyUsed ] = func;
    }
    void
    ModelAlgebraicFactory::addFunctionNewtonInitialGuess( function_newton_initial_guess_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionNewtonInitialGuess.size()).str() : key;
        M_addFunctionNewtonInitialGuess[ keyUsed ] = func;
    }
    void
    ModelAlgebraicFactory::addFunctionJacobianAssembly( function_assembly_jacobian_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionJacobianAssembly.size()).str() : key;
        M_addFunctionJacobianAssembly[ keyUsed ] = func;
    }
    void
    ModelAlgebraicFactory::addFunctionResidualAssembly( function_assembly_residual_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionResidualAssembly.size()).str() : key;
        M_addFunctionResidualAssembly[ keyUsed ] = func;
    }
    void
    ModelAlgebraicFactory::addFunctionJacobianPostAssembly( function_assembly_jacobian_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionJacobianPostAssembly.size()).str() : key;
        M_addFunctionJacobianPostAssembly[ keyUsed ] = func;
    }
    void
    ModelAlgebraicFactory::addFunctionResidualPostAssembly( function_assembly_residual_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionResidualPostAssembly.size()).str() : key;
        M_addFunctionResidualPostAssembly[ keyUsed ] = func;
    }


    void
    ModelAlgebraicFactory::addVectorLinearRhsAssembly( vector_ptrtype const& vec, double scaling, std::string const& key, bool cstPart )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addVectorLinearRhsAssembly.size()).str() : key;
        M_addVectorLinearRhsAssembly[ keyUsed ] = std::make_tuple( vec, scaling, cstPart );
    }
    void
    ModelAlgebraicFactory::addVectorResidualAssembly( vector_ptrtype const& vec, double scaling, std::string const& key, bool cstPart )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addVectorResidualAssembly.size()).str() : key;
        M_addVectorResidualAssembly[ keyUsed ] = std::make_tuple( vec, scaling, cstPart );
    }

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::updateInformationObject( pt::ptree & p )
    {
        pt::ptree subPt;
        subPt.put( "prefix",this->backend()->prefix() );
        subPt.put( "type", "BACKEND_PETSC" );
        p.put_child( "Backend", subPt );
        // KSP
        subPt.clear();
        subPt.put( "type", this->backend()->kspType() );
        subPt.put( "rtol", this->backend()->rTolerance() );
        subPt.put( "dtol", this->backend()->dTolerance() );
        subPt.put( "atol", this->backend()->aTolerance() );
        subPt.put( "maxit", this->backend()->maxIterationsKSP() );
        subPt.put( "reuse-prec", this->backend()->reusePrec() );
        if ( this->backend()->reusePrec() )
            subPt.put( "maxit reuse-prec", this->backend()->maxIterationsKSPReuse() );
        p.put_child( "KSP", subPt );

        //SNES
        subPt.clear();
        subPt.put( "rtol",this->backend()->rToleranceSNES() );
        subPt.put( "stol",this->backend()->sToleranceSNES() );
        subPt.put( "atol",this->backend()->aToleranceSNES() );
        subPt.put( "maxit", this->backend()->maxIterationsSNES() );
        subPt.put( "reuse-jac", this->backend()->reuseJac() );
        if ( this->backend()->reuseJac() )
        {
            subPt.put( "reuse-jac rebuild at first Newton step", this->backend()->reuseJacRebuildAtFirstNewtonStep() );
            subPt.put( "maxit with reuse", this->backend()->maxIterationsSNESReuse() );
        }
        p.put_child( "SNES", subPt );

        // KSP in SNES
        subPt.clear();
        subPt.put( "rtol",this->backend()->rtoleranceKSPinSNES() );
        subPt.put( "maxit",this->backend()->maxIterationsKSPinSNES() );
        subPt.put( "reuse-prec", this->backend()->reusePrec() );
        if ( this->backend()->reusePrec() )
        {
            subPt.put( "reuse-prec maxit", this->backend()->maxIterationsKSPinSNESReuse() );
            subPt.put( "reuse-prec rebuild at first Newton step", this->backend()->reusePrecRebuildAtFirstNewtonStep() );
        }
        p.put_child( "KSP in SNES", subPt );

        // PC
        subPt.clear();
        subPt.put( "type", this->backend()->pcType() );
        if ( this->backend()->pcType() == "lu" )
            subPt.put( "mat-solver-package", this->backend()->pcFactorMatSolverPackageType() );
        p.put_child( "PC", subPt );
    }

    std::shared_ptr<std::ostringstream>
    ModelAlgebraicFactory::getInfo() const
    {
        std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );

        *_ostr << "\n||==============================================||"
               << "\n||---------------Info : SolverNum---------------||"
               << "\n||==============================================||";

        *_ostr << "\n   Backend "
               << "\n     -- prefix : " << this->backend()->prefix()
               << "\n     -- type : " << "BACKEND_PETSC"
               << "\n   KSP "
               << "\n     -- type : " << this->backend()->kspType()
               << "\n     -- tolerance : res=" << this->backend()->rTolerance() << ", diff=" << this->backend()->dTolerance() << ", abs=" << this->backend()->aTolerance()
               << "\n     -- maxit : " << this->backend()->maxIterationsKSP()
               << "\n     -- maxit : " << this->backend()->maxIterationsKSPReuse()
               << "\n     -- reuse prec : " << (this->backend()->reusePrec()?"true":"false")
               << "\n   SNES "
               << "\n     -- reuse jac : " << (this->backend()->reuseJac()?"true":"false")
               << "\n     -- reuse jac rebuild at first Newton step : " << (this->backend()->reuseJacRebuildAtFirstNewtonStep()?"true":"false")
               << "\n     -- tolerance : res=" << this->backend()->rToleranceSNES() << ", step=" << this->backend()->sToleranceSNES() << ", abs=" << this->backend()->aToleranceSNES()
               << "\n     -- maxit : " << this->backend()->maxIterationsSNES()
               << "\n     -- maxit with reuse : " << this->backend()->maxIterationsSNESReuse()
               << "\n   KSP in SNES "
               << "\n     -- tolerance : res=" << this->backend()->rtoleranceKSPinSNES()
               << "\n     -- maxit : " << this->backend()->maxIterationsKSPinSNES()
               << "\n     -- maxit with reuse : " << this->backend()->maxIterationsKSPinSNESReuse()
               << "\n     -- reuse prec : " << (this->backend()->reusePrec()?"true":"false")
               << "\n     -- reuse prec rebuild at first Newton step : " << (this->backend()->reusePrecRebuildAtFirstNewtonStep()?"true":"false")
               << "\n   PC "
               << "\n     -- type : " << this->backend()->pcType();
        if ( this->backend()->pcType() == "lu" )
            *_ostr << "\n     -- mat solver package : " << this->backend()->pcFactorMatSolverPackageType();

        return _ostr;
    }

    void
    ModelAlgebraicFactory::printInfo() const
    {
        if ( this->model()->verboseAllProc() ) std::cout << this->getInfo()->str();
        else if (this->model()->worldComm().isMasterRank() )
            std::cout << this->getInfo()->str();
    }


    void
    ModelAlgebraicFactory::solve( std::string const& type,vector_ptrtype& sol )
    {
        if ( type == "Linear" || type == "LinearSystem" )
        {
            this->solveLinear( sol );
        }
        else if ( type == "Picard" || type == "FixPoint" )
        {
            this->solvePicard( sol );
        }
        else if ( type == "Newton" )
        {
            this->solveNewton( sol );
        }
        else
            CHECK( false ) << "invalid solver type " << type;
    }

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::solveLinear( vector_ptrtype& U )
    {
        auto model = this->model();
        if (this->model()->verbose()) Feel::FeelModels::Log(this->model()->prefix()+".ModelAlgebraicFactory","linearSolver", "start",
                                                            this->model()->worldComm(),this->model()->verboseAllProc());

        this->model()->timerTool("Solve").start();

        // assembling cst part
        if (!M_hasBuildLinearSystemCst)
        {
            M_CstJ->zero();
            M_CstR->zero();
            //this->model()->updateLinearPDE(U,M_CstJ,M_CstR,true,M_Extended,false);
            ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true,M_Extended,false);
            this->model()->updateLinearPDE( dataLinearCst );
            for ( auto const& func : M_addFunctionLinearAssembly )
                func.second( dataLinearCst );
            M_CstR->close();
            for ( auto const& av : M_addVectorLinearRhsAssembly )
                if ( std::get<2>( av.second ) )
                    M_CstR->add( std::get<1>( av.second ), std::get<0>( av.second ) );
            M_hasBuildLinearSystemCst = true;
        }
        else if ( this->model()->rebuildCstPartInLinearSystem() || this->model()->needToRebuildCstPart() ||
                  !this->model()->useCstMatrix() || !this->model()->useCstVector() )
        {
            M_CstJ->zero();
            M_CstR->zero();
            //this->model()->updateLinearPDE(U,M_CstJ,M_CstR,true,M_Extended,false);
            ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true,M_Extended,false);
            this->model()->updateLinearPDE( dataLinearCst );
            for ( auto const& func : M_addFunctionLinearAssembly )
                func.second( dataLinearCst );
            M_CstR->close();
            for ( auto const& av : M_addVectorLinearRhsAssembly )
                if ( std::get<2>( av.second ) )
                    M_CstR->add( std::get<1>( av.second ), std::get<0>( av.second ) );
        }
        this->model()->setNeedToRebuildCstPart(false);

        // copying cst matrix/vector
        if (this->model()->useCstMatrix())
        {
            M_J->zero();
            M_J->addMatrix(1.0, M_CstJ );
        }
        if (this->model()->useCstVector())
        {
            M_R->zero();
            M_R->add(1.0, M_CstR );
        }

        // set M_Extended to zero if really used and build
        if (this->model()->hasExtendedPattern() && this->model()->buildMatrixPrecond() )
            M_Extended->zero();

        // assembling non cst part
        ModelAlgebraic::DataUpdateLinear dataLinearNonCst(U,M_J,M_R,false,M_Extended,true);
        // apply before addFunctionLinearAssembly because due to Strong dirichlet
        for ( auto const& func : M_addFunctionLinearAssembly )
            func.second( dataLinearNonCst );
        this->model()->updateLinearPDE( dataLinearNonCst );

        M_R->close();
        // add maybe vector to rhs
        for ( auto const& av : M_addVectorLinearRhsAssembly )
            if ( !std::get<2>( av.second ) )
                M_R->add( std::get<1>( av.second ), std::get<0>( av.second ) );

        // dof elimination
        this->model()->updateLinearPDEDofElimination( dataLinearNonCst );

        // post-assembly
        for ( auto const& func : M_addFunctionLinearPostAssembly )
            func.second( dataLinearNonCst );

        // assembling matrix used for preconditioner
        this->model()->updatePreconditioner(U,M_J,M_Extended,M_Prec);

        // add extended term (important : after updatePreconditioner !)
        // and only do if shared_ptr M_J != M_Extended
        if (this->model()->hasExtendedPattern() && this->model()->buildMatrixPrecond() )
            M_J->addMatrix(1.0, M_Extended );

        double mpiTimerAssembly = this->model()->timerTool("Solve").elapsed("algebraic-assembly");

        if (this->model()->verboseSolverTimer())
          Feel::FeelModels::Log(this->model()->prefix()+".ModelAlgebraicFactory","linearSolver",
                         (boost::format("finish assembling in %1% s") % mpiTimerAssembly ).str(),
                         this->model()->worldComm(),this->model()->verboseSolverTimerAllProc());

        this->model()->timerTool("Solve").restart();

        // set preconditioner
        //M_PrecondManage->setMatrix(M_Prec);

        this->model()->updateInHousePreconditioner( M_J, U );

        pre_solve_type pre_solve = std::bind(&model_type::preSolveLinear, model, std::placeholders::_1, std::placeholders::_2);
        post_solve_type post_solve = std::bind(&model_type::postSolveLinear, model, std::placeholders::_1, std::placeholders::_2);

        // solve linear system
        auto const solveStat = M_backend->solve( _matrix=M_J,
                                                 _solution=U,
                                                 _rhs=M_R,
                                                 _prec=M_PrecondManage,
                                                 _pre=pre_solve,
                                                 _post=post_solve );

        if ( this->model()->errorIfSolverNotConverged() )
            CHECK( solveStat.isConverged() ) << "the linear solver has not converged in "
                                             << solveStat.nIterations() << " iterations with norm "
                                             << solveStat.residual() << "\n";

        double tElapsed = this->model()->timerTool("Solve").stop("algebraic-solve");
        this->model()->timerTool("Solve").setAdditionalParameter("ksp-niter",int(solveStat.nIterations()) );
#if 0
        if ( option(_name="clear-preconditioner-after-use",_prefix=this->model()->prefix()).as<bool>() )
        {
            Feel::FeelModels::Log(this->model()->prefix()+".ModelAlgebraicFactory","linearSolver","clear-preconditioner-after-use",
                           this->model()->worldComm(),this->model()->verboseSolverTimerAllProc());
            M_PrecondManage->clear();
            MatrixPetsc<double> * pmatrix = dynamic_cast<MatrixPetsc<double>*>( M_J.get() );
            pmatrix->mapSplitPC().clear();
        }
#endif
        if (this->model()->verboseSolverTimer())
          Feel::FeelModels::Log(this->model()->prefix()+".ModelAlgebraicFactory","linearSolver",
                         (boost::format("finish solve in %1% s")%tElapsed ).str(),
                         this->model()->worldComm(),this->model()->verboseSolverTimerAllProc());
    }



    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//


    void
    ModelAlgebraicFactory::updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J/*, vector_ptrtype& R*/ )
    {
        auto model = this->model();
        model->timerTool("Solve").start();

        auto R = this->backend()->newVector(X->mapPtr());

        J->zero();

        if ( model->useCstMatrix())
        {
            J->addMatrix(1.0, M_CstJ );
        }
        else
        {
            ModelAlgebraic::DataUpdateJacobian dataJacobianCst(X, J, R, true, M_Extended,false);
            model->updateJacobian( dataJacobianCst );
            for ( auto const& func : M_addFunctionJacobianAssembly )
                func.second( dataJacobianCst );
        }

        ModelAlgebraic::DataUpdateJacobian dataJacobianNonCst(X,J,R,false, M_Extended,false);
        if ( M_usePseudoTransientContinuation )
        {
            dataJacobianNonCst.addInfo( "use-pseudo-transient-continuation" );
            CHECK( !M_pseudoTransientContinuationDeltaAndResidual.empty() ) << "must have at least one value";
            dataJacobianNonCst.addDoubleInfo( "pseudo-transient-continuation.delta", M_pseudoTransientContinuationDeltaAndResidual.back().first );
        }
        // apply before addFunctionJacobianAssembly because due to Strong dirichlet
        for ( auto const& func : M_addFunctionJacobianAssembly )
            func.second( dataJacobianNonCst );
        model->updateJacobian( dataJacobianNonCst );

        // dof elimination
        model->updateJacobianDofElimination( dataJacobianNonCst );

        for ( auto const& func : M_addFunctionJacobianPostAssembly )
            func.second( dataJacobianNonCst );

        model->updateInHousePreconditioner( J, X );

        double tElapsed = model->timerTool("Solve").stop();
        model->timerTool("Solve").addDataValue("algebraic-jacobian",tElapsed);
    }
    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::updateResidual(const vector_ptrtype& X, vector_ptrtype& R)
    {
        auto model = this->model();
        model->timerTool("Solve").start();

        R->zero();

        if ( model->useCstVector())
        {
            R->add(1.0, M_CstR );
        }
        else
        {
            ModelAlgebraic::DataUpdateResidual dataResidualCst( X, R, true, true );
            model->updateResidual( dataResidualCst );
            for ( auto const& func : M_addFunctionResidualAssembly )
                func.second( dataResidualCst );
            R->close();
            for ( auto const& av : M_addVectorResidualAssembly )
                if ( std::get<2>( av.second ) )
                    R->add( std::get<1>( av.second ), std::get<0>( av.second ) );
        }

        bool doOptimization = this->model()->useLinearJacobianInResidual() && this->model()->useCstMatrix();
        // add linear contribution from jacobian terms
        if (doOptimization)
            R->addVector(*X, *M_CstJ );

        ModelAlgebraic::DataUpdateResidual dataResidualNonCst( X, R, false, doOptimization );
        // apply before addFunctionResidualAssembly because due to Strong dirichlet
        for ( auto const& func : M_addFunctionResidualAssembly )
            func.second( dataResidualNonCst );
        model->updateResidual( dataResidualNonCst );

        R->close();

        for ( auto const& av : M_addVectorResidualAssembly )
            if ( !std::get<2>( av.second ) )
                R->add( std::get<1>( av.second ), std::get<0>( av.second ) );

        // dof elimination
        model->updateResidualDofElimination( dataResidualNonCst );

        for ( auto const& func : M_addFunctionResidualPostAssembly )
            func.second( dataResidualNonCst );

        double tElapsed = model->timerTool("Solve").stop();
        model->timerTool("Solve").addDataValue("algebraic-residual",tElapsed);
    }

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::solveNewton( vector_ptrtype& U )
    {
        auto model = this->model();
        if ( model->verbose() )
            Feel::FeelModels::Log( model->prefix()+".ModelAlgebraicFactory","NonLinearSolverNewton", "start",
                                   model->worldComm(),model->verboseAllProc() );

        if ( M_usePseudoTransientContinuation )
            M_pseudoTransientContinuationDeltaAndResidual.clear();
        //---------------------------------------------------------------------//
        //---------------------------------------------------------------------//
        //---------------------------------------------------------------------//
        model->timerTool("Solve").start();
        model->updateNewtonInitialGuess(U);
        for ( auto const& func : M_addFunctionNewtonInitialGuess )
            func.second( U );
        //U->close();
        model->timerTool("Solve").elapsed("algebraic-newton-bc");
        model->timerTool("Solve").restart();
        //---------------------------------------------------------------------//
        if (model->useCstMatrix())
        {
            if (!M_hasBuildLinearJacobian)
            {
                M_CstJ->zero();
                ModelAlgebraic::DataUpdateJacobian dataJacobianCst(U, M_CstJ, M_R, true, M_Extended,false );
                model->updateJacobian( dataJacobianCst );
                for ( auto const& func : M_addFunctionJacobianAssembly )
                    func.second( dataJacobianCst );
                M_CstJ->close();
                M_hasBuildLinearJacobian = true;
            }
            else if (model->rebuildLinearPartInJacobian() || model->needToRebuildCstPart())
            {
                M_CstJ->zero();
                ModelAlgebraic::DataUpdateJacobian dataJacobianCst(U, M_CstJ, M_R, true, M_Extended,false );
                model->updateJacobian( dataJacobianCst );
                for ( auto const& func : M_addFunctionJacobianAssembly )
                    func.second( dataJacobianCst );
                M_CstJ->close();
            }
        }
        //model->updatePreconditioner(U,M_Prec);
        //---------------------------------------------------------------------//
        model->timerTool("Solve").elapsed("algebraic-jacobian");
        model->timerTool("Solve").restart();
        //---------------------------------------------------------------------//
        if (model->useCstVector())
        {
            if (!M_hasBuildResidualCst)
            {
                M_CstR->zero();
                // Warning : the second true is very important in order to build M_CstR!!!!!!
                ModelAlgebraic::DataUpdateResidual dataResidualCst( U, M_CstR, true, true );
                model->updateResidual( dataResidualCst );
                for ( auto const& func : M_addFunctionResidualAssembly )
                    func.second( dataResidualCst );
                M_CstR->close();
                for ( auto const& av : M_addVectorResidualAssembly )
                    if ( std::get<2>( av.second ) )
                        M_CstR->add( std::get<1>( av.second ), std::get<0>( av.second ) );
                M_hasBuildResidualCst = true;
            }
            else if (model->rebuildCstPartInResidual() || model->needToRebuildCstPart() ) // first if is not an option (just optimisation with fsi semi-implicit)
            {
                M_CstR->zero();
                // Warning : the second true is very important in order to build M_CstR!!!!!!
                ModelAlgebraic::DataUpdateResidual dataResidualCst( U, M_CstR, true, true );
                model->updateResidual( dataResidualCst );
                for ( auto const& func : M_addFunctionResidualAssembly )
                    func.second( dataResidualCst );
                M_CstR->close();
                for ( auto const& av : M_addVectorResidualAssembly )
                    if ( std::get<2>( av.second ) )
                        M_CstR->add( std::get<1>( av.second ), std::get<0>( av.second ) );
            }
        }
        model->setNeedToRebuildCstPart(false);
        //---------------------------------------------------------------------//
        model->timerTool("Solve").elapsed("algebraic-residual");
        model->timerTool("Solve").restart();
        //---------------------------------------------------------------------//

        pre_solve_type pre_solve = std::bind(&model_type::preSolveNewton, model, std::placeholders::_1, std::placeholders::_2);
        post_solve_type post_solve = std::bind(&model_type::postSolveNewton, model, std::placeholders::_1, std::placeholders::_2);
        update_nlsolve_type update_nlsolve = std::bind(&self_type::updateNewtonIteration, std::ref( *this ), std::placeholders::_1, std::placeholders::_2,std::placeholders::_3,std::placeholders::_4 );

        auto const solveStat = M_backend->nlSolve( _jacobian=M_J,
                                                   _solution=U,
                                                   _residual=M_R,
                                                   _prec=M_PrecondManage,
                                                   _pre=pre_solve,
                                                   _post=post_solve,
                                                   _update=update_nlsolve );
        if ( false )
            Feel::FeelModels::Log(model->prefix()+".ModelAlgebraicFactory","NonLinearSolverNewton",
                                  "solver stat :\n" +
                                  (boost::format("   - isConverged : %1%\n") %solveStat.isConverged() ).str() +
                                  (boost::format("   - nIterations : %1%\n") %solveStat.nIterations() ).str() +
                                  (boost::format("   - residual : %1%\n") %solveStat.residual() ).str(),
                                  model->worldComm(),model->verboseSolverTimerAllProc());


        if ( model->errorIfSolverNotConverged() )
            CHECK( solveStat.isConverged() ) << "the non-linear solver has not converged in "
                                             << solveStat.nIterations() << " iterations with norm "
                                             << solveStat.residual() << "\n";

        //---------------------------------------------------------------------//
        //---------------------------------------------------------------------//
        //---------------------------------------------------------------------//
        model->timerTool("Solve").elapsed();
        double tElapsed = model->timerTool("Solve").accumulateTime();
        model->timerTool("Solve").setDataValue("algebraic-nlsolve",tElapsed);
        model->timerTool("Solve").setAdditionalParameter("snes-niter",int(solveStat.nIterations()) );
        model->timerTool("Solve").stop();

        if (model->verboseSolverTimer()) Feel::FeelModels::Log( model->prefix()+".ModelAlgebraicFactory","NonLinearSolverNewton",
                                                                (boost::format("finish in %1% s")%tElapsed ).str(),
                                                                model->worldComm(),model->verboseSolverTimerAllProc());
    }

    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::rebuildCstJacobian( vector_ptrtype U )
    {
        M_CstJ->zero();
        ModelAlgebraic::DataUpdateJacobian dataJacobianCst(U, M_CstJ, M_R, true, M_Extended,false );
        this->model()->updateJacobian( dataJacobianCst );
        for ( auto const& func : M_addFunctionJacobianAssembly )
            func.second( dataJacobianCst );
        M_CstJ->close();
    }


    void
    ModelAlgebraicFactory::rebuildCstLinearPDE( vector_ptrtype U )
    {
        M_CstJ->zero();
        M_CstR->zero();
        ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true,M_Extended,false);
        this->model()->updateLinearPDE(dataLinearCst);
        for ( auto const& func : M_addFunctionLinearAssembly )
            func.second( dataLinearCst );
        M_CstJ->close();
        M_CstR->close();
    }





    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//


    void
    ModelAlgebraicFactory::solvePicard( vector_ptrtype& U )
    {
        if (this->model()->verbose()) Feel::FeelModels::Log(this->model()->prefix()+".ModelAlgebraicFactory","AlgoPicard", "start",
                                               this->model()->worldComm(),this->model()->verboseAllProc());

        this->model()->timerTool("Solve").start();

        this->model()->timerTool("Solve").setDataValue("algebraic-assembly",0.);
        this->model()->timerTool("Solve").setDataValue("algebraic-solve",0.);

        bool useConvergenceAlgebraic = true;

        // assembling cst part
        if ( !M_hasBuildLinearSystemCst ||
             this->model()->rebuildCstPartInLinearSystem() || this->model()->needToRebuildCstPart() ||
             !this->model()->useCstMatrix() || !this->model()->useCstVector() )
        {
            this->model()->timerTool("Solve").start();

            M_CstJ->zero();
            M_CstR->zero();
            ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true,M_Extended,false);
            this->model()->updateLinearPDE( dataLinearCst );
            M_hasBuildLinearSystemCst = true;

            double tAssemblyElapsed = this->model()->timerTool("Solve").stop();
            this->model()->timerTool("Solve").addDataValue("algebraic-assembly",tAssemblyElapsed);
        }
        this->model()->setNeedToRebuildCstPart(false);

        vector_ptrtype Uold;
        vector_ptrtype residual;
        if ( !useConvergenceAlgebraic )
            Uold = M_backend->newVector( M_R->mapPtr() );
        else
            residual = M_backend->newVector( M_R->mapPtr() );

        bool hasConverged = false;

        double convergenceRate=1;
        double rtol = this->backend()->rToleranceSNES();//1e-6;
        double atol = this->backend()->aToleranceSNES();//1e-50;
        int fixPointMaxIt = this->backend()->maxIterationsSNES();
        int cptIteration=0;
        for ( ; cptIteration < fixPointMaxIt ; ++cptIteration )
        {
            if ( !useConvergenceAlgebraic )
            {
                *Uold = *U;
                Uold->close();
            }

            this->model()->timerTool("Solve").start();

            // copying cst matrix/vector
            if (this->model()->useCstMatrix() && this->model()->useCstVector() )
            {
                M_J->zero();
                M_J->addMatrix(1.0, M_CstJ );
                M_R->zero();
                M_R->add(1.0, M_CstR );
            }
            else if ( cptIteration > 0 )
            {
                M_J->zero();
                M_R->zero();
                ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_J,M_R,true,M_Extended,false);
                this->model()->updateLinearPDE( dataLinearCst );
            }

            // assembling non cst part
            ModelAlgebraic::DataUpdateLinear dataLinearNonCst(U,M_J,M_R,false,M_Extended,true);
            this->model()->updateLinearPDE( dataLinearNonCst );

            // dof elimination
            this->model()->updateLinearPDEDofElimination( dataLinearNonCst );

            // post-assembly
            for ( auto const& func : M_addFunctionLinearPostAssembly )
                func.second( dataLinearNonCst );

            double tAssemblyElapsed = this->model()->timerTool("Solve").elapsed();
            this->model()->timerTool("Solve").addDataValue("algebraic-assembly",tAssemblyElapsed);

            if (this->model()->verboseSolverTimer())
                Feel::FeelModels::Log(this->model()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                      (boost::format("picard iteration[%1%] finish assembling in %2% s") %cptIteration % tAssemblyElapsed ).str(),
                                      this->model()->worldComm(),this->model()->verboseSolverTimerAllProc());

            this->model()->timerTool("Solve").restart();

            if ( useConvergenceAlgebraic )
            {
                M_J->close();
                M_R->close();
                residual->zero();
                residual->addVector( U,M_J );
                residual->add( -1., M_R );
                residual->close();
                double nomResidual = residual->l2Norm();
                double normRhs = M_R->l2Norm();
                if (this->model()->verboseSolverTimer())
                    Feel::FeelModels::Log( this->model()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                           (boost::format("picard iteration[%1%] residual norm : %2%")%cptIteration %nomResidual ).str(),
                                           this->model()->worldComm(),this->model()->verboseSolverTimerAllProc() );

                if ( nomResidual < std::max(rtol*normRhs,atol) )
                {
                    hasConverged = true;
                    this->model()->timerTool("Solve").stop();
                    break;
                }
            }

            // assembling matrix used for preconditioner
            this->model()->updatePreconditioner(U,M_J,M_Extended,M_Prec);

            // update in-house preconditioners
            this->model()->updateInHousePreconditioner( M_J, U );

            // set pre/post solve functions
            pre_solve_type pre_solve = std::bind(&model_type::preSolvePicard, this->model(), std::placeholders::_1, std::placeholders::_2);
            post_solve_type post_solve = std::bind(&model_type::postSolvePicard, this->model(), std::placeholders::_1, std::placeholders::_2);
            // solve linear system
            auto const solveStat = M_backend->solve( _matrix=M_J,
                                                     _solution=U,
                                                     _rhs=M_R,
                                                     _prec=M_PrecondManage,
                                                     _pre=pre_solve,
                                                     _post=post_solve );

            if ( this->model()->errorIfSolverNotConverged() )
                CHECK( solveStat.isConverged() ) << "the linear solver has not converged in "
                                                 << solveStat.nIterations() << " iterations with norm "
                                                 << solveStat.residual() << "\n";

            double tSolveElapsed = this->model()->timerTool("Solve").stop();
            this->model()->timerTool("Solve").addDataValue("algebraic-solve",tSolveElapsed);

            //this->model()->timerTool("Solve").setAdditionalParameter("ksp-niter",int(solveStat.nIterations()) );
            if (this->model()->verboseSolverTimer())
                Feel::FeelModels::Log(this->model()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                      (boost::format("picard iteration[%1%] finish sub solve in %2% s")%cptIteration %tSolveElapsed ).str(),
                                      this->model()->worldComm(),this->model()->verboseSolverTimerAllProc());

            if ( !useConvergenceAlgebraic )
            {
                convergenceRate = this->model()->updatePicardConvergence( U,Uold );
                if (this->model()->verboseSolverTimer())
                    Feel::FeelModels::Log( this->model()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                           (boost::format("fix point convergence : %1%")%convergenceRate ).str(),
                                           this->model()->worldComm(),this->model()->verboseSolverTimerAllProc() );

                if ( convergenceRate < rtol )
                {
                    hasConverged = true;
                    break;
                }

            }
        } // for ( ; cptIteration < fixPointMaxIt ; ++cptIteration )

        double tFixPointElapsed = this->model()->timerTool("Solve").stop();
        //this->model()->timerTool("Solve").setAdditionalParameter("ksp-niter",int(solveStat.nIterations()) );

        if (this->model()->verboseSolverTimer())
            Feel::FeelModels::Log(this->model()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                  (boost::format("finish solve in %1% s")%tFixPointElapsed ).str(),
                                  this->model()->worldComm(),this->model()->verboseSolverTimerAllProc());

    }

void
ModelAlgebraicFactory::updateNewtonIteration( int step, vector_ptrtype residual, vector_ptrtype sol, typename backend_type::solvernonlinear_type::UpdateIterationData const& data )
{
    this->model()->updateNewtonIteration( step, residual, sol, data );

    if ( M_usePseudoTransientContinuation )
    {
        /**
         * some references used :
         * - Coffey, T. S., Kelley, C. T., & Keyes, D. E. (2003). Pseudotransient continuation and differential-algebraic equations. SIAM Journal on Scientific Computing, 25(2), 553-569.
         * - Kelley, C. T., & Keyes, D. E. (1998). Convergence analysis of pseudo-transient continuation. SIAM Journal on Numerical Analysis, 35(2), 508-523.
         * - Ceze, M., & Fidkowski, K. (2013). Pseudo-transient continuation, solution update methods, and CFL strategies for DG discretizations of the RANS-SA equations. In 21st AIAA computational fluid dynamics conference (p. 2686).
         * - Mavriplis, D. (2018). A Residual Smoothing Strategy for Accelerating Newton Method Continuation. arXiv preprint arXiv:1805.03756.
         */
        if ( M_pseudoTransientContinuationEvolutionMethod == "SER" )
        {
            if ( M_pseudoTransientContinuationSerVariant == "residual" )
            {
                double resNorm = residual->l2Norm();
                if ( step == 0 )
                    M_pseudoTransientContinuationDeltaAndResidual.push_back( std::make_pair(M_pseudoTransientContinuationDelta0,resNorm ) );
                else
                {
                    double curDelta = M_pseudoTransientContinuationDeltaAndResidual.back().first*M_pseudoTransientContinuationDeltaAndResidual.back().second/resNorm;
                    double newDelta = std::min( curDelta, M_pseudoTransientContinuationDeltaMax );
                    M_pseudoTransientContinuationDeltaAndResidual.push_back( std::make_pair(newDelta,resNorm ) );
                }
            }
            else if ( M_pseudoTransientContinuationSerVariant == "solution" )
            {
                if ( step == 0 )
                    M_pseudoTransientContinuationDeltaAndResidual.push_back( std::make_pair(M_pseudoTransientContinuationDelta0,0. ) );
                else
                {
                    M_pseudoTransientContinuationPreviousSolution->add( -1.,sol );
                    double diffNorm = M_pseudoTransientContinuationPreviousSolution->l2Norm();
                    double curDelta = M_pseudoTransientContinuationDeltaAndResidual.back().first/diffNorm;
                    double newDelta = std::min( curDelta, M_pseudoTransientContinuationDeltaMax );
                }
                *M_pseudoTransientContinuationPreviousSolution = *sol;
            }
        }
        else if ( M_pseudoTransientContinuationEvolutionMethod == "EXPur" )
        {
            double lambda = data.doubleInfos("linesearch.lambda" );
            //std::cout << "QQ lambda="<< lambda << "\n";
            if ( step == 0 )
                M_pseudoTransientContinuationDeltaAndResidual.push_back( std::make_pair(M_pseudoTransientContinuationDelta0,0. ) );
            else
            {
                double beta1 = M_pseudoTransientContinuationExpurBetaHigh;
                double beta2 = M_pseudoTransientContinuationExpurBetaLow;
                double lastDelta = M_pseudoTransientContinuationDeltaAndResidual.back().first;
                if ( lambda >= M_pseudoTransientContinuationExpurThresholdHigh )
                    M_pseudoTransientContinuationDeltaAndResidual.push_back( std::make_pair( beta1*lastDelta,0. ) );
                else if ( lambda <= M_pseudoTransientContinuationExpurThresholdLow )
                    M_pseudoTransientContinuationDeltaAndResidual.push_back( std::make_pair( beta2*lastDelta,0. ) );
                else
                    M_pseudoTransientContinuationDeltaAndResidual.push_back( std::make_pair( lastDelta,0. ) );
            }
        }
        //std::cout << "CFL="<<M_pseudoTransientContinuationDeltaAndResidual.back().first<<"\n";
    }
}


} // namespace FeelModels
} // namespace Feel
