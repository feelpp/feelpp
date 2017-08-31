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

    ModelAlgebraicFactory::ModelAlgebraicFactory(appli_ptrtype const& __app, backend_ptrtype const& __backend)
        :
        M_appli(__app),
        M_backend( __backend ),
        M_hasBuildLinearJacobian(false),
        M_hasBuildResidualCst(false),
        M_hasBuildLinearSystemCst(false)
    {
        this->application()->timerTool("Constructor").start();
        auto graph = M_appli->buildMatrixGraph();
        this->application()->timerTool("Constructor").elapsed("graph");
        if ( graph )
        {
            this->application()->timerTool("Constructor").restart();
            this->buildMatrixVector( graph,graph->mapRow().indexSplit() );
            this->application()->timerTool("Constructor").elapsed("matrixVector");

            this->application()->timerTool("Constructor").restart();
            this->buildOthers();
            this->application()->timerTool("Constructor").elapsed("algebraicOthers");
        }
        this->application()->timerTool("Constructor").stop();
    }

    //---------------------------------------------------------------------------------------------------------------//

    ModelAlgebraicFactory::ModelAlgebraicFactory(appli_ptrtype const& __app,
                           backend_ptrtype const& __backend,
                           graph_ptrtype const& graph,
                           indexsplit_ptrtype const& indexSplit )
        :
        M_appli(__app),
        M_backend(__backend ),
        M_hasBuildLinearJacobian(false),
        M_hasBuildResidualCst(false),
        M_hasBuildLinearSystemCst(false)
    {
        if (this->application()->verbose()) Feel::FeelModels::Log(this->application()->prefix()+".MethodNum","constructor1", "start",
                                               this->application()->worldComm(),this->application()->verboseAllProc());

        this->init(graph,indexSplit);

        if (this->application()->verbose()) Feel::FeelModels::Log(this->application()->prefix()+".MethodNum","constructor", "finish",
                                               this->application()->worldComm(),this->application()->verboseAllProc());
    }

    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::init(graph_ptrtype const& graph,
                     indexsplit_ptrtype const& indexSplit)
    {
        this->application()->timerTool("Constructor").start();
        this->buildMatrixVector(graph,indexSplit);
        this->application()->timerTool("Constructor").elapsed("matrixVector");

        this->application()->timerTool("Constructor").restart();
        this->buildOthers();
        this->application()->timerTool("Constructor").stop("algebraicOthers");
    }

    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::buildMatrixVector(graph_ptrtype const& graph,
                                  indexsplit_ptrtype const& indexSplit)
    {
        //vectors
        M_R = M_backend->newVector(graph->mapRowPtr()); //  size1,size1 ;
        if (this->application()->useCstVector())
            M_CstR = M_backend->newVector(graph->mapRowPtr());//size1,size1 );
        else M_CstR = M_R;
        if (this->application()->verbose()) Feel::FeelModels::Log(this->application()->prefix()+".MethodNum","buildMatrixVector","build all vectors",
                                                           this->application()->worldComm(),this->application()->verboseAllProc());

        size_type size1=0,size2=0;
        //matrix with standart pattern
        M_J = M_backend->newMatrix(size1,size2,size1,size2,graph,indexSplit);
        if (this->application()->useCstMatrix())
            M_CstJ = M_backend->newMatrix(size1,size2,size1,size2,graph,indexSplit);
        else M_CstJ = M_J;
        if (this->application()->buildMatrixPrecond())
            M_Prec = M_backend->newMatrix( size1,size2,size1,size2,graph,indexSplit );
        else M_Prec = M_J;
        //matrix with extended pattern
        if (this->application()->hasExtendedPattern() && this->application()->buildMatrixPrecond())
            M_Extended = M_backend->newMatrix( size1,size2,size1,size2,graph,indexSplit );
        else
            M_Extended = M_J;

        if (this->application()->verbose()) Feel::FeelModels::Log(this->application()->prefix()+".MethodNum","buildMatrixVector","build all matrix",
                                                           this->application()->worldComm(),this->application()->verboseAllProc());

        // M_J = M_backend->newMatrix( _trial=M_Xh, _test=M_Xh, _pattern_block=myblockpattern );
        // M_Extended = M_backend->newZeroMatrix( size1,size2,size1,size2 );

        if (this->application()->printGraph())
            graph->printPython( this->application()->printGraphFileName() );
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
        M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian,
                                                       boost::ref( *this ), _1, _2/*, M_R*/ );
        M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual,
                                                       boost::ref( *this ), _1, _2 );

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
        boost::shared_ptr<NullSpace<value_type> > mynullspace( new NullSpace<value_type>(this->backend(),nullSpace) );
        this->backend()->attachNullSpace( mynullspace );
    }
    void
    ModelAlgebraicFactory::attachNearNullSpace( NullSpace<value_type> const& nearNullSpace )
    {
        CHECK( this->backend() ) << "backend not init\n";
        boost::shared_ptr<NullSpace<value_type> > myNearNullSpace( new NullSpace<value_type>(this->backend(),nearNullSpace) );
        this->backend()->attachNearNullSpace( myNearNullSpace );
    }
    void
    ModelAlgebraicFactory::attachNearNullSpace( int k, NullSpace<value_type> const& nearNullSpace )
    {
        CHECK( this->backend() ) << "backend not init\n";
        CHECK( M_PrecondManage ) << "preconditioner not init\n";
        boost::shared_ptr<NullSpace<value_type> > myNearNullSpace( new NullSpace<value_type>(this->backend(),nearNullSpace) );
        M_PrecondManage->attachNearNullSpace( k, myNearNullSpace );
    }


    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//

    boost::shared_ptr<std::ostringstream>
    ModelAlgebraicFactory::getInfo() const
    {
        boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );

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
        if ( this->application()->verboseAllProc() ) std::cout << this->getInfo()->str();
        else if (this->application()->worldComm().isMasterRank() )
            std::cout << this->getInfo()->str();
    }


    void
    ModelAlgebraicFactory::solve( std::string const& type,vector_ptrtype& sol )
    {
        if ( type == "LinearSystem" )
        {
            this->linearSolver( sol );
        }
        else if ( type == "Picard" || type == "FixPoint" )
        {
            this->AlgoPicard( sol );
        }
        else if ( type == "Newton" )
        {
            this->AlgoNewton2( sol );
        }
        else
            CHECK( false ) << "invalid solver type " << type;
    }

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::linearSolver( vector_ptrtype& U )
    {
        if (this->application()->verbose()) Feel::FeelModels::Log(this->application()->prefix()+".ModelAlgebraicFactory","linearSolver", "start",
                                               this->application()->worldComm(),this->application()->verboseAllProc());

        this->application()->timerTool("Solve").start();

        // assembling cst part
        if (!M_hasBuildLinearSystemCst)
        {
            M_CstJ->zero();
            M_CstR->zero();
            //this->application()->updateLinearPDE(U,M_CstJ,M_CstR,true,M_Extended,false);
            ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true,M_Extended,false);
            this->application()->updateLinearPDE( dataLinearCst );
            M_hasBuildLinearSystemCst = true;
        }
        else if ( this->application()->rebuildCstPartInLinearSystem() || this->application()->needToRebuildCstPart() ||
                  !this->application()->useCstMatrix() || !this->application()->useCstVector() )
        {
            M_CstJ->zero();
            M_CstR->zero();
            //this->application()->updateLinearPDE(U,M_CstJ,M_CstR,true,M_Extended,false);
            ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true,M_Extended,false);
            this->application()->updateLinearPDE( dataLinearCst );
        }
        this->application()->setNeedToRebuildCstPart(false);

        // copying cst matrix/vector
        if (this->application()->useCstMatrix())
        {
            M_J->zero();
            M_J->addMatrix(1.0, M_CstJ );
        }
        if (this->application()->useCstVector())
        {
            M_R->zero();
            M_R->add(1.0, M_CstR );
        }

        // set M_Extended to zero if really used and build
        if (this->application()->hasExtendedPattern() && this->application()->buildMatrixPrecond() )
            M_Extended->zero();

        // pre-assembly (optional)
        if ( this->addFunctionLinearPreAssemblyNonCst != NULL )
            this->addFunctionLinearPreAssemblyNonCst( M_J,M_R );

        // assembling non cst part
        //this->application()->updateLinearPDE(U,M_J,M_R,false,M_Extended,true);
        ModelAlgebraic::DataUpdateLinear dataLinearNonCst(U,M_J,M_R,false,M_Extended,true);
        this->application()->updateLinearPDE( dataLinearNonCst );

        // post-assembly (optional)
        if ( this->addFunctionLinearPostAssembly != NULL )
            this->addFunctionLinearPostAssembly(M_J,M_R);

        // assembling matrix used for preconditioner
        this->application()->updatePreconditioner(U,M_J,M_Extended,M_Prec);

        // add extended term (important : after updatePreconditioner !)
        // and only do if shared_ptr M_J != M_Extended
        if (this->application()->hasExtendedPattern() && this->application()->buildMatrixPrecond() )
            M_J->addMatrix(1.0, M_Extended );

        double mpiTimerAssembly = this->application()->timerTool("Solve").elapsed("algebraic-assembly");

        if (this->application()->verboseSolverTimer())
          Feel::FeelModels::Log(this->application()->prefix()+".ModelAlgebraicFactory","linearSolver",
                         (boost::format("finish assembling in %1% s") % mpiTimerAssembly ).str(),
                         this->application()->worldComm(),this->application()->verboseSolverTimerAllProc());

        this->application()->timerTool("Solve").restart();

        // set preconditioner
        //M_PrecondManage->setMatrix(M_Prec);

        this->application()->updateInHousePreconditioner( M_J, U );

        // solve linear system
        auto const solveStat = M_backend->solve( _matrix=M_J,
                                                 _solution=U,
                                                 _rhs=M_R,
                                                 _prec=M_PrecondManage );

        if ( this->application()->errorIfSolverNotConverged() )
            CHECK( solveStat.isConverged() ) << "the linear solver has not converged in "
                                             << solveStat.nIterations() << " iterations with norm "
                                             << solveStat.residual() << "\n";

        double tElapsed = this->application()->timerTool("Solve").stop("algebraic-solve");
        this->application()->timerTool("Solve").setAdditionalParameter("ksp-niter",int(solveStat.nIterations()) );
#if 0
        if ( option(_name="clear-preconditioner-after-use",_prefix=this->application()->prefix()).as<bool>() )
        {
            Feel::FeelModels::Log(this->application()->prefix()+".ModelAlgebraicFactory","linearSolver","clear-preconditioner-after-use",
                           this->application()->worldComm(),this->application()->verboseSolverTimerAllProc());
            M_PrecondManage->clear();
            MatrixPetsc<double> * pmatrix = dynamic_cast<MatrixPetsc<double>*>( M_J.get() );
            pmatrix->mapSplitPC().clear();
        }
#endif
        if (this->application()->verboseSolverTimer())
          Feel::FeelModels::Log(this->application()->prefix()+".ModelAlgebraicFactory","linearSolver",
                         (boost::format("finish solve in %1% s")%tElapsed ).str(),
                         this->application()->worldComm(),this->application()->verboseSolverTimerAllProc());
    }



    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//


    void
    ModelAlgebraicFactory::updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J/*, vector_ptrtype& R*/ )
    {
        this->application()->timerTool("Solve").start();

        auto R = this->backend()->newVector(X->mapPtr());

        J->zero();

        if ( this->addFunctionJacobianPreAssembly != NULL )
            this->addFunctionJacobianPreAssembly( X, J );

        if (this->application()->useCstMatrix())
            {
                J->addMatrix(1.0, M_CstJ );
            }
        else
            {
                //this->application()->updateJacobian( X, J, R, true, M_Extended,false );
                ModelAlgebraic::DataUpdateJacobian dataJacobianCst(X, J, R, true, M_Extended,false);
                this->application()->updateJacobian( dataJacobianCst );
            }

        //M_appli->updateJacobian(X,J,R,false, M_Extended,false);
        ModelAlgebraic::DataUpdateJacobian dataJacobianNonCst(X,J,R,false, M_Extended,false);
        M_appli->updateJacobian( dataJacobianNonCst );

        this->application()->updateInHousePreconditioner( J, X );

        double tElapsed = this->application()->timerTool("Solve").stop();
        this->application()->timerTool("Solve").addDataValue("algebraic-jacobian",tElapsed);
    }
    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::updateResidual(const vector_ptrtype& X, vector_ptrtype& R)
    {
        this->application()->timerTool("Solve").start();

        R->zero();

        if ( this->addFunctionResidualPreAssembly != NULL )
            this->addFunctionResidualPreAssembly( X, R );

        if (this->application()->useCstVector())
            {
                R->add(1.0, M_CstR );
            }
        else
            {
                //M_appli->updateResidual( X, R, true, true );
                ModelAlgebraic::DataUpdateResidual dataResidualCst( X, R, true, true );
                M_appli->updateResidual( dataResidualCst );
            }

        bool doOptimization = this->application()->useLinearJacobianInResidual() && this->application()->useCstMatrix();
        // add linear contribution from jacobian terms
        if (doOptimization)
            R->addVector(*X, *M_CstJ );

        //M_appli->updateResidual(X,R,false,doOptimization);
        ModelAlgebraic::DataUpdateResidual dataResidualNonCst( X, R, false, doOptimization );
        M_appli->updateResidual( dataResidualNonCst );

        double tElapsed = this->application()->timerTool("Solve").stop();
        this->application()->timerTool("Solve").addDataValue("algebraic-residual",tElapsed);
    }

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::AlgoNewton2(vector_ptrtype& U)
    {
        if (this->application()->verbose()) Feel::FeelModels::Log(this->application()->prefix()+".ModelAlgebraicFactory","NonLinearSolverNewton", "start",
                                               this->application()->worldComm(),this->application()->verboseAllProc());

        //---------------------------------------------------------------------//
        //---------------------------------------------------------------------//
        //---------------------------------------------------------------------//
        this->application()->timerTool("Solve").start();
        M_appli->updateNewtonInitialGuess(U);
        //U->close();
        this->application()->timerTool("Solve").elapsed("algebraic-newton-bc");
        this->application()->timerTool("Solve").restart();
        //---------------------------------------------------------------------//
        if (this->application()->useCstMatrix())
            {
                if (!M_hasBuildLinearJacobian)
                    {
                        M_CstJ->zero();
                        //M_appli->updateJacobian( U, M_CstJ, M_R, true, M_Extended,false );
                        ModelAlgebraic::DataUpdateJacobian dataJacobianCst(U, M_CstJ, M_R, true, M_Extended,false );
                        M_appli->updateJacobian( dataJacobianCst );
                        M_CstJ->close();
                        M_hasBuildLinearJacobian = true;
                    }
                else if (this->application()->rebuildLinearPartInJacobian() || this->application()->needToRebuildCstPart())
                    {
                        M_CstJ->zero();
                        //M_appli->updateJacobian( U, M_CstJ, M_R, true, M_Extended,false );
                        ModelAlgebraic::DataUpdateJacobian dataJacobianCst(U, M_CstJ, M_R, true, M_Extended,false );
                        M_appli->updateJacobian( dataJacobianCst );
                        M_CstJ->close();
                    }
            }
        //M_appli->updatePreconditioner(U,M_Prec);
        //---------------------------------------------------------------------//
        this->application()->timerTool("Solve").elapsed("algebraic-jacobian");
        this->application()->timerTool("Solve").restart();
        //---------------------------------------------------------------------//
        if (this->application()->useCstVector())
            {
                if (!M_hasBuildResidualCst)
                    {
                        M_CstR->zero();
                        // Warning : the second true is very important in order to build M_CstR!!!!!!
                        //M_appli->updateResidual( U, M_CstR, true, true );
                        ModelAlgebraic::DataUpdateResidual dataResidualCst( U, M_CstR, true, true );
                        M_appli->updateResidual( dataResidualCst );
                        M_CstR->close();
                        M_hasBuildResidualCst = true;
                    }
                else if (this->application()->rebuildCstPartInResidual() || this->application()->needToRebuildCstPart() ) // first if is not an option (just optimisation with fsi semi-implicit)
                    {
                        M_CstR->zero();
                        // Warning : the second true is very important in order to build M_CstR!!!!!!
                        //M_appli->updateResidual( U, M_CstR, true, true );
                        ModelAlgebraic::DataUpdateResidual dataResidualCst( U, M_CstR, true, true );
                        M_appli->updateResidual( dataResidualCst );
                        M_CstR->close();
                    }
            }
        this->application()->setNeedToRebuildCstPart(false);
        //---------------------------------------------------------------------//
        this->application()->timerTool("Solve").elapsed("algebraic-residual");
        this->application()->timerTool("Solve").restart();
        //---------------------------------------------------------------------//
        
        pre_solve_type pre_solve = std::bind(&appli_type::preSolveNewton, M_appli, std::placeholders::_1, std::placeholders::_2);
        post_solve_type post_solve = std::bind(&appli_type::postSolveNewton, M_appli, std::placeholders::_1, std::placeholders::_2);
        
        auto const solveStat = M_backend->nlSolve( _jacobian=M_J,
                                                   _solution=U,
                                                   _residual=M_R,
                                                   _prec=M_PrecondManage,
                                                   _pre=pre_solve,
                                                   _post=post_solve);
        if ( false )
            Feel::FeelModels::Log(this->application()->prefix()+".ModelAlgebraicFactory","NonLinearSolverNewton",
                           "solver stat :\n" +
                           (boost::format("   - isConverged : %1%\n") %solveStat.isConverged() ).str() +
                           (boost::format("   - nIterations : %1%\n") %solveStat.nIterations() ).str() +
                           (boost::format("   - residual : %1%\n") %solveStat.residual() ).str(),
                           this->application()->worldComm(),this->application()->verboseSolverTimerAllProc());


        if ( this->application()->errorIfSolverNotConverged() )
            CHECK( solveStat.isConverged() ) << "the non-linear solver has not converged in "
                                             << solveStat.nIterations() << " iterations with norm "
                                             << solveStat.residual() << "\n";

        //---------------------------------------------------------------------//
        //---------------------------------------------------------------------//
        //---------------------------------------------------------------------//
        this->application()->timerTool("Solve").elapsed();
        double tElapsed = this->application()->timerTool("Solve").accumulateTime();
        this->application()->timerTool("Solve").setDataValue("algebraic-nlsolve",tElapsed);
        this->application()->timerTool("Solve").setAdditionalParameter("snes-niter",int(solveStat.nIterations()) );
        this->application()->timerTool("Solve").stop();

        if (this->application()->verboseSolverTimer()) Feel::FeelModels::Log( this->application()->prefix()+".ModelAlgebraicFactory","NonLinearSolverNewton",
                                                                       (boost::format("finish in %1% s")%tElapsed ).str(),
                                                                       this->application()->worldComm(),this->application()->verboseSolverTimerAllProc());
    }

    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::rebuildCstJacobian( vector_ptrtype U )
    {
        M_CstJ->zero();
        //M_appli->updateJacobian( U, M_CstJ, M_R, true, M_Extended,false );
        ModelAlgebraic::DataUpdateJacobian dataJacobianCst(U, M_CstJ, M_R, true, M_Extended,false );
        M_appli->updateJacobian( dataJacobianCst );
        M_CstJ->close();
    }


    void
    ModelAlgebraicFactory::rebuildCstLinearPDE( vector_ptrtype U )
    {
        M_CstJ->zero();
        M_CstR->zero();
        ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true,M_Extended,false);
        M_appli->updateLinearPDE(dataLinearCst);
        M_CstJ->close();
        M_CstR->close();
    }





    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//


    void
    ModelAlgebraicFactory::AlgoPicard( vector_ptrtype U )
    {
        if (this->application()->verbose()) Feel::FeelModels::Log(this->application()->prefix()+".ModelAlgebraicFactory","AlgoPicard", "start",
                                               this->application()->worldComm(),this->application()->verboseAllProc());

        this->application()->timerTool("Solve").start();

        this->application()->timerTool("Solve").setDataValue("algebraic-assembly",0.);
        this->application()->timerTool("Solve").setDataValue("algebraic-solve",0.);

        bool useConvergenceAlgebraic = true;

        // assembling cst part
        if ( !M_hasBuildLinearSystemCst ||
             this->application()->rebuildCstPartInLinearSystem() || this->application()->needToRebuildCstPart() ||
             !this->application()->useCstMatrix() || !this->application()->useCstVector() )
        {
            this->application()->timerTool("Solve").start();

            M_CstJ->zero();
            M_CstR->zero();
            ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true,M_Extended,false);
            this->application()->updateLinearPDE( dataLinearCst );
            M_hasBuildLinearSystemCst = true;

            double tAssemblyElapsed = this->application()->timerTool("Solve").stop();
            this->application()->timerTool("Solve").addDataValue("algebraic-assembly",tAssemblyElapsed);
        }
        this->application()->setNeedToRebuildCstPart(false);

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

            this->application()->timerTool("Solve").start();

            // copying cst matrix/vector
            if (this->application()->useCstMatrix() && this->application()->useCstVector() )
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
                this->application()->updateLinearPDE( dataLinearCst );
            }

            // pre-assembly (optional)
            if ( this->addFunctionLinearPreAssemblyNonCst != NULL )
                this->addFunctionLinearPreAssemblyNonCst( M_J,M_R );

            // assembling non cst part
            ModelAlgebraic::DataUpdateLinear dataLinearNonCst(U,M_J,M_R,false,M_Extended,true);
            this->application()->updateLinearPDE( dataLinearNonCst );

            // post-assembly (optional)
            if ( this->addFunctionLinearPostAssembly != NULL )
                this->addFunctionLinearPostAssembly(M_J,M_R);

            double tAssemblyElapsed = this->application()->timerTool("Solve").elapsed();
            this->application()->timerTool("Solve").addDataValue("algebraic-assembly",tAssemblyElapsed);

            if (this->application()->verboseSolverTimer())
                Feel::FeelModels::Log(this->application()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                      (boost::format("picard iteration[%1%] finish assembling in %2% s") %cptIteration % tAssemblyElapsed ).str(),
                                      this->application()->worldComm(),this->application()->verboseSolverTimerAllProc());

            this->application()->timerTool("Solve").restart();

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
                if (this->application()->verboseSolverTimer())
                    Feel::FeelModels::Log( this->application()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                           (boost::format("picard iteration[%1%] residual norm : %2%")%cptIteration %nomResidual ).str(),
                                           this->application()->worldComm(),this->application()->verboseSolverTimerAllProc() );

                if ( nomResidual < std::max(rtol*normRhs,atol) )
                {
                    hasConverged = true;
                    this->application()->timerTool("Solve").stop();
                    break;
                }
            }

            // assembling matrix used for preconditioner
            this->application()->updatePreconditioner(U,M_J,M_Extended,M_Prec);

            // update in-house preconditioners
            this->application()->updateInHousePreconditioner( M_J, U );

            // set pre/post solve functions
            pre_solve_type pre_solve = std::bind(&appli_type::preSolvePicard, M_appli, std::placeholders::_1, std::placeholders::_2);
            post_solve_type post_solve = std::bind(&appli_type::postSolvePicard, M_appli, std::placeholders::_1, std::placeholders::_2);
            // solve linear system
            auto const solveStat = M_backend->solve( _matrix=M_J,
                                                     _solution=U,
                                                     _rhs=M_R,
                                                     _prec=M_PrecondManage,
                                                     _pre=pre_solve,
                                                     _post=post_solve );

            if ( this->application()->errorIfSolverNotConverged() )
                CHECK( solveStat.isConverged() ) << "the linear solver has not converged in "
                                                 << solveStat.nIterations() << " iterations with norm "
                                                 << solveStat.residual() << "\n";

            double tSolveElapsed = this->application()->timerTool("Solve").stop();
            this->application()->timerTool("Solve").addDataValue("algebraic-solve",tSolveElapsed);

            //this->application()->timerTool("Solve").setAdditionalParameter("ksp-niter",int(solveStat.nIterations()) );
            if (this->application()->verboseSolverTimer())
                Feel::FeelModels::Log(this->application()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                      (boost::format("picard iteration[%1%] finish sub solve in %2% s")%cptIteration %tSolveElapsed ).str(),
                                      this->application()->worldComm(),this->application()->verboseSolverTimerAllProc());

            if ( !useConvergenceAlgebraic )
            {
                convergenceRate = this->application()->updatePicardConvergence( U,Uold );
                if (this->application()->verboseSolverTimer())
                    Feel::FeelModels::Log( this->application()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                           (boost::format("fix point convergence : %1%")%convergenceRate ).str(),
                                           this->application()->worldComm(),this->application()->verboseSolverTimerAllProc() );

                if ( convergenceRate < rtol )
                {
                    hasConverged = true;
                    break;
                }

            }
        } // for ( ; cptIteration < fixPointMaxIt ; ++cptIteration )

        double tFixPointElapsed = this->application()->timerTool("Solve").stop();
        //this->application()->timerTool("Solve").setAdditionalParameter("ksp-niter",int(solveStat.nIterations()) );

        if (this->application()->verboseSolverTimer())
            Feel::FeelModels::Log(this->application()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                  (boost::format("finish solve in %1% s")%tFixPointElapsed ).str(),
                                  this->application()->worldComm(),this->application()->verboseSolverTimerAllProc());

    }


    //---------------------------------------------------------------------------------------------------------------//
    //OLD version (without petsc) : not up
    void
    ModelAlgebraicFactory::AlgoNewton(vector_ptrtype U)
    {
#if 0
        std::cout << "[ModelAlgebraicFactory] : Newton Algo start\n";

        //M_appli->updateNewtonInitialGuess(U);


        //vector_ptrtype R( M_backend->newVector( M_Xh ) );
        //sparse_matrix_ptrtype J( M_backend->newMatrix(M_Xh,M_Xh) );
        vector_ptrtype DeltaU( M_backend->newVector( M_Xh ) );

        vector_ptrtype Uold( M_backend->newVector( M_Xh ) );
        vector_ptrtype Rold( M_backend->newVector( M_Xh ) );

        double Err=1;double ErrDeltaU=1;double ErrDiff=1;
        double ErrTol=1e-8;//1e-7//1e-10;//8;//1e-5;
        double ErrTolDiff=1e-8;
        uint cpt=0;

        //compute residual
        M_appli->updateResidual( U, M_R );
        auto nResidu=M_R->l2Norm();
        while ( Err>ErrTol && ErrDiff>ErrTolDiff && cpt<50)
            {

                vector_ptrtype BBB( M_backend->newVector( M_Xh ) );
                sparse_matrix_ptrtype JJJJ( M_backend->newMatrix(M_Xh,M_Xh) );

                *Uold=*U;

                auto nDeltaUOld=DeltaU->l2Norm();

                M_R->scale(-1.0);

                M_appli->updateJacobian( U, M_J, M_R );

                if (cpt>0) {*Rold=*M_R;if (nResidu!=0) Rold->scale(1./nResidu);ErrDiff=Rold->l2Norm();}
                std::cout << "\n Residu " << ErrDiff << "\n";

                std::cout << "[ModelAlgebraicFactory] : Newton solve start\n";
                /*M_backend->solve( _matrix=J,
                  _solution=DeltaU,
                  _rhs=R,
                  _maxit=1000
                  //_atolerance=1e-50,
                  //_dtolerance=1e7
                  );*/
                M_backend->solve( M_J, DeltaU, M_R );

                std::cout << "[ModelAlgebraicFactory] : Newton solve finish\n";

                U->add(1.,DeltaU); //*U += *DeltaU;

                M_appli->updateResidual( U, M_R );

                ///////////
                if (nDeltaUOld!=0) DeltaU->scale(1./nDeltaUOld);
                if (cpt>0) ErrDeltaU = DeltaU->l2Norm();
                //if (nDeltaUOld!=0) ErrDeltaU/=nDeltaUOld;
                Err=ErrDeltaU;



                std::cout<<"[ModelAlgebraicFactory] : Newton  : iter " << cpt << " Delta Norm L2 : " << Err << "\n";
#if 0
                auto nUUUold=Uold->l2Norm();
                Uold->add(-1.,U);
                if (nUUUold!=0) Uold->scale(1./nUUUold);
                if (cpt>0) ErrDiff=Uold->l2Norm();
                std::cout << "\n norm difff ="<< ErrDiff <<"\n";
#endif

                //*(M_appli->getSolution()) = *U;
                //M_appli->exportResults();

                ++cpt;

            }

        std::cout << "[ModelAlgebraicFactory] : Newton Algo finish\n";
#endif
    }




} // namespace FeelModels
} // namespace Feel
