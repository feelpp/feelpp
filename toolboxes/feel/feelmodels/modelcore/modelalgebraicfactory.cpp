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

#include <feel/feelmodels/modelcore/modelalgebraicfactory.hpp>
#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/enums.hpp>

namespace Feel
{
namespace FeelModels
{

ModelAlgebraicFactory::ModelAlgebraicFactory( std::string const& prefix, po::variables_map const& vm )
    :
    M_useSolverPtAP( false ),
    M_dofElimination_strategy( Feel::ContextOnMap[soption(_prefix=prefix,_name="on.type",_vm=vm)] ),
    M_dofElimination_valueOnDiagonal( doption(_prefix=prefix,_name="on.value_on_diagonal",_vm=vm) ),
    M_hasBuildLinearJacobian(false),
    M_hasBuildResidualCst(false),
    M_hasBuildLinearSystemCst(false),
    M_usePseudoTransientContinuation( boption(_prefix=prefix,_name="pseudo-transient-continuation",_vm=vm) ),
    M_pseudoTransientContinuationEvolutionMethod( soption(_prefix=prefix,_name="pseudo-transient-continuation.evolution",_vm=vm) ),
    M_pseudoTransientContinuationDelta0( doption(_prefix=prefix,_name="pseudo-transient-continuation.delta0",_vm=vm) ),
    M_pseudoTransientContinuationDeltaMax( doption(_prefix=prefix,_name="pseudo-transient-continuation.delta-max",_vm=vm) ),
    M_pseudoTransientContinuationSerVariant( soption(_prefix=prefix,_name="pseudo-transient-continuation.ser-variant",_vm=vm) ),
    M_pseudoTransientContinuationExpurThresholdHigh( doption(_prefix=prefix,_name="pseudo-transient-continuation.expur.threshold-high",_vm=vm) ),
    M_pseudoTransientContinuationExpurThresholdLow( doption(_prefix=prefix,_name="pseudo-transient-continuation.expur.threshold-low",_vm=vm) ),
    M_pseudoTransientContinuationExpurBetaHigh( doption(_prefix=prefix,_name="pseudo-transient-continuation.expur.beta-high",_vm=vm) ),
    M_pseudoTransientContinuationExpurBetaLow( doption(_prefix=prefix,_name="pseudo-transient-continuation.expur.beta-low",_vm=vm) )
{}

ModelAlgebraicFactory::ModelAlgebraicFactory( model_ptrtype const& model, backend_ptrtype const& backend )
    :
    ModelAlgebraicFactory( model->prefix(), model->clovm() )
{
    model->log( model->prefix()+".MethodNum","constructor1", "start" );
    this->init( model,backend );
    model->log( model->prefix()+".MethodNum","constructor1", "finish" );
}

ModelAlgebraicFactory::ModelAlgebraicFactory( model_ptrtype const& model, backend_ptrtype const& backend,
                                              graph_ptrtype const& graph, indexsplit_ptrtype const& indexSplit )
    :
    ModelAlgebraicFactory( model->prefix() )
{
    model->log( model->prefix()+".MethodNum","constructor2", "start" );
    this->init( model,backend,graph,indexSplit );
    model->log( model->prefix()+".MethodNum","constructor2", "finish" );
}

//---------------------------------------------------------------------------------------------------------------//

void
ModelAlgebraicFactory::init( model_ptrtype const& model, backend_ptrtype const& backend )
{
    model->timerTool("Constructor").start();
    auto graph = model->buildMatrixGraph();
    model->timerTool("Constructor").elapsed("graph");
    auto indexSplit = ( graph )? graph->mapRow().indexSplit() : indexsplit_ptrtype();
    this->init( model,backend,graph,graph->mapRow().indexSplit() );
    model->timerTool("Constructor").stop();
}
void
ModelAlgebraicFactory::init( model_ptrtype const& model, backend_ptrtype const& backend,
                             graph_ptrtype const& graph, indexsplit_ptrtype const& indexSplit )
{
    M_model = model;
    this->setFunctionLinearAssembly( std::bind( &model_type::updateLinearPDE,
                                                std::ref( *model ), std::placeholders::_1 ) );
    this->setFunctionLinearDofElimination( std::bind( &model_type::updateLinearPDEDofElimination,
                                                      std::ref( *model ), std::placeholders::_1 ) );
    this->setFunctionNewtonInitialGuess( std::bind( &model_type::updateNewtonInitialGuess,
                                                    std::ref( *model ), std::placeholders::_1 ) );
    this->setFunctionJacobianAssembly( std::bind( &model_type::updateJacobian,
                                                  std::ref( *model ), std::placeholders::_1 ) );
    this->setFunctionResidualAssembly( std::bind( &model_type::updateResidual,
                                                  std::ref( *model ), std::placeholders::_1 ) );

    this->setFunctionUpdateInHousePreconditionerLinear( std::bind( static_cast<void(model_type::*)(ModelAlgebraic::DataUpdateLinear&) const>(&model_type::updateInHousePreconditioner),
                                                                   std::ref( *model ), std::placeholders::_1 ) );
    this->setFunctionUpdateInHousePreconditionerJacobian( std::bind( static_cast<void(model_type::*)(ModelAlgebraic::DataUpdateJacobian&) const>(&model_type::updateInHousePreconditioner),
                                                                     std::ref( *model ), std::placeholders::_1 ) );

    this->init( backend, graph, indexSplit );
}
void
ModelAlgebraicFactory::init( backend_ptrtype const& backend, graph_ptrtype const& graph, indexsplit_ptrtype const& indexSplit )
{
    M_backend = backend;

    if ( graph )
    {
        this->model()->timerTool("Constructor").start();
        this->buildMatrixVector(graph,indexSplit);
        this->model()->timerTool("Constructor").elapsed("matrixVector");

        this->model()->timerTool("Constructor").restart();
        this->buildOthers();
        this->model()->timerTool("Constructor").stop("algebraicOthers");
    }
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
        // matrix used for compute the preconditioner (the same)
        M_Prec = M_J;

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

void ModelAlgebraicFactory::initExplictPartOfSolution()
{
    if ( !M_explictPartOfSolution )
        M_explictPartOfSolution = M_backend->newVector( M_R->mapPtr() );
}

void ModelAlgebraicFactory::initSolverPtAP( sparse_matrix_ptrtype matP, sparse_matrix_ptrtype matQ )
    {
        CHECK ( matP ) << "invalid matP";
        // TODO CHECK compatibility with A

        M_useSolverPtAP = true;
        M_solverPtAP_matP = matP;
        M_solverPtAP_matQ = matQ;

        // if already built we assume that the new matP as the same stencil
        if ( !M_solverPtAP_backend )
        {
            M_solverPtAP_backend = backend_type::build( soption( _name="backend" ), this->model()->prefix(), this->model()->worldCommPtr() );

            M_solverPtAP_matPtAP = M_solverPtAP_backend->newZeroMatrix( M_solverPtAP_matP->mapColPtr(), M_solverPtAP_matP->mapColPtr() );
            M_solverPtAP_matPtAP->clear();
            M_J->PtAP( *M_solverPtAP_matP, *M_solverPtAP_matPtAP );

            M_solverPtAP_matPtAP->setIndexSplit( M_J->indexSplit() );


            M_solverPtAP_solution = M_solverPtAP_backend->newVector( M_solverPtAP_matPtAP->mapColPtr() );
            M_solverPtAP_PtF = M_solverPtAP_backend->newVector( M_solverPtAP_matP->mapColPtr() );

            M_solverPtAP_Psolution = M_solverPtAP_backend->newVector( M_solverPtAP_matP->mapRowPtr() );

            M_solverPtAP_prec = preconditioner(_pc=(PreconditionerType) M_solverPtAP_backend->pcEnumType() /*LU_PRECOND*/,
                                               _matrix=M_solverPtAP_matPtAP,
                                               _backend= M_solverPtAP_backend,
                                               _pcfactormatsolverpackage=(MatSolverPackageType) M_solverPtAP_backend->matSolverPackageEnumType(),
                                               _worldcomm=M_solverPtAP_backend->comm(),
                                               _prefix=M_solverPtAP_backend->prefix() );

            // necessary to use the nonlinear solver
            M_solverPtAP_backend->nlSolver()->jacobian = std::bind( &self_type::updateJacobian,
                                                                    std::ref( *this ), std::placeholders::_1, std::placeholders::_2 );
            M_solverPtAP_backend->nlSolver()->residual = std::bind( &self_type::updateResidual,
                                                                    std::ref( *this ), std::placeholders::_1, std::placeholders::_2 );
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
    ModelAlgebraicFactory::reset( backend_ptrtype backend,
                                  graph_ptrtype const& graph,
                                  indexsplit_ptrtype const& indexSplit)
    {
        M_hasBuildLinearJacobian=false;
        M_hasBuildResidualCst=false;
        M_hasBuildLinearSystemCst=false;

        this->init( backend, graph, indexSplit );
    }

    void
    ModelAlgebraicFactory::attachNullSpace( NullSpace<value_type> const& nullSpace )
    {
        CHECK( this->backend() ) << "backend not init\n";
        std::shared_ptr<NullSpace<value_type> > mynullspace( new NullSpace<value_type>(this->backend(),nullSpace) );
        this->backend()->attachNullSpace( mynullspace );
        if ( M_solverPtAP_backend )
            M_solverPtAP_backend->attachNullSpace( mynullspace );
    }
    void
    ModelAlgebraicFactory::attachNearNullSpace( NullSpace<value_type> const& nearNullSpace )
    {
        CHECK( this->backend() ) << "backend not init\n";
        std::shared_ptr<NullSpace<value_type> > myNearNullSpace( new NullSpace<value_type>(this->backend(),nearNullSpace) );
        this->backend()->attachNearNullSpace( myNearNullSpace );
        if ( M_solverPtAP_backend )
            M_solverPtAP_backend->attachNearNullSpace( myNearNullSpace );
    }
    void
    ModelAlgebraicFactory::attachNearNullSpace( int k, NullSpace<value_type> const& nearNullSpace, std::string const& prefix )
    {
        CHECK( this->backend() ) << "backend not init\n";
        CHECK( M_PrecondManage ) << "preconditioner not init\n";
        std::shared_ptr<NullSpace<value_type> > myNearNullSpace( new NullSpace<value_type>(this->backend(),nearNullSpace) );
        M_PrecondManage->attachNearNullSpace( k, myNearNullSpace, prefix );
        if ( M_solverPtAP_prec )
            M_solverPtAP_prec->attachNearNullSpace( k, myNearNullSpace, prefix );
    }
    void
    ModelAlgebraicFactory::attachNearNullSpace( int k, NullSpace<value_type> const& nearNullSpace )
    {
        this->attachNearNullSpace( k, nearNullSpace, this->preconditionerTool()->name() );
    }


    void
    ModelAlgebraicFactory::attachAuxiliarySparseMatrix( std::string const& key,sparse_matrix_ptrtype const& mat )
    {
        M_PrecondManage->attachAuxiliarySparseMatrix( key, mat );
        if ( M_solverPtAP_prec )
            M_solverPtAP_prec->attachAuxiliarySparseMatrix( key, mat );
    }

    bool
    ModelAlgebraicFactory::hasAuxiliarySparseMatrix( std::string const& key ) const
    {
        return M_PrecondManage->hasAuxiliarySparseMatrix( key );
    }

    sparse_matrix_ptrtype const&
    ModelAlgebraicFactory::auxiliarySparseMatrix( std::string const& key ) const
    {
        return M_PrecondManage->auxiliarySparseMatrix( key );
    }

    void
    ModelAlgebraicFactory::attachOperatorPCD( std::string const& key, typename preconditioner_type::operator_pcdbase_ptrtype const& opPCD )
    {
        M_PrecondManage->attachOperatorPCD( key, opPCD );
        if ( M_solverPtAP_prec )
            M_solverPtAP_prec->attachOperatorPCD( key, opPCD );
    }


    void
    ModelAlgebraicFactory::addFunctionLinearAssembly( function_assembly_linear_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionLinearAssembly.size()).str() : key;
        M_addFunctionLinearAssembly[ keyUsed ] = func;
    }
    void
    ModelAlgebraicFactory::addFunctionLinearDofElimination( function_assembly_linear_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionLinearDofElimination.size()).str() : key;
        M_addFunctionLinearDofElimination[ keyUsed ] = func;
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
    ModelAlgebraicFactory::addFunctionJacobianDofElimination( function_assembly_jacobian_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionJacobianDofElimination.size()).str() : key;
        M_addFunctionJacobianDofElimination[ keyUsed ] = func;
    }
    void
    ModelAlgebraicFactory::addFunctionResidualDofElimination( function_assembly_residual_type const& func, std::string const& key )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addFunctionResidualDofElimination.size()).str() : key;
        M_addFunctionResidualDofElimination[ keyUsed ] = func;
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
        M_addVectorLinearRhsAssembly[ keyUsed ] = std::make_tuple( vec, scaling, cstPart, true );
    }
    void
    ModelAlgebraicFactory::addVectorResidualAssembly( vector_ptrtype const& vec, double scaling, std::string const& key, bool cstPart )
    {
        std::string keyUsed = ( key.empty() )? (boost::format("FEELPP_DEFAULT_%1%")%M_addVectorResidualAssembly.size()).str() : key;
        M_addVectorResidualAssembly[ keyUsed ] = std::make_tuple( vec, scaling, cstPart, true );
    }


    void
    ModelAlgebraicFactory::setActivationAddVectorLinearRhsAssembly( std::string const& key, bool b )
    {
        if ( M_addVectorLinearRhsAssembly.find( key ) == M_addVectorLinearRhsAssembly.end() )
            return;
        std::get<3>( M_addVectorLinearRhsAssembly[key] ) = b;
    }
    void
    ModelAlgebraicFactory::setActivationAddVectorResidualAssembly( std::string const& key, bool b )
    {
        if ( M_addVectorResidualAssembly.find( key ) == M_addVectorResidualAssembly.end() )
            return;
        std::get<3>( M_addVectorResidualAssembly[key] ) = b;

    }

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::updateInformationObject( nl::json & p ) const
    {
        if ( p.contains( "Backend" ) )
            return;
        nl::json subPt;
        subPt.emplace( "prefix",this->backend()->prefix() );
        subPt.emplace( "type", backend_type::enumToKind( this->backend()->type() ) /*"BACKEND_PETSC"*/ );
        p["Backend"] = subPt;

        // KSP
        subPt.clear();
        subPt.emplace( "type", this->backend()->kspType() );
        subPt.emplace( "rtol", this->backend()->rTolerance() );
        subPt.emplace( "dtol", this->backend()->dTolerance() );
        subPt.emplace( "atol", this->backend()->aTolerance() );
        subPt.emplace( "maxit", this->backend()->maxIterationsKSP() );
        subPt.emplace( "reuse-prec", this->backend()->reusePrec() );
        if ( this->backend()->reusePrec() )
            subPt.emplace( "maxit reuse-prec", this->backend()->maxIterationsKSPReuse() );
        p["KSP"] = subPt;

        //SNES
        subPt.clear();
        subPt.emplace( "rtol",this->backend()->rToleranceSNES() );
        subPt.emplace( "stol",this->backend()->sToleranceSNES() );
        subPt.emplace( "atol",this->backend()->aToleranceSNES() );
        subPt.emplace( "maxit", this->backend()->maxIterationsSNES() );
        subPt.emplace( "reuse-jac", this->backend()->reuseJac() );
        if ( this->backend()->reuseJac() )
        {
            subPt.emplace( "reuse-jac rebuild at first Newton step", this->backend()->reuseJacRebuildAtFirstNewtonStep() );
            subPt.emplace( "maxit with reuse", this->backend()->maxIterationsSNESReuse() );
        }
        p["SNES"] = subPt;

        // KSP in SNES
        subPt.clear();
        subPt.emplace( "rtol",this->backend()->rtoleranceKSPinSNES() );
        subPt.emplace( "maxit",this->backend()->maxIterationsKSPinSNES() );
        subPt.emplace( "reuse-prec", this->backend()->reusePrec() );
        if ( this->backend()->reusePrec() )
        {
            subPt.emplace( "reuse-prec maxit", this->backend()->maxIterationsKSPinSNESReuse() );
            subPt.emplace( "reuse-prec rebuild at first Newton step", this->backend()->reusePrecRebuildAtFirstNewtonStep() );
        }
        p["KSP in SNES"] = subPt;

        // PC
        subPt.clear();
        subPt.emplace( "type", this->backend()->pcType() );
        if ( this->backend()->pcType() == "lu" )
            subPt.emplace( "mat-solver-package", this->backend()->pcFactorMatSolverPackageType() );
        p["PC"] = subPt;
    }

tabulate_informations_ptr_t
ModelAlgebraicFactory::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New();

    for ( std::string const& section : std::vector<std::string>({"Backend","KSP","SNES","KSP in SNES", "PC" }) )
        {
            if ( !jsonInfo.contains( section ) )
                continue;
            Feel::Table tabInfoSectionEntries;
            TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoSectionEntries, jsonInfo.at( section ), tabInfoProp );
            tabInfoSectionEntries.format()
                .setShowAllBorders( false )
                .setColumnSeparator(":")
                .setHasRowSeparator( false );
            tabInfo->add( section, TabulateInformations::New( tabInfoSectionEntries ) );
        }

    return tabInfo;
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
            ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true);
            dataLinearCst.copyInfos( this->dataInfos() );
            M_functionLinearAssembly( dataLinearCst );
            for ( auto const& func : M_addFunctionLinearAssembly )
                func.second( dataLinearCst );
            M_CstR->close();
            for ( auto const& av : M_addVectorLinearRhsAssembly )
                if ( std::get<2>( av.second ) && std::get<3>( av.second ) )
                    M_CstR->add( std::get<1>( av.second ), std::get<0>( av.second ) );
            M_hasBuildLinearSystemCst = true;
        }
        else if ( this->model()->rebuildCstPartInLinearSystem() || this->model()->needToRebuildCstPart() ||
                  !this->model()->useCstMatrix() || !this->model()->useCstVector() )
        {
            M_CstJ->zero();
            M_CstR->zero();
            ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true);
            dataLinearCst.copyInfos( this->dataInfos() );
            M_functionLinearAssembly( dataLinearCst );
            for ( auto const& func : M_addFunctionLinearAssembly )
                func.second( dataLinearCst );
            M_CstR->close();
            for ( auto const& av : M_addVectorLinearRhsAssembly )
                if ( std::get<2>( av.second ) && std::get<3>( av.second ) )
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

        // assembling non cst part
        ModelAlgebraic::DataUpdateLinear dataLinearNonCst(U,M_J,M_R,false);
        dataLinearNonCst.copyInfos( this->dataInfos() );
        M_functionLinearAssembly( dataLinearNonCst );
        for ( auto const& func : M_addFunctionLinearAssembly )
            func.second( dataLinearNonCst );

        M_R->close();
        // add maybe vector to rhs
        for ( auto const& av : M_addVectorLinearRhsAssembly )
            if ( !std::get<2>( av.second ) && std::get<3>( av.second ) )
                M_R->add( std::get<1>( av.second ), std::get<0>( av.second ) );

        if ( M_explictPartOfSolution )
        {
            M_explictPartOfSolution->scale( -1. );
            M_R->addVector(*M_explictPartOfSolution, *M_J );
            M_explictPartOfSolution->scale( -1. );
        }

        pre_solve_type pre_solve = std::bind(&self_type::preSolveLinear, std::ref( *this ), std::placeholders::_1, std::placeholders::_2);
        post_solve_type post_solve = std::bind(&self_type::postSolveLinear, std::ref( *this ), std::placeholders::_1, std::placeholders::_2);

        if ( M_useSolverPtAP )
        {
            M_solverPtAP_matQ->multVector( *U, *M_solverPtAP_solution );
            //*M_solverPtAP_solution = *U; // need to have the same size/datamap

            // PtAP = P^T A P
            M_J->PtAP( *M_solverPtAP_matP, *M_solverPtAP_matPtAP );
            // PtF = P^T F
            M_solverPtAP_matP->multVector( *M_R, *M_solverPtAP_PtF, true );

            ModelAlgebraic::DataUpdateLinear dataLinearPtAP(M_solverPtAP_solution,M_solverPtAP_matPtAP,M_solverPtAP_PtF,false);
            // dof elimination
            M_functionLinearDofElimination( dataLinearPtAP );
            for ( auto const& func : M_addFunctionLinearDofElimination )
                func.second( dataLinearPtAP );
            // others dof eliminations (not need an expression,just for fix maybe the well posed system)
            if ( M_solverPtAP_dofEliminationIds )
            {
                // we assume that all shared dofs are present, not need to appy a sync
                std::vector<int> _dofs;_dofs.assign( M_solverPtAP_dofEliminationIds->begin(), M_solverPtAP_dofEliminationIds->end());
                auto tmp = M_solverPtAP_backend->newVector( U->mapPtr() );
                M_solverPtAP_matPtAP->zeroRows( _dofs, *tmp, *M_solverPtAP_PtF, M_dofElimination_strategy, M_dofElimination_valueOnDiagonal );
            }

            // post-assembly
            for ( auto const& func : M_addFunctionLinearPostAssembly )
                func.second( dataLinearPtAP );
        }
        else
        {
            // dof elimination
            M_functionLinearDofElimination( dataLinearNonCst );
            for ( auto const& func : M_addFunctionLinearDofElimination )
                func.second( dataLinearNonCst );

            // post-assembly
            for ( auto const& func : M_addFunctionLinearPostAssembly )
                func.second( dataLinearNonCst );
        }

        double mpiTimerAssembly = this->model()->timerTool("Solve").elapsed("algebraic-assembly");
        if (this->model()->verboseSolverTimer())
            Feel::FeelModels::Log(this->model()->prefix()+".ModelAlgebraicFactory","linearSolver",
                                  (boost::format("finish assembling in %1% s") % mpiTimerAssembly ).str(),
                                  this->model()->worldComm(),this->model()->verboseSolverTimerAllProc());
        this->model()->timerTool("Solve").restart();

        //-------------------------------------------------//
        // solve linear system
        //-------------------------------------------------//
        typename backend_type::solve_return_type solveStat;
        if ( M_useSolverPtAP )
        {
            ModelAlgebraic::DataUpdateLinear dataLinearSolverPtAP(U,M_solverPtAP_matPtAP,M_solverPtAP_PtF,false);

            std::invoke( M_functionUpdateInHousePreconditionerLinear, dataLinearSolverPtAP );

            solveStat = M_solverPtAP_backend->solve( _matrix=M_solverPtAP_matPtAP,
                                                     _solution=M_solverPtAP_solution,
                                                     _rhs=M_solverPtAP_PtF,//M_R,
                                                     _prec=M_solverPtAP_prec,
                                                     _pre=pre_solve,
                                                     _post=post_solve );

            M_solverPtAP_matP->multVector( M_solverPtAP_solution, U );
        }
        else
        {
            // set preconditioner
            //M_PrecondManage->setMatrix(M_Prec);
            std::invoke( M_functionUpdateInHousePreconditionerLinear, dataLinearNonCst );

            // solve linear system
            solveStat = M_backend->solve( _matrix=M_J,
                                          _solution=U,
                                          _rhs=M_R,
                                          _prec=M_PrecondManage,
                                          _pre=pre_solve,
                                          _post=post_solve );
        }

        if ( M_explictPartOfSolution )
        {
            U->add( 1.0, M_explictPartOfSolution );
        }

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
    ModelAlgebraicFactory::updateJacobian(const vector_ptrtype& XX, sparse_matrix_ptrtype& JJ/*, vector_ptrtype& R*/ )
    {
        auto model = this->model();
        model->timerTool("Solve").start();

        auto R = this->backend()->newVector(XX->mapPtr()); // TO REMOVE

        vector_ptrtype currentSolution = XX;
        sparse_matrix_ptrtype currentJacobian = JJ;
        if ( M_useSolverPtAP )
        {
            currentSolution = M_solverPtAP_Psolution;
            currentJacobian = M_J;
            M_solverPtAP_matP->multVector( XX, currentSolution );
        }

        currentJacobian->zero();

        if ( model->useCstMatrix())
        {
            currentJacobian->addMatrix(1.0, M_CstJ );
        }
        else
        {
            ModelAlgebraic::DataUpdateJacobian dataJacobianCst(currentSolution, currentJacobian, R, true);
            dataJacobianCst.copyInfos( this->dataInfos() );
            M_functionJacobianAssembly( dataJacobianCst );
            for ( auto const& func : M_addFunctionJacobianAssembly )
                func.second( dataJacobianCst );
        }

        ModelAlgebraic::DataUpdateJacobian dataJacobianNonCst(currentSolution, currentJacobian,R,false);
        dataJacobianNonCst.copyInfos( this->dataInfos() );
        if ( M_usePseudoTransientContinuation )
        {
            dataJacobianNonCst.addInfo( "use-pseudo-transient-continuation" );
            CHECK( !M_pseudoTransientContinuationDeltaAndResidual.empty() ) << "must have at least one value";
            dataJacobianNonCst.addDoubleInfo( "pseudo-transient-continuation.delta", M_pseudoTransientContinuationDeltaAndResidual.back().first );
        }

        M_functionJacobianAssembly( dataJacobianNonCst );
        for ( auto const& func : M_addFunctionJacobianAssembly )
            func.second( dataJacobianNonCst );

        currentJacobian->close();
        if ( M_useSolverPtAP )
        {
            // PtAP = P^T A P
            M_J->PtAP( *M_solverPtAP_matP, *JJ );
        }

        ModelAlgebraic::DataUpdateJacobian dataJacobianDofElimination(XX,JJ,R,false);
        dataJacobianDofElimination.copyInfos( this->dataInfos() );

        // dof elimination
        model->updateJacobianDofElimination( dataJacobianDofElimination );
        for ( auto const& func : M_addFunctionJacobianDofElimination )
            func.second( dataJacobianDofElimination );

        if ( dataJacobianDofElimination.hasDofEliminationIds() )
        {
            // we assume that all shared dofs are present, not need to appy a sync
            std::vector<int> _dofs;_dofs.assign( dataJacobianDofElimination.dofEliminationIds().begin(), dataJacobianDofElimination.dofEliminationIds().end());
            auto tmp = dataJacobianDofElimination.vectorUsedInStrongDirichlet();
            JJ->zeroRows( _dofs, *tmp, *tmp, M_dofElimination_strategy, M_dofElimination_valueOnDiagonal );
        }

        if ( M_useSolverPtAP && M_solverPtAP_dofEliminationIds )
        {
            // we assume that all shared dofs are present, not need to appy a sync
            std::vector<int> _dofs;_dofs.assign( M_solverPtAP_dofEliminationIds->begin(), M_solverPtAP_dofEliminationIds->end());
            auto tmp = dataJacobianDofElimination.vectorUsedInStrongDirichlet();
            JJ->zeroRows( _dofs, *tmp, *tmp, M_dofElimination_strategy, M_dofElimination_valueOnDiagonal );
        }

        for ( auto const& func : M_addFunctionJacobianPostAssembly )
            func.second( dataJacobianDofElimination );

        std::invoke( M_functionUpdateInHousePreconditionerJacobian, dataJacobianDofElimination );

        double tElapsed = model->timerTool("Solve").stop();
        model->timerTool("Solve").addDataValue("algebraic-jacobian",tElapsed);
    }
    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::updateResidual(const vector_ptrtype& XX, vector_ptrtype& RR)
    {
        auto model = this->model();
        model->timerTool("Solve").start();

        vector_ptrtype currentSolution = XX;
        vector_ptrtype currentResidual = RR;
        if ( M_useSolverPtAP )
        {
            currentSolution = M_solverPtAP_Psolution;
            currentResidual = M_R;
            M_solverPtAP_matP->multVector( XX, currentSolution );
        }

        currentResidual->zero();

        if ( model->useCstVector())
        {
            currentResidual->add(1.0, M_CstR );
        }
        else
        {
            ModelAlgebraic::DataUpdateResidual dataResidualCst( currentSolution, currentResidual, true, true );
            dataResidualCst.copyInfos( this->dataInfos() );
            M_functionResidualAssembly( dataResidualCst );
            for ( auto const& func : M_addFunctionResidualAssembly )
                func.second( dataResidualCst );
            currentResidual->close();
            for ( auto const& av : M_addVectorResidualAssembly )
                if ( std::get<2>( av.second ) && std::get<3>( av.second ) )
                    currentResidual->add( std::get<1>( av.second ), std::get<0>( av.second ) );
            if ( M_explictPartOfSolution )
                currentResidual->add(1.0,M_contributionsExplictPartOfSolutionWithNewton);
        }

        bool doOptimization = this->model()->useLinearJacobianInResidual() && this->model()->useCstMatrix();
        // add linear contribution from jacobian terms
        if (doOptimization)
            currentResidual->addVector(*currentSolution, *M_CstJ );

        ModelAlgebraic::DataUpdateResidual dataResidualNonCst( currentSolution, currentResidual, false, doOptimization );
        dataResidualNonCst.copyInfos( this->dataInfos() );
        if ( M_explictPartOfSolution )
            dataResidualNonCst.addVectorInfo( "explicit-part-of-solution", M_explictPartOfSolution );
        M_functionResidualAssembly( dataResidualNonCst );
        for ( auto const& func : M_addFunctionResidualAssembly )
            func.second( dataResidualNonCst );

        currentResidual->close();

        for ( auto const& av : M_addVectorResidualAssembly )
            if ( !std::get<2>( av.second ) && std::get<3>( av.second ) )
                currentResidual->add( std::get<1>( av.second ), std::get<0>( av.second ) );


        if ( M_useSolverPtAP )
        {
            // PtF = P^T F
            M_solverPtAP_matP->multVector( *M_R, *RR, true );
        }

        ModelAlgebraic::DataUpdateResidual dataResidualDofElimination( XX,RR, false, doOptimization );
        dataResidualDofElimination.copyInfos( this->dataInfos() );
        // dof elimination
        model->updateResidualDofElimination( dataResidualDofElimination );
        for ( auto const& func : M_addFunctionResidualDofElimination )
            func.second( dataResidualDofElimination );

        if ( dataResidualDofElimination.hasDofEliminationIds() )
        {
            for ( size_type k : dataResidualDofElimination.dofEliminationIds() )
                RR->set( k, 0. );
            // we assume that all shared dofs are present, not need to appy a sync
        }
        if ( M_useSolverPtAP && M_solverPtAP_dofEliminationIds )
        {
            for ( size_type k : *M_solverPtAP_dofEliminationIds )
                RR->set( k, 0. );
        }

        for ( auto const& func : M_addFunctionResidualPostAssembly )
            func.second( dataResidualDofElimination );

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
        ModelAlgebraic::DataNewtonInitialGuess dataInitialGuess( U );
        M_functionNewtonInitialGuess( dataInitialGuess );
        for ( auto const& func : M_addFunctionNewtonInitialGuess )
            func.second( dataInitialGuess );

        if ( dataInitialGuess.hasDofEliminationIds() )
        {
            // create view in order to avoid mpi comm inside petsc
            auto UView = VectorUblas<value_type>::createView( *U );
            // imposed values by ordering the entities
            std::vector<ElementsType> fromEntities = { MESH_ELEMENTS, MESH_FACES, MESH_EDGES, MESH_POINTS };
            for ( ElementsType entity : fromEntities )
            {
                if ( dataInitialGuess.hasDofEliminationIds( entity ) )
                    sync( UView, "=", dataInitialGuess.dofEliminationIds( entity ) );
            }
        }

        model->timerTool("Solve").elapsed("algebraic-newton-initial-guess");
        model->timerTool("Solve").restart();

#if 0
        if ( M_explictPartOfSolution )
            M_solverPtAP_solution->add( -1.0, M_explictPartOfSolution );
#endif

        //---------------------------------------------------------------------//
        if (model->useCstMatrix())
        {
            if (!M_hasBuildLinearJacobian || model->rebuildLinearPartInJacobian() || model->needToRebuildCstPart() )
            {
                M_CstJ->zero();
                ModelAlgebraic::DataUpdateJacobian dataJacobianCst(U, M_CstJ, M_R, true );
                dataJacobianCst.copyInfos( this->dataInfos() );
                if ( M_explictPartOfSolution )
                    dataJacobianCst.addVectorInfo( "explicit-part-of-solution", M_explictPartOfSolution );
                M_functionJacobianAssembly( dataJacobianCst );
                for ( auto const& func : M_addFunctionJacobianAssembly )
                    func.second( dataJacobianCst );
                M_CstJ->close();
                M_hasBuildLinearJacobian = true;
            }
        }

        //---------------------------------------------------------------------//
        model->timerTool("Solve").elapsed("algebraic-jacobian");
        model->timerTool("Solve").restart();
        //---------------------------------------------------------------------//
        if (model->useCstVector())
        {
            // model->rebuildCstPartInResidual() is not an option (just optimisation with fsi semi-implicit) TODO REMOVE THIS OPTION!!!
            if (!M_hasBuildResidualCst || model->rebuildCstPartInResidual() || model->needToRebuildCstPart())
            {
                M_CstR->zero();
                // Warning : the second true is very important in order to build M_CstR!!!!!!
                ModelAlgebraic::DataUpdateResidual dataResidualCst( U, M_CstR, true, true );
                dataResidualCst.copyInfos( this->dataInfos() );
                M_functionResidualAssembly( dataResidualCst );
                for ( auto const& func : M_addFunctionResidualAssembly )
                    func.second( dataResidualCst );
                M_CstR->close();
                for ( auto const& av : M_addVectorResidualAssembly )
                    if ( std::get<2>( av.second ) && std::get<3>( av.second ) )
                        M_CstR->add( std::get<1>( av.second ), std::get<0>( av.second ) );
                M_hasBuildResidualCst = true;
            }
        }
        model->setNeedToRebuildCstPart(false);
        //---------------------------------------------------------------------//
        model->timerTool("Solve").elapsed("algebraic-residual");
        model->timerTool("Solve").restart();
        //---------------------------------------------------------------------//

        if ( M_explictPartOfSolution )
        {
            std::vector<std::string> infoAdded = { "ignore-assembly.rhs" };
            if (model->useCstVector())
            {
                this->evaluateResidual( M_explictPartOfSolution, M_CstR, infoAdded, false );
                M_CstR->close();
            }
            else
            {
                if ( !M_contributionsExplictPartOfSolutionWithNewton )
                    M_contributionsExplictPartOfSolutionWithNewton = M_backend->newVector( M_R->mapPtr() );
                M_contributionsExplictPartOfSolutionWithNewton->zero();
                this->evaluateResidual( M_explictPartOfSolution, M_contributionsExplictPartOfSolutionWithNewton, infoAdded, false );
                M_contributionsExplictPartOfSolutionWithNewton->close();
            }
        }

        pre_solve_type pre_solve = std::bind(&self_type::preSolveNewton, std::ref( *this ), std::placeholders::_1, std::placeholders::_2);
        post_solve_type post_solve = std::bind(&self_type::postSolveNewton, std::ref( *this ), std::placeholders::_1, std::placeholders::_2);
        update_nlsolve_type update_nlsolve = std::bind(&self_type::updateNewtonIteration, std::ref( *this ), std::placeholders::_1, std::placeholders::_2,std::placeholders::_3,std::placeholders::_4 );

        typename backend_type::nl_solve_return_type solveStat;
        if ( M_useSolverPtAP )
        {
            // TODO REMOVE EXPLICIT PART
            M_solverPtAP_matQ->multVector( *U, *M_solverPtAP_solution );
            //*M_solverPtAP_solution = *U; // need to have the same size/datamap

            solveStat = M_solverPtAP_backend->nlSolve( _jacobian=M_solverPtAP_matPtAP,
                                                       _solution=M_solverPtAP_solution,
                                                       _residual=M_solverPtAP_PtF,
                                                       _prec=M_solverPtAP_prec,
                                                       _pre=pre_solve,
                                                       _post=post_solve,
                                                       _update=update_nlsolve );

            M_solverPtAP_matP->multVector( M_solverPtAP_solution, U );
        }
        else
        {
            solveStat = M_backend->nlSolve( _jacobian=M_J,
                                            _solution=U,
                                            _residual=M_R,
                                            _prec=M_PrecondManage,
                                            _pre=pre_solve,
                                            _post=post_solve,
                                            _update=update_nlsolve );
        }

        if ( M_explictPartOfSolution )
            U->add( 1.0, M_explictPartOfSolution );

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
    ModelAlgebraicFactory::applyAssemblyLinear(const vector_ptrtype& U, sparse_matrix_ptrtype& lhs, vector_ptrtype& rhs, std::vector<std::string> const& infos, bool applyDofElimination ) const
    {
        ModelAlgebraic::DataUpdateLinear dataLinear(U,lhs,rhs,true);
        dataLinear.copyInfos( this->dataInfos() );
        for ( std::string const& info : infos )
            dataLinear.addInfo( info );
        this->applyAssemblyLinear( dataLinear, applyDofElimination );
    }
    void
    ModelAlgebraicFactory::applyAssemblyLinear( ModelAlgebraic::DataUpdateLinear & dataLinear, bool applyDofElimination ) const
    {
        // cst part
        dataLinear.setBuildCstPart( true );
        M_functionLinearAssembly( dataLinear );
        for ( auto const& func : M_addFunctionLinearAssembly )
            func.second( dataLinear );

        // non cst part
        dataLinear.setBuildCstPart( false );
        M_functionLinearAssembly( dataLinear );
        for ( auto const& func : M_addFunctionLinearAssembly )
            func.second( dataLinear );

        // add maybe vector to rhs
        for ( auto const& av : M_addVectorLinearRhsAssembly )
            if ( /*std::get<2>( av.second ) &&*/ std::get<3>( av.second ) )
                dataLinear.rhs()->add( std::get<1>( av.second ), std::get<0>( av.second ) );

        // dof elimination
        if ( applyDofElimination )
        {
            M_functionLinearDofElimination( dataLinear );
            for ( auto const& func : M_addFunctionLinearDofElimination )
                func.second( dataLinear );
        }

        // post-assembly
        for ( auto const& func : M_addFunctionLinearPostAssembly )
            func.second( dataLinear );
    }

    void
    ModelAlgebraicFactory::evaluateResidual(const vector_ptrtype& U, vector_ptrtype& R, std::vector<std::string> const& infos, bool applyDofElimination ) const
    {
        ModelAlgebraic::DataUpdateResidual dataResidual( U, R, true, false );
        for ( std::string const& info : infos )
            dataResidual.addInfo( info );
        this->evaluateResidual( dataResidual, applyDofElimination );
    }

    void
    ModelAlgebraicFactory::evaluateResidual( ModelAlgebraic::DataUpdateResidual & dataResidual, bool applyDofElimination ) const
    {
        auto model = this->model();
        dataResidual.setBuildCstPart( true );
        dataResidual.setUseJacobianLinearTerms( false );
        vector_ptrtype& R = dataResidual.residual();
        //R->zero();
        M_functionResidualAssembly( dataResidual );
        for ( auto const& func : M_addFunctionResidualAssembly )
            func.second( dataResidual );
        dataResidual.setBuildCstPart( false );
        M_functionResidualAssembly( dataResidual );
        for ( auto const& func : M_addFunctionResidualAssembly )
            func.second( dataResidual );
        R->close();
        for ( auto const& av : M_addVectorResidualAssembly )
            if ( std::get<3>( av.second ) )
                R->add( std::get<1>( av.second ), std::get<0>( av.second ) );

        // dof elimination
        if ( applyDofElimination )
        {
            model->updateResidualDofElimination( dataResidual );
            for ( auto const& func : M_addFunctionResidualDofElimination )
                func.second( dataResidual );
            if ( dataResidual.hasDofEliminationIds() )
            {
                for ( size_type k : dataResidual.dofEliminationIds() )
                    R->set( k, 0. );
            }
        }

        for ( auto const& func : M_addFunctionResidualPostAssembly )
            func.second( dataResidual );
    }

    //---------------------------------------------------------------------------------------------------------------//

    void
    ModelAlgebraicFactory::rebuildCstJacobian( vector_ptrtype U )
    {
        M_CstJ->zero();
        ModelAlgebraic::DataUpdateJacobian dataJacobianCst(U, M_CstJ, M_R, true );
        M_functionJacobianAssembly( dataJacobianCst );
        for ( auto const& func : M_addFunctionJacobianAssembly )
            func.second( dataJacobianCst );
        M_CstJ->close();
    }


    void
    ModelAlgebraicFactory::rebuildCstLinearPDE( vector_ptrtype U )
    {
        M_CstJ->zero();
        M_CstR->zero();
        ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true);
        M_functionLinearAssembly(dataLinearCst);
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

        auto dataLinear = std::make_shared<ModelAlgebraic::DataUpdateLinear>(U,M_J,M_R,false);
        dataLinear->copyInfos( this->dataInfos() );
        std::shared_ptr<ModelAlgebraic::DataUpdateLinear> dataLinearPtAP;
        if ( M_useSolverPtAP )
            dataLinearPtAP = std::make_shared<ModelAlgebraic::DataUpdateLinear>(M_solverPtAP_solution,M_solverPtAP_matPtAP,M_solverPtAP_PtF,true);

        this->updatePicardIteration( 0, U );

        // assembling cst part
        if ( !M_hasBuildLinearSystemCst ||
             this->model()->rebuildCstPartInLinearSystem() || this->model()->needToRebuildCstPart() ||
             !this->model()->useCstMatrix() || !this->model()->useCstVector() )
        {
            this->model()->timerTool("Solve").start();

            M_CstJ->zero();
            M_CstR->zero();
            ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_CstJ,M_CstR,true);
            dataLinearCst.copyInfos( this->dataInfos() );
            M_functionLinearAssembly( dataLinearCst );
            M_hasBuildLinearSystemCst = true;

            double tAssemblyElapsed = this->model()->timerTool("Solve").stop();
            this->model()->timerTool("Solve").addDataValue("algebraic-assembly",tAssemblyElapsed);
        }
        this->model()->setNeedToRebuildCstPart(false);

        vector_ptrtype Uold = M_backend->newVector( M_R->mapPtr() );
        *Uold = *U;
        vector_ptrtype residual;
        if ( useConvergenceAlgebraic )
            residual = M_backend->newVector( M_R->mapPtr() );

        bool hasConverged = false;

        double convergenceRate=1;
        double rtol = this->backend()->rToleranceSNES();//1e-6;
        double atol = this->backend()->aToleranceSNES();//1e-50;
        double stol = this->backend()->sToleranceSNES();
        double initialResiduValue = 0;
        int fixPointMaxIt = this->backend()->maxIterationsSNES();
        int cptIteration=0;
        for ( ; cptIteration < fixPointMaxIt ; ++cptIteration )
        {
            if ( cptIteration > 0 )
                this->updatePicardIteration( cptIteration, U );

            if ( !useConvergenceAlgebraic )
            {
                *Uold = *U;
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
                dataLinear->setBuildCstPart( true );
                //ModelAlgebraic::DataUpdateLinear dataLinearCst(U,M_J,M_R,true);
                M_functionLinearAssembly( *dataLinear );
            }

            // assembling non cst part
            dataLinear->setBuildCstPart( false );
            //ModelAlgebraic::DataUpdateLinear dataLinearNonCst(U,M_J,M_R,false);
            M_functionLinearAssembly( *dataLinear );

            if ( M_explictPartOfSolution )
            {
                M_explictPartOfSolution->scale( -1. );
                M_R->addVector(*M_explictPartOfSolution, *M_J );
                M_explictPartOfSolution->scale( -1. );
            }

            std::shared_ptr<ModelAlgebraic::DataUpdateLinear> dataDofElimination = dataLinear;
            if ( M_useSolverPtAP )
            {
                *M_solverPtAP_solution = *U;
                // PtAP = P^T A P
                M_J->PtAP( *M_solverPtAP_matP, *M_solverPtAP_matPtAP );
                // PtF = P^T F
                M_solverPtAP_matP->multVector( *M_R, *M_solverPtAP_PtF, true );

                dataDofElimination = dataLinearPtAP;
            }

            // dof elimination
            M_functionLinearDofElimination( *dataDofElimination );
            for ( auto const& func : M_addFunctionLinearDofElimination )
                func.second( *dataDofElimination );

            // others dof eliminations (not need an expression)
            if ( M_useSolverPtAP && M_solverPtAP_dofEliminationIds )
            {
                // we assume that all shared dofs are present, not need to appy a sync
                std::vector<int> _dofs;_dofs.assign( M_solverPtAP_dofEliminationIds->begin(), M_solverPtAP_dofEliminationIds->end());
                auto tmp = M_solverPtAP_backend->newVector( U->mapPtr() );
                M_solverPtAP_matPtAP->zeroRows( _dofs, *tmp, *tmp, M_dofElimination_strategy, M_dofElimination_valueOnDiagonal );
            }

            // post-assembly
            for ( auto const& func : M_addFunctionLinearPostAssembly )
                func.second( *dataDofElimination );


            double tAssemblyElapsed = this->model()->timerTool("Solve").elapsed();
            this->model()->timerTool("Solve").addDataValue("algebraic-assembly",tAssemblyElapsed);

            if (this->model()->verboseSolverTimer())
                Feel::FeelModels::Log(this->model()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                      (boost::format("picard iteration[%1%] finish assembling in %2% s") %cptIteration % tAssemblyElapsed ).str(),
                                      this->model()->worldComm(),this->model()->verboseSolverTimerAllProc());

            this->model()->timerTool("Solve").restart();

            auto themat = dataDofElimination->matrix();
            auto therhs = dataDofElimination->rhs();
            auto thesol = dataDofElimination->currentSolution();

            if ( useConvergenceAlgebraic )
            {
                themat->close();
                therhs->close();
                residual->zero();
                residual->addVector( thesol, themat );//U,M_J );
                residual->add( -1., therhs /*M_R*/ );
                if ( M_useSolverPtAP && M_solverPtAP_dofEliminationIds )
                {
                    for ( size_type k : *M_solverPtAP_dofEliminationIds )
                        residual->set( k, 0. );
                    residual->close();
                }
                //residual->close();
                double nomResidual = residual->l2Norm();
                //double normRhs = M_R->l2Norm();
                if (this->model()->verboseSolverTimer())
                    Feel::FeelModels::Log( this->model()->prefix()+".ModelAlgebraicFactory","AlgoPicard",
                                           (boost::format("picard iteration[%1%] residual norm : %2%")%cptIteration %nomResidual ).str(),
                                           this->model()->worldComm(),this->model()->verboseSolverTimerAllProc() );

                if ( cptIteration == 0 )
                    initialResiduValue = nomResidual;
                else
                {
                    //if ( nomResidual < std::max(rtol*normRhs,atol) )
                    if ( nomResidual < std::max(rtol*initialResiduValue,atol) || nomResidual < atol )
                    {
                        hasConverged = true;
                        this->model()->timerTool("Solve").stop();
                        break;
                    }
                }
            }

            // update in-house preconditioners
            std::invoke( M_functionUpdateInHousePreconditionerLinear, *dataDofElimination );

            // set pre/post solve functions
            pre_solve_type pre_solve = std::bind(&model_type::preSolvePicard, this->model(), std::placeholders::_1, std::placeholders::_2);
            post_solve_type post_solve = std::bind(&model_type::postSolvePicard, this->model(), std::placeholders::_1, std::placeholders::_2);

            backend_ptrtype thebackend = M_useSolverPtAP? M_solverPtAP_backend : M_backend;
            preconditioner_ptrtype theprec = M_useSolverPtAP? M_solverPtAP_prec : M_PrecondManage;

            // solve linear system
            auto const solveStat = thebackend->solve( _matrix=themat,
                                                      _solution=thesol,
                                                      _rhs=therhs,
                                                      _prec=theprec,
                                                      _pre=pre_solve,
                                                      _post=post_solve );

            if ( M_useSolverPtAP )
            {
                M_solverPtAP_matP->multVector( thesol, U );
            }

            if ( M_explictPartOfSolution )
            {
                U->add( 1.0, M_explictPartOfSolution );
            }

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


            if ( useConvergenceAlgebraic )
            {
                Uold->add( -1., U );
                double normDelta = Uold->l2Norm();
                double normU = U->l2Norm();
                if ( normDelta < std::max(stol*normU,atol) )
                {
                    hasConverged = true;
                    break;
                }
            }
            *Uold = *U;
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

void
ModelAlgebraicFactory::updatePicardIteration( int step, vector_ptrtype sol )
{
    this->model()->updatePicardIteration( step, sol );
}


void
ModelAlgebraicFactory::preSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    this->model()->preSolveNewton( rhs, sol );
}
void
ModelAlgebraicFactory::postSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    this->model()->postSolveNewton( rhs, sol );
}
void
ModelAlgebraicFactory::preSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    this->model()->preSolvePicard( rhs, sol );
}
void
ModelAlgebraicFactory::postSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    this->model()->postSolvePicard( rhs, sol );
}
void
ModelAlgebraicFactory::preSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    this->model()->preSolveLinear( rhs, sol );
}
void
ModelAlgebraicFactory::postSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    this->model()->postSolveLinear( rhs, sol );
}


} // namespace FeelModels
} // namespace Feel
