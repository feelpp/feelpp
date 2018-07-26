/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2018-03-06

  Copyright (C) 2018 Feel++ Consortium

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

#include <feel/feelmodels/heatfluid/heatfluid.hpp>

#include <feel/feelvf/vf.hpp>
/*#include <feel/feelvf/form.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/operators.hpp>
 #include <feel/feelvf/operations.hpp>*/

#include <feel/feelmodels/modelmesh/createmesh.hpp>

namespace Feel
{
namespace FeelModels
{

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
HEATFLUID_CLASS_TEMPLATE_TYPE::HeatFluid( std::string const& prefix,
                                          bool buildMesh,
                                          WorldComm const& worldComm,
                                          std::string const& subPrefix,
                                          ModelBaseRepository const& modelRep )
    :
    super_type( prefix, worldComm, subPrefix, modelRep )
{
    this->log("HeatFluid","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatFluidConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatFluidSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatFluidPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatFluidTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // build mesh
    if ( buildMesh )
        this->initMesh();
    //-----------------------------------------------------------------------------//
    this->log("HeatFluid","constructor", "finish");
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_useNaturalConvection = boption(_prefix=this->prefix(),_name="use-natural-convection");//"Forced-Convection";

    M_BoussinesqRefTemperature = doption(_name="Boussinesq.ref-temperature",_prefix=this->prefix());
    std::string gravityStr;
    if ( Environment::vm().count(prefixvm(this->prefix(),"gravity-force").c_str()) )
        gravityStr = soption(_name="gravity-force",_prefix=this->prefix());
    else if (nDim == 2 )
        gravityStr = "{0,-9.80665}";
    else if (nDim == 3 )
        gravityStr = "{0,0,-9.80665}";
    M_gravityForce = expr<nDim,1,2>( gravityStr,"",this->worldComm(),this->repository().expr() );

    M_useSemiImplicitTimeScheme = boption(_name="use-semi-implicit-time-scheme",_prefix=this->prefix());
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("HeatFluid","initMesh", "start");
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("initMesh");
    this->log("HeatFluid","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()



HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
int
HEATFLUID_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = M_heatModel->nBlockMatrixGraph() + M_fluidModel->nBlockMatrixGraph();
    return nBlock;
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
HEATFLUID_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);

    int nBlockFluid = M_fluidModel->nBlockMatrixGraph();
    int nBlockHeat = M_heatModel->nBlockMatrixGraph();

    int startIndexBlockFluid = 0;
    int startIndexBlockHeat = nBlockFluid;

    auto blockMatFluid = M_fluidModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockFluid ;++tk1 )
        for (int tk2=0;tk2<nBlockFluid ;++tk2 )
            myblockGraph(startIndexBlockFluid+tk1,startIndexBlockFluid+tk2) = blockMatFluid(tk1,tk2);

    BlocksStencilPattern patCoupling1(1,M_fluidModel->functionSpace()->nSpaces,size_type(Pattern::ZERO));
    patCoupling1(0,0) = size_type(Pattern::COUPLED);
    myblockGraph(startIndexBlockHeat,startIndexBlockFluid) = stencil(_test=M_heatModel->spaceTemperature(),
                                                                     _trial=M_fluidModel->functionSpace(),
                                                                     _pattern_block=patCoupling1,
                                                                     _diag_is_nonzero=false,_close=false)->graph();

    BlocksStencilPattern patCoupling2(M_fluidModel->functionSpace()->nSpaces,1,size_type(Pattern::ZERO));
    patCoupling2(0,0) = size_type(Pattern::COUPLED);
    if ( M_fluidModel->stabilizationGLS() )
        patCoupling2(1,0) = size_type(Pattern::COUPLED);
    myblockGraph(startIndexBlockFluid,startIndexBlockHeat) = stencil(_test=M_fluidModel->functionSpace(),
                                                                     _trial=M_heatModel->spaceTemperature(),
                                                                     _pattern_block=patCoupling2,
                                                                     _diag_is_nonzero=false,_close=false)->graph();

    auto blockMatHeat = M_heatModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockHeat ;++tk1 )
        for (int tk2=0;tk2<nBlockHeat ;++tk2 )
            myblockGraph(startIndexBlockHeat+tk1,startIndexBlockHeat+tk2) = blockMatHeat(tk1,tk2);

    myblockGraph.close();

    return myblockGraph;
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("HeatFluid","init", "start" );
    this->timerTool("Constructor").start();

    M_heatModel.reset( new heat_model_type(prefixvm(this->prefix(),"heat"), false, this->worldComm(),
                                                           this->subPrefix(), this->repository() ) );
    if ( !M_heatModel->modelPropertiesPtr() )
        M_heatModel->setModelProperties( this->modelPropertiesPtr() );
    M_heatModel->setMesh( this->mesh() );
    M_heatModel->init( false );

    M_fluidModel.reset( new fluid_model_type(prefixvm(this->prefix(),"fluid"), false, this->worldComm(),
                                             this->subPrefix(), this->repository() ) );
    if ( !M_fluidModel->modelPropertiesPtr() )
        M_fluidModel->setModelProperties( this->modelPropertiesPtr() );
    M_fluidModel->setMesh( this->mesh() );
    M_fluidModel->init( false );

    if ( !this->isStationary() )
    {
        this->updateTime( this->timeStepBase()->time() );
        this->setTimeInitial( this->timeStepBase()->timeInitial() );
    }
    else
        M_useSemiImplicitTimeScheme = false;

    if ( !M_useNaturalConvection )
    {
        M_heatModel->initAlgebraicFactory();
        M_fluidModel->initAlgebraicFactory();
    }
    else
    {
        M_fluidModel->setStabilizationGLSDoAssembly( false );
        if ( M_useSemiImplicitTimeScheme )
            M_fluidModel->setSolverName("Oseen");
    }

    for ( auto const& rangeData : M_fluidModel->materialProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        if ( !M_heatModel->thermalProperties()->hasMaterial( matName ) )
            continue;
        M_rangeMeshElementsByMaterial[matName] = rangeData.second;
    }

    // post-process
    this->initPostProcess();

    // backend
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldComm() );

    // block vector solution
    auto blockVectorSolutionFluid = M_fluidModel->blockVectorSolution();
    auto blockVectorSolutionHeat = M_heatModel->blockVectorSolution();
    int nBlockFluid = blockVectorSolutionFluid.size();
    int nBlockHeat = blockVectorSolutionHeat.size();
    int nBlock = nBlockFluid + nBlockHeat;
    M_blockVectorSolution.resize( nBlock );
    int indexBlock=0;
    for ( int k=0;k<nBlockFluid ;++k )
        M_blockVectorSolution(indexBlock+k) = blockVectorSolutionFluid(k);
    indexBlock += nBlockFluid;
    for ( int k=0;k<nBlockHeat ;++k )
        M_blockVectorSolution(indexBlock+k) = blockVectorSolutionHeat(k);
    indexBlock += nBlockHeat;
    // init monolithic vector associated to the block vector
    M_blockVectorSolution.buildVector( this->backend() );

    size_type currentStartBlockSpaceIndex = 0;
    this->setStartSubBlockSpaceIndex( "fluid", currentStartBlockSpaceIndex );
    currentStartBlockSpaceIndex += blockVectorSolutionFluid.vectorMonolithic()->map().numberOfDofIdToContainerId();
    this->setStartSubBlockSpaceIndex( "heat", currentStartBlockSpaceIndex );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
    {
        if ( M_useNaturalConvection /*M_solverName == "Newton"*/ )
        {
            M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );
        }
    }

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("HeatFluid","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("HeatFluid","initPostProcess", "start");
    this->timerTool("Constructor").start();

    std::string modelName = "heat-fluid";
    auto const& exportsFields = this->modelProperties().postProcess().exports( modelName ).fields();
    M_postProcessFieldExportedHeatt = M_heatModel->postProcessFieldExported( exportsFields, "heat" );
    M_postProcessFieldExportedFluid = M_fluidModel->postProcessFieldExported( exportsFields, "fluid" );

    if ( !M_postProcessFieldExportedHeatt.empty() || !M_postProcessFieldExportedFluid.empty() )
    {
        std::string geoExportType="static";//change_coords_only, change, static
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export",
                               _geo=geoExportType,
                               _path=this->exporterPath() );

        if ( this->doRestart() && this->restartPath().empty() )
        {
            if ( M_exporter->doExport() )
                M_exporter->restart(this->timeInitial());
        }
    }

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("HeatFluid","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}



HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
HEATFLUID_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << M_heatModel->getInfo()->str();
    *_ostr << M_fluidModel->getInfo()->str();

    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||---------------Info : HeatFluid---------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository();
    *_ostr << "\n   Physical Model"
           << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient" )
           << "\n     -- natural convection  : " << std::string( (M_useNaturalConvection)? "ON":"OFF" );
    //*_ostr << "\n   Numerical Solver"
    //       << "\n     -- solver   : " << M_solverName;
    if ( M_exporter )
    {
        *_ostr << "\n   Exporter"
               << "\n     -- type            : " << M_exporter->type()
               << "\n     -- freq save       : " << M_exporter->freq();
        std::string fieldExportedHeat;
        for ( std::string const& fieldName : M_postProcessFieldExportedHeatt )
            fieldExportedHeat=(fieldExportedHeat.empty())? fieldName : fieldExportedHeat + " - " + fieldName;
        std::string fieldExportedFluid;
        for ( std::string const& fieldName : M_postProcessFieldExportedFluid )
            fieldExportedFluid=(fieldExportedFluid.empty())? fieldName : fieldExportedFluid + " - " + fieldName;
        *_ostr << "\n     -- fields [heat] : " << fieldExportedHeat
               << "\n     -- fields [fluid] : " << fieldExportedFluid;
    }

    *_ostr << "\n   Processors"
           << "\n     -- number of proc : " << this->worldComm().globalSize()
           << "\n     -- current rank : " << this->worldComm().globalRank();

    if ( M_algebraicFactory )
        *_ostr << M_algebraicFactory->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->heatModel()->updateTimeStep();
    this->fluidModel()->updateTimeStep();
    this->updateTime( this->timeStepBase()->time() );
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("HeatFluid","exportResults", "start");
    this->timerTool("PostProcessing").start();

    if ( M_exporter && M_exporter->doExport() )
    {
        bool hasFieldToExportHeat = M_heatModel->updateExportedFields( M_exporter,M_postProcessFieldExportedHeatt,time );
        bool hasFieldToExportFluid = M_fluidModel->updateExportedFields( M_exporter,M_postProcessFieldExportedFluid,time );
        if ( hasFieldToExportHeat || hasFieldToExportFluid )
        {
            this->upload( M_exporter->path() );
            M_exporter->save();
        }
    }

    M_heatModel->exportResults( time );
    M_fluidModel->exportResults( time );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("HeatFluid","exportResults", "finish");
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    M_heatModel->updateParameterValues();
    M_fluidModel->updateParameterValues();
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("HeatFluid","solve", "start");
    this->timerTool("Solve").start();

    if ( !M_useNaturalConvection )
    {
        M_fluidModel->solve();
        M_heatModel->setFieldVelocityConvectionIsUsed( true );
        for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
        {
            auto const& range = rangeData.second;
            M_heatModel->updateFieldVelocityConvection( range,idv(M_fluidModel->fieldVelocity()) );
        }
        M_heatModel->solve();
    }
    else
    {
        this->updateParameterValues();

        M_fluidModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("fluid") );
        M_heatModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("heat") );

        M_blockVectorSolution.updateVectorFromSubVectors();

        if ( M_useSemiImplicitTimeScheme )
            M_algebraicFactory->solve( "Linear", M_blockVectorSolution.vectorMonolithic() );
        else
            M_algebraicFactory->solve( "Newton", M_blockVectorSolution.vectorMonolithic() );

        M_blockVectorSolution.localize();
    }

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("HeatFluid","solve", (boost::format("finish in %1% s")%tElapsed).str() );

}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::postSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    M_heatModel->postSolveNewton( rhs, sol );
    M_fluidModel->postSolveNewton( rhs, sol );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::postSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    M_heatModel->postSolvePicard( rhs, sol );
    M_fluidModel->postSolvePicard( rhs, sol );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::postSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    M_heatModel->postSolveLinear( rhs, sol );
    M_fluidModel->postSolveLinear( rhs, sol );
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("HeatFluid","updateLinearPDE", "start"+sc);

    M_heatModel->updateLinearPDE( data );
    M_fluidModel->updateLinearPDE( data );

    if ( buildNonCstPart )
    {
        auto XhVP = M_fluidModel->spaceVelocityPressure();
        auto const& U = M_fluidModel->fieldVelocityPressure();
        auto u = U.template element<0>();
        auto XhT = M_heatModel->spaceTemperature();
        auto const& t = M_heatModel->fieldTemperature();

        auto mylfVP = form1( _test=XhVP, _vector=F,
                             _rowstart=M_fluidModel->rowStartInVector()+0 );

        auto mybfTT = form2( _test=XhT,_trial=XhT,_matrix=A,
                             _rowstart=M_heatModel->rowStartInMatrix(),
                             _colstart=M_heatModel->colStartInMatrix() );
        auto mybfVPT = form2( _test=XhVP,_trial=XhT,_matrix=A,
                             _rowstart=M_fluidModel->rowStartInMatrix(),
                             _colstart=M_heatModel->colStartInMatrix() );

        auto UConvection = M_fluidModel->timeStepBDF()->poly();
        auto uConvection = UConvection.template element<0>();

        for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            auto const& rhoHeatCapacity = M_heatModel->thermalProperties()->rhoHeatCapacity( matName );

            if ( rhoHeatCapacity.isConstant() )
            {
                double rhoHeatCapacityValue = rhoHeatCapacity.value();
                mybfTT +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityValue*(gradt(t)*idv(uConvection))*id(t),
                               _geomap=this->geomap() );
            }
            else
            {
                auto rhoHeatCapacityExpr = rhoHeatCapacity.expr();
                mybfTT +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityExpr*(gradt(t)*idv(uConvection))*id(t),
                               _geomap=this->geomap() );
            }
            auto const& rho = M_heatModel->thermalProperties()->rho( matName );
            auto const& thermalExpansion = M_heatModel->thermalProperties()->thermalExpansion( matName );
            CHECK( rhoHeatCapacity.isConstant() && thermalExpansion.isConstant() ) << "TODO";
            double rhoValue = rho.value();
            double beta = thermalExpansion.value();
            double T0 = M_BoussinesqRefTemperature;
            mybfVPT +=
                integrate( _range=range,
                           _expr= rhoValue*(beta*idt(t))*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );
            mylfVP +=
                integrate( _range=range,
                           _expr= rhoValue*(beta*T0)*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );

            if ( M_heatModel->stabilizationGLS() )
            {
                auto const& thermalConductivity = M_heatModel->thermalProperties()->thermalConductivity( matName );
                if ( thermalConductivity.isMatrix() )
                    CHECK( false ) << "TODO";
                else if ( thermalConductivity.isConstant() )
                    M_heatModel->updateLinearPDEStabilizationGLS( cst(rhoHeatCapacity.value()),cst(thermalConductivity.value()),idv(uConvection),range,data );
                else
                    CHECK( false ) << "TODO";
            }
            if ( M_fluidModel->stabilizationGLS() )
            {
                auto rhoF = idv(M_fluidModel->materialProperties()->fieldRho());
                //auto mu = Feel::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
                auto mu = idv(M_fluidModel->materialProperties()->fieldMu());
                auto expraddedInGLSResidualLF = rhoValue*beta*T0*M_gravityForce;
                auto exprAddedInGLSResidualBF = rhoValue*beta*idt(t)*M_gravityForce;
                M_fluidModel->updateLinearPDEStabilisationGLS( data, rhoF, mu, matName, hana::make_tuple(expraddedInGLSResidualLF),hana::make_tuple(std::make_pair(mybfVPT, exprAddedInGLSResidualBF)) );
            }

        }

    }

    this->log("HeatFluid","updateLinearPDE", "finish");

}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
    this->log("HeatFluid","updateNewtonInitialGuess","start" );
    M_heatModel->updateNewtonInitialGuess( U );
    M_fluidModel->updateNewtonInitialGuess( U );
    this->log("HeatFluid","updateNewtonInitialGuess","finish" );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("HeatFluid","updateJacobian", "start"+sc);

    auto mesh = this->mesh();

    M_heatModel->updateJacobian( data );
    M_fluidModel->updateJacobian( data );

    if ( buildNonCstPart )
    {
        auto XhVP = M_fluidModel->spaceVelocityPressure();
        auto const U = XhVP->element(XVec, M_fluidModel->rowStartInVector()+0 );
        auto u = U.template element<0>();

        auto XhT = M_heatModel->spaceTemperature();
        auto t = XhT->element(XVec, M_heatModel->rowStartInVector() );
        auto const& thermalProperties = M_heatModel->thermalProperties();

        auto bfVPT = form2( _test=XhVP,_trial=XhT,_matrix=J,
                            _rowstart=M_fluidModel->rowStartInMatrix(),
                            _colstart=M_heatModel->colStartInMatrix() );

        for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            auto const& rho = M_heatModel->thermalProperties()->rho( matName );
            auto const& rhoHeatCapacity = M_heatModel->thermalProperties()->rhoHeatCapacity( matName );
            auto const& thermalExpansion = M_heatModel->thermalProperties()->thermalExpansion( matName );
            CHECK( rho.isConstant() && rhoHeatCapacity.isConstant() && thermalExpansion.isConstant() ) << "TODO";

            double rhoHeatCapacityValue = rhoHeatCapacity.value();
            double rhoValue = rho.value();
            double beta = thermalExpansion.value();

            form2( _test=XhT,_trial=XhT,_matrix=J,
                   _rowstart=M_heatModel->rowStartInMatrix(),
                   _colstart=M_heatModel->colStartInMatrix() ) +=
                integrate( _range=range,
                           _expr= rhoHeatCapacityValue*(gradt(t)*idv(u))*id(t),
                           _geomap=this->geomap() );

            form2( _test=XhT,_trial=XhVP,_matrix=J,
                   _rowstart=M_heatModel->rowStartInMatrix(),
                   _colstart=M_fluidModel->colStartInMatrix() ) +=
                integrate( _range=range,
                           _expr= rhoHeatCapacityValue*(gradv(t)*idt(u))*id(t),
                           _geomap=this->geomap() );

            bfVPT +=
                integrate( _range=range,
                           _expr= rhoValue*beta*idt(t)*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );

            if ( M_heatModel->stabilizationGLS() )
            {
                auto const& thermalConductivity = M_heatModel->thermalProperties()->thermalConductivity( matName );
                if ( thermalConductivity.isMatrix() )
                    CHECK( false ) << "TODO";
                else if ( thermalConductivity.isConstant() )
                    M_heatModel->updateJacobianStabilizationGLS( cst(rhoHeatCapacity.value()),cst(thermalConductivity.value()),idv(u),range,data );
                else
                    CHECK( false ) << "TODO";
            }

            if ( M_fluidModel->stabilizationGLS() )
            {
                auto rhoF = idv(M_fluidModel->materialProperties()->fieldRho());
                //auto mu = Feel::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
                auto mu = idv(M_fluidModel->materialProperties()->fieldMu());
                auto exprAddedInGLSResidual = rhoValue*beta*idt(t)*M_gravityForce;
                M_fluidModel->updateJacobianStabilisationGLS( data, U, rhoF, mu, matName, std::make_pair(bfVPT, exprAddedInGLSResidual) );
            }

        }
    }

    this->log("HeatFluid","updateJacobian", "finish"+sc);

}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool useJacobianLinearTerms = data.useJacobianLinearTerms();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("HeatFluid","updateResidual", "start"+sc);

    auto mesh = this->mesh();

    M_heatModel->updateResidual( data );
    M_fluidModel->updateResidual( data );

    if ( buildNonCstPart )
    {
        auto XhVP = M_fluidModel->spaceVelocityPressure();
        auto const U = XhVP->element(XVec, M_fluidModel->rowStartInVector()+0 );
        auto u = U.template element<0>();

        auto XhT = M_heatModel->spaceTemperature();
        auto const t = XhT->element(XVec, M_heatModel->rowStartInVector()+0 );
        auto mylfT = form1( _test=XhT, _vector=R,
                            _rowstart=M_heatModel->rowStartInVector()+0 );
        auto mylfVP = form1( _test=XhVP, _vector=R,
                             _rowstart=M_fluidModel->rowStartInVector()+0 );


        for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            auto const& rhoHeatCapacity = M_heatModel->thermalProperties()->rhoHeatCapacity( matName );
            if ( rhoHeatCapacity.isConstant() )
            {
                double rhoHeatCapacityValue = rhoHeatCapacity.value();
                mylfT +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityValue*(gradv(t)*idv(u))*id(t),
                               _geomap=this->geomap() );
            }
            else
            {
                auto rhoHeatCapacityExpr = rhoHeatCapacity.expr();
                mylfT +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityExpr*(gradv(t)*idv(u))*id(t),
                               _geomap=this->geomap() );
            }
            auto const& rho = M_heatModel->thermalProperties()->rho( matName );
            auto const& thermalExpansion = M_heatModel->thermalProperties()->thermalExpansion( matName );
            CHECK( rhoHeatCapacity.isConstant() && thermalExpansion.isConstant() ) << "TODO";
            double rhoValue = rho.value();
            double beta = thermalExpansion.value();
            double T0 = M_BoussinesqRefTemperature;
            mylfVP +=
                integrate( _range=range,
                           _expr= rhoValue*(beta*(idv(t)-T0))*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );


            if ( M_heatModel->stabilizationGLS() )
            {
                auto const& thermalConductivity = M_heatModel->thermalProperties()->thermalConductivity( matName );
                if ( thermalConductivity.isMatrix() )
                    CHECK( false ) << "TODO";
                else if ( thermalConductivity.isConstant() )
                    M_heatModel->updateResidualStabilizationGLS( cst(rhoHeatCapacity.value()),cst(thermalConductivity.value()),idv(u),range,data );
                else
                    CHECK( false ) << "TODO";
            }

            if ( M_fluidModel->stabilizationGLS() )
            {
                auto rhoF = idv(M_fluidModel->materialProperties()->fieldRho());
                //auto mu = Feel::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
                auto mu = idv(M_fluidModel->materialProperties()->fieldMu());
                auto expraddedInGLSResidual = rhoValue*(beta*(idv(t)-T0))*M_gravityForce;
                M_fluidModel->updateResidualStabilisationGLS( data, U, rhoF, mu, matName, expraddedInGLSResidual );
            }

        }

    }

    this->log("HeatFluid","updateResidual", "finish"+sc);
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    M_heatModel->updateLinearPDEDofElimination( data );
    M_fluidModel->updateLinearPDEDofElimination( data );
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    M_heatModel->updateJacobianDofElimination( data );
    M_fluidModel->updateJacobianDofElimination( data );
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    M_heatModel->updateResidualDofElimination( data );
    M_fluidModel->updateResidualDofElimination( data );
}


} // end namespace FeelModels
} // end namespace Feel
