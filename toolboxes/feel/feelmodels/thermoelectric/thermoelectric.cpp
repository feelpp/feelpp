/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2016-12-12

  Copyright (C) 2016 Feel++ Consortium

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
   \file thermoelectric.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-12-12
 */

#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>

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

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::ThermoElectric( std::string const& prefix,
                                                    bool buildMesh,
                                                    worldcomm_ptr_t const& worldComm,
                                                    std::string const& subPrefix,
                                                    ModelBaseRepository const& modelRep )
    :
    super_type( prefix, worldComm, subPrefix, modelRep )
{
    this->log("ThermoElectric","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoElectricConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoElectricSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoElectricPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoElectricTimeStepping.data";
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
        this->createMesh();
    //-----------------------------------------------------------------------------//
    this->log("ThermoElectric","constructor", "finish");
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_solverNewtonInitialGuessUseLinearThermoElectric = boption(_prefix=this->prefix(),_name="solver-newton.initial-guess.use-linear-thermo-electric");
    M_solverNewtonInitialGuessUseLinearHeat = boption(_prefix=this->prefix(),_name="solver-newton.initial-guess.use-linear-heat");
    M_solverNewtonInitialGuessUseLinearElectric = boption(_prefix=this->prefix(),_name="solver-newton.initial-guess.use-linear-electric");
    if ( M_solverNewtonInitialGuessUseLinearThermoElectric )
    {
        M_solverNewtonInitialGuessUseLinearHeat = true;
        M_solverNewtonInitialGuessUseLinearElectric = true;
    }
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::createMesh()
{
    this->log("ThermoElectric","createMesh", "start");
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("createMesh");
    this->log("ThermoElectric","createMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = M_heatModel->nBlockMatrixGraph() + M_electricModel->nBlockMatrixGraph();
    return nBlock;
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);

    int indexBlock=0;

    int nBlockHeat = M_heatModel->nBlockMatrixGraph();
    auto blockMatHeat = M_heatModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockHeat ;++tk1 )
        for (int tk2=0;tk2<nBlockHeat ;++tk2 )
            myblockGraph(indexBlock+tk1,indexBlock+tk2) = blockMatHeat(tk1,tk2);

    BlocksStencilPattern patCoupling1(1,nBlockHeat,size_type(Pattern::ZERO));
    patCoupling1(0,0) = size_type(Pattern::COUPLED);
    myblockGraph(indexBlock,indexBlock+nBlockHeat) = stencil(_test=M_heatModel->spaceTemperature(),
                                                                     _trial=M_electricModel->spaceElectricPotential(),
                                                                     _pattern_block=patCoupling1,
                                                                     _diag_is_nonzero=false,_close=false)->graph();

    if ( true )
    {
        BlocksStencilPattern patCoupling2(nBlockHeat,1,size_type(Pattern::ZERO));
        patCoupling2(0,0) = size_type(Pattern::COUPLED);
        myblockGraph(indexBlock+nBlockHeat,indexBlock) = stencil(_test=M_electricModel->spaceElectricPotential(),
                                                                         _trial=M_heatModel->spaceTemperature(),
                                                                         _pattern_block=patCoupling2,
                                                                         _diag_is_nonzero=false,_close=false)->graph();
    }

    indexBlock += nBlockHeat;

    int nBlockElectric = M_electricModel->nBlockMatrixGraph();
    auto blockMatElectric = M_electricModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockElectric ;++tk1 )
        for (int tk2=0;tk2<nBlockElectric ;++tk2 )
            myblockGraph(indexBlock+tk1,indexBlock+tk2) = blockMatElectric(tk1,tk2);

    myblockGraph.close();

    return myblockGraph;
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("ThermoElectric","init", "start" );
    this->timerTool("Constructor").start();

#if 0
    M_modelName = this->modelProperties().models().model("thermo-electric").equations();
    if ( M_modelName.empty() )
        M_modelName = "ThermoElectric";

    M_solverName = "Linear";
    if ( M_modelName == "ThermoElectric-linear" )
    {
        M_modelName = "ThermoElectric";
        M_solverName = "Linear";
    }
    else if ( M_modelName == "ThermoElectric-nonlinear" )
    {
        M_modelName = "ThermoElectric";
        M_solverName = "Newton";
    }
    CHECK( M_modelName == "ThermoElectric" ) << "invalid model name : " << M_modelName << "\n";
#endif
 

    M_heatModel.reset( new heat_model_type(prefixvm(this->prefix(),"heat"), false, this->worldCommPtr(),
                                                           this->subPrefix(), this->repository() ) );
    if ( !M_heatModel->modelPropertiesPtr() )
        M_heatModel->setModelProperties( this->modelPropertiesPtr() );
    M_heatModel->setMesh( this->mesh() );
    M_heatModel->init( false );

    M_electricModel.reset( new electric_model_type(prefixvm(this->prefix(),"electric"), false, this->worldCommPtr(),
                                                   this->subPrefix(), this->repository() ) );
    if ( !M_electricModel->modelPropertiesPtr() )
        M_electricModel->setModelProperties( this->modelPropertiesPtr() );
    M_electricModel->setMesh( this->mesh() );
    M_electricModel->init( false );

    M_modelName = "ThermoElectric";
    M_solverName = "Linear";
    if ( M_electricModel->electricProperties()->hasElectricConductivityDependingOnSymbol( "heat_T" ) )
        M_solverName = "Newton";

    M_modelUseJouleEffect = true;

    if ( M_solverName == "Linear" || M_solverNewtonInitialGuessUseLinearHeat )
    {
        M_heatModel->initAlgebraicFactory();
        M_heatModel->algebraicFactory()->addFunctionLinearAssembly( boost::bind( &self_type::updateLinearPreAssemblyJouleLaw,
                                                                                 boost::ref( *this ), _1 ) );
        M_heatModel->algebraicFactory()->addFunctionResidualAssembly( boost::bind( &self_type::updateResidualPreAssemblyJouleLaw,
                                                                                   boost::ref( *this ), _1 ) );
    }
    if ( M_solverName == "Linear" || M_solverNewtonInitialGuessUseLinearElectric )
    {
        M_electricModel->initAlgebraicFactory();
        M_electricModel->algebraicFactory()->addFunctionLinearAssembly( boost::bind( &self_type::updateLinearElectricDependingOnTemperature,
                                                                                     boost::ref( *this ), _1 ) );
    }


#if 0
    M_rangeMeshElements = ( M_heatModel->thermalProperties()->isDefinedOnWholeMesh() && M_electricModel->electricProperties()->isDefinedOnWholeMesh() )?
        elements(this->mesh() ) :
        intersect( M_heatModel->rangeMeshElements(), M_electricModel->rangeMeshElements() );
#endif

    for ( auto const& rangeData : M_electricModel->electricProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        if ( !M_heatModel->thermalProperties()->hasMaterial( matName ) )
            continue;
        M_rangeMeshElementsByMaterial[matName] = rangeData.second;
    }
    // post-process
    this->initPostProcess();

    // backend
    M_backendMonolithic = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() );

    // block vector solution
    auto blockVectorSolutionHeat = M_heatModel->blockVectorSolution();
    auto blockVectorSolutionElectric = M_electricModel->blockVectorSolution();
    int nBlockHeat = blockVectorSolutionHeat.size();
    int nBlockElectric = blockVectorSolutionElectric.size();
    int nBlock = nBlockHeat + nBlockElectric;
    M_blockVectorSolutionMonolithic.resize( nBlock );
    int indexBlock=0;
    for ( int k=0;k<nBlockHeat ;++k )
        M_blockVectorSolutionMonolithic(indexBlock+k) = blockVectorSolutionHeat(k);
    indexBlock += nBlockHeat;
    for ( int k=0;k<nBlockElectric ;++k )
        M_blockVectorSolutionMonolithic(indexBlock+k) = blockVectorSolutionElectric(k);
    indexBlock += nBlockElectric;
    // init monolithic vector associated to the block vector
    M_blockVectorSolutionMonolithic.buildVector( this->backend() );

    size_type currentStartBlockSpaceIndex = 0;
    this->setStartSubBlockSpaceIndex( "heat", currentStartBlockSpaceIndex );
    currentStartBlockSpaceIndex += blockVectorSolutionHeat.vectorMonolithic()->map().numberOfDofIdToContainerId();
    this->setStartSubBlockSpaceIndex( "electric", currentStartBlockSpaceIndex );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
    {
        if ( M_solverName == "Newton" )
        {
            M_algebraicFactoryMonolithic.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );
        }
    }

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("ThermoElectric","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("ThermoElectric","initPostProcess", "start");
    this->timerTool("Constructor").start();

    std::string modelName = "thermo-electric";
    auto const& exportsFields = this->modelProperties().postProcess().exports( modelName ).fields();
    M_postProcessFieldExportedHeat = M_heatModel->postProcessFieldExported( exportsFields, "heat" );
    M_postProcessFieldExportedElectric = M_electricModel->postProcessFieldExported( exportsFields, "electric" );

    if ( !M_postProcessFieldExportedHeat.empty() || !M_postProcessFieldExportedElectric.empty() )
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
    this->log("ThermoElectric","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << M_heatModel->getInfo()->str();
    *_ostr << M_electricModel->getInfo()->str();

    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||----------Info : ThermoElectric---------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository();
    //*_ostr << "\n   Physical Model"
    //       << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient");
    *_ostr << "\n   Numerical Solver"
           << "\n     -- solver   : " << M_solverName;
    if ( M_exporter )
    {
        *_ostr << "\n   Exporter"
               << "\n     -- type            : " << M_exporter->type()
               << "\n     -- freq save       : " << M_exporter->freq();
        std::string fieldExportedHeat;
        for ( std::string const& fieldName : M_postProcessFieldExportedHeat )
            fieldExportedHeat=(fieldExportedHeat.empty())? fieldName : fieldExportedHeat + " - " + fieldName;
        std::string fieldExportedElectric;
        for ( std::string const& fieldName : M_postProcessFieldExportedElectric )
            fieldExportedElectric=(fieldExportedElectric.empty())? fieldName : fieldExportedElectric + " - " + fieldName;
        *_ostr << "\n     -- fields [heat] : " << fieldExportedHeat
               << "\n     -- fields [electric] : " << fieldExportedElectric;
    }
    *_ostr << "\n   Processors"
           << "\n     -- number of proc : " << this->worldComm().globalSize()
           << "\n     -- current rank : " << this->worldComm().globalRank();

    if ( M_algebraicFactoryMonolithic )
        *_ostr << M_algebraicFactoryMonolithic->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("ThermoElectric","exportResults", "start");
    this->timerTool("PostProcessing").start();

    if ( M_exporter && M_exporter->doExport() )
    {
        bool hasFieldToExportHeat = M_heatModel->updateExportedFields( M_exporter,M_postProcessFieldExportedHeat,time );
        bool hasFieldToExportElectric = M_electricModel->updateExportedFields( M_exporter,M_postProcessFieldExportedElectric,time );
        if ( hasFieldToExportHeat || hasFieldToExportElectric )
        {
            M_exporter->save();
            this->upload( M_exporter->path() );
        }
    }

    M_heatModel->exportResults( time );
    M_electricModel->exportResults( time );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("ThermoElectric","exportResults", "finish");
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    M_heatModel->updateParameterValues();
    M_electricModel->updateParameterValues();
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("ThermoElectric","solve", "start");
    this->timerTool("Solve").start();

    this->updateParameterValues();

    this->setStartBlockSpaceIndex( 0 );

    if ( M_solverName == "Linear" )
    {
        M_electricModel->solve();
        M_heatModel->solve();
    }
    else if ( M_solverName == "Newton" )
    {
        // initial guess
        if ( M_solverNewtonInitialGuessUseLinearElectric )
            M_electricModel->solve();
        if ( M_solverNewtonInitialGuessUseLinearHeat )
            M_heatModel->solve();

        // solve non linear monolithic system
        M_heatModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("heat") );
        M_electricModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("electric") );
        M_blockVectorSolutionMonolithic.updateVectorFromSubVectors();
        M_algebraicFactoryMonolithic->solve( "Newton", M_blockVectorSolutionMonolithic.vectorMonolithic() );
        M_blockVectorSolutionMonolithic.localize();

        M_electricModel->updateElectricField();
        this->updateCurrentDensity();
    }

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("ThermoElectric","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearElectricDependingOnTemperature( DataUpdateLinear & data ) const
{
    this->log("ThermoElectric","updateLinearElectricDependingOnTemperature","start" );

    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    if ( !buildNonCstPart )
        return;
    sparse_matrix_ptrtype& A = data.matrix();

    auto mesh = this->mesh();
    auto XhV = M_electricModel->spaceElectricPotential();
    auto const& v = M_electricModel->fieldElectricPotential();
    auto XhT = M_heatModel->spaceTemperature();
    auto const& t = M_heatModel->fieldTemperature();

    auto mybf = form2( _test=XhV,_trial=XhV,_matrix=A,
                       _pattern=size_type(Pattern::COUPLED),
                       _rowstart=this->rowStartInMatrix() ,
                       _colstart=this->colStartInMatrix() );

    for ( auto const& rangeData : M_rangeMeshElementsByMaterial )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& electricConductivity = M_electricModel->electricProperties()->electricConductivity( matName );
        if ( electricConductivity.isConstant() )
            continue;
        std::string symbolStr = "heat_T";
        if ( !electricConductivity.expr().expression().hasSymbol( symbolStr ) )
            continue;
        this->log("ThermoElectric","updateLinearElectricDependingOnTemperature","update material : "+matName );
        auto sigma = electricConductivity.expr( symbolStr, idv(t) );
        mybf +=
            integrate( _range=range,
                       _expr= sigma*inner(gradt(v),grad(v)),
                       _geomap=this->geomap() );
    }

    this->log("ThermoElectric","updateLinearElectricDependingOnTemperature","finish" );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearPreAssemblyJouleLaw( DataUpdateLinear & data ) const
{
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    if ( !buildNonCstPart )
        return;
    vector_ptrtype& F = data.rhs();
    this->updateGenericPreAssemblyJouleLaw( F, false );
}
THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateResidualPreAssemblyJouleLaw( DataUpdateResidual & data ) const
{
    bool buildCstPart = data.buildCstPart();
    if  ( !buildCstPart )
        return;
    vector_ptrtype& R = data.residual();
    this->updateGenericPreAssemblyJouleLaw( R, true );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateGenericPreAssemblyJouleLaw( vector_ptrtype& F, bool applyOnResidual ) const
{
    if ( !M_modelUseJouleEffect )
        return;
    this->log("ThermoElectric","updateGenericPreAssemblyJouleLaw","start" );

    auto mesh = this->mesh();
    auto XhV = M_electricModel->spaceElectricPotential();
    auto const& v = M_electricModel->fieldElectricPotential();
    auto XhT = M_heatModel->spaceTemperature();
    auto const& t = M_heatModel->fieldTemperature();

    auto myLinearForm = form1( _test=XhT,_vector=F,
                               _rowstart=M_heatModel->rowStartInVector() );

    int sign = ( applyOnResidual )? -1 : 1;
    for ( auto const& rangeData : M_rangeMeshElementsByMaterial )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& electricConductivity = M_electricModel->electricProperties()->electricConductivity( matName );
        if ( electricConductivity.isConstant() )
        {
            double sigma = sign*electricConductivity.value();
            myLinearForm +=
                integrate( _range=range,
                           _expr= sigma*inner(gradv(v)/*,gradv(v)*/)*id(t),
                           _geomap=this->geomap() );
        }
        else
        {
            std::string symbolStr = "heat_T";
            auto sigma = sign*electricConductivity.expr( symbolStr, idv(t) );
            //auto sigma = sign*electricConductivity.expr();
            //auto sigma = idv(M_electricModel->electricProperties()->fieldElectricConductivity());
            myLinearForm +=
                integrate( _range=range,
                           _expr= sigma*inner(gradv(v)/*,gradv(v)*/)*id(t),
                           _geomap=this->geomap() );
        }
    }

    this->log("ThermoElectric","updateGenericPreAssemblyJouleLaw","finish" );
}



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
    this->log("ThermoElectric","updateNewtonInitialGuess","start" );
    M_heatModel->updateNewtonInitialGuess( U );
    M_electricModel->updateNewtonInitialGuess( U );
    this->log("ThermoElectric","updateNewtonInitialGuess","finish" );
}
THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("ThermoElectric","updateJacobian", "start"+sc);
    size_type startBlockIndexTemperature = M_heatModel->startBlockSpaceIndexVector()+0;
    size_type startBlockIndexElectricPotential = M_electricModel->startBlockSpaceIndexVector()+0;

    auto mesh = this->mesh();

    if ( !buildCstPart )
    {
        auto XhV = M_electricModel->spaceElectricPotential();
        auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );
        auto XhT = M_heatModel->spaceTemperature();
        //auto const& t = M_heatModel->fieldTemperature();
        auto const t = XhT->element(XVec, this->rowStartInVector()+startBlockIndexTemperature );

        auto mybfTT = form2( _test=XhT,_trial=XhT,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexTemperature,
                             _colstart=this->colStartInMatrix()+startBlockIndexTemperature );
        auto mybfTV = form2( _test=XhT,_trial=XhV,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexTemperature,
                             _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );
        auto mybfVV = form2( _test=XhV,_trial=XhV,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                             _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );
        auto mybfVT = form2( _test=XhV,_trial=XhT,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                             _colstart=this->colStartInMatrix()+startBlockIndexTemperature );

        for ( auto const& rangeData : M_rangeMeshElementsByMaterial )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            auto const& electricConductivity = M_electricModel->electricProperties()->electricConductivity( matName );
            if ( electricConductivity.isConstant() )
            {
                if ( M_modelUseJouleEffect )
                {
                    double sigma = electricConductivity.value();
                    mybfTV +=
                        integrate( _range=range,
                                   _expr= -sigma*2*inner(gradt(v),gradv(v))*id( t ),
                                   _geomap=this->geomap() );
                }
            }
            else
            {
                std::string symbolStr = "heat_T";
                auto sigma = electricConductivity.expr( symbolStr, idv(t) );
                if ( M_modelUseJouleEffect )
                {
                    mybfTV +=
                        integrate( _range=range,
                                   _expr= -sigma*2*inner(gradt(v),gradv(v))*id( t ),
                                   _geomap=this->geomap() );
                }

                if ( sigma.expression().hasSymbol( symbolStr ) )
                {
                    auto sigmaDiff = diff( electricConductivity.expr(),symbolStr,1,"",this->worldComm(),this->repository().expr());
                    auto sigmaDiffEval = expr( sigmaDiff, symbolExpr( symbolStr, idv(t) ) );

                    if ( M_modelUseJouleEffect )
                    {
                        mybfTT +=
                            integrate( _range=range,
                                       _expr= -sigmaDiffEval*idt(t)*inner(gradv(v)/*,gradv(v)*/)*id( t ),
                                       _geomap=this->geomap() );
                    }

                    mybfVV +=
                        integrate( _range=range,
                                   _expr= sigma*inner(gradt(v),grad(v)),
                                   _geomap=this->geomap() );
                    mybfVT +=
                        integrate( _range=range,
                                   _expr= sigmaDiffEval*idt(t)*inner(gradv(v),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }
    }

    DataUpdateJacobian dataSubPhysics( data );
    M_heatModel->updateJacobian( dataSubPhysics );
    M_electricModel->updateJacobian( dataSubPhysics );

    this->log("ThermoElectric","updateJacobian", "finish"+sc);
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool useJacobianLinearTerms = data.useJacobianLinearTerms();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("ThermoElectric","updateResidual", "start"+sc);

    size_type startBlockIndexTemperature = M_heatModel->startBlockSpaceIndexVector()+0;
    size_type startBlockIndexElectricPotential = M_electricModel->startBlockSpaceIndexVector()+0;

    auto mesh = this->mesh();

    DataUpdateResidual dataSubPhysics( data );
    M_heatModel->updateResidual( dataSubPhysics );
    M_electricModel->updateResidual( dataSubPhysics );

    if ( !buildCstPart )
    {
        auto XhV = M_electricModel->spaceElectricPotential();
        auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );
        auto XhT = M_heatModel->spaceTemperature();
        auto const t = XhT->element(XVec, this->rowStartInVector()+startBlockIndexTemperature );
        auto mylfT = form1( _test=XhT, _vector=R,
                            _rowstart=this->rowStartInVector()+startBlockIndexTemperature );
        auto mylfV = form1( _test=XhV, _vector=R,
                            _rowstart=this->rowStartInVector()+startBlockIndexElectricPotential );

        for ( auto const& rangeData : M_rangeMeshElementsByMaterial )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            auto const& electricConductivity = M_electricModel->electricProperties()->electricConductivity( matName );
            if ( electricConductivity.isConstant() )
            {
                if ( M_modelUseJouleEffect )
                {
                    double sigma = electricConductivity.value();
                    mylfT +=
                        integrate( _range=range,
                                   _expr= -sigma*inner(gradv(v)/*,gradv(v)*/)*id( t ),
                                   _geomap=this->geomap() );
                }
            }
            else
            {
                std::string symbolStr = "heat_T";
                auto sigma = electricConductivity.expr( symbolStr, idv(t) );
                if ( M_modelUseJouleEffect )
                {
                    mylfT +=
                        integrate( _range=range,
                                   _expr= -sigma*inner(gradv(v)/*,gradv(v)*/)*id( t ),
                                   _geomap=this->geomap() );
                }
                if ( sigma.expression().hasSymbol( symbolStr ) )
                {
                    mylfV +=
                        integrate( _range=range,
                                   _expr= sigma*inner(gradv(v),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }
    }

    this->log("ThermoElectric","updateResidual", "finish"+sc);
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    M_heatModel->updateJacobianDofElimination( data );
    M_electricModel->updateJacobianDofElimination( data );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    M_heatModel->updateResidualDofElimination( data );
    M_electricModel->updateResidualDofElimination( data );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void  
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateCurrentDensity()
{
    auto const& v = M_electricModel->fieldElectricPotential();
    auto const& t = M_heatModel->fieldTemperature();

    for ( auto const& rangeData : M_rangeMeshElementsByMaterial )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& electricConductivity = M_electricModel->electricProperties()->electricConductivity( matName );
        if ( electricConductivity.isConstant() )
        {
            double sigma = electricConductivity.value();
            auto cd = -sigma*trans(gradv(v));
            M_electricModel->updateCurrentDensity( cd, range );
        }
        else
        {
            auto sigma = electricConductivity.expr();
            std::string symbolStr = "heat_T";
            if ( sigma.expression().hasSymbol( symbolStr ) )
            {
                auto sigma = electricConductivity.expr( symbolStr, idv(t) );
                auto cd = -sigma*trans(gradv(v));
                M_electricModel->updateCurrentDensity( cd, range );
            }
            else
            {
                auto cd = -sigma*trans(gradv(v));
                M_electricModel->updateCurrentDensity( cd, range );
            }
        }
    }
}

} // end namespace FeelModels
} // end namespace Feel
