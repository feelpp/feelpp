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
#include <feel/feelmodels/modelcore/utils.hpp>

namespace Feel
{
namespace FeelModels
{

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::ThermoElectric( std::string const& prefix,
                                                    std::string const& keyword,
                                                    worldcomm_ptr_t const& worldComm,
                                                    std::string const& subPrefix,
                                                    ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelPhysics<mesh_type::nDim>( "thermo-electric" ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep )
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
    this->log("ThermoElectric","constructor", "finish");
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_solverName = soption(_prefix=this->prefix(),_name="solver");
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
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("ThermoElectric","initMesh", "start");
    this->timerTool("Constructor").start();

    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( this->mesh() ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("createMesh");
    this->log("ThermoElectric","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

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

    M_heatModel.reset( new heat_model_type(prefixvm(this->prefix(),"heat"), "heat", this->worldCommPtr(),
                                           this->subPrefix(), this->repository() ) );
    M_electricModel.reset( new electric_model_type(prefixvm(this->prefix(),"electric"), "electric", this->worldCommPtr(),
                                                   this->subPrefix(), this->repository() ) );

    if ( this->physics().empty() )
    {
        typename ModelPhysics<mesh_type::nDim>::subphysic_description_type subPhyicsDesc;
        subPhyicsDesc[M_heatModel->physicType()] = std::make_tuple( M_heatModel->keyword(), M_heatModel );
        subPhyicsDesc[M_electricModel->physicType()] = std::make_tuple( M_electricModel->keyword(), M_electricModel );
        this->initPhysics( this->keyword(), this->modelProperties().models(), subPhyicsDesc );
    }

    // physical properties
    if ( !M_materialsProperties )
    {
        //auto paramValues = this->modelProperties().parameters().toParameterValues();
        //this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    if ( !this->mesh() )
        this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    // init heat toolbox
    M_heatModel->setPhysics( this->physics( M_heatModel->physicType() ), M_heatModel->keyword() );
    M_heatModel->setManageParameterValues( false );
    if ( !M_heatModel->modelPropertiesPtr() )
    {
        M_heatModel->setModelProperties( this->modelPropertiesPtr() );
        M_heatModel->setManageParameterValuesOfModelProperties( false );
    }
    M_heatModel->setMesh( this->mesh() );
    M_heatModel->setMaterialsProperties( M_materialsProperties );
    M_heatModel->init( false );

    // init electric toolbox
    M_electricModel->setPhysics( this->physics( M_electricModel->physicType() ), M_electricModel->keyword() );
    M_electricModel->setManageParameterValues( false );
    if ( !M_electricModel->modelPropertiesPtr() )
    {
        M_electricModel->setModelProperties( this->modelPropertiesPtr() );
        M_electricModel->setManageParameterValuesOfModelProperties( false );
    }
    M_electricModel->setMesh( this->mesh() );
    M_electricModel->setMaterialsProperties( M_materialsProperties );
    M_electricModel->init( false );

    M_modelName = "ThermoElectric";
    if ( M_solverName == "automatic" )
    {
        if ( this->materialsProperties()->hasElectricConductivityDependingOnSymbol( "heat_T" ) )
            M_solverName = "Newton";
        else
            M_solverName = "Linear";
    }
    M_modelUseJouleEffect = true;

    if ( M_solverName == "Linear" || M_solverNewtonInitialGuessUseLinearHeat )
    {
        M_heatModel->initAlgebraicFactory();
        M_heatModel->algebraicFactory()->setFunctionLinearAssembly( boost::bind( &self_type::updateLinear_Heat,
                                                                                 boost::ref( *this ), _1 ) );
        M_heatModel->algebraicFactory()->setFunctionResidualAssembly( boost::bind( &self_type::updateResidual_Heat,
                                                                                   boost::ref( *this ), _1 ) );
    }
    if ( M_solverName == "Linear" || M_solverNewtonInitialGuessUseLinearElectric )
    {
        M_electricModel->initAlgebraicFactory();
        M_electricModel->algebraicFactory()->setFunctionLinearAssembly( boost::bind( &self_type::updateLinear_Electric,
                                                                                     boost::ref( *this ), _1 ) );
    }


#if 0
    M_rangeMeshElements = ( M_heatModel->thermalProperties()->isDefinedOnWholeMesh() && M_electricModel->electricProperties()->isDefinedOnWholeMesh() )?
        elements(this->mesh() ) :
        intersect( M_heatModel->rangeMeshElements(), M_electricModel->rangeMeshElements() );
#endif

    // post-process
    this->initPostProcess();

    // update constant parameters into
    this->updateParameterValues();

    // backend
    this->initAlgebraicBackend();

    // block vector solution
    auto const& blockVectorSolutionHeat = *M_heatModel->algebraicBlockVectorSolution();
    auto const& blockVectorSolutionElectric = *M_electricModel->algebraicBlockVectorSolution();
    int nBlockHeat = blockVectorSolutionHeat.size();
    int nBlockElectric = blockVectorSolutionElectric.size();
    int nBlock = nBlockHeat + nBlockElectric;
    auto bvs = this->initAlgebraicBlockVectorSolution( nBlock );
    int indexBlock=0;
    int numberOfBlockSpaceHeat = 0;
    for ( int k=0;k<nBlockHeat ;++k )
    {
        bvs->operator()(indexBlock+k) = blockVectorSolutionHeat(k);
        numberOfBlockSpaceHeat += blockVectorSolutionHeat(k)->map().numberOfDofIdToContainerId();
    }
    indexBlock += nBlockHeat;
    for ( int k=0;k<nBlockElectric ;++k )
        bvs->operator()(indexBlock+k) = blockVectorSolutionElectric(k);
    indexBlock += nBlockElectric;
    // init monolithic vector associated to the block vector
    bvs->buildVector( this->backend() );

    size_type currentStartBlockSpaceIndex = 0;
    this->setStartSubBlockSpaceIndex( "heat", currentStartBlockSpaceIndex );
    currentStartBlockSpaceIndex += numberOfBlockSpaceHeat;
    this->setStartSubBlockSpaceIndex( "electric", currentStartBlockSpaceIndex );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
    {
        if ( M_solverName == "Newton" || M_solverName == "Picard" )
        {
            auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
            this->setAlgebraicFactory( algebraicFactory );
        }
    }

    this->setIsUpdatedForUse( true );

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


    // need to not include export fields of material of subphysics
    std::set<std::string> ppExportsAllFieldsAvailableHeat = Feel::FeelModels::detail::set_difference( this->heatModel()->postProcessExportsAllFieldsAvailable(),
                                                                                                      this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->heatModel()->physicsAvailable() ) );
    std::set<std::string> ppExportsAllFieldsAvailableElectric = Feel::FeelModels::detail::set_difference( this->electricModel()->postProcessExportsAllFieldsAvailable(),
                                                                                                          this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->electricModel()->physicsAvailable() ) );
    std::set<std::string> ppExportsAllFieldsAvailable;
    for ( auto const& s : ppExportsAllFieldsAvailableHeat )
        ppExportsAllFieldsAvailable.insert( prefixvm( this->heatModel()->keyword(), s) );
    for ( auto const& s : ppExportsAllFieldsAvailableElectric )
        ppExportsAllFieldsAvailable.insert( prefixvm( this->electricModel()->keyword(), s) );

    this->setPostProcessExportsAllFieldsAvailable( ppExportsAllFieldsAvailable );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    super_type::initPostProcess();

    if ( !this->postProcessExportsFields().empty() )
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
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    // p.put( "toolbox-heat", M_heatModel->journalSectionName() );
    // p.put( "toolbox-electric", M_electricModel->journalSectionName() );

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );

    // Numerical Solver
    nl::json subPt;
    subPt.emplace( "solver", M_solverName );
    p["Numerical Solver"] = subPt;

    // Exporter
#if 0
    if ( M_exporter )
    {
        subPt.clear();
        subPt.put( "type",M_exporter->type() );
        subPt.put( "freq save",M_exporter->freq() );
        pt::ptree subPt2;
        for ( std::string const& fieldName : this->postProcessExportsFields() )
            subPt2.push_back( std::make_pair("", pt::ptree( fieldName ) ) );
        subPt.put_child( "fields", subPt2 );
        p.put_child( "Exporter", subPt );
    }
#endif

    // Algebraic Solver
    if ( this->algebraicFactory() )
        this->algebraicFactory()->updateInformationObject( p["Algebraic Solver"] );

    p["Toolbox Heat"] = M_heatModel->journalSection().to_string();
    p["Toolbox Electric"] = M_electricModel->journalSection().to_string();
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfo->add( "Materials Properties", this->materialsProperties()->tabulateInformations(jsonInfo.at("Materials Properties"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    // Numerical Solver
    if ( jsonInfo.contains( "Numerical Solver" ) )
    {
        Feel::Table tabInfoNumSolver;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoNumSolver, jsonInfo.at("Numerical Solver"), tabInfoProp );
        tabInfo->add( "Numerical Solver",  TabulateInformations::New( tabInfoNumSolver, tabInfoProp ) );
    }

    if ( jsonInfo.contains( "Algebraic Solver" ) )
        tabInfo->add( "Algebraic Solver", model_algebraic_factory_type::tabulateInformations( jsonInfo.at("Algebraic Solver"), tabInfoProp ) );

    // generate sub toolboxes info
    if ( M_heatModel && jsonInfo.contains( "Toolbox Heat" ) )
    {
        nl::json::json_pointer jsonPointerHeat( jsonInfo.at( "Toolbox Heat" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerHeat ) )
        {
            auto tabInfos_heat = M_heatModel->tabulateInformations( JournalManager::journalData().at( jsonPointerHeat ), tabInfoProp );
            TabulateInformationsSections::cast( tabInfos_heat )->erase( "Materials Properties" );
            tabInfo->add( "Toolbox Heat", tabInfos_heat );
        }
    }
    if ( M_electricModel && jsonInfo.contains( "Toolbox Electric" ) )
    {
        nl::json::json_pointer jsonPointerElectric( jsonInfo.at( "Toolbox Electric" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerElectric ) )
        {
            auto tabInfos_electric = M_electricModel->tabulateInformations( JournalManager::journalData().at( jsonPointerElectric ), tabInfoProp );
            TabulateInformationsSections::cast( tabInfos_electric )->erase( "Materials Properties" );
            tabInfo->add( "Toolbox Electric", tabInfos_electric );
        }
    }

    return tabInfo;
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
#if 0
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
        std::string fieldExported;
        for ( std::string const& fieldName : this->postProcessExportsFields() )
            fieldExported=(fieldExported.empty())? fieldName : fieldExported + " - " + fieldName;
        *_ostr << "\n     -- fields [heat] : " << fieldExported;
    }
    *_ostr << "\n   Processors"
           << "\n     -- number of proc : " << this->worldComm().globalSize()
           << "\n     -- current rank : " << this->worldComm().globalRank();

    if ( this-<algebraicFactory() )
        *_ostr << this->algebraicFactory()->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";
#endif
    return _ostr;
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    this->heatModel()->startTimeStep();
    this->updateTime( this->heatModel()->time() );
    this->updateParameterValues();
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->heatModel()->updateTimeStep();
    this->updateTime( this->heatModel()->time() );
    this->updateParameterValues();
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("ThermoElectric","exportResults", "start");
    this->timerTool("PostProcessing").start();

    auto mfields = this->modelFields();
    auto symbolExpr = this->symbolsExpr( mfields );
    //std::cout << "holalla \n "<< symbolExpr.names() << std::endl;
    M_heatModel->exportResults( time, symbolExpr );
    M_electricModel->exportResults( time, symbolExpr );

    auto exprExport =  hana::concat( M_materialsProperties->exprPostProcessExports( this->mesh(),this->physicsAvailable(),symbolExpr ),
                                     hana::concat( M_heatModel->exprPostProcessExportsToolbox( symbolExpr,M_heatModel->keyword() ),
                                                   M_electricModel->exprPostProcessExportsToolbox( symbolExpr,M_electricModel->keyword() ) ) );
    this->executePostProcessExports( M_exporter, time, mfields, symbolExpr, exprExport );

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
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );

    this->setParameterValues( paramValues );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }
    M_heatModel->setParameterValues( paramValues );
    M_electricModel->setParameterValues( paramValues );
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("ThermoElectric","solve", "start");
    this->timerTool("Solve").start();

    this->setStartBlockSpaceIndex( 0 );

    if ( M_solverName == "Linear" )
    {
        M_electricModel->solve();
        M_heatModel->solve();
        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
    }
    else if ( M_solverName == "Newton" || M_solverName == "Picard" )
    {
        // initial guess
        if ( M_solverNewtonInitialGuessUseLinearElectric )
            M_electricModel->solve();
        if ( M_solverNewtonInitialGuessUseLinearHeat )
            M_heatModel->solve();

        // solve non linear monolithic system
        M_heatModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("heat") );
        M_electricModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("electric") );
        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
        this->algebraicFactory()->solve( M_solverName, this->algebraicBlockVectorSolution()->vectorMonolithic() );
        this->algebraicBlockVectorSolution()->localize();
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

} // end namespace FeelModels
} // end namespace Feel
