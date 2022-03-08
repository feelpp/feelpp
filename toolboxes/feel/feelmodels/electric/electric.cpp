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

#include <feel/feelmodels/electric/electric.hpp>

#include <feel/feelvf/vf.hpp>

namespace Feel
{
namespace FeelModels
{

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
ELECTRIC_CLASS_TEMPLATE_TYPE::Electric( std::string const& prefix,
                                        std::string const& keyword,
                                        worldcomm_ptr_t const& worldComm,
                                        std::string const& subPrefix,
                                        ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelPhysics<nDim>( "electric" ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep )
{
    this->log("Electric","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ElectricConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ElectricSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ElectricPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ElectricTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    this->log("Electric","constructor", "finish");
}


ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{

}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("Electric","initMesh", "start");
    this->timerTool("Constructor").start();

    if ( auto ptMeshes = this->modelProperties().pTree().get_child_optional("Meshes") )
        super_type::super_model_meshes_type::setup( *ptMeshes, {this->keyword()} );
    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( this->mesh() ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("createMesh");
    this->log("Electric","initMesh",(boost::format("finish in %1% s")%tElpased).str() );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
int
ELECTRIC_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = 1;
    return nBlock;
}


ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
ELECTRIC_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    myblockGraph(0,0) = stencil(_test=this->spaceElectricPotential(),
                                _trial=this->spaceElectricPotential() )->graph();
    return myblockGraph;
}


ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("Electric","init", "start" );
    this->timerTool("Constructor").start();

    // physics
    this->initPhysics( this->keyword(), this->modelProperties().models() );

    // physical properties
    if ( !M_materialsProperties )
    {
        auto paramValues = this->modelProperties().parameters().toParameterValues();
        this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    //if ( !this->mesh() )
    this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    auto mom = this->materialsProperties()->materialsOnMesh(this->mesh());
    // functionspace
    if ( mom->isDefinedOnWholeMesh( this->physicsAvailableFromCurrentType() ) )
    {
        M_rangeMeshElements = elements(this->mesh());
        M_XhElectricPotential = space_electricpotential_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), mom->markers( this->physicsAvailableFromCurrentType() ));
        M_XhElectricPotential = space_electricpotential_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
    }
    M_fieldElectricPotential.reset( new element_electricpotential_type(M_XhElectricPotential,"V"));

    this->initBoundaryConditions();

    this->initInitialConditions();

    // post-process
    this->initPostProcess();

    // update constant parameters
    this->updateParameterValues();

    // backend
    this->initAlgebraicBackend();

    size_type currentStartIndex = 0;
    this->setStartSubBlockSpaceIndex( "electric-potential", currentStartIndex );

    // vector solution
    int nBlock = this->nBlockMatrixGraph();
    auto bvs = this->initAlgebraicBlockVectorSolution( nBlock );
    bvs->operator()(0) = this->fieldElectricPotentialPtr();
    // init petsc vector associated to the block
    bvs->buildVector( this->backend() );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
        this->initAlgebraicFactory();

    this->setIsUpdatedForUse( true );

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("Electric","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::initInitialConditions()
{
    if ( !this->doRestart() )
    {
        std::vector<element_electricpotential_ptrtype> icElectricPotentialFields;
        CHECK( this->isStationary() ) << "TODO";
        icElectricPotentialFields = { this->fieldElectricPotentialPtr() };

        auto paramValues = this->modelProperties().parameters().toParameterValues();
        this->modelProperties().initialConditions().setParameterValues( paramValues );

        this->updateInitialConditions( "electric-potential", M_rangeMeshElements, this->symbolsExpr(), icElectricPotentialFields );
    }
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    if ( !this->modelProperties().boundaryConditions().hasSection( this->keyword() ) )
        return;
    M_boundaryConditions.setup( *this,this->modelProperties().boundaryConditions().section( this->keyword() ) );

    for ( auto const& [bcName,bcData] : M_boundaryConditions.electricPotentialImposed() )
        bcData->updateDofEliminationIds( *this, "electric-potential", this->spaceElectricPotential() );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("Electric","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->setPostProcessExportsAllFieldsAvailable( {"electric-potential","electric-field","current-density","joules-losses"} );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    this->setPostProcessSaveAllFieldsAvailable( {"electric-potential","electric-field","electric-conductivity","current-density","joules-losses"} );
    super_type::initPostProcess();

    if ( !this->postProcessExportsFields().empty() )
    {
        std::string geoExportType="static";//change_coords_only, change, static
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export",
                               _geo=geoExportType,
                               _path=this->exporterPath() );

        if (this->doRestart() && this->restartPath().empty() )
        {
            if ( M_exporter->doExport() )
                M_exporter->restart(this->timeInitial());
        }
    }

    auto se = this->symbolsExpr();
    this->template initPostProcessMeshes<mesh_type>( se );

    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasures().restart( this->timeInitial() );
    }

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("Electric","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
    this->setAlgebraicFactory( algebraicFactory );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    super_physics_type::updateInformationObjectFromCurrentType( p["Physics"] );

    // Physics
    nl::json subPt;
    subPt.emplace( "time mode", std::string( (this->isStationary())?"Stationary":"Transient") );
    p["Physics2"] = subPt;

    // Boundary Conditions
    M_boundaryConditions.updateInformationObject( p["Boundary Conditions"] );

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );

    // Function Spaces
    subPt.clear();
    subPt["Electric Potential"] = M_XhElectricPotential->journalSection().to_string();
    p.emplace( "Function Spaces", subPt );

    // fields
    this->modelFields().updateInformationObject( p["Fields"] );

    // Algebraic Solver
    if ( this->algebraicFactory() )
        this->algebraicFactory()->updateInformationObject( p["Algebraic Solver"] );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
ELECTRIC_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics") )
        tabInfo->add( "Physics", super_physics_type::tabulateInformations( jsonInfo.at("Physics"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics2") )
    {
        Feel::Table tabInfoPhysics;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoPhysics, jsonInfo.at("Physics2"), tabInfoProp );
        tabInfo->add( "Physics2", TabulateInformations::New( tabInfoPhysics, tabInfoProp ) );
    }

    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfo->add( "Materials Properties", this->materialsProperties()->tabulateInformations(jsonInfo.at("Materials Properties"), tabInfoProp ) );

    if ( jsonInfo.contains("Boundary Conditions") )
        tabInfo->add( "Boundary Conditions", ElectricBoundaryConditions::tabulateInformations( jsonInfo.at("Boundary Conditions"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    if ( jsonInfo.contains("Function Spaces") )
    {
        auto const& jsonInfoFunctionSpaces = jsonInfo.at("Function Spaces");
        auto tabInfoFunctionSpaces = TabulateInformationsSections::New( tabInfoProp );
        nl::json::json_pointer jsonPointerSpaceElectricPotential( jsonInfoFunctionSpaces.at( "Electric Potential" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerSpaceElectricPotential ) )
            tabInfoFunctionSpaces->add( "Electric Potential", TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( JournalManager::journalData().at( jsonPointerSpaceElectricPotential ), tabInfoProp ) );
        tabInfo->add( "Function Spaces", tabInfoFunctionSpaces );
    }

    // fields
    if ( jsonInfo.contains("Fields") )
        tabInfo->add( "Fields", TabulateInformationTools::FromJSON::tabulateInformationsModelFields( jsonInfo.at("Fields"), tabInfoProp ) );

    if ( jsonInfo.contains( "Algebraic Solver" ) )
        tabInfo->add( "Algebraic Solver", model_algebraic_factory_type::tabulateInformations( jsonInfo.at("Algebraic Solver"), tabInfoProp ) );

    return tabInfo;
}


ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
ELECTRIC_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
#if 0
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||-----------Info : Electric--------------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository();
    *_ostr << "\n   Physical Model"
           << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient");
    *_ostr << "\n   Boundary conditions"
           << M_bcDirichletMarkerManagement.getInfoDirichletBC()
           << M_bcNeumannMarkerManagement.getInfoNeumannBC()
           << M_bcRobinMarkerManagement.getInfoRobinBC();
#if 0
    *_ostr << this->materialsProperties()->getInfoMaterialParameters()->str();

    *_ostr << "\n   Mesh Discretization"
           << "\n     -- mesh filename      : " << this->meshFile()
           << "\n     -- number of element : " << M_mesh->numGlobalElements()
           << "\n     -- order             : " << nOrderGeo;
    *_ostr << "\n   Space ElectricPotential Discretization"
           << "\n     -- order         : " << nOrderPolyElectricPotential
           << "\n     -- number of dof : " << M_XhElectricPotential->nDof() << " (" << M_XhElectricPotential->nLocalDof() << ")";

    if ( this->algebraicFactory() )
        *_ostr << this->algebraicFactory()->getInfo()->str();
#endif
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";
#endif
    return _ostr;
}



ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    this->exportResults( time, mfields, se, this->exprPostProcessExports( se ) );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );
    for ( auto [physicName,physicData] : this->physics/*FromCurrentType*/() )
        physicData->updateParameterValues( paramValues );

    this->updateParameterValues_postProcess( paramValues, prefixvm("postprocess",this->keyword(),"_" ) );

    this->setParameterValues( paramValues );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    this->log("Electric","setParameterValues", "start");

    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->modelProperties().initialConditions().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        physicData->setParameterValues( paramValues );

    M_boundaryConditions.setParameterValues( paramValues );

    this->log("Electric","setParameterValues", "finish");
}


ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Electric","solve", "start");
    this->timerTool("Solve").start();

    this->setStartBlockSpaceIndex( 0 );

    this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
    this->algebraicFactory()->solve( "LinearSystem", this->algebraicBlockVectorSolution()->vectorMonolithic() );
    this->algebraicBlockVectorSolution()->localize();

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("Electric","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}


} // end namespace FeelModels
} // end namespace Feel
