/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Romain Hild <romain.hild@cemosis.fr>
       Date: 2021-05-04

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
   \file mixedpoisson.hpp
   \author Romain Hild <romain.hild@cemosis.fr>
   \date 2021-05-04
 */

#include <feel/feelmodels/hdg/mixedpoisson.hpp>
#include <feel/feelmesh/complement.hpp>
#include <numeric>
#include <algorithm>

namespace Feel
{
namespace FeelModels
{

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::MixedPoisson( std::string const& prefix,
                                                MixedPoissonPhysics const& physic,
                                                worldcomm_ptr_t const& _worldComm,
                                                std::string const& subPrefix,
                                                ModelBaseRepository const& modelRep )
    : super_type( prefix, MixedPoissonPhysicsMap[physic]["keyword"], _worldComm, subPrefix, modelRep),
      ModelPhysics<nDim>( MixedPoissonPhysicsMap[physic]["keyword"] ),
      ModelBase( prefix, MixedPoissonPhysicsMap[physic]["keyword"], _worldComm, subPrefix, modelRep),
      M_physic(physic),
      M_physicMap(MixedPoissonPhysicsMap[physic]),
      M_potentialKey(MixedPoissonPhysicsMap[physic]["potentialK"]),
      M_fluxKey(MixedPoissonPhysicsMap[physic]["fluxK"]),
      M_tauCst(doption( prefixvm(this->prefix(), "tau_constant") )),
      M_useSC(boption( prefixvm(this->prefix(), "use-sc")) ),
      M_postMatrixInit(false),
      M_solverName(soption(_prefix=this->prefix(), _name="solver")),
      M_useNearNullSpace(boption(_prefix=this->prefix(), _name="use-near-null-space"))
{
    this->log("MixedPoisson","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedPoissonConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedPoissonSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedPoissonPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedPoissonTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    this->log("MixedPoisson","constructor", "finish");
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    // time stepping
    M_timeStepping = soption(_name="time-stepping",_prefix=this->prefix());
    M_timeStepThetaValue = doption(_name="time-stepping.theta.value",_prefix=this->prefix());
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
int
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = 1;
    return nBlock;
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    myblockGraph(0,0) = stencil(_test=this->spacePotential(),
                                _trial=this->spacePotential() )->graph();
    return myblockGraph;
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
    this->setAlgebraicFactory( algebraicFactory );
}



MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("MixedPoisson","init", "start" );
    this->timerTool("Constructor").start();

    this->initModelProperties();

    this->initPhysics( this->shared_from_this(), this->modelProperties().models() );

    // physical properties
    if ( !M_materialsProperties )
    {
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    this->initBoundaryConditions();

    this->initFunctionSpaces();

    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();

    // update constant parameters
    this->updateParameterValues();

    // update initial conditions
    this->updateInitialConditions( this->symbolsExpr() );

    // post-process
    this->initPostProcess();

    // backend
    this->initAlgebraicBackend();

    // algebraic solver
    if ( buildModelAlgebraicFactory )
        this->initAlgebraicFactory();

    solve::strategy s = M_useSC ? solve::strategy::static_condensation : solve::strategy::monolithic;
    auto pps = product( M_Whp );

    M_A = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(*M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), *this->backend(), (s>=solve::strategy::static_condensation)?false:true );
    M_F = makeSharedVectorCondensed<value_type>(s, blockVector(*M_ps), *this->backend(), false);
    M_App = makeSharedMatrixCondensed<value_type>(solve::strategy::local,  csrGraphBlocks(pps, Pattern::ZERO), this->backend(), false );
    M_Fpp = makeSharedVectorCondensed<value_type>(solve::strategy::local, blockVector(pps), this->backend(), false);

    this->setIsUpdatedForUse( true );

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("MixedPoisson","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    auto mom = this->materialsProperties()->materialsOnMesh(this->mesh());
    // functionspace
    if ( mom->isDefinedOnWholeMesh( this->physicsAvailableFromCurrentType() ) )
    {
        M_rangeMeshElements = elements(this->mesh());
        M_Vh = space_flux_type::New( _mesh=this->mesh(), _extended_doftable=true, _worldscomm=this->worldsComm() );
        M_Wh = space_potential_type::New( _mesh=this->mesh(), _extended_doftable=true, _worldscomm=this->worldsComm() );
        M_Whp = space_postpotential_type::New( _mesh=this->mesh(), _extended_doftable=true, _worldscomm=this->worldsComm() );
        M_P0dh = space_p0dh_type::New( _mesh=this->mesh(), _extended_doftable=true, _worldscomm=this->worldsComm() );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), mom->markers( this->physicsAvailableFromCurrentType() ));
        M_Vh = space_flux_type::New( _mesh=this->mesh(), _extended_doftable=true, _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
        M_Wh = space_potential_type::New( _mesh=this->mesh(), _extended_doftable=true, _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
        M_Whp = space_postpotential_type::New( _mesh=this->mesh(), _extended_doftable=true, _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
        M_P0dh = space_p0dh_type::New( _mesh=this->mesh(), _extended_doftable=true, _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
    }
    M_up = std::make_shared<element_flux_type>(M_Vh, M_physicMap["fluxSymbol"]);
    M_pp = std::make_shared<element_potential_type>(M_Wh, M_physicMap["potentialSymbol"]);
    M_ppp = std::make_shared<element_postpotential_type>(M_Whp, "post"+M_physicMap["potentialSymbol"]);

#if 0
    auto ibcMarkers = std::accumulate(M_boundaryConditions->integral().begin(),
                                      M_boundaryConditions->integral().end(),
                                      std::set<std::string>(),
                                      [](auto& prev, auto const& pair) {
                                          prev.insert(pair.second->markers().begin(), pair.second->markers().end());
                                          return prev;
                                      });
#else
    std::set<std::string> ibcMarkers;
    for ( auto const& [bcName,bcData] : M_boundaryConditions->integral() )
        ibcMarkers.insert( bcData->markers().begin(),bcData->markers().end() );
    for ( auto const& [bcName,bcData] : M_boundaryConditions->couplingODEs() )
        ibcMarkers.insert( bcData->markers().begin(),bcData->markers().end() );
#endif
    std::set<int> ibcMeshMarkers;
    std::for_each(ibcMarkers.begin(), ibcMarkers.end(),
                  [this,&ibcMeshMarkers](auto const& x) {
                      ibcMeshMarkers.insert(this->mesh()->markerName(x));
                  });
    auto complement_integral_faces = complement(faces(support(M_Wh)),
                                                [ibcMeshMarkers]( auto const& ewrap ) {
                                                    auto const& e = unwrap_ref( ewrap );
                                                    return ( e.hasMarker() && ibcMeshMarkers.count(e.marker().value()) ); // WARNING now we can have multiple marker values by entity
                                                });
    M_gammaMinusIntegral = complement(boundaryfaces(support(M_Wh)),
                                      [ibcMeshMarkers]( auto const& ewrap ) {
                                          auto const& e = unwrap_ref( ewrap );
                                          return ( e.hasMarker() && ibcMeshMarkers.count(e.marker().value()) );
                                      });
    auto face_mesh = createSubmesh( _mesh=this->mesh(), _range=complement_integral_faces, _update=0 );
    M_Mh = space_trace_type::New( _mesh=face_mesh, _extended_doftable=true, _worldscomm=this->worldsComm() );
    M_phat = std::make_shared<element_trace_type>(M_Mh, "phat");

    auto ibc_mesh = createSubmesh( _mesh=this->mesh(), _range=markedfaces(this->mesh(), ibcMarkers), _update=0 );
    M_Ch = space_traceibc_type::New( _mesh=ibc_mesh, _extended_doftable=true, _worldscomm=this->worldsComm() );
    M_mup = element_traceibc_vector_type(this->constantSpacesSize(),
                                         std::make_shared<element_traceibc_type>(M_Ch, "mup"));

    auto ibcSpaces = std::make_shared<ProductSpace<space_traceibc_ptrtype, true> >( this->constantSpacesSize(), M_Ch);
    this->setSpaceProperties(ibcSpaces);
    M_ps = std::make_shared<product2_space_type>(ibcSpaces, M_Vh, M_Wh, M_Mh);
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::setSpaceProperties(product_space_ptrtype const& ibcSpaces)
{
    std::vector<std::string> props(this->constantSpacesSize(), "Ibc");
    ibcSpaces->setProperties( props );
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("MixedPoisson","initMesh", "start");
    this->timerTool("Constructor").start();

    if ( this->modelProperties().jsonData().contains("Meshes") )
        super_type::super_model_meshes_type::setup( this->modelProperties().jsonData().at("Meshes"), {this->keyword()} );
    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( this->mesh() ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("createMesh");
    this->log("MixedPoisson","initMesh",(boost::format("finish in %1% s")%tElpased).str() );
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    M_boundaryConditions = std::make_shared<boundary_conditions_type>( this->shared_from_this() );
    if ( !this->modelProperties().boundaryConditions().hasSection( this->keyword() ) )
        return;

    M_boundaryConditions->setup( this->modelProperties().boundaryConditions().section( this->keyword() ) );
#if 0
    for ( auto const& [bcName,bcData] : M_boundaryConditions->dirichlet() )
        bcData->updateDofEliminationIds( *this, this->unknownName(), this->spaceUnknown() );
#endif
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    this->log("MixedPoisson","initTimeStep", "start" );
    this->timerTool("Constructor").start();

    double ti = this->timeInitial();
    double tf = this->timeFinal();
    double dt = this->timeStep();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();

    if ( M_timeStepping == "BDF" || M_timeStepping == "Theta" )
    {
        int bdfOrder = 1;
        if ( M_timeStepping == "BDF" )
            bdfOrder = ioption(_prefix=this->prefix(),_name="bdf.order");
        int nConsecutiveSave = std::max( 3, bdfOrder ); // at least 3 is required when restart with theta scheme
        M_bdfPotential = this->createBdf( this->spacePotential(),M_potentialKey, bdfOrder, nConsecutiveSave, myFileFormat );
    }
    else if( M_timeStepping == "Newmark" )
    {
        M_newmarkPotential = newmark( _space=this->spacePotential(),
                                      _name=M_potentialKey+suffixName,
                                      _prefix=this->prefix(),
                                      _initial_time=ti, _final_time=tf, _time_step=dt,
                                      _restart=this->doRestart(), _restart_path=this->restartPath(), _restart_at_last_save=this->restartAtLastSave(),
                                      _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
    }

    if (!this->doRestart())
    {
        // up current time
        if ( M_timeStepping == "Newmark" )
            this->updateTime( M_newmarkPotential->timeInitial() );
        else
            this->updateTime( M_bdfPotential->timeInitial() );
    }
    else
    {
        double tir;
        if ( M_timeStepping == "Newmark" )
        {
            // restart time step
            tir = M_newmarkPotential->restart();
            // load a previous solution as current solution
            *this->fieldPotentialPtr() = M_newmarkPotential->previousUnknown();
        }
        else
        {
            // start time step
            tir = M_bdfPotential->restart();
            // load a previous solution as current solution
            *this->fieldPotentialPtr() = M_bdfPotential->unknown(0);
        }
        // up initial time
        this->setTimeInitial( tir );
        // up current time
        this->updateTime( tir );
    }

    double tElapsed = this->timerTool("Constructor").stop("initTimeStep");
    this->log("MixedPoisson","initTimeStep", (boost::format("finish in %1% s") %tElapsed).str() );
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("MixedPoisson","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->setPostProcessExportsAllFieldsAvailable( {M_potentialKey, M_fluxKey, "post"+M_potentialKey} );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    this->setPostProcessSaveAllFieldsAvailable( {M_potentialKey } );
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
    this->log("MixedPoisson","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );

    this->setParameterValues( paramValues );
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    this->log("MixedPoisson","setParameterValues", "start");

    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        physicData->setParameterValues( paramValues );

    M_boundaryConditions->setParameterValues( paramValues );

    this->log("MixedPoisson","setParameterValues", "finish");
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    this->exportResults( time, mfields, se, this->exprPostProcessExports( se ) );
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    super_physics_type::updateInformationObjectFromCurrentType( p["Physics"] );

    M_boundaryConditions->updateInformationObject( p["Boundary Conditions"] );
    // // Physics
    nl::json subPt;
    // subPt.emplace( "time mode", std::string( (this->isStationary())?"Stationary":"Transient") );
    // p["Physics"] = subPt;

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );


    // Function Spaces
    subPt.clear();
    subPt["Potential"] = M_Wh->journalSection().to_string();
    subPt["Flux"] = M_Vh->journalSection().to_string();
    subPt["Trace"] = M_Mh->journalSection().to_string();
    p.emplace( "Function Spaces", subPt );

    if ( !this->isStationary() )
    {
        subPt.clear();
        subPt.emplace( "initial time", this->timeStepBase()->timeInitial() );
        subPt.emplace( "final time", this->timeStepBase()->timeFinal() );
        subPt.emplace( "time step", this->timeStepBase()->timeStep() );
        subPt.emplace( "type", M_timeStepping );
        p["Time Discretization"] = subPt;
    }

    // // Algebraic Solver
    // if ( M_algebraicFactory )
    //     M_algebraicFactory->updateInformationObject( p["Algebraic Solver"] );
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics") )
        tabInfo->add( "Physics", super_physics_type::tabulateInformations( jsonInfo.at("Physics"), tabInfoProp ) );

    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfo->add( "Materials Properties", this->materialsProperties()->tabulateInformations(jsonInfo.at("Materials Properties"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    if ( jsonInfo.contains("Boundary Conditions") )
        tabInfo->add( "Boundary Conditions", boundary_conditions_type::tabulateInformations( jsonInfo.at("Boundary Conditions"), tabInfoProp ) );

    //tabInfoSections.push_back( std::make_pair( "Boundary conditions",  tabulate::Table{} ) );


    if ( jsonInfo.contains("Function Spaces") )
    {
        auto const& jsonInfoFunctionSpaces = jsonInfo.at("Function Spaces");
        auto tabInfoFunctionSpaces = TabulateInformationsSections::New( tabInfoProp );
        nl::json::json_pointer jsonPointerSpacePotential( jsonInfoFunctionSpaces.at( "Potential" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerSpacePotential ) )
            tabInfoFunctionSpaces->add( "Potential", TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( JournalManager::journalData().at( jsonPointerSpacePotential ), tabInfoProp ) );
        nl::json::json_pointer jsonPointerSpaceFlux( jsonInfoFunctionSpaces.at( "Flux" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerSpaceFlux ) )
            tabInfoFunctionSpaces->add( "Flux", TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( JournalManager::journalData().at( jsonPointerSpaceFlux ), tabInfoProp ) );
        nl::json::json_pointer jsonPointerSpaceTrace( jsonInfoFunctionSpaces.at( "Trace" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerSpaceTrace ) )
            tabInfoFunctionSpaces->add( "Trace", TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( JournalManager::journalData().at( jsonPointerSpaceTrace ), tabInfoProp ) );
        tabInfo->add( "Function Spaces", tabInfoFunctionSpaces );
    }

    if ( jsonInfo.contains("Time Discretization") )
    {
        Feel::Table tabInfoTimeDiscr;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoTimeDiscr, jsonInfo.at("Time Discretization"), tabInfoProp );
        tabInfo->add( "Time Discretization", TabulateInformations::New( tabInfoTimeDiscr, tabInfoProp ) );
    }

    // if ( jsonInfo.contains( "Algebraic Solver" ) )
    //     tabInfo->add( "Algebraic Solver", model_algebraic_factory_type::tabulateInformations( jsonInfo.at("Algebraic Solver"), tabInfoProp ) );

    return tabInfo;

}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<NullSpace<double>>
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::nullSpace(std::string const& name)
{
    if constexpr( is_tensor2symm )
    {
        auto mode1 = this->M_Mh->element( oneX() );
        auto mode2 = this->M_Mh->element( oneY() );
        if constexpr( nDim == 2 )
        {
            auto mode3 = this->M_Mh->element( vec(Py(),-Px()) );
            auto ns = NullSpace<double>({ mode1, mode2, mode3 });
            return std::make_shared<NullSpace<double>>(Feel::backend(_name=name), ns);
        } else if constexpr( nDim == 3 )
        {
            auto mode3 = this->M_Mh->element( oneZ() );
            auto mode4 = this->M_Mh->element( vec(Py(),-Px(),cst(0.)) );
            auto mode5 = this->M_Mh->element( vec(-Pz(),cst(0.),Px()) );
            auto mode6 = this->M_Mh->element( vec(cst(0.),Pz(),-Py()) );
            auto ns = NullSpace<double>({ mode1, mode2, mode3, mode4, mode5, mode6 });
            return std::make_shared<NullSpace<double>>(Feel::backend(_name=name), ns);
        }
        else
            return nullptr;
    }
    else
        return nullptr;
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::solve()
{
    if( M_solverName == "Linear" )
        this->solveLinear();
    else
        this->solvePicard();
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::solveLinear()
{
    solve::strategy s = M_useSC ? solve::strategy::static_condensation : solve::strategy::monolithic;
    M_A = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(*M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), *this->backend(), (s>=solve::strategy::static_condensation)?false:true );
    M_F = makeSharedVectorCondensed<value_type>(s, blockVector(*M_ps), *this->backend(), false);

    // M_A->zero();
    // M_F->zero();

    auto U = M_ps->element();
    U.buildVector(this->backend());

    auto A = std::dynamic_pointer_cast<typename super_type::backend_type::sparse_matrix_type>(M_A);
    auto F = std::dynamic_pointer_cast<typename super_type::backend_type::vector_type>(M_F);
    CHECK(A) << "Dynamic cast for M_A not ok !!!";
    CHECK(F) << "Dynamic cast for M_F not ok !!!";
    this->algebraicFactory()->applyAssemblyLinear(U.vectorMonolithic(), A, F);

    auto bbf = blockform2( *M_ps, M_A);
    auto blf = blockform1( *M_ps, M_F);

    if constexpr( is_tensor2symm )
    {
        if( M_useNearNullSpace )
        {
            auto name = this->prefix()+".sc";
            Feel::backend( _name=name )->attachNearNullSpace( this->nullSpace(name) );
        }
    }

    bbf.solve(_solution=U, _rhs=blf, _condense=M_useSC, _name=this->prefix());

    M_up = std::make_shared<element_flux_type>( U(0_c) );
    M_pp = std::make_shared<element_potential_type>( U(1_c));
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::solvePicard()
{
    CHECK(false) << "Solve Picard not yet implemented";
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::solvePostProcess()
{
    Feel::cout << "solving post process" << std::endl;;
    auto pps = product( M_Whp );
    auto PP = pps.element();
    PP.buildVector(this->backend());
    auto App = std::dynamic_pointer_cast<typename super_type::backend_type::sparse_matrix_type>(M_App);
    auto Fpp = std::dynamic_pointer_cast<typename super_type::backend_type::vector_type>(M_Fpp);
    DataUpdateLinear dataPost(PP.vectorMonolithic(), App, Fpp, !M_postMatrixInit);
    this->updatePostPDE(dataPost);

    auto bbf = blockform2( pps, M_App );
    auto blf = blockform1( pps, M_Fpp );
    bbf.solve( _solution=PP, _rhs=blf, _name="sc.post", _local=true);

    M_ppp = std::make_shared<element_postpotential_type>( PP(0_c) );
    M_ppp->plusAssign(M_ppp->ewiseMean(M_P0dh), -1);
    M_ppp->plusAssign(M_pp->ewiseMean(M_P0dh), 1);
    // M_ppp -= M_ppp->ewiseMean(M_P0dh);
    // M_ppp += M_pp->ewiseMean(M_P0dh);
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    this->log("MixedPoisson","startTimeStep", "start");

    // // some time stepping require to compute residual without time derivative
    // this->updateTimeStepCurrentResidual();

    if ( M_timeStepping == "Newmark" )
    {
        // start time step
        if ( !this->doRestart() )
            M_newmarkPotential->start(this->fieldPotential());
    }
    else if ( M_timeStepping == "BDF" || M_timeStepping == "Theta" )
    {
        // start time step
        if (!this->doRestart())
            M_bdfPotential->start( M_bdfPotential->unknowns() );
    }
    // up current time
    this->updateTime( M_bdfPotential->time() );

    // update all expressions in bc or in house prec
    this->updateParameterValues();

    this->log("MixedPoisson","startTimeStep", "finish");
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->log("MixedPoisson","updateTimeStep", "start");
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    // // some time stepping require to compute residual without time derivative
    // this->updateTimeStepCurrentResidual();

    bool rebuildCstAssembly = false;
    if ( M_timeStepping == "BDF" )
    {
        int previousTimeOrder = this->timeStepBdfPotential()->timeOrder();
        M_bdfPotential->next( this->fieldPotential() );
        int currentTimeOrder = this->timeStepBdfPotential()->timeOrder();
        rebuildCstAssembly = previousTimeOrder != currentTimeOrder && this->timeStepBase()->strategy() == TS_STRATEGY_DT_CONSTANT;
    }
    else if ( M_timeStepping == "Theta" )
    {
        M_bdfPotential->next( this->fieldPotential() );
    }
    else if ( M_timeStepping == "Newmark" )
    {
        M_newmarkPotential->next( this->fieldPotential() );
    }
    // up current time
    this->updateTime( this->timeStepBase()->time() );

    if ( rebuildCstAssembly )
        this->setNeedToRebuildCstPart(true);

    this->updateParameterValues();

    this->timerTool("TimeStepping").stop("updateTimeStep");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("MixedPoisson","updateTimeStep", "finish");
}

} // namespace FeelModels
} // namespace Feel
