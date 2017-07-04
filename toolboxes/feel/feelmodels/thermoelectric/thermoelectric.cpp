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
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    std::string const& rootRepository )
    :
    super_type( prefix, worldComm, subPrefix, rootRepository ),
    M_electricProperties( new electricproperties_type( prefixvm(prefix,"electric") ) )
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

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"ThermoElectric.info") );
    //-----------------------------------------------------------------------------//
    // load config (from json)
    this->loadConfigBCFile();
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
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerRobinBC();

    this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC("elimination", marker(d) );
    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "Neumann" );
    for( auto const& d : this->M_bcNeumann )
        this->addMarkerNeumannBC(NeumannBCShape::SCALAR,marker(d));

    this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "electric-potential", "Robin" );
    for( auto const& d : this->M_bcRobin )
        this->addMarkerRobinBC( marker(d) );

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "VolumicForces" );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::loadConfigMeshFile(std::string const& geofilename)
{
    CHECK( false ) << "not allow";
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{

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
    int nBlock = 1;
    if ( M_modelName == "ThermoElectric" )
        nBlock += M_thermodynModel->nBlockMatrixGraph();
    return nBlock;
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    int indexBlock=0;
    if ( M_modelName == "ThermoElectric" )
    {
        int nBlockThermoDyn = M_thermodynModel->nBlockMatrixGraph();
        for (int tk1=0;tk1<nBlockThermoDyn ;++tk1 )
            for (int tk2=0;tk2<nBlockThermoDyn ;++tk2 )
                myblockGraph(indexBlock+tk1,indexBlock+tk2) = M_thermodynModel->buildBlockMatrixGraph()(tk1,tk2);
        indexBlock += nBlockThermoDyn;

        BlocksStencilPattern patCoupling1(1,nBlockThermoDyn,size_type(Pattern::ZERO));
        patCoupling1(0,0) = size_type(Pattern::COUPLED);
        myblockGraph(0,indexBlock) = stencil(_test=M_thermodynModel->spaceTemperature(),
                                             _trial=this->spaceElectricPotential(),
                                             _pattern_block=patCoupling1,
                                             _diag_is_nonzero=false,_close=false)->graph();
    }
    myblockGraph(indexBlock,indexBlock) = stencil(_test=this->spaceElectricPotential(),
                                                  _trial=this->spaceElectricPotential() )->graph();
    ++indexBlock;
    myblockGraph.close();

    return myblockGraph;
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("ThermoElectric","init", "start" );
    this->timerTool("Constructor").start();

    //std::string theThermoElectricModel = this->modelProperties().model();
    M_modelName = this->modelProperties().model();
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
    CHECK( ( M_modelName == "Electric" ) || ( M_modelName == "ThermoElectric" ) ) << "invalid model name : " << M_modelName << "\n";

    bool doSolveOnlyElectricPotential =
        ( ( M_modelName == "ThermoElectric" ) && ( M_solverName == "Linear" ) ) ||
        ( M_modelName == "Electric" );

    if ( M_modelName == "ThermoElectric" )
    {
        M_thermodynModel.reset( new thermodyn_model_type(prefixvm(this->prefix(),"thermo"), false, this->worldComm(),
                                                         this->subPrefix(), this->rootRepositoryWithoutNumProc() ) );
        M_thermodynModel->loadMesh( this->mesh() );
        // disable thermo exporter if we use fluid exporter
        M_thermodynModel->setDoExportResults( false );
        M_thermodynModel->init( doSolveOnlyElectricPotential /*false*/ );
        M_thermodynModel->setRowStartInMatrix( 0 );
        M_thermodynModel->setColStartInMatrix( 0 );
        M_thermodynModel->setRowStartInVector( 0 );
        if ( M_solverName == "Linear" )
        {
            M_thermodynModel->algebraicFactory()->addFunctionLinearPreAssemblyNonCst = boost::bind( &self_type::updateLinearPreAssemblyJouleLaw,
                                                                                                    boost::ref( *this ), _1, _2 );
        }
    }

    // physical properties
    M_electricProperties->updateForUse( M_mesh, this->modelProperties().materials(),  this->localNonCompositeWorldsComm());

    // functionspace
    if ( M_electricProperties->isDefinedOnWholeMesh() )
    {
        M_rangeMeshElements = elements(M_mesh);
        M_XhElectricPotential = space_electricpotential_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
        M_XhElectricField = space_electricfield_type::New(_mesh=M_mesh, _worldscomm=this->worldsComm() );
    }
    else
    {
        M_rangeMeshElements = markedelements(M_mesh, M_electricProperties->markers());
        M_XhElectricPotential = space_electricpotential_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
        M_XhElectricField = space_electricfield_type::New(_mesh=M_mesh, _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
    }
    M_fieldElectricPotential.reset( new element_electricpotential_type(M_XhElectricPotential,"V"));
    M_fieldElectricField.reset( new element_electricfield_type(M_XhElectricField,"E"));

    // backend : use worldComm of Xh
    M_backendMonolithic = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldComm() );

    size_type currentStartIndex = 0;// velocity and pressure before
    if ( M_modelName == "ThermoElectric" )
    {
        M_startBlockIndexFieldsInMatrix["temperature"] = currentStartIndex;
        currentStartIndex += M_thermodynModel->nBlockMatrixGraph();
    }
    M_startBlockIndexFieldsInMatrix["potential-electric"] = currentStartIndex;

    // vector solution
    int nBlock = this->nBlockMatrixGraph();
    M_blockVectorSolutionMonolithic.resize( nBlock );
    int indexBlock=0;
    if ( M_modelName == "ThermoElectric" )
    {
        int nBlockThermoDyn = M_thermodynModel->nBlockMatrixGraph();
        for (int tk1=0;tk1<nBlockThermoDyn ;++tk1 )
            M_blockVectorSolutionMonolithic(indexBlock+tk1) = M_thermodynModel->blockVectorSolution()(tk1);
        indexBlock+=nBlockThermoDyn;
    }
    M_blockVectorSolutionMonolithic(indexBlock) = this->fieldElectricPotentialPtr();

    // init petsc vector associated to the block
    M_blockVectorSolutionMonolithic.buildVector( this->backend() );


    // start or restart time step scheme and exporter
#if 0
    if (!this->doRestart())
    {
        // start time step
        M_bdfTemperature->start(this->fieldTemperature());
        // up current time
        this->updateTime( M_bdfTemperature->time() );
    }
    else
    {
        // start time step
        M_bdfTemperature->restart();
        // load a previous solution as current solution
        *this->fieldTemperaturePtr() = M_bdfTemperature->unknown(0);
        // up initial time
        this->setTimeInitial( M_bdfTemperature->timeInitial() );
        // up current time
        this->updateTime( M_bdfTemperature->time() );
    }
#endif
    this->updateBoundaryConditionsForUse();

    // post-process
    this->initPostProcess();

    // algebraic solver
    if ( buildModelAlgebraicFactory )
    {
        if ( !doSolveOnlyElectricPotential )//monolithic
        {
            // matrix graph of non zero
            typename model_algebraic_factory_type::graph_ptrtype graph( new typename model_algebraic_factory_type::graph_type( this->buildBlockMatrixGraph() ) );
            // tool for assembly and solver
            M_algebraicFactoryMonolithic.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend(),
                                                                                  graph, graph->mapRow().indexSplit() ) );
        }
        if ( doSolveOnlyElectricPotential ) // only
        {
            M_backendElectricModel = backend_type::build( soption( _name="backend" ), prefixvm(this->prefix(),"electric"), this->worldComm() );
            auto graph = stencil(_test=this->spaceElectricPotential(),
                                 _trial=this->spaceElectricPotential() )->graph();
            // tool for assembly and solver
            M_algebraicFactoryElectricModel.reset( new model_algebraic_factory_type( this->shared_from_this(),M_backendElectricModel,
                                                                                     graph, graph->mapRow().indexSplit() ) );
        }
    }

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("ThermoElectric","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateBoundaryConditionsForUse()
{
    auto mesh = this->mesh();
    auto XhElectricPotential = this->spaceElectricPotential();

    auto & dofsWithValueImposedElectricPotential = M_dofsWithValueImposed["electric-potential"];
    dofsWithValueImposedElectricPotential.clear();
    std::set<std::string> electricPotentialMarkers;

    // strong Dirichlet bc on temperature from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = this->markerDirichletBCByNameId( "elimination",marker(d) );
        electricPotentialMarkers.insert( listMark.begin(), listMark.end() );
    }
    auto meshMarkersElectricPotentialByEntities = detail::distributeMarkerListOnSubEntity( mesh, electricPotentialMarkers );

    // on topological faces
    auto const& listMarkedFacesElectricPotential = std::get<0>( meshMarkersElectricPotentialByEntities );
    for ( auto const& faceWrap : markedfaces(mesh,listMarkedFacesElectricPotential ) )
    {
        auto const& face = unwrap_ref( faceWrap );
        auto facedof = XhElectricPotential->dof()->faceLocalDof( face.id() );
        for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
            dofsWithValueImposedElectricPotential.insert( it->index() );
    }
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("ThermoElectric","initPostProcess", "start");
    this->timerTool("Constructor").start();

    M_postProcessFieldExported.clear();
    if ( this->modelProperties().postProcess().find("Fields") != this->modelProperties().postProcess().end() )
        for ( auto const& o : this->modelProperties().postProcess().find("Fields")->second )
        {
            if ( M_modelName == "ThermoElectric" )
                if ( o == "temperature" || o == "all" ) M_postProcessFieldExported.insert( ThermoElectricPostProcessFieldExported::Temperature );
            if ( o == "electric-potential" || o == "all" ) M_postProcessFieldExported.insert( ThermoElectricPostProcessFieldExported::ElectricPotential );
            if ( o == "electric-field" || o == "all" ) M_postProcessFieldExported.insert( ThermoElectricPostProcessFieldExported::ElectricField );
            if ( o == "pid" || o == "all" ) M_postProcessFieldExported.insert( ThermoElectricPostProcessFieldExported::Pid );
        }


    std::string geoExportType="static";//change_coords_only, change, static
    M_exporter = exporter( _mesh=this->mesh(),
                           _name="Export",
                           _geo=geoExportType,
                           _path=this->exporterPath() );

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("ThermoElectric","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::restartExporters()
{
    if (this->doRestart() && this->restartPath().empty() )
    {
        if ( M_exporter->doExport() )
            M_exporter->restart(this->timeInitial());
    }
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
boost::shared_ptr<std::ostringstream>
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    if ( M_modelName == "ThermoElectric" )
        *_ostr << M_thermodynModel->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||-----------Info : Electric--------------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << prefixvm(this->prefix(),"electric")
           << "\n   Root Repository : " << this->rootRepository();
    *_ostr << "\n   Physical Model"
           << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient");
    *_ostr << "\n   Boundary conditions"
           << this->getInfoDirichletBC()
           << this->getInfoNeumannBC()
           << this->getInfoRobinBC();
    *_ostr << M_electricProperties->getInfoMaterialParameters()->str();
    *_ostr << "\n   Mesh Discretization"
           << "\n     -- msh filename      : " << this->mshfileStr()
           << "\n     -- number of element : " << M_mesh->numGlobalElements()
           << "\n     -- order             : " << nOrderGeo;
    *_ostr << "\n   Space ElectricPotential Discretization"
           << "\n     -- order         : " << nOrderPolyElectricPotential
           << "\n     -- number of dof : " << M_XhElectricPotential->nDof() << " (" << M_XhElectricPotential->nLocalDof() << ")";
    if ( M_algebraicFactoryElectricModel )
        *_ostr << M_algebraicFactoryElectricModel->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";



    if ( M_modelName == "ThermoElectric" )
    {
        std::string myexporterType = M_exporter->type();
        int myexporterFreq = M_exporter->freq();
        std::string doExport_str;
        if ( this->hasPostProcessFieldExported( ThermoElectricPostProcessFieldExported::Temperature ) )
            doExport_str=(doExport_str.empty())?"temperature":doExport_str+" - temperature";
        if ( this->hasPostProcessFieldExported( ThermoElectricPostProcessFieldExported::ElectricPotential ) )
            doExport_str=(doExport_str.empty())?"electric-potential":doExport_str+" - electric-potential";
        if ( this->hasPostProcessFieldExported( ThermoElectricPostProcessFieldExported::ElectricField ) )
            doExport_str=(doExport_str.empty())?"electric-field":doExport_str+" - electric-field";
        if ( this->hasPostProcessFieldExported( ThermoElectricPostProcessFieldExported::Pid ) )
            doExport_str=(doExport_str.empty())?"pid":doExport_str+" - pid";


        *_ostr << "\n||==============================================||"
               << "\n||==============================================||"
               << "\n||==============================================||"
               << "\n||----------Info : ThermoElectric---------------||"
               << "\n||==============================================||"
               << "\n||==============================================||"
               << "\n||==============================================||"
               << "\n   Prefix : " << this->prefix()
               << "\n   Root Repository : " << this->rootRepository();
        *_ostr << "\n   Physical Model"
               << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient");
        *_ostr << "\n   Exporter"
               << "\n     -- type            : " << myexporterType
               << "\n     -- freq save       : " << myexporterFreq
               << "\n     -- fields exported : " << doExport_str
               << "\n   Processors"
               << "\n     -- number of proc : " << this->worldComm().globalSize()
               << "\n     -- current rank : " << this->worldComm().globalRank();

        if ( M_algebraicFactoryMonolithic )
            *_ostr << M_algebraicFactoryMonolithic->getInfo()->str();
        *_ostr << "\n||==============================================||"
               << "\n||==============================================||"
               << "\n||==============================================||"
               << "\n";
    }
    return _ostr;
}



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    if ( !M_exporter->doExport() ) return;

    this->log("ThermoElectric","exportResults", "start");
    this->timerTool("PostProcessing").start();

    bool hasFieldToExport = false;

    if ( this->hasPostProcessFieldExported( ThermoElectricPostProcessFieldExported::Temperature ) &&
         M_thermodynModel->mesh()->isSameMesh( this->mesh() ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"temperature"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"temperature")),
                                       M_thermodynModel->fieldTemperature() );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( ThermoElectricPostProcessFieldExported::ElectricPotential ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"electric-potential"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"electric-potential")),
                                       this->fieldElectricPotential() );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( ThermoElectricPostProcessFieldExported::ElectricField ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"electric-fields"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"electric-fields")),
                                       *M_fieldElectricField );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( ThermoElectricPostProcessFieldExported::Pid ) )
    {
        M_exporter->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
        hasFieldToExport = true;
    }

    if ( hasFieldToExport )
        M_exporter->save();

    // this->exportMeasures( time );

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
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateElectricField()
{
    std::string M_computeElectricFieldProjType = "nodal";
    /*if ( M_computeElectricFieldProjType == "L2" )
     *M_fieldElectricFieldContinuous = M_l2proj->operator()( -trans(gradv(this->fieldElectricPotential() ) ) );
     else*/ if ( M_computeElectricFieldProjType == "nodal" )
        M_fieldElectricField->on(_range=M_rangeMeshElements,
                                 _expr=-trans(gradv( this->fieldElectricPotential() ) ) );
    else
        CHECK( false ) << "invalid M_computeElectricFieldProjType " << M_computeElectricFieldProjType << "\n";
}



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    this->modelProperties().parameters().updateParameterValues();

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    M_bcDirichlet.setParameterValues( paramValues );
    M_bcNeumann.setParameterValues( paramValues );
    M_bcRobin.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("ThermoElectric","solve", "start");
    this->timerTool("Solve").start();

    if ( M_modelName == "ThermoElectric" )
        M_thermodynModel->updateParameterValues();
    this->updateParameterValues();

    //std::string theThermoElectricModel = this->modelProperties().model();
    //if ( theThermoElectricModel == "ThermoElectric-linear" )
    if ( M_solverName == "Linear" )
    {
        auto mySolutionVectorLinear = M_backendElectricModel->toBackendVectorPtr( this->fieldElectricPotential() );
        M_algebraicFactoryElectricModel->solve( "LinearSystem", mySolutionVectorLinear );
        if ( M_modelName == "ThermoElectric" )
            M_thermodynModel->solve();
    }
    else if ( M_solverName == "Newton" )
    {
        M_blockVectorSolutionMonolithic.updateVectorFromSubVectors();
        M_algebraicFactoryMonolithic->solve( "Newton", M_blockVectorSolutionMonolithic.vector() );
        M_blockVectorSolutionMonolithic.localize();
    }

    this->updateElectricField();

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
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("ThermoElectric","updateLinearPDE", "start"+sc);
    boost::mpi::timer thetimer;

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    auto const& v = this->fieldElectricPotential();

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix() ,
                                              _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=XhV, _vector=F,
                               _rowstart=this->rowStartInVector() );

    //--------------------------------------------------------------------------------------------------//

    auto sigma = idv(M_electricProperties->fieldElectricConductivity());
    if ( buildCstPart )
    {
        bilinearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= sigma*inner(gradt(v),grad(v)),
                       _geomap=this->geomap() );
    }

    // update source term
    this->updateSourceTermLinearPDE(F, buildCstPart);

    // update bc
    this->updateWeakBCLinearPDE(A,F,buildCstPart);

    if ( !buildCstPart && _doBCStrongDirichlet)
    {
        this->updateBCStrongDirichletLinearPDE(A,F);
    }
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearPreAssemblyJouleLaw( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const
{
    this->log("ThermoElectric","updateLinearPreAssemblyJouleLaw","start" );

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    auto const& v = this->fieldElectricPotential();
    auto XhT = M_thermodynModel->spaceTemperature();
    auto const& t = M_thermodynModel->fieldTemperature();
    auto sigma = idv(M_electricProperties->fieldElectricConductivity());
    form1( _test=XhT,_vector=F,
           _rowstart=M_thermodynModel->rowStartInVector() ) +=
        integrate( _range=M_rangeMeshElements,
                       _expr= sigma*inner(gradv(v),gradv(v))*id(t),
                       _geomap=this->geomap() );

    this->log("ThermoElectric","updateLinearPreAssemblyJouleLaw","finish" );
}



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
    M_thermodynModel->updateNewtonInitialGuess( U );

    if ( M_bcDirichlet.empty() ) return;

    this->log("ThermoElectric","updateNewtonInitialGuess","start" );

    auto mesh = this->mesh();
    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;
    auto v = this->spaceElectricPotential()->element( U, this->rowStartInVector()+startBlockIndexElectricPotential );
    for( auto const& d : M_bcDirichlet )
    {
        v.on(_range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
             _expr=expression(d) );
    }
    // synchronize electric potential dof on interprocess
    auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("electric-potential");
    if ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )
        sync( v, "=", itFindDofsWithValueImposed->second );

    this->log("ThermoElectric","updateNewtonInitialGuess","finish" );
}
THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool _BuildCstPart = data.buildCstPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("ThermoElectric","updateJacobian", "start"+sc);
    size_type startBlockIndexTemperature = M_startBlockIndexFieldsInMatrix.find( "temperature" )->second;
    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    // auto const& v = this->fieldElectricPotential();
    auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );
    auto XhT = M_thermodynModel->spaceTemperature();
    auto const& t = M_thermodynModel->fieldTemperature();

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                                              _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );

    auto sigma = idv(M_electricProperties->fieldElectricConductivity());
    if ( buildCstPart )
    {
        bilinearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= sigma*inner(gradt(v),grad(v)),
                       _geomap=this->geomap() );
    }
    if ( !buildCstPart )
    {
        form2( _test=XhT,_trial=XhV,_matrix=J,
               _pattern=size_type(Pattern::COUPLED),
               _rowstart=this->rowStartInMatrix()+startBlockIndexTemperature,
               _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential  ) +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -sigma*2*inner(gradt(v),gradv(v))*id( t ),
                       _geomap=this->geomap() );
    }

    DataUpdateJacobian dataThermo( data );
    dataThermo.setDoBCStrongDirichlet( false );
    M_thermodynModel->updateJacobian( dataThermo );

    this->updateBCWeakJacobian( v,J,buildCstPart );

    if ( buildNonCstPart && _doBCStrongDirichlet )
    {
        this->updateBCStrongDirichletJacobian( J,RBis );
        M_thermodynModel->updateBCStrongDirichletJacobian( J,RBis );
    }
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletJacobian(sparse_matrix_ptrtype& J,vector_ptrtype& RBis ) const
{
    if ( this->M_bcDirichlet.empty() ) return;

    this->log("ThermoElectric","updateBCStrongDirichletJacobian","start" );

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    auto const& v = this->fieldElectricPotential();
    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;
    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                                              _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );

    for( auto const& d : this->M_bcDirichlet )
    {
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                _element=v,_rhs=RBis,_expr=cst(0.) );
    }

    this->log("ThermoElectric","updateBCStrongDirichletJacobian","finish" );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateBCWeakJacobian( element_electricpotential_external_storage_type const& v, sparse_matrix_ptrtype& J, bool buildCstPart ) const
{
    if ( this->M_bcRobin.empty() ) return;

    if ( !buildCstPart )
    {
        auto XhV = this->spaceElectricPotential();
        auto mesh = XhV->mesh();
        size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;

        auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                                                  _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );
        for( auto const& d : this->M_bcRobin )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( marker(d) ) ),
                           _expr= expression1(d)*idt(v)*id(v),
                           _geomap=this->geomap() );
        }
    }
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool _BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("ThermoElectric","updateResidual", "start"+sc);

    size_type startBlockIndexTemperature = M_startBlockIndexFieldsInMatrix.find( "temperature" )->second;
    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    // auto const& v = this->fieldElectricPotential();
    auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );

    auto XhT = M_thermodynModel->spaceTemperature();
    // auto const& t = M_thermodynModel->fieldTemperature();
    auto const t = XhT->element(XVec, this->rowStartInVector()+startBlockIndexTemperature );

    auto myLinearFormThermo = form1( _test=XhV, _vector=R,
                                     _rowstart=this->rowStartInVector() );
    auto myLinearFormElectric = form1( _test=XhV, _vector=R,
                               _rowstart=this->rowStartInVector() + startBlockIndexElectricPotential );

    auto sigma = idv(M_electricProperties->fieldElectricConductivity());
    if (!buildCstPart && !UseJacobianLinearTerms )
    {
        myLinearFormElectric +=
            integrate( _range=M_rangeMeshElements,
                       _expr= sigma*inner(gradv(v),grad(v)),
                       _geomap=this->geomap() );
    }
    if ( !buildCstPart )
    {
        myLinearFormThermo +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -sigma*inner(gradv(v),gradv(v))*id( t ),
                       _geomap=this->geomap() );
    }

    DataUpdateResidual dataThermo( data );
    dataThermo.setDoBCStrongDirichlet( false );
    M_thermodynModel->updateResidual( dataThermo );

    this->updateSourceTermResidual( R,buildCstPart ) ;

    this->updateBCWeakResidual( v,R,buildCstPart );

    if ( !buildCstPart && _doBCStrongDirichlet &&
         (this->hasMarkerDirichletBCelimination() || M_thermodynModel->hasMarkerDirichletBCelimination() ) )
    {
        R->close();
        this->updateBCDirichletStrongResidual( R );
        M_thermodynModel->updateBCDirichletStrongResidual( R );
    }
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateBCDirichletStrongResidual( vector_ptrtype& R ) const
{
    if ( this->M_bcDirichlet.empty() ) return;

    this->log("ThermoElectric","updateBCDirichletStrongResidual","start" );

    auto XhV = this->spaceElectricPotential();
    auto mesh = XhV->mesh();
    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;
    auto v = this->spaceElectricPotential()->element( R,this->rowStartInVector()+startBlockIndexElectricPotential );
    auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("electric-potential");
    auto const& dofsWithValueImposedElectricPotential = ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )? itFindDofsWithValueImposed->second : std::set<size_type>();
    for ( size_type thedof : dofsWithValueImposedElectricPotential )
        v.set( thedof,0. );
    sync( v, "=", dofsWithValueImposedElectricPotential );

    this->log("ThermoElectric","updateBCDirichletStrongResidual","finish" );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateBCWeakResidual( element_electricpotential_external_storage_type const& v, vector_ptrtype& R, bool buildCstPart ) const
{
    if ( this->M_bcNeumann.empty() && this->M_bcRobin.empty() ) return;

    auto XhV = this->spaceElectricPotential();
    auto mesh = XhV->mesh();
    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;

    auto myLinearForm = form1( _test=XhV, _vector=R,
                               _rowstart=this->rowStartInVector()+startBlockIndexElectricPotential );
    if ( buildCstPart )
    {
        for( auto const& d : this->M_bcNeumann )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerNeumannBC(NeumannBCShape::SCALAR,marker(d)) ),
                           _expr= -expression(d)*id(v),
                           _geomap=this->geomap() );
        }
    }
    for( auto const& d : this->M_bcRobin )
    {
        if ( !buildCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( marker(d) ) ),
                           _expr= expression1(d)*idv(v)*id(v),
                           _geomap=this->geomap() );
        }
        if ( buildCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( marker(d) ) ),
                           _expr= -expression1(d)*expression2(d)*id(v),
                           _geomap=this->geomap() );
        }
    }
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateSourceTermResidual( vector_ptrtype& R, bool buildCstPart ) const
{
    if ( this->M_volumicForcesProperties.empty() ) return;

    if ( buildCstPart )
    {
        auto XhV = this->spaceElectricPotential();
        auto mesh = XhV->mesh();
        size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;
        auto myLinearForm = form1( _test=XhV, _vector=R,
                                   _rowstart=this->rowStartInVector()+startBlockIndexElectricPotential );
        auto const& v = this->fieldElectricPotential();

        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto rangeEltUsed = (marker(d).empty())? M_rangeMeshElements : markedelements(this->mesh(),marker(d));
            myLinearForm +=
                integrate( _range=rangeEltUsed,
                           _expr= -expression(d)*id(v),
                           _geomap=this->geomap() );
        }
    }
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    if ( this->M_bcDirichlet.empty() ) return;

    this->log("ThermoElectric","updateBCStrongDirichletLinearPDE","start" );

    auto XhV = this->spaceElectricPotential();
    auto const& v = this->fieldElectricPotential();
    auto mesh = XhV->mesh();

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );
    for( auto const& d : this->M_bcDirichlet )
    {
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                _element=v,_rhs=F,_expr=expression(d) );
    }

    this->log("ThermoElectric","updateBCStrongDirichletLinearPDE","finish" );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateSourceTermLinearPDE( vector_ptrtype& F, bool buildCstPart ) const
{
#if 0
    if ( this->M_overwritemethod_updateSourceTermLinearPDE != NULL )
    {
        this->M_overwritemethod_updateSourceTermLinearPDE(F,buildCstPart);
        return;
    }
#endif

    if ( this->M_volumicForcesProperties.empty() ) return;

    if ( !buildCstPart )
    {
        auto XhV = this->spaceElectricPotential();
        auto const& v = this->fieldElectricPotential();
        auto mesh = XhV->mesh();

        auto myLinearForm = form1( _test=XhV, _vector=F,
                                   _rowstart=this->rowStartInVector() );

        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto rangeEltUsed = (marker(d).empty())? M_rangeMeshElements : markedelements(this->mesh(),marker(d));
            myLinearForm +=
                integrate( _range=rangeEltUsed,
                           _expr= expression(d)*id(v),
                           _geomap=this->geomap() );
        }
    }

}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const
{
    if ( this->M_bcNeumann.empty() && this->M_bcRobin.empty() ) return;

    if ( !buildCstPart )
    {
        auto XhV = this->spaceElectricPotential();
        auto const& v = this->fieldElectricPotential();
        auto mesh = XhV->mesh();

        auto myLinearForm = form1( _test=XhV, _vector=F,
                                   _rowstart=this->rowStartInVector() );
        for( auto const& d : this->M_bcNeumann )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerNeumannBC(NeumannBCShape::SCALAR,marker(d)) ),
                           _expr= expression(d)*id(v),
                           _geomap=this->geomap() );
        }

        auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=this->rowStartInMatrix(),
                                                  _colstart=this->colStartInMatrix() );
        for( auto const& d : this->M_bcRobin )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( marker(d) ) ),
                           _expr= expression1(d)*idt(v)*id(v),
                           _geomap=this->geomap() );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( marker(d) ) ),
                           _expr= expression1(d)*expression2(d)*id(v),
                           _geomap=this->geomap() );
        }

    }
}


} // end namespace FeelModels
} // end namespace Feel
