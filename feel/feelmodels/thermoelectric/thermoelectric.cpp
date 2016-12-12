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
    // load info from .bc file
    this->loadConfigBCFile();
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // build mesh, space, exporter,...
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

    this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC("elimination", marker(d) );
    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "Neumann" );
    for( auto const& d : this->M_bcNeumann )
        this->addMarkerNeumannBC(NeumannBCShape::SCALAR,marker(d));

    this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "electric-potential", "Robin" );
    for( auto const& d : this->M_bcRobin )
        this->addMarkerRobinBC( marker(d) );

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "VolumicForces" );
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
    //TODO
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
    return M_thermodynModel->nBlockMatrixGraph() + 1;
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    int indexBlock=0;
    int nBlockThermoDyn = M_thermodynModel->nBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockThermoDyn ;++tk1 )
        for (int tk2=0;tk2<nBlockThermoDyn ;++tk2 )
            myblockGraph(indexBlock+tk1,indexBlock+tk2) = M_thermodynModel->buildBlockMatrixGraph()(tk1,tk2);
    indexBlock += nBlockThermoDyn;

    myblockGraph(indexBlock,indexBlock) = stencil(_test=this->spaceElectricPotential(),
                                                  _trial=this->spaceElectricPotential() )->graph();

#if 0
    BlocksStencilPattern patCoupling1(1,space_fluid_type::nSpaces,size_type(Pattern::ZERO));
    patCoupling1(0,0) = size_type(Pattern::COUPLED);
    myblockGraph(indexBlock+1,0) = stencil(_test=M_thermodynModel->spaceTemperature(),
                                           _trial=this->functionSpace(),
                                           _pattern_block=patCoupling1,
                                           _diag_is_nonzero=false,_close=false)->graph();

    BlocksStencilPattern patCoupling2(space_fluid_type::nSpaces,1,size_type(Pattern::ZERO));
    patCoupling2(0,0) = size_type(Pattern::COUPLED);
    myblockGraph(0,indexBlock+1) = stencil(_test=this->functionSpace(),
                                           _trial=M_thermodynModel->spaceTemperature(),
                                           _pattern_block=patCoupling2,
                                           _diag_is_nonzero=false,_close=false)->graph();
    ++indexBlock;
#endif

    myblockGraph.close();

    return myblockGraph;
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("ThermoElectric","init", "start" );
    this->timerTool("Constructor").start();

    std::string theThermoElectricModel = this->modelProperties().model();
    bool doSolveOnlyElectricPotential = ( theThermoElectricModel == "ThermoElectric-linear" );

    
    M_thermodynModel.reset( new thermodyn_model_type(prefixvm(this->prefix(),"thermo"), false, this->worldComm(),
                                                     this->subPrefix(), this->rootRepositoryWithoutNumProc() ) );
    M_thermodynModel->loadMesh( this->mesh() );
    // disable thermo exporter if we use fluid exporter
    M_thermodynModel->setDoExportResults( false );
    M_thermodynModel->init( doSolveOnlyElectricPotential /*false*/ );
    M_thermodynModel->setRowStartInMatrix( 0 );
    M_thermodynModel->setColStartInMatrix( 0 );
    M_thermodynModel->setRowStartInVector( 0 );


    // functionspace
    M_XhElectricPotential = space_electricpotential_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    M_fieldElectricPotential.reset( new element_electricpotential_type(M_XhElectricPotential,"V"));

    // physical properties
    M_XhScalarP0 = space_scalar_P0_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    M_electricProperties->initFromSpace( M_XhScalarP0 );
    M_electricProperties->updateFromModelMaterials( this->modelProperties().materials() );


    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldComm() );

    size_type currentStartIndex = 0;// velocity and pressure before
    M_startBlockIndexFieldsInMatrix["temperature"] = 0;
    M_startBlockIndexFieldsInMatrix["potential-electric"] = M_thermodynModel->nBlockMatrixGraph();

    // vector solution
    int nBlock = this->nBlockMatrixGraph();
    M_blockVectorSolution.resize( nBlock );
    int nBlockThermoDyn = M_thermodynModel->nBlockMatrixGraph();
    int indexBlock=0;
    for (int tk1=0;tk1<nBlockThermoDyn ;++tk1 )
        M_blockVectorSolution(indexBlock+tk1) = M_thermodynModel->blockVectorSolution()(tk1);
    indexBlock+=nBlockThermoDyn;
    M_blockVectorSolution(indexBlock) = this->fieldElectricPotentialPtr();

    // init petsc vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );


    
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
    // post-process
    this->initPostProcess();

    // algebraic solver
    if ( buildModelAlgebraicFactory )
    {
        if ( false )//monolithic
        {
            // matrix graph of non zero
            typename model_algebraic_factory_type::graph_ptrtype graph( new typename model_algebraic_factory_type::graph_type( this->buildBlockMatrixGraph() ) );

            // tool for assembly and solver
            M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend(),
                                                                        graph, graph->mapRow().indexSplit() ) );
        }
        if ( true ) // only
        {
            // auto M_backend2 = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldComm() );
            // matrix graph of non zero
            auto graph = this->buildBlockMatrixGraph()(1,1);
            // auto graph = stencil(_test=this->spaceElectricPotential(),
            //                      _trial=this->spaceElectricPotential() )->graph();
            // tool for assembly and solver
            M_algebraicFactory.reset( new model_algebraic_factory_type(this->shared_from_this(),this->backend(),
                                                                       graph, graph->mapRow().indexSplit() ) );
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
    if ( M_algebraicFactory )
        *_ostr << M_algebraicFactory->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    *_ostr << M_thermodynModel->getInfo()->str();

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

    if ( true )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"electric-potential"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"electric-potential")),
                                       this->fieldElectricPotential() );
        hasFieldToExport = true;
    }
    if ( M_thermodynModel->mesh()->isSameMesh( this->mesh() ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"temperature"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"temperature")),
                                       M_thermodynModel->fieldTemperature() );
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
    
    M_thermodynModel->updateParameterValues();
    this->updateParameterValues();

    std::string theThermoElectricModel = this->modelProperties().model();
    if ( theThermoElectricModel == "ThermoElectric-linear" )
    {
        // auto mySolutionVectorLinear = M_backend()->newVector( this->spaceElectricPotential()
        auto mySolutionVectorLinear = M_backend->toBackendVectorPtr( this->fieldElectricPotential() );
        M_algebraicFactory->solve( "LinearSystem", mySolutionVectorLinear/*this->blockVectorSolution().vector()*/ );

        M_thermodynModel->algebraicFactory()->addFunctionLinearPreAssemblyNonCst = boost::bind( &self_type::updateLinearPreAssemblyJouleLaw,
                                                                                                boost::ref( *this ), _1, _2 );
        M_thermodynModel->solve();
    }
    else
    {
        CHECK( false ) << "TODO";
        // M_thermodynModel->algebraicFactory()->addFunctionLinearPreAssemblyNonCst = nullptr;
        //M_algebraicFactory->linearSolver(this->blockVectorSolution().vector());
        M_algebraicFactory->solve( "LinearSystem", this->blockVectorSolution().vector() );
        //M_algebraicFactory->solve( "Newton", this->blockVectorSolution().vector() );

        M_blockVectorSolution.localize();
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
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("ThermoElectric","updateLinearPDE", "start"+sc);
    boost::mpi::timer thetimer;

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    auto XhT = this->thermodynModel()->spaceTemperature();
    auto const& t = this->thermodynModel()->fieldTemperature();
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
            integrate( _range=elements(mesh),
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
    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    auto const& v = this->fieldElectricPotential();
    auto XhT = M_thermodynModel->spaceTemperature();
    auto const& t = M_thermodynModel->fieldTemperature();

    auto sigma = idv(M_electricProperties->fieldElectricConductivity());
    form1( _test=XhT,_vector=F,
           _rowstart=M_thermodynModel->rowStartInMatrix() ) +=
        integrate( _range=elements(mesh),
                       _expr= sigma*inner(gradv(v),gradv(v))*id(t),
                       _geomap=this->geomap() );
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

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    auto const& v = this->fieldElectricPotential();
    auto XhT = M_thermodynModel->spaceTemperature();
    auto const& t = M_thermodynModel->fieldTemperature();

    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;
    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                                              _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );

    if ( buildCstPart )
    {
        auto sigma = idv(M_electricProperties->fieldElectricConductivity());
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= sigma*inner(gradt(v),grad(v)),
                       _geomap=this->geomap() );
    }

    if ( !buildNonCstPart )
    {
#if 0
        form2( _test=XhT,_trial=XhV,_matrix=A,
               _pattern=size_type(Pattern::COUPLED),
               _rowstart=this->rowStartInMatrix(),
               _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential  ) +=
            integrate( _range=elements(mesh),
                       _expr= sigma*inner(gradt(v),gradt(v))*id( t ),
                       _geomap=this->geomap() );
#endif
    }

    DataUpdateJacobian dataThermo( data );
    dataThermo.setDoBCStrongDirichlet( false );
    M_thermodynModel->updateJacobian( dataThermo );

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
#if 0
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

#endif
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
#if 0
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool _BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("ThermoDynamics","updateResidual", "start"+sc);

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    auto const& v = this->fieldElectricPotential();
    auto XhT = M_thermodynModel->spaceTemperature();
    auto const& t = M_thermodynModel->fieldTemperature();

    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;
    auto myLinearFormThermo = form1( _test=XhV, _vector=R,
                                     _rowstart=this->rowStartInVector() );
    auto myLinearFormElectric = form1( _test=XhV, _vector=R,
                               _rowstart=this->rowStartInVector() + startBlockIndexElectricPotential );

    auto sigma = idv(M_electricProperties->fieldElectricConductivity());
    if (!buildCstPart && !UseJacobianLinearTerms )
    {
        myLinearFormElectric +=
            integrate( _range=elements(mesh),
                       _expr= sigma*inner(gradv(v),grad(v)),
                       _geomap=this->geomap() );
    }
    if ( !buildCstPart )
    {
        myLinearFormThermo +=
            integrate( _range=elements(mesh),
                       _expr= -sigma*inner(gradv(v),gradv(v))*id( t ),
                       _geomap=this->geomap() );
    }

    if ( !buildCstPart && _doBCStrongDirichlet /*&& this->hasMarkerDirichletBCelimination()*/ )
    {
        R->close();
        this->updateBCDirichletStrongResidual( R );
        M_thermodynModel->updateBCDirichletStrongResidual( R );
    }

#endif
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateBCDirichletStrongResidual( vector_ptrtype& R ) const
{
#if 0
    if ( this->M_bcDirichlet.empty() ) return;

    this->log("ThermoDynamics","updateBCDirichletStrongResidual","start" );

    auto XhV = this->spaceElectricPotential();
    auto const& v = this->fieldElectricPotential();
    auto mesh = XhV->mesh();

    auto mesh = this->mesh();
    auto u = this->spaceTemperature()->element( R,this->rowStartInVector() );

    for( auto const& d : this->M_bcDirichlet )
    {
        u.on(_range=markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
             _expr=cst(0.) );
    }

    this->log("ThermoDynamics","updateBCDirichletStrongResidual","finish" );

#endif
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
            if ( marker(d).empty() )
                myLinearForm +=
                    integrate( _range=elements(this->mesh()),
                               _expr= expression(d)*id(v),
                               _geomap=this->geomap() );
            else
                myLinearForm +=
                    integrate( _range=markedelements(this->mesh(),marker(d)),
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
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,marker(d)) ),
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
                integrate( _range=markedfaces(this->mesh(),marker(d) ),
                           _expr= expression1(d)*idt(v)*id(v),
                           _geomap=this->geomap() );
            myLinearForm +=
                integrate( _range=markedfaces(this->mesh(),marker(d) ),
                           _expr= expression1(d)*expression2(d)*id(v),
                           _geomap=this->geomap() );
        }

    }
}


} // end namespace FeelModels
} // end namespace Feel
