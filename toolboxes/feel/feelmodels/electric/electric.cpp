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
/*#include <feel/feelvf/form.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/operators.hpp>
 #include <feel/feelvf/operations.hpp>*/

#include <feel/feelmodels/modelmesh/createmesh.hpp>

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
    ModelPhysics<nDim>( "electric" )
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

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("createMesh");
    this->log("Electric","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()



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

    if ( !M_mesh )
        this->initMesh();

    // physical properties
    if ( !M_materialsProperties )
    {
        auto paramValues = this->modelProperties().parameters().toParameterValues();
        this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->prefix(), this->repository().expr() ) );
        M_materialsProperties->updateForUse( M_mesh, this->modelProperties().materials(), *this );
    }

    // functionspace
    if ( this->materialsProperties()->isDefinedOnWholeMesh( this->physic() ) )
    {
        M_rangeMeshElements = elements(M_mesh);
        M_XhElectricPotential = space_electricpotential_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    }
    else
    {
        M_rangeMeshElements = markedelements(M_mesh, this->materialsProperties()->markers( this->physic() ));
        M_XhElectricPotential = space_electricpotential_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
    }
    M_fieldElectricPotential.reset( new element_electricpotential_type(M_XhElectricPotential,"V"));

    this->initBoundaryConditions();

    this->initInitialConditions();

    // post-process
    this->initPostProcess();

    // update constant parameters
    this->updateParameterValues();

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() );

    size_type currentStartIndex = 0;// velocity and pressure before
    this->setStartSubBlockSpaceIndex( "potential-electric", currentStartIndex );

    // vector solution
    int nBlock = this->nBlockMatrixGraph();
    M_blockVectorSolution.resize( nBlock );
    int indexBlock=0;
    M_blockVectorSolution(indexBlock) = this->fieldElectricPotentialPtr();

    // init petsc vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );

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
    M_bcDirichletMarkerManagement.clearMarkerDirichletBC();
    M_bcNeumannMarkerManagement.clearMarkerNeumannBC();
    M_bcRobinMarkerManagement.clearMarkerRobinBC();

    this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
        M_bcDirichletMarkerManagement.addMarkerDirichletBC("elimination", name(d), markers(d) );
    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "Neumann" );
    for( auto const& d : this->M_bcNeumann )
        M_bcNeumannMarkerManagement.addMarkerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d),markers(d));

    this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "electric-potential", "Robin" );
    for( auto const& d : this->M_bcRobin )
         M_bcRobinMarkerManagement.addMarkerRobinBC( name(d),markers(d) );

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "VolumicForces" );

    auto mesh = this->mesh();
    auto XhElectricPotential = this->spaceElectricPotential();

    std::set<std::string> electricPotentialMarkers;

    // strong Dirichlet bc on electric-potential from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) );
        electricPotentialMarkers.insert( listMark.begin(), listMark.end() );
    }
    auto meshMarkersElectricPotentialByEntities = detail::distributeMarkerListOnSubEntity( mesh, electricPotentialMarkers );

    // on topological faces
    auto const& listMarkedFacesElectricPotential = std::get<0>( meshMarkersElectricPotentialByEntities );
    if ( !listMarkedFacesElectricPotential.empty() )
        this->updateDofEliminationIds( "potential-electric", XhElectricPotential, markedfaces( mesh,listMarkedFacesElectricPotential ) );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("Electric","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->setPostProcessExportsAllFieldsAvailable( {"electric-potential","electric-field","current-density","joules-losses"} );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->physics() ) );
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

    // point measures
    auto fieldNamesWithSpaceElectricPotential = std::make_pair( std::set<std::string>({"electric-potential"}), this->spaceElectricPotential() );
    //auto fieldNamesWithSpaceElectricField = std::make_pair( std::set<std::string>({"electric-field"}), this->spaceElectricField() );
    auto fieldNamesWithSpaces = hana::make_tuple( fieldNamesWithSpaceElectricPotential/*, fieldNamesWithSpaceElectricField*/ );
    M_measurePointsEvaluation = std::make_shared<measure_points_evaluation_type>( fieldNamesWithSpaces );
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( this->keyword() ) )
    {
       M_measurePointsEvaluation->init( evalPoints );
    }

    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasuresIO().restart( "time", this->timeInitial() );
        else
            this->postProcessMeasuresIO().setMeasure( "time", this->timeInitial() ); //just for have time in first column
    }

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("Electric","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateInformationObject( pt::ptree & p )
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.get_child_optional( "Prefix" ) )
        return;

    p.put( "Prefix", this->prefix() );
    p.put( "Root Repository", this->rootRepository() );

    // Physical Model
    pt::ptree subPt, subPt2;
    subPt.put( "time mode", std::string( (this->isStationary())?"Stationary":"Transient") );
    p.put_child( "Physical Model", subPt );

    // Boundary Conditions
    subPt.clear();
    subPt2.clear();
    M_bcDirichletMarkerManagement.updateInformationObjectDirichletBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    subPt2.clear();
    M_bcNeumannMarkerManagement.updateInformationObjectNeumannBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    subPt2.clear();
    M_bcRobinMarkerManagement.updateInformationObjectRobinBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    p.put_child( "Boundary Conditions",subPt );

#if 0
    // Materials parameters
    subPt.clear();
    this->thermalProperties()->updateInformationObject( subPt );
    p.put_child( "Materials parameters", subPt );
#endif

    // Mesh and FunctionSpace
    subPt.clear();
    subPt.put("filename", this->meshFile());
    M_mesh->putInformationObject( subPt );
    p.put( "Mesh",  M_mesh->journalSectionName() );
    p.put( "FunctionSpace ElectricPotential",  M_XhElectricPotential->journalSectionName() );

    // Algebraic Solver
    if ( M_algebraicFactory )
    {
        subPt.clear();
        M_algebraicFactory->updateInformationObject( subPt );
        p.put_child( "Algebraic Solver", subPt );
    }

}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
ELECTRIC_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
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
    *_ostr << this->materialsProperties()->getInfoMaterialParameters()->str();
    *_ostr << "\n   Mesh Discretization"
           << "\n     -- mesh filename      : " << this->meshFile()
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

    return _ostr;
}



ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->exportResults( time, this->symbolsExpr() );
}
#if 0
ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateElectricField()
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
#endif

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );

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
        this->materialsProperties()->setParameterValues( paramValues );
    }
    M_bcDirichlet.setParameterValues( paramValues );
    M_bcNeumann.setParameterValues( paramValues );
    M_bcRobin.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );

    this->log("Electric","setParameterValues", "finish");
}


ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Electric","solve", "start");
    this->timerTool("Solve").start();

    this->setStartBlockSpaceIndex( 0 );

    M_blockVectorSolution.updateVectorFromSubVectors();
    M_algebraicFactory->solve( "LinearSystem", M_blockVectorSolution.vectorMonolithic() );
    M_blockVectorSolution.localize();

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("Electric","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}



ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    this->updateLinearPDE( data, this->symbolsExpr() );
}


ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    if ( M_bcDirichlet.empty() ) return;

    this->log("Electric","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto mesh = this->mesh();
    size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "potential-electric" );
    auto v = this->spaceElectricPotential()->element( U, this->rowStartInVector()+startBlockIndexElectricPotential );
    for( auto const& d : M_bcDirichlet )
    {
        v.on(_range=markedfaces(mesh, M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) ) ),
             _expr=expression(d) );
    }

    // update info for synchronization
    this->updateDofEliminationIds( "potential-electric", data );

    this->log("Electric","updateNewtonInitialGuess","finish" );
}
ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateJacobian( data, this->symbolsExpr( this->modelFields( XVec, this->rowStartInVector() ) ) );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    if ( this->M_bcDirichlet.empty() ) return;

    this->log("Electric","updateJacobianDofElimination","start" );

    this->updateDofEliminationIds( "potential-electric", data );

    this->log("Electric","updateJacobianDofElimination","finish" );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateResidual( data, this->symbolsExpr( this->modelFields( XVec, this->rowStartInVector() ) ) );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( this->M_bcDirichlet.empty() ) return;

    this->log("Electric","updateResidualDofElimination","start" );

    this->updateDofEliminationIds( "potential-electric", data );

    this->log("Electric","updateResidualDofElimination","finish" );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    if ( this->M_bcDirichlet.empty() ) return;

    this->log("Electric","updateLinearPDEDofElimination","start" );

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
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
            on( _range=markedfaces(mesh, M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) ) ),
                _element=v,_rhs=F,_expr=expression(d) );
    }

    this->log("Electric","updateLinearPDEDofElimination","finish" );
}

} // end namespace FeelModels
} // end namespace Feel
