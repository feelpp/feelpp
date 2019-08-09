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
    M_electricProperties( new electricproperties_type( prefix ) )
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
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().materials().setParameterValues( paramValues );
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
    M_fieldCurrentDensity.reset( new element_electricfield_type(M_XhElectricField,"j"));

    this->initBoundaryConditions();

    // post-process
    this->initPostProcess();

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
ELECTRIC_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerRobinBC();

    this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC("elimination", name(d), markers(d) );
    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "Neumann" );
    for( auto const& d : this->M_bcNeumann )
        this->addMarkerNeumannBC(NeumannBCShape::SCALAR,name(d),markers(d));

    this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "electric-potential", "Robin" );
    for( auto const& d : this->M_bcRobin )
        this->addMarkerRobinBC( name(d),markers(d) );

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().getScalarFields( "electric-potential", "VolumicForces" );

    auto mesh = this->mesh();
    auto XhElectricPotential = this->spaceElectricPotential();

    std::set<std::string> electricPotentialMarkers;

    // strong Dirichlet bc on electric-potential from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = this->markerDirichletBCByNameId( "elimination",name(d) );
        electricPotentialMarkers.insert( listMark.begin(), listMark.end() );
    }
    auto meshMarkersElectricPotentialByEntities = detail::distributeMarkerListOnSubEntity( mesh, electricPotentialMarkers );

    // on topological faces
    auto const& listMarkedFacesElectricPotential = std::get<0>( meshMarkersElectricPotentialByEntities );
    if ( !listMarkedFacesElectricPotential.empty() )
    {
        auto therange = markedfaces( mesh,listMarkedFacesElectricPotential );
        auto dofsToAdd = XhElectricPotential->dofs( therange );
        XhElectricPotential->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
        this->dofEliminationIdsAll("potential-electric",MESH_FACES).insert( dofsToAdd.begin(), dofsToAdd.end() );
        auto dofsMultiProcessToAdd = XhElectricPotential->dofs( therange, ComponentType::NO_COMPONENT, true );
        this->dofEliminationIdsMultiProcess("potential-electric",MESH_FACES).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
    }
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
std::set<std::string>
ELECTRIC_CLASS_TEMPLATE_TYPE::postProcessFieldExported( std::set<std::string> const& ifields, std::string const& prefix ) const
{
    std::set<std::string> res;
    for ( auto const& o : ifields )
    {
        if ( o == prefixvm(prefix,"electric-potential") || o == prefixvm(prefix,"all") )
            res.insert( "electric-potential" );
        if ( o == prefixvm(prefix,"electric-field") || o == prefixvm(prefix,"all") )
            res.insert( "electric-field" );
        if ( o == prefixvm(prefix,"conductivity") || o == prefixvm(prefix,"all") )
            res.insert( "conductivity" );
        if ( o == prefixvm(prefix,"current-density") || o == prefixvm(prefix,"all") )
            res.insert( "current-density" );
        if ( o == prefixvm(prefix,"pid") || o == prefixvm(prefix,"all") )
            res.insert( "pid" );
    }
    return res;
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("Electric","initPostProcess", "start");
    this->timerTool("Constructor").start();

    M_postProcessFieldExported.clear();
    M_postProcessFieldExported = this->postProcessFieldExported( this->modelProperties().postProcess().exports( this->keyword() ).fields() );

    if ( !M_postProcessFieldExported.empty() )
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
    this->updateInformationObjectDirichletBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    subPt2.clear();
    this->updateInformationObjectNeumannBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    subPt2.clear();
    this->updateInformationObjectRobinBC( subPt2 );
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
           << this->getInfoDirichletBC()
           << this->getInfoNeumannBC()
           << this->getInfoRobinBC();
    *_ostr << M_electricProperties->getInfoMaterialParameters()->str();
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

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::exportFields( double time )
{
    bool hasFieldToExport = this->updateExportedFields( M_exporter, M_postProcessFieldExported, time );
    if ( hasFieldToExport )
    {
        M_exporter->save();
        this->upload( M_exporter->path() );
    }
}
ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
bool
ELECTRIC_CLASS_TEMPLATE_TYPE::updateExportedFields( export_ptrtype exporter, std::set<std::string> const& fields, double time )
{
    if ( !exporter ) return false;
    if ( !exporter->doExport() ) return false;

    bool hasFieldToExport = false;
    if ( fields.find( "electric-potential" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"electric-potential"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"electric-potential")),
                                     this->fieldElectricPotential() );
        hasFieldToExport = true;
    }
    if ( fields.find( "electric-field" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"electric-fields"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"electric-fields")),
                                     *M_fieldElectricField );
        hasFieldToExport = true;
    }
    if ( fields.find( "conductivity" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"conductivity"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"conductivity")),
                                     M_electricProperties->fieldElectricConductivity() );
        hasFieldToExport = true;
    }
    if ( fields.find( "current-density" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"current-density"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"current-density")),
                                     *M_fieldCurrentDensity );
        hasFieldToExport = true;
    }
    if ( fields.find( "pid" ) != fields.end() )
    {
        exporter->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
        hasFieldToExport = true;
    }
    return hasFieldToExport;
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::exportMeasures( double time )
{
    this->exportMeasures( time, this->symbolsExpr() );
}

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

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateCurrentDensity()
{
    for ( auto const& rangeData : M_electricProperties->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& electricConductivity = M_electricProperties->electricConductivity( matName );
        if ( electricConductivity.isConstant() )
        {
            double sigma = electricConductivity.value();
            auto cd = -sigma*trans(gradv( this->fieldElectricPotential() ));
            this->updateCurrentDensity( cd, range );
        }
        else
        {
            auto sigma = electricConductivity.expr();
            if ( sigma.expression().hasSymbol( "heat_T" ) )
                continue;
            auto cd = -sigma*trans(gradv( this->fieldElectricPotential() ));
            this->updateCurrentDensity( cd, range );
        }
    }
}



ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();

    this->electricProperties()->setParameterValues( paramValues );
    M_bcDirichlet.setParameterValues( paramValues );
    M_bcNeumann.setParameterValues( paramValues );
    M_bcRobin.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );
}


ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Electric","solve", "start");
    this->timerTool("Solve").start();

    this->updateParameterValues();

    this->setStartBlockSpaceIndex( 0 );

    M_blockVectorSolution.updateVectorFromSubVectors();
    M_algebraicFactory->solve( "LinearSystem", M_blockVectorSolution.vectorMonolithic() );
    M_blockVectorSolution.localize();

    this->updateElectricField();
    this->updateCurrentDensity();

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
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Electric","updateLinearPDE", "start"+sc);
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

    if ( buildCstPart )
    {
        for ( auto const& rangeData : M_electricProperties->rangeMeshElementsByMaterial() )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            auto const& electricConductivity = M_electricProperties->electricConductivity( matName );
            if ( electricConductivity.isConstant() )
            {
                double sigma = electricConductivity.value();
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= sigma*inner(gradt(v),grad(v)),
                               _geomap=this->geomap() );
            }
            else
            {
                auto sigma = electricConductivity.expr();
                if ( sigma.expression().hasSymbol( "heat_T" ) )
                    continue;
                //auto sigma = idv(M_electricProperties->fieldElectricConductivity());
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= sigma*inner(gradt(v),grad(v)),
                               _geomap=this->geomap() );
            }
        }
    }

    // update source term
    if ( !buildCstPart )
    {
        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto rangeEltUsed = (markers(d).empty())? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeEltUsed,
                           _expr= expression(d)*id(v),
                           _geomap=this->geomap() );
        }
    }


    // update bc
    this->updateLinearPDEWeakBC(A,F,buildCstPart);
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
        v.on(_range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
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
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool _BuildCstPart = data.buildCstPart();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Electric","updateJacobian", "start"+sc);
    size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "potential-electric" );

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    // auto const& v = this->fieldElectricPotential();
    auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                                              _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );

    for ( auto const& rangeData : M_electricProperties->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& electricConductivity = M_electricProperties->electricConductivity( matName );
        if ( electricConductivity.isConstant() )
        {
            if ( buildCstPart )
            {
                double sigma = electricConductivity.value();
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= sigma*inner(gradt(v),grad(v)),
                               _geomap=this->geomap() );
            }
        }
        else
        {
            auto sigma = electricConductivity.expr();
            if ( sigma.expression().hasSymbol( "heat_T" ) )
                continue;
            if ( buildCstPart )
            {
                //auto sigma = idv(M_electricProperties->fieldElectricConductivity());
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= sigma*inner(gradt(v),grad(v)),
                               _geomap=this->geomap() );
            }
        }
    }

    this->updateJacobianWeakBC( v,J,buildCstPart );
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
ELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobianWeakBC( element_electricpotential_external_storage_type const& v, sparse_matrix_ptrtype& J, bool buildCstPart ) const
{
    if ( this->M_bcRobin.empty() ) return;

    if ( !buildCstPart )
    {
        auto XhV = this->spaceElectricPotential();
        auto mesh = XhV->mesh();
        size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "potential-electric" );

        auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                                                  _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );
        for( auto const& d : this->M_bcRobin )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= expression1(d)*idt(v)*id(v),
                           _geomap=this->geomap() );
        }
    }
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool _BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Electric","updateResidual", "start"+sc);

    size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "potential-electric" );

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    // auto const& v = this->fieldElectricPotential();
    auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );


    auto myLinearForm = form1( _test=XhV, _vector=R,
                               _rowstart=this->rowStartInVector() + startBlockIndexElectricPotential );


    for ( auto const& rangeData : M_electricProperties->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& electricConductivity = M_electricProperties->electricConductivity( matName );
        if ( electricConductivity.isConstant() )
        {
            if (!buildCstPart && !UseJacobianLinearTerms )
            {
                double sigma = electricConductivity.value();
                myLinearForm +=
                    integrate( _range=range,
                               _expr= sigma*inner(gradv(v),grad(v)),
                               _geomap=this->geomap() );
            }
        }
        else
        {
            auto sigma = electricConductivity.expr();
            if ( sigma.expression().hasSymbol( "heat_T" ) )
                continue;
            if (!buildCstPart && !UseJacobianLinearTerms )
            {
                //auto sigma = idv(M_electricProperties->fieldElectricConductivity());
                myLinearForm +=
                    integrate( _range=range,
                               _expr= sigma*inner(gradv(v),grad(v)),
                               _geomap=this->geomap() );
            }
        }
    }
    // source term
    if ( buildCstPart )
    {
        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto rangeEltUsed = (markers(d).empty())? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeEltUsed,
                           _expr= -expression(d)*id(v),
                           _geomap=this->geomap() );
        }
    }

    // weak bc
    this->updateResidualWeakBC( v,R,buildCstPart );

    this->log("Electric","updateResidual", "finish"+sc);
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
ELECTRIC_CLASS_TEMPLATE_TYPE::updateResidualWeakBC( element_electricpotential_external_storage_type const& v, vector_ptrtype& R, bool buildCstPart ) const
{
    if ( this->M_bcNeumann.empty() && this->M_bcRobin.empty() ) return;

    auto XhV = this->spaceElectricPotential();
    auto mesh = XhV->mesh();
    size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "potential-electric" );

    auto myLinearForm = form1( _test=XhV, _vector=R,
                               _rowstart=this->rowStartInVector()+startBlockIndexElectricPotential );
    if ( buildCstPart )
    {
        for( auto const& d : this->M_bcNeumann )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                           _expr= -expression(d)*id(v),
                           _geomap=this->geomap() );
        }
    }
    for( auto const& d : this->M_bcRobin )
    {
        if ( !buildCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= expression1(d)*idv(v)*id(v),
                           _geomap=this->geomap() );
        }
        if ( buildCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= -expression1(d)*expression2(d)*id(v),
                           _geomap=this->geomap() );
        }
    }
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
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
                _element=v,_rhs=F,_expr=expression(d) );
    }

    this->log("Electric","updateLinearPDEDofElimination","finish" );
}


ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearPDEWeakBC( sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart ) const
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
                integrate( _range=markedfaces(mesh,this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
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
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= expression1(d)*idt(v)*id(v),
                           _geomap=this->geomap() );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= expression1(d)*expression2(d)*id(v),
                           _geomap=this->geomap() );
        }

    }
}


} // end namespace FeelModels
} // end namespace Feel
