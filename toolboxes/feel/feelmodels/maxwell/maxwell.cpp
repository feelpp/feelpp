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
   \file maxwell.cpp
   \author Romain Hild <romain.hild@unistra.fr>
   \date 2018-05-03
 */

#include <feel/feelmodels/maxwell/maxwell.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelmodels/modelmesh/createmesh.hpp>

namespace Feel
{
namespace FeelModels
{

MAXWELL_CLASS_TEMPLATE_DECLARATIONS
MAXWELL_CLASS_TEMPLATE_TYPE::Maxwell( std::string const& prefix,
                                      std::string const& keyword,
                                      worldcomm_ptr_t const& worldComm,
                                      std::string const& subPrefix,
                                      ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    M_maxwellProperties( std::make_shared<maxwellproperties_type>( prefix ) ),
    M_epsilon( doption(_name="regularization-epsilon", _prefix=prefix) )
{
    this->log("Maxwell","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MaxwellConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MaxwellSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MaxwellPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MaxwellTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    // this->setFilenameSaveInfo( prefixvm(this->prefix(),"Maxwell.info") );
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // build mesh
    // if ( buildMesh )
    //     this->createMesh();
    //-----------------------------------------------------------------------------//
    this->log("Maxwell","constructor", "finish");

}

MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{

}

MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("Maxwell","initMesh", "start");
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("initMesh");
    this->log("Maxwell","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // initMesh()


MAXWELL_CLASS_TEMPLATE_DECLARATIONS
int
MAXWELL_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = 1;
    return nBlock;
}


MAXWELL_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
MAXWELL_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    myblockGraph(0,0) = stencil(_test=this->spaceMagneticPotential(),
                                _trial=this->spaceMagneticPotential() )->graph();
    return myblockGraph;
}

MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("Maxwell","init", "start" );
    this->timerTool("Constructor").start();

    if ( !M_mesh )
        this->initMesh();

    // physical properties
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().materials().setParameterValues( paramValues );
    M_maxwellProperties->updateForUse( M_mesh, this->modelProperties().materials(),  this->localNonCompositeWorldsComm());

    // functionspace
    if ( M_maxwellProperties->isDefinedOnWholeMesh() )
    {
        M_rangeMeshElements = elements(M_mesh);
        M_XhMagneticPotential = space_magneticpotential_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
        M_XhMagneticField = space_magneticfield_type::New(_mesh=M_mesh, _worldscomm=this->worldsComm() );
    }
    else
    {
        M_rangeMeshElements = markedelements(M_mesh, M_maxwellProperties->markers());
        M_XhMagneticPotential = space_magneticpotential_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
        M_XhMagneticField = space_magneticfield_type::New(_mesh=M_mesh, _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
    }
    M_fieldMagneticPotential.reset( new element_magneticpotential_type(M_XhMagneticPotential,"A"));
    M_fieldMagneticField.reset( new element_magneticfield_type(M_XhMagneticField,"B"));

    this->initBoundaryConditions();

    // post-process
    this->initPostProcess();

    // update fields
    this->updateFields( this->symbolsExpr() );

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() );

    size_type currentStartIndex = 0;// velocity and pressure before
    this->setStartSubBlockSpaceIndex( "potential-maxwell", currentStartIndex );

    // vector solution
    int nBlock = this->nBlockMatrixGraph();
    M_blockVectorSolution.resize( nBlock );
    int indexBlock=0;
    M_blockVectorSolution(indexBlock) = this->fieldMagneticPotentialPtr();

    // init petsc vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );


    // algebraic solver
    if ( buildModelAlgebraicFactory )
        this->initAlgebraicFactory();

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("Maxwell","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerRobinBC();

    this->initDirichlet();
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC("nitsche", name(d), markers(d) );
        // this->addMarkerDirichletBC("elimination", marker(d) );
    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "magnetic-potential", "Neumann" );
    for( auto const& d : this->M_bcNeumann )
        this->addMarkerNeumannBC(NeumannBCShape::SCALAR,name(d),markers(d));

    this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "magnetic-potential", "Robin" );
    for( auto const& d : this->M_bcRobin )
        this->addMarkerRobinBC( name(d),markers(d) );

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "magnetic-potential", "VolumicForces" );

    auto mesh = this->mesh();
    auto XhMagneticPotential = this->spaceMagneticPotential();

    // auto & dofsWithValueImposedMagneticPotential = M_dofsWithValueImposed["magnetic-potential"];
    // dofsWithValueImposedMagneticPotential.clear();
    std::set<std::string> magneticPotentialMarkers;

    // strong Dirichlet bc on magnetic-potential from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = this->markerDirichletBCByNameId( "elimination",name(d) );
        magneticPotentialMarkers.insert( listMark.begin(), listMark.end() );
    }
    auto meshMarkersMagneticPotentialByEntities = detail::distributeMarkerListOnSubEntity( mesh, magneticPotentialMarkers );

    // // on topological faces
    // auto const& listMarkedFacesMagneticPotential = std::get<0>( meshMarkersMagneticPotentialByEntities );
    // for ( auto const& faceWrap : markedfaces(mesh,listMarkedFacesMagneticPotential ) )
    // {
    //     auto const& face = unwrap_ref( faceWrap );
    //     auto facedof = XhMagneticPotential->dof()->faceLocalDof( face.id() );
    //     for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
    //         dofsWithValueImposedMagneticPotential.insert( it->index() );
    // }
    // on topological faces
    auto const& listMarkedFacesMagneticPotential = std::get<0>( meshMarkersMagneticPotentialByEntities );
    if ( !listMarkedFacesMagneticPotential.empty() )
    {
        auto therange = markedfaces( mesh,listMarkedFacesMagneticPotential );
        auto dofsToAdd = XhMagneticPotential->dofs( therange );
        XhMagneticPotential->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
        this->dofEliminationIdsAll("potential-magnetic",MESH_FACES).insert( dofsToAdd.begin(), dofsToAdd.end() );
        auto dofsMultiProcessToAdd = XhMagneticPotential->dofs( therange, ComponentType::NO_COMPONENT, true );
        this->dofEliminationIdsMultiProcess("potential-magnetic",MESH_FACES).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
    }
}



// MAXWELL_CLASS_TEMPLATE_DECLARATIONS
// std::set<std::string>
// MAXWELL_CLASS_TEMPLATE_TYPE::postProcessFieldExported( std::set<std::string> const& ifields, std::string const& prefix ) const
// {
//     std::set<std::string> res;
//     for ( auto const& o : ifields )
//     {
//         if ( o == prefixvm(prefix,"magnetic-potential") || o == prefixvm(prefix,"all") )
//             res.insert( "magnetic-potential" );
//         if ( o == prefixvm(prefix,"magnetic-field") || o == prefixvm(prefix,"all") )
//             res.insert( "magnetic-field" );
//         // if ( o == prefixvm(prefix,"conductivity") || o == prefixvm(prefix,"all") )
//         //     res.insert( "conductivity" );
//         if ( o == prefixvm(prefix,"pid") || o == prefixvm(prefix,"all") )
//             res.insert( "pid" );
//     }
//     return res;
// }

MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("Maxwell","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->setPostProcessExportsAllFieldsAvailable( {"magnetic-potential","magnetic-field","magnetic-permeability","pid"} );
    this->setPostProcessSaveAllFieldsAvailable( {"magnetic-potential","magnetic-field","magnetic-permeability"} );
    super_type::initPostProcess();

    // M_postProcessFieldExported.clear();
    // std::string modelName = "maxwell";
    // M_postProcessFieldExported = this->postProcessFieldExported( this->modelProperties().postProcess().exports( modelName ).fields() );

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
    auto fieldNamesWithSpaceMagneticPotential = std::make_pair( std::set<std::string>({"magnetic-potential"}), this->spaceMagneticPotential() );
    auto fieldNamesWithSpaceMagneticField = std::make_pair( std::set<std::string>({"magnetic-field"}), this->spaceMagneticField() );
    auto fieldNamesWithSpaces = hana::make_tuple( fieldNamesWithSpaceMagneticPotential, fieldNamesWithSpaceMagneticField );
    M_measurePointsEvaluation = std::make_shared<measure_points_evaluation_type>( fieldNamesWithSpaces );
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( this->keyword() ) )
    {
       M_measurePointsEvaluation->init( evalPoints );
    }

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("Maxwell","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );
}

MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::updateInformationObject( pt::ptree & p )
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
    p.put( "FunctionSpace MagneticPotential",  M_XhMagneticPotential->journalSectionName() );

    // Algebraic Solver
    if ( M_algebraicFactory )
    {
        subPt.clear();
        M_algebraicFactory->updateInformationObject( subPt );
        p.put_child( "Algebraic Solver", subPt );
    }

}

MAXWELL_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
MAXWELL_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||-----------Info : Maxwell--------------------||"
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
    *_ostr << M_maxwellProperties->getInfoMaterialParameters()->str();
    *_ostr << "\n   Mesh Discretization"
           << "\n     -- mesh filename      : " << this->meshFile()
           << "\n     -- number of element : " << M_mesh->numGlobalElements()
           << "\n     -- order             : " << nOrderGeo;
    *_ostr << "\n   Space MagneticPotential Discretization"
#if MAXWELL_H1
           << "\n     -- family        : " << "Lagrange"
#else
           << "\n     -- family        : " << "Nedelec"
#endif
           << "\n     -- order         : " << nOrderPolyMagneticPotential
           << "\n     -- number of dof : " << M_XhMagneticPotential->nDof() << " (" << M_XhMagneticPotential->nLocalDof() << ")";
    if ( M_algebraicFactory )
        *_ostr << M_algebraicFactory->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}


MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    // this->log("Maxwell","exportResults", "start");
    // this->timerTool("PostProcessing").start();

    // this->exportFields( time );

    // this->exportMeasures( time );

    // this->timerTool("PostProcessing").stop("exportResults");
    // if ( this->scalabilitySave() )
    // {
    //     if ( !this->isStationary() )
    //         this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
    //     this->timerTool("PostProcessing").save();
    // }
    // this->log("Maxwell","exportResults", "finish");
    this->exportResults( time, this->symbolsExpr() );
}

// MAXWELL_CLASS_TEMPLATE_DECLARATIONS
// void
// MAXWELL_CLASS_TEMPLATE_TYPE::exportFields( double time )
// {
//     bool hasFieldToExport = this->updateExportedFields( M_exporter, M_postProcessFieldExported, time );
//     if ( hasFieldToExport )
//         M_exporter->save();
// }

// MAXWELL_CLASS_TEMPLATE_DECLARATIONS
// bool
// MAXWELL_CLASS_TEMPLATE_TYPE::updateExportedFields( export_ptrtype exporter, std::set<std::string> const& fields, double time )
// {
//     if ( !exporter ) return false;
//     if ( !exporter->doExport() ) return false;

//     bool hasFieldToExport = false;
//     if ( fields.find( "magnetic-potential" ) != fields.end() )
//     {
//         exporter->step( time )->add( prefixvm(this->prefix(),"magnetic-potential"),
//                                      prefixvm(this->prefix(),prefixvm(this->subPrefix(),"magnetic-potential")),
//                                      this->fieldMagneticPotential() );
//         hasFieldToExport = true;
//     }
//     if ( fields.find( "magnetic-field" ) != fields.end() )
//     {
//         exporter->step( time )->add( prefixvm(this->prefix(),"magnetic-fields"),
//                                      prefixvm(this->prefix(),prefixvm(this->subPrefix(),"magnetic-fields")),
//                                      *M_fieldMagneticField );
//         hasFieldToExport = true;
//     }
//     // if ( fields.find( "conductivity" ) != fields.end() )
//     // {
//     //     exporter->step( time )->add( prefixvm(this->prefix(),"conductivity"),
//     //                                  prefixvm(this->prefix(),prefixvm(this->subPrefix(),"conductivity")),
//     //                                  M_maxwellProperties->fieldElectricConductivity() );
//     //     hasFieldToExport = true;
//     // }
//     if ( fields.find( "pid" ) != fields.end() )
//     {
//         exporter->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
//         hasFieldToExport = true;
//     }
//     return hasFieldToExport;
// }

// MAXWELL_CLASS_TEMPLATE_DECLARATIONS
// void
// MAXWELL_CLASS_TEMPLATE_TYPE::exportMeasures( double time )
// {

// }

MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::updateMagneticField()
{
    std::string M_computeMagneticFieldProjType = "nodal";
    /*if ( M_computeMagneticFieldProjType == "L2" )
     *M_fieldMagneticFieldContinuous = M_l2proj->operator()( -trans(gradv(this->fieldMagneticPotential() ) ) );
     else*/ if ( M_computeMagneticFieldProjType == "nodal" )
        M_fieldMagneticField->on(_range=M_rangeMeshElements,
// #if FEELPP_DIM==3
                                 _expr=curlv( this->fieldMagneticPotential() )
// #else
//                                  _expr=curlxv( this->fieldMagneticPotential() )
// #endif
                                 );
    else
        CHECK( false ) << "invalid M_computeMagneticFieldProjType " << M_computeMagneticFieldProjType << "\n";
}


MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();

    this->maxwellProperties()->setParameterValues( paramValues );
    M_bcDirichlet.setParameterValues( paramValues );
    M_bcNeumann.setParameterValues( paramValues );
    M_bcRobin.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );
}


MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Maxwell","solve", "start");
    this->timerTool("Solve").start();

    this->updateParameterValues();

    this->setStartBlockSpaceIndex( 0 );

    M_blockVectorSolution.updateVectorFromSubVectors();
    M_algebraicFactory->solve( "LinearSystem", M_blockVectorSolution.vectorMonolithic() );
    M_blockVectorSolution.localize();

    this->updateMagneticField();

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("Maxwell","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}


MAXWELL_CLASS_TEMPLATE_DECLARATIONS
void
MAXWELL_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    this->updateLinearPDE( data, this->symbolsExpr() );
}


// MAXWELL_CLASS_TEMPLATE_DECLARATIONS
// void
// MAXWELL_CLASS_TEMPLATE_TYPE::updateLinearPDEStrongDirichletBC( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const
// {
//     if ( this->M_bcDirichlet.empty() ) return;

//     this->log("Maxwell","updateBCStrongDirichletLinearPDE","start" );

//     auto XhV = this->spaceMagneticPotential();
//     auto const& v = this->fieldMagneticPotential();
//     auto mesh = XhV->mesh();

//     auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
//                                               _pattern=size_type(Pattern::COUPLED),
//                                               _rowstart=this->rowStartInMatrix(),
//                                               _colstart=this->colStartInMatrix() );
//     for( auto const& d : this->M_bcDirichlet )
//     {
//         bilinearForm_PatternCoupled +=
//             on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
//                 _element=v,_rhs=F,_expr=expression(d) );
//     }

//     this->log("Maxwell","updateBCStrongDirichletLinearPDE","finish" );
// }


// MAXWELL_CLASS_TEMPLATE_DECLARATIONS
// void
// MAXWELL_CLASS_TEMPLATE_TYPE::updateLinearPDEWeakBC( sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart ) const
// {
//     if ( this->M_bcNeumann.empty() && this->M_bcRobin.empty() && this->M_bcDirichlet.empty() ) return;

//     if ( !buildCstPart )
//     {
//         auto XhA = this->spaceMagneticPotential();
//         auto const& v = this->fieldMagneticPotential();
//         auto mesh = XhA->mesh();

//         auto myLinearForm = form1( _test=XhA, _vector=F,
//                                    _rowstart=this->rowStartInVector() );
//         auto bilinearForm_PatternCoupled = form2( _test=XhA,_trial=XhA,_matrix=A,
//                                                   _pattern=size_type(Pattern::COUPLED),
//                                                   _rowstart=this->rowStartInMatrix(),
//                                                   _colstart=this->colStartInMatrix() );

//         for( auto const& d : this->M_bcDirichlet )
//         {
//             myLinearForm +=
//                 integrate( _range=markedfaces(mesh,this->markerDirichletBCByNameId("nitsche",name(d)) ),
// // #if FEELPP_DIM==3
//                            _expr= trans(expression(d))*curl(v) + 1e5*trans(expression(d))*cross(id(v),N())/hFace(),
// // #else
// //                            _expr=expression(d)*curlx(v) + 1e5*expression(d)*cross(id(v),N())/hFace(),
// // #endif
//                            _geomap=this->geomap() );

//             bilinearForm_PatternCoupled +=
//                 integrate( _range=markedfaces(mesh,this->markerDirichletBCByNameId("nitsche",name(d)) ),
// // // #if FEELPP_DIM==3
//                            _expr=trans(curlt(v))*cross(id(v),N()) + trans(curl(v))*cross(idt(v),N())
//                            + 1e5*trans(cross(idt(v),N()))*cross(id(v),N())/hFace(),
// // // #else
// // //                            _expr=curlxt(v)*cross(id(v),N()) + curlx(v)*cross(idt(v),N())
// // //                            + 1e5*trans(cross(idt(v),N()))*cross(id(v),N())/hFace(),
// // // #endif
//                            _geomap=this->geomap() );
//         }
//     }
// }

}// namespace FeelModels
} // namespace Feel
