/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2014-06-04

  Copyright (C) 2014 Université Joseph Fourier (Grenoble I)

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
   \file solidmechanics.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-06-04
 */

#include <feel/feelmodels/solid/solidmechanics.hpp>

#include <feel/feelmodels/modelvf/solidmecgeomapeulerian.hpp>
#include <feel/feelmodels/modelalg/functionSup.cpp>

namespace Feel
{
namespace FeelModels
{

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::SolidMechanics( std::string __prefix,
                                                    bool __buildMesh,
                                                    WorldComm const& __worldComm,
                                                    std::string __subPrefix,
                                                    std::string __appliShortRepository )
    :
    super_type(__prefix,__buildMesh,__worldComm,__subPrefix, __appliShortRepository)
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());

    //-----------------------------------------------------------------------------//
    // load info from .bc file
    this->loadConfigBCFile();
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // set of worldComm for the function spaces
    this->createWorldsComm();
    //-----------------------------------------------------------------------------//
    // build  mesh, space,exporter,...
    if (__buildMesh) this->build();
    //-----------------------------------------------------------------------------//

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerFluidStructureInterfaceBC();
    this->clearMarkerRobinBC();


    fs::path curPath=fs::current_path();
    bool hasChangedRep=false;
    if ( curPath != fs::path(this->ginacExprCompilationDirectory()) )
    {
        this->log("SolidMechanics","loadConfigBCFile", "change repository (temporary) for build ginac expr with default name : "+ this->appliRepository() );
        bool hasChangedRep=true;
        Environment::changeRepository( _directory=boost::format(this->ginacExprCompilationDirectory()), _subdir=false );
    }

    std::string dirichletbcType = "elimination";//soption(_name="dirichletbc.type",_prefix=this->prefix());
    M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "Dirichlet" );
    for( auto const& d : M_bcDirichlet )
        this->addMarkerDirichletBC( dirichletbcType, marker(d), ComponentType::NO_COMPONENT );
    M_bcDirichletX = this->modelProperties().boundaryConditions().getScalarFields( "displacement_x", "Dirichlet" );
    for( auto const& d : M_bcDirichletX )
        this->addMarkerDirichletBC( dirichletbcType, marker(d), ComponentType::X );
    M_bcDirichletY = this->modelProperties().boundaryConditions().getScalarFields( "displacement_y", "Dirichlet" );
    for( auto const& d : M_bcDirichletY )
        this->addMarkerDirichletBC( dirichletbcType, marker(d), ComponentType::Y );
    M_bcDirichletZ = this->modelProperties().boundaryConditions().getScalarFields( "displacement_z", "Dirichlet" );
    for( auto const& d : M_bcDirichletZ )
        this->addMarkerDirichletBC( dirichletbcType, marker(d), ComponentType::Z );
    M_bcNeumannScalar = this->modelProperties().boundaryConditions().getScalarFields( "displacement", "Neumann_scalar" );
    for( auto const& d : M_bcNeumannScalar )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d));
    M_bcNeumannVectorial = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "Neumann_vectorial" );
    for( auto const& d : M_bcNeumannVectorial )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::VECTORIAL,marker(d));
    M_bcInterfaceFSI = this->modelProperties().boundaryConditions().getScalarFields( "displacement", "interface_fsi" );
    for( auto const& d : M_bcInterfaceFSI )
        this->addMarkerFluidStructureInterfaceBC( marker(d) );
    M_bcRobin = this->modelProperties().boundaryConditions().template getVectorFieldsList<super_type::nDim>( "displacement", "robin" );
    for( auto const& d : M_bcRobin )
        this->addMarkerRobinBC( marker(d) );

    M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "VolumicForces" );

    // go back to previous repository
    if ( hasChangedRep )
        Environment::changeRepository( _directory=boost::format(curPath.string()), _subdir=false );

    if ( this->modelProperties().postProcess().find("Fields") != this->modelProperties().postProcess().end() )
        for ( auto const& o : this->modelProperties().postProcess().find("Fields")->second )
        {
            if ( o == "displacement" || o == "all" ) this->M_doExportDisplacement = true;
            if ( o == "velocity" || o == "all" ) this->M_doExportVelocity = true;
            if ( o == "acceleration" || o == "all" ) this->M_doExportAcceleration = true;
            if ( o == "stress" || o == "normal-stress" || o == "all" ) this->M_doExportNormalStress = true;
            if ( o == "pressure" || o == "all" ) this->M_doExportPressure = true;
            if ( o == "all" ) this->M_doExportVelocityInterfaceFromFluid = true;
        }


}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::loadConfigMeshFile( std::string const& geofilename )
{
    CHECK( false ) << "not allow";
}

//---------------------------------------------------------------------------------------------------//
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::loadConfigMeshFile1dReduced( std::string const& geofilename )
{
    CHECK( false ) << "not allow";
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::solve( bool upVelAcc )
{
    M_bcDirichlet.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcDirichletX.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcDirichletY.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcDirichletZ.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcNeumannScalar.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcNeumannVectorial.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcRobin.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_volumicForcesProperties.setParameterValues( this->modelProperties().parameters().toParameterValues() );

    super_type::solve( upVelAcc );
}

//---------------------------------------------------------------------------------------------------//
namespace detail
{
template <typename MeshType>
boost::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string> >
distributeMarkerListOnSubEntity( boost::shared_ptr<MeshType> const& mesh, std::list<std::string> const& listMarker )
{
    boost::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string> > res;
    for ( std::string const& marker : listMarker )
    {
        if ( !mesh->hasMarker( marker ) ) continue;

        if ( /*mesh->markerDim( marker ) == (MeshType::nDim-1)*/mesh->hasFaceMarker( marker ) )
        {
            res.template get<0>().push_back( marker );
            //std::cout << "has face marker " << marker << "\n";
        }
        else if ( /*MeshType::nDim == 3 && mesh->markerDim( marker ) == (MeshType::nDim-2)*/ mesh->hasEdgeMarker( marker ) )
        {
            res.template get<1>().push_back( marker );
            //std::cout << "has edge marker " << marker << "\n";
        }
        else if ( /*mesh->markerDim( marker ) == 0*/ mesh->hasPointMarker( marker ) )
        {
            res.template get<2>().push_back( marker );
            //std::cout << "has point marker " << marker << "\n";
        }
        else
        {
            //std::cout << "unknow marker " << marker << " with dim " << mesh->markerDim( marker ) << "\n";
        }
    }
    return res;
}
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
    if ( M_bcDirichlet.empty() && M_bcDirichletX.empty() && M_bcDirichletY.empty() && M_bcDirichletZ.empty() ) return;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateNewtonInitialGuess", "start",
                                               this->worldComm(),this->verboseAllProc());

    auto const& u = this->fieldDisplacement();
    auto Xh = this->functionSpace();
    auto mesh = this->mesh();
    size_type rowStartInVector = this->rowStartInVector();

    if (Xh->worldsComm()[0].isActive()) // only on Displacement Proc
    {
        for( auto const& d : M_bcDirichlet )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d) ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                modifVec(markedfaces(mesh,listMarkerFaces ),
                         u, U, expression(d), rowStartInVector );
            if ( !listMarkerEdges.empty() )
                modifVec(markededges(mesh,listMarkerEdges),
                         u, U, expression(d), rowStartInVector );
            if ( !listMarkerPoints.empty() )
                modifVec(markedpoints(mesh,listMarkerPoints),
                         u, U, expression(d), rowStartInVector );
        }
        for( auto const& d : M_bcDirichletX )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d), Component::X ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                modifVec(markedfaces(mesh,listMarkerFaces),
                         u[Component::X], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerEdges.empty() )
                modifVec(markededges(mesh,listMarkerEdges),
                         u[Component::X], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerPoints.empty() )
                modifVec(markedpoints(mesh,listMarkerPoints),
                         u[Component::X], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
        }
        for( auto const& d : M_bcDirichletY )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d), Component::Y ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                modifVec(markedfaces(mesh,listMarkerFaces),
                         u[Component::Y], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerEdges.empty() )
                modifVec(markededges(mesh,listMarkerEdges),
                         u[Component::Y], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerPoints.empty() )
                modifVec(markedpoints(mesh,listMarkerPoints),
                         u[Component::Y], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
        }
        for( auto const& d : M_bcDirichletZ )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d), Component::Z ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                modifVec(markedfaces(mesh,listMarkerFaces),
                         u[Component::Z], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerEdges.empty() )
                modifVec(markededges(mesh,listMarkerEdges),
                         u[Component::Z], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerPoints.empty() )
                modifVec(markedpoints(mesh,listMarkerPoints),
                         u[Component::Z], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
        }
    }
    U->close();

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateNewtonInitialGuess", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletStrongResidual(vector_ptrtype& R) const
{
    if ( M_bcDirichlet.empty() && M_bcDirichletX.empty() && M_bcDirichletY.empty() && M_bcDirichletZ.empty() ) return;

    auto const& u = this->fieldDisplacement();
    size_type rowStartInVector = this->rowStartInVector();

    R->close();
    if ( this->functionSpaceDisplacement()->worldsComm()[0].isActive() ) // only on Displacement Proc
    {
        for( auto const& d : M_bcDirichlet )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                modifVec( markedfaces(this->mesh(),listMarkerFaces),
                          u, R, 0*vf::one(), rowStartInVector );
            if ( !listMarkerEdges.empty() )
                modifVec( markededges(this->mesh(),listMarkerEdges),
                          u, R, 0*vf::one(), rowStartInVector );
            if ( !listMarkerPoints.empty() )
                modifVec( markedpoints(this->mesh(),listMarkerPoints),
                          u, R, 0*vf::one(), rowStartInVector );
        }
        for( auto const& d : M_bcDirichletX )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d), Component::X ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                modifVec( markedfaces(this->mesh(),listMarkerFaces),
                          u[Component::X], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerEdges.empty() )
                modifVec( markededges(this->mesh(),listMarkerEdges),
                          u[Component::X], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerPoints.empty() )
                modifVec( markedpoints(this->mesh(),listMarkerPoints),
                          u[Component::X], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
        }
        for( auto const& d : M_bcDirichletY )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d), Component::Y ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                modifVec( markedfaces(this->mesh(),listMarkerFaces),
                          u[Component::Y], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerEdges.empty() )
                modifVec( markededges(this->mesh(),listMarkerEdges),
                          u[Component::Y], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerPoints.empty() )
                modifVec( markedpoints(this->mesh(),listMarkerPoints),
                          u[Component::Y], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
        }
        for( auto const& d : M_bcDirichletZ )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d), Component::Z ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                modifVec( markedfaces(this->mesh(),listMarkerFaces),
                          u[Component::Z], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerEdges.empty() )
                modifVec( markededges(this->mesh(),listMarkerEdges),
                          u[Component::Z], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
            if ( !listMarkerPoints.empty() )
                modifVec( markedpoints(this->mesh(),listMarkerPoints),
                          u[Component::Z], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
        }
    }
}
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletStrongJacobian( sparse_matrix_ptrtype& J ) const
{
    if ( M_bcDirichlet.empty() && M_bcDirichletX.empty() && M_bcDirichletY.empty() && M_bcDirichletZ.empty() ) return;

    auto RBis = this->backend()->newVector( J->mapRowPtr() );
    auto bilinearForm_PatternCoupled = form2( _test=this->functionSpace(),_trial=this->functionSpace(),_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldDisplacement();
    for( auto const& d : M_bcDirichlet )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) );
        auto const& listMarkerFaces = ret.template get<0>();
        auto const& listMarkerEdges = ret.template get<1>();
        auto const& listMarkerPoints = ret.template get<2>();
        if ( !listMarkerFaces.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedfaces(this->mesh(), listMarkerFaces),
                    _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/ );
        if ( !listMarkerEdges.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markededges(this->mesh(), listMarkerEdges),
                    _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/ );
        if ( !listMarkerPoints.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedpoints(this->mesh(), listMarkerPoints),
                    _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/ );
    }
    for( auto const& d : M_bcDirichletX )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),Component::X ) );
        auto const& listMarkerFaces = ret.template get<0>();
        auto const& listMarkerEdges = ret.template get<1>();
        auto const& listMarkerPoints = ret.template get<2>();
        if ( !listMarkerFaces.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedfaces(this->mesh(), listMarkerFaces),
                    _element=u[Component::X],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
        if ( !listMarkerEdges.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markededges(this->mesh(), listMarkerEdges),
                    _element=u[Component::X],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
        if ( !listMarkerPoints.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedpoints(this->mesh(), listMarkerPoints),
                    _element=u[Component::X],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
    }
    for( auto const& d : M_bcDirichletY )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),Component::Y ) );
        auto const& listMarkerFaces = ret.template get<0>();
        auto const& listMarkerEdges = ret.template get<1>();
        auto const& listMarkerPoints = ret.template get<2>();
        if ( !listMarkerFaces.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedfaces(this->mesh(), listMarkerFaces),
                    _element=u[Component::Y],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
        if ( !listMarkerEdges.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markededges(this->mesh(), listMarkerEdges),
                    _element=u[Component::Y],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
        if ( !listMarkerPoints.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedpoints(this->mesh(), listMarkerPoints),
                    _element=u[Component::Y],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
    }
    for( auto const& d : M_bcDirichletZ )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),Component::Z ) );
        auto const& listMarkerFaces = ret.template get<0>();
        auto const& listMarkerEdges = ret.template get<1>();
        auto const& listMarkerPoints = ret.template get<2>();
        if ( !listMarkerFaces.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedfaces(this->mesh(), listMarkerFaces),
                    _element=u[Component::Z],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
        if ( !listMarkerEdges.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markededges(this->mesh(), listMarkerEdges),
                    _element=u[Component::Z],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
        if ( !listMarkerPoints.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedpoints(this->mesh(), listMarkerPoints),
                    _element=u[Component::Z],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
    }
}
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletStrongLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    if ( this->isStandardModel() )
    {
        if ( M_bcDirichlet.empty() && M_bcDirichletX.empty() && M_bcDirichletY.empty() && M_bcDirichletZ.empty() ) return;

        auto Xh = this->functionSpace();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->fieldDisplacement();
        for( auto const& d : M_bcDirichlet )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                bilinearForm +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=u,_rhs=F,_expr=expression(d) );
            if ( !listMarkerEdges.empty() )
                bilinearForm +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=u,_rhs=F,_expr=expression(d) );
            if ( !listMarkerPoints.empty() )
                bilinearForm +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=u,_rhs=F,_expr=expression(d) );
        }
        for( auto const& d : M_bcDirichletX )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),Component::X ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                bilinearForm +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=u[Component::X],_rhs=F,_expr=expression(d) );
            if ( !listMarkerEdges.empty() )
                bilinearForm +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=u[Component::X],_rhs=F,_expr=expression(d) );
            if ( !listMarkerPoints.empty() )
                bilinearForm +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=u[Component::X],_rhs=F,_expr=expression(d) );
        }
        for( auto const& d : M_bcDirichletY )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),Component::Y ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                bilinearForm +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=u[Component::Y],_rhs=F,_expr=expression(d) );
            if ( !listMarkerEdges.empty() )
                bilinearForm +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=u[Component::Y],_rhs=F,_expr=expression(d) );
            if ( !listMarkerPoints.empty() )
                bilinearForm +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=u[Component::Y],_rhs=F,_expr=expression(d) );
        }
        for( auto const& d : M_bcDirichletZ )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),Component::Z ) );
            auto const& listMarkerFaces = ret.template get<0>();
            auto const& listMarkerEdges = ret.template get<1>();
            auto const& listMarkerPoints = ret.template get<2>();
            if ( !listMarkerFaces.empty() )
                bilinearForm +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=u[Component::Z],_rhs=F,_expr=expression(d) );
            if ( !listMarkerEdges.empty() )
                bilinearForm +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=u[Component::Z],_rhs=F,_expr=expression(d) );
            if ( !listMarkerPoints.empty() )
                bilinearForm +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=u[Component::Z],_rhs=F,_expr=expression(d) );
        }
    }
    else if ( this->is1dReducedModel() )
    {
        if ( M_bcDirichlet.empty() ) return;

        auto Xh = this->functionSpace1dReduced();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->fieldDisplacementScal1dReduced();
        //WARNING : fixed at zero
        for( auto const& d : M_bcDirichlet )
            bilinearForm +=
                on( _range=markedfaces(Xh->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                    _element=u, _rhs=F,
                    _expr=cst(0.) );
    }
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannResidual(vector_ptrtype& R) const
{
    if ( M_bcNeumannScalar.empty() && M_bcNeumannVectorial.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : M_bcNeumannScalar )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d)) ),
                       _expr= -expression(d)*inner( N(),id(v) ),
                       _geomap=this->geomap() );

    for( auto const& d : M_bcNeumannVectorial )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::VECTORIAL,marker(d)) ),
                       _expr= -inner( expression(d),id(v) ),
                       _geomap=this->geomap() );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
    if ( M_bcNeumannScalar.empty() && M_bcNeumannVectorial.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : M_bcNeumannScalar )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d)) ),
                       _expr= expression(d)*inner( N(),id(v) ),
                        _geomap=this->geomap() );

    for( auto const& d : M_bcNeumannVectorial )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::VECTORIAL,marker(d)) ),
                       _expr= inner( expression(d),id(v) ),
                       _geomap=this->geomap() );
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinResidual(element_displacement_type const& u, vector_ptrtype& R) const
{

    if ( M_bcRobin.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );

    // Warning : take only first component of expression1
    for( auto const& d : M_bcRobin )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)/*this->markerRobinBC()*/),
                       _expr= inner( expression1(d)(0,0)*idv(u) - expression2(d) ,id(u) ),
                       _geomap=this->geomap() );

}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinJacobian( sparse_matrix_ptrtype& J) const
{
    if ( M_bcRobin.empty() ) return;

    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldDisplacement();

    // Warning : take only first component of expression1
    for( auto const& d : M_bcRobin )
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)/*this->markerRobinBC()*/),
                       _expr= expression1(d)(0,0)*inner( idt(u) ,id(u) ),
                       _geomap=this->geomap() );

}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinLinearPDE( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const
{
    if ( M_bcRobin.empty() ) return;

    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& u = this->fieldDisplacement();

    // Warning : take only first component of expression1
    for( auto const& d : M_bcRobin )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)/*this->markerRobinBC()*/),
                       _expr= expression1(d)(0,0)*inner( idt(u) ,id(u) ),
                       _geomap=this->geomap() );
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)/*this->markerRobinBC()*/),
                       _expr= inner( expression2(d) , id(u) ),
                       _geomap=this->geomap() );

    }

}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateSourceTermResidual( vector_ptrtype& R ) const
{

    if ( M_volumicForcesProperties.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : M_volumicForcesProperties )
    {
        if ( marker(d).empty() )
            myLinearForm +=
                integrate( _range=elements(this->mesh()),
                           _expr= -inner( expression(d),id(v) ),
                           _geomap=this->geomap() );
        else
            myLinearForm +=
                integrate( _range=markedelements(this->mesh(),marker(d)),
                           _expr= -inner( expression(d),id(v) ),
                           _geomap=this->geomap() );
    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateSourceTermLinearPDE( vector_ptrtype& F ) const
{
    if ( M_volumicForcesProperties.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : M_volumicForcesProperties )
    {
        if ( marker(d).empty() )
            myLinearForm +=
                integrate( _range=elements(this->mesh()),
                           _expr= inner( expression(d),id(v) ),
                           _geomap=this->geomap() );
        else
            myLinearForm +=
                integrate( _range=markedelements(this->mesh(),marker(d)),
                           _expr= inner( expression(d),id(v) ),
                           _geomap=this->geomap() );
    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCFollowerPressureResidual( typename super_type::element_displacement_type const& u, vector_ptrtype& R ) const
{
#if 0 // TODO
    auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
    if ( !bcDef.hasFollowerPressure() ) return;
    //CHECK( false ) << "TODO";
    auto Xh = this->functionSpaceDisplacement();
    auto linearForm = form1( _test=Xh, _vector=R,_rowstart=this->rowStartInVector() );
    ForEachBC( bcDef,cl::follower_pressure,
               linearForm +=
               /**/ integrate( _range=markedfaces(this->mesh(),PhysicalName),
                               _expr= -Expression*inner(Feel::vf::FeelModels::solidMecGeomapEulerian(u)*N(),id(u) ),
                               _geomap=this->geomap() ); );
#endif
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCFollowerPressureJacobian( typename super_type::element_displacement_type const& u, sparse_matrix_ptrtype& J) const
{
#if 0 // TODO
    auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
    if ( !bcDef.hasFollowerPressure() ) return;
    //CHECK( false ) << "TODO";
    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
   ForEachBC( bcDef,cl::follower_pressure,
               bilinearForm +=
              /**/ integrate( _range=markedfaces(this->mesh(),PhysicalName) ,
                              _expr= -Expression*inner(Feel::vf::FeelModels::solidMecGeomapEulerianJacobian(u)*N(),id(u) ),
                              _geomap=this->geomap() ) );
#endif
}


} // namespace FeelModels
} // namespace Feel
