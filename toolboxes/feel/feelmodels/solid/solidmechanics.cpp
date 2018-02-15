/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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

namespace Feel
{
namespace FeelModels
{

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::SolidMechanics( std::string const& prefix,
                                                    bool buildMesh,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    ModelBaseRepository const& modelRep )
    :
    super_type( prefix, buildMesh, worldComm, subPrefix, modelRep )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());

    //-----------------------------------------------------------------------------//
    // load info from .json file
    this->loadConfigBCFile();
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // set of worldComm for the function spaces
    this->createWorldsComm();
    //-----------------------------------------------------------------------------//
    // build  mesh, space,exporter,...
    if (buildMesh) this->build();
    //-----------------------------------------------------------------------------//

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::self_ptrtype
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::New( std::string const& prefix,
                                         bool buildMesh,
                                         WorldComm const& worldComm,
                                         std::string const& subPrefix,
                                         ModelBaseRepository const& modelRep )
{
    return boost::make_shared<self_type>( prefix, buildMesh, worldComm, subPrefix, modelRep );
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerFluidStructureInterfaceBC();
    this->clearMarkerRobinBC();

    std::string dirichletbcType = "elimination";//soption(_name="dirichletbc.type",_prefix=this->prefix());
    this->M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC( dirichletbcType, marker(d), ComponentType::NO_COMPONENT );
    for ( ComponentType comp : std::vector<ComponentType>( { ComponentType::X, ComponentType::Y, ComponentType::Z } ) )
    {
        std::string compTag = ( comp ==ComponentType::X )? "x" : (comp == ComponentType::Y )? "y" : "z";
        this->M_bcDirichletComponents[comp] = this->modelProperties().boundaryConditions().getScalarFields( (boost::format("displacement_%1%")%compTag).str(), "Dirichlet" );
        for( auto const& d : this->M_bcDirichletComponents.find(comp)->second )
            this->addMarkerDirichletBC( dirichletbcType, marker(d), comp );
    }

    this->M_bcNeumannScalar = this->modelProperties().boundaryConditions().getScalarFields( "displacement", "Neumann_scalar" );
    for( auto const& d : this->M_bcNeumannScalar )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d));
    this->M_bcNeumannVectorial = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "Neumann_vectorial" );
    for( auto const& d : this->M_bcNeumannVectorial )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::VECTORIAL,marker(d));
    this->M_bcNeumannTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<super_type::nDim>( "displacement", "Neumann_tensor2" );
    for( auto const& d : this->M_bcNeumannTensor2 )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::TENSOR2,marker(d));

    this->M_bcInterfaceFSI = this->modelProperties().boundaryConditions().getScalarFields( "displacement", "interface_fsi" );
    for( auto const& d : this->M_bcInterfaceFSI )
        this->addMarkerFluidStructureInterfaceBC( marker(d) );

    this->M_bcRobin = this->modelProperties().boundaryConditions().template getVectorFieldsList<super_type::nDim>( "displacement", "robin" );
    for( auto const& d : this->M_bcRobin )
        this->addMarkerRobinBC( marker(d) );

    this->M_bcNeumannEulerianFrameScalar = this->modelProperties().boundaryConditions().getScalarFields( { { "displacement", "Neumann_eulerian_scalar" },{ "displacement", "FollowerPressure" } } );
    for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
        this->addMarkerNeumannEulerianFrameBC(super_type::NeumannEulerianFrameBCShape::SCALAR,marker(d));
    this->M_bcNeumannEulerianFrameVectorial = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "Neumann_eulerian_vectorial" );
    for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
        this->addMarkerNeumannEulerianFrameBC(super_type::NeumannEulerianFrameBCShape::VECTORIAL,marker(d));
    this->M_bcNeumannEulerianFrameTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<super_type::nDim>( "displacement", "Neumann_eulerian_tensor2" );
    for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
        this->addMarkerNeumannEulerianFrameBC(super_type::NeumannEulerianFrameBCShape::TENSOR2,marker(d));

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "VolumicForces" );

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
    this->modelProperties().parameters().updateParameterValues();

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->M_bcDirichlet.setParameterValues( paramValues );
    for ( auto & bcDirComp : this->M_bcDirichletComponents )
        bcDirComp.second.setParameterValues( paramValues );
    this->M_bcNeumannScalar.setParameterValues( paramValues );
    this->M_bcNeumannVectorial.setParameterValues( paramValues );
    this->M_bcNeumannTensor2.setParameterValues( paramValues );
    this->M_bcNeumannEulerianFrameScalar.setParameterValues( paramValues );
    this->M_bcNeumannEulerianFrameVectorial.setParameterValues( paramValues );
    this->M_bcNeumannEulerianFrameTensor2.setParameterValues( paramValues );
    this->M_bcRobin.setParameterValues( paramValues );
    this->M_volumicForcesProperties.setParameterValues( paramValues );

    super_type::solve( upVelAcc );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletStrongJacobian( sparse_matrix_ptrtype& J, vector_ptrtype& RBis ) const
{
    if ( !this->hasDirichletBC() ) return;

    //auto RBis = this->backend()->newVector( J->mapRowPtr() );
    auto bilinearForm_PatternCoupled = form2( _test=this->functionSpace(),_trial=this->functionSpace(),_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldDisplacement();
    for( auto const& d : this->M_bcDirichlet )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) );
        auto const& listMarkerFaces = std::get<0>( ret );
        auto const& listMarkerEdges = std::get<1>( ret );
        auto const& listMarkerPoints = std::get<2>( ret );
        if ( !listMarkerFaces.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedfaces(this->mesh(), listMarkerFaces),
                    _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/,
                    _prefix=this->prefix() );
        if ( !listMarkerEdges.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markededges(this->mesh(), listMarkerEdges),
                    _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/,
                    _prefix=this->prefix() );
        if ( !listMarkerPoints.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedpoints(this->mesh(), listMarkerPoints),
                    _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/,
                    _prefix=this->prefix() );
    }
    for ( auto const& bcDirComp : this->M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),comp ) );
            auto const& listMarkerFaces = std::get<0>( ret );
            auto const& listMarkerEdges = std::get<1>( ret );
            auto const& listMarkerPoints = std::get<2>( ret );
            if ( !listMarkerFaces.empty() )
                bilinearForm_PatternCoupled +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=u[comp],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/,
                        _prefix=this->prefix() );
            if ( !listMarkerEdges.empty() )
                bilinearForm_PatternCoupled +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=u[comp],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/,
                        _prefix=this->prefix() );
            if ( !listMarkerPoints.empty() )
                bilinearForm_PatternCoupled +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=u[comp],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/,
                        _prefix=this->prefix() );

        }
    }
}
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletStrongLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    if ( this->isStandardModel() )
    {
        if ( !this->hasDirichletBC() ) return;

        auto Xh = this->functionSpace();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->fieldDisplacement();
        for( auto const& d : this->M_bcDirichlet )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) );
            auto const& listMarkerFaces = std::get<0>( ret );
            auto const& listMarkerEdges = std::get<1>( ret );
            auto const& listMarkerPoints = std::get<2>( ret );
            if ( !listMarkerFaces.empty() )
                bilinearForm +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=u,_rhs=F,_expr=expression(d),
                        _prefix=this->prefix() );
            if ( !listMarkerEdges.empty() )
                bilinearForm +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=u,_rhs=F,_expr=expression(d),
                        _prefix=this->prefix() );
            if ( !listMarkerPoints.empty() )
                bilinearForm +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=u,_rhs=F,_expr=expression(d),
                        _prefix=this->prefix() );
        }
        for ( auto const& bcDirComp : this->M_bcDirichletComponents )
        {
            ComponentType comp = bcDirComp.first;
            for( auto const& d : bcDirComp.second )
            {
                auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),comp ) );
                auto const& listMarkerFaces = std::get<0>( ret );
                auto const& listMarkerEdges = std::get<1>( ret );
                auto const& listMarkerPoints = std::get<2>( ret );
                if ( !listMarkerFaces.empty() )
                    bilinearForm +=
                        on( _range=markedfaces(this->mesh(), listMarkerFaces),
                            _element=u[comp],_rhs=F,_expr=expression(d),
                            _prefix=this->prefix() );
                if ( !listMarkerEdges.empty() )
                    bilinearForm +=
                        on( _range=markededges(this->mesh(), listMarkerEdges),
                            _element=u[comp],_rhs=F,_expr=expression(d),
                            _prefix=this->prefix() );
                if ( !listMarkerPoints.empty() )
                    bilinearForm +=
                        on( _range=markedpoints(this->mesh(), listMarkerPoints),
                            _element=u[comp],_rhs=F,_expr=expression(d),
                            _prefix=this->prefix() );
            }
        }
    }
    else if ( this->is1dReducedModel() )
    {
        if ( this->M_bcDirichlet.empty() ) return;

        auto Xh = this->functionSpace1dReduced();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->fieldDisplacementScal1dReduced();
        //WARNING : fixed at zero
        for( auto const& d : this->M_bcDirichlet )
            bilinearForm +=
                on( _range=markedfaces(Xh->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                    _element=u, _rhs=F, _expr=cst(0.),
                    _prefix=this->prefix() );
    }
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannResidual(vector_ptrtype& R) const
{
    if ( this->M_bcNeumannScalar.empty() && this->M_bcNeumannVectorial.empty() && this->M_bcNeumannTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_bcNeumannScalar )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d)) ),
                       _expr= -expression(d)*inner( N(),id(v) ),
                       _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannVectorial )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::VECTORIAL,marker(d)) ),
                       _expr= -inner( expression(d),id(v) ),
                       _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannTensor2 )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::TENSOR2,marker(d)) ),
                       _expr= -inner( expression(d)*N(),id(v) ),
                       _geomap=this->geomap() );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
    if ( this->M_bcNeumannScalar.empty() && this->M_bcNeumannVectorial.empty() && this->M_bcNeumannTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_bcNeumannScalar )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d)) ),
                       _expr= expression(d)*inner( N(),id(v) ),
                        _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannVectorial )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::VECTORIAL,marker(d)) ),
                       _expr= inner( expression(d),id(v) ),
                       _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannTensor2 )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::TENSOR2,marker(d)) ),
                       _expr= inner( expression(d)*N(),id(v) ),
                       _geomap=this->geomap() );
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinResidual(element_displacement_external_storage_type const& u, vector_ptrtype& R) const
{

    if ( this->M_bcRobin.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );

    // Warning : take only first component of expression1
    for( auto const& d : this->M_bcRobin )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)/*this->markerRobinBC()*/),
                       _expr= inner( expression1(d)(0,0)*idv(u) - expression2(d) ,id(u) ),
                       _geomap=this->geomap() );

}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinJacobian( sparse_matrix_ptrtype& J) const
{
    if ( this->M_bcRobin.empty() ) return;

    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldDisplacement();

    // Warning : take only first component of expression1
    for( auto const& d : this->M_bcRobin )
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)/*this->markerRobinBC()*/),
                       _expr= expression1(d)(0,0)*inner( idt(u) ,id(u) ),
                       _geomap=this->geomap() );

}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinLinearPDE( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const
{
    if ( this->M_bcRobin.empty() ) return;

    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& u = this->fieldDisplacement();

    // Warning : take only first component of expression1
    for( auto const& d : this->M_bcRobin )
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

    if ( this->M_volumicForcesProperties.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_volumicForcesProperties )
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
    if ( this->M_volumicForcesProperties.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_volumicForcesProperties )
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
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCFollowerPressureResidual( typename super_type::element_displacement_external_storage_type const& u, vector_ptrtype& R ) const
{
    if ( this->M_bcNeumannEulerianFrameScalar.empty() && this->M_bcNeumannEulerianFrameVectorial.empty() && this->M_bcNeumannEulerianFrameTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)),
                       _expr= -expression(d)*inner( Feel::vf::FeelModels::solidMecGeomapEulerian(u)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)),
                       _expr= -inner( Feel::vf::FeelModels::solidMecGeomapEulerian(u)*expression(d),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)),
                       _expr= -inner( Feel::vf::FeelModels::solidMecGeomapEulerian(u)*expression(d)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCFollowerPressureJacobian( typename super_type::element_displacement_external_storage_type const& u, sparse_matrix_ptrtype& J) const
{
    if ( this->M_bcNeumannEulerianFrameScalar.empty() && this->M_bcNeumannEulerianFrameVectorial.empty() && this->M_bcNeumannEulerianFrameTensor2.empty() ) return;

    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );

    for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)) ,
                       _expr= -expression(d)*inner(Feel::vf::FeelModels::solidMecGeomapEulerianJacobian(u)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)) ,
                       _expr= -inner(Feel::vf::FeelModels::solidMecGeomapEulerianJacobian(u)*expression(d),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)) ,
                       _expr= -inner(Feel::vf::FeelModels::solidMecGeomapEulerianJacobian(u)*expression(d)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
}


} // namespace FeelModels
} // namespace Feel
