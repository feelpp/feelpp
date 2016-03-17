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
#include <feel/feelmodels/modelalg/functionSup.cpp>

namespace Feel
{
namespace FeelModels
{

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::SolidMechanics( std::string const& prefix,
                                                    bool buildMesh,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    std::string const& rootRepository )
    :
    super_type( prefix, buildMesh, worldComm, subPrefix, rootRepository )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());

    //-----------------------------------------------------------------------------//
    // load info from .json file
    this->loadConfigBCFile();
    this->loadConfigPostProcess();
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
                                         std::string const& appliShortRepository )
{
    return boost::make_shared<self_type>( prefix, buildMesh, worldComm, subPrefix, appliShortRepository );
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
    M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "Dirichlet" );
    for( auto const& d : M_bcDirichlet )
        this->addMarkerDirichletBC( dirichletbcType, marker(d), ComponentType::NO_COMPONENT );
    for ( ComponentType comp : std::vector<ComponentType>( { ComponentType::X, ComponentType::Y, ComponentType::Z } ) )
    {
        std::string compTag = ( comp ==ComponentType::X )? "x" : (comp == ComponentType::Y )? "y" : "z";
        M_bcDirichletComponents[comp] = this->modelProperties().boundaryConditions().getScalarFields( (boost::format("displacement_%1%")%compTag).str(), "Dirichlet" );
        for( auto const& d : M_bcDirichletComponents.find(comp)->second )
            this->addMarkerDirichletBC( dirichletbcType, marker(d), comp );
    }

    M_bcNeumannScalar = this->modelProperties().boundaryConditions().getScalarFields( "displacement", "Neumann_scalar" );
    for( auto const& d : M_bcNeumannScalar )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d));
    M_bcNeumannVectorial = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "Neumann_vectorial" );
    for( auto const& d : M_bcNeumannVectorial )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::VECTORIAL,marker(d));
    M_bcNeumannTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<super_type::nDim>( "displacement", "Neumann_tensor2" );
    for( auto const& d : M_bcNeumannTensor2 )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::TENSOR2,marker(d));

    M_bcInterfaceFSI = this->modelProperties().boundaryConditions().getScalarFields( "displacement", "interface_fsi" );
    for( auto const& d : M_bcInterfaceFSI )
        this->addMarkerFluidStructureInterfaceBC( marker(d) );

    M_bcRobin = this->modelProperties().boundaryConditions().template getVectorFieldsList<super_type::nDim>( "displacement", "robin" );
    for( auto const& d : M_bcRobin )
        this->addMarkerRobinBC( marker(d) );

    M_bcNeumannEulerianFrameScalar = this->modelProperties().boundaryConditions().getScalarFields( { { "displacement", "Neumann_eulerian_scalar" },{ "displacement", "FollowerPressure" } } );
    for( auto const& d : M_bcNeumannEulerianFrameScalar )
        this->addMarkerNeumannEulerianFrameBC(super_type::NeumannEulerianFrameBCShape::SCALAR,marker(d));
    M_bcNeumannEulerianFrameVectorial = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "Neumann_eulerian_vectorial" );
    for( auto const& d : M_bcNeumannEulerianFrameVectorial )
        this->addMarkerNeumannEulerianFrameBC(super_type::NeumannEulerianFrameBCShape::VECTORIAL,marker(d));
    M_bcNeumannEulerianFrameTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<super_type::nDim>( "displacement", "Neumann_eulerian_tensor2" );
    for( auto const& d : M_bcNeumannEulerianFrameTensor2 )
        this->addMarkerNeumannEulerianFrameBC(super_type::NeumannEulerianFrameBCShape::TENSOR2,marker(d));

    M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "VolumicForces" );

}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::loadConfigPostProcess()
{
    if ( this->modelProperties().postProcess().find("Fields") != this->modelProperties().postProcess().end() )
        for ( auto const& o : this->modelProperties().postProcess().find("Fields")->second )
        {
            if ( o == "displacement" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Displacement );
            if ( o == "velocity" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Velocity );
            if ( o == "acceleration" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Acceleration );
            if ( o == "stress" || o == "normal-stress" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::NormalStress );
            if ( o == "pressure" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Pressure );
            if ( o == "material-properties" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::MaterialProperties );
            if ( o == "pid" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Pid );
            if ( o == "fsi" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::FSI );
            if ( o == "Von-Mises" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::VonMises );
            if ( o == "Tresca" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::Tresca );
            if ( o == "principal-stresses" || o == "all" ) this->M_postProcessFieldExported.insert( SolidMechanicsPostProcessFieldExported::PrincipalStresses );
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
    this->modelProperties().parameters().updateParameterValues();

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    M_bcDirichlet.setParameterValues( paramValues );
    for ( auto & bcDirComp : M_bcDirichletComponents )
        bcDirComp.second.setParameterValues( paramValues );
    M_bcNeumannScalar.setParameterValues( paramValues );
    M_bcNeumannVectorial.setParameterValues( paramValues );
    M_bcNeumannTensor2.setParameterValues( paramValues );
    M_bcNeumannEulerianFrameScalar.setParameterValues( paramValues );
    M_bcNeumannEulerianFrameVectorial.setParameterValues( paramValues );
    M_bcNeumannEulerianFrameTensor2.setParameterValues( paramValues );
    M_bcRobin.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );

    super_type::solve( upVelAcc );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
    if ( !this->hasDirichletBC() ) return;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateNewtonInitialGuess", "start",
                                               this->worldComm(),this->verboseAllProc());

    auto Xh = this->functionSpace();
    auto u = Xh->element( U,0 );
    auto mesh = this->mesh();
    size_type rowStartInVector = this->rowStartInVector();

    if ( !Xh->worldsComm()[0].isActive()) // only on Displacement Proc
        return;

    for( auto const& d : M_bcDirichlet )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d) ) );
        auto const& listMarkerFaces = std::get<0>( ret );
        auto const& listMarkerEdges = std::get<1>( ret );
        auto const& listMarkerPoints = std::get<2>( ret );
        if ( !listMarkerFaces.empty() )
            u.on(_range=markedfaces(mesh,listMarkerFaces ),
                 _expr=expression(d) );
        if ( !listMarkerEdges.empty() )
            u.on(_range=markededges(mesh,listMarkerEdges),
                 _expr=expression(d) );
        if ( !listMarkerPoints.empty() )
            u.on(_range=markedpoints(mesh,listMarkerPoints),
                 _expr=expression(d) );
    }
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d), comp ) );
            auto const& listMarkerFaces = std::get<0>( ret );
            auto const& listMarkerEdges = std::get<1>( ret );
            auto const& listMarkerPoints = std::get<2>( ret );
            if ( !listMarkerFaces.empty() )
                u[comp].on(_range=markedfaces(mesh,listMarkerFaces ),
                           _expr=expression(d) );
            if ( !listMarkerEdges.empty() )
                u[comp].on(_range=markededges(mesh,listMarkerEdges),
                           _expr=expression(d) );
            if ( !listMarkerPoints.empty() )
                u[comp].on(_range=markedpoints(mesh,listMarkerPoints),
                           _expr=expression(d) );
        }
    }

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateNewtonInitialGuess", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletStrongResidual(vector_ptrtype& R) const
{
    if ( !this->hasDirichletBC() ) return;

    size_type rowStartInVector = this->rowStartInVector();
    auto u = this->functionSpace()->element( R,0 );

    R->close();
    if ( !this->functionSpaceDisplacement()->worldsComm()[0].isActive() ) // only on Displacement Proc
        return;
    for( auto const& d : M_bcDirichlet )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) );
        auto const& listMarkerFaces = std::get<0>( ret );
        auto const& listMarkerEdges = std::get<1>( ret );
        auto const& listMarkerPoints = std::get<2>( ret );

        auto exprUsed = vf::zero<super_type::nDim,1>();// 0*vf::one();
        if ( !listMarkerFaces.empty() )
            u.on(_range=markedfaces(this->mesh(),listMarkerFaces),
                 _expr=exprUsed );
        if ( !listMarkerEdges.empty() )
            u.on(_range=markededges(this->mesh(),listMarkerEdges),
                 _expr=exprUsed );
        if ( !listMarkerPoints.empty() )
            u.on(_range=markedpoints(this->mesh(),listMarkerPoints),
                 _expr=exprUsed );
    }
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d), comp ) );
            auto const& listMarkerFaces = std::get<0>( ret );
            auto const& listMarkerEdges = std::get<1>( ret );
            auto const& listMarkerPoints = std::get<2>( ret );
            auto exprUsed = vf::cst(0.);
            if ( !listMarkerFaces.empty() )
                u[comp].on(_range=markedfaces(this->mesh(),listMarkerFaces),
                           _expr=exprUsed );
            if ( !listMarkerEdges.empty() )
                u[comp].on(_range=markededges(this->mesh(),listMarkerEdges),
                           _expr=exprUsed );
            if ( !listMarkerPoints.empty() )
                u[comp].on(_range=markedpoints(this->mesh(),listMarkerPoints),
                           _expr=exprUsed );
        }
    }

}
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
    for( auto const& d : M_bcDirichlet )
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
    for ( auto const& bcDirComp : M_bcDirichletComponents )
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
        for( auto const& d : M_bcDirichlet )
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
        for ( auto const& bcDirComp : M_bcDirichletComponents )
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
                    _element=u, _rhs=F, _expr=cst(0.),
                    _prefix=this->prefix() );
    }
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannResidual(vector_ptrtype& R) const
{
    if ( M_bcNeumannScalar.empty() && M_bcNeumannVectorial.empty() && M_bcNeumannTensor2.empty() ) return;

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
    for( auto const& d : M_bcNeumannTensor2 )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::TENSOR2,marker(d)) ),
                       _expr= -inner( expression(d)*N(),id(v) ),
                       _geomap=this->geomap() );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
    if ( M_bcNeumannScalar.empty() && M_bcNeumannVectorial.empty() && M_bcNeumannTensor2.empty() ) return;

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
    for( auto const& d : M_bcNeumannTensor2 )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::TENSOR2,marker(d)) ),
                       _expr= inner( expression(d)*N(),id(v) ),
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
    if ( M_bcNeumannEulerianFrameScalar.empty() && M_bcNeumannEulerianFrameVectorial.empty() && M_bcNeumannEulerianFrameTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    for( auto const& d : M_bcNeumannEulerianFrameScalar )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)),
                       _expr= -expression(d)*inner( Feel::vf::FeelModels::solidMecGeomapEulerian(u)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : M_bcNeumannEulerianFrameVectorial )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)),
                       _expr= -inner( Feel::vf::FeelModels::solidMecGeomapEulerian(u)*expression(d),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : M_bcNeumannEulerianFrameTensor2 )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)),
                       _expr= -inner( Feel::vf::FeelModels::solidMecGeomapEulerian(u)*expression(d)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCFollowerPressureJacobian( typename super_type::element_displacement_type const& u, sparse_matrix_ptrtype& J) const
{
    if ( M_bcNeumannEulerianFrameScalar.empty() && M_bcNeumannEulerianFrameVectorial.empty() && M_bcNeumannEulerianFrameTensor2.empty() ) return;

    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );

    for( auto const& d : M_bcNeumannEulerianFrameScalar )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)) ,
                       _expr= -expression(d)*inner(Feel::vf::FeelModels::solidMecGeomapEulerianJacobian(u)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : M_bcNeumannEulerianFrameVectorial )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)) ,
                       _expr= -inner(Feel::vf::FeelModels::solidMecGeomapEulerianJacobian(u)*expression(d),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : M_bcNeumannEulerianFrameTensor2 )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)) ,
                       _expr= -inner(Feel::vf::FeelModels::solidMecGeomapEulerianJacobian(u)*expression(d)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
}


} // namespace FeelModels
} // namespace Feel
