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

#include <feel/feelmodels2/solid/solidmechanics.hpp>

#include <feel/feelmodels2/modelvf/solidmecgeomapeulerian.hpp>
#include <feel/feelmodels2/modelalg/functionSup.cpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::SolidMechanics( bool __isStationary,
                                                                     std::string __prefix,
                                                                     WorldComm const& __worldComm,
                                                                     bool __buildMesh,
                                                                     std::string __subPrefix,
                                                                     std::string __appliShortRepository )
    :
    super_type(__isStationary,__prefix,__buildMesh,__worldComm,__subPrefix, __appliShortRepository)
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


template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->M_markerNameFSI.clear();
    this->M_markerNameBCRobin.clear();

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

    M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "displacement", "VolumicForces" );

#if 0 // TODO
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::paroi_mobile>() )
        M_markerNameFSI.push_back(PhysicalName);
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::robin_vec>() )
        M_markerNameBCRobin.push_back(PhysicalName);
#endif
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::loadConfigMeshFile( std::string const& geofilename )
{
    CHECK( false ) << "not allow";
}

//---------------------------------------------------------------------------------------------------//
template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::loadConfigMeshFile1dReduced( std::string const& geofilename )
{
    CHECK( false ) << "not allow";
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::init( bool buildModelAlgebraicFactory )
{
    super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::solve( bool upVelAcc )
{
    M_bcDirichlet.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcDirichletX.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcDirichletY.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcDirichletZ.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcNeumannScalar.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcNeumannVectorial.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_volumicForcesProperties.setParameterValues( this->modelProperties().parameters().toParameterValues() );

    super_type::solve( upVelAcc );
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateCLDirichlet( vector_ptrtype& U ) const
{
    if ( M_bcDirichlet.empty() && M_bcDirichletX.empty() && M_bcDirichletY.empty() && M_bcDirichletZ.empty() ) return;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateCLDirichlet", "start",
                                               this->worldComm(),this->verboseAllProc());

    auto const& u = this->fieldDisplacement();
    auto Xh = this->functionSpace();
    auto mesh = this->mesh();
    size_type rowStartInVector = this->rowStartInVector();

    if (Xh->worldsComm()[0].isActive()) // only on Displacement Proc
    {
        for( auto const& d : M_bcDirichlet )
            modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                     u, U, expression(d), rowStartInVector );
        for( auto const& d : M_bcDirichletX )
            modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",marker(d), Component::X ) ),
                     u[Component::X], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
        for( auto const& d : M_bcDirichletY )
            modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",marker(d), Component::Y ) ),
                     u[Component::Y], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
        for( auto const& d : M_bcDirichletZ )
            modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",marker(d), Component::Z ) ),
                     u[Component::Z], U, expression(d),rowStartInVector, element_displacement_type::nComponents );
    }
    U->close();

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateCLDirichlet", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateBCDirichletStrongResidual(vector_ptrtype& R) const
{
    if ( M_bcDirichlet.empty() && M_bcDirichletX.empty() && M_bcDirichletY.empty() && M_bcDirichletZ.empty() ) return;

    auto const& u = this->fieldDisplacement();
    size_type rowStartInVector = this->rowStartInVector();

    R->close();
    if ( this->functionSpaceDisplacement()->worldsComm()[0].isActive() ) // only on Displacement Proc
    {
        for( auto const& d : M_bcDirichlet )
            modifVec( markedfaces(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                      u, R, 0*vf::one(), rowStartInVector );
        for( auto const& d : M_bcDirichletX )
            modifVec( markedfaces(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d), Component::X ) ),
                      u[Component::X], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
        for( auto const& d : M_bcDirichletY )
            modifVec( markedfaces(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d), Component::Y ) ),
                      u[Component::Y], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
        for( auto const& d : M_bcDirichletZ )
            modifVec( markedfaces(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d), Component::Z ) ),
                      u[Component::Z], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
    }
}
template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateBCDirichletStrongJacobian( sparse_matrix_ptrtype& J ) const
{
    if ( M_bcDirichlet.empty() && M_bcDirichletX.empty() && M_bcDirichletY.empty() && M_bcDirichletZ.empty() ) return;

    auto RBis = this->backend()->newVector( J->mapRowPtr() );
    auto bilinearForm_PatternCoupled = form2( _test=this->functionSpace(),_trial=this->functionSpace(),_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldDisplacement();
    for( auto const& d : M_bcDirichlet )
        bilinearForm_PatternCoupled +=
            /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                     _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/ );
    for( auto const& d : M_bcDirichletX )
        bilinearForm_PatternCoupled +=
            /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",marker(d),Component::X ) ),
                     _element=u[Component::X],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
    for( auto const& d : M_bcDirichletY )
        bilinearForm_PatternCoupled +=
            /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",marker(d),Component::Y ) ),
                     _element=u[Component::Y],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
    for( auto const& d : M_bcDirichletZ )
        bilinearForm_PatternCoupled +=
            /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",marker(d),Component::Z ) ),
                     _element=u[Component::Z],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
}
template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateBCDirichletStrongLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
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
            bilinearForm +=
                on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                    _element=u,_rhs=F,_expr=expression(d) );
        for( auto const& d : M_bcDirichletX )
            bilinearForm +=
                on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",marker(d),Component::X ) ),
                    _element=u[Component::X],_rhs=F,_expr=expression(d) );
        for( auto const& d : M_bcDirichletY )
            bilinearForm +=
                on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",marker(d),Component::Y ) ),
                    _element=u[Component::Y],_rhs=F,_expr=expression(d) );
        for( auto const& d : M_bcDirichletZ )
            bilinearForm +=
                on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",marker(d),Component::Z ) ),
                    _element=u[Component::Z],_rhs=F,_expr=expression(d) );
    }
    else if ( this->is1dReducedModel() )
    {
        if ( M_bcDirichlet.empty() ) return;

        auto Xh = this->functionSpace1dReduced();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->getDisplacementScal1dReduced();
        //WARNING : fixed at zero
        for( auto const& d : M_bcDirichlet )
            bilinearForm +=
                on( _range=markedfaces(Xh->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                    _element=u, _rhs=F,
                    _expr=cst(0.) );
    }
}


template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateBCNeumannResidual(vector_ptrtype& R) const
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

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
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


template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateBCRobinResidual(vector_ptrtype& R) const
{
#if 0 // TODO
#if 0
    auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
    if ( bcDef.hasRobinVec() )
    {
        auto const& v = this->fieldDisplacement();
        double alpha_robin = 1e4;
        ForEachBC( bcDef,cl::robin_vec,
                   form1( _test=this->functionSpace(), _vector=R) +=
                   /**/ integrate( _range=markedfaces(this->mesh(),PhysicalName),
                                   _expr= alpha_robin*trans(idv(u))*id(v),
                                   _geomap=this->geomap() ) );
    }
#endif
#endif
}


template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateSourceTermResidual( vector_ptrtype& R ) const
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

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateSourceTermLinearPDE( vector_ptrtype& F ) const
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

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateBCFollowerPressureResidual( typename super_type::element_displacement_type const& u, vector_ptrtype& R ) const
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

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanics<ConvexType,OrderDisp,UseCstMechProp>::updateBCFollowerPressureJacobian( typename super_type::element_displacement_type const& u, sparse_matrix_ptrtype& J) const
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
