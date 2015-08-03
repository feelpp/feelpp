/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-07-17

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file codegen_solidmec.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-07-17
 */

#include "codegen_solidmec.hpp"

#include <feel/feelfilters/geotool.hpp>

#include <feel/feelmodels/modelvf/solidmecgeomapeulerian.hpp>
#include <feel/feelmodels/modelalg/functionSup.cpp>

#undef SOLIDMECHANICS_MESH
#undef SOLIDMECHANICS_1D_REDUCED_MESH
#include "solid.mesh"

namespace Feel
{
namespace FeelModels
{


SOLIDMECHANICS_CLASS_NAME::SOLIDMECHANICS_CLASS_NAME( std::string __prefix,
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


void
SOLIDMECHANICS_CLASS_NAME::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerFluidStructureInterfaceBC();
    this->clearMarkerRobinBC();

    auto const bcDef = SOLIDMECHANICS_BC(this/*this->shared_from_this()*/);
    std::string dirichletbcType = "elimination";//soption(_name="dirichletbc.type",_prefix=this->prefix());
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_vec>() )
        this->addMarkerDirichletBC( dirichletbcType, PhysicalName, ComponentType::NO_COMPONENT );
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_x>() )
        this->addMarkerDirichletBC( dirichletbcType, PhysicalName, ComponentType::X );
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_y>() )
        this->addMarkerDirichletBC( dirichletbcType, PhysicalName, ComponentType::Y );
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_z>() )
        this->addMarkerDirichletBC( dirichletbcType, PhysicalName, ComponentType::Z );
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::neumann_scal>() )
        this->addMarkerNeumannBC(NeumannBCShape::SCALAR,PhysicalName);
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::neumann_vec>() )
        this->addMarkerNeumannBC(NeumannBCShape::VECTORIAL,PhysicalName);
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::paroi_mobile>() )
        this->addMarkerFluidStructureInterfaceBC( PhysicalName );
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::robin_vec>() )
        this->addMarkerRobinBC( PhysicalName );
}

//---------------------------------------------------------------------------------------------------//

void
SOLIDMECHANICS_CLASS_NAME::loadConfigMeshFile( std::string const& geofilename )
{

#if !defined(SOLIDMECHANICS_MESH) && !defined(SOLIDMECHANICS_MESH1) && !defined(SOLIDMECHANICS_MESH2) && !defined(SOLIDMECHANICS_MESH3)
        CHECK( false ) << " SOLIDMECHANICS_MESH is not defined : probably .mesh is wrong\n";
#endif

        switch ( this->geotoolMeshIndex() )
        {
        default :
        case 0 :
        {
#if defined(SOLIDMECHANICS_MESH)
            SOLIDMECHANICS_MESH(M_meshSize, geofilename );
            M_mesh=mesh;
#else
            CHECK( false ) << "geometry not define with MeshIndex=" << this->geotoolMeshIndex() << "\n";
#endif
        }
        break;
        case 1 :
        {
#if defined(SOLIDMECHANICS_MESH1)
            SOLIDMECHANICS_MESH1(M_meshSize, geofilename );
            M_mesh=mesh;
#else
            CHECK( false ) << "geometry not define with MeshIndex=" << this->geotoolMeshIndex() << "\n";
#endif
        }
        break;
        case 2 :
        {
#if defined(SOLIDMECHANICS_MESH2)
            SOLIDMECHANICS_MESH2(M_meshSize, geofilename );
            M_mesh=mesh;
#else
            CHECK( false ) << "geometry not define with MeshIndex=" << this->geotoolMeshIndex() << "\n";
#endif
        }
        break;
        case 3 :
        {
#if defined(SOLIDMECHANICS_MESH3)
            SOLIDMECHANICS_MESH3(M_meshSize, geofilename );
            M_mesh=mesh;
#else
            CHECK( false ) << "geometry not define with MeshIndex=" << this->geotoolMeshIndex() << "\n";
#endif
        }
        break;

        } // switch

}

//---------------------------------------------------------------------------------------------------//
void
SOLIDMECHANICS_CLASS_NAME::loadConfigMeshFile1dReduced( std::string const& geofilename )
{
#if defined( SOLIDMECHANICS_1D_REDUCED_MESH )
    SOLIDMECHANICS_1D_REDUCED_MESH(M_meshSize,geofilename);
    M_mesh_1d_reduced = mesh;
#else
    CHECK( false ) << " SOLIDMECHANICS_1D_REDUCED_MESH not define\n";
#endif

}

//---------------------------------------------------------------------------------------------------//

void
SOLIDMECHANICS_CLASS_NAME::init( bool buildMethodNum )
{
    super_type::init( buildMethodNum, this->shared_from_this() );
}

//---------------------------------------------------------------------------------------------------//

void
SOLIDMECHANICS_CLASS_NAME::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
    auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
    if ( !bcDef.hasDirichletVec() && !bcDef.hasDirichletX() && !bcDef.hasDirichletY() && !bcDef.hasDirichletZ() )
        return;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateNewtonInitialGuess", "start",
                                        this->worldComm(),this->verboseAllProc());

    auto const& u = this->fieldDisplacement();
    auto Xh = this->functionSpace();
    auto mesh = this->mesh();
    auto const rowStartInVector = this->rowStartInVector();

    if (Xh->worldsComm()[0].isActive()) // only on Displacement Proc
    {
        ForEachBC( bcDef,cl::dirichlet_vec,
                   modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",PhysicalName ) ),
                            u, U, Expression,rowStartInVector ) );
        ForEachBC( bcDef,cl::dirichlet_x,
                   modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",PhysicalName, Component::X ) ),
                            u[Component::X], U, Expression,rowStartInVector, element_displacement_type::nComponents ) );
        ForEachBC( bcDef,cl::dirichlet_y,
                   modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",PhysicalName, Component::Y ) ),
                            u[Component::Y], U, Expression,rowStartInVector, element_displacement_type::nComponents ) );
        ForEachBC( bcDef,cl::dirichlet_z,
                   modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",PhysicalName, Component::Z ) ),
                            u[Component::Z], U, Expression,rowStartInVector, element_displacement_type::nComponents ) );
    }
    U->close();

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateNewtonInitialGuess", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

void
SOLIDMECHANICS_CLASS_NAME::updateBCDirichletStrongResidual(vector_ptrtype& R) const
{
    auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
    if ( !bcDef.hasDirichletVec() && !bcDef.hasDirichletX() && !bcDef.hasDirichletY() && !bcDef.hasDirichletZ() )
        return;

    auto const& u = this->fieldDisplacement();
    auto const rowStartInVector = this->rowStartInVector();

    R->close();
    if (M_Xh->worldsComm()[0].isActive()) // only on Displacement Proc
    {
        for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_vec>() )
            modifVec( markedfaces(this->mesh(),this->markerDirichletBCByNameId( "elimination",PhysicalName ) ),
                      u, R, 0*vf::one(), rowStartInVector );
        for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_x>() )
            modifVec( markedfaces(this->mesh(),this->markerDirichletBCByNameId( "elimination",PhysicalName, Component::X ) ),
                      u[Component::X], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
        for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_y>() )
            modifVec( markedfaces(this->mesh(),this->markerDirichletBCByNameId( "elimination",PhysicalName, Component::Y ) ),
                      u[Component::Y], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
        for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_z>() )
            modifVec( markedfaces(this->mesh(),this->markerDirichletBCByNameId( "elimination",PhysicalName, Component::Z ) ),
                      u[Component::Z], R, vf::cst(0.), rowStartInVector, element_displacement_type::nComponents );
    }
}
void
SOLIDMECHANICS_CLASS_NAME::updateBCDirichletStrongJacobian( sparse_matrix_ptrtype& J ) const
{
    auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
    if ( !bcDef.hasDirichletVec() && !bcDef.hasDirichletX() && !bcDef.hasDirichletY() && !bcDef.hasDirichletZ() )
        return;

    auto RBis = M_backend->newVector( J->mapRowPtr() );
    auto bilinearForm_PatternCoupled = form2( _test=this->functionSpace(),_trial=this->functionSpace(),_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldDisplacement();
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_vec>() )
        bilinearForm_PatternCoupled +=
            /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",PhysicalName ) ),
                     _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/ );
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_x>() )
        bilinearForm_PatternCoupled +=
            /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",PhysicalName,Component::X ) ),
                     _element=u[Component::X],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_y>() )
        bilinearForm_PatternCoupled +=
            /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",PhysicalName,Component::Y ) ),
                     _element=u[Component::Y],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_z>() )
        bilinearForm_PatternCoupled +=
            /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",PhysicalName,Component::Z ) ),
                     _element=u[Component::Z],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/ );

}
void
SOLIDMECHANICS_CLASS_NAME::updateBCDirichletStrongLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
    if ( !bcDef.hasDirichletVec() && !bcDef.hasDirichletX() && !bcDef.hasDirichletY() && !bcDef.hasDirichletZ() )
        return;

    if ( this->isStandardModel() )
    {
        auto Xh = this->functionSpace();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->fieldDisplacement();
        ForEachBC( bcDef,cl::dirichlet_vec,
                   bilinearForm +=
                   /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",PhysicalName ) ),
                            _element=u,_rhs=F,_expr=Expression ); );
        ForEachBC( bcDef,cl::dirichlet_x,
                   bilinearForm +=
                   /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",PhysicalName,Component::X ) ),
                            _element=u[Component::X],_rhs=F,_expr=Expression ); );
        ForEachBC( bcDef,cl::dirichlet_y,
                   bilinearForm +=
                   /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",PhysicalName,Component::Y ) ),
                            _element=u[Component::Y],_rhs=F,_expr=Expression ); );
        ForEachBC( bcDef,cl::dirichlet_z,
                   bilinearForm +=
                   /**/ on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",PhysicalName,Component::Z ) ),
                            _element=u[Component::Z],_rhs=F,_expr=Expression ); );
    }
    else if ( this->is1dReducedModel() )
    {
        auto Xh = this->functionSpace1dReduced();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->getDisplacementScal1dReduced();
        //WARNING : fixed at zero
        ForEachBC( bcDef,cl::dirichlet_vec,
                   bilinearForm +=
                   /**/ on( _range=markedfaces(Xh->mesh(),this->markerDirichletBCByNameId( "elimination",PhysicalName ) ),
                            _element=u, _rhs=F,
                            _expr=cst(0.) );
                   );
    }

}


void
SOLIDMECHANICS_CLASS_NAME::updateBCNeumannResidual(vector_ptrtype& R) const
{
    auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
    if ( !bcDef.hasNeumannScal() && !bcDef.hasNeumannVec() )
        return;

    auto const& v = this->fieldDisplacement();

    double alpha_f=M_genAlpha_alpha_f;

    if ( bcDef.hasNeumannScal() )
    {
        ForEachBC( bcDef,cl::neumann_scal,
                   form1( _test=this->functionSpace(), _vector=R) +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::SCALAR,PhysicalName) ),
                                   _expr= -alpha_f*Expression*inner( N(),id(v) ),
                                   _geomap=this->geomap() ) );
    }
    if ( bcDef.hasNeumannVec() )
    {
        ForEachBC( bcDef,cl::neumann_vec,
                   form1( _test=this->functionSpace(), _vector=R) +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::VECTORIAL,PhysicalName) ),
                                   _expr= -alpha_f*inner( Expression, id(v) ),
                                   _geomap=this->geomap() ) );
    }
}

void
SOLIDMECHANICS_CLASS_NAME::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
   auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
    if ( !bcDef.hasNeumannScal() && !bcDef.hasNeumannVec() )
        return;

    auto const& v = this->fieldDisplacement();

    double alpha_f=M_genAlpha_alpha_f;

    if ( bcDef.hasNeumannScal() )
    {
        ForEachBC( bcDef,cl::neumann_scal,
                   form1( _test=this->functionSpace(), _vector=F) +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::SCALAR,PhysicalName) ),
                                   _expr= alpha_f*Expression*inner( N(),id(v) ),
                                   _geomap=this->geomap() ) );
    }
    if ( bcDef.hasNeumannVec() )
    {
        ForEachBC( bcDef,cl::neumann_vec,
                   form1( _test=this->functionSpace(), _vector=F) +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::VECTORIAL,PhysicalName) ),
                                   _expr= alpha_f*inner( Expression,id(v) ),
                                   _geomap=this->geomap() ) );
    }
}


void
SOLIDMECHANICS_CLASS_NAME::updateBCRobinResidual(element_displacement_type const& u, vector_ptrtype& R) const
{
    CHECK( false ) << "TODO";
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
}

void
SOLIDMECHANICS_CLASS_NAME::updateBCRobinJacobian( sparse_matrix_ptrtype& J) const
{
    CHECK( false ) << "TODO";
}

void
SOLIDMECHANICS_CLASS_NAME::updateBCRobinLinearPDE( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const
{
    CHECK( false ) << "TODO";
}


void
SOLIDMECHANICS_CLASS_NAME::updateSourceTermResidual( vector_ptrtype& R ) const
{
#if defined(SOLIDMECHANICS_VOLUME_FORCE)
    auto f = SOLIDMECHANICS_VOLUME_FORCE(this->shared_from_this());
    auto const& v = this->fieldDisplacement();

    form1( _test=this->functionSpace(), _vector=R ) +=
        integrate( _range=elements(this->mesh()),
                   _expr= -trans(f)*id(v),
                   _geomap=this->geomap() );
#endif
}

void
SOLIDMECHANICS_CLASS_NAME::updateSourceTermLinearPDE( vector_ptrtype& F ) const
{
#if defined(SOLIDMECHANICS_VOLUME_FORCE)
    auto f = SOLIDMECHANICS_VOLUME_FORCE(this->shared_from_this());
    auto const& v = this->fieldDisplacement();

    form1( _test=this->functionSpace(), _vector=F ) +=
        integrate( _range=elements(this->mesh()),
                   _expr= trans(f)*id(v),
                   _geomap=this->geomap() );
#endif
}

void
SOLIDMECHANICS_CLASS_NAME::updateBCFollowerPressureResidual( element_displacement_type const& u, vector_ptrtype& R ) const
{
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
}

void
SOLIDMECHANICS_CLASS_NAME::updateBCFollowerPressureJacobian( element_displacement_type const& u, sparse_matrix_ptrtype& J) const
{
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
}




} // namespace FeelModels

} // namespace Feel
