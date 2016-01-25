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
 \file codegen_fluidmec.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-17
 */

#include "codegen_fluidmec.hpp"

#include <feel/feelfilters/geotool.hpp>
#include <feel/feelmodels/modelalg/functionSup.cpp>


#undef FLUIDMECHANICS_MESH
#undef FLUIDMECHANICS_MESH1
#undef FLUIDMECHANICS_MESH2
#undef FLUIDMECHANICS_MESH3
#include "fluid.mesh"

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_NAME::FLUIDMECHANICS_CLASS_NAME( std::string __prefix,
                                                      bool __buildMesh,
                                                      WorldComm const& __worldComm,
                                                      std::string __subPrefix,
                                                      std::string __appliShortRepository )
    :
    super_type( __prefix,__worldComm, __buildMesh, __subPrefix,__appliShortRepository)
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"FluidMechanics.info") );
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

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

void
FLUIDMECHANICS_CLASS_NAME::loadConfigBCFile()
{
    // use ALE mode or not
    auto const bcDef = FLUIDMECHANICS_BC(this/*->shared_from_this()*/); // not shared_from_this because call in constructor
    M_isMoveDomain = bcDef.hasMovingBoundary();
    M_hasFluidOutlet = bcDef.hasFluidOutlet();
    // clear
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerALEMeshBC();
    this->clearMarkerSlipBC();
    this->clearMarkerPressureBC();
    M_fluidOutletsBCType.clear();
    // boundary conditions
    std::string dirichletbcType = soption(_name="dirichletbc.type",_prefix=this->prefix());
    CHECK( dirichletbcType=="elimination" || dirichletbcType=="nitsche" || dirichletbcType=="lm" ) << "invalid dirichletbc.type " << dirichletbcType;
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::dirichlet_vec>() )
    {
        this->addMarkerDirichletBC( dirichletbcType, PhysicalName, ComponentType::NO_COMPONENT );
        this->addMarkerALEMeshBC("fixed",PhysicalName);
    }
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::paroi_mobile>() )
    {
        this->addMarkerDirichletBC( dirichletbcType, PhysicalName, ComponentType::NO_COMPONENT );
        this->addMarkerALEMeshBC("moving",PhysicalName);
    }
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::neumann_scal>() )
    {
        this->addMarkerNeumannBC(NeumannBCShape::SCALAR,PhysicalName);
        this->addMarkerALEMeshBC("fixed",PhysicalName);
    }
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::pressure>() )
    {
        this->addMarkerPressureBC(PhysicalName);
        this->addMarkerALEMeshBC("fixed",PhysicalName);
    }
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::slip>() )
    {
        this->addMarkerSlipBC(PhysicalName);
        this->addMarkerALEMeshBC("fixed",PhysicalName);
    }

    M_fluidOutletType = soption(_name="fluid-outlet.type", _prefix=this->prefix());
    CHECK( M_fluidOutletType == "free" || M_fluidOutletType == "windkessel" ) << "invalid fluid-outlet.type " << M_fluidOutletType;
    for ( std::string const& PhysicalName : bcDef.getMarkerNameList<cl::fluid_outlet>() )
        M_fluidOutletsBCType[M_fluidOutletType].push_back(PhysicalName);
}



void
FLUIDMECHANICS_CLASS_NAME::loadConfigMeshFile( std::string const& geofilename )
{
#if !defined(FLUIDMECHANICS_MESH) && !defined(FLUIDMECHANICS_MESH1) && !defined(FLUIDMECHANICS_MESH2) && !defined(FLUIDMECHANICS_MESH3)
    CHECK( false ) << " FLUIDMECHANICS_MESH is not defined : probably .mesh is wrong\n";
#endif

    switch ( this->geotoolMeshIndex() )
    {
    default :
    case 0 :
    {
#if defined(FLUIDMECHANICS_MESH)
        FLUIDMECHANICS_MESH(M_meshSize, geofilename );
        M_mesh=mesh;
#else
        CHECK( false ) << "geometry not define with MeshIndex=" << this->geotoolMeshIndex() << "\n";
#endif
    }
    break;
    case 1 :
    {
#if defined(FLUIDMECHANICS_MESH1)
        FLUIDMECHANICS_MESH1(M_meshSize, geofilename );
        M_mesh=mesh;
#else
        CHECK( false ) << "geometry not define with MeshIndex=" << this->geotoolMeshIndex() << "\n";
#endif
    }
    break;
    case 2 :
    {
#if defined(FLUIDMECHANICS_MESH2)
        FLUIDMECHANICS_MESH2(M_meshSize, geofilename );
        M_mesh=mesh;
#else
        CHECK( false ) << "geometry not define with MeshIndex=" << this->geotoolMeshIndex() << "\n";
#endif
    }
    break;
    case 3 :
    {
#if defined(FLUIDMECHANICS_MESH3)
        FLUIDMECHANICS_MESH3(M_meshSize, geofilename );
        M_mesh=mesh;
#else
        CHECK( false ) << "geometry not define with MeshIndex=" << this->geotoolMeshIndex() << "\n";
#endif
    }
    break;

    } // switch

}


void
FLUIDMECHANICS_CLASS_NAME::init(bool buildMethodNum)
{
    super_type::init( buildMethodNum, this->shared_from_this() );
}



//---------------------------------------------------------------------------------------------------------//

void
FLUIDMECHANICS_CLASS_NAME::updateInitialNewtonSolutionBCDirichlet(vector_ptrtype& U) const
{
    this->log("FluidMechanics","updateInitialNewtonSolutionBCDirichlet", "start");

    if ( !this->hasMarkerDirichletBCelimination() && !this->hasMarkerDirichletBClm() ) return;

    boost::mpi::timer timerBCnewton;

    auto condlim = FLUIDMECHANICS_BC(this->shared_from_this());
    auto const& u = this->fieldVelocity();
    auto Xh = this->functionSpace();
    auto mesh = this->mesh();
    auto const rowStartInVector = this->rowStartInVector();

    if (Xh->worldsComm()[0].isActive()) // only on Velocity Proc
    {
        // modif vector with BC
        ForEachBC( condlim,cl::dirichlet_vec,
                   modifVec(markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",PhysicalName ) ), u, U, Expression, rowStartInVector );
                   modifVec(markedfaces(mesh, this->markerDirichletBCByNameId( "lm",PhysicalName ) ), u, U, Expression, rowStartInVector );
                   );
    }
    //U->close();

    double t1=timerBCnewton.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateInitialNewtonSolutionBCDirichlet",
                                               "finish in "+(boost::format("%1% s") % t1).str(),
                                               this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

void
FLUIDMECHANICS_CLASS_NAME::updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    using namespace Feel::vf;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletLinearPDE", "start",
                                               this->worldComm(),this->verboseAllProc());
    boost::mpi::timer timerBClinear;

    auto const& bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    auto const& u = this->fieldVelocity();

    ForEachBC( bcDef, cl::dirichlet_vec,
               //if (this->worldComm().globalRank()==0) std::cout << "\n integrator on =" << PhysicalName << "\n" << std::endl;
               form2( _test=M_Xh, _trial=M_Xh, _matrix=A,
                      _rowstart=this->rowStartInMatrix(),
                      _colstart=this->colStartInMatrix() ) +=
               on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",PhysicalName ) /*PhysicalName*/),
                   _element=u,
                   _rhs=F,
                   _expr=Expression ) );

               double t1=timerBClinear.elapsed();
                      if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletLinearPDE",
                                                                 "finish in "+(boost::format("%1% s") % t1).str(),
                                                                 this->worldComm(),this->verboseAllProc());
                                                                 }

//---------------------------------------------------------------------------------------------------------//

void
FLUIDMECHANICS_CLASS_NAME::updateBCStrongDirichletJacobian(sparse_matrix_ptrtype& J,vector_ptrtype& RBis) const
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletJacobian", "start",
                                               this->worldComm(),this->verboseAllProc());
    using namespace Feel::vf;

    boost::timer btimeStrongCL;

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();

    auto const rowStartInMatrix = this->rowStartInMatrix();
    auto const colStartInMatrix = this->colStartInMatrix();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );

    auto const& u = this->fieldVelocity();

    //auto RBis = M_backend->newVector( J->mapRowPtr() );

    //boundaries conditions
    auto const& bcDef = FLUIDMECHANICS_BC(this->shared_from_this());

    ForEachBC( bcDef,cl::dirichlet_vec,
               bilinearForm_PatternCoupled +=
               /**/ on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",PhysicalName )/*PhysicalName*/),
                        _element=u,
                        _rhs=RBis,
                        _expr= 0*one()/*Expression-idv(u)*/ ) );

    std::ostringstream ostr3;ostr3<<btimeStrongCL.elapsed()<<"s";
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletJacobian", "finish in "+ostr3.str(),
                                               this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

void
FLUIDMECHANICS_CLASS_NAME::updateBCStrongDirichletResidual(vector_ptrtype& R) const
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletResidual", "start",
                                               this->worldComm(),this->verboseAllProc());
    boost::mpi::timer timerBCresidu;

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto const rowStartInVector = this->rowStartInVector();

    auto const& u = this->fieldVelocity();

    //boundaries conditions
    auto const& bcDef = FLUIDMECHANICS_BC(this->shared_from_this());

    R->close();
    if (Xh->worldsComm()[0].isActive()) // only on Velocity Proc
    {
        // Zero because is know
        ForEachBC( bcDef,cl::dirichlet_vec,
                   modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",PhysicalName )/*PhysicalName*/), u, R, 0.*vf::one(),rowStartInVector ) );
        //ForEachBC( bcDef,cl::fbm_dirichlet,
        //           modifVec(markedfaces(mesh,PhysicalName), u, R, vf::one()-vf::one(),rowStartInVector ) );
    }

    double t1=timerBCresidu.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletResidual",
                                               "finish in "+(boost::format("%1% s") % t1).str(),
                                               this->worldComm(),this->verboseAllProc());
}

void
FLUIDMECHANICS_CLASS_NAME::updateSourceTermLinearPDE( vector_ptrtype& F, bool BuildCstPart ) const
{
    if ( M_overwritemethod_updateSourceTermLinearPDE != NULL )
    {
        M_overwritemethod_updateSourceTermLinearPDE(F,BuildCstPart);
        return;
    }

#if defined(FLUIDMECHANICS_VOLUME_FORCE)
    bool BuildSourceTerm = !BuildCstPart;
    if (this->useFSISemiImplicitScheme())
        BuildSourceTerm=BuildCstPart;

    if ( BuildSourceTerm )
    {
        auto myLinearForm =form1( _test=this->functionSpace(), _vector=F,
                                  _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldVelocity();
        auto f = FLUIDMECHANICS_VOLUME_FORCE(this->shared_from_this()) ;
        myLinearForm +=
            integrate( _range=elements(this->mesh()),
                       _expr= trans(f)*id(v),
#if defined(FLUIDMECHANICS_VOLUME_FORCE_QUADORDER)
                       _quad=_Q<FLUIDMECHANICS_VOLUME_FORCE_QUADORDER>(),
                       _quad1=_Q<FLUIDMECHANICS_VOLUME_FORCE_QUADORDER>(),
#endif
                       _geomap=this->geomap() );
    }
#endif
}
void
FLUIDMECHANICS_CLASS_NAME::updateSourceTermResidual( vector_ptrtype& F ) const
{
    if ( M_overwritemethod_updateSourceTermResidual != NULL )
    {
        M_overwritemethod_updateSourceTermResidual(F);
        return;
    }

#if defined(FLUIDMECHANICS_VOLUME_FORCE)
    auto myLinearForm =form1( _test=this->functionSpace(), _vector=F,
                              _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldVelocity();
    auto f = FLUIDMECHANICS_VOLUME_FORCE(this->shared_from_this()) ;
    myLinearForm +=
        integrate( _range=elements(this->mesh()),
                   _expr= -trans(f)*id(v),
#if defined(FLUIDMECHANICS_VOLUME_FORCE_QUADORDER)
                   _quad=_Q<FLUIDMECHANICS_VOLUME_FORCE_QUADORDER>(),
                   _quad1=_Q<FLUIDMECHANICS_VOLUME_FORCE_QUADORDER>(),
#endif
                   _geomap=this->geomap() );
#endif
}

void
FLUIDMECHANICS_CLASS_NAME::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
    auto const bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    if ( bcDef.hasNeumannScal() )
    {
        auto myLinearForm = form1( _test=this->functionSpace(), _vector=F,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldVelocity();
        ForEachBC( bcDef, cl::neumann_scal,
                   myLinearForm +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,PhysicalName) ),
                                   _expr= Expression*inner( N(),id(v) ),
                                   _geomap=this->geomap() ) );
    }
}

void
FLUIDMECHANICS_CLASS_NAME::updateBCNeumannResidual( vector_ptrtype& R ) const
{
    auto const bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    if ( bcDef.hasNeumannScal() )
    {
        auto myLinearForm = form1( _test=this->functionSpace(), _vector=R,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldVelocity();
        ForEachBC( bcDef, cl::neumann_scal,
                   myLinearForm +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,PhysicalName) ),
                                   _expr= -Expression*inner( N(),id(v) ),
                                   _geomap=this->geomap() ) );
    }
}

void
FLUIDMECHANICS_CLASS_NAME::updateBCPressureLinearPDE( vector_ptrtype& F ) const
{
    auto const bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    if ( bcDef.hasPressure() )
    {
        auto myLinearForm = form1( _test=this->functionSpace(), _vector=F,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldVelocity();
        ForEachBC( bcDef,cl::pressure,
                   myLinearForm +=
                   /**/ integrate( _range=markedfaces(this->mesh(),PhysicalName),
                                   _expr= trans(-Expression*N())*id(v),
                                   _geomap=this->geomap() ) );
    }
}

void
FLUIDMECHANICS_CLASS_NAME::updateBCPressureResidual( vector_ptrtype& R ) const
{
    auto const bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    if ( bcDef.hasPressure() )
    {
        auto myLinearForm = form1( _test=this->functionSpace(), _vector=R,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldVelocity();
        ForEachBC( bcDef,cl::pressure,
                   myLinearForm +=
                   /**/ integrate( _range=markedfaces(this->mesh(),PhysicalName),
                                   _expr= -trans(-Expression*N())*id(v),
                                   _geomap=this->geomap() ) );
    }
}

void
FLUIDMECHANICS_CLASS_NAME::updateBCDirichletLagMultLinearPDE( vector_ptrtype& F ) const
{
    auto const bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    if ( bcDef.hasDirichletVec() )
    {
        auto lambdaBC = this->XhDirichletLM()->element();
        size_type startDofIndexDirichletLM = this->startDofIndexFieldsInMatrix().find("dirichletlm")->second;
        ForEachBC( bcDef, cl::dirichlet_vec,
                   form1( _test=this->XhDirichletLM(),_vector=F,
                          _rowstart=this->rowStartInVector()+startDofIndexDirichletLM ) +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "lm",PhysicalName ) ),
                                   //_range=markedelements(this->meshDirichletLM(),PhysicalName),
                                   _expr= inner( Expression,id(lambdaBC) ),
                                   _geomap=this->geomap() ) );
    }
}

void
FLUIDMECHANICS_CLASS_NAME::updateBCDirichletLagMultResidual( vector_ptrtype& R ) const
{
    auto const bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    if ( bcDef.hasDirichletVec() )
    {
        auto lambdaBC = this->XhDirichletLM()->element();
        size_type startDofIndexDirichletLM = this->startDofIndexFieldsInMatrix().find("dirichletlm")->second;
        ForEachBC( bcDef, cl::dirichlet_vec,
                   form1( _test=this->XhDirichletLM(),_vector=R,
                          _rowstart=this->rowStartInVector()+startDofIndexDirichletLM ) +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "lm",PhysicalName ) ),
                                   //_range=markedelements(this->meshDirichletLM(),PhysicalName),
                                   _expr= -inner( Expression,id(lambdaBC) ),
                                   _geomap=this->geomap() ) );
    }
}

void
FLUIDMECHANICS_CLASS_NAME::updateBCDirichletNitscheLinearPDE( vector_ptrtype& F ) const
{
    auto const bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    if ( bcDef.hasDirichletVec() )
    {
        auto myLinearForm = form1( _test=this->functionSpace(), _vector=F,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldVelocity();
        ForEachBC( bcDef, cl::dirichlet_vec,
                   myLinearForm +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "nitsche",PhysicalName ) ),
                                   _expr= this->dirichletBCnitscheGamma()*trans(Expression)*id(v)/hFace(),
                                   _geomap=this->geomap() ) );
    }
}

void
FLUIDMECHANICS_CLASS_NAME::updateBCDirichletNitscheResidual( vector_ptrtype& R ) const
{
    auto const bcDef = FLUIDMECHANICS_BC(this->shared_from_this());
    if ( bcDef.hasDirichletVec() )
    {
        auto myLinearForm = form1( _test=this->functionSpace(), _vector=R,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldVelocity();
        ForEachBC( bcDef, cl::dirichlet_vec,
                   myLinearForm +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "nitsche",PhysicalName ) ),
                                   _expr= - this->dirichletBCnitscheGamma()*inner( Expression, id(v) )/hFace(),
                                   _geomap=this->geomap() ) );
    }
}


} // namespace FeelModels
} // namespace Feel
