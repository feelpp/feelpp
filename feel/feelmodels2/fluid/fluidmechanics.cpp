/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
 \file fluidmechanics.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-17
 */

#include <feel/feelmodels2/fluid/fluidmechanics.hpp>

#include <feel/feelmodels2/modelalg/functionSup.cpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::FluidMechanics( //bool __isStationary,
                                                    std::string __prefix,
                                                    bool __buildMesh,
                                                    WorldComm const& __worldComm,
                                                    std::string __subPrefix,
                                                    std::string __appliShortRepository )
    :
    super_type( __prefix, __buildMesh,__worldComm, __subPrefix,__appliShortRepository)
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    // clear
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerALEMeshBC();
    this->clearMarkerSlipBC();
    this->clearMarkerPressureBC();
    this->M_fluidOutletsBCType.clear();

    // boundary conditions
    this->M_isMoveDomain = false;
    std::string dirichletbcType = soption(_name="dirichletbc.type",_prefix=this->prefix());
    CHECK( dirichletbcType=="elimination" || dirichletbcType=="nitsche" || dirichletbcType=="lm" ) << "invalid dirichletbc.type " << dirichletbcType;
    M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "velocity", "Dirichlet" );
    for( auto const& d : M_bcDirichlet )
    {
        this->addMarkerDirichletBC( dirichletbcType, marker(d), ComponentType::NO_COMPONENT );
        this->addMarkerALEMeshBC("fixed",marker(d));
    }
    M_bcMovingBoundary = this->modelProperties().boundaryConditions().getScalarFields( "velocity", "interface_fsi"/*"moving_boundary"*/ );
    for( auto const& d : M_bcMovingBoundary )
    {
        //this->addMarkerDirichletBC( dirichletbcType, marker(d), ComponentType::NO_COMPONENT ); // ??
        this->addMarkerALEMeshBC("moving",marker(d));
        this->M_isMoveDomain=true;
    }
    M_bcNeumannScalar = this->modelProperties().boundaryConditions().getScalarFields( "velocity", "Neumann_scalar" );
    for( auto const& d : M_bcNeumannScalar )
    {
        this->addMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d));
        this->addMarkerALEMeshBC("fixed",marker(d));
    }
    M_bcPressure = this->modelProperties().boundaryConditions().getScalarFields( "pressure", "weak" );
    for( auto const& d : M_bcPressure )
    {
        this->addMarkerPressureBC(marker(d));
        this->addMarkerALEMeshBC("fixed",marker(d));
    }
    M_bcSlip = this->modelProperties().boundaryConditions().getScalarFields( "velocity", "slip" );
    for( auto const& d : M_bcSlip )
    {
        this->addMarkerSlipBC(marker(d));
        this->addMarkerALEMeshBC("fixed",marker(d));
    }

    this->M_hasFluidOutlet=false;
    this->M_fluidOutletType = soption(_name="fluid-outlet.type", _prefix=this->prefix());
    CHECK( this->M_fluidOutletType == "free" || this->M_fluidOutletType == "windkessel" ) << "invalid fluid-outlet.type " << this->M_fluidOutletType;
    M_bcFluidOutlets = this->modelProperties().boundaryConditions().getScalarFields( "fluid", "outlet" );
    for( auto const& d : M_bcFluidOutlets )
    {
        this->M_fluidOutletsBCType[this->M_fluidOutletType].push_back(marker(d));
        this->M_hasFluidOutlet=true;
    }

    M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "fluid", "VolumicForces" );

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::loadConfigMeshFile( std::string const& geofilename )
{
    CHECK( false ) << "not allow";
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::solve()
{
    std::map<std::string,double> mySymbolsValues = { {"t",this->currentTime()} };
    M_bcDirichlet.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcDirichlet.setParameterValues( mySymbolsValues );
    M_bcNeumannScalar.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcNeumannScalar.setParameterValues( mySymbolsValues );
    M_volumicForcesProperties.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_volumicForcesProperties.setParameterValues( mySymbolsValues );
    super_type::solve();
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateCLDirichlet(vector_ptrtype& U) const
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateCLDirichlet", "start",
                                               this->worldComm(),this->verboseAllProc());

    if ( !this->hasMarkerDirichletBCelimination() && !this->hasMarkerDirichletBClm() ) return;

    boost::mpi::timer timerBCnewton;

    auto const& u = this->fieldVelocity();
    auto Xh = this->functionSpace();
    auto mesh = this->mesh();
    size_type rowStartInVector = this->rowStartInVector();

    if (Xh->worldsComm()[0].isActive()) // only on Velocity Proc
    {
        // modif vector with BC
        for( auto const& d : M_bcDirichlet )
        {
            modifVec(markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",marker(d) ) ), u, U, expression(d), rowStartInVector );
            modifVec(markedfaces(mesh, this->markerDirichletBCByNameId( "lm",marker(d) ) ), u, U, expression(d), rowStartInVector );
        }
#if defined( FEELPP_MODELS_HAS_MESHALE ) // must be move in base class
        if (this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet")
        {
#if 0
            std::list<std::string> movingBCmarkers = detail::intersectionList( this->markersNameMovingBoundary(),
                                                                               this->markerDirichletBCelimination() );
            // modif vector with BC
            for ( std::string const& marker : movingBCmarkers )
                modifVec(markedfaces(mesh, marker), u, U, vf::idv(this->meshVelocity2()), rowStartInVector );
#else
            modifVec(markedfaces(mesh, this->markersNameMovingBoundary()), u, U, vf::idv(this->meshVelocity2()), rowStartInVector );
#endif
        }
#endif
    }
    U->close();

    double t1=timerBCnewton.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateCLDirichlet",
                                               "finish in "+(boost::format("%1% s") % t1).str(),
                                               this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    using namespace Feel::vf;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletLinearPDE", "start",
                                               this->worldComm(),this->verboseAllProc());
    boost::mpi::timer timerBClinear;

    auto Xh = this->functionSpace();
    auto const& u = this->fieldVelocity();

#if defined( FEELPP_MODELS_HAS_MESHALE ) // must be move in base class
    if (this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet")
    {
        form2( _test=Xh, _trial=Xh, _matrix=A,
               _rowstart=this->rowStartInMatrix(),
               _colstart=this->colStartInMatrix() ) +=
            on( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                _element=u,
                _rhs=F,
                _expr=idv(this->meshVelocity2()) );
    }
#endif

    for( auto const& d : M_bcDirichlet )
        form2( _test=Xh, _trial=Xh, _matrix=A,
               _rowstart=this->rowStartInMatrix(),
               _colstart=this->colStartInMatrix() ) +=
            on( _range=markedfaces(this->mesh(), this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                _element=u,
                _rhs=F,
                _expr=expression(d) );

    double t1=timerBClinear.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletLinearPDE",
                                               "finish in "+(boost::format("%1% s") % t1).str(),
                                               this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletJacobian(sparse_matrix_ptrtype& J) const
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletJacobian", "start",
                                               this->worldComm(),this->verboseAllProc());
    using namespace Feel::vf;

    boost::timer btimeStrongCL;

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldVelocity();
    auto RBis = this->backend()->newVector( J->mapRowPtr() );

#if defined( FEELPP_MODELS_HAS_MESHALE ) // must be move in base class
    if (this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet")
    {
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markersNameMovingBoundary()),
                _element=u,
                _rhs=RBis,
                _expr= 0*one()/*idv(this->meshVelocity2()) - idv(u)*/ );
    }
#endif
    for( auto const& d : M_bcDirichlet )
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                _element=u,
                _rhs=RBis,
                _expr= 0*one()/*Expression-idv(u)*/ );

    std::ostringstream ostr3;ostr3<<btimeStrongCL.elapsed()<<"s";
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletJacobian", "finish in "+ostr3.str(),
                                               this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletResidual(vector_ptrtype& R) const
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletResidual", "start",
                                               this->worldComm(),this->verboseAllProc());
    boost::mpi::timer timerBCresidu;

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    size_type rowStartInVector = this->rowStartInVector();
    auto const& u = this->fieldVelocity();

    R->close();
    if (Xh->worldsComm()[0].isActive()) // only on Velocity Proc
    {
        // Zero because is know
        for( auto const& d : M_bcDirichlet )
            modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",marker(d) )),
                     u, R, 0.*vf::one(),rowStartInVector );
#if defined( FEELPP_MODELS_HAS_MESHALE ) // must be move in base class
        if (this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet")
        {
            modifVec(markedfaces(mesh,this->markersNameMovingBoundary()), u, R, 0*vf::one(),rowStartInVector );
        }
#endif
    }

    double t1=timerBCresidu.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletResidual",
                                               "finish in "+(boost::format("%1% s") % t1).str(),
                                               this->worldComm(),this->verboseAllProc());
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateSourceTermLinearPDE( vector_ptrtype& F, bool BuildCstPart ) const
{
    if ( this->M_overwritemethod_updateSourceTermLinearPDE != NULL )
    {
        this->M_overwritemethod_updateSourceTermLinearPDE(F,BuildCstPart);
        return;
    }

    if ( M_volumicForcesProperties.empty() ) return;

    bool BuildSourceTerm = !BuildCstPart;
    if (this->useFSISemiImplicitScheme())
        BuildSourceTerm=BuildCstPart;

    if ( BuildSourceTerm )
    {
        auto myLinearForm =form1( _test=this->functionSpace(), _vector=F,
                                  _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldVelocity();
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
}
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateSourceTermResidual( vector_ptrtype& F ) const
{
    if ( this->M_overwritemethod_updateSourceTermResidual != NULL )
    {
        this->M_overwritemethod_updateSourceTermResidual(F);
        return;
    }

    if ( M_volumicForcesProperties.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpace(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldVelocity();
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
    if ( M_bcNeumannScalar.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpace(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldVelocity();
    for( auto const& d : M_bcNeumannScalar )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d)) ),
                       _expr= expression(d)*inner( N(),id(v) ),
                       _geomap=this->geomap() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannResidual( vector_ptrtype& R ) const
{
    if ( M_bcNeumannScalar.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpace(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldVelocity();
    for( auto const& d : M_bcNeumannScalar )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d)) ),
                       _expr= -expression(d)*inner( N(),id(v) ),
                       _geomap=this->geomap() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCPressureLinearPDE( vector_ptrtype& F ) const
{
#if 0 //TODO
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
#endif // TODO
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCPressureResidual( vector_ptrtype& R ) const
{
#if 0 //TODO
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
#endif // TODO
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletLagMultLinearPDE( vector_ptrtype& F ) const
{
    if ( M_bcDirichlet.empty() ) return;

    auto lambdaBC = this->XhDirichletLM()->element();
    size_type startDofIndexDirichletLM = this->startDofIndexFieldsInMatrix().find("dirichletlm")->second;
    for( auto const& d : M_bcDirichlet )
        form1( _test=this->XhDirichletLM(),_vector=F,
               _rowstart=this->rowStartInVector()+startDofIndexDirichletLM ) +=
            integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "lm",marker(d) ) ),
                       //_range=markedelements(this->meshDirichletLM(),PhysicalName),
                       _expr= inner( expression(d),id(lambdaBC) ),
                       _geomap=this->geomap() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletLagMultResidual( vector_ptrtype& R ) const
{
    if ( M_bcDirichlet.empty() ) return;

    auto lambdaBC = this->XhDirichletLM()->element();
    size_type startDofIndexDirichletLM = this->startDofIndexFieldsInMatrix().find("dirichletlm")->second;
    for( auto const& d : M_bcDirichlet )
        form1( _test=this->XhDirichletLM(),_vector=R,
               _rowstart=this->rowStartInVector()+startDofIndexDirichletLM ) +=
            integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "lm",marker(d) ) ),
                       //_range=markedelements(this->meshDirichletLM(),PhysicalName),
                       _expr= -inner( expression(d),id(lambdaBC) ),
                       _geomap=this->geomap() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletNitscheLinearPDE( vector_ptrtype& F ) const
{
    if ( M_bcDirichlet.empty() ) return;
    auto myLinearForm = form1( _test=this->functionSpace(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldVelocity();
    for( auto const& d : M_bcDirichlet )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "nitsche",marker(d) ) ),
                       _expr= this->dirichletBCnitscheGamma()*inner( expression(d),id(v) )/hFace(),
                       _geomap=this->geomap() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletNitscheResidual( vector_ptrtype& R ) const
{
    if ( M_bcDirichlet.empty() ) return;
    auto myLinearForm = form1( _test=this->functionSpace(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldVelocity();
    for( auto const& d : M_bcDirichlet )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerDirichletBCByNameId( "nitsche",marker(d) ) ),
                       _expr= -this->dirichletBCnitscheGamma()*inner( expression(d),id(v) )/hFace(),
                       _geomap=this->geomap() );
}


} // namespace FeelModels

} // namespace Feel
