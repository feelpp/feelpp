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
 \file fluidmechanics.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-17
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelmodels/modelalg/functionSup.cpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelpde/preconditionerblockns.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::FluidMechanics( std::string __prefix,
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
    if (__buildMesh) this->build();
    //-----------------------------------------------------------------------------//

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

namespace detailbc
{
std::list<std::string>
generateMarkerBCList(BoundaryConditions const bc,std::string const& field,std::string const& bcname, std::string const& marker, std::string const& param="number" )
{
    std::list<std::string> markerList;

    std::pair<bool,int> numberOfMarkerRead = bc.iparam( field/*"velocity"*/, bcname/*"Dirichlet"*/, marker/*(d)*/, param/*"number"*/ );
    int numberOfMarker = ( numberOfMarkerRead.first )? numberOfMarkerRead.second : 1;
    for (int k=0 ; k<numberOfMarker ; ++k )
    {
        std::string currentMarker = ( numberOfMarker == 1 )? marker : (boost::format("%1%%2%")%marker %k).str();
        markerList.push_back( currentMarker );
    }

    return markerList;
}
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

    // change path
    fs::path curPath=fs::current_path();
    bool hasChangedRep=false;
    if ( curPath != fs::path(this->ginacExprCompilationDirectory()) )
    {
        this->log("FluidMechanics","loadConfigBCFile", "change repository (temporary) for build ginac expr with default name : "+ this->appliRepository() );
        bool hasChangedRep=true;
        Environment::changeRepository( _directory=boost::format(this->ginacExprCompilationDirectory()), _subdir=false );
    }

    // boundary conditions
    this->M_isMoveDomain = false;
    M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "velocity", "Dirichlet" );
    for( auto const& d : M_bcDirichlet )
    {
        std::pair<bool,std::string> dirichletbcTypeRead = this->modelProperties().boundaryConditions().sparam( "velocity", "Dirichlet", marker(d), "type" );
        std::string dirichletbcType = ( dirichletbcTypeRead.first )? dirichletbcTypeRead.second : soption(_name="dirichletbc.type",_prefix=this->prefix());
        CHECK( dirichletbcType=="elimination" || dirichletbcType=="nitsche" || dirichletbcType=="lm" ) << "invalid dirichletbc.type " << dirichletbcType;

        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "velocity", "Dirichlet", marker(d) );
        this->setMarkerDirichletBCByNameId( dirichletbcType, marker(d), markerList,ComponentType::NO_COMPONENT );

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Dirichlet", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
    }
    M_bcMovingBoundary = this->modelProperties().boundaryConditions().getScalarFields( { { "velocity", "interface_fsi" }, { "velocity","moving_boundary"} } );
    for( auto const& d : M_bcMovingBoundary )
    {
        this->addMarkerALEMeshBC("moving",marker(d));
        this->M_isMoveDomain=true;
    }
    M_bcNeumannScalar = this->modelProperties().boundaryConditions().getScalarFields( "velocity", "Neumann_scalar" );
    for( auto const& d : M_bcNeumannScalar )
    {
        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "velocity", "Neumann_scalar", marker(d) );
        this->setMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d),markerList);

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Neumann_scalar", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
    }
    M_bcPressure = this->modelProperties().boundaryConditions().getScalarFields( "pressure", "weak" );
    for( auto const& d : M_bcPressure )
    {
        this->addMarkerPressureBC(marker(d));
        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "pressure", "weak", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        this->addMarkerALEMeshBC(bcTypeMeshALE,marker(d));
    }
    M_bcSlip = this->modelProperties().boundaryConditions().getScalarFields( "velocity", "slip" );
    for( auto const& d : M_bcSlip )
    {
        this->addMarkerSlipBC(marker(d));
        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "slip", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        this->addMarkerALEMeshBC(bcTypeMeshALE,marker(d));
    }
    M_bcFluidOutlets = this->modelProperties().boundaryConditions().getScalarFields( "fluid", "outlet" );
    for( auto const& d : M_bcFluidOutlets )
    {
        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "fluid", "outlet", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");

        std::string typeOutlet = soption(_name="fluid-outlet.type", _prefix=this->prefix());//"free";
        std::pair<bool,std::string> typeOutletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "outlet", marker(d), "type" );
        if ( typeOutletRead.first )
        {
            typeOutlet = typeOutletRead.second;
            CHECK( typeOutlet == "free" || typeOutlet == "windkessel" ) << "invalid outlet type " << typeOutlet;
        }
        std::string typeCouplingWindkesselOutlet = soption(_name="fluid-outlet.windkessel.coupling", _prefix=this->prefix());
        std::pair<bool,std::string> typeCouplingWindkesselOutletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "outlet", marker(d), "windkessel_coupling" );
        if ( typeCouplingWindkesselOutletRead.first )
        {
            typeCouplingWindkesselOutlet = typeCouplingWindkesselOutletRead.second;
            CHECK( typeCouplingWindkesselOutlet == "implicit" || typeCouplingWindkesselOutlet == "explicit" ) << "invalid windkessel coupling type " << typeCouplingWindkesselOutlet;
        }
        std::pair<bool,double> WindkesselRdRead = this->modelProperties().boundaryConditions().dparam( "fluid", "outlet", marker(d), "windkessel_Rd" );
        std::pair<bool,double> WindkesselRpRead = this->modelProperties().boundaryConditions().dparam( "fluid", "outlet", marker(d), "windkessel_Rp" );
        std::pair<bool,double> WindkesselCdRead = this->modelProperties().boundaryConditions().dparam( "fluid", "outlet", marker(d), "windkessel_Cd" );
        double WindkesselRd = ( WindkesselRdRead.first )? WindkesselRdRead.second : 1.;
        double WindkesselRp = ( WindkesselRpRead.first )? WindkesselRpRead.second : 1.;
        double WindkesselCd = ( WindkesselCdRead.first )? WindkesselCdRead.second : 1.;

        std::tuple<std::string,double,double,double> windkesselParam = std::make_tuple(typeCouplingWindkesselOutlet,WindkesselRd,WindkesselRp,WindkesselCd);

        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "fluid", "outlet", marker(d) );
        for (std::string const& currentMarker : markerList )
        {
            this->M_fluidOutletsBCType.push_back(std::make_tuple(currentMarker,typeOutlet, windkesselParam ));
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
        }
    }

    M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "fluid", "VolumicForces" );

    // go back to previous repository
    if ( hasChangedRep )
        Environment::changeRepository( _directory=boost::format(curPath.string()), _subdir=false );

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::loadConfigPostProcess()
{
    if ( this->modelProperties().postProcess().find("Fields") != this->modelProperties().postProcess().end() )
        for ( auto const& o : this->modelProperties().postProcess().find("Fields")->second )
        {
            if ( o == "velocity" || o == "all" ) this->M_doExportVelocity = true;
            if ( o == "pressure" || o == "all" ) this->M_doExportPressure = true;
            if ( o == "displacement" || o == "all" ) this->M_doExportMeshDisplacement = true;
            if ( o == "vorticity" || o == "all" ) this->M_doExportVorticity = true;
            if ( o == "stress" || o == "normal-stress" || o == "all" ) this->M_doExportNormalStress = true;
            if ( o == "wall-shear-stress" || o == "all" ) this->M_doExportWallShearStress = true;
            if ( o == "viscosity" || o == "all" ) this->M_doExportViscosity = true;
        }
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

    if ( buildModelAlgebraicFactory && ( soption(_prefix=this->prefix(),_name="pc-type" ) == "blockns" ) )
    {
        BoundaryConditions bcPrecPCD;
        bcPrecPCD.clear();

        auto itFindFieldVelocity = this->modelProperties().boundaryConditions().find("velocity");
        bool hasFindFieldVelocity = itFindFieldVelocity != this->modelProperties().boundaryConditions().end();
        if ( hasFindFieldVelocity )
        {
            auto itFindDirichletType = itFindFieldVelocity->second.find("Dirichlet");
            if ( itFindDirichletType != itFindFieldVelocity->second.end() )
            {
                for ( auto const& myBcDesc : itFindDirichletType->second )
                    bcPrecPCD["velocity"]["Dirichlet"].push_back( myBcDesc );
            }
            auto itFindNeumannScalType = itFindFieldVelocity->second.find("Neumann_scalar");
            if ( itFindNeumannScalType != itFindFieldVelocity->second.end() )
            {
                for ( auto const& myBcDesc : itFindNeumannScalType->second )
                    bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc );
            }
        }
        auto itFindFieldFluid = this->modelProperties().boundaryConditions().find("fluid");
        if ( itFindFieldFluid != this->modelProperties().boundaryConditions().end() )
        {
            auto itFindOutletType = itFindFieldFluid->second.find("outlet");
            if ( itFindOutletType != itFindFieldFluid->second.end() )
            {
                for ( auto const& myBcDesc : itFindOutletType->second )
                    bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc );
            }
        }
        // TODO other bc (fsi,...)
#if 1
        if ( Environment::isMasterRank() )
        {
            for( auto const& s : bcPrecPCD )
            {
                std::cout << "field " << s.first << "\n";
                for( auto const& t : s.second )
                {
                    std::cout << " - type " << t.first << "\n";
                    for( auto const& c : t.second )
                    {
                        if ( c.hasExpression2() )
                            std::cout << "  . boundary  " << c.marker() << " expr : " << c.expression1() << " expr2:" << c.expression2() << "\n";
                        else
                            std::cout << "  . boundary  " << c.marker() << " expr : " << c.expression() << "\n";
                    }
                }
            }
        }
#endif
        CHECK( this->algebraicFactory()->preconditionerTool()->matrix() ) << "no matrix define in preconditionerTool";
        double myalpha = (this->isStationary())? 0 : this->densityViscosityModel()->cstRho()*this->timeStepBDF()->polyDerivCoefficient(0);
        auto a_blockns = Feel::blockns( _space=this->functionSpace(),
                                        _type=soption(_prefix=this->prefix(),_name="blockns.type"),//"PCD",
                                        _bc=bcPrecPCD,
                                        _matrix=this->algebraicFactory()->preconditionerTool()->matrix(),
                                        _prefix="velocity",
                                        _mu=this->densityViscosityModel()->cstMu(),
                                        _rho=this->densityViscosityModel()->cstRho(),
                                        _alpha=myalpha );
        this->algebraicFactory()->preconditionerTool()->attachInHousePreconditioners("blockns",a_blockns);
    }

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateInHousePreconditionerPCD( sparse_matrix_ptrtype const& mat,vector_ptrtype const& vecSol) const
{
    if ( this->algebraicFactory()->preconditionerTool()->hasInHousePreconditioners( "blockns" ) )
    {
        this->log("FluidMechanics","updateInHousePreconditionerPCD", "start");

        boost::shared_ptr< PreconditionerBlockNS<typename super_type::space_fluid_type> > myPrecBlockNs =
            boost::dynamic_pointer_cast< PreconditionerBlockNS<typename super_type::space_fluid_type> >( this->algebraicFactory()->preconditionerTool()->inHousePreconditioners( "blockns" ) );

        if ( this->isStationary() )
            myPrecBlockNs->setAlpha( 0. );
        else
            myPrecBlockNs->setAlpha( this->densityViscosityModel()->cstRho()*this->timeStepBDF()->polyDerivCoefficient(0) );

        if ( this->pdeType() == "Stokes" )
        {
            myPrecBlockNs->update( mat );
        }
        else if ( this->pdeType() == "Oseen" )
        {
            auto BetaU = this->timeStepBDF()->poly();
            auto betaU = BetaU.template element<0>();
            auto const& rho = this->densityViscosityModel()->fieldRho();

            if (this->isMoveDomain() )
            {
#if defined( FEELPP_MODELS_HAS_MESHALE )
                myPrecBlockNs->update( mat, idv(rho)*( idv(betaU)-idv(this->meshVelocity()) ) );
#endif
            }
            else
            {
                myPrecBlockNs->update( mat, idv(rho)*idv(betaU) );
            }
        }
        else if ( this->pdeType() == "Navier-Stokes" )
        {
            auto U = this->functionSpace()->element();
            // copy vector values in fluid element
            for ( size_type k=0;k<this->functionSpace()->nLocalDofWithGhost();++k )
                U(k) = vecSol->operator()(/*rowStartInVector+*/k);
            auto u = U.template element<0>();
            auto const& rho = this->densityViscosityModel()->fieldRho();

            if (this->isMoveDomain() )
            {
#if defined( FEELPP_MODELS_HAS_MESHALE )
                myPrecBlockNs->update( mat, idv(rho)*( idv(u)-idv(this->meshVelocity()) ) );
#endif
            }
            else
            {
                myPrecBlockNs->update( mat, idv(rho)*idv(u) );
            }
        }

        this->log("FluidMechanics","updateInHousePreconditionerPCD", "finish");
    }
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::solve()
{
    M_bcDirichlet.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_bcNeumannScalar.setParameterValues( this->modelProperties().parameters().toParameterValues() );
    M_volumicForcesProperties.setParameterValues( this->modelProperties().parameters().toParameterValues() );

    if ( this->algebraicFactory() && this->algebraicFactory()->preconditionerTool()->hasInHousePreconditioners( "blockns" ) )
    {
        boost::shared_ptr< PreconditionerBlockNS<typename super_type::space_fluid_type> > myPrecBlockNs =
            boost::dynamic_pointer_cast< PreconditionerBlockNS<typename super_type::space_fluid_type> >( this->algebraicFactory()->preconditionerTool()->inHousePreconditioners( "blockns" ) );
        myPrecBlockNs->setParameterValues( this->modelProperties().parameters().toParameterValues() );
    }

    super_type::solve();
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateInitialNewtonSolutionBCDirichlet(vector_ptrtype& U) const
{
    if ( !this->hasMarkerDirichletBCelimination() && !this->hasMarkerDirichletBClm() ) return;

    this->log("FluidMechanics","updateInitialNewtonSolutionBCDirichlet", "start");

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
    }
    //U->close();

    double t1=timerBCnewton.elapsed();
    this->log("FluidMechanics","updateInitialNewtonSolutionBCDirichlet","finish in "+(boost::format("%1% s") % t1).str() );
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
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletJacobian(sparse_matrix_ptrtype& J,vector_ptrtype& RBis) const
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
    //auto RBis = this->backend()->newVector( J->mapRowPtr() );

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

    //R->close();
    if (Xh->worldsComm()[0].isActive()) // only on Velocity Proc
    {
        // Zero because is know
        for( auto const& d : M_bcDirichlet )
            modifVec(markedfaces(mesh,this->markerDirichletBCByNameId( "elimination",marker(d) )),
                     u, R, 0.*vf::one(),rowStartInVector );
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
