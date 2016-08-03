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
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::FluidMechanics( std::string const& prefix,
                                                    bool buildMesh,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    std::string const& rootRepository )
    :
    super_type( prefix, buildMesh, worldComm, subPrefix, rootRepository )
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
    if ( buildMesh ) this->build();
    //-----------------------------------------------------------------------------//

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::self_ptrtype
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::New( std::string const& prefix, bool buildMesh,
                                         WorldComm const& worldComm, std::string const& subPrefix,
                                         std::string const& rootRepository )
{
    return boost::make_shared<self_type>( prefix, buildMesh, worldComm, subPrefix, rootRepository );

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

    // boundary conditions
    this->M_isMoveDomain = false;
    M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "velocity", "Dirichlet" );
    for( auto const& d : M_bcDirichlet )
    {
        std::pair<bool,std::string> dirichletbcTypeRead = this->modelProperties().boundaryConditions().sparam( "velocity", "Dirichlet", marker(d), "method" );
        std::string dirichletbcType = ( dirichletbcTypeRead.first )? dirichletbcTypeRead.second : soption(_name="dirichletbc.type",_prefix=this->prefix());
        CHECK( dirichletbcType=="elimination" || dirichletbcType=="nitsche" || dirichletbcType=="lm" ) << "invalid dirichletbc.type " << dirichletbcType;

        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "velocity", "Dirichlet", marker(d) );
        this->setMarkerDirichletBCByNameId( dirichletbcType, marker(d), markerList,ComponentType::NO_COMPONENT );

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Dirichlet", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
    }
    for ( ComponentType comp : std::vector<ComponentType>( { ComponentType::X, ComponentType::Y, ComponentType::Z } ) )
    {
        std::string compTag = ( comp ==ComponentType::X )? "x" : (comp == ComponentType::Y )? "y" : "z";
        std::string bcDirichletCompField = (boost::format("velocity_%1%")%compTag).str();
        std::string bcDirichletCompKeyword = "Dirichlet";
        M_bcDirichletComponents[comp] = this->modelProperties().boundaryConditions().getScalarFields( { { bcDirichletCompField, bcDirichletCompKeyword } } );
        for( auto const& d : M_bcDirichletComponents.find(comp)->second )
        {
            std::pair<bool,std::string> dirichletbcTypeRead = this->modelProperties().boundaryConditions().sparam( bcDirichletCompField, bcDirichletCompKeyword, marker(d), "method" );
            std::string dirichletbcType = ( dirichletbcTypeRead.first )? dirichletbcTypeRead.second : soption(_name="dirichletbc.type",_prefix=this->prefix());
            CHECK( dirichletbcType=="elimination" || dirichletbcType=="nitsche" || dirichletbcType=="lm" ) << "invalid dirichletbc.type " << dirichletbcType;

            std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), bcDirichletCompField, bcDirichletCompKeyword, marker(d) );
            this->setMarkerDirichletBCByNameId( dirichletbcType, marker(d), markerList, comp );

            std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( bcDirichletCompField, bcDirichletCompKeyword, marker(d), "alemesh_bc" );
            std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
            for (std::string const& currentMarker : markerList )
                this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
        }
    }

    for( std::string const& bcMarker : this->modelProperties().boundaryConditions().markers( { { "velocity", "interface_fsi" }, { "velocity","moving_boundary"} } ) )
    {
        this->addMarkerALEMeshBC("moving",bcMarker);
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
    M_bcNeumannVectorial = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "velocity", "Neumann_vectorial" );
    for( auto const& d : M_bcNeumannVectorial )
    {
        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "velocity", "Neumann_vectorial", marker(d) );
        this->setMarkerNeumannBC(super_type::NeumannBCShape::VECTORIAL,marker(d),markerList);

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Neumann_vectorial", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
    }
    M_bcNeumannTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<super_type::nDim>( "velocity", "Neumann_tensor2" );
    for( auto const& d : M_bcNeumannTensor2 )
    {
        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "velocity", "Neumann_tensor2", marker(d) );
        this->setMarkerNeumannBC(super_type::NeumannBCShape::TENSOR2,marker(d),markerList);

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Neumann_tensor2", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
    }

    M_bcPressure = this->modelProperties().boundaryConditions().getScalarFields( "pressure", "Dirichlet" );
    for( auto const& d : M_bcPressure )
    {
        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "pressure", "Dirichlet", marker(d) );
        this->setMarkerPressureBC(marker(d),markerList);

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "pressure", "Dirichlet", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
    }
    for( std::string const& bcMarker : this->modelProperties().boundaryConditions().markers("fluid", "slip") )
    {
        this->addMarkerSlipBC( bcMarker );
        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "slip", bcMarker, "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        this->addMarkerALEMeshBC(bcTypeMeshALE,bcMarker);
    }
    for( std::string const& bcMarker : this->modelProperties().boundaryConditions().markers("fluid", "outlet") )
    {
        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "fluid", "outlet", bcMarker, "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");

        std::string typeOutlet = soption(_name="fluid-outlet.type", _prefix=this->prefix());//"free";
        std::pair<bool,std::string> typeOutletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "outlet", bcMarker, "model" );
        if ( typeOutletRead.first )
        {
            typeOutlet = typeOutletRead.second;
            CHECK( typeOutlet == "free" || typeOutlet == "windkessel" ) << "invalid outlet model " << typeOutlet;
        }
        std::string typeCouplingWindkesselOutlet = soption(_name="fluid-outlet.windkessel.coupling", _prefix=this->prefix());
        std::pair<bool,std::string> typeCouplingWindkesselOutletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "outlet", bcMarker, "windkessel_coupling" );
        if ( typeCouplingWindkesselOutletRead.first )
        {
            typeCouplingWindkesselOutlet = typeCouplingWindkesselOutletRead.second;
            CHECK( typeCouplingWindkesselOutlet == "implicit" || typeCouplingWindkesselOutlet == "explicit" ) << "invalid windkessel coupling type " << typeCouplingWindkesselOutlet;
        }
        std::pair<bool,double> WindkesselRdRead = this->modelProperties().boundaryConditions().dparam( "fluid", "outlet", bcMarker, "windkessel_Rd" );
        std::pair<bool,double> WindkesselRpRead = this->modelProperties().boundaryConditions().dparam( "fluid", "outlet", bcMarker, "windkessel_Rp" );
        std::pair<bool,double> WindkesselCdRead = this->modelProperties().boundaryConditions().dparam( "fluid", "outlet", bcMarker, "windkessel_Cd" );
        double WindkesselRd = ( WindkesselRdRead.first )? WindkesselRdRead.second : 1.;
        double WindkesselRp = ( WindkesselRpRead.first )? WindkesselRpRead.second : 1.;
        double WindkesselCd = ( WindkesselCdRead.first )? WindkesselCdRead.second : 1.;

        std::tuple<std::string,double,double,double> windkesselParam = std::make_tuple(typeCouplingWindkesselOutlet,WindkesselRd,WindkesselRp,WindkesselCd);

        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "fluid", "outlet", bcMarker );
        for (std::string const& currentMarker : markerList )
        {
            this->M_fluidOutletsBCType.push_back(std::make_tuple(currentMarker,typeOutlet, windkesselParam ));
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
        }
    }
    for( std::string const& bcMarker : this->modelProperties().boundaryConditions().markers("fluid", "inlet") )
    {
        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "fluid", "inlet", bcMarker, "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");

        std::string shapeInlet;
        std::pair<bool,std::string> shapeInletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "inlet", bcMarker, "shape" );
        if ( shapeInletRead.first )
        {
            shapeInlet = shapeInletRead.second;
            CHECK( shapeInlet == "constant" || shapeInlet == "parabolic" ) << "invalid inlet shape " << shapeInlet;
        }
        else
            CHECK( false ) << "inlet shape not given";

        std::string constraintInlet;
        std::pair<bool,std::string> constraintInletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "inlet", bcMarker, "constraint" );
        if ( constraintInletRead.first )
        {
            constraintInlet = constraintInletRead.second;
            CHECK( constraintInlet == "velocity_max" || constraintInlet == "flow_rate" ) << "invalid inlet constraint " << constraintInlet;
        }
        else
            CHECK( false ) << "inlet constraint not given";

        std::string fullTypeInlet = (boost::format("%1%_%2%")%constraintInlet %shapeInlet).str();

        std::string exprFluidInlet;
        std::pair<bool,std::string> exprFluidInletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "inlet", bcMarker, "expr" );
        if ( exprFluidInletRead.first )
            exprFluidInlet = exprFluidInletRead.second;
        else
            CHECK( false ) << "inlet expr not given";

        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "fluid", "inlet", bcMarker );
        for (std::string const& currentMarker : markerList )
        {
            this->M_fluidInletDesc.push_back(std::make_tuple(currentMarker,fullTypeInlet, expr<2>( exprFluidInlet,"",this->worldComm(),this->directoryLibSymbExpr() )) );
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
        }
    }

    M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<super_type::nDim>( "fluid", "VolumicForces" );

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::loadConfigPostProcess()
{
    if ( this->modelProperties().postProcess().find("Fields") != this->modelProperties().postProcess().end() )
        for ( auto const& o : this->modelProperties().postProcess().find("Fields")->second )
        {
            if ( o == "velocity" || o == "all" ) this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Velocity );
            if ( o == "pressure" || o == "all" ) this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Pressure );
            if ( o == "displacement" || o == "all" ) this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Displacement );
            if ( o == "vorticity" || o == "all" ) this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Vorticity );
            if ( o == "stress" || o == "normal-stress" || o == "all" ) this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::NormalStress );
            if ( o == "wall-shear-stress" || o == "all" ) this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::WallShearStress );
            if ( o == "density" || o == "all" ) this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Density );
            if ( o == "viscosity" || o == "all" ) this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Viscosity );
            if ( o == "pid" || o == "all" ) this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Pid );

            if ( o == "alemesh" /*|| o == "all"*/ ) this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::ALEMesh );
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
                {
                    auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",myBcDesc.marker() ) );
                    auto const& listMarkerFaces = std::get<0>( ret );
                    ExpressionStringAtMarker myBcDesc2( myBcDesc );
                    myBcDesc2.setMeshMarkers( listMarkerFaces );
                    bcPrecPCD["velocity"]["Dirichlet"].push_back( myBcDesc2 );
                }
            }
            // For weak Dirichlet (Nitche,Magrange Multiplier ) ???
            // TODO Dirchlet component

            auto itFindNeumannScalType = itFindFieldVelocity->second.find("Neumann_scalar");
            if ( itFindNeumannScalType != itFindFieldVelocity->second.end() )
            {
                for ( auto const& myBcDesc : itFindNeumannScalType->second )
                {
                    auto markList = this->markerNeumannBC( super_type::NeumannBCShape::SCALAR,myBcDesc.marker() );
                    if ( markList.empty() ) continue;
                    ExpressionStringAtMarker myBcDesc2( myBcDesc );
                    myBcDesc2.setMeshMarkers( markList );
                    bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc2 );
                }
            }
            auto itFindNeumannVecType = itFindFieldVelocity->second.find("Neumann_vectorial");
            if ( itFindNeumannVecType != itFindFieldVelocity->second.end() )
            {
                for ( auto const& myBcDesc : itFindNeumannVecType->second )
                {
                    auto markList = this->markerNeumannBC( super_type::NeumannBCShape::VECTORIAL,myBcDesc.marker() );
                    if ( markList.empty() ) continue;
                    ExpressionStringAtMarker myBcDesc2( myBcDesc );
                    myBcDesc2.setMeshMarkers( markList );
                    bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc2 );
                }
            }
            auto itFindNeumannTensor2Type = itFindFieldVelocity->second.find("Neumann_tensor2");
            if ( itFindNeumannTensor2Type != itFindFieldVelocity->second.end() )
            {
                for ( auto const& myBcDesc : itFindNeumannTensor2Type->second )
                {
                    auto markList = this->markerNeumannBC( super_type::NeumannBCShape::TENSOR2,myBcDesc.marker() );
                    if ( markList.empty() ) continue;
                    ExpressionStringAtMarker myBcDesc2( myBcDesc );
                    myBcDesc2.setMeshMarkers( markList );
                    bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc2 );
                }
            }
        }
#if 0
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
#else
        if ( !this->M_fluidOutletsBCType.empty() )
        {
            std::list<std::string> markList;
            for ( auto const& bcOutlet : this->M_fluidOutletsBCType )
                markList.push_back( std::get<0>(bcOutlet) );
            ExpressionStringAtMarker myBcDesc2( std::make_tuple( "expression","wind","0","","" ) );
            myBcDesc2.setMeshMarkers( markList );
            bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc2 );
        }
#endif

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
                        std::ostringstream ostrMarkers;
                        ostrMarkers << "(";
                        for ( std::string const& mark : c.meshMarkers() )
                            ostrMarkers << mark << " ";
                        ostrMarkers << ")";
                        if ( c.hasExpression2() )
                            std::cout << "  . boundary  " << c.marker() << " " << ostrMarkers.str() << " expr : " << c.expression1() << " expr2:" << c.expression2() << "\n";
                        else
                            std::cout << "  . boundary  " << c.marker() << " " << ostrMarkers.str() << " expr : " << c.expression() << "\n";
                    }
                }
            }
        }
#endif
        CHECK( this->algebraicFactory()->preconditionerTool()->matrix() ) << "no matrix define in preconditionerTool";
        // auto myalpha = (this->isStationary())? 0 : this->densityViscosityModel()->cstRho()*this->timeStepBDF()->polyDerivCoefficient(0);
        auto myalpha = (!this->isStationary())*idv(this->densityViscosityModel()->fieldRho())*this->timeStepBDF()->polyDerivCoefficient(0);

        typedef typename super_type::space_fluid_type space_type;
        typedef typename super_type::space_densityviscosity_type properties_space_type;

        boost::shared_ptr< PreconditionerBlockNS<space_type, properties_space_type> > a_blockns = Feel::blockns( _space=this->functionSpace(),
                                        _properties_space=this->densityViscosityModel()->fieldDensityPtr()->functionSpace(),
                                        _type=soption(_prefix=this->prefix(),_name="blockns.type"),//"PCD",
                                        _bc=bcPrecPCD,
                                        _matrix=this->algebraicFactory()->preconditionerTool()->matrix(),
                                        _prefix="velocity",
                                        _mu=idv(this->densityViscosityModel()->fieldMu()),
                                        _rho=idv(this->densityViscosityModel()->fieldRho()),
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

        typedef typename super_type::space_fluid_type space_type;
        typedef typename super_type::space_densityviscosity_type properties_space_type;

        boost::shared_ptr< PreconditionerBlockNS<space_type, properties_space_type> > myPrecBlockNs =
            boost::dynamic_pointer_cast< PreconditionerBlockNS<space_type, properties_space_type> >( this->algebraicFactory()->preconditionerTool()->inHousePreconditioners( "blockns" ) );

        auto myalpha = (!this->isStationary())*idv(this->densityViscosityModel()->fieldRho())*this->timeStepBDF()->polyDerivCoefficient(0);
        myPrecBlockNs->setAlpha( myalpha );
        myPrecBlockNs->setMu( idv(this->densityViscosityModel()->fieldMu()) );
        myPrecBlockNs->setRho( idv(this->densityViscosityModel()->fieldRho()) );

        if ( this->modelName() == "Stokes" )
        {
            myPrecBlockNs->update( mat );
        }
        else if ( ( this->modelName() == "Navier-Stokes" && this->solverName() == "Oseen" ) || this->modelName() == "Oseen" )
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
        else if ( this->modelName() == "Navier-Stokes" )
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
    this->modelProperties().parameters().updateParameterValues();

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    M_bcDirichlet.setParameterValues( paramValues );
    for ( auto & bcDirComp : M_bcDirichletComponents )
        bcDirComp.second.setParameterValues( paramValues );
    M_bcNeumannScalar.setParameterValues( paramValues );
    M_bcNeumannVectorial.setParameterValues( paramValues );
    M_bcNeumannTensor2.setParameterValues( paramValues );
    M_bcPressure.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );
    this->updateFluidInletVelocity();

    if ( this->algebraicFactory() && this->algebraicFactory()->preconditionerTool()->hasInHousePreconditioners( "blockns" ) )
    {
        typedef typename super_type::space_fluid_type space_type;
        typedef typename super_type::space_densityviscosity_type properties_space_type;

        boost::shared_ptr< PreconditionerBlockNS<space_type, properties_space_type> > myPrecBlockNs =
            boost::dynamic_pointer_cast< PreconditionerBlockNS<space_type, properties_space_type> >( this->algebraicFactory()->preconditionerTool()->inHousePreconditioners( "blockns" ) );
        myPrecBlockNs->setParameterValues( paramValues );
    }

    if ( this->M_useThermodynModel && this->M_useGravityForce )
        this->M_thermodynModel->updateParameterValues();

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

    size_type rowStartInVector = this->rowStartInVector();
    //auto const& u = this->fieldVelocity();
    auto Xh = this->functionSpace();
    auto up = Xh->element( U, rowStartInVector );
    auto u = up.template element<0>();
    auto mesh = this->mesh();

    if (!Xh->worldsComm()[0].isActive()) // only on Velocity Proc
        return;

    // store markers for each entities in order to apply strong bc with priority (points erase edges erace faces)
    std::map<std::string, std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string> > > mapMarkerBCToEntitiesMeshMarker;
    for( auto const& d : M_bcDirichlet )
    {
        mapMarkerBCToEntitiesMeshMarker[marker(d)] =
            detail::distributeMarkerListOnSubEntity(mesh,{
                    /**/this->markerDirichletBCByNameId( "elimination",marker(d) ),
                        this->markerDirichletBCByNameId( "lm",marker(d) ) } );
    }
    std::map<std::pair<std::string,ComponentType>, std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string> > > mapCompMarkerBCToEntitiesMeshMarker;
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            mapCompMarkerBCToEntitiesMeshMarker[std::make_pair(marker(d),comp)] =
                detail::distributeMarkerListOnSubEntity(mesh, {
                        /**/this->markerDirichletBCByNameId( "elimination",marker(d), comp ),
                            this->markerDirichletBCByNameId( "lm", marker(d), comp )  } );
        }
    }

    // modif vector with bc define on topological faces
    for( auto const& d : M_bcDirichlet )
    {
        auto const& listMarkerFaces = std::get<0>( mapMarkerBCToEntitiesMeshMarker.find( marker(d) )->second );
        if ( !listMarkerFaces.empty() )
            u.on(_range=markedfaces(mesh,listMarkerFaces ),
                 _expr=expression(d) );
    }
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto const& listMarkerFaces = std::get<0>(  mapCompMarkerBCToEntitiesMeshMarker.find( std::make_pair(marker(d),comp) )->second );
            if ( !listMarkerFaces.empty() )
                u[comp].on(_range=markedfaces(mesh,listMarkerFaces ),
                           _expr=expression(d) );
        }
    }

    // modif vector with bc define on edges (only 3d)
    for( auto const& d : M_bcDirichlet )
    {
        auto const& listMarkerEdges = std::get<1>( mapMarkerBCToEntitiesMeshMarker.find( marker(d) )->second );
        if ( !listMarkerEdges.empty() )
            u.on(_range=markededges(mesh,listMarkerEdges ),
                 _expr=expression(d) );
    }
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto const& listMarkerEdges = std::get<1>(  mapCompMarkerBCToEntitiesMeshMarker.find( std::make_pair(marker(d),comp) )->second );
            if ( !listMarkerEdges.empty() )
                u[comp].on(_range=markededges(mesh,listMarkerEdges ),
                           _expr=expression(d) );
        }
    }

    // modif vector with bc define on points
    for( auto const& d : M_bcDirichlet )
    {
        auto const& listMarkerPoints = std::get<2>( mapMarkerBCToEntitiesMeshMarker.find( marker(d) )->second );
        if ( !listMarkerPoints.empty() )
            u.on(_range=markedfaces(mesh,listMarkerPoints ),
                 _expr=expression(d) );
    }
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto const& listMarkerPoints = std::get<2>(  mapCompMarkerBCToEntitiesMeshMarker.find( std::make_pair(marker(d),comp) )->second );
            if ( !listMarkerPoints.empty() )
                u[comp].on(_range=markedfaces(mesh,listMarkerPoints ),
                           _expr=expression(d) );
        }
    }

    double t1=timerBCnewton.elapsed();
    this->log("FluidMechanics","updateInitialNewtonSolutionBCDirichlet","finish in "+(boost::format("%1% s") % t1).str() );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    using namespace Feel::vf;

    if ( !this->hasDirichletBC() ) return;
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletLinearPDE", "start",
                                               this->worldComm(),this->verboseAllProc());
    boost::mpi::timer timerBClinear;

    auto Xh = this->functionSpace();
    auto mesh = this->mesh();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    // special hack : use bilinearform only on velocity
    // use operotor on(..) doesnt work on composite space (velocity+pressure) ( see bilinearform.hpp l731 no support component)
    auto bilinearFormComp = form2( _test=this->functionSpaceVelocity(),_trial=this->functionSpaceVelocity(),_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldVelocity();

    // store markers for each entities in order to apply strong bc with priority (points erase edges erace faces)
    std::map<std::string, std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string> > > mapMarkerBCToEntitiesMeshMarker;
    for( auto const& d : M_bcDirichlet )
    {
        mapMarkerBCToEntitiesMeshMarker[marker(d)] =
            detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d) ) );
    }
    std::map<std::pair<std::string,ComponentType>, std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string> > > mapCompMarkerBCToEntitiesMeshMarker;
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            mapCompMarkerBCToEntitiesMeshMarker[std::make_pair(marker(d),comp)] =
                detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d),comp ) );
        }
    }

    // apply on() with bc define on topological faces
    for( auto const& d : M_bcDirichlet )
    {
        auto const& listMarkerFaces = std::get<0>( mapMarkerBCToEntitiesMeshMarker.find( marker(d) )->second );
        if ( !listMarkerFaces.empty() )
            bilinearForm +=
                on( _range=markedfaces( mesh, listMarkerFaces ),
                    _element=u, _rhs=F, _expr=expression(d) );
    }
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto const& listMarkerFaces = std::get<0>(  mapCompMarkerBCToEntitiesMeshMarker.find( std::make_pair(marker(d),comp) )->second );
            if ( !listMarkerFaces.empty() )
                bilinearFormComp +=
                    on( _range=markedfaces( mesh, listMarkerFaces ),
                        _element=this->M_Solution->template element<0>()[comp], //u[comp],
                        _rhs=F, _expr=expression(d) );
        }
    }

    // apply on() with bc define on edges (only 3d)
    for( auto const& d : M_bcDirichlet )
    {
        auto const& listMarkerEdges = std::get<1>( mapMarkerBCToEntitiesMeshMarker.find( marker(d) )->second );
        if ( !listMarkerEdges.empty() )
            bilinearForm +=
                on( _range=markedfaces( mesh, listMarkerEdges ),
                    _element=u, _rhs=F, _expr=expression(d) );
    }
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto const& listMarkerEdges = std::get<1>(  mapCompMarkerBCToEntitiesMeshMarker.find( std::make_pair(marker(d),comp) )->second );
            if ( !listMarkerEdges.empty() )
                bilinearFormComp +=
                    on( _range=markedfaces( this->mesh(), listMarkerEdges ),
                        _element=this->M_Solution->template element<0>()[comp], //u[comp],
                        _rhs=F, _expr=expression(d) );
        }
    }

    // apply on() with bc define on points
    for( auto const& d : M_bcDirichlet )
    {
        auto const& listMarkerPoints = std::get<2>( mapMarkerBCToEntitiesMeshMarker.find( marker(d) )->second );
        if ( !listMarkerPoints.empty() )
            bilinearForm +=
                on( _range=markedfaces( mesh, listMarkerPoints ),
                    _element=u,
                    _rhs=F, _expr=expression(d) );
    }
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto const& listMarkerPoints = std::get<2>(  mapCompMarkerBCToEntitiesMeshMarker.find( std::make_pair(marker(d),comp) )->second );
            if ( !listMarkerPoints.empty() )
                bilinearFormComp +=
                    on( _range=markedfaces( mesh, listMarkerPoints ),
                        _element=this->M_Solution->template element<0>()[comp], //u[comp],
                        _rhs=F, _expr=expression(d) );
        }
    }


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
    if ( !this->hasDirichletBC() ) return;
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletJacobian", "start",
                                               this->worldComm(),this->verboseAllProc());
    using namespace Feel::vf;

    boost::timer btimeStrongCL;

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldVelocity();
    //auto RBis = this->backend()->newVector( J->mapRowPtr() );

    for( auto const& d : M_bcDirichlet )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d) ) );
        auto const& listMarkerFaces = std::get<0>( ret );
        auto const& listMarkerEdges = std::get<1>( ret );
        auto const& listMarkerPoints = std::get<2>( ret );
        auto exprUsed = vf::zero<super_type::nDim,1>();// 0*vf::one();
        if ( !listMarkerFaces.empty() )
            bilinearForm +=
                on( _range=markedfaces( mesh, listMarkerFaces ),
                    _element=u, _rhs=RBis, _expr=exprUsed );
        if ( !listMarkerEdges.empty() )
            bilinearForm +=
                on( _range=markededges(mesh, listMarkerEdges),
                    _element=u,_rhs=RBis,_expr=exprUsed );
        if ( !listMarkerPoints.empty() )
            bilinearForm +=
                on( _range=markedpoints(mesh, listMarkerPoints),
                    _element=u,_rhs=RBis,_expr=exprUsed );
    }
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto bilinearFormComp = form2( _test=this->functionSpaceVelocity(),_trial=this->functionSpaceVelocity(),_matrix=J,
                                           _rowstart=this->rowStartInMatrix(),
                                           _colstart=this->colStartInMatrix() );
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),comp ) );
            auto const& listMarkerFaces = std::get<0>( ret );
            auto const& listMarkerEdges = std::get<1>( ret );
            auto const& listMarkerPoints = std::get<2>( ret );
            auto exprUsed = cst(0.);
            if ( !listMarkerFaces.empty() )
                bilinearFormComp +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=this->M_Solution->template element<0>()[comp], //u[comp],
                        _rhs=RBis,_expr=exprUsed );
            if ( !listMarkerEdges.empty() )
                bilinearFormComp +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=this->M_Solution->template element<0>()[comp], //u[comp],
                        _rhs=RBis,_expr=exprUsed );
            if ( !listMarkerPoints.empty() )
                bilinearFormComp +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=this->M_Solution->template element<0>()[comp], //u[comp],
                        _rhs=RBis,_expr=exprUsed );
        }
    }

    std::ostringstream ostr3;ostr3<<btimeStrongCL.elapsed()<<"s";
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletJacobian", "finish in "+ostr3.str(),
                                               this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletResidual(vector_ptrtype& R) const
{
    if ( !this->hasDirichletBC() ) return;
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateBCStrongDirichletResidual", "start",
                                               this->worldComm(),this->verboseAllProc());
    boost::mpi::timer timerBCresidu;

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    size_type rowStartInVector = this->rowStartInVector();
    //auto const& u = this->fieldVelocity();
    auto up = Xh->element(R,rowStartInVector);
    auto u = up.template element<0>();

    //R->close();
    if (!Xh->worldsComm()[0].isActive()) // only on Velocity Proc
        return;


    // Zero because is know
    for( auto const& d : M_bcDirichlet )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",marker(d) ) );
        auto const& listMarkerFaces = std::get<0>( ret );
        auto const& listMarkerEdges = std::get<1>( ret );
        auto const& listMarkerPoints = std::get<2>( ret );
        auto exprUsed = vf::zero<super_type::nDim,1>();// 0*vf::one();
        if ( !listMarkerFaces.empty() )
            u.on(_range=markedfaces(mesh,listMarkerFaces ),
                 _expr=exprUsed );
        if ( !listMarkerEdges.empty() )
            u.on(_range=markededges(mesh,listMarkerEdges),
                 _expr=exprUsed );
        if ( !listMarkerPoints.empty() )
            u.on(_range=markedpoints(mesh,listMarkerPoints),
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
            auto exprUsed = vf::zero<1,1>();// cst(0.);
            if ( !listMarkerFaces.empty() )
                u[comp].on(_range=markedfaces(mesh,listMarkerFaces ),
                           _expr=exprUsed );
            if ( !listMarkerEdges.empty() )
                u[comp].on(_range=markededges(mesh,listMarkerEdges),
                           _expr=exprUsed );
            if ( !listMarkerPoints.empty() )
                u[comp].on(_range=markedpoints(mesh,listMarkerPoints),
                           _expr=exprUsed );
        }
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
    if ( M_bcNeumannScalar.empty() && M_bcNeumannVectorial.empty() && M_bcNeumannTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpace(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldVelocity();
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannResidual( vector_ptrtype& R ) const
{
    if ( M_bcNeumannScalar.empty() && M_bcNeumannVectorial.empty() && M_bcNeumannTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpace(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldVelocity();
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCPressureLinearPDE( vector_ptrtype& F ) const
{
    if ( M_bcPressure.empty() ) return;
    auto myLinearForm = form1( _test=this->functionSpace(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldVelocity();
    for( auto const& d : M_bcPressure )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerPressureBC(marker(d)) ),
                       _expr= -expression(d)*trans(N())*id(v),
                       _geomap=this->geomap() );
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCPressureResidual( vector_ptrtype& R ) const
{
    if ( M_bcPressure.empty() ) return;
    auto myLinearForm = form1( _test=this->functionSpace(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldVelocity();
    for( auto const& d : M_bcPressure )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerPressureBC(marker(d)) ),
                       _expr= expression(d)*trans(N())*id(v),
                       _geomap=this->geomap() );
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletLagMultLinearPDE( vector_ptrtype& F ) const
{
    if ( M_bcDirichlet.empty() ) return;

    auto lambdaBC = this->XhDirichletLM()->element();
    size_type startBlockIndexDirichletLM = this->startBlockIndexFieldsInMatrix().find("dirichletlm")->second;
    for( auto const& d : M_bcDirichlet )
        form1( _test=this->XhDirichletLM(),_vector=F,
               _rowstart=this->rowStartInVector()+startBlockIndexDirichletLM ) +=
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
    size_type startBlockIndexDirichletLM = this->startBlockIndexFieldsInMatrix().find("dirichletlm")->second;
    for( auto const& d : M_bcDirichlet )
        form1( _test=this->XhDirichletLM(),_vector=R,
               _rowstart=this->rowStartInVector()+startBlockIndexDirichletLM ) +=
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
