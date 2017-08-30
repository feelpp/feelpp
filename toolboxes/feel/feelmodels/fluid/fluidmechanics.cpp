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
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::self_ptrtype
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::New( std::string const& prefix, bool buildMesh,
                                         WorldComm const& worldComm, std::string const& subPrefix,
                                         std::string const& rootRepository )
{
    return boost::make_shared<self_type>( prefix, buildMesh, worldComm, subPrefix, rootRepository );

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
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::solve()
{
    this->modelProperties().parameters().updateParameterValues();

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    //this->modelProperties().materials().setParameterValues( paramValues );
    this->M_bcDirichlet.setParameterValues( paramValues );
    for ( auto & bcDirComp : this->M_bcDirichletComponents )
        bcDirComp.second.setParameterValues( paramValues );
    this->M_bcNeumannScalar.setParameterValues( paramValues );
    this->M_bcNeumannVectorial.setParameterValues( paramValues );
    this->M_bcNeumannTensor2.setParameterValues( paramValues );
    this->M_bcPressure.setParameterValues( paramValues );
    this->M_volumicForcesProperties.setParameterValues( paramValues );
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


} // namespace FeelModels

} // namespace Feel
