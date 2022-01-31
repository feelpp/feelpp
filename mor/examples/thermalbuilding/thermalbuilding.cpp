//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Vincent Chabannes <vincent.chabannes@feelpp.org>
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 12 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!

#include <feel/feelmor/crbplugin.hpp>
#include "thermalbuilding.hpp"

namespace Feel
{

po::options_description
makeThermalBuildingOptions()
{
    po::options_description myoptions( "ThermalBuilding options" );
    myoptions.add_options()
        ( "thermal-building.gamma", po::value<double>()->default_value( 10 ), "penalisation term" )
        ( "thermal-building.conductivity.air", po::value<double>()->default_value( /*0.0262*4*/ /*2.9*/ 1. ), "conductivity.air" )
        ( "thermal-building.conductivity.internal-walls", po::value<double>()->default_value( 0.25 ), "conductivity.internal-walls" )
        ( "thermal-building.conductivity.internal-doors", po::value<double>()->default_value( 0.13 ), "conductivity.internal-doors" )
        ;
    return myoptions;
}
AboutData
makeThermalBuildingAbout( std::string const& str )
{
    Feel::AboutData about( /*AppName  */ str.c_str(),
                           /*ProgName */ str.c_str(),
                           /*Version  */ "0.1",
                           /*ShortDesc*/ "Thermal Building Application",
                           /*Licence  */ Feel::AboutData::License_GPL,
                           /*Copyright*/ "Copyright (c) 2016 Feel++ Consortium" );
    return about;
}

ThermalBuilding::ThermalBuilding()
    :
    super_type( "ThermalBuilding" )
{}

void
ThermalBuilding::initBetaQ()
{
    M_Qa=2;//1;
    M_Nl=2;
    M_Ql.resize( 2 );
    M_Ql[0]=7;//6;
    M_Ql[1]=1;
    this->M_betaAq.resize( M_Qa );
    this->M_betaFq.resize( M_Nl );
    for(int i=0; i<M_Nl; i++)
    {
        int ql=M_Ql[i];
        this->M_betaFq[i].resize( ql );
    }
}

void
ThermalBuilding::initDataStructureAffineDecomposition()
{
    this->initBetaQ();
#if 0
    this->M_betaAq.resize( M_Qa );
    this->M_betaFq.resize( M_Nl );
    for(int i=0; i<M_Nl; i++)
    {
        int ql=M_Ql[i];
        this->M_betaFq[i].resize( ql );
    }
#endif
    this->M_Aq.resize( M_Qa );
    this->M_Fq.resize( M_Nl );
    for(int l=0; l<M_Nl; l++)
    {
        this->M_Fq[l].resize( M_Ql[l] );
    }

}

ThermalBuilding::super_type::betaq_type
ThermalBuilding::computeBetaQ( parameter_type const& mu , double time , bool only_terms_time_dependent )
{
    if ( this->M_betaAq.empty() )
        this->initBetaQ();

    this->M_betaAq[0] = 1;
    this->M_betaAq[1] = mu(6);
    for (int k = 0;k<5;++k)
        this->M_betaFq[0][k] = mu( k );
    this->M_betaFq[0][5] = mu(5)*mu(6);
    this->M_betaFq[0][6] = mu(5);
    this->M_betaFq[1][0] = 1.;
    return boost::make_tuple( this->M_betaAq, this->M_betaFq );

}
ThermalBuilding::super_type::betaq_type
ThermalBuilding::computeBetaQ( parameter_type const& mu )
{
    return this->computeBetaQ( mu,0, false );
}

void
ThermalBuilding::initModel()
{
    CHECK( is_linear && !is_time_dependent ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&NonLinear) == NonLinear) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&TimeDependent) == TimeDependent) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";

    auto mesh = loadMesh( _mesh=new ThermalBuildingConfig::mesh_type,
                          _update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES),
                          _savehdf5=0 );

    this->setFunctionSpaces( ThermalBuildingConfig::space_type::New( mesh ) );
    //this->setSymbolicExpressionBuildDir("$repository/crb/thermalbuilding/symbolicexpr/");

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";
    }

    M_Qa=2;//1;
    M_Nl=2;
    M_Ql.resize( 2 );
    M_Ql[0]=7;//6;
    M_Ql[1]=1;
    this->initDataStructureAffineDecomposition();

    Dmu->setDimension( 7 );
#if 1
    // 0.24275
    /// [parameters]
    auto mu_min = Dmu->element();
    //mu_min << 310,310,310,310,310,280,0.14;
    mu_min << 300,300,300,300,300,270,0.14;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 340,340,340,340,340,290,0.29;
    Dmu->setMax( mu_max );
#else
    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    CHECK( mu_min.size() == this->crbModelParameters().nCrbParameters() ) << "invalid size compatibility";
    for ( int k=0; k<this->crbModelParameters().nCrbParameters(); ++k )
    {
        mu_min(k) = this->crbModelParameters().crbParameters(k).min();
        mu_max(k) = this->crbModelParameters().crbParameters(k).max();
    }
    Dmu->setMin( mu_min );
    Dmu->setMax( mu_max );
#endif

    this->assembleData();
}
void
ThermalBuilding::assembleData()
{
    auto mesh = this->Xh->mesh();
    auto u = this->Xh->element();
    auto v = this->Xh->element();

    std::vector<std::string> volumeMarkers = { "air", "internal-walls", "internal-doors" };
    for ( std::string const& markerVol : volumeMarkers )
        CHECK( mesh->hasMarker( markerVol ) ) << "mesh does not have the volume marker : " <<  markerVol;

    double penaldir= doption(_name="thermal-building.gamma");

    // init matrix/vectors
    for (int k = 0 ; k<this->M_Aq.size() ; ++k)
        this->M_Aq[k] = backend()->newMatrix( this->Xh, this->Xh );
    for (int k = 0 ; k<this->M_Fq[0].size() ; ++k)
        this->M_Fq[0][k] = backend()->newVector( this->Xh );
    this->M_Fq[1][0] = backend()->newVector( this->Xh );

    std::map<std::string,double> conductivity;
    conductivity["air"] = doption(_name="thermal-building.conductivity.air");
    conductivity["internal-walls"] = doption(_name="thermal-building.conductivity.internal-walls");
    conductivity["internal-doors"] = doption(_name="thermal-building.conductivity.internal-doors");

    for ( std::string const& marker : volumeMarkers )
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[0]) +=
            integrate( _range= markedelements( mesh,marker ),
                       _expr= conductivity[marker]*gradt( u )*trans( grad( v ) ) );


    // weak Dirichlet condition on heaters
    std::vector<std::string> heatersSurfMarkers = { "heater-livingroom","heater-kitchen","heater-bedroom1", "heater-bedroom2","heater-bathroom" };
    form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[0]) +=
        integrate( _range=markedfaces( mesh,  std::list<std::string>( heatersSurfMarkers.begin(),heatersSurfMarkers.end() ) ),
                   _expr= conductivity["air"]*( -(gradt(u)*id(v)+grad(v)*idt(u))*N() + penaldir*idt(u)*id(v)/hFace() ) );

    for (int k=0;k<heatersSurfMarkers.size();++k )
        form1(_test=this->Xh, _vector=this->M_Fq[0][k]) +=
            integrate( _range=markedfaces( mesh,heatersSurfMarkers[k] ),
                       _expr=conductivity["air"]*( -grad(v)*N() + penaldir*id(v)/hFace() ) );


    // rsi + enduit + parpain + polystirene + ba13 + rse
    // 1.0/(0.06+0.01/0.5 + 0.3/0.8 + 0.20/0.032 +0.016/0.313 +0.14)

    // Robin condition on contatc exterior
    std::map<std::string,double> heatTransferCoefficient;
    heatTransferCoefficient["front-door"] = 1.0/(0.06+0.06/0.150+0.1/0.029+0.14);
#if 0
    std::vector<std::string> markerContactExterior = { "exterior-walls","front-door" };
    for ( std::string const& mark : markerContactExterior )
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[0]) +=
            integrate( _range=markedfaces(this->mesh(), mark ),
                       _expr= heatTransferCoefficient[mark]*idt(v)*id(v) );

    int idFq0ExternalTemp = heatersSurfMarkers.size();
    for ( std::string const& mark : markerContactExterior )
        form1(_test=this->Xh, _vector=this->M_Fq[0][idFq0ExternalTemp]) +=
            integrate( _range=markedfaces(this->mesh(),mark),
                       _expr= heatTransferCoefficient[mark]*id(v) );
#else
    form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[0]) +=
        integrate( _range=markedfaces(this->mesh(), "front-door" ),
                   _expr= heatTransferCoefficient["front-door"]*idt(v)*id(v) );

    form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[1]) +=
        integrate( _range=markedfaces(this->mesh(), "exterior-walls" ),
                   _expr= idt(v)*id(v) );

    //int idFq0ExternalTemp = heatersSurfMarkers.size();
    form1(_test=this->Xh, _vector=this->M_Fq[0][5]) +=
        integrate( _range=markedfaces(this->mesh(),"exterior-walls" ),
                   _expr= /*heatTransferCoefficient[mark]**/id(v) );

    form1(_test=this->Xh, _vector=this->M_Fq[0][6]) +=
        integrate( _range=markedfaces(this->mesh(),"front-door" ),
                   _expr= heatTransferCoefficient["front-door"]*id(v) );

#endif
    // define inner product
    auto energy = backend()->newMatrix(_test=Xh,_trial=Xh);
    energy->addMatrix(1.,this->M_Aq[0] );
    energy->addMatrix( 0.18,this->M_Aq[1] );
    energy->close();
    this->addEnergyMatrix( energy );


    // output
    double areaAir = integrate( _range=markedelements( mesh, "air" ), _expr= cst(1.) ).evaluate()(0,0);
    form1(_test=this->Xh, _vector=this->M_Fq[1][0])
        += integrate( _range=markedelements( mesh, "air" ), _expr= id( v )/areaAir );

}
double
ThermalBuilding::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{
    //CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    auto mesh = Xh->mesh();
    double output=0;
    // right hand side (compliant)
    if ( output_index == 0 )
    {
        //output  = integrate( markedfaces( mesh,"BR" ), -mu(0)*expr(soption("functions.f"))*(gradv(u)*N()+doption("gamma")*idv(u)) ).evaluate()(0,0)
        //    + integrate( markedfaces( mesh,"BL" ), -mu(1)*expr(soption("functions.g"))*(gradv(u)*N()+doption("gamma")*idv(u)) ).evaluate()(0,0);
    }
    else if ( output_index == 1 )
    {
        output = mean(markedelements(mesh,"air" ),idv(u))(0,0);
        std::cout << " ThermalBuilding::output " << output << "\n";
    }
    // else if ( output_index == 2 )
    // {
    //     output = mean(elements(mesh),idv(u)).evaluat()(0,0);
    // }
    else
        throw std::logic_error( "[ThermalBuilding::output] error with output_index : only 0 or 1 " );
    return output;

}
FEELPP_CRB_PLUGIN( ThermalBuilding, thermalbuilding )
} // namespace Feel
