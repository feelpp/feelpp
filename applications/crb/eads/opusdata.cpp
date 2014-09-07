/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-11-15

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file opusdata.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-11-15
 */
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>

#include <feel/options.hpp>
#include <opusdefs.hpp>
#include <opusdata.hpp>
#include <opuscomponent.hpp>
#include <feel/feeldiscr/bdf2.hpp>

namespace Feel
{
Feel::po::options_description
OpusData::makeOptions()
{
    Feel::po::options_description opusoptions( "Opus benchmark #1 options" );
    opusoptions.add_options()
    ( "steady", Feel::po::value<bool>()->default_value( 1 ), "compute the steady state" )
    ( "stab", Feel::po::value<bool>()->default_value( 0 ), "compute with stabilisation" )

    ( "d", Feel::po::value<int>()->default_value( 2 ), "time step value" )

    ( "order-geo", Feel::po::value<int>()->default_value( 2 ), "order of geometry" )
    ( "order-time", Feel::po::value<int>()->default_value( 1 ), "order of time discretisation (time)" )
    ( "order-temp", Feel::po::value<int>()->default_value( 2 ), "order of spatial discretisation (temperature)" )
    ( "order-u", Feel::po::value<int>()->default_value( 2 ), "order of spatial discretisation (velocity)" )
    ( "order-p", Feel::po::value<int>()->default_value( 1 ), "order of spatial discretisation (pressure)" )
    ( "gamma-bc", Feel::po::value<double>()->default_value( 20 ), "penalisation parameter" )
    ( "gamma-u", Feel::po::value<double>()->default_value( 10 ), "stabilisation parameter for velocity" )
    ( "gamma-p", Feel::po::value<double>()->default_value( 10 ), "stabilisation parameter for velocity" )
    ( "gamma-t", Feel::po::value<double>()->default_value( 0.01 ), "stabilisation parameter for heat convection" )
    ( "gamma-divdiv", Feel::po::value<double>()->default_value( 0.0 ), "stabilisation parameter for divergence jumps" )

    ( "delta-divdiv", Feel::po::value<double>()->default_value( 0.0 ), "divergence penalty term" )

    ( "eps-pseudo-compress", Feel::po::value<double>()->default_value( 0.0 ), "pseudo compressibility term (pressure coefficient)" )


    ( "linalg-same-prec", Feel::po::value<int>()->default_value( 1 ), "use same preconditioner" )
    ( "init", Feel::po::value<int>()->default_value( 1 ), "initialize Navier-Stokes solver (0=zero, 1=Stokes)" )

    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )

    ( "exportresults", Feel::po::value<int>()->default_value( 0 ), "export strategy (0=no export, 1=same mesh, 2=finest mesh)" )
    ( "export-profiles", "export 1D profiles" )

    ( "export-matlab", "export matrix and vectors in matlab" )
    ( "rb-nmax", Feel::po::value<int>()->default_value( 2 ), "Sampling dimension" )
    ( "rb-taille", Feel::po::value<int>()->default_value( 2 ), "taille" )
    ( "rb-epsilon1", Feel::po::value<double>()->default_value( 1e-3 ), "Estimator accuracy" )
    ( "rb-epsilon2", Feel::po::value<double>()->default_value( 1e-3 ), "Estimator accuracy" )
    ( "rb-tau", Feel::po::value<double>()->default_value( 1 ), "trust constant" )
    ( "mufile", Feel::po::value<std::string>()->default_value( "Echanti.txt" ),"file of parameters" )
    ( "pmu", "many parameters" )
    ( "rb", "reduced basis mode" )
    ( "adap", "construction adaptatif" )
    ( "orth", "with process Gramm Schmidt Orthonormalisation" )
    ;

    return opusoptions
           .add( makeComponentOptions() )
           .add( bdf_options( "temperature" ) );
}

OpusData::OpusData( int d )
    :
    M_is_steady( true ),
    M_dimension( d ),
    M_order_g( 1 ),
    M_order_time( 1 ),
    M_order_temp( 2 ),
    M_order_u( 2 ),
    M_order_p( 1 ),
    //M_h( 0.5 ),
    M_h( 1. ),

#if 0
    M_gamma_p( 0.25 ),
    M_gamma_u( 0.25 ),
    M_gamma_conv_T( 0.1 ),
#else
    M_gamma_bc( 10 ),
    M_gamma_p( 10 ),
    M_gamma_u( 10 ),
    M_gamma_conv_T( 0.01 ),
#endif
    M_gamma_divdiv( 0.0 ),
    M_delta_divdiv( 0.0 ),
    M_eps_compress( 0.0 ),

    M_use_same_prec( 1 ),
    M_init( INIT_FLOW_WITH_ZERO ),
    M_export( EXPORT_LAGP1_MESH ),

    M_components(),

    M_dirichlet_temp(),
    M_dirichlet_velocity(),
    M_dirichlet_pressure()


{
    M_components.insert( std::make_pair( "PCB", OpusComponent( "pcb", /*k*/0.2,    /*rhoC*/2*1e6,   /*Q*/0,   /*h*/13*1e-2, /*e*/2*1e-3 ) ) );
    M_components.insert( std::make_pair( "IC1", OpusComponent( "ic1", /*k*/10,     /*rhoC*/1.4*1e6, /*Q*/1e6, /*h*/2*1e-2,  /*e*/2*1e-3 ) ) );
    M_components.insert( std::make_pair( "IC2", OpusComponent( "ic2", /*k*/10,     /*rhoC*/1.4*1e6, /*Q*/1e6, /*h*/2*1e-2,  /*e*/2*1e-3 ) ) );
    //M_components.insert( std::make_pair("AIR", OpusComponent( "air", /*k*/3*1e-2, /*rhoC*/1100,     /*Q*/0,   /*h*/13*1e-2, /*e*/5*1e-2 )) );
    M_components.insert( std::make_pair( "AIR", OpusComponent( "air", /*k*/3*1e-2, /*rhoC*/1100,     /*Q*/0,   /*h*/13*1e-2, /*e*/4*1e-3, /*flow rate*/5e-3 ) ) );

    M_dirichlet_velocity.push_back( "Gamma_4_AIR1" );
    M_dirichlet_velocity.push_back( "Gamma_4_AIR4" );
    M_dirichlet_velocity.push_back( "Gamma_4_AIR" );
    M_dirichlet_velocity.push_back( "Gamma_1" );
    M_dirichlet_velocity.push_back( "Gamma_2" );

    print();
}

OpusData::OpusData( OpusData const& opusdata )
    :
    M_is_steady( opusdata.M_is_steady ),
    M_dimension( opusdata.M_dimension ),
    M_order_g( opusdata.M_order_g ),
    M_order_time( opusdata.M_order_time ),
    M_order_temp( opusdata.M_order_temp ),
    M_order_u( opusdata.M_order_u ),
    M_order_p( opusdata.M_order_p ),
    M_h( opusdata.M_h ),

    M_gamma_p( opusdata.M_gamma_p ),
    M_gamma_u( opusdata.M_gamma_u ),
    M_gamma_conv_T( opusdata.M_gamma_conv_T ),

    M_gamma_divdiv( opusdata.M_gamma_divdiv ),

    M_delta_divdiv( opusdata.M_delta_divdiv ),

    M_eps_compress( 0.0 ),

    M_use_same_prec( 1 ),

    M_init( INIT_FLOW_WITH_ZERO ),

    M_export( EXPORT_LAGP1_MESH ),

    M_components( opusdata.M_components ),
    M_dirichlet_temp(),
    M_dirichlet_velocity(),
    M_dirichlet_pressure()

{
    M_dirichlet_velocity.push_back( "Gamma_4_AIR1" );
    M_dirichlet_velocity.push_back( "Gamma_4_AIR4" );
    M_dirichlet_velocity.push_back( "Gamma_4_AIR" );
    M_dirichlet_velocity.push_back( "Gamma_1" );
    M_dirichlet_velocity.push_back( "Gamma_2" );

    print();
}
OpusData::OpusData( int d, Feel::po::variables_map const& vm )
    :
    M_vm ( vm ),

    M_components(),
    M_dirichlet_velocity(),
    M_dirichlet_pressure()

{

    M_is_steady = vm["steady"].as<bool>();

    M_dimension = d;

    M_order_g = vm["order-geo"].as<int>();
    M_order_time = vm["order-time"].as<int>();
    M_order_temp = vm["order-temp"].as<int>();
    M_order_u = vm["order-u"].as<int>();
    M_order_p = vm["order-u"].as<int>();

    M_h = vm["hsize"].as<double>();

    M_dirichlet_velocity.push_back( "Gamma_4_AIR1" );
    M_dirichlet_velocity.push_back( "Gamma_4_AIR4" );
    M_dirichlet_velocity.push_back( "Gamma_4_AIR" );
    M_dirichlet_velocity.push_back( "Gamma_1" );
    M_dirichlet_velocity.push_back( "Gamma_2" );
    M_dirichlet_pressure.push_back( "Gamma_3" );

    M_gamma_bc = vm["gamma-bc"].as<double>();
    M_gamma_p = vm["gamma-p"].as<double>();
    M_gamma_u = vm["gamma-u"].as<double>();
    M_gamma_conv_T = vm["gamma-t"].as<double>();
    M_gamma_divdiv = vm["gamma-divdiv"].as<double>();
    M_delta_divdiv = vm["delta-divdiv"].as<double>();
    M_eps_compress = vm["eps-pseudo-compress"].as<double>();

    M_use_same_prec = vm["linalg-same-prec"].as<int>();
    M_init = vm["init"].as<int>();
    M_export = vm["exportresults"].as<int>();

    M_components.insert( std::make_pair( "PCB", OpusComponent( "pcb", vm ) ) );
    M_components.insert( std::make_pair( "IC1", OpusComponent( "ic1", vm ) ) );
    M_components.insert( std::make_pair( "IC2", OpusComponent( "ic2", vm ) ) );
    M_components.insert( std::make_pair( "AIR", OpusComponent( "air", vm ) ) );

    print();


}
OpusData::~OpusData()
{}
void
OpusData::print() const
{

    LOG(INFO) << "Simulation parameters\n";
    LOG(INFO) << "=====================\n";
    LOG(INFO) << "h = " << this->h() << "\n";

    LOG(INFO) << "Stabilisation/Penalisation parameters\n";
    LOG(INFO) << "=====================================\n";
    LOG(INFO) << "gamma-bc     = " << this->gammaBc() << "\n";
    LOG(INFO) << "gamma-u      = " << this->gammaU() << "\n";
    LOG(INFO) << "gamma-t      = " << this->gammaTemp() << "\n";
    LOG(INFO) << "gamma-p      = " << this->gammaP() << "\n";
    LOG(INFO) << "gamma-divdiv = " << this->gammaDivDiv() << "\n";
    LOG(INFO) << "delta-divdiv = " << this->deltaDivDiv() << "\n";

    LOG(INFO) << "eps-peusdo-compress = " << this->epsPseudoCompressibility() << "\n";

    LOG(INFO) << "linalg-same-prec = " << this->useSamePreconditioner() << "\n";

    LOG(INFO) << "init = " << this->init() << "\n";

    LOG(INFO) << "export = " << this->doExport() << "\n";

    std::ostringstream ostr;
    ostr << this->component( "PCB" ) << "\n"
         << this->component( "IC1" ) << "\n"
         << this->component( "IC2" ) << "\n"
         << this->component( "AIR" ) << "\n";
    LOG(INFO) << ostr.str();
}

void
OpusData::load( const std::string &filename )
{
#if 0
    // Create an empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    // Load the XML file into the property tree. If reading fails
    // (cannot open file, parse error), an exception is thrown.
    read_xml( filename, pt );

    M_is_steady = pt.get<bool>( "opus.eads.fem.discretisation.time.steady", true );

    M_order_time = pt.get<int>( "opus.eads.fem.discretisation.time.order", 2 );
    M_order_temp = pt.get<int>( "opus.eads.fem.discretisation.temperature.order", 2 );
    M_order_u = pt.get<int>( "opus.eads.fem.discretisation.velocity.order", 2 );
    M_order_p = pt.get<int>( "opus.eads.fem.discretisation.pressure.order", 2 );

    M_h = pt.get<double>( "opus.eads.fem.mesh.h", 1 );
    M_order_g = pt.get<int>( "opus.eads.fem.mesh.order", 1 );

    M_gamma_bc = pt.get<double>( "opus.eads.fem.discretisation.dirichlet.gamma.bc", 10 );
    M_gamma_p = pt.get<double>( "opus.eads.fem.discretisation.dirichlet.gamma.p", 10 );
    M_gamma_u = pt.get<double>( "opus.eads.fem.discretisation.dirichlet.gamma.u", 10 );
    //M_gamma_divdiv = vm["gamma-divdiv"].as<double>();
    //M_delta_divdiv = vm["delta-divdiv"].as<double>();
    //M_eps_compress = vm["eps-pseudo-compress"].as<double>();

    //M_use_same_prec = vm["linalg-same-prec"].as<int>();
    M_init = pt.get<double>( "opus.eads.fem.discretisation.time.init", 0 );
    M_export = pt.get<double>( "opus.eads.fem.export", 0 );

    M_components["PCB"].load( pt );
    M_components["IC1"].load( pt );
    M_components["IC2"].load( pt );
    M_components["AIR"].load( pt );
#endif
}

void
OpusData::save( const std::string &filename )
{
#if 0
    // Create an empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    pt.put( "opus.eads.fem.discretisation.time.steady", M_is_steady );
    pt.put( "opus.eads.fem.discretisation.time.order", M_order_time );
    pt.put( "opus.eads.fem.discretisation.temperature.order", M_order_temp );
    pt.put( "opus.eads.fem.discretisation.velocity.order", M_order_u );
    pt.put( "opus.eads.fem.discretisation.pressure.order", M_order_p );

    pt.put( "opus.eads.fem.mesh.h", M_h );
    pt.put( "opus.eads.fem.mesh.order", M_order_g );

    pt.put( "opus.eads.fem.discretisation.dirichlet.gamma.bc", M_gamma_bc );
    pt.put( "opus.eads.fem.discretisation.dirichlet.gamma.p", M_gamma_p );
    pt.put( "opus.eads.fem.discretisation.dirichlet.gamma.u", M_gamma_u );

    pt.put( "opus.eads.fem.discretisation.time.init", M_init );
    pt.put( "opus.eads.fem.export", M_export );

    M_components["PCB"].save( pt );
    M_components["IC1"].save( pt );
    M_components["IC2"].save( pt );
    M_components["AIR"].save( pt );


    // Write the property tree to the XML file.
    write_xml( filename, pt );
#endif
}

void
OpusData::createMatlabScript()
{
    std::ofstream M_matlab( "opusdata.m" );

    M_matlab <<"function opusdata(Tspan) \n"
             <<"% % \n"
             <<"% % Results :  \n"
             <<"% %     * On figure 1, 4 graphics appear :\n"
             <<"% %         - the first with the drag coefficient Cd ; \n"
             <<"% %         - the second with the lift coeficient Cl ; \n"
             <<"% %         - the third with Delta P ; \n"
             <<"% %         - the last one with the constant 'Strouhal Number'. \n"
             <<"% %     * On figure 2 deux graphics appear :\n"
             <<"% %         - the first with frequencies used to calculate the Strouhal Number ; \n"
             <<"% %         - the other with the last coblumn of file 'opusdata.txt'.\n"
             <<"% %     * All graphics have the same representation : \n"
             <<"% %         - in blue the theorical results of the Benchmark ; \n"
             <<"% %         - in red the experimental results.\n"
             <<"% % \n"
             <<"% % This script is quite useful to see the results for the case test2D-2.  \n"
             <<"\n"
             <<"% Initialisation \n"
             <<"    Matrix=load(['opusdata.txt']); \n"
             <<"    [L,C]=size(Matrix); \n"
             <<"\n"
             <<"    % OpusData : \n"
             <<"disp('opusdata :');\n"
             <<"%    [Tzero,Tfinal,dt,time_order]=bdf();\n"
             //<<"    Re=" << M_Re << "; \n"
             //<<"    h=" << M_h << "; \n"
             //<<"    hcyl_scale=" << M_hcyl_scale << "; \n"
             <<"    dimension=" << M_dimension << "; \n"
             <<"    dt=" << M_vm["bdf.time-step"].as<double>() << "; \n"
             <<"IndexL=Matrix(:,1)>=Tspan(1) & Matrix(:,1)<= Tspan(2);\n"
             <<"% Benchmarks limits \n"
             <<"    % Drag coefficient :  \n"
             <<"        minCD=-3.22*ones(size(Matrix(IndexL,1),1),1); \n"
             <<"        maxCD=-3.24*ones(size(Matrix(IndexL,1)),1); \n"
             <<"    % Lift coefficient \n"
             <<"        minCL=0.99*ones(size(Matrix(IndexL,1),1),1); \n"
             <<"        maxCL=1.01*ones(size(Matrix(IndexL,1),1),1); \n"
             <<"    % Delta P \n"
             <<"        minDP=2.46*ones(size(Matrix(IndexL,1),1),1); \n"
             <<"        maxDP=2.5*ones(size(Matrix(IndexL,1),1),1); \n"
             <<"    % Strouhal Number \n"
             <<"        Minstrou=0.295*ones(size(Matrix(IndexL,1),1),1); \n"
             <<"        Maxstrou=0.305*ones(size(Matrix(IndexL,1),1),1); \n"
             <<"\n"
             <<"% Plotting graphs \n"
             <<"\n"
             //<<"            title(['Simplex_" << M_dimension << "_" << M_order_g << "/P2P1/h_" <<M_h << "_" << M_hcyl_scale << "/Re_" << M_Re << "']); \n"
             <<"    % First value \n"
             <<"        % Drag coefficient \n"
             <<"            subplot(3,2,1) \n"
             <<"            plot(Matrix(IndexL,1),Matrix(IndexL,2),'r',Matrix(IndexL,1),minCD,'b',Matrix(IndexL,1),maxCD,'b'); \n"
             <<"            ylabel('C_D'); \n"
             <<"            xlabel('time (s)'); \n"
             <<"            disp('Cd plotted'); \n"
             <<"\n"
             <<"    % Second value \n"
             <<"        % Lift coefficient \n"
             <<"            subplot(3,2,2) \n"
             <<"            plot(Matrix(IndexL,1),Matrix(IndexL,3),'r',Matrix(IndexL,1),minCL,'b',Matrix(IndexL,1),maxCL,'b'); \n"
             <<"            ylabel('C_L'); \n"
             <<"            xlabel('time (s)'); \n"
             <<"            disp('Cl plotted'); \n"
             <<"    % 3rd value \n"
             <<"        % Delta P \n"
             <<"            subplot(3,2,3) \n"
             <<"            plot(Matrix(IndexL,1),Matrix(IndexL,4),'r',Matrix(IndexL,1),minDP,'b',Matrix(IndexL,1),maxDP,'b'); \n"
             <<"            ylabel('\\Delta P'); \n"
             <<"            xlabel('time (s)'); \n"
             <<"            disp('DP plotted'); \n"
             <<"    % Last value \n"
             <<"        subplot(3,2,4); \n"
             <<"        plot(Matrix(IndexL,1),Matrix(IndexL,5),'r'); \n"
             <<"        ylabel('\\|u\\|_\\infty'); \n"
             <<"        xlabel('time (s)'); \n"
             <<"\n"
             <<"    % 5th value \n"
             <<"        % Frequency: f  \n"
             <<"            F=fft(Matrix(IndexL,3)); \n"
             <<"            N=length(F); \n"
             <<"            P=F(1:N/2).*conj(F(1:N/2))/(N/2); \n"
             <<"            freq=(1:N/2)/(N/2)*0.5/dt; \n"
             <<"        % Fréquencies plotting \n"
             <<"            subplot(3,2,5); \n"
             <<"            plot(freq,P,'-ro'); \n"
             <<"            ylabel('f'); \n"
             <<"            xlabel('N)'); \n"
             <<"            disp('Frequencies plotted'); \n"
             <<"\n"
             <<"    % Calculating Strouhal Number \n"
             <<"        [maxf,imaxf]=max(P); \n"
             <<"        %max_CL_frequency=freq(imaxf); \n"
             <<"        %DT=1/freq(imaxf); \n"
             <<"        strouhal=0.1*freq(imaxf)/(2*1.5/3); \n"
             <<"        Ystrou=strouhal*ones(size(Matrix(IndexL,1),1),1); \n"
             <<"    % Plotting Strouhal Number \n"
             <<"        subplot(3,2,6) \n"
             <<"        plot(Matrix(IndexL,1),Ystrou,'r',Matrix(IndexL,1),Minstrou,'b',Matrix(IndexL,1),Maxstrou,'b'); \n"
             <<"        xlabel('time (s)'); \n"
             <<"        ylabel('S_t'); \n"
             <<"        disp('Strouhal Number plotted'); \n"
             <<"\n"
             <<"    hold off; \n"
             <<"% end \n" ;

}

gmsh_ptrtype
OpusData::createMesh( double h, bool ref )
{
    gmsh_ptrtype gmshp( new Gmsh );

    //h = this->h();
    std::cout << "createMesh h=" << h << "\n";

    //std::ofstream costr( "constants.geo" );
    std::ostringstream costr;
    costr << gmshp->preamble()<<"\n";

    costr << "m=1;mm=10^-3;"
          << "\n"
          << "//\n"
          << "// Integrated circuit : IC\n"
          << "//\n"
          << "// thickness\n"
          << "e_IC  = " << this->component( "IC1" ).e() << "*m;\n"
          << "// length\n"
          << "L_IC  = " << this->component( "IC1" ).h() << "*m;\n"
          << "// position of the first IC\n"
          << "h_1   = 20*mm;\n"
          << "// position of the second IC\n"
          << "h_2   = 70*mm;\n"
          << "\n"
          << "\n"
          << "//\n"
          << "// PCB\n"
          << "//\n"
          << "// thickness\n"
          << "e_PCB = " << this->component( "PCB" ).e() << "*m;\n"
          << "// height\n"
          << "h_PCB = " << this->component( "PCB" ).h() << "*m;\n"
          << "\n"
          << "\n"
          << "//\n"
          << "// Air\n"
          << "//\n"
          << "// thickness\n";

    if ( ref )
        costr << "e_A = " << 5e-2 << "*m;\n";

    else
        costr << "e_A = " << this->component( "AIR" ).e() << "*m;\n";

    //<< "e_A = " << 5e-2 << "*m;\n"
    //<< "e_A = " << 2e-2 << "*m;\n"
    //<< "e_A = " << 4e-3 << "*m;\n"
    costr << "h=" << h << "*mm - 1e-8;\n"
          //<< "h=" << 1 << "*mm - 1e-8;\n"
          << "// Surface\n"
          << "Include \"geometry_heat.geo\";\n";


    std::ostringstream nameStr;
    std::ofstream ostr( "geometry_heat.geo" );
    ostr << gmshp->preamble()<<"\n";
    ostr
            << "/**\n"
            << " * Geometry for the Opus model\n"
            << " */\n"
            << "hs = "<<h<<";\n"
            << "p1=newp;Point(p1) = {0,0,0,hs};\n"
            << "p2=newp;Point(p2) = {e_PCB,0,0,hs};\n"
            << "\n"
            << "p3=newp;Point(p3) = {e_PCB,h_1,0,hs};\n"
            << "p4=newp;Point(p4) = {e_PCB,h_1+L_IC,0,hs};\n"
            << "\n"
            << "p5=newp;Point(p5) = {e_PCB,h_2,0,h};\n"
            << "p6=newp;Point(p6) = {e_PCB,h_2+L_IC,0,hs};\n"
            << "\n"
            << "p7=newp;Point(p7) = {e_PCB,h_PCB,0,hs};\n"
            << "p8=newp;Point(p8) = {0,h_PCB,0,hs};\n"
            << "\n"
            << "air_p1=p2;\n"
            << "air_p21=newp;Point(air_p21) = {e_PCB+e_IC,0,0,hs};\n"
            << "air_p22=newp;Point(air_p22) = {e_PCB+e_A,0,0,hs};\n"
            << "air_p31=newp;Point(air_p31) = {e_PCB+e_A,h_PCB,0,hs};\n"
            << "air_p32=newp;Point(air_p32) = {e_PCB+e_IC,h_PCB,0,hs};\n"
            << "air_p4=p7;\n"
            << "\n"
            << "air_p5=p3;\n"
            << "air_p51=newp;Point(air_p51) = {e_PCB+e_IC,h_1,0,hs};\n"
            << "air_p6=p4;\n"
            << "air_p61=newp;Point(air_p61) = {e_PCB+e_IC,h_1+L_IC,0,hs};\n"
            << "\n"
            << "air_p7=p5;\n"
            << "air_p71=newp;Point(air_p71) = {e_PCB+e_IC,h_2,0,hs};\n"
            << "air_p8=p6;\n"
            << "air_p81=newp;Point(air_p81) = {e_PCB+e_IC,h_2+L_IC,0,hs};\n"
            << "\n"
            << "ic1_p1=p3;\n"
            << "ic1_p2=air_p51;\n"
            << "ic1_p3=air_p61;\n"
            << "ic1_p4=p4;\n"
            << "\n"
            << "ic2_p1=p5;\n"
            << "ic2_p2=air_p71;\n"
            << "ic2_p3=air_p81;\n"
            << "ic2_p4=p6;\n"
            << "\n"
            << "Line(1) = {1, 2};\n"
            << "Line(2) = {2, 9};\n"
            << "Line(3) = {9, 10};\n"
            << "Line(4) = {10, 11};\n"
            << "Line(5) = {11, 12};\n"
            << "Line(6) = {12, 7};\n"
            << "Line(7) = {7, 8};\n"
            << "Line(8) = {7, 8};\n"
            << "Line(9) = {8, 1};\n"
            << "Line(10) = {2, 3};\n"
            << "Line(11) = {3, 4};\n"
            << "Line(12) = {4, 5};\n"
            << "Line(13) = {5, 6};\n"
            << "Line(14) = {6, 7};\n"
            << "Line(15) = {9, 13};\n"
            << "Line(16) = {13, 14};\n"
            << "Line(17) = {14, 15};\n"
            << "Line(18) = {15, 16};\n"
            << "Line(19) = {16, 12};\n"
            << "Line(20) = {3, 13};\n"
            << "Line(21) = {4, 14};\n"
            << "Line(22) = {5, 15};\n"
            << "Line(23) = {6, 16};\n"
            << "Line Loop(24) = {4, 5, -19, -18, -17, -16, -15, 3};\n"
            << "Plane Surface(25) = {24};\n"
            << "Line Loop(26) = {19, 6, -14, 23};\n"
            << "Plane Surface(27) = {26};\n"
            << "Line Loop(28) = {18, -23, -13, 22};\n"
            << "Plane Surface(29) = {28};\n"
            << "Line Loop(30) = {17, -22, -12, 21};\n"
            << "Plane Surface(31) = {30};\n"
            << "Line Loop(32) = {16, -21, -11, 20};\n"
            << "Plane Surface(33) = {32};\n"
            << "Line Loop(34) = {10, 20, -15, -2};\n"
            << "Plane Surface(35) = {34};\n"
            << "Line Loop(36) = {9, 1, 10, 11, 12, 13, 14, 7};\n"
            << "Plane Surface(37) = {36};\n"
            << "\n"
            << "//Physical Line(38) = {10, 11, 12, 13, 14};\n"
            << "Physical Line(\"Gamma_1\") = {9};\n"
            << "Physical Line(\"Gamma_2\") = {4};\n"
            << "Physical Line(\"Gamma_IC1_PCB\") = {11};\n"
            << "Physical Line(\"Gamma_IC2_PCB\") = {13};\n"
            << "\n"
            << "Physical Line(\"Gamma_4_AIR1\") = {2};\n"
            << "Physical Line(\"Gamma_4_AIR4\") = {3};\n"
            << "Physical Line(\"Gamma_3_AIR4\") = {5};\n"
            << "Physical Line(\"Gamma_3_AIR3\") = {6};\n"
            << "Physical Line(\"Gamma_3_PCB\") = {7};\n"
            << "Physical Line(\"Gamma_4_PCB\") = {1};\n"
            << "Physical Surface(\"PCB\") = {37};\n"
            << "Physical Surface(\"AIR123\") = {35, 31, 27};\n"
            << "Physical Surface(\"IC1\") = {33};\n"
            << "Physical Surface(\"IC2\") = {29};\n"
            << "Physical Surface(\"AIR4\") = {25};\n";

    nameStr << "opusthermal";
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

gmsh_ptrtype
OpusData::createMeshLine( double h )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;

    gmsh_ptrtype gmshp( new Gmsh );

    std::ofstream costr( "constantsline.geo" );
    costr<< gmshp->preamble()<<"\n";
    costr << "m=1;mm=10^-3;"
          << "\n"
          << "//\n"
          << "// Integrated circuit : IC\n"
          << "//\n"
          << "// thickness\n"
          << "e_IC  = " << this->component( "IC1" ).e() << "*m;\n"
          << "// length\n"
          << "L_IC  = " << this->component( "IC1" ).h() << "*m;\n"
          << "// position of the first IC\n"
          << "h_1   = 20*mm;\n"
          << "// position of the second IC\n"
          << "h_2   = 70*mm;\n"
          << "\n"
          << "\n"
          << "//\n"
          << "// PCB\n"
          << "//\n"
          << "// thickness\n"
          << "e_PCB = " << this->component( "PCB" ).e() << "*m;\n"
          << "// height\n"
          << "h_PCB = " << this->component( "PCB" ).h() << "*m;\n"
          << "\n"
          << "\n"
          << "//\n"
          << "// Air\n"
          << "//\n"
          << "// thickness\n"
          << "e_A = " << this->component( "AIR" ).e() << "*m;\n";

    ostr << gmshp->preamble()<<"\n";
    ostr << "Include \"constantsline.geo\";"
         << "\n"
         << "h=" << h << "*mm-1e-8;\n"
         << "p1=newp;Point(p1) = {e_PCB,0,0,h};\n"
         << "p2=newp;Point(p2) = {e_PCB,h_PCB,0,h};\n"
         << "Line(1) = {p1, p2};\n"
         << "\n"
         << "\n"
         << "Physical Line(\"Profile\") = {1};\n";


    nameStr << "opusline";

    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;

}

gmsh_ptrtype
OpusData::createMeshCrossSection2( double h )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;

    gmsh_ptrtype gmshp( new Gmsh );

    std::ofstream costr( "constantscs2.geo" );
    costr << gmshp->preamble()<<"\n";
    costr << "m=1;mm=10^-3;"
          << "\n"
          << "//\n"
          << "// Integrated circuit : IC\n"
          << "//\n"
          << "// thickness\n"
          << "e_IC  = " << this->component( "IC1" ).e() << "*m;\n"
          << "// length\n"
          << "L_IC  = " << this->component( "IC1" ).h() << "*m;\n"
          << "// position of the first IC\n"
          << "h_1   = 20*mm;\n"
          << "// position of the second IC\n"
          << "h_2   = 70*mm;\n"
          << "\n"
          << "\n"
          << "//\n"
          << "// PCB\n"
          << "//\n"
          << "// thickness\n"
          << "e_PCB = " << this->component( "PCB" ).e() << "*m;\n"
          << "// height\n"
          << "h_PCB = " << this->component( "PCB" ).h() << "*m;\n"
          << "\n"
          << "\n"
          << "//\n"
          << "// Air\n"
          << "//\n"
          << "// thickness\n"
          << "e_A = " << this->component( "AIR" ).e() << "*m;\n";

    ostr << "Include \"constantscs2.geo\";"
         << "\n"
         << "h=" << 0.05 << "*mm-1e-8;\n"
         << "p1=newp;Point(p1) = {0,0.08,0,h};\n"
         << "p2=newp;Point(p2) = {e_PCB+e_A,0.08,0,h};\n"
         << "Line(1) = {p1, p2};\n"
         << "\n"
         << "\n"
         << "Physical Line(\"Profile\") = {1};\n";


    nameStr << "opus-cross-section-2";

    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}

gmsh_ptrtype
OpusData::createMeshAir( double h )
{
    gmsh_ptrtype gmshp( new Gmsh );

    std::ostringstream costr;
    costr << gmshp->preamble()<<"\n";
    costr << "m=1;mm=10^-3;"
          << "\n"


          << "//\n"
          << "// Integrated circuit : IC\n"
          << "//\n"
          << "// thickness\n"
          << "e_IC  = " << this->component( "IC1" ).e() << "*m;\n"
          << "// length\n"
          << "L_IC  = " << this->component( "IC1" ).h() << "*m;\n"
          << "// position of the first IC\n"
          << "h_1   = 20*mm;\n"
          << "// position of the second IC\n"
          << "h_2   = 70*mm;\n"
          << "\n"
          << "\n"
          << "//\n"
          << "// PCB\n"
          << "//\n"
          << "// thickness\n"
          << "e_PCB = " << this->component( "PCB" ).e() << "*m;\n"
          << "// height\n"
          << "h_PCB = " << this->component( "PCB" ).h() << "*m;\n"
          << "\n"
          << "\n"
          << "//\n"
          << "// Air\n"
          << "//\n"
          << "// thickness\n"
          << "e_A = " << this->component( "AIR" ).e() << "*m;\n"
          << "\n"
          << "h=" << M_vm["fluid.hsize"].as<double>() << "*mm - 1e-8;\n"
          //<< "h=" << 1 << "*mm - 1e-8;\n"
          << "Include \"geometry_air.geo\";\n"
          << "Physical Surface(\"AIR\") = {25};\n\n";

    std::ofstream ostr( "geometry_air.geo" );;
    std::ostringstream nameStr;

    ostr << gmshp->preamble()<<"\n";
    ostr << "p1=newp;Point(p1) = {0,0,0,h};\n"
         << "p2=newp;Point(p2) = {e_PCB,0,0,h};\n"
         << "\n"
         << "p3=newp;Point(p3) = {e_PCB,h_1,0,h};\n"
         << "p4=newp;Point(p4) = {e_PCB,h_1+L_IC,0,h};\n"
         << "\n"
         << "p5=newp;Point(p5) = {e_PCB,h_2,0,h};\n"
         << "p6=newp;Point(p6) = {e_PCB,h_2+L_IC,0,h};\n"
         << "\n"
         << "p7=newp;Point(p7) = {e_PCB,h_PCB,0,h};\n"
         << "p8=newp;Point(p8) = {0,h_PCB,0,h};\n"
         << "\n"
         << "air_p1=p2;\n"
         << "air_p2=newp;Point(air_p2) = {e_PCB+e_A,0,0,h};\n"
         << "air_p3=newp;Point(air_p3) = {e_PCB+e_A,h_PCB,0,h};\n"
         << "air_p4=p7;\n"
         << "\n"
         << "air_p5=p3;\n"
         << "air_p51=newp;Point(air_p51) = {e_PCB+e_IC,h_1,0,h};\n"
         << "air_p6=p4;\n"
         << "air_p61=newp;Point(air_p61) = {e_PCB+e_IC,h_1+L_IC,0,h};\n"
         << "\n"
         << "air_p7=p5;\n"
         << "air_p71=newp;Point(air_p71) = {e_PCB+e_IC,h_2,0,h};\n"
         << "air_p8=p6;\n"
         << "air_p81=newp;Point(air_p81) = {e_PCB+e_IC,h_2+L_IC,0,h};\n"
         << "\n"
         << "ic1_p1=p3;\n"
         << "ic1_p2=air_p51;\n"
         << "ic1_p3=air_p61;\n"
         << "ic1_p4=p4;\n"
         << "\n"
         << "ic2_p1=p5;\n"
         << "ic2_p2=air_p71;\n"
         << "ic2_p3=air_p81;\n"
         << "ic2_p4=p6;\n"
         << "\n"
         << "Line(1) = {1, 2};\n"
         << "Line(2) = {2, 9};\n"
         << "Line(3) = {9, 10};\n"
         << "Line(4) = {10, 7};\n"
         << "Line(5) = {7, 6};\n"
         << "Line(6) = {6, 14};\n"
         << "Line(7) = {14, 13};\n"
         << "Line(8) = {13, 5};\n"
         << "Line(9) = {5, 6};\n"
         << "Line(10) = {5, 4};\n"
         << "Line(11) = {4, 12};\n"
         << "Line(12) = {12, 11};\n"
         << "Line(13) = {11, 3};\n"
         << "Line(14) = {3, 4};\n"
         << "Line(15) = {3, 2};\n"
         << "Line(16) = {7, 8};\n"
         << "Line(17) = {8, 1};\n"
         << "Line Loop(18) = {17, 1, -15, 14, -10, 9, -5, 16};\n"
         << "Plane Surface(19) = {18};\n"
         << "Line Loop(20) = {9, 6, 7, 8};\n"
         << "Plane Surface(21) = {20};\n"
         << "Line Loop(22) = {14, 11, 12, 13};\n"
         << "Plane Surface(23) = {22};\n"
         << "Line Loop(24) = {5, 6, 7, 8, 10, 11, 12, 13, 15, 2, 3, 4};\n"
         << "Plane Surface(25) = {24};\n"
         << "\n"
         << "\n"

         << "Physical Line(\"Gamma_2\") = {3};\n"
         << "Physical Line(\"Gamma_3_AIR\") = {4};\n"

         << "Physical Line(\"Gamma_1\") = {5, 6, 7, 8, 10, 11, 12, 13, 15};\n"
         << "Physical Line(\"Gamma_4_AIR\") = {2};\n"
         << "//Physical Surface(\"AIR\") = {25};\n\n";

    nameStr << "opusair";

    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}


OpusData::Inflow::Inflow( OpusData const& opusdata, double time )
    :
    M_opusdata( opusdata ),
    M_time( time )
{}


double
OpusData::Inflow::operator()( uint16_type c1, uint16_type c2, node_type const& p, node_type const& /*n*/ ) const
{
#if 0
    double time_term = 1.0;

    if ( M_opusdata.inflowType() == INFLOW_UNSTEADY )
        time_term = sin( M_PI*M_time/8 );

    if ( M_opusdata.d() == 2 )
    {
        switch ( c1 )
        {
        case 0:
            return 4*M_opusdata.Um()*p[1]*( M_opusdata.H()-p[1] )*time_term/pow( M_opusdata.H(),2 );

        case 1:
        default:
            return 0;
        }
    }

    else
    {

        // 3D
        switch ( c1 )
        {
        case 0:
            return 16*M_opusdata.Um()*p[1]*( M_opusdata.H()-p[1] )*p[2]*( M_opusdata.H()-p[2] )*time_term/pow( M_opusdata.H(),4 );

        case 1:
        case 2:
        default:
            return 0;
        }
    }

#endif // 0
}


}
