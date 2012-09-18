/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2008-05-02

   Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file data.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-02
*/
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>

#include <feel/options.hpp>

#include <data.hpp>
#include <turek.hpp>


Data::data_ptrtype
Data::New( Feel::po::variables_map const& vm )
{
    using namespace Feel;

    if ( vm["d"].as<int>() == 2 )
    {
        if ( vm["order-u"].as<int>() == 2 )
        {
            if ( vm["order-geo"].as<int>() == 1 )
                return data_ptrtype( new Turek<2,2,1>( vm ) );

            //else
            //return data_ptrtype( new Turek<2,2,2>( vm ) );
        }

#if 1
    }

#else

        else if ( vm["order-u"].as<int>() == 3 )
        {
            if ( vm["order-geo"].as<int>() == 1 )
                //return data_ptrtype();//data_ptrtype( new Turek<2,3,1>( vm ) );
                return data_ptrtype( new Turek<2,3,1>( vm ) );

            else
                return data_ptrtype( new Turek<2,3,2>( vm ) );
        }

        else if ( vm["order-u"].as<int>() == 4 )
        {
            if ( vm["order-geo"].as<int>() == 1 )
                return data_ptrtype( new Turek<2,4,1>( vm ) );

            else
                return data_ptrtype( new Turek<2,4,2>( vm ) );
        }
    }

    else if ( vm["d"].as<int>() == 3 )
    {
        if ( vm["order-u"].as<int>() == 2 )
        {
            if ( vm["order-geo"].as<int>() == 1 )
                //return data_ptrtype();
                return data_ptrtype( new Turek<3,2,1>( vm ) );

            else
                //return data_ptrtype();
                return data_ptrtype( new Turek<3,2,2>( vm ) );
        }

        else if ( vm["order-u"].as<int>() == 3 )
        {
            if ( vm["order-geo"].as<int>() == 1 )
                return data_ptrtype( new Turek<3,3,1>( vm ) );

            else
                //return data_ptrtype();
                return data_ptrtype( new Turek<3,3,2>( vm ) );
        }

        else if ( vm["order-u"].as<int>() == 4 )
        {
            if ( vm["order-geo"].as<int>() == 1 )
                return data_ptrtype( new Turek<3,4,1>( vm ) );

            else
                return data_ptrtype( new Turek<3,4,2>( vm ) );
        }
    }

#endif
    throw std::logic_error( "invalid solver specifications" );
}

Feel::po::options_description
Data::makeOptions()
{
    Feel::po::options_description turekoptions( "Turek benchmark options" );
    turekoptions.add_options()
    ( "d", Feel::po::value<int>()->default_value( 2 ), "time step value" )

    ( "Re", Feel::po::value<double>()->default_value( 20 ), "Reynolds number" )
    ( "rho", Feel::po::value<double>()->default_value( 1 ), "density (kg/m^3)" )
    ( "nu", Feel::po::value<double>()->default_value( 1e-3 ), "dynamic viscosity (Pa s)" )
    ( "inflow-type", Feel::po::value<int>()->default_value( 0 ), "inflow type : 0=steady Poiseuille, 1=unsteady Poiseuille" )
    ( "cross-section-type", Feel::po::value<int>()->default_value( 0 ), "cross section type : 0=circular, 1=square" )

    ( "order-geo", Feel::po::value<int>()->default_value( 2 ), "order of geometry" )
    ( "order-u", Feel::po::value<int>()->default_value( 2 ), "order of spatial discretisation (velocity)" )
    ( "order-p", Feel::po::value<int>()->default_value( 1 ), "order of spatial discretisation (pressure)" )

    ( "gamma-bc", Feel::po::value<double>()->default_value( 100 ), "penalisation parameter" )
    ( "gamma-u", Feel::po::value<double>()->default_value( 10 ), "stabilisation parameter for velocity" )
    ( "gamma-p", Feel::po::value<double>()->default_value( 10 ), "stabilisation parameter for velocity" )

    ( "gamma-divdiv", Feel::po::value<double>()->default_value( 0.0 ), "stabilisation parameter for divergence jumps" )

    ( "delta-divdiv", Feel::po::value<double>()->default_value( 0.0 ), "divergence penalty term" )

    ( "eps-pseudo-compress", Feel::po::value<double>()->default_value( 0.0 ), "pseudo compressibility term (pressure coefficient)" )


    ( "linalg-same-prec", Feel::po::value<int>()->default_value( 1 ), "use same preconditioner" )
    ( "init", Feel::po::value<int>()->default_value( 1 ), "initialize Navier-Stokes solver (0=zero, 1=Stokes)" )


    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ( "h-cyl-scale", Feel::po::value<double>()->default_value( 5 ), "scale by which hsize is divided on the cylinder" )

    ( "export", Feel::po::value<int>()->default_value( 0 ), "export strategy (0=no export, 1=same mesh, 2=finest mesh)" )

    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return turekoptions.add( Feel::feel_options() );
}

Feel::AboutData
Data::makeAbout()
{
    Feel::AboutData about( "turek" ,
                           "turek" ,
                           "0.1",
                           "nD(n=2,3) Turek  benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}
Data::Data( int d )
    :
    M_dimension( d ),
    M_order_g( 1 ),
    M_order_u( 2 ),
    M_order_p( 1 ),
    M_h( 0.1 ),
    M_hcyl_scale( M_h/5 ),

    M_gamma_p( 0.25 ),
    M_gamma_u( 0.25 ),
    M_gamma_divdiv( 0.0 ),
    M_delta_divdiv( 0.0 ),
    M_eps_compress( 0.0 ),

    M_use_same_prec( 1 ),
    M_init( INIT_WITH_ZERO ),
    M_export( EXPORT_LAGP1_MESH ),

    M_H( 0.41 ),
    M_D( 0.1 ),
    M_Re( 20 ),
    M_nu( 1e-3 ),
    M_rho( 1e3 ),
    M_inflow_type( INFLOW_STEADY ),
    M_cross_section_type( CROSS_SECTION_CIRCULAR ),


    M_dirichlet_velocity(),
    M_dirichlet_pressure()

{

    M_dirichlet_velocity.push_back( "inflow" );
    M_dirichlet_velocity.push_back( "wall" );
    M_dirichlet_velocity.push_back( "cylinder" );

    M_dirichlet_pressure.push_back( "outflow" );
    print();
}

Data::Data( Data const& data )
    :
    M_dimension( data.M_dimension ),
    M_order_g( data.M_order_g ),
    M_order_u( data.M_order_u ),
    M_order_p( data.M_order_p ),
    M_h( data.M_h ),
    M_hcyl_scale( data.M_hcyl_scale ),




    M_gamma_p( data.M_gamma_p ),
    M_gamma_u( data.M_gamma_u ),

    M_gamma_divdiv( data.M_gamma_divdiv ),

    M_delta_divdiv( data.M_delta_divdiv ),

    M_eps_compress( 0.0 ),

    M_use_same_prec( 1 ),

    M_init( INIT_WITH_ZERO ),

    M_export( EXPORT_LAGP1_MESH ),

    M_H( data.M_H ),
    M_D( data.M_D ),
    M_Re( data.M_Re ),
    M_nu( data.M_nu ),
    M_rho( data.M_rho ),
    M_inflow_type( data.M_inflow_type ),
    M_cross_section_type( data.M_cross_section_type ),

    M_dirichlet_velocity(),
    M_dirichlet_pressure()

{
    M_dirichlet_velocity.push_back( "inflow" );
    M_dirichlet_velocity.push_back( "wall" );
    M_dirichlet_velocity.push_back( "cylinder" );

    M_dirichlet_pressure.push_back( "outflow" );

    print();
}
Data::Data( int d, Feel::po::variables_map const& vm )
    :
    M_vm ( vm ),
    M_dirichlet_velocity(),
    M_dirichlet_pressure()
{
    M_dimension = d;

    M_order_g = vm["order-geo"].as<int>();
    M_order_u = vm["order-u"].as<int>();
    M_order_p = vm["order-u"].as<int>();

    M_h = vm["hsize"].as<double>();
    M_hcyl_scale = vm["h-cyl-scale"].as<double>();

    M_H = 0.41;
    M_D = 0.1;

    M_Re = vm["Re"].as<double>();
    M_nu = vm["nu"].as<double>();
    M_rho = vm["rho"].as<double>();
    M_inflow_type = vm["inflow-type"].as<int>();
    M_cross_section_type = vm["cross-section-type"].as<int>();
    M_dirichlet_velocity.push_back( "inflow" );
    M_dirichlet_velocity.push_back( "wall" );
    M_dirichlet_velocity.push_back( "cylinder" );
    M_dirichlet_pressure.push_back( "outflow" );

    M_gamma_bc = vm["gamma-bc"].as<double>();
    M_gamma_p = vm["gamma-p"].as<double>();
    M_gamma_u = vm["gamma-u"].as<double>();
    M_gamma_divdiv = vm["gamma-divdiv"].as<double>();
    M_delta_divdiv = vm["delta-divdiv"].as<double>();
    M_eps_compress = vm["eps-pseudo-compress"].as<double>();

    M_use_same_prec = vm["linalg-same-prec"].as<int>();
    M_init = vm["init"].as<int>();
    M_export = vm["export"].as<int>();



    print();

}
Data::~Data()
{}
void
Data::print() const
{

    LOG(INFO) << "Simulation parameters\n";
    LOG(INFO) << "=====================\n";
    LOG(INFO) << "D = " << this->D() << "\n";
    LOG(INFO) << "H = " << this->H() << "\n";
    LOG(INFO) << "h = " << this->h() << "\n";
    LOG(INFO) << "h-cyl-scale = " << this->hCylinderScale() << "\n";
    LOG(INFO) << "Re = " << this->Re() << "\n";
    LOG(INFO) << "rho = " << this->rho() << "\n";
    LOG(INFO) << "nu = " << this->nu() << "\n";
    LOG(INFO) << "mu = " << this->mu() << "\n";
    LOG(INFO) << "inflow type = " << this->inflowType() << " (0: steady Poiseuille, 1: unsteady Poiseuille)\n";
    LOG(INFO) << "cross section type = " << this->crossSectionType() << " (0: circilar, 1: square)\n";

    LOG(INFO) << "Ubar = " << this->Ubar() << "\n";
    LOG(INFO) << "Um = " << this->Um() << "\n";

    LOG(INFO) << "Stabilisation/Penalisation parameters\n";
    LOG(INFO) << "=====================================\n";
    LOG(INFO) << "gamma-bc = " << this->gammaBc() << "\n";
    LOG(INFO) << "gamma-u = " << this->gammaU() << "\n";
    LOG(INFO) << "gamma-p = " << this->gammaP() << "\n";
    LOG(INFO) << "gamma-divdiv = " << this->gammaDivDiv() << "\n";
    LOG(INFO) << "delta-divdiv = " << this->deltaDivDiv() << "\n";

    LOG(INFO) << "eps-peusdo-compress = " << this->epsPseudoCompressibility() << "\n";

    LOG(INFO) << "linalg-same-prec = " << this->useSamePreconditioner() << "\n";

    LOG(INFO) << "init = " << this->init() << "\n";

    LOG(INFO) << "export = " << this->doExport() << "\n";



}
double
Data::scalingForce() const
{
    return  2./( this->rho()*Feel::math::pow( this->Ubar(),2 )*this->D()*Feel::math::pow( this->H(),M_dimension-2 ) );
}
Data::node_type
Data::xa() const
{
    switch ( M_dimension )
    {
    case 2:
    {
        node_type xa( 2 );
        xa[0]=0.15;
        xa[1]=0.2;
        return xa;
    }
    break;

    case 3:
    {
        node_type xa( 3 );
        xa[0]=0.45;
        xa[1]=0.20;
        xa[2]=0.205;
        return xa;
    }
    break;
    }

    return node_type();
}
Data::node_type
Data::xe() const
{
    switch ( M_dimension )
    {
    case 2:
    {
        node_type xe( 2 );
        xe[0]=0.25;
        xe[1]=0.2;
        return xe;
    }
    break;

    case 3:
    {
        node_type xe( 3 );
        xe[0]=0.55;
        xe[1]=0.20;
        xe[2]=0.205;
        return xe;
    }
    break;
    }

    return node_type();
}

void
Data::createMatlabScript()
{
    std::ofstream M_matlab( "data.m" );

    M_matlab <<"function data(Tspan) \n"
             <<"% % \n"
             <<"% % Results :  \n"
             <<"% %     * On figure 1, 4 graphics appear :\n"
             <<"% %         - the first with the drag coefficient Cd ; \n"
             <<"% %         - the second with the lift coeficient Cl ; \n"
             <<"% %         - the third with Delta P ; \n"
             <<"% %         - the last one with the constant 'Strouhal Number'. \n"
             <<"% %     * On figure 2 deux graphics appear :\n"
             <<"% %         - the first with frequencies used to calculate the Strouhal Number ; \n"
             <<"% %         - the other with the last column of file 'data.txt'.\n"
             <<"% %     * All graphics have the same representation : \n"
             <<"% %         - in blue the theorical results of the Benchmark ; \n"
             <<"% %         - in red the experimental results.\n"
             <<"% % \n"
             <<"% % This script is quite useful to see the results for the case test2D-2.  \n"
             <<"\n"
             <<"% Initialisation \n"
             <<"    Matrix=load(['data.txt']); \n"
             <<"    [L,C]=size(Matrix); \n"
             <<"\n"
             <<"    % Data : \n"
             <<"disp('data :');\n"
             <<"%    [Tzero,Tfinal,dt,time_order]=bdf();\n"
             <<"    Re=" << M_Re << "; \n"
             <<"    h=" << M_h << "; \n"
             <<"    hcyl_scale=" << M_hcyl_scale << "; \n"
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
             <<"            title(['Simplex_" << M_dimension << "_" << M_order_g << "/P2P1/h_" <<M_h << "_" << M_hcyl_scale << "/Re_" << M_Re << "']); \n"
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

std::pair<std::string,std::string>
Data::createCylinder()
{
    std::ostringstream ostr;
    std::ostringstream nameStr;

    switch ( M_dimension )
    {
    case 2:
    {
        ostr << "h=" << M_h << ";\n"
             << "\n"
             << "H=0.41;\n"
             << "D=0.1;\n"
             << "R=0.05;\n"
             << "L=2.2;\n"
             << "xc=0.2;\n"
             << "yc=0.2;\n"
             << "\n"
             << "Point(1) = {0,0,0,h};\n"
             << "Point(2) = {L,0,0,h};\n"
             << "Point(3) = {L,H,0,h};\n"
             << "Point(4) = {0,H,0,h};\n"
             << "\n"
             << "Line(1) = {1,2};\n"
             << "Line(2) = {2,3};\n"
             << "Line(3) = {3,4};\n"
             << "Line(4) = {4,1};\n"
             << "\n"
             << "Line Loop(20)={1,2,3,4};\n";

        if ( crossSectionType() == CROSS_SECTION_CIRCULAR )
        {
            ostr << "hcyl=h/" << M_hcyl_scale << ";\n"
                 << "Point(5) = {xc,yc,0,hcyl};\n"
                 << "Point(6) = {xc+R,yc,0,hcyl};\n"
                 << "Point(7) = {xc-R,yc,0,hcyl};\n"
                 << "Point(8) = {xc,yc+R,0,hcyl};\n"
                 << "Point(9) = {xc,yc-R,0,hcyl};\n"
                 << "\n"
                 << "Point( 11 ) = {0.1, 0.2, 0, hcyl};\n"
                 << "Point( 12 ) = {0.2, 0.3, 0, hcyl};\n"
                 << "Point( 13 ) = {0.2, 0.1, 0, hcyl};\n"
                 << "Point( 14 ) = {L, 0.1, 0, hcyl};\n"
                 << "Point( 15 ) = {L, 0.3, 0, hcyl};\n"
                 << "\n"
                 << "Circle(5) = {6,5,8};\n"
                 << "Circle(6) = {8,5,7};\n"
                 << "Circle(7) = {7,5,9};\n"
                 << "Circle(8) = {9,5,6};\n";
        }

        else // CROSS_SECTION_SQUARE
        {
            ostr << "// square : \n"
                 << "X=0.15; \n"
                 << "Y=0.15; \n"
                 << "Point(5) = {X,Y,0,hcyl}; \n"
                 << "Point(6) = {X,Y+D,0,hcyl}; \n"
                 << "Point(7) = {X+D,Y+D,0,hcyl}; \n"
                 << "Point(8) = {X+D,Y,0,hcyl}; \n"
                 << "Line(5) = {8,7}; \n"
                 << "Line(6) = {7,6}; \n"
                 << "Line(7) = {6,5}; \n"
                 << "Line(8) = {5,8}; \n";
        }

        ostr << "\n"
             << "\n"
             << "Line Loop(21) = {5,6,7,8};\n"
             << "\n"
             << "Plane Surface(30) = {20,-21};\n"
             << "\n"
             << "Physical Line(\"inflow\") = {4};\n"
             << "Physical Line(\"wall\") = {3,1};\n"
             << "Physical Line(\"cylinder\") = {5,6,7,8};\n"
             << "Physical Line(\"outflow\") = {2};\n"
             << "\n"
             << "Physical Surface(\"fluid\")={30};\n";

        nameStr << "cylinder";
    }
    break;

    case 3:
        if ( crossSectionType() == CROSS_SECTION_CIRCULAR )
        {
            ostr << "H=0.41;\n"
                 << "D=0.1;\n"
                 << "R=0.05;\n"
                 << "L=2.5;\n"
                 << "\n"
                 << "//h=0.015625;\n"
                 << "h= " << M_h << ";\n"
                 << "\n"
                 << "/*\n"
                 << " *\n"
                 << " */\n"
                 << "Point(1) = {0,0,0,h};\n"
                 << "Point(2) = {L,0,0,h};\n"
                 << "Point(3) = {L,H,0,h};\n"
                 << "Point(4) = {0,H,0,h};\n"
                 << "\n"
                 << "Line(1) = {1,2};\n"
                 << "Line(2) = {2,3};\n"
                 << "Line(3) = {3,4};\n"
                 << "Line(4) = {4,1};\n"
                 << "\n"
                 << "Line Loop(1)={1,2,3,4};\n"
                 << "\n"
                 << "\n"
                 << "/*\n"
                 << " * cylinder points\n"
                 << " */\n"
                 << "hcyl=h/4;\n"
                 << "Point(5) = {0.5,0.2,0,hcyl};\n"
                 << "Point(6) = {0.55,0.2,0,hcyl};\n"
                 << "Point(7) = {0.45,0.2,0,hcyl};\n"
                 << "Point(8) = {0.5,0.25,0,hcyl};\n"
                 << "Point(9) = {0.5,0.15,0,hcyl};\n"
                 << "Point(10) = {0.535355339059327,0.235355339059327,0,hcyl};\n"
                 << "Circle(5) = {7,5,10};\n"
                 << "Circle(6) = {10,5,9};\n"
                 << "Circle(7) = {9,5,7};\n"
                 << "Line Loop(8) = {5,6,7};\n"
                 << "Plane Surface(9) = {1,8};\n"
                 << "\n"
                 << "\n"
                 << "Extrude Surface {9, {0.0,0.0,H}};\n"
                 << "\n"
                 << "\n"
                 << "\n"
                 << "inlet=40;\n"
                 << "outlet=50;\n"
                 << "wall=60;\n"
                 << "cylinder=70;\n"
                 << "Physical Surface(\"inflow\") = {33};\n"
                 << "Physical Surface(\"wall\") = {29,9,21,46};\n"
                 << "Physical Surface(\"outflow\") = {25};\n"
                 << "Physical Surface(\"cylinder\") = {45,37,41};\n"
                 << "\n"
                 << "\n"
                 << "Physical Volume( 1 ) = { 1 };\n";
            nameStr << "cylinder";
        }

        else
        {
            ostr <<"//Variables :\n"
                 <<"H=0.41;\n"
                 <<"D=0.1;\n"
                 <<"L=2.5;\n"
                 <<"X=0.45;\n"
                 <<"Y=0.15;\n"
                 <<"h="<< M_h <<";\n"
                 <<"hcyl=h/"<< M_hcyl_scale <<";\n"
                 <<"// Boundary :\n"
                 <<"Point(1) = {0,0,0,h};\n"
                 <<"Point(2) = {0,H,0,h};\n"
                 <<"Point(3) = {L,H,0,h};\n"
                 <<"Point(4) = {L,0,0,h};\n"
                 <<"// Cylinder :\n"
                 <<"Point(5) = {X,Y,0,hcyl};\n"
                 <<"Point(6) = {X,Y+D,0,hcyl};\n"
                 <<"Point(7) = {X+D,Y+D,0,hcyl};\n"
                 <<"Point(8) = {X+D,Y,0,hcyl};\n"
                 <<"Line(1) = {8,7};\n"
                 <<"Line(2) = {7,6};\n"
                 <<"Line(3) = {6,5};\n"
                 <<"Line(4) = {5,8};\n"
                 <<"Line(5) = {1,4};\n"
                 <<"Line(6) = {4,3};\n"
                 <<"Line(7) = {3,2};\n"
                 <<"Line(8) = {2,1};\n"
                 <<"Line Loop(9) = {5,6,7,8};\n"
                 <<"Line Loop(11) = {4,1,2,3};\n"
                 <<"Plane Surface(12) = {-9,11};\n"
                 <<"Extrude {0,0,H} {\n"
                 <<"  Surface{12};\n"
                 <<"}\n"
                 <<"Surface Loop(55) = {54,25,12,29,33,37,41,45,49,53};\n"
                 <<"Volume(56) = {55};\n"
                 <<"Physical Surface(\"inflow\") = {25};\n"
                 <<"Physical Surface(\"outflow\") = {33};\n"
                 <<"Physical Surface(\"wall\") = {54,29,12,37};\n"
                 <<"Physical Surface(\"cylinder\") = {53,41,45,49};\n";

            nameStr << "cylinder";
        }

        break;

    default:
        std::ostringstream os;
        os << "invalid dimension: " << M_dimension;
        throw std::logic_error( os.str() );
    }

    return std::make_pair( nameStr.str(), ostr.str() );
}


Data::Inflow::Inflow( Data const& data, double time )
    :
    M_data( data ),
    M_time( time )
{}


double
Data::Inflow::operator()( uint16_type c1, uint16_type c2, node_type const& p, node_type const& /*n*/ ) const
{

    double time_term = 1.0;

    if ( M_data.inflowType() == INFLOW_UNSTEADY )
        time_term = sin( M_PI*M_time/8 );

    if ( M_data.d() == 2 )
    {
        switch ( c1 )
        {
        case 0:
            return 4*M_data.Um()*p[1]*( M_data.H()-p[1] )*time_term/pow( M_data.H(),2 );

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
            return 16*M_data.Um()*p[1]*( M_data.H()-p[1] )*p[2]*( M_data.H()-p[2] )*time_term/pow( M_data.H(),4 );

        case 1:
        case 2:
        default:
            return 0;
        }
    }
}

#if 0
void createOctave()
{
    std::ofstream ostr( "cylinder.m" );

    ostr << "function cylinder (Ti,Tf,dt, file,color)\n"
         << "  %%  Copyright (C) 2007 Université Joseph Fourier\n"
         << "  %% \n"
         << "  %%  This file is part of feel_cylinder.\n"
         << "  %% \n"
         << "  %%  feel_cylinder is free software; you can redistribute it and/or modify\n"
         << "  %%  it under the terms of the GNU General Public License as published by\n"
         << "  %%  the Free Software Foundation; either version 2, or (at your option)\n"
         << "  %%  any later version.\n"
         << "  %% \n"
         << "  %%  feel_cylinder is distributed in the hope that it will be useful, but\n"
         << "  %%  WITHOUT ANY WARRANTY; without even the implied warranty of\n"
         << "  %%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n"
         << "  %%  General Public License for more details.\n"
         << "  %% \n"
         << "  %%  You should have received a copy of the GNU General Public License\n"
         << "  %%  along with Octave; see the file COPYING.  If not, write to the Free\n"
         << "  %%  Software Foundation, 59 Temple Place - Suite 330, Boston, MA\n"
         << "  %%  02111-1307, USA.\N"
         << "  \n"
         << "  %% usage:  cylinder (dt)\n"
         << "  %%\n"
         << "  %%\n"
         << "\n"
         << "  %%  Author: Christophe Prud'homme <christophe.prudhomme@feelpp.org>\n"
         << "  %%  Keywords: CD, CL, DP and Frequency calculation\n"
         << "  %%  Maintainer: Christophe Prud'homme <christophe.prudhomme@feelpp.org>\n"
         << "\n"
         << "  \n"
         << "  if ( file == "" )\n"
         << "    file='quantities.dat';\n"
         << "  endif \n"
         << "  if ( color == "" )\n"
         << "    color='r*-';\n"
         << "  endif \n"
         << "  \n"
         << "  x=load(file);\n"
         << "  \n"
         << "  Nsteps=(Tf-Ti)/dt;\n"
         << "  NTi=Ti/dt;\n"
         << "  NT=Tf/dt;\n"
         << "  \n"
         << "  nrows=3;\n"
         << "  ncols=2;\n"
         << "  \n"
         << "  %%top_title(\"Cylinder Test Case Re=100, dt=0.0125\");\n"
         << "  \n"
         << "  %% CD\N"
         << "  %%minCD=5.57*ones(size(x(NTi:NT,1),1),1);\n"
         << "  %%maxCD=5.59*ones(size(x(NTi:NT,1),1),1);\n"
         << "  minCD=3.22*ones(size(x(NTi:NT,1),1),1);\n"
         << "  maxCD=3.24*ones(size(x(NTi:NT,1),1),1);\n"
         << "\n"
         << "  subplot(nrows,ncols,1);\n"
         << "  plot(x(NTi:NT,1),x(NTi:NT,2),[ color, ';CD(t);'],\n"
         << "       x(NTi:NT,1),minCD,'@11;min CD(t);',\n"
         << "       x(NTi:NT,1),maxCD,'@12;max CD(t);');\n"
         << "  disp(['CD plotted']);\n"
         << "  \n"
         << "  %% CL\N"
         << "  %minCL=0.0104*ones(size(x(NTi:NT,1),1),1);\n"
         << "  %maxCL=0.0110*ones(size(x(NTi:NT,1),1),1);\n"
         << "  minCL=0.99*ones(size(x(NTi:NT,1),1),1);\n"
         << "  maxCL=1.01*ones(size(x(NTi:NT,1),1),1);\n"
         << "\n"
         << "  subplot(nrows,ncols,2);\n"
         << "  plot(x(NTi:NT,1),x(NTi:NT,3), [color, ';CL(t);'],\n"
         << "       x(NTi:NT,1),minCL,'@11;min CL(t);',\n"
         << "       x(NTi:NT,1),maxCL,'@12;max CL(t);');\n"
         << "  \n"
         << "  disp(['CL plotted']);\n"
         << "  \n"
         << "  %% D P\n"
         << "  subplot(nrows,ncols,3);\n"
         << "  plot(x(NTi:NT,1),x(NTi:NT,5), [color, ';DP(t);']);\n"
         << "  \n"
         << "  disp(['DP plotted']);\n"
         << "  \n"
         << "  %%\n"
         << "  F=fft(x(NTi:NT,3));\n"
         << "  N=length(F);\n"
         << "  P=F(1:N/2).*conj(F(1:N/2))/(N/2);\n"
         << "  freq=(1:N/2)/(N/2)*0.5/dt;\n"
         << "  subplot(nrows,ncols,4);\n"
         << "  plot(freq,P,[color, ';P=f(CL(t));']);\n"
         << "  \n"
         << "  disp(['Frequencies plotted']);\n"
         << "  \n"
         << "  %% strouhal\n"
         << "  [maxf,imaxf]=max(P)\n"
         << "  max_CL_frequency=freq(imaxf)\n"
         << "  DT=1/freq(imaxf)\n"
         << "  strouhal=0.1*freq(imaxf)/(2*1.5/3)\n"
         << "  \n"
         << "  %%P=polyfit(\n"
         << "  \n"
         << "\n";
}

#endif
