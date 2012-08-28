/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-03-04

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)

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
   \file convection_run.hpp
   \author Christophe Prud'homme <christophe.nprudhomme@ujf-grenoble.fr>
   \date 2009-03-04
 */
#include "convection.hpp"

// <int Order_s, int Order_p, int Order_t>
void
Convection::run()
{
    std::cout << "start run()\n";
    std::cout << "gr=" << this->vm()["gr"]. as<double>() << std::endl;
    std::cout << "pr=" << this->vm()["pr"]. as<double>() << std::endl;
    std::cout << "h=" << this->vm()["hsize"]. as<double>() << std::endl;
    std::cout << "parameter defined\n";


    using namespace Feel::vf;

    // All together timer
    timers["all"].first.restart();

    std::ofstream timings( "runtime.txt" );
    std::cout << "1\n";

    //
    // --- MESH ---
    //
    timers["mesh"].first.restart();
    mesh_ptrtype mesh( new mesh_type );

    if (this->vm()["readMesh"]. as<int>()){
        std::string repository = this->vm()["input_dir"]. as<std::string>() ;
        std::string file_mesh = this->vm()["mesh_name"]. as<std::string>() ;;
        std::string complete_name = repository + file_mesh;
        std::cout << "Meshes read in file : " << complete_name <<std::endl;
        
        mesh  =  loadGMSHMesh( _mesh=new mesh_type,
                              _filename=complete_name,
                              _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );        
    }
    else{
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=createMesh(),
                               _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    }

    Log() << "Tfixed: " << mesh->markerName( "Tfixed" ) << ": " << integrate(markedfaces(mesh,"Tfixed"), cst(1.) ).evaluate()(0,0) << "\n";
    Log() << "Tflux: " << mesh->markerName( "Tflux" ) << ": " << integrate(markedfaces(mesh,"Tflux"), cst(1.) ).evaluate()(0,0)  << "\n";
    //Log() << "Fflux: " << mesh->markerName( "Fflux" ) << ": " << integrate(markedfaces(mesh,"Fflux"), cst(1.) ).evaluate()(0,0)  << "\n";
    Log() << "Tinsulated: " << mesh->markerName( "Tinsulated" ) << ": " << integrate(markedfaces(mesh,"Tinsulated"), cst(1.) ).evaluate()(0,0)  << "\n";
    timers["mesh"].second=timers["mesh"].first.elapsed();
    timings << "[Mesh] Time : " << timers["mesh"].second << std::endl;
    //
    // --- END MESH SECTION ---


    //
    // --- FUNCTION SPACE ---
    //
    timers["fspace"].first.restart();
    // Espace des fonctions et elements
    Xh = space_type::New( mesh );

    element_type U( Xh, "u" );
    element_type Un( Xh, "un" );
    element_type V( Xh, "v" );
    element_type W( Xh, "v" );

    element_0_type u = U. element<0>(); // fonction vitesse
    element_0_type un = Un. element<0>(); // fonction vitesse
    element_0_type v = V. element<0>(); // fonction test vitesse
    element_1_type p = U. element<1>(); // fonction pression


    element_1_type pn = Un. element<1>(); // fonction pression
    element_1_type q = V. element<1>(); // fonction test pression
    element_2_type t = U. element<2>(); // fonction temperature
    element_2_type tn = Un. element<2>(); // fonction temperature
    element_2_type s = V. element<2>(); // fonction test temperature
#if defined( FEELPP_USE_LM )
    element_3_type xi = U. element<3>(); // fonction multipliers
    element_3_type eta = V. element<3>(); // fonction test multipliers
#endif

    Log() << "[convection::run()] u.size() = " << u.size() << " u.start() = " << u.start() << "\n";
    Log() << "[convection::run()] p.size() = " << p.size() << " p.start() = " << p.start() << "\n";
    Log() << "[convection::run()] t.size() = " << t.size() << " p.start() = " << t.start() << "\n";
#if defined( FEELPP_USE_LM )
    Log() << "[convection::run()] xi.size() = " << xi.size() << " p.start() = " << xi.start() << "\n";
#endif
    Log() << "[convection::run()] U.size() = " << U.size() << " Xh ndof = " << Xh->nDof() << "\n";

#if CONVECTION_DIM == 2
    u = vf::project( Xh-> functionSpace<0>(), elements( mesh ), vec( Px()*Py(),Py()*Px() ) );
    un = vf::project( Xh-> functionSpace<0>(), elements( mesh ), vec( Px()*Py(),Py()*Px() ) );
    std::cout << integrate( elements( mesh ), idv( u ) , _Q<2>() ).evaluate() << "\n";
    std::cout << integrate( boundaryfaces( mesh ), gradv( u )*N() , _Q<1>() ).evaluate() << "\n";
#endif

    p = vf::project( Xh->  functionSpace<1>(), elements( mesh ), exp( Px() ) );
    pn = vf::project( Xh->  functionSpace<1>(), elements( mesh ), exp( Px() ) );
    t = vf::project( Xh->  functionSpace<2>(), elements( mesh ), sin( Py() ) );
    tn = vf::project( Xh->  functionSpace<2>(), elements( mesh ), sin( Py() ) );


    std::cout << integrate( elements( mesh ), idv( p ) , _Q<3>() ).evaluate() << "\n";
    std::cout << integrate( elements( mesh ), idv( t ) , _Q<6>() ).evaluate() << "\n";


    std::cout << integrate( boundaryfaces( mesh ), gradv( p )*N() , _Q<3>() ).evaluate() << "\n";
    std::cout << integrate( boundaryfaces( mesh ), gradv( t )*N() , _Q<6>() ).evaluate() << "\n";

#if defined( FEELPP_USE_LM )
    xi = vf::project( Xh->  functionSpace<3>(), elements( mesh ), constant( 1.0 ) );
    std::cout << integrate( elements( mesh ), idv( xi ), _Q<1>() ).evaluate() << "\n";
#endif

    std::cout<< "----1----"<<std::endl;
    int adim=this->vm()["adim"]. as<int>();
    timers["fspace"].second=timers["fspace"].first.elapsed();
    timings << "[F spaces] Time : " << timers["fspace"].second << std::endl;
    std::cout<< "----2----"<<std::endl;
    // init to 0 and then later reuse previous grashof results to
    // initialize the solver

#if CONVECTION_DIM==2
    u = vf::project( Xh-> functionSpace<0>(), elements( mesh ), vec( constant( 0.0 ),constant( 0.0 ) ) );
    un = vf::project( Xh-> functionSpace<0>(), elements( mesh ), vec( constant( 0.0 ),constant( 0.0 ) ) );
#else
    u = vf::project( Xh-> functionSpace<0>(), elements( mesh ), vec( constant( 0.0 ),constant( 0.0 ), cst(0.) ) );
    un = vf::project( Xh-> functionSpace<0>(), elements( mesh ), vec( constant( 0.0 ),constant( 0.0 ), cst(0.) ) );
#endif
    p = vf::project( Xh->  functionSpace<1>(), elements( mesh ), constant( 0.0 ) );
    pn = vf::project( Xh->  functionSpace<1>(), elements( mesh ), constant( 0.0 ) );

    if ( adim==1 )
    {
        t = vf::project( Xh->  functionSpace<2>(), elements( mesh ), constant( 0.0 ) );
        tn = vf::project( Xh->  functionSpace<2>(), elements( mesh ), constant( 0.0 ) );
    }

    else
    {
        t = vf::project( Xh->  functionSpace<2>(), elements( mesh ), constant( 300 ) );
        tn = vf::project( Xh->  functionSpace<2>(), elements( mesh ), constant( 300 ) );
    }
#if defined( FEELPP_USE_LM )
    xi = vf::project( Xh->  functionSpace<3>(), elements( mesh ), constant( 0.0 ) );
#endif
    //initialisation timer :
    double titi = 0;
    boost::shared_ptr<export_type> exporter( export_type::New( this->vm(),
            ( boost::format( "%1%" )
              % this->about().appName() ).str() ) );

    // set up the non linear solver

    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual,
                                                   boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian,
                                                   boost::ref( *this ), _1, _2 );

    std::cout<< "----3----"<<std::endl;
    // Output for the benchmark data for each grashof number
    std::ofstream benchOut( "benchmark.dat" );



    std::cout<< "----4----"<<std::endl;

    vector_ptrtype R( M_backend->newVector( Xh ) );
    sparse_matrix_ptrtype J( M_backend->newMatrix( Xh,Xh ) );

    Log() << "============================================================\n";
    std::cout << "============================================================\n";
    double gr( this->vm()["gr"]. as<double>() );
    M_current_Grashofs = gr;
    double pr = this->vm()["pr"]. as<double>();
    M_current_Prandtl = pr;
    Log() << "Grashof = " << M_current_Grashofs << "\n";
    Log() << "Prandtl = " << M_current_Prandtl << "\n";
    std::cout << "Grashof = " << M_current_Grashofs << "\n";
    std::cout << "Prandtl = " << M_current_Prandtl << "\n";

    int N=std::max( 1.0,std::max( std::ceil( std::log( gr ) ),std::ceil( std::log( pr )-std::log( 1e-2 ) ) ) );

    for ( int i = 0; i < N; ++i )
    {
        int denom = ( N==1 )?1:N-1;
        M_current_Grashofs = math::exp( math::log( 1. )+i*( math::log( gr )-math::log( 1. ) )/denom );
        M_current_Prandtl = math::exp( math::log( 1e-2 )+i*( math::log( pr )-math::log( 1e-2 ) )/denom );
        std::cout << "i/N = " << i << "/" << N <<std::endl;
        std::cout << " intermediary Grashof = " << M_current_Grashofs<<std::endl;
        std::cout<< " and Prandtl = " << M_current_Prandtl << "\n"<<std::endl;
        M_backend->nlSolve( _solution = U );

        if ( exporter->doExport() )
        {
            Log() << "exportResults starts\n";

            exporter->step( i )->setMesh( mesh );

            exporter->step( i )->add( "u", u );
            exporter->step( i )->add( "t", t );
            exporter->step( i )->add( "p", p );

            exporter->save();
            Log() << "exportResults done\n";
        }

    }

    // value mean-pressure
    double meas = integrate( elements( mesh ),constant( 1.0 )  ).evaluate()( 0, 0 );
    std::cout << "measure(Omega)=" << meas << " (should be equal to 1)\n";
    std::cout << "mean pressure = "
              << integrate( elements( mesh ) ,idv( p ) ).evaluate()( 0,0 )/meas << "\n";

#if defined( FEELPP_USE_LM )
    Log() << "value of the Lagrange multiplier xi= " << xi( 0 ) << "\n";
    std::cout << "value of the Lagrange multiplier xi= " << xi( 0 ) << "\n";
#endif

    double mean_div_u = integrate( elements( mesh ),
                                   divv( u ) ).evaluate()( 0, 0 );
    std::cout << "mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = integrate( elements( mesh ),
                                       divv( u )*divv( u ) ).evaluate()( 0, 0 );
    std::cout << "||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";

    // calcul le nombre de Nusselt
    double AverageT = integrate( markedfaces( mesh,"Tflux" ) ,
                                 idv( t ) ).evaluate()( 0,0 ) ;
    std::cout << "AverageT = " << AverageT << std::endl;

#if 0
#if CONVECTION_DIM==2
    double Flux = integrate( markedfaces( mesh, "Fflux" ) ,
                             trans( idv( u ) )*vec( constant( -1.0 ),constant( 0.0 ) ) ).evaluate()( 0,0 ) ;
#else
    double Flux = integrate( markedfaces( mesh, "Fflux" ) ,
                            trans( idv( u ) )*vec( constant( -1.0 ),constant( 0.0 ), cst(0.0) ) ).evaluate()( 0,0 ) ;
#endif
    std::cout << "Flux = " << Flux << std::endl;
#endif
    // benchOut << M_current_Grashofs << " " << AverageT << " " << Flux << std::endl;

    //this->exportResults( boost::format("") , U, 0 );



    //benchOut.close();

    timers["all"].second=timers["all"].first.elapsed();
    timings << "[Run] Total execution time : " << timers["all"].second << std::endl;

    timings.close();

}

// instantiation
// class Convection<2,1,2>;
