/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-03-04

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-03-04
 */
#include "convection.hpp"

template <int Order_s, int Order_p, int Order_t>
void
Convection<Order_s, Order_p, Order_t>::run()
{
    std::cout << "start run()\n";
    std::cout << "gr=" << this->vm()["gr"].template as<double>() << std::endl;
    std::cout << "pr=" << this->vm()["pr"].template as<double>() << std::endl;
    std::cout << "h=" << this->vm()["hsize"].template as<double>() << std::endl;
    this->addParameterValue( Order_s )
        .addParameterValue( Order_p )
        .addParameterValue( Order_t )
        .addParameterValue( this->vm()["gr"].template as<double>() )
        .addParameterValue( this->vm()["pr"].template as<double>() )
        .addParameterValue( this->vm()["hsize"].template as<double>() );
    std::cout << "parameter defined\n";
    RunStatus ierr = this->preProcessing();
    if ( ierr == RUN_EXIT )
        return;

	using namespace Feel::vf;

	// All together timer
	timers["all"].first.restart();

    std::ofstream timings("runtime.txt");

    //
	// --- MESH ---
    //
	timers["mesh"].first.restart();
	mesh_ptrtype mesh( new mesh_type );
    mesh = createMesh();
	timers["mesh"].second=timers["mesh"].first.elapsed();
	timings << "[Mesh] Time : " << timers["mesh"].second << std::endl;
    //
	// --- END MESH SECTION ---
    //


    //
	// --- FUNCTION SPACE ---
    //
	timers["fspace"].first.restart();
	// Espace des fonctions et elements
    Xh = space_type::New( mesh );

    element_type U( Xh, "u" );
    element_type V( Xh, "v" );
	element_type W( Xh, "v" );
	element_0_type u = U.template element<0>(); // fonction vitesse
    element_0_type v = V.template element<0>(); // fonction test vitesse
	element_1_type p = U.template element<1>(); // fonction pression
	element_1_type q = V.template element<1>(); // fonction test pression
	element_2_type t = U.template element<2>(); // fonction temperature
	element_2_type s = V.template element<2>(); // fonction test temperature
    element_3_type xi = U.template element<3>(); // fonction multipliers
	element_3_type eta = V.template element<3>(); // fonction test multipliers

    Log() << "[convection::run()] u.size() = " << u.size() << " u.start() = " << u.start() << "\n";
    Log() << "[convection::run()] p.size() = " << p.size() << " p.start() = " << p.start() << "\n";
    Log() << "[convection::run()] t.size() = " << t.size() << " p.start() = " << t.start() << "\n";
    Log() << "[convection::run()] xi.size() = " << xi.size() << " p.start() = " << xi.start() << "\n";
    Log() << "[convection::run()] U.size() = " << U.size() << " Xh ndof = " << Xh->nDof() << "\n";

    u = vf::project( Xh->template functionSpace<0>(), elements(mesh), vec(Px()*Py(),Py()*Px()));
    p = vf::project( Xh->template  functionSpace<1>(), elements(mesh), exp(Px()) );
    t = vf::project( Xh->template  functionSpace<2>(), elements(mesh), sin(Py()) );
    xi = vf::project( Xh->template  functionSpace<3>(), elements(mesh), constant(1.0) );

    std::cout << integrate( elements(mesh), _Q<2>(), idv(u) ).evaluate() << "\n";
    std::cout << integrate( elements(mesh), _Q<3>(), idv(p) ).evaluate() << "\n";
    std::cout << integrate( elements(mesh), _Q<6>(), idv(t) ).evaluate() << "\n";
    std::cout << integrate( elements(mesh), _Q<1>(), idv(xi) ).evaluate() << "\n";

    std::cout << integrate( boundaryfaces(mesh), _Q<1>(), gradv(u)*N() ).evaluate() << "\n";
    std::cout << integrate( boundaryfaces(mesh), _Q<3>(), gradv(p)*N() ).evaluate() << "\n";
    std::cout << integrate( boundaryfaces(mesh), _Q<6>(), gradv(t)*N() ).evaluate() << "\n";
    std::cout << integrate( boundaryfaces(mesh), _Q<1>(), gradv(xi)*N() ).evaluate() << "\n";

	timers["fspace"].second=timers["fspace"].first.elapsed();
	timings << "[F spaces] Time : " << timers["fspace"].second << std::endl;

    // set up the non linear solver
    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual,
                                                   boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian,
                                                   boost::ref( *this ), _1, _2 );

    // Output for the benchmark data for each grashof number
    std::ofstream benchOut("benchmark.dat");


    D = sparse_matrix_ptrtype( M_backend->newMatrix(Xh,Xh) );
    form2( Xh,Xh, D,_init=true);

    F = vector_ptrtype( M_backend->newVector( Xh ) );
    form1( Xh, F,_init=true);

    // init to 0 and then later reuse previous grashof results to
    // initialize the solver
    u = vf::project( Xh->template functionSpace<0>(), elements(mesh), vec(constant(0.0),constant(0.0)));
    p = vf::project( Xh->template  functionSpace<1>(), elements(mesh), constant(0.0) );
    t = vf::project( Xh->template  functionSpace<2>(), elements(mesh), constant(0.0) );
    xi = vf::project( Xh->template  functionSpace<3>(), elements(mesh), constant(0.0) );


    //M_oplin->close();
    //M_oplin->mat().printMatlab( "L.m" );

    vector_ptrtype R( M_backend->newVector( Xh ) );
    sparse_matrix_ptrtype J( M_backend->newMatrix(Xh,Xh) );

    Log() << "============================================================\n";
    std::cout << "============================================================\n";
    double gr(this->vm()["gr"].template as<double>());
    M_current_Grashofs = gr;
    double pr = this->vm()["pr"].template as<double>();
    M_current_Prandtl = pr;
    Log() << "Grashof = " << M_current_Grashofs << "\n";
    Log() << "Prandtl = " << M_current_Prandtl << "\n";
    std::cout << "Grashof = " << M_current_Grashofs << "\n";
    std::cout << "Prandtl = " << M_current_Prandtl << "\n";

    int N=std::max(1.0,std::max(std::ceil(std::log(gr)),std::ceil(std::log(pr)-std::log(1e-2))));
    for( int i = 0;i < N; ++i )
        {
            int denom = (N==1)?1:N-1;
            M_current_Grashofs = math::exp( math::log(1.)+i*(math::log(gr)-math::log(1.))/denom );
            M_current_Prandtl = math::exp( math::log(1e-2)+i*(math::log(pr)-math::log(1e-2))/denom );
            std::cout << "i/N = " << i << "/" << N
                      << " intermediary Grashof = " << M_current_Grashofs
                      << " and Prandtl = " << M_current_Prandtl << "\n";
            this->solve( J, U, R );
        }

    // value mean-pressure
    double meas = integrate( elements(mesh), _Q<0>(),
                             constant(1.0) ).evaluate()( 0, 0);
    std::cout << "measure(Omega)=" << meas << " (should be equal to 1)\n";
    std::cout << "mean pressure = " << integrate(elements(mesh) ,
                                                 _Q<Order_p>() ,
                                                 idv(p) ).evaluate()(0,0)/meas << "\n";

    Log() << "value of the Lagrange multiplier xi= " << xi(0) << "\n";
    std::cout << "value of the Lagrange multiplier xi= " << xi(0) << "\n";

    double mean_div_u = integrate( elements(mesh), _Q<Order_s-1>(),
                                   divv(u) ).evaluate()( 0, 0 );
    std::cout << "mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = integrate( elements(mesh), _Q<2*(Order_s-1)>(),
                                       divv(u)*divv(u) ).evaluate()( 0, 0 );
    std::cout << "||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";

    // calcul le nombre de Nusselt
    double AverageT = integrate(markedfaces(mesh,mesh->markerName( "Tflux" )) ,
                                _Q<Order_t>() ,
                                idv(t) ).evaluate()(0,0) ;
    std::cout << "AverageT = " << AverageT << std::endl;


    double Flux = integrate(markedfaces(mesh,mesh->markerName( "Fflux" )) ,
                            _Q<Order_s>() ,
                            trans(idv(u))*vec(constant(-1.0),constant(0.0)) ).evaluate()(0,0) ;
    std::cout << "Flux = " << Flux << std::endl;

    benchOut << M_current_Grashofs << " " << AverageT << " " << Flux << std::endl;

    this->exportResults( boost::format("") , U, 0 );

	benchOut.close();

	timers["all"].second=timers["all"].first.elapsed();
	timings << "[Run] Total execution time : " << timers["all"].second << std::endl;

	timings.close();

    this->addOutputValue( AverageT ).addOutputValue( Flux ).addOutputValue( math::sqrt( div_u_error_L2 ) );
    this->postProcessing();
}

// instantiation
template class Convection<2,1,2>;
