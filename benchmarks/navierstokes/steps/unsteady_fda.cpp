/* -*- Mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date     : Tue Feb 25 12:13:15 2014

   Copyright (C) 2014-2015 Feel++ Consortium

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
#include <feel/feel.hpp>
#include <feel/feelpde/preconditionerblockns.hpp>

int main(int argc, char**argv )
{
    //! [marker1]
    using namespace Feel;
    po::options_description unsteadyfdaoptions( "Unsteady FDA options" );
    unsteadyfdaoptions.add_options()
        ( "H", po::value<double>()->default_value( 0.41 ), "height of the channel" )
        ( "Um", po::value<double>()->default_value( 0.3 ), "max velocity at inflow" )
        ( "mu", po::value<double>()->default_value( 1.0 ), "viscosity" )
        ( "rho", po::value<double>()->default_value( 1.0 ), "coeff" )
        ( "ns.preconditioner", po::value<std::string>()->default_value( "petsc" ), "Navier-Stokes preconditioner: petsc, PCD, PMM" )
        ;
    unsteadyfdaoptions.add( backend_options( "ns" ) );
    Environment env( _argc=argc, _argv=argv,
                     _desc=unsteadyfdaoptions,
                     _about=about(_name="unsteady_fda_3d",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    constexpr int dim = FEELPP_DIM;
    constexpr int p_order = FEELPP_ORDER_P;
    tic();
    double Di=0.012;
    auto mesh = loadMesh( new Mesh<Simplex<dim>> );
    size_type nbdyfaces = nelements(boundaryfaces(mesh));
    if ( Environment::isMasterRank() )
    {
        std::cout << " - mesh entities" << std::endl;
        std::cout << "      number of elements : " << mesh->numGlobalElements() << std::endl;
        std::cout << "         number of faces : " << mesh->numGlobalFaces() << std::endl;
        std::cout << "number of boundary faces : " << nbdyfaces << std::endl;
        if ( FEELPP_DIM > 2 )
            std::cout << "      number of edges : " << mesh->numGlobalEdges() << std::endl;
        std::cout << "      number of points : " << mesh->numGlobalPoints() << std::endl;
        std::cout << "    number of vertices : " << mesh->numGlobalVertices() << std::endl;
        std::cout << " - mesh sizes" << std::endl;
        std::cout << "                h max : " << mesh->hMax() << std::endl;
        std::cout << "                h min : " << mesh->hMin() << std::endl;
        std::cout << "                h avg : " << mesh->hAverage() << std::endl;
        std::cout << "              measure : " << mesh->measure() << std::endl;
    }

    //CHECK( mesh->hasMarkers( {"wall","inlet"} ) ) << "Mesh markers wall or inlet are not set properly in "  << soption("gmsh.filename");
    toc("mesh");tic();
    // Taylor Hood P_N+1 velocity P_N  pressure space (N==p_order)
    auto Vh = THch<p_order>( mesh );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();
    
    double mu = doption(_name="mu");
    double rho = doption(_name="rho");
    
    if ( Environment::isMasterRank() )
    {
        std::cout << "Re\t\tU-order\t\tP-order\t\tHsize\tFunctionSpace\tLocalDOF\tVelocity\tPressure\n";
        std::cout.width(16);
        std::cout << std::left << 0.0461*0.012*rho/mu;
        std::cout.width(16);
        std::cout << std::left << p_order+1;
        std::cout.width(16);
        std::cout << std::left << p_order;
        std::cout << std::left << doption( "gmsh.hsize" )<<"\t";
        std::cout.width(16);
        std::cout << std::left << Vh->nDof();
        std::cout.width(16);
        std::cout << std::left << Vh->nLocalDof();
        std::cout.width(16);
        std::cout << std::left << Vh->functionSpace<0>()->nDof();
        std::cout.width(16);
        std::cout << std::left << Vh->functionSpace<1>()->nDof() << "\n";
        if( soption("ns.preconditioner") != "petsc" )
        {
            std::cout << "[blockns]\n";
            std::cout << " - cd: " << boption( "blockns.cd" ) << "\n";
            std::cout << " - pcd: " << boption( "blockns.pcd" ) << "\n";
            std::cout << " - pcd.inflow: " << soption( "blockns.pcd.inflow" ) << "\n";
            std::cout << " - pcd.outflow: " << soption( "blockns.pcd.outflow" ) << "\n";
            std::cout << " - pcd.order: " << ioption( "blockns.pcd.order" ) << "\n\n";
            
            /*std::cout << " - Stokes rtol: " << doption( "stokes.ksp-rtol" ) << "\n";
            std::cout << " - Picard rtol: " << doption( "picard.ksp-rtol" ) << "\n\n";
            
            std::cout << " - Ap preconditioner: " << soption( "Ap.pc-type" ) << "\n";
            std::cout << " - Ap relative tolerence: " << doption( "Ap.ksp-rtol" ) << "\n";
            std::cout << " - Ap reuse-prec: " << boption( "Ap.reuse-prec" ) << "\n\n";
            
            std::cout << " - Mp preconditioner: " << soption( "Mp.pc-type" ) << "\n";
            std::cout << " - Mp relative tolerence: " << doption( "Mp.ksp-rtol" ) << "\n";
            std::cout << " - Mp reuse-prec: " << boption( "Mp.reuse-prec" ) << "\n\n";
            
            std::cout << " - Fu preconditioner: " << soption( "Fu.pc-type" ) << "\n";
            std::cout << " - Fu relative tolerence: " << doption( "Fu.ksp-rtol" ) << "\n";
            std::cout << " - Fu reuse-prec: " << boption( "Fu.reuse-prec" ) << "\n\n";*/

        }
        
    }

    toc("space");tic();
    
    auto g = expr<FEELPP_DIM,1>( soption(_name="functions.g"), "g" );
    
    auto deft = gradt( u );
    auto def = grad( v );

    auto mybdf = bdf( _space=Vh, _name="mybdf" );
    U = mybdf->unknown(0);
    
    auto ft = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh), at = form2( _trial=Vh, _test=Vh);

    a = integrate( _range=elements( mesh ), _expr=mu*inner( deft, grad(v) ) + rho*mybdf->polyDerivCoefficient(0)*trans(idt(u))*id(u) );
    a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) - divt( u )*id( q ) );
    auto e = exporter( _mesh=mesh );
    auto w = Vh->functionSpace<0>()->element( curlv(u), "w" );

    /*
     * retrieve vector fields from boundary condition factory
     */
    //auto dirichlet_conditions = BoundaryConditionFactory::instance().getVectorFields<dim> ( "velocity", "Dirichlet" );
    BoundaryConditions bcs;
    map_vector_field<dim,1,2> dirichlet_conditions { bcs.getVectorFields<dim> ( "velocity", "Dirichlet" ) };
  
    /*
     * Navier-Stokes block preconditioners
     */
    auto a_blockns = blockns( _space=Vh,
                              _type=soption("ns.preconditioner"),
                              _bc=BoundaryConditionFactory::instance(),
                              _matrix= a.matrixPtr(),
                              _alpha=rho*mybdf->polyDerivCoefficient(0),
                              _mu=mu,
                              _rho=rho,
                              _prefix="velocity" );
    
    toc("bdf, forms,...");

    for ( mybdf->restart();  mybdf->isFinished() == false; mybdf->next(U) )
    {
        if ( Environment::isMasterRank() )
        {
            std::cout << "------------------------------------------------------------\n";
            std::cout << "Time " << mybdf->time() << "s\n";
        }
        tic();
        tic();
        auto bdf_poly = mybdf->polyDeriv();
        auto rhsu =  bdf_poly.element<0>();
        auto extrap = mybdf->poly();
        auto extrapu = extrap.element<0>();
        // add BDF term to the right hand side from previous time steps
        ft = integrate( _range=elements(mesh), _expr=rho*(trans(idv(rhsu))*id(u) ) );
        toc("update rhs");tic();

        at.zero();
        at += a;
        at += integrate( _range=elements( mesh ), _expr= rho*trans(gradt(u)*idv(extrapu))*id(v) );
        for( auto const& condition : dirichlet_conditions )
        {
            dirichlet_conditions.setParameterValues( { {"t",mybdf->time()}}  );
            at+=on(_range=markedfaces(mesh,marker(condition)), _rhs=ft, _element=u,_expr=expression(condition));
            /*at+=on(_range=markedfaces(mesh,"wall"), _rhs=ft, _element=u,
            _expr=zero<FEELPP_DIM,1>() );
             at+=on(_range=markedfaces(mesh,"inlet"), _rhs=ft, _element=u, _expr=-g*N() );*/
        }
        toc("update lhs");tic();

        if ( soption("ns.preconditioner") != "petsc" )
        {
            a_blockns->update( at.matrixPtr(), idv(extrapu), dirichlet_conditions );
            at.solveb(_rhs=ft,_solution=U,_backend=backend(_name="ns",_rebuild=false),_prec=a_blockns);
        }
        else
        {
            // use petsc preconditioner
            at.solveb(_rhs=ft,_solution=U,_backend=backend(_name="ns"));
        }
        toc("Solve");tic();
        Environment::saveTimers( true );

        /*double areaIn = integrate(_range=markedfaces(mesh,"inlet"), _expr=vf::cst(1.) ).evaluate()(0,0);
        double areaOut = integrate(_range=markedfaces(mesh,"outlet"), _expr=vf::cst(1.) ).evaluate()(0,0);

        auto Uz= minmax(_range=markedfaces(mesh,"inlet"),_pset=_Q<2>(),_expr=idv(u.comp(Component::X)));
        double ReynoldsIn =rho*Uz.get<1>()*Di/mu;

        auto intUz = integrate(_range=markedfaces(mesh,"inlet"), _expr=idv(u.comp(Component::X)) ).evaluate()(0,0) ;
        auto meanU = intUz/areaIn;

        double nElt = nelements(markedfaces(mesh, "face1"),true);
        double area1 = integrate(_range=markedfaces(mesh,"face1"), _expr=cst(1.) ).evaluate()(0,0);*/

        auto flowIn = integrate(_range=markedfaces(mesh,"inlet"), _expr=inner(idv(u),N())).evaluate()(0,0) ;
        auto flowOut = integrate(_range=markedfaces(mesh,"outlet"), _expr=inner(idv(u),N())).evaluate()(0,0) ;

        auto flow1 = integrate(_range=markedfaces(mesh,"face1"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow2 = integrate(_range=markedfaces(mesh,"face2"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow3 = integrate(_range=markedfaces(mesh,"face3"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow4 = integrate(_range=markedfaces(mesh,"face4"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow5 = integrate(_range=markedfaces(mesh,"face5"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow6 = integrate(_range=markedfaces(mesh,"face6"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow7 = integrate(_range=markedfaces(mesh,"face7"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow8 = integrate(_range=markedfaces(mesh,"face8"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow9 = integrate(_range=markedfaces(mesh,"face9"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow10 = integrate(_range=markedfaces(mesh,"face10"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow11 = integrate(_range=markedfaces(mesh,"face11"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;
        auto flow12 = integrate(_range=markedfaces(mesh,"face12"), _expr=inner(idv(u),oneZ() )).evaluate()(0,0) ;

        //double Reynolds =rho*meanU*Di/mu;

        if (Environment::isMasterRank())
        {

            /*std::cout<<  "     Time  "
                     << "flowIn  " << "flowOut "  << "sum(flowIn+flowOut) "
                     << "areaIn "
                     << "areaOut "
                     << "Reynolds"
                     << "\n";

            std::cout <<  mybdf->time() << " "
                      << std::scientific
                      << flowIn << " " << flowOut << " " << std::abs(flowIn+flowOut) << " "
                      << areaIn << " "
                      << areaOut << " "
                      <<ReynoldsIn<<" "
                      << "\n";*/

            std::cout<< "\n \n \n ============= FLOW OVER RADIAL SECTIONS =================\n";
            std::cout<< "Flow in: "<<flowIn<<"\n";
            std::cout<< "Flow out: "<<flowOut<<"\n";
            /*std::cout<< "Nb elem face 1: "<<nElt<<"\n";
            std::cout<< "Area of face 1: "<<area1<<"\n";*/
            std::cout<< "Flow 1: "<<flow1<<"  ( "<< ((flow1-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 2: "<<flow2<<"  ( "<< ((flow2-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 3: "<<flow3<<"  ( "<< ((flow3-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 4: "<<flow4<<"  ( "<< ((flow4-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 5: "<<flow5<<"  ( "<< ((flow5-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 6: "<<flow6<<"  ( "<< ((flow6-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 7: "<<flow7<<"  ( "<< ((flow7-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 8: "<<flow8<<"  ( "<< ((flow8-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 9: "<<flow9<<"  ( "<< ((flow9-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 10: "<<flow10<<"  ( "<< ((flow10-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 11: "<<flow11<<"  ( "<< ((flow11-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "Flow 12: "<<flow12<<"  ( "<< ((flow12-5.216650303e-6)/5.216650303e-6)*100<<"  ) \n";
            std::cout<< "\n \n \n ==========================================================\n";

            }


        w.on( _range=elements(mesh), _expr=curlv(u) );
        e->step(mybdf->time())->add( "u", u );
        //e->step(mybdf->time())->add( "w", w );
        e->step(mybdf->time())->add( "p", p );
        auto mean_p = mean( _range=elements(mesh), _expr=idv(p) )( 0, 0 );
        e->step(mybdf->time())->addScalar( "mean_p", mean_p );
        e->save();
        //e->restart( mybdf->time() );
        toc("export");
        toc("time step");


    }
    //! [marker1]



    return 0;
}
