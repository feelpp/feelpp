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

bool display_flowrates( std::ostream& os, double t, std::map<std::string, double> const& flowrates, double Q );
bool display_flowrates_header( std::ostream& os, std::map<std::string, double> const& flowrates );

int main(int argc, char**argv )
{
    //! [marker1]
    using namespace Feel;
    po::options_description nsphantomoptions( "Unsteady FDA options" );
    nsphantomoptions.add_options()
        ( "H", po::value<double>()->default_value( 0.41 ), "height of the channel" )
        ( "Um", po::value<double>()->default_value( 0.3 ), "max velocity at inflow" )
        ( "mu", po::value<double>()->default_value( 1.0 ), "viscosity" )
        ( "rho", po::value<double>()->default_value( 1.0 ), "coeff" )
        ( "Q", po::value<double>()->default_value( 0.0000052166503023447670669 ), "Flow rate" )
        ( "pulsatile", po::value<bool>()->default_value( 1.0 ), "Equal to 0 for a cst flow and equal to 1 for a pulsatile flow" )
        ( "filename", po::value<std::string>()->default_value( "data.dat" ), "Data filename" )
        ( "ns.preconditioner", po::value<std::string>()->default_value( "petsc" ), "Navier-Stokes preconditioner: petsc, PCD, PMM" )
        ;
    nsphantomoptions.add( backend_options( "ns" ) );
    Environment env( _argc=argc, _argv=argv,
                     _desc=nsphantomoptions,
                     _about=about(_name="ns_phantom_3d",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    constexpr int dim = FEELPP_DIM;
    constexpr int p_order = FEELPP_ORDER_P;
    tic();
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
    double Q = doption(_name="Q");
    auto filename=soption("filename");
    bool pulsatile = boption(_name="pulsatile");

    if ( Environment::isMasterRank() )
    {
        std::cout << "U-order\t\tP-order\t\tHsize\tFunctionSpace\tLocalDOF\tVelocity\tPressure\n";
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

        }

    }
    toc("space");tic();

    auto g = expr<FEELPP_DIM,1>( soption(_name="functions.g"), "g" );

    auto deft = gradt( u );
    auto def = grad( v );

    auto mybdf = bdf( _space=Vh, _name="mybdf" );
    U = mybdf->unknown(0);

    std::vector<std::string> flowrates_str={"inlet","outlet"};
    auto f_flowIn = form1(_test=Vh->functionSpace<0>() );
    std::map<std::string,form1_type<functionspace_type<decltype(Vh->functionSpace<0>())>>> flowrates_f;
    for( auto const& f : flowrates_str )
    {
        flowrates_f[f] = form1(_test=Vh->functionSpace<0>() );
        flowrates_f[f] = integrate(_range=markedfaces(mesh,f), _expr=inner(id(u),N()));

    }
    std::ofstream ofs( "flowrates.md" );
    std::map<std::string,double> flowrates;
    for( auto & f : flowrates_f )
    {
        flowrates[f.first]= 0;;
    }
    display_flowrates_header( ofs, flowrates );

    auto ft = form1( _test=Vh );
    auto a = form2( _trial=Vh, _test=Vh), at = form2( _trial=Vh, _test=Vh);
    a = integrate( _range=elements( mesh ), _expr=mu*inner( deft, grad(v) ) + rho*mybdf->polyDerivCoefficient(0)*trans(idt(u))*id(u) );
    a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) - divt( u )*id( q ) );
    auto e = exporter( _mesh=mesh );

    //*********************                poiseuille                *********************
    auto area_inlet = integrate(_range=markedfaces(mesh,"inlet"), _expr=cst(1.)).evaluate()( 0, 0 );
    auto poiseuille = 1-4*(Py()*Py()+Px()*Px())/25;
    auto int_poiseuille = integrate(_range=markedfaces(mesh,"inlet"), _expr=poiseuille).evaluate()( 0, 0 );
    double coeff = area_inlet/int_poiseuille;
    if ( Environment::isMasterRank() )
    {
        std::cout<<"[NAvier-Stokes]Area inlet = "<<area_inlet<<"\n";
        std::cout<<"[NAvier-Stokes]Coefficient = "<<coeff<<"\n";
    }
    //*************************************************************************************
    /*
     * retrieve vector fields from boundary condition factory
     */
    //auto dirichlet_conditions = BoundaryConditionFactory::instance().getVectorFields<dim> ( "velocity", "Dirichlet" );
    BoundaryConditions bcs;
    map_vector_field<dim,1,2> dirichlet_conditions { bcs.getVectorFields<dim> ( "velocity", "Dirichlet" ) };
    
    double temps;
    double v_mean;
    
    std::ifstream fichier(filename, std::ios::in);
    if(!fichier)
    {
        if ( Environment::isMasterRank() )
        {
            std::cout<<"Echec d'ouverture du fichier \n";
        }
    }
    else
    {
        if ( Environment::isMasterRank() )
        {
            std::cout<<"Lecture du fichier "<<filename<<".dat \n\n";
        }
    }
    if (pulsatile==1)
    {
        fichier >> temps >> v_mean;
        dirichlet_conditions.setParameterValues( { {"v",v_mean }});
        dirichlet_conditions.setParameterValues( { {"c",coeff }});
        if ( Environment::isMasterRank() )
        {
            std::cout<<"\t [NAvier-Stokes]Temps = "<<temps<<"\n";
            std::cout<<"\t [NAvier-Stokes]v_mean = "<<v_mean<<"\n";

        }
    }
    
    /*
     * Navier-Stokes block preconditioners
     */
    auto a_blockns = blockns( _space=Vh,
                              _type=soption("ns.preconditioner"),
                              _bc=BoundaryConditionFactory::instance(),
                              _matrix= a.matrixPtr(),
                              _properties_space = Pdh<0>( mesh ),
                              _alpha=cst(rho*mybdf->polyDerivCoefficient(0)),
                              _mu=cst(mu),
                              _rho=cst(rho),
                              _prefix="velocity" );

    toc("bdf, forms,...");

    for ( mybdf->restart();  mybdf->isFinished() == false; mybdf->next(U) )
    {
        if ( Environment::isMasterRank() )
        {
            std::cout << "------------------------------------------------------------\n";
            std::cout << "Time " << mybdf->time() << "s\n";
            LOG(INFO) << "------------------------------------------------------------\n";
            LOG(INFO) << "Time " << mybdf->time() << "s\n";
        }
        tic();
        tic();
        auto bdf_poly = mybdf->polyDeriv();
        auto rhsu =  bdf_poly.element<0>();
        auto extrap = mybdf->poly();
        auto extrapu = extrap.element<0>();
        // add BDF term to the right hand side from previous time steps
        ft = integrate( _range=elements(mesh), _expr=rho*(trans(idv(rhsu))*id(u) ) );
        toc("update rhs");
        tic();

        tic();
        at.zero();
        toc("update lhs zero");
        tic();
        at += a;
        toc("update lhs at += a ");
        tic();
        at += integrate( _range=elements( mesh ), _expr= rho*trans(gradt(u)*idv(extrapu))*id(v) );
        toc("update lhs at += convection ");
        tic();
        for( auto const& condition : dirichlet_conditions )
        {
            if (pulsatile==1)
            {
                dirichlet_conditions.setParameterValues( { {"v",v_mean }});
                dirichlet_conditions.setParameterValues( { {"c",coeff }});
            }
            at+=on(_range=markedfaces(mesh,marker(condition)), _rhs=ft, _element=u,_expr=expression(condition));
        }
        toc("update lhs dirichlet ");
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
        toc("Solve");

        tic();

        for( auto & f : flowrates_f )
        {
            if (f.first=="inlet" || f.first=="outlet")
                flowrates[f.first]= f.second(u);
            else
                flowrates[f.first]= 0.5*f.second(u);
        }

        if (Environment::isMasterRank())
        {
            std::cout<< "\n \n \n ============= FLOW OVER RADIAL SECTIONS =================\n";
            bool ok = display_flowrates( std::cout, mybdf->time(), flowrates, Q );
            if ( !ok )
                std::cout << "WARNING ! flow rate issue at time " << mybdf->time() << "\n";
            display_flowrates( ofs, mybdf->time(), flowrates, Q );
            std::cout<< "\n \n \n ==========================================================\n";
        }

        toc("Postprocess");
        tic();
        e->step(mybdf->time())->add( "u", u );
        e->step(mybdf->time())->add( "p", p );
        auto mean_p = mean( _range=elements(mesh), _expr=idv(p) )( 0, 0 );
        e->step(mybdf->time())->addScalar( "mean_p", mean_p );
        e->save();
        //e->restart( mybdf->time() );
        toc("Exporter");
        toc("time step");
        if( pulsatile == 1)
        {
            fichier >> temps >> v_mean;
            if ( Environment::isMasterRank() )
            {
                std::cout<<"[NAvier-Stokes]Temps = "<<temps<<"\n";
                std::cout<<"[NAvier-Stokes]v_mean = "<<v_mean<<"\n";
            }
        }
    }
    //! [marker1]

    return 0 ;
}
