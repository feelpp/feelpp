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
    po::options_description unsteadyfdaoptions( "Unsteady FDA options" );
    unsteadyfdaoptions.add_options()
        ( "H", po::value<double>()->default_value( 0.41 ), "height of the channel" )
        ( "Um", po::value<double>()->default_value( 0.3 ), "max velocity at inflow" )
        ( "mu", po::value<double>()->default_value( 1.0 ), "viscosity" )
        ( "rho", po::value<double>()->default_value( 1.0 ), "coeff" )
        ( "Q", po::value<double>()->default_value( 0.0000052166503023447670669 ), "Flow rate" )
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
    constexpr int geo_order = FEELPP_ORDER_GEO;
    
    tic();
    auto mesh = loadMesh( new Mesh<Simplex<dim,geo_order,dim>> );
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
    Environment::saveTimers(true);
    
    //CHECK( mesh->hasMarkers( {"wall","inlet"} ) ) << "Mesh markers wall or inlet are not set properly in "  << soption("gmsh.filename");
    toc("mesh");tic();
    // Taylor Hood P_N+1 velocity P_N  pressure space (N==p_order)
    auto Vh = THch<p_order>( mesh );
    //auto P1h = lagrangeP1( _space=Vh->functionSpace<0>() );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();
    
    double mu = doption(_name="mu");
    double rho = doption(_name="rho");
    double Q = doption(_name="Q");
    
    if ( Environment::isMasterRank() )
    {
        std::cout << "Re\t\tU-order\t\tP-order\t\tHsize\tFunctionSpace\tLocalDOF\tVelocity\tPressure\n";
        std::cout.width(16);
        std::cout << std::left << Q*0.012*rho/(pi*0.006*0.006*mu);
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
    }

    toc("space");tic();
    
    auto g = expr<FEELPP_DIM,1>( soption(_name="functions.g"), "g" );
    
    auto deft = gradt( u );
    auto def = grad( v );

    auto mybdf = bdf( _space=Vh, _name="mybdf" );
    U = mybdf->unknown(0);
 
    std::vector<std::string> flowrates_str={"inlet","outlet","face1","face2","face3","face4","face5","face6","face7","face8","face9","face10","face11","face12"};
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
    //auto e = exporter( _mesh=mesh );
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
                              _properties_space = Pdh<0>( mesh ),
                              _alpha=cst(rho*mybdf->polyDerivCoefficient(0)),
                              _mu=cst(mu),
                              _rho=cst(rho),
                              _prefix="velocity" );
    
    toc("bdf, forms,...");
    auto b = backend(_prefix="ns",_name="ns");
    
    auto precPetsc = preconditioner( _prefix="ns",_matrix=a.matrixPtr(),_pc=b->pcEnumType(),
                                    _pcfactormatsolverpackage=b->matSolverPackageEnumType(), _backend=b->shared_from_this(),
                                    _worldcomm=b->comm() );


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
            dirichlet_conditions.setParameterValues( { {"t",mybdf->time()}}  );
            at+=on(_range=markedfaces(mesh,marker(condition)), _rhs=ft, _element=u,_expr=expression(condition));
            /*at+=on(_range=markedfaces(mesh,"wall"), _rhs=ft, _element=u,
             _expr=zero<FEELPP_DIM,1>() );
             at+=on(_range=markedfaces(mesh,"inlet"), _rhs=ft, _element=u, _expr=-g*N() );*/
        }
        toc("update lhs dirichlet ");
        toc("update lhs");tic();

        if ( soption("ns.preconditioner") != "petsc" )
        {
            a_blockns->update( at.matrixPtr(), rho*idv(extrapu), dirichlet_conditions );
            at.solveb(_rhs=ft,_solution=U,_backend=backend(_name="ns",_rebuild=false),_prec=a_blockns);
        }
        else
        {
            // use petsc preconditioner
            //at.solveb(_rhs=ft,_solution=U,_backend=backend(_name="ns"));
            backend(_name="ns")->solve(_matrix=at.matrixPtr(),_solution=U,_rhs=ft.vectorPtr(),_prec=precPetsc );
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
        
        auto e = exporter( _mesh=mesh );
        e->step(mybdf->time())->add( "u", u );
        w.on( _range=elements(mesh), _expr=curlv(u) );
        e->step(mybdf->time())->add( "w", w );
        e->step(mybdf->time())->add( "p", p );
        auto mean_p = mean( _range=elements(mesh), _expr=idv(p) )( 0, 0 );
        e->step(mybdf->time())->addScalar( "mean_p", mean_p );
        e->save();
        Environment::saveTimers( true );
            
    }
        toc("Exporter");
        toc("time step");
    
    tic();
    
    if(geo_order>1 || p_order>1)
    {
        auto meshP1 = lagrangeP1(_space=Vh->template functionSpace<0>())->mesh();
        auto XhVisuU = Pchv<1>(meshP1,true);
        auto XhVisuP = Pch<1>(meshP1,true);
        auto opIVisuU = opInterpolation(_domainSpace=Vh->template functionSpace<0>(),
                                        _imageSpace=XhVisuU,
                                        _type=InterpolationNonConforme(false,true,false) );
        auto opIVisuP = opInterpolation(_domainSpace=Vh->template functionSpace<1>(),
                                        _imageSpace=XhVisuP,
                                        _type=InterpolationNonConforme(false,true,false) );
        
        auto uVisu = opIVisuU->operator()(u);
        auto pVisu = opIVisuP->operator()(p);
        
        // exporter mesh and harmonic extension
        auto e2 = exporter( _mesh=meshP1, _name="VisuP1" );
        e2->step(mybdf->time())->setMesh( meshP1 );
        e2->step(mybdf->time())->add( "u", uVisu );
        e2->step(mybdf->time())->add( "p", pVisu );
        e2->save();

    }
    toc("ExporterHO");
    Environment::saveTimers( true );



    return 0;
}
