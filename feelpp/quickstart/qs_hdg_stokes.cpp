/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel++ library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Daniele Prada <daniele.prada85@gmail.com>
			 Lorenzo Sala <sala@unistra.fr>
       Date: 2016-02-10

  Copyright (C) 2016-present Feel++ Consortium

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
#include <fmt/core.h>
#include <feel/feelcore/environment.hpp>
//#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vonmises.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelvf/eig.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpython/pyexpr.hpp>

namespace Feel {

inline
po::options_description
makeOptions()
{
    po::options_description hdgoptions( "test qs_hdg_stokes options" );
    hdgoptions.add_options()
        ( "pyexpr.filename", po::value<std::string>()->default_value("$cfgdir/../python/stokes.py"), "python file to evaluate" )
        ( "hsize", po::value<double>()->default_value( 0.8 ), "mesh size" )
        ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ( "mu", po::value<std::string>()->default_value( "1" ), "viscosity" )
        ( "hdg.robin", po::value<double>()->default_value( 1.0 ), "robin left hand side factor")
        ( "hface", po::value<int>()->default_value( 0 ), "hface" )
        ( "hdg.tau.constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "hdg.tau.order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
#if FEELPP_DIM==2
        ( "velocity", po::value<std::string>()->default_value( "Array([0,0])" ), "velocity is given" )
        ( "f", po::value<std::string>()->default_value( "Array([0,0])" ), "external volumic load" )
        ( "stressn", po::value<std::string>()->default_value( "Array([0,0])" ), "external surfacic load" )
#else
        ( "velocity", po::value<std::string>()->default_value( "Array([0,0,0])" ), "velocity is given" )
        ( "f", po::value<std::string>()->default_value( "Array([0,0,0])" ), "external volumic load" )
        ( "stressn", po::value<std::string>()->default_value( "Array([0,0,0])" ), "external surfacic load" )
#endif
        ( "potential", po::value<std::string>()->default_value( "0" ), "pressure is given" )
        ( "order", po::value<int>()->default_value( 1 ), "approximation order"  )
        ( "quad", po::value<int>()->default_value( 4 ), "quadrature order for rhs and norms"  )
        ( "exact", po::value<bool>()->default_value( true ), "velocity is give and provides the exact solution or an approximation"  )
        ( "use-near-null-space", po::value<bool>()->default_value( false ), "use near-null-space for AMG"  )
        ( "use-null-space", po::value<bool>()->default_value( false ), "use null-space for AMG"  )
        ;
    return hdgoptions;
}

inline
AboutData
makeAbout()
{
    AboutData about( "qs_hdg_stokes" ,
                     "qs_hdg_stokes" ,
                     "0.1",
                     "Quickstart for HDG method for Stokes problem",
                     AboutData::License_GPL,
                     "Copyright (c) 2021 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Philippe Ricka", "developer", "pricka@math.unistra.fr", "" );

    return about;

}



template<int Dim, int OrderP, int OrderG = 1>
int hdg_stokes( std::map<std::string,std::string>& locals )
{

	typedef Simplex<Dim,OrderG> convex_type;
	//! mesh type
    typedef Mesh<convex_type> mesh_type;
    using Feel::cout;

    auto tau_constant =  cst(doption("hdg.tau.constant"));
    int tau_order =  ioption("hdg.tau.order");
    auto mu = expr("1");//locals.at("mu"));
    auto r = doption("hdg.robin");

#if 0
    auto velocity_exact = expr<Dim,1>( locals.at("velocity") );
    auto grad_velocity_exact = expr<Dim,1>( locals.at("grad_velocity") );
    auto delta_exact = expr<Dim,Dim>( locals.at("stress") );
#else
    auto stressn = expr<Dim,1>(locals.at("stressn"));
    auto rhs_f = expr<Dim,1>(locals.at("f"));
    auto velocity = expr<Dim,1>(locals.at("velocity"));
    auto potential = expr(locals.at("potential"));
#endif
    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;
    
    tic();
    auto mesh = loadMesh( new Mesh<Simplex<Dim>> );
    toc("mesh",true);

    // ****** Hybrid-mixed formulation ******
    // We treat Vh, Wh, Ph, and Mh separately
    tic();

    auto Vh = Pdhms<OrderP>( mesh, true );
    auto Wh = Pdhv<OrderP>( mesh, true );
    auto Ph = Pdh<OrderP>( mesh, true );
    auto Phm = Pch<0>( mesh, true );
    auto face_mesh = createSubmesh( _mesh=mesh, _range=faces(mesh), _update=0 );
    auto Mh = Pdhv<OrderP>( face_mesh,true );

    toc("spaces",true);

    cout << "Vh<" << OrderP << "> : " << Vh->nDof() << std::endl
         << "Wh<" << OrderP << "> : " << Wh->nDof() << std::endl
         << "Ph<" << OrderP << "> : " << Ph->nDof() << std::endl
         << "Phm<" << OrderP << "> : " << Phm->nDof() << std::endl
         << "Mh<" << OrderP << "> : " << Mh->nDof() << std::endl;

    auto delta = Vh->element( "delta" );
    auto gamma = Vh->element( "gamma" );
    auto u     = Wh->element( "u" );
    auto v     = Wh->element( "v" );
    auto p     = Ph->element( "p" );
    auto q     = Ph->element( "q" );
    auto qm    = Phm->element( "qm" );
    auto uhat  = Mh->element( "uhat" );
    auto m     = Mh->element( "m" );

    // Number of dofs associated with each space
    auto nDofdelta = delta.functionSpace()->nDof();
    auto nDofu     = u.functionSpace()->nDof();
    auto nDofp     = p.functionSpace()->nDof();
    auto nDofuhat  = uhat.functionSpace()->nDof();

    solve::strategy strategy = boption("sc.condense")?solve::strategy::static_condensation:solve::strategy::monolithic;
	double sc_param = 1;
    if( boption("sc.condense") )
    {
        sc_param = 0.5;
    }
	tic();
    auto ps = product( Vh, Wh, Ph, Mh, Phm );
	auto a = blockform2( ps, strategy , backend() );
	auto rhs = blockform1( ps, strategy , backend() );


    // Building the RHS
    auto M0h = Pdh<0>( face_mesh,true );
    auto H     = M0h->element( "H" );
    if ( ioption("hface" ) == 0 )
        H.on( _range=elements(face_mesh), _expr=pow(mesh->hMax(),tau_order) );
    else if ( ioption("hface" ) == 1 )
        H.on( _range=elements(face_mesh), _expr=pow(mesh->hMin(),tau_order) );
    else if ( ioption("hface" ) == 2 )
        H.on( _range=elements(face_mesh), _expr=pow(mesh->hAverage(),tau_order) );
    else
        H.on( _range=elements(face_mesh), _expr=pow(h(),tau_order) );

    bool bc_only_dirichlet = ( nelements(markedfaces(mesh, "Dirichlet")) == nelements(boundaryfaces(mesh)) );
    tic();
    rhs( 1_c ) += integrate( _range = elements( mesh ), _quad = ioption( "quad" ),
                             _expr = trans( rhs_f ) * id( v ) );
    rhs(3_c) += integrate(_range=markedfaces(mesh, {"Neumann", "Robin"}),_quad=ioption("quad"),
                          _expr= (-1)*inner(stressn,id(m)) );
    rhs(3_c) += integrate(_range=markedfaces(mesh, {"Dirichlet", "Robin"}),_quad=ioption("quad"),
                          _expr=trans(velocity)*id(m) );
    rhs( 4_c ) += integrate( _range = elements(mesh),
                             _expr = 0 * id( qm ) );
    toc("rhs",true);

    // Building the matrix
    tic();
    a( 0_c, 0_c ) +=  integrate(_range=elements(mesh),
                                _expr=inner(idt(delta),id(gamma)) );
    toc("a(0,0)", true);

    tic();
    a( 0_c, 1_c ) += integrate(_range=elements(mesh),_expr=(trans(idt(u))*div(gamma)) );
    toc("a(0,1)", true);
    tic();
    a( 0_c, 3_c) += integrate(_range=internalfaces(mesh),
                              _expr=-(trans(idt(uhat))*leftface((id(gamma)*N()))+
                                      trans(idt(uhat))*rightface((id(gamma)*N())) ));
    a( 0_c, 3_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=-trans(idt(uhat))*(id(gamma)*N()) );
    toc("a(0,3)", true);
    tic();
    a( 1_c, 0_c) += integrate(_range=elements(mesh),
                              _expr=2*mu*inner(idt(delta),grad(v)));
    a( 1_c, 0_c) += integrate(_range=internalfaces(mesh),
                              _expr=-2*mu*(trans(leftface(id(v)))*leftfacet((idt(delta)*N()))+
                                         trans(rightface(id(v)))*rightfacet((idt(delta)*N())) ));
    a( 1_c, 0_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=-2*mu*trans(id(v))*(idt(delta)*N()) );                 
    toc("a(1,0)", true);

    a( 1_c, 1_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=mu*tau_constant*trans(id(v))*idt(u) );
    a( 1_c, 1_c) += integrate(_range=internalfaces(mesh),
                              _expr=mu*tau_constant*(trans(leftface(id(v)))*leftfacet(idt(u))+
                                                     trans(rightface(id(v)))*rightfacet(idt(u))) );

    tic();
    a( 1_c, 2_c) += integrate(_range=elements(mesh),
                              _expr=(-1)*idt(p)*div(v) );
    a( 1_c, 2_c) += integrate(_range=internalfaces(mesh),
                              _expr=inner(jumpt(idt(p)),id(v)) );
    a( 1_c, 2_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=normal(v)*idt(p) );
    toc("a(1,2)", true);
    tic();
    a( 1_c, 3_c) += integrate(_range=internalfaces(mesh),
                              _expr=-mu*tau_constant*(trans(leftface(id(v)))*idt(uhat)+
                                                   trans(rightface(id(v)))*idt(uhat)) );
    a( 1_c, 3_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=-mu*tau_constant*(trans(id(v))*idt(uhat)) );
    toc("a(1,3)", true);
    tic();
    a( 2_c, 1_c) += integrate(_range=elements(mesh),
                              _expr=(-1)*grad(q)*idt(u) );
    toc("a(2,1)", true);
    tic();
    a( 2_c, 3_c) += integrate(_range=internalfaces(mesh),
                              _expr=inner(jump(id(q)),idt(uhat)) );
    a( 2_c, 3_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=id(q)*(trans(idt(uhat))*N()) );
    if ( bc_only_dirichlet )
    {
        a( 2_c, 4_c) += integrate(_range=elements(mesh),
                                  _expr=idt(qm)*id(q) );       
        a( 4_c, 2_c) += integrate(_range=elements(mesh),
                                  _expr=id(qm)*idt(q) ); 
    }
    else
    {
        a( 4_c, 4_c) += integrate(_range=elements(mesh),
                              _expr=id(qm)*idt(qm) );
    } 
    toc("a(2,3)", true);
    tic();
    a( 3_c, 0_c) += integrate(_range=internalfaces(mesh),
                              _expr=-2*mu*(trans(id(m))*leftfacet(idt(delta)*N())+
                                         trans(id(m))*rightfacet(idt(delta)*N())) );
//    a( 3_c, 0_c) += integrate(_range=boundaryfaces(mesh),
//                              _expr=-mu*trans(id(m))*idt(delta)*N());
    a( 3_c, 0_c ) += integrate(_range = markedfaces( mesh, {"Neumann", "Robin"} ),
                               _expr = -2*mu * trans(id(m))*(idt(delta)*N() ));
    toc("a(3,0)", true);
    tic();
    a( 3_c, 1_c) += integrate(_range=internalfaces(mesh),
                              _expr=mu*tau_constant*(trans(leftfacet(idt(u)))*id(m)+
                                                   trans(rightfacet(idt(u)))*id(m)) );
    a( 3_c, 1_c ) += integrate( _range = markedfaces( mesh, {"Neumann", "Robin"} ),
                                _expr = mu * tau_constant * ( trans( idt( u ) ) * id( m ) ) );
    toc("a(3,1)", true);
    tic();
    a( 3_c, 2_c) += integrate(_range=internalfaces(mesh),
                              _expr=inner(id(m),jumpt(idt(p))) ); //normalt(m)*(rightface(id(p))-leftface(id(p))) );
    a( 3_c, 2_c) += integrate(_range=markedfaces(mesh,{"Neumann", "Robin"}),
                              _expr=idt(p)*trans(id(m))*N() );
    toc("a(3,2)", true);
    tic();
    a( 3_c, 3_c) += integrate(_range=internalfaces(mesh),                  
                              _expr=-sc_param*mu*tau_constant*trans(idt(uhat))*id(m) );
    a( 3_c, 3_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                              _expr=trans(idt(uhat))*id(m) );
    a( 3_c, 3_c) += integrate(_range=markedfaces(mesh,"Robin"),
                              _expr=r*trans(idt(uhat))*id(m) );
    a( 3_c, 3_c ) += integrate( _range = markedfaces( mesh, {"Neumann", "Robin"} ),
                                _expr = -mu*tau_constant*trans( idt( uhat ) ) * id( m ) );
    toc("a(3,3)", true);
    tic();
#if 0
    if ( boption( "use-near-null-space" ) ||  boption( "use-null-space" ) )
    {
        auto b = backend( _name="sc" );
        std::shared_ptr<NullSpace<double> > myNullSpace( new NullSpace<double>(b,qsNullSpace(Mh,mpl::int_<FEELPP_DIM>())) );
        if ( boption( "use-near-null-space" ) )
            b->attachNearNullSpace( myNullSpace );
        if ( boption( "use-null-space" ) )
            b->attachNullSpace( myNullSpace );
    }
#endif

    auto U = ps.element();
    auto Ue = ps.element();
    //a.solve( _solution=U, _rhs=rhs, _rebuild=true, _condense=boption("sc.condense"));
    a.solve( _solution=U, _rhs=rhs, _condense=boption("sc.condense"), _condenser=condenser_stokes() );
    toc("solve",true);
    cout << "[Hdg] solve done" << std::endl;

    auto deltap = U(0_c);
    auto up = U(1_c);
    auto pp = U(2_c);
    auto uhatp = U(3_c);
    if ( Environment::isSequential() && boption("exporter.matlab") )
    {
        deltap.printMatlab("s");
        up.printMatlab("u");
        pp.printMatlab("p");
        uhatp.printMatlab("uhat");
    }
    int status_velocity = 1, status_stress = 1;
    if ( boption( "exact" ) )
    {
        auto pressure_exact = potential;
        auto velocity_exact = velocity;
        auto delta_exact = expr<Dim,Dim>(locals.at("strain"));
        auto grad_velocity_exact = expr<Dim,Dim>(locals.at("grad_velocity"));
        Ue(0_c).on( _range=elements(mesh), _expr=delta_exact );
        Ue(1_c).on( _range=elements(mesh), _expr=velocity_exact );
        Ue(2_c).on( _range=elements(mesh), _expr=pressure_exact );
        Ue(3_c).on( _range=faces(mesh), _expr=velocity_exact );
        if ( Environment::isSequential() && boption("exporter.matlab") )
        {
            Ue(0_c).printMatlab("se");
            Ue(1_c).printMatlab("ue");
            Ue(2_c).printMatlab("pe");
            Ue(3_c).printMatlab("uhate");
        }
        
        auto l2err_delta = normL2( _range=elements(mesh), _expr=delta_exact - idv(deltap),_quad=ioption("quad") );
        auto l2err_vel = normL2( _range=elements(mesh), _expr=velocity_exact - idv(up),_quad=ioption("quad") );
        auto mean_pe = mean( _range = elements( mesh ), _expr = pressure_exact,_quad=ioption("quad") )(0,0);
        auto mean_p = mean( _range = elements( mesh ), _expr = idv(pp) )(0,0);
        auto l2err_pres = normL2( _range = elements( mesh ), _expr = ( pressure_exact-(bc_only_dirichlet*mean_pe)) - idv( pp ),_quad=ioption("quad") );
        auto l2vel = normL2( _range = elements( mesh ), _expr = velocity_exact, _quad = ioption( "quad" ) );
        Feel::cout << fmt::format( "{:<30}: {: .4e}", "L2 Error strain", l2err_delta ) << std::endl;
        Feel::cout << fmt::format( "{:<30}: {: .4e}", "L2 Error velocity", l2err_vel ) << std::endl;
        Feel::cout << fmt::format( "{:<30}: {: .4e}", "L2 velocity", l2vel ) << std::endl;
        if ( std::abs(l2vel) > 1e-10 )
            Feel::cout << fmt::format( "{:<30}: {: .4e}", "L2 relative error velocity", l2err_vel/l2vel ) << std::endl;
        Feel::cout << fmt::format( "{:<30}: {: .4e}", "mean pressure exact", mean_pe ) << std::endl;
        Feel::cout << fmt::format( "{:<30}: {: .4e}", "mean pressure", mean_p ) << std::endl;
        Feel::cout << fmt::format( "{:<30}: {: .4e}", "L2 error pressure", l2err_pres ) << std::endl;
        
        auto h = doption("hsize");
        ofstream errors;
        errors.open ("hdg_stokes_errors.csv");
        errors << h << l2err_delta << l2err_vel << l2err_pres <<"\n";
        errors.close();
        toc("error");
                    
        // CHECKER
        auto norms_stress = [&]( std::string const& solution ) ->std::map<std::string,double>
            {
                tic();
                double l2 = normL2( _range=elements(mesh), _expr=delta_exact - idv(deltap) );
                toc("L2 stress error norm");
                
                return { { "L2", l2 } };
            };
        // compute l2 and h1 norm of u-u_h where u=solution
        auto norms_velocity = [&]( std::string const& solution ) ->std::map<std::string,double>
            {
			tic();
			double l2 = normL2( _range=elements(mesh), _expr=velocity_exact - idv(up) );
			toc("L2 velocity error norm");

			tic();
			double h1 = normH1(_range=elements(mesh), _expr=idv(up)- velocity_exact, _grad_expr=gradv(up)-grad_velocity_exact );
			toc("H1 velocity error norm");

			return { { "L2", l2 } , {  "H1", h1 } };
            };
#if 0
        status_velocity = checker("L2/H1 velocity norms",velocity_exact).runOnce( norms_velocity, rate::hp( mesh->hMax(), Wh->fe()->order() ) );
        status_stress = checker("L2 stress norms",velocity_exact).runOnce( norms_stress, rate::hp( mesh->hMax(), Vh->fe()->order() ) );
#endif        
        delta.on( _range=elements(mesh), _expr=delta_exact );
        u.on( _range=elements(mesh), _expr=velocity_exact );
        p.on( _range=elements(mesh), _expr=pressure_exact );
    }

    tic();
    std::string exportName =  "hdg_stokes";
    std::string deltaName = "stress";
    std::string delta_exName = "stress-ex";
    std::string uName = "velocity";
    std::string u_exName = "velocity-ex";
    std::string pName = "pressure";
    std::string p_exName = "pressure-ex";
    auto e = exporter( _mesh=mesh, _name=exportName );
    e->setMesh( mesh );
    e->addRegions();
    e->add( deltaName, deltap, "nodal" );
    e->add( uName, up, "nodal" );
    e->add( uName, pp, "nodal" );
    std::set<std::string> reps({ "nodal", "element" });
    e->add( "vonmises", vonmises(idv(deltap)), reps );
    e->add( "principal_stress", eig(idv(deltap)), reps );
    e->add( "magnitude_stress", sqrt(inner(idv(deltap))), reps );
    
    if ( boption("exact" ) )
    {
        e->add( delta_exName, delta, "nodal" );
        e->add( u_exName, u, "nodal" );
        e->add( p_exName, p, "nodal" );
    }

    e->save();

    toc("export");
#if 0
    return status_stress || status_velocity;
#endif
    return 1;
}

} // Feel

int main( int argc, char** argv )
{
    // tag::env[]
    using namespace Feel;

	Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="qs_hdg_stokes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    // end::env[]

    // Exact solutions
    std::map<std::string,std::string> locals{
        {"dim",std::to_string(FEELPP_DIM)},
        {"exact",std::to_string(boption("exact"))}, 
        {"mu",soption("mu")},
        {"potential",soption("potential")},
        {"velocity", soption("velocity")},
        {"grad_velocity",""},
        {"strain",""},
        {"stress",""},
        {"stressn",soption("stressn")},
        {"f",soption("f")}};
    Feel::pyexprFromFile( Environment::expand(soption("pyexpr.filename")), locals  );


    for( auto d: locals )
        Feel::cout << fmt::format( " -- symbol {} expression {}\n", d.first, d.second) << std::endl;

    if ( ioption( "order" ) == 1 )
        return !hdg_stokes<FEELPP_DIM,1>( locals );
#if 1        
    if ( ioption( "order" ) == 2 )
        return !hdg_stokes<FEELPP_DIM,2>( locals );
#endif
    return 0;
}
