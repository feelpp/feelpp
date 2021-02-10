/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel++ library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Daniele Prada <daniele.prada85@gmail.com>
			 Lorenzo Sala <sala@unistra.fr>
       Date: 2016-02-10

  Copyright (C) 2016 Feel++ Consortium

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
#include <feel/feel.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelvf/vf.hpp>
namespace Feel {

inline
po::options_description
makeOptions()
{
    po::options_description testhdivoptions( "test qs_hdg_stokes options" );
    testhdivoptions.add_options()
        ( "pyexpr.filename", po::value<std::string>()->default_value("${top_srcdir}/quickstart/stokes.py"), "python file to evaluate" )
        ( "hsize", po::value<double>()->default_value( 0.8 ), "mesh size" )
        ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ( "mu", po::value<double>()->default_value( 1 ), "viscosity" )
        ( "hface", po::value<int>()->default_value( 0 ), "hface" )
        ( "hdg.tau.constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "hdg.tau.order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
#if FEELPP_DIM==2
        ( "u", po::value<std::string>()->default_value( "Array([0,0])" ), "velocity is given" )
        ( "f", po::value<std::string>()->default_value( "Array([0,0])" ), "external volumic load" )
        ( "stressn", po::value<std::string>()->default_value( "Array([0,0])" ), "external surfacic load" )
#else
        ( "u", po::value<std::string>()->default_value( "Array([0,0,0])" ), "velocity is given" )
        ( "f", po::value<std::string>()->default_value( "Array([0,0,0])" ), "external volumic load" )
        ( "stressn", po::value<std::string>()->default_value( "Array([0,0,0])" ), "external surfacic load" )
#endif
        ;
    return testhdivoptions.add( Feel::feel_options() ).add( backend_options("sc"));
}

inline
AboutData
makeAbout()
{
    AboutData about( "qs_hdg_stokes" ,
                     "qs_hdg_stokes" ,
                     "0.1",
                     "Quickstart for HDG method for Stokes",
                     AboutData::License_GPL,
                     "Copyright (c) 2016-2017 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Daniele Prada", "developer", "daniele.prada85@gmail.com", "" );
    about.addAuthor( "Lorenzo Sala", "developer", "sala@unistra.fr", "" );
    about.addAuthor( "Philippe Ricka", "developer", "pricka@math.unistra.fr", "" );

    return about;

}


template<int Dim, int OrderP, int OrderG=1>
int hdg_stokes( std::map<std::string,std::map<std::string,std::string>>& locals )
{

	typedef Simplex<Dim,OrderG> convex_type;
	//! mesh type
    typedef Mesh<convex_type> mesh_type;
    using Feel::cout;

    auto tau_constant =  cst(doption("hdg.tau.constant"));
    int tau_order =  ioption("hdg.tau.order");
    auto mu = locals.at("mu");
    auto stressn = expr<Dim,1>(locals.at("stressn"));
    auto u_d = expr<Dim,1>(locals.at("u"));
    auto f = expr<Dim,1>(locals.at("f"));

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
    auto face_mesh = createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );
    auto Mh = Pdhv<OrderP>( face_mesh,true );

    toc("spaces",true);

    cout << "Vh<" << OrderP   << "> : " << Vh->nDof() << std::endl
         << "Wh<" << OrderP   << "> : " << Wh->nDof() << std::endl
         << "Ph<" << OrderP   << "> : " << Ph->nDof() << std::endl
         << "Mh<" << OrderP   << "> : " << Mh->nDof() << std::endl;

    auto delta = Vh->element( "sigma" );
    auto gamma = Vh->element( "v" );
    auto u     = Wh->element( "u" );
    auto v     = Wh->element( "w" );
    auto p     = Ph->element( "p" );
    auto q     = Ph->element( "q" );
    auto uhat  = Mh->element( "uhat" );
    auto m     = Mh->element( "m" );

    // Number of dofs associated with each space
    auto nDofsigma = delta.functionSpace()->nDof();
    auto nDofu     = u.functionSpace()->nDof();
    auto nDofp     = p.functionSpace()->nDof();
    auto nDofuhat  = uhat.functionSpace()->nDof();


	bool condense = 0;
	double sc_param = 1;
    if( condense )
        sc_param = 0.5;
	

	tic();
    auto ps = product( Vh, Wh, Ph, Mh );
	auto a = blockform2( ps , backend() );
	auto rhs = blockform1( ps , backend() );


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

    tic();
    rhs(1_c) += integrate(_range=elements(mesh),
                          _expr=trans(f)*idt(v) );
    rhs(3_c) -= integrate(_range=markedfaces(mesh, "Neumann"),
                          _expr=inner(stressn,idt(m)) );
    rhs(3_c) += integrate(_range=markedfaces(mesh, "Dirichlet"),
                          _expr=trans(u_d)*idt(m) );
    toc("rhs",true);

    // Building the matrix
    tic();
    a( 0_c, 0_c ) +=  integrate(_range=elements(mesh),
                                _expr=inner(idt(delta),id(gamma)) );
    toc("a(0,0)", true);

    tic();
    a( 0_c, 1_c ) += integrate(_range=elements(mesh),_expr=(-2*trans(idt(u))*div(gamma)));
    toc("a(0,1)", true);
    tic();
    a( 0_c, 3_c) += integrate(_range=internalfaces(mesh),
                              _expr=2*(trans(idt(uhat))*leftface(normal(gamma))+
                                      trans(idt(uhat))*rightface(normal(gamma))) );
    a( 0_c, 3_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=2*trans(idt(uhat))*(normal(gamma)));
    toc("a(0,3)", true);
    tic();
    a( 1_c, 0_c) += integrate(_range=elements(mesh),
                              _expr=mu*trans(id(delta))*gradt(v));
    toc("a(1,0)", true);
    tic();
    a( 1_c, 2_c) += integrate(_range=elements(mesh),
                              _expr=-inner(id(p)*eye<Dim>(),gradt(v)) );
    toc("a(1,2)", true);
    tic();
    a( 2_c, 1_c) += integrate(_range=elements(mesh),
                              _expr=-trans(id(u))*gradt(q) );
    toc("a(2,1)", true);
    tic();
    a( 2_c, 3_c) += integrate(_range=internalfaces(mesh),
                              _expr=normal(uhat)*jumpt(q) );
    a( 2_c, 3_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=normal(uhat)*idt(q) );
    toc("a(2,3)", true);
    tic();
    a( 3_c, 0_c) -= integrate(_range=internalfaces(mesh),
                              _expr=mu*(leftfacet(normal(delta))*id(m)+
                                         rightfacet(normal(delta)*id(m))) );
    a( 3_c, 0_c) -= integrate(_range=boundaryfaces(mesh),
                              _expr=mu*trans(normal(delta))*idt(m) );
    toc("a(3,0)", true);
    tic();
    a( 3_c, 1_c) += integrate(_range=internalfaces(mesh),
                              _expr=-tau_constant*(trans(leftface(id(u)))*idt(m)+
                                                   trans(rightface(id(u)))*idt(m)) );
    a( 3_c, 1_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=-tau_constant*(trans(id(u))*idt(m)) );
    toc("a(3,1)", true);
    tic();
    a( 3_c, 2_c) += integrate(_range=internalfaces(mesh),
                              _expr=trans(jump(p))*idt(m) );
    a( 3_c, 2_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=trans(p*N())*idt(m) );
    toc("a(3,2)", true);
    tic();
    a( 3_c, 3_c) += integrate(_range=internalfaces(mesh),                  //
                              _expr=tau_constant*trans(id(uhat))*idt(m) ); // these terms should be counted twice because of the sum on all the elements
    a( 3_c, 3_c) += integrate(_range=boundaryfaces(mesh),                  // each face apprearing twice
                              _expr=tau_constant*trans(id(uhat))*idt(m) ); //

    tic();
    auto U = ps.element();
    auto Ue = ps.element();
    a.solve( _solution=U, _rhs=rhs, _rebuild=true, _condense=condense);
    toc("solve",true);
    cout << "[Hdg] solve done" << std::endl;

    auto sigmap = U(0_c);
    auto up = U(1_c);
    auto pp = U(2_c);
    auto uhatp = U(3_c);

    int status = 1;
#if 0
    if ( displ.count("exact" ) )
    {
        auto displ_exact = displ.at("exact");
        auto sigma_exact = locals.at("stress").at("exact");
        auto grad_displ_exact = locals.at("grad_displ").at("exact");
        Ue(0_c).on( _range=elements(mesh), _expr=expr<Dim,Dim>( sigma_exact ) );
        Ue(1_c).on( _range=elements(mesh), _expr=expr<Dim,1>( displ_exact ) );

        auto l2err_sigma = normL2( _range=elements(mesh), _expr=expr<Dim,Dim>(sigma_exact) - idv(sigmap) );
        toc("error");
                    
        // CHECKER
        // compute l2 and h1 norm of u-u_h where u=solution
        auto norms_displ = [&]( std::string const& solution ) ->std::map<std::string,double>
            {
			tic();
			double l2 = normL2( _range=elements(mesh), _expr=expr<Dim,1>(displ_exact) - idv(up) );
			toc("L2 error norm");

			tic();
			double h1 = normH1(_range=elements(mesh), _expr=idv(up)- expr<Dim,1>(displ_exact), _grad_expr=gradv(up)-expr<Dim,Dim>(grad_displ_exact) );
			toc("H1 error norm");

			return { { "L2", l2 } , {  "H1", h1 } };
            };

        status = checker(displ_exact).runOnce( norms_displ, rate::hp( mesh->hMax(), Wh->fe()->order() ) );
        v.on( _range=elements(mesh), _expr=expr<Dim,Dim>(sigma_exact) );
        w.on( _range=elements(mesh), _expr=expr<Dim,1>(displ_exact) );
    }
#endif

    tic();
    std::string exportName =  "hdg_stokes";
    std::string sigmaName = "stress";
    std::string sigma_exName = "stress-ex";
    std::string uName = "u";
    std::string u_exName = "u-ex";

    auto e = exporter( _mesh=mesh, _name=exportName );
    e->setMesh( mesh );
    e->add( sigmaName, sigmap );
    e->add( uName, up );
#if 0
    if ( displ.count("exact" ) )
    {
        e->add( sigma_exName, v );
        e->add( u_exName, w );
    }
#endif
    e->save();

    toc("export");

    return !status;
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
    std::map<std::string,std::map<std::string,std::string>> locals{{"dim",{{"dim",std::to_string(FEELPP_DIM)}}},
        {"u", {{ "exact", soption("u")}} },
        {"grad_u",{}},
        {"strain",{}},
        {"stress",{}},
        {"stressn",{}},
        {"f",{}} };
    Feel::pyexprFromFile( Environment::expand(soption("pyexpr.filename")), locals  );

    for( auto d: locals )
        std::cout << d.first << ":" << d.second << std::endl;
    
    int status = hdg_stokes<FEELPP_DIM,2>( locals );
    return !status;

}
