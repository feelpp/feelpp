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
namespace Feel {




inline
po::options_description
makeOptions()
{
    po::options_description testhdivoptions( "test qs_hdg_elasticity options" );
    testhdivoptions.add_options()
        ( "pyexpr.filename", po::value<std::string>(), "filename" )
        ( "hsize", po::value<double>()->default_value( 0.8 ), "mesh size" )
        ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ( "lambda", po::value<std::string>()->default_value( "1" ), "lambda" )
        ( "mu", po::value<std::string>()->default_value( "1" ), "mu" )
        ( "hface", po::value<int>()->default_value( 0 ), "hface" )
        ( "use_hypercube", po::value<bool>()->default_value( false ), "use hypercube or a given geometry" )
        ( "tau_constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "tau_order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( "hdg.tau.constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "hdg.tau.order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ;
    return testhdivoptions.add( Feel::feel_options() ).add( backend_options("sc"));
}

inline
AboutData
makeAbout()
{
    AboutData about( "qs_hdg_elasticity" ,
                     "qs_hdg_elasticity" ,
                     "0.1",
                     "Quickstart for HDG method for linear elasticity",
                     AboutData::License_GPL,
                     "Copyright (c) 2016-2017 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Daniele Prada", "developer", "daniele.prada85@gmail.com", "" );
    about.addAuthor( "Lorenzo Sala", "developer", "sala@unistra.fr", "" );

    return about;

}


template<int Dim, int OrderP, int OrderG=1>
int hdg_elasticity()
{

	typedef Simplex<Dim,OrderG> convex_type;
	//! mesh type
    typedef Mesh<convex_type> mesh_type;
	//! the exporter factory type
	typedef Exporter<mesh_type> export_type;
	//! the exporter factory (shared_ptr<> type)
	typedef boost::shared_ptr<export_type> export_ptrtype;

    using Feel::cout;

    auto tau_constant =  cst(doption("hdg.tau.constant"));
    int tau_order =  ioption("hdg.tau.order");

    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;
    
#if defined(FEELPP_HAS_SYMPY )
    // Exact solutions
    auto dict = Feel::pyexprFromFile( soption("pyexpr.filename"), {"displ", "grad_displ", "strain", "stress", "stressn", "f", "c1", "c2"}  );

    cout << "displ : " << dict.at("displ") << std::endl;
    cout << "grad_displ : " << dict.at("grad_displ") << std::endl;
    cout << "strain : " << dict.at("strain") << std::endl;
    cout << "stress : " << dict.at("stress") << std::endl;
    cout << "f : " << dict.at("f") << std::endl;
    cout << "c1 : " << dict.at("c1") << std::endl;
    cout << "c2 : " << dict.at("c2") << std::endl;

    std::string displ_exact_str = dict.at("displ");
    std::string grad_displ_exact_str = dict.at("grad_displ");
    std::string stress_exact_str = dict.at("stress");
    auto displ_exact = expr<Dim,1>( displ_exact_str );
    auto grad_displ_exact = expr<Dim,Dim>( grad_displ_exact_str );
    auto sigma_exact = expr<Dim,Dim>( stress_exact_str );
    auto stressn_exact = expr<Dim,1>( dict.at("stressn") );
    auto f_exact = expr<Dim,1>( dict.at("f") );
    auto c1 = expr( dict.at("c1") );
    auto c2 = expr( dict.at("c2") );
#else
    std::string	sigma_exact_str = soption("solution.sigma");
    std::string displ_exact_str = soption("solution.u");
    auto sigma_exact = expr<Dim,Dim>(sigma_exact_str);
    auto displ_exact = expr<Dim,1>(displ_exact_str);
    auto f_exact = expr<Dim,1>( soption( "functions.f") );
    auto lambda = expr(soption("lambda"));
    auto mu     = expr(soption("mu"));
	auto c1     = cst(0.5)/mu;
 	auto c2     = -lambda/(cst(2.) * mu * (cst(Dim)*lambda + cst(2.)*mu));
#endif
    tic();
    auto mesh = loadMesh( new Mesh<Simplex<Dim>> );
    toc("mesh",true);

    // ****** Hybrid-mixed formulation ******
    // We treat Vh, Wh, and Mh separately
    tic();

    auto Vh = Pdhms<OrderP>( mesh, true );
    auto Wh = Pdhv<OrderP>( mesh, true );
    auto face_mesh = createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );
    auto Mh = Pdhv<OrderP>( face_mesh,true );

    toc("spaces",true);

    cout << "Vh<" << OrderP   << "> : " << Vh->nDof() << std::endl
         << "Wh<" << OrderP+1 << "> : " << Wh->nDof() << std::endl
         << "Mh<" << OrderP   << "> : " << Mh->nDof() << std::endl;

    auto sigma = Vh->element( "sigma" );
    auto v     = Vh->element( "v" );
    auto u     = Wh->element( "u" );
    auto w     = Wh->element( "w" );
    auto uhat  = Mh->element( "uhat" );
    auto m     = Mh->element( "m" );

    // Number of dofs associated with each space
    auto nDofsigma = sigma.functionSpace()->nDof();
    auto nDofu     = u.functionSpace()->nDof();
    auto nDofuhat  = uhat.functionSpace()->nDof();


	bool condense = 1;
	double sc_param = 1;
    if( condense )
        sc_param = 0.5;
	

	tic();
    auto ps = product( Vh, Wh, Mh );
	auto a = blockform2( ps , backend() );
	auto rhs = blockform1( ps , backend() );


    // Building the RHS
    auto M0h = Pdh<0>( face_mesh,true );
    auto H     = M0h->element( "H" );
    if ( ioption("hface" ) == 0 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hMax()) );
    else if ( ioption("hface" ) == 1 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hMin()) );
    else if ( ioption("hface" ) == 2 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hAverage()) );
    else
        H.on( _range=elements(face_mesh), _expr=h() );

    tic();
    rhs(1_c) += integrate(_range=elements(mesh), _expr=trans(f_exact)*id(w));

    // in convergence test Neumann condition is given from the displacement and
    // constitutive law
    rhs(2_c) += integrate(_range=markedfaces(mesh,"Neumann"), _expr=trans(id(m))*(stressn_exact));
    rhs(2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"), _expr=trans(id(m))*displ_exact);
    toc("rhs",true);

    // Building the matrix
    tic();
    a( 0_c, 0_c ) +=  integrate(_range=elements(mesh),_expr=(c1*inner(idt(sigma),id(v))) );
    toc("a(0,0).1", true);
    tic();
    a( 0_c, 0_c ) += integrate(_range=elements(mesh),_expr=c2*tracet(sigma)*trace(v) );
    toc("a(0,0).2", true);

    tic();
    a( 0_c, 1_c ) += integrate(_range=elements(mesh),_expr=(trans(idt(u))*div(v)));
    toc("a(0,1)", true);
    tic();
    a( 0_c, 2_c) += integrate(_range=internalfaces(mesh),
                              _expr=-( trans(idt(uhat))*leftface(normal(v))+
                                       trans(idt(uhat))*rightface(normal(v))) );
    a( 0_c, 2_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=-trans(idt(uhat))*(normal(v)));
    toc("a(0,2)", true);
    tic();
    a( 1_c, 0_c) += integrate(_range=elements(mesh),
                              _expr=(trans(id(w))*divt(sigma)));
    toc("a(1,0)", true);
    tic();
    // begin dp: here we need to put the projection of u on the faces
    a( 1_c, 1_c) += integrate(_range=internalfaces(mesh),_expr=-tau_constant *
                              ( leftfacet( pow(idv(H),tau_order)*trans(idt(u)))*leftface(id(w)) +
                                rightfacet( pow(idv(H),tau_order)*trans(idt(u)))*rightface(id(w) )));

    a( 1_c, 1_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=-(tau_constant * pow(idv(H),tau_order)*trans(idt(u))*id(w)));
    toc("a(1,1)", true);
    tic();
    a( 1_c, 2_c) += integrate(_range=internalfaces(mesh), _expr=tau_constant *
                              ( trans(idt(uhat))*(leftface( pow(idv(H),tau_order)*id(w))+
                                                  rightface( pow(idv(H),tau_order)*id(w) ))));

    a( 1_c, 2_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=tau_constant * trans(idt(uhat)) * pow(idv(H),tau_order)*id(w) );
    toc("a(1,2)", true);
    tic();
    a( 2_c, 0_c) += integrate(_range=internalfaces(mesh),
                              _expr=( trans(id(m))*(leftfacet(normalt(sigma))+
                                                    rightfacet(normalt(sigma))) ) );
    toc("a(2,0)", true);
    tic();
    a( 2_c, 1_c) += integrate(_range=internalfaces(mesh),
                              _expr=-tau_constant * trans(id(m)) * (leftfacet( pow(idv(H),tau_order)*idt(u) )+
                                                                    rightfacet( pow(idv(H),tau_order)*idt(u) )));
    toc("a(2,1)", true);
    tic();
    a( 2_c, 2_c) += integrate(_range=internalfaces(mesh),
                              _expr=sc_param*tau_constant * trans(idt(uhat)) * id(m) * ( leftface( pow(idv(H),tau_order) )+
                                                                                          rightface( pow(idv(H),tau_order) )));

    a( 2_c, 2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                              _expr=trans(idt(uhat)) * id(m) );

    toc("a(2,2)", true);
    tic();
    a( 2_c, 0_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                              _expr=( trans(id(m))*(normalt(sigma)) ));

    a( 2_c, 1_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                              _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),tau_order)*idt(u) ) );

    a( 2_c, 2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                              _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),tau_order) ) );
    toc("a(2,{0,1,2})", true);
    toc("matrices",true);

    tic();
    auto U = ps.element();
    auto Ue = ps.element();
    a.solve( _solution=U, _rhs=rhs, _rebuild=true, _condense=condense);
    toc("solve",true);
    cout << "[Hdg] solve done" << std::endl;

    auto sigmap = U(0_c);
    auto up = U(1_c);
    auto uhatp = U(2_c);

    Ue(0_c).on( _range=elements(mesh), _expr=sigma_exact );
    Ue(1_c).on( _range=elements(mesh), _expr=displ_exact );

#if 0
        Feel::cout << "sigma exact: \t" << Ue(0_c) << std::endl;
        Feel::cout << "sigma: \t" << sigmap << std::endl;
        Feel::cout << "u exact: \t" << Ue(1_c) << std::endl;
        Feel::cout << "u: \t" << up << std::endl;
        Feel::cout << "uhat: \t" << uhatp << std::endl;
#endif


    // ****** Compute error ******
    tic();
    bool has_dirichlet = nelements(markedfaces(mesh,"Dirichlet"),true) >= 1;
    BOOST_ASSERT(has_dirichlet);

   	/*
    How does Feel++ handle BC on single components? Say, Dirichlet on u_x and
    Neumann on dot(sigma*n, e_y)?
    */

    auto l2err_sigma = normL2( _range=elements(mesh), _expr=sigma_exact - idv(sigmap) );
   	toc("error");

	// CHECKER
    // compute l2 and h1 norm of u-u_h where u=solution
    auto norms_displ = [&]( std::string const& solution ) ->std::map<std::string,double>
		{
			tic();
			double l2 = normL2( _range=elements(mesh), _expr=displ_exact - idv(up) );
			toc("L2 error norm");

			tic();
			double h1 = normH1(_range=elements(mesh), _expr=idv(up)- displ_exact, _grad_expr=gradv(up)-grad_displ_exact );
			toc("H1 error norm");

			return { { "L2", l2 } , {  "H1", h1 } };
		};

    int status = checker(displ_exact_str).runOnce( norms_displ, rate::hp( mesh->hMax(), Wh->fe()->order() ) );


    tic();
    std::string exportName =  "hdg_elasticity";
    std::string sigmaName = "stress";
    std::string sigma_exName = "stress-ex";
    std::string uName = "displacement";
    std::string u_exName = "displacement-ex";

    v.on( _range=elements(mesh), _expr=sigma_exact );
    w.on( _range=elements(mesh), _expr=displ_exact );
    auto e = exporter( _mesh=mesh, _name=exportName );
    e->setMesh( mesh );
    e->add( sigmaName, sigmap );
    e->add( uName, up );
    e->add( sigma_exName, v );
    e->add( u_exName, w );
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
                     _about=about(_name="qs_hdg_elasticity",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    // end::env[]
   
    int status = hdg_elasticity<FEELPP_DIM,2>();
    return !status;

}
