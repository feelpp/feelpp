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
        ( "hsize", po::value<double>()->default_value( 0.8 ), "mesh size" )
        ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ( "lambda", po::value<std::string>()->default_value( "1" ), "lambda" )
        ( "mu", po::value<std::string>()->default_value( "1" ), "mu" )
        ( "u_exact", po::value<std::string>()->default_value("empty"), "u exact" )
		( "f", po::value<std::string>()->default_value("empty"), "divergence of the stress tensor")
        ( "hface", po::value<int>()->default_value( 0 ), "hface" )
        ( "use_hypercube", po::value<bool>()->default_value( false ), "use hypercube or a given geometry" )
        ( "tau_constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "tau_order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
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

        using Feel::cout;

    auto tau_constant =  cst(doption("hdg.tau.constant"));
    int tau_order =  ioption("hdg.tau.order");

    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;

    auto K = expr(soption("k"));
    auto lambda = cst(1.)/K;
    
#if defined(FEELPP_HAS_SYMPY )
    // Exact solutions
    std::ostringstream ostr;
    ostr <<  "from sympy2ginac import *\n"
         << "s=syms(" << Dim << ");\n"
         << "ns=nsyms(" << Dim << ");\n"
         << "p=sympify("<< soption("solution.sympy.p") << ");\n"
         << "k=sympify("<< soption("k") << ");\n"
         << "grad_p=grad(p,s);\n"
         << "flux=-k*grad(p,s);\n"
         << "u=flux;\n"
         << "un=n(flux,1,ns);\n"
         << "f=div(flux,s);\n";
    cout << ostr.str() << std::endl;

    auto dict = Feel::pyexpr( ostr.str().c_str(), {"p", "grad_p", "u", "un", "f"} );

    cout << "k : " << K /*<< "\tlambda : " << lambda*/ << std::endl;
    cout << "p : " << dict.at("p") << std::endl;
    cout << "grad_p : " << dict.at("grad_p") << std::endl;
    cout << "u : " << dict.at("u") << std::endl;
    cout << "un : " << dict.at("un") << std::endl;
    cout << "f : " << dict.at("f") << std::endl;

    std::string p_exact_str = dict.at("p");
    std::string u_exact_str = dict.at("u");
    auto p_exact = expr( p_exact_str );
    auto u_exact = expr<Dim,1>( u_exact_str );
    auto un_exact = expr( dict.at("un") );
    auto f_exact = expr( dict.at("f") );
#else
    std::string p_exact_str = soption("solution.p");
    std::string u_exact_str = soption("solution.u");
    auto p_exact = expr(p_exact_str);
    auto u_exact = expr<Dim,1>(u_exact_str);
    auto un_exact = trans(u_exact)*N();
    auto f_exact = expr( soption( "functions.f") );
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



	tic();
    auto ps = product( Vh, Wh, Mh );
	auto a = blockform2( ps , backend() );
	auto rhs = blockform1( ps , backend() );


    // Building the RHS
    M0h_ptr_t M0h = Pdh<0>( face_mesh,true );
    auto H     = M0h->element( "H" );
    if ( ioption("hface" ) == 0 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hMax()) );
    else if ( ioption("hface" ) == 1 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hMin()) );
    else if ( ioption("hface" ) == 2 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hAverage()) );
    else
        H.on( _range=elements(face_mesh), _expr=h() );

    rhs(1_c) += integrate(_range=elements(mesh), _expr=trans(f)*id(w));
    // rhs(1_c) += integrate(_range=elements(mesh), _expr= divv(sigma_exact)*id(w));

    cout << "rhs2 works fine" << std::endl;

    // in convergence test Neumann condition is given from the displacement and
    // constitutive law
    rhs(2_c) += integrate(_range=markedfaces(mesh,"Neumann"), _expr=trans(id(m))*(sigma_exact*N()));

    rhs(2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                          _expr=trans(id(m))*u_exact);


    cout << "rhs3 works fine" << std::endl;

    // Building the matrix
    a( 0_c, 0_c ) +=  integrate(_range=elements(mesh),_expr=(c1*inner(idt(sigma),id(v))) );
    a( 0_c, 0_c ) += integrate(_range=elements(mesh),_expr=(c2*trace(idt(sigma))*trace(id(v))) );

    a( 0_c, 1_c ) += integrate(_range=elements(mesh),_expr=(trans(idt(u))*div(v)));

    a( 0_c, 2_c) += integrate(_range=internalfaces(mesh),
                              _expr=-( trans(idt(uhat))*leftface(id(v)*N())+
                                       trans(idt(uhat))*rightface(id(v)*N())) );
    a( 0_c, 2_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=-trans(idt(uhat))*(id(v)*N()));

    a( 1_c, 0_c) += integrate(_range=elements(mesh),
                              _expr=(trans(id(w))*divt(sigma)));
    // begin dp: here we need to put the projection of u on the faces
    a( 1_c, 1_c) += integrate(_range=internalfaces(mesh),_expr=-tau_constant *
                              ( leftfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*leftface(id(w)) +
                                rightfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*rightface(id(w) )));

    a( 1_c, 1_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=-(tau_constant * pow(idv(H),M_tau_order)*trans(idt(u))*id(w)));

    a( 1_c, 2_c) += integrate(_range=internalfaces(mesh), _expr=tau_constant *
                              ( leftfacet(trans(idt(uhat)))*leftface( pow(idv(H),M_tau_order)*id(w))+
                                rightfacet(trans(idt(uhat)))*rightface( pow(idv(H),M_tau_order)*id(w) )));

    a( 1_c, 2_c) += integrate(_range=boundaryfaces(mesh),
                              _expr=tau_constant * trans(idt(uhat)) * pow(idv(H),M_tau_order)*id(w) );


    a( 2_c, 0_c) += integrate(_range=internalfaces(mesh),
                              _expr=( trans(id(m))*(leftfacet(idt(sigma)*N())+
                                                    rightfacet(idt(sigma)*N())) ) );
    a( 2_c, 1_c) += integrate(_range=internalfaces(mesh),
                              _expr=-tau_constant * trans(id(m)) * (leftfacet( pow(idv(H),M_tau_order)*idt(u) )+
                                                                    rightfacet( pow(idv(H),M_tau_order)*idt(u) )));

    a( 2_c, 2_c) += integrate(_range=internalfaces(mesh),
                              _expr=sc_param*tau_constant * trans(idt(uhat)) * id(m) * ( leftface( pow(idv(H),M_tau_order) )+
                                                                                          rightface( pow(idv(H),M_tau_order) )));

    a( 2_c, 2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                              _expr=trans(idt(uhat)) * id(m) );


    a( 2_c, 0_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                              _expr=( trans(id(m))*(idt(sigma)*N()) ));

    a( 2_c, 1_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                              _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tau_order)*idt(u) ) );

    a( 2_c, 2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                              _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tau_order) ) );

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
    Ue(1_c).on( _range=elements(mesh), _expr=u_exact );

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
  	auto l2err_u     = normL2( _range=elements(mesh), _expr=u_exact - idv(up) );
	auto h1err_u 	 = normH1(_range=elements(mesh), _expr=idv(up)- u_exact, _grad_expr=gradv(up)-grad(u_exact) );
   	toc("error");

	if ( !checker().check() )
	{
    	cout << "||sigma_exact - sigma||_L2 = " << l2err_sigma << std::endl;
    	cout << "||u_exact - u||_L2 = " << l2err_u << std::endl;
	}

	// CHECKER
	if ( checker().check() )
	{

		// compute l2 and h1 norm of u-u_h where u=solution
		auto norms = [=]( std::string const& u_exact ) ->std::map<std::string,double>
		{
			tic();
			double l2 = l2err_u; 
			toc("L2 error norm");
			
			tic();
			double h1 = h1err_u; 
			toc("H1 error norm");
			
			return { { "L2", l2 } , {  "H1", h1 } };
		};
		
		status = checker().runOnce( norms, rate::hp( mesh->hMax(), Wh->fe()->order() ) );
		
	}

    tic();
    std::string exportName =  ( boost::format( "%1%" ) % this->about().appName() ).str();
    std::string sigmaName = "stress";
    std::string sigma_exName = "stress-ex";
    std::string uName = "displacement";
    std::string u_exName = "displacement-ex";

    v.on( _range=elements(mesh), _expr=sigma_exact , _quad=_Q<expr_order>());
    w.on( _range=elements(mesh), _expr=u_exact , _quad=_Q<expr_order>());
    export_ptrtype exporter_cvg( export_type::New( exportName ) );

    exporter_cvg->setMesh( mesh );
    exporter_cvg->add( sigmaName, sigmap );
    exporter_cvg->add( uName, up );
    exporter_cvg->add( sigma_exName, v );
    exporter_cvg->add( u_exName, w );
    exporter_cvg->save();

    toc("export");

	// return !status;
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
