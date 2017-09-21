//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 26 Aug 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <feel/feel.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/time.hpp>
#include <feel/feelpython/pyexpr.hpp>
namespace Feel {


inline
po::options_description
makeOptions()
{
    po::options_description hdgoptions( "HDG options" );
    hdgoptions.add_options()
        ( "k", po::value<std::string>()->default_value( "-1" ), "diffusion coefficient" )
        ( "solution.p", po::value<std::string>()->default_value( "1" ), "solution p exact" )
        ( "hdg.tau.constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "hdg.tau.order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ;
    return hdgoptions;
}

inline
AboutData
makeAbout()
{
    AboutData about( "qs_hdg_laplacian" ,
                     "qs_hdg_laplacian" ,
                     "0.1",
                     "Quickstart HDG Laplacian",
                     AboutData::License_GPL,
                     "Copyright (c) 2017 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

template<int Dim, int OrderP>
int hdg_laplacian()
{
    using Feel::cout;

    auto tau_constant =  cst(doption("hdg.tau.constant"));
    int tau_order =  ioption("hdg.tau.order");

    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;

    auto K = expr(soption("k"));
    auto lambda = cst(1.)/K;
    

    // Exact solutions
    std::ostringstream ostr;
    ostr <<  "from sympy2ginac import *\n"
         << "s=syms(" << Dim << ");\n"
         << "ns=nsyms(" << Dim << ");\n"
         << "p=sympify("<< soption("solution.p") << ");\n"
         << "k=sympify("<< soption("k") << ");\n"
         << "grad_p=grad(p,s);\n"
         << "flux=-k*grad(p,s);\n"
         << "u=flux;\n"
         << "un=n(flux,1,ns);\n"
         << "f=div(flux,s);\n";
    cout << ostr.str() << std::endl;

    auto m = Feel::pyexpr( ostr.str().c_str(), {"p", "grad_p", "u", "un", "f"} );

    cout << "k : " << K /*<< "\tlambda : " << lambda*/ << std::endl;
    cout << "p : " << m.at("p") << std::endl;
    cout << "grad_p : " << m.at("grad_p") << std::endl;
    cout << "u : " << m.at("u") << std::endl;
    cout << "un : " << m.at("un") << std::endl;
    cout << "f : " << m.at("f") << std::endl;

    auto p_exact = expr(m.at("p"));
    auto u_exact = expr<Dim,1>( m.at("u") );
    auto un_exact = expr( m.at("un") );
    auto f_exact = expr( m.at("f") );

    tic();
    auto mesh = loadMesh( new Mesh<Simplex<Dim>> );
    toc("mesh",true);

    // ****** Hybrid-mixed formulation ******
    // We treat Vh, Wh, and Mh separately
    tic();

    auto Vh = Pdhv<OrderP>( mesh, true );
    auto Wh = Pdh<OrderP>( mesh, true );
    auto face_mesh = createSubmesh( mesh, faces(mesh ), EXTRACTION_KEEP_MESH_RELATION, 0 );
    auto Mh = Pdh<OrderP>( face_mesh,true );

    toc("spaces",true);

    auto Xh = Pdh<0>(face_mesh);
    auto uf = Xh->element(cst(1.));

    cout << "Vh<" << OrderP << "> : " << Vh->nDof() << std::endl
         << "Wh<" << OrderP << "> : " << Wh->nDof() << std::endl
         << "Mh<" << OrderP << "> : " << Mh->nDof() << std::endl;

    auto u = Vh->element( "u" );
    auto v = Vh->element( "v" );
    auto p = Wh->element( "p" );
    auto q = Wh->element( "q" );
    auto w = Wh->element( "w" );
    auto phat = Mh->element( "phat" );
    auto l = Mh->element( "lambda" );

    tic();
    auto ps = product( Vh, Wh, Mh );
    toc("space",true);
    tic();
    solve::strategy strategy = boption("sc.condense")?solve::strategy::static_condensation:solve::strategy::monolithic;
    auto a = blockform2( ps, strategy ,backend() );
    auto rhs = blockform1( ps, strategy, backend() );
    toc("forms",true);

    tic();
    // Building the RHS
    //
    // This is only a part of the RHS - how to build the whole RHS? Is it right to
    // imagine we moved it to the left? SKIPPING boundary conditions for the moment.
    // How to identify Dirichlet/Neumann boundaries?
    rhs(1_c) += integrate(_range=elements(mesh),
                          _expr=-f_exact*id(w));

    rhs(2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                          _expr=id(l)*un_exact );
    rhs(2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                          _expr=id(l)*p_exact);

    toc("rhs",true);
    tic();
    //
    // First row a(0_c,:)
    //
    tic();
    a(0_c,0_c) += integrate(_range=elements(mesh),_expr=(trans(lambda*idt(u))*id(v)) );
    toc("a(0,0)",FLAGS_v>0);

    tic();
    a(0_c,1_c) += integrate(_range=elements(mesh),_expr=-(idt(p)*div(v)));
    toc("a(0,1)",FLAGS_v>0);

    tic();
    a(0_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=( idt(phat)*(leftface(trans(id(v))*N())+
                                               rightface(trans(id(v))*N()))), _verbose=true );

    a(0_c,2_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=idt(phat)*(trans(id(v))*N()));
    toc("a(0,2)",FLAGS_v>0);

    //
    // Second row a(1_c,:)
    //
    tic();
    a(1_c,0_c) += integrate(_range=elements(mesh),_expr=-(id(w)*divt(u)));
    toc("a(1,0)",FLAGS_v>0);

    tic();
    a(1_c,1_c) += integrate(_range=internalfaces(mesh),
                            _expr=-tau_constant *
                            ( leftfacet( pow(h(),tau_order)*idt(p))*leftface(id(w)) +
                              rightfacet( pow(h(),tau_order)*idt(p))*rightface(id(w) )));
    a(1_c,1_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=-(tau_constant * pow(h(),tau_order)*id(w)*idt(p)));
    toc("a(1,1)",FLAGS_v>0);

    tic();
    a(1_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=tau_constant * idt(phat) *
                            ( leftface( pow(h(),tau_order)*id(w) )+
                              rightface( pow(h(),tau_order)*id(w) )));

    a(1_c,2_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=tau_constant * idt(phat) * pow(h(),tau_order)*id(w) );
    toc("a(1,2)",FLAGS_v>0);

    //
    // Third row a(2_c,:)
    // 
    tic();
    a(2_c,0_c) += integrate(_range=internalfaces(mesh),
                            _expr=( id(l)*(leftfacet(trans(idt(u))*N())+rightfacet(trans(idt(u))*N()))),
                            //_expr=( cst(2.)*(leftfacet(trans(idt(u))*N())+rightfacet(trans(idt(u))*N())) ),
                            _verbose=true);
        
    toc("a(2,0).1",FLAGS_v>0);
        
    tic();
    // BC
    a(2_c,0_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=( id(l)*(trans(idt(u))*N())));
    toc("a(2,0).3",FLAGS_v>0);

    tic();
    a(2_c,1_c) += integrate(_range=internalfaces(mesh),
                            _expr=tau_constant * id(l) * ( leftfacet( pow(h(),tau_order)*idt(p) )+
                                                           rightfacet( pow(h(),tau_order)*idt(p) )),_verbose=true);
    a(2_c,1_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=tau_constant * id(l) * ( pow(h(),tau_order)*idt(p) ),_verbose=true );
    toc("a(2,1)",FLAGS_v>0);

    tic();
    a(2_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=-0.5*tau_constant * idt(phat) * id(l) * ( leftface( pow(h(),tau_order) )+
                                                                            rightface( pow(h(),tau_order) )));
    a(2_c,2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=-tau_constant * idt(phat) * id(l) * ( pow(h(),tau_order) ) );
    a(2_c,2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                            _expr=idt(phat) * id(l) );
    toc("a(2,2)",FLAGS_v>0);

    toc("matrices",true);

    tic();
    auto U=ps.element();
    a.solve( _solution=U, _rhs=rhs, _condense=boption("sc.condense"));
    toc("solve",true);

    // ****** Compute error ******
    auto up = U(0_c);
    auto pp = U(1_c);



    tic();
    v.on( _range=elements(mesh), _expr=u_exact );
    q.on( _range=elements(mesh), _expr=p_exact );
    auto e = exporter( _mesh=mesh );
    e->setMesh( mesh );
    e->add( "flux", U(0_c) );
    e->add( "potential", U(1_c) );
    e->add( "flux.exact", v );
    e->add( "potential.exact", q );
    e->save();
    toc("export");

    tic();

    bool has_dirichlet = nelements(markedfaces(mesh,"Dirichlet"),true) >= 1;


    // tag::check[]
    // compute l2 and h1 norm of u-u_h where u=solution
    auto norms_u = [=]( std::string const& solution ) ->std::map<std::string,double>
        {
            tic();
            auto l2err_u = normL2( _range=elements(mesh), _expr=u_exact - idv(up) );
            toc("L2 error norm u");
            return { { "L2", l2err_u }};
        };

    auto norms_p = [=]( std::string const& solution ) ->std::map<std::string,double>
        {
            tic();
            double l2err_p = 1e+30;
            if ( has_dirichlet )
            {
                l2err_p = normL2( _range=elements(mesh), _expr=expr(solution) - idv(pp) );
            }
            else
            {
                auto mean_p_exact = mean( elements(mesh), expr(solution) )(0,0);
                auto mean_p = mean( elements(mesh), idv(pp) )(0,0);
                l2err_p = normL2( elements(mesh),
                                  (expr(solution) - cst(mean_p_exact)) - (idv(pp) - cst(mean_p)) );
            }
            toc("L2 error norm p");
            return { { "L2", l2err_p }};
        };
#if 0
    int status = checker().runOnce( {norms_p,norms_u},
                                    { rate::hp( mesh->hMax(), Vh->fe()->order() ), rate::hp( mesh->hMax(), Vh->fe()->order() ) } );
#else
    int status1 = checker().runOnce( norms_p, rate::hp( mesh->hMax(), Vh->fe()->order() ) );
    int status2 = checker().runOnce( norms_u, rate::hp( mesh->hMax(), Vh->fe()->order() ) );
#endif
    // end::check[]

    return status1 || status2;
}



} // Feel

int main( int argc, char** argv )
{
    // tag::env[]
    using namespace Feel;

	Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="qs_hdg_laplacian",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    // end::env[]
   
    int status = hdg_laplacian<FEELPP_DIM,2>();
    return !status;
}
