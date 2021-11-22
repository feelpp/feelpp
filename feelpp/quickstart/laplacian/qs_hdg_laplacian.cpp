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
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feeldiscr/check.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <tabulate/table.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/print.hpp>

#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/operators2.hpp>
#include <feel/feelpython/pyexpr.hpp>

#include <feel/feelpde/cg_laplacian.hpp>

namespace Feel {


inline
po::options_description
makeOptions()
{
    po::options_description hdgoptions( "HDG options" );
    hdgoptions.add_options()
        ( "k", po::value<std::string>()->default_value( "1" ), "diffusion coefficient" )
        ( "r_1", po::value<std::string>()->default_value( "1" ), "Robin lhs coefficient" )
        ( "r_2", po::value<std::string>()->default_value( "" ), "Robin rhs coefficient" )
        ( "pyexpr.filename", po::value<std::string>()->default_value( "${top_srcdir}/quickstart/laplacian.py" ), "python filename to execute" )
        ( "solution.p", po::value<std::string>()->default_value( "1" ), "solution p exact" )
        ( "solution.sympy.p", po::value<std::string>()->default_value( "1" ), "solution p exact (if we use sympy)" )        
#if (FEELPP_DIM==2)
        ( "solution.u", po::value<std::string>()->default_value( "{0,0}" ), "solution u exact" )
#else
        ( "solution.u", po::value<std::string>()->default_value( "{0,0,0}" ), "solution u exact" )
#endif
        ( "hdg.tau.constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "hdg.tau.order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( "solvecg", po::value<bool>()->default_value( false ), "solve corresponding problem with CG"  )
        ( "order", po::value<int>()->default_value( 1 ), "approximation order"  )
        ( "rhs_quad", po::value<int>()->default_value( 4 ), "quadrature order"  )
        ( "mass.quad", po::value<bool>()->default_value( true ), "use quadrature for mass matrix assembly"  )
        ( "sc.transpose", po::value<bool>()->default_value( false ), "transpose block in SC assembly"  )
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
                     "Copyright (c) 2017-2020 Feel++ Consortium" );
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
    
    
#if defined(FEELPP_HAS_SYMPY)

    std::map<std::string,std::string> inputs{{"dim",std::to_string(Dim)},{"k",soption("k")},{"p",soption("checker.solution")},{"grad_p",""}, {"u",""}, {"un",""}, {"f",""}, {"g",""}, {"r_1",soption("r_1")}, {"r_2",soption("r_2")}};
    // if we do not check the results with a manufactured solution,
    // the right hand side is given by functions.f otherwise it is computed by the python script
    auto thechecker = checker( _name= "L2/H1 convergence", 
                               _solution_key="p",
                               _gradient_key="grad_p",
                               _inputs=inputs
                               );
    auto locals = thechecker.runScript();

    std::string p_exact_str = locals.at("p");
    std::string u_exact_str = locals.at("u");
    auto p_exact = expr( p_exact_str );
    auto u_exact = expr<FEELPP_DIM,1>( u_exact_str );
    auto k = expr( locals.at("k") );
    auto lambda = cst(1.)/k;
    auto un = expr( locals.at("un") );
    auto f = expr( locals.at("f") );
    auto g = expr( locals.at("g") );
    auto r_1 = expr( locals.at("r_1") );
    auto r_2 = expr( locals.at("r_2") ); 
#else
    std::string p_exact_str = soption("solution.p");
    std::string u_exact_str = soption("solution.u");
    auto p_exact = expr(p_exact_str);
    auto u_exact = expr<Dim,1>(u_exact_str);
    auto k = expr(soption("k"));
    auto lambda = cst(1.)/k;
    auto un = trans(u_exact)*N();
    auto f = expr( soption( "functions.f") );
    auto g = p_exact;
    auto r_1 = cst(0.);
    auto r_2 = un;
#endif
    tic();
    auto mesh = loadMesh( new Mesh<Simplex<Dim>> );
    toc("mesh",true);

    // ****** Hybrid-mixed formulation ******
    // We treat Vh, Wh, and Mh separately
    tic();

    auto Vh = Pdhv<OrderP>( mesh, true );
    auto Wh = Pdh<OrderP>( mesh, true );
    auto face_mesh = createSubmesh( _mesh=mesh, _range=faces(mesh ), _update=0 );
    auto Mh = Pdh<OrderP>( face_mesh,true );

    toc("spaces",true);
    auto P0dh = Pdh<0>(mesh);
    auto Xh = Pdh<0>(face_mesh);
    auto uf = Xh->element(cst(1.));

    cout << "Exact potential if applicable: " << p_exact_str << "\n"
         << "Exact flux if applicable: " << u_exact_str << "\n";

    cout << "#elts: " << mesh->numGlobalElements() << std::endl
         << "#faces: " << mesh->numGlobalFaces() << std::endl
         << "#facesMh: " << face_mesh->numGlobalElements() << std::endl
         << "Vh<" << OrderP << "> : " << Vh->nDof() << std::endl
         << "Wh<" << OrderP << "> : " << Wh->nDof() << std::endl
         << "Mh<" << OrderP << "> : " << Mh->nDof() << std::endl;
    cout << mesh->numGlobalElements()  << " " << mesh->numGlobalFaces() << " "
         << Vh->nDof() << " " << Wh->nDof() << " " << Mh->nDof() << std::endl;

    int status_cg = 0;
    if ( boption( "solvecg" ) == true )
    {
        Feel::cout << "-- CG<" << OrderP+1 << "> starts ----------------------------------------------------------\n";
        auto cgXh = Pch<OrderP+1>(mesh);
        Feel::cout << "cgXh<" << OrderP+1 << "> : " << cgXh->nDof() << std::endl;
        auto u = cgLaplacian( cgXh, std::tuple{k,f,p_exact,un,r_1,r_2} );
#if defined(FEELPP_HAS_SYMPY)
        if ( u )        
            status_cg = check( checker( _name= "L2/H1 convergence cG", 
                                        _solution_key="p",
                                        _gradient_key="grad_p",
                                        _inputs=locals
                                       ), *u );
#endif
        Feel::cout << "-- CG<" << OrderP+1 << "> done ----------------------------------------------------------\n";
    }
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
    bool sctrans = boption( "sc.transpose" );
    auto a = blockform2( ps, strategy ,backend() );
    auto rhs = blockform1( ps, strategy, backend() );
    toc("forms",true);

    tic(); // total
    tic(); // assembly
    tic();
    // Building the RHS
    //
    // This is only a part of the RHS - how to build the whole RHS? Is it right to
    // imagine we moved it to the left? SKIPPING boundary conditions for the moment.
    // How to identify Dirichlet/Neumann boundaries?
    rhs(1_c) += integrate(_range=elements(mesh),
                          _expr=f*id(w), _quad=ioption("rhs_quad") );

    rhs(2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                          _expr=id(l)*un, _quad=ioption("rhs_quad")  );
    rhs(2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                          _expr=id(l)*p_exact, _quad=ioption("rhs_quad") );
    rhs(2_c) += integrate( _range=markedfaces(mesh, "Robin"),
                           _expr=id(l)*r_2);
    
    toc("rhs",true);
    tic();
    //
    // First row a(0_c,:)
    //
    tic();
    if ( boption( "mass.quad" ) )
        a(0_c,0_c) += integrate(_range=elements(mesh),_expr=(trans(lambda*idt(u))*id(v)) );
    else
        a(0_c,0_c) += integrate(_range=elements(mesh),_expr=mass(u,v) );

    toc("a(0,0)",FLAGS_v>0);

    tic();
    a(0_c,1_c) += integrate(_range=elements(mesh),_expr=-(idt(p)*div(v)));
    toc("a(0,1)",FLAGS_v>0);

    tic();
    a(0_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=( idt(phat)*(leftface(normal(v))+
                                               rightface(normal(v)))) );
    a(0_c,2_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=idt(phat)*(normal(v)));
    toc("a(0,2)",FLAGS_v>0);

    //
    // Second row a(1_c,:)
    //
    tic();
    a(1_c,0_c) += integrate(_range=elements(mesh),_expr=(id(w)*divt(u)));
    toc("a(1,0)",FLAGS_v>0);

    tic();
    a(1_c,1_c) += integrate(_range=internalfaces(mesh),
                            _expr=tau_constant *
                            ( leftfacet( idt(p))*leftface(id(w)) +
                              rightfacet( idt(p))*rightface(id(w) )));
    a(1_c,1_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=(tau_constant * id(w)*idt(p)));
    toc("a(1,1)",FLAGS_v>0);

    tic();
    a(1_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=-tau_constant * idt(phat) *
                            ( leftface( id(w) )+
                              rightface( id(w) )));
    a(1_c,2_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=-tau_constant * idt(phat) * id(w) );
    toc("a(1,2)",FLAGS_v>0);

    //
    // Third row a(2_c,:)
    //
    
    tic();
    a(2_c,0_c) += integrate(_range=internalfaces(mesh),
                            _expr=( id(l)*(leftfacet(normalt(u))+rightfacet(normalt(u))))
                            //_expr=( cst(2.)*(leftfacet(trans(idt(u))*N())+rightfacet(trans(idt(u))*N())) ),
                            );
    toc("a(2,0).1",FLAGS_v>0);

    tic();
    // BC
    a(2_c,0_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=( id(l)*(normalt(u))));
    toc("a(2,0).3",FLAGS_v>0);

    tic();
    a(2_c,1_c) += integrate(_range=internalfaces(mesh),
                            _expr=tau_constant * id(l) * ( leftfacet( idt(p) )+
                                                           rightfacet( idt(p) )));

    a(2_c,1_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=tau_constant * id(l) * ( idt(p) ) );
    toc("a(2,1)",FLAGS_v>0);

    tic();
    a(2_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=-(1.-0.5*boption("sc.condense"))*tau_constant * idt(phat) * id(l) );
    a(2_c,2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=-tau_constant * idt(phat) * id(l)  );
    a(2_c,2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                            _expr=idt(phat) * id(l) );
    // Robin
    a( 2_c, 0_c ) += integrate(_range=markedfaces(mesh,"Robin"),
                               _expr=id(l)*normalt(u) );
    a( 2_c, 1_c ) += integrate(_range=markedfaces(mesh,"Robin"),
                               _expr=tau_constant * id(l) * idt(p)  );
    a( 2_c, 2_c ) += integrate(_range=markedfaces(mesh,"Robin"),
                               _expr=-tau_constant * idt(phat) * id(l) );
    a( 2_c, 2_c ) += integrate(_range=markedfaces(mesh,"Robin"),
                               _expr=-r_1*idt(phat) * id(l) );

    toc("a(2,2)",FLAGS_v>0);


    toc("matrices",true);
    toc("assembly",true);
    
    tic(); // solver+postpro time
    tic();
    auto U=ps.element();
    a.solve( _solution=U, _rhs=rhs, _condense=boption("sc.condense"));
    toc("solve",true);

    
    // ****** Compute error ******
    auto up = U(0_c);
    auto pp = U(1_c);

    tic();
    tic();
    auto Whp = Pdh<OrderP+1>( mesh, true );
    auto pps = product( Whp );
    auto PP = pps.element();
    auto ppp = PP(0_c);
    toc("postproceSsing.space",FLAGS_v>0);
    tic();
    tic();
    auto b = blockform2( pps, solve::strategy::local, backend() );
    b( 0_c, 0_c ) = integrate( _range=elements(mesh), _expr=inner(gradt(ppp),grad(ppp)));
    toc("postprocessing.assembly.a",FLAGS_v>0);
    tic();
    auto ell = blockform1( pps, solve::strategy::local, backend() );
    ell(0_c) = integrate( _range=elements(mesh), _expr=-lambda*grad(ppp)*idv(up));
    toc("postprocessing.assembly.l",FLAGS_v>0);
    toc("postprocessing.assembly",FLAGS_v>0);
    
    tic();
    tic();
    b.solve( _solution=PP, _rhs=ell, _name="sc.post", _local=true);
    toc("postprocessing.solve.local",FLAGS_v>0);
    ppp=PP(0_c);
    tic();
    tic();
    ppp -= ppp.ewiseMean(P0dh);
    toc("postprocessing.solve.correction.ppp",FLAGS_v>0);
    tic();
    ppp += pp.ewiseMean(P0dh);
    toc("postprocessing.solve.correction.pp",FLAGS_v>0);
    toc("postprocessing.solve.correction",FLAGS_v>0);
    toc("postprocessing.solve");
    toc("postprocessing");

    toc("solver+postprocessing");
    toc("assembly+solver+postprocessing");
    

    tic();
    v.on( _range=elements(mesh), _expr=u_exact );
    q.on( _range=elements(mesh), _expr=p_exact );
    auto e = exporter( _mesh=mesh );
    e->setMesh( mesh );
    e->add( "flux", U(0_c) );
    e->add( "potential", U(1_c) );
    e->add( "potentialpp", PP(0_c) );
    e->add( "flux.exact", v );
    e->add( "potential.exact", q );
    e->save();
    toc("export");

    tic();

    int status1 = 0, status2 = 0, status3 = 0;
#if defined(FEELPP_HAS_SYMPY)
    bool has_dirichlet = nelements(markedfaces(mesh,"Dirichlet"),true) >= 1;
    solution_t s_t = has_dirichlet?solution_t::unique:solution_t::up_to_a_constant;
    status1 = check( checker( _name= "L2/H1 convergence of potential", 
                              _solution_key="p",
                              _gradient_key="grad_p",
                              _inputs=locals
                              ), pp, s_t );
    status2 = check( checker( _name= "L2 convergence of the flux", 
                              _solution_key="u",
                              _inputs=locals
                              ), up );
    status3 = check( checker( _name= "L2/H1 convergence of postprocessed potential", 
                              _solution_key="p",
                              _gradient_key="grad_p",
                              _inputs=locals
                              ), ppp, s_t );    
    // end::check[]
#endif

    return status_cg || status1 || status2 || status3;
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
    if ( ioption( "order" ) == 1 )
        return hdg_laplacian<FEELPP_DIM,1>();
    if ( ioption( "order" ) == 2 )
        return hdg_laplacian<FEELPP_DIM,2>();

 #if 0   
    if ( ioption( "order" ) == 3 )
        return hdg_laplacian<FEELPP_DIM,3>();

    if ( ioption( "order" ) == 4 )
        return hdg_laplacian<FEELPP_DIM,4>();
#endif
    return 1;
}
