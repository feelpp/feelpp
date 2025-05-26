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
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/print.hpp>

#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/operators2.hpp>
#include <feel/feelpython/pyexpr.hpp>

#include <feel/feelpde/cg_laplacian.hpp>
#include <fmt/format.h>
#include <fmt/ranges.h>    // for fmt::join

namespace Feel {

// somewhere in your .cpp or in a header:
static const std::string zero_vec = [](){
    // default‐initialize an array of FEELPP_DIM zeros
    std::array<int,FEELPP_DIM> a{};  // all elements == 0
    // join with commas and wrap in braces:
    return fmt::format("{{{}}}", fmt::join(a, ","));
}();


inline
po::options_description
makeOptions()
{
    po::options_description hdgoptions( "HDG options" );
    hdgoptions.add_options()
        ( "k", po::value<std::string>()->default_value( "1" ), "diffusion coefficient" )
        ( "beta", po::value<std::string>()->default_value( zero_vec ), "advection coefficient" )
        ( "r_1", po::value<std::string>()->default_value( "1" ), "Robin lhs coefficient" )
        ( "r_2", po::value<std::string>()->default_value( "" ), "Robin rhs coefficient" )
        ( "pyexpr.filename", po::value<std::string>()->default_value( "${top_srcdir}/quickstart/laplacian.py" ), "python filename to execute" )
        ( "solution.p", po::value<std::string>()->default_value( "1" ), "solution p exact" )
        ( "solution.sympy.p", po::value<std::string>()->default_value( "1" ), "solution p exact (if we use sympy)" )
        ( "solution.u", po::value<std::string>()->default_value( zero_vec ), "solution u exact" )
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


    int tau_order =  ioption("hdg.tau.order");

    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;


#if defined(FEELPP_HAS_SYMPY)

    // 1) Start the interpreter
    //py::scoped_interpreter guard{}; 

    // 2) Import our API module
    auto sympy_api = py::module_::import("feelpp.sympy.api");
    auto get_coeffs = sympy_api.attr("get_coefficients");

    // 3) Build the inputs map
    std::map<std::string,std::string> inputs{
        {"dim",std::to_string(Dim)},
        {"k",soption("k")},
        {"p",soption("checker.solution")},
        {"r_1",soption("r_1")}, 
        {"r_2",soption("r_2")}, 
        {"beta","{1,1}"}};
    // 4) Call into Python: mode "adr" for advection–diffusion–reaction, or "laplacian", etc.
    //    All kwargs are passed as strings and sympified on the Python side.
    py::dict py_kwargs;
    for ( auto const& [k,v] : inputs )
        py_kwargs[k.c_str()] = v;

    // e.g. "adr" — change to "laplacian" or "wave" or your new "stokes"/"darcy"
    py::object result = get_coeffs(py::str("adr"), py_kwargs);

    // 5) Cast back to a C++ map<string,string>
    auto locals = result.cast<std::map<std::string,std::string>>();
    // 6) Extract exactly as before, but now it's coming straight from Python objects:
    auto p_exact_str = locals.at("p");
    auto u_exact_str = locals.at("u");
    auto f_str       = locals.at("f");
    auto un_str      = locals.at("un");
    auto g_str       = locals.at("g");
    auto r1_str      = locals.at("r_1");
    auto r2_str      = locals.at("r_2");

    std::cout << fmt::format("p_exact_str: {}\nu_exact_str: {}\nf_str: {}\nun_str: {}\ng_str: {}\nr1_str: {}\nr2_str: {}",
                     p_exact_str, u_exact_str, f_str, un_str, g_str, r1_str, r2_str) << std::endl;
    // 7) Convert to Feel++ expressions
    auto p_exact = expr(p_exact_str);
    auto u_exact = expr<FEELPP_DIM,1>(u_exact_str);
    auto k        = expr(locals.at("k"));
    auto lambda   = cst(1.)/k;
    auto un       = expr(un_str);
    auto f        = expr(f_str);
    auto g        = expr(g_str);
    auto r_1      = expr(r1_str);
    auto r_2      = expr(r2_str);
    auto beta     = expr<FEELPP_DIM,1>(locals.at("beta"));
#if 0    
    std::map<std::string,std::string> inputs{
            {"dim",std::to_string(Dim)},
            {"k",soption("k")},
            {"beta", soption("beta")},
            {"p",soption("checker.solution")},
            {"grad_p",""},
            {"u",""},
            {"un",""},
            {"f",""},
            {"g",""},
            {"r_1",soption("r_1")},
            {"r_2",soption("r_2")}
    };
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
    auto beta = expr<FEELPP_DIM,1>( locals.at("beta") );
    auto lambda = cst(1.)/k;
    auto un = expr( locals.at("un") );
    auto f = expr( locals.at("f") );
    auto g = expr( locals.at("g") );
    auto r_1 = expr( locals.at("r_1") );
    auto r_2 = expr( locals.at("r_2") );
#endif

    // 8) Print the coefficients
    //std::cout << fmt::format("p_exact: {}\nu_exact: {}\nk: {}\nbeta: {}\nlambda: {}\nun: {}\nf: {}\ng: {}\nr_1: {}\nr_2: {}",
    //                 p_exact_str, u_exact_str, k, beta, lambda, un, f, g, r_1, r_2) << std::endl;    
#else
    std::string p_exact_str = soption("solution.p");
    std::string u_exact_str = soption("solution.u");
    auto p_exact = expr(p_exact_str);
    auto u_exact = expr<Dim,1>(u_exact_str);
    auto k = expr(soption("k"));
    auto beta = expr<Dim,1>(soption("beta"));
    auto lambda = cst(1.)/k;
    auto un = trans(u_exact)*N();
    auto f = expr( soption( "functions.f") );
    auto g = p_exact;
    auto r_1 = cst(0.);
    auto r_2 = un;
#endif
    auto beta_n = trans(beta) * N();
    tic();
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<Dim>> );
    toc("mesh",true);

    // ****** Hybrid-mixed formulation ******
    // We treat Vh, Wh, and Mh separately
    tic();

    auto Vh = Pdhv<OrderP>( mesh );
    auto Wh = Pdh<OrderP>( mesh );
    auto face_mesh = createSubmesh( _mesh=mesh, _range=faces(mesh ), _update=0 );
    auto Mh = Pdh<OrderP>( face_mesh );

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
        auto u = cgLaplacian( _space=cgXh, _data=std::tuple{k,f,p_exact,un,r_1,r_2} );
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
    auto inflow = chi(beta_n < 0.); // inflow: trace comes from p_hat
    auto outflow = chi(beta_n >= 0.); // outflow: trace comes from p
    auto tau_D =  cst(doption("hdg.tau.constant"));
    auto tau_S = max(beta_n,0)+tau_D/h();
    a(0_c,1_c) += integrate(_range=elements(mesh),_expr=-(idt(p)*div(v)) + (trans(beta)*idt(p))*id(v) );
    a(0_c,1_c) += integrate(_range=internalfaces(mesh),
                            _expr=( leftfacet( outflow*idt(p) )*leftface(normal(v)) +
                                                    rightfacet( outflow*idt(p) )*rightface(normal(v)) ));
    a(0_c,1_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=outflow*idt(p)*normal(v));
    toc("a(0,1)",FLAGS_v>0);

    tic();
    a(0_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=( inflow*idt(phat)*(leftface(normal(v))+
                                               rightface(normal(v)))) );
    a(0_c,2_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=inflow*idt(phat)*(normal(v)));
    toc("a(0,2)",FLAGS_v>0);

    //
    // Second row a(1_c,:)
    //
    tic();
    
    a(1_c,0_c) += integrate(_range=elements(mesh),_expr=(grad(w)*idt(u)) );
    toc("a(1,0)",FLAGS_v>0);

    tic();
    a(1_c,1_c) += integrate(_range=internalfaces(mesh),
                            _expr=tau_S *
                            ( leftfacet( idt(p))*leftface(id(w)) +
                              rightfacet( idt(p))*rightface(id(w) )));
    a(1_c,1_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=(tau_S * id(w)*idt(p)));
    toc("a(1,1)",FLAGS_v>0);

    tic();


    // Use λ on inflow, u on outflow
    auto upwind_trace = inflow*idt(phat) + outflow*idt(p);

    a(1_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=-tau_S * idt(phat) *
                            ( leftface( id(w) )+
                              rightface( id(w) )));
    a(1_c,2_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=-tau_S * idt(phat) * id(w) );
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
                            _expr=tau_S * id(l) * ( leftfacet( idt(p) )+
                                                           rightfacet( idt(p) )));

    a(2_c,1_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=tau_S * id(l) * ( idt(p) ) );
    toc("a(2,1)",FLAGS_v>0);

    tic();
    a(2_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=-(1.-0.5*boption("sc.condense"))*tau_S * idt(phat) * id(l) );
    a(2_c,2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=-tau_S * idt(phat) * id(l)  );
    a(2_c,2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                            _expr=idt(phat) * id(l) );
    // Robin
    a( 2_c, 0_c ) += integrate(_range=markedfaces(mesh,"Robin"),
                               _expr=id(l)*normalt(u) );
    a( 2_c, 1_c ) += integrate(_range=markedfaces(mesh,"Robin"),
                               _expr=tau_S * id(l) * idt(p)  );
    a( 2_c, 2_c ) += integrate(_range=markedfaces(mesh,"Robin"),
                               _expr=-tau_S * idt(phat) * id(l) );
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
    auto Whp = Pdh<OrderP+1>( mesh );
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

    try
    {
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
    }
    catch( ... )
    {
        handleExceptions();
    }
    return 1;
}
