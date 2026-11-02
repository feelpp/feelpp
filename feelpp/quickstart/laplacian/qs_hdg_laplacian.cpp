/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2017-2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
//! @file
//! @date 26 Aug 2017
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
#include <feel/feelvf/vf_eval.hpp>
#include <feel/feelvf/print.hpp>

#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/operators2.hpp>
#include <feel/feelpython/pyexpr.hpp>

#include <feel/feelpde/cg_laplacian.hpp>
#include "../discretization.hpp"
#include "../local_mean.hpp"
#include <algorithm>

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
        ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" )
        ( "rhs_quad", po::value<int>()->default_value( 0 ), "RHS quadrature order (0: expression order with a 2*k floor, scaled by geometry order); increase explicitly for nonpolynomial data" )
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

template<int Dim>
int hdg_laplacian()
{
    using Feel::cout;

    auto tau_constant =  cst(doption("hdg.tau.constant"));
    int tau_order =  ioption("hdg.tau.order");
    int quadratureOrder = ioption( "rhs_quad" );
    CHECK_GE( quadratureOrder, 0 ) << "rhs_quad must be non-negative";
    if ( quadratureOrder == 0 )
        quadratureOrder = quad_order_from_expression;
    int order = ioption( "order" );
    CHECK( order >= 1 ) << "HDG approximation order must be >= 1";
    auto runtimeOrder = RuntimeOrder::checked( order );
    auto runtimePostOrder = RuntimeOrder::checked( order + 1 );
    int geometryOrder = ioption( "geometry.order" );
    CHECK( geometryOrder >= 1 ) << "Mesh geometry order must be >= 1";
    auto runtimeGeometryOrder = RuntimeOrder::checked( geometryOrder );

    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;


#if defined(FEELPP_HAS_SYMPY)

    std::string sympyPotential = soption( "solution.sympy.p" );
    if ( sympyPotential == "1" && soption( "solution.p" ) != "1" )
        sympyPotential = soption( "checker.solution" );

    std::map<std::string,std::string> inputs{{"dim",std::to_string(Dim)},{"k",soption("k")},{"p",sympyPotential},{"grad_p",""}, {"u",""}, {"un",""}, {"f",""}, {"g",""}, {"r_1",soption("r_1")}, {"r_2",soption("r_2")}};
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
    auto un = trans(u_exact)*N();
    auto f = expr( locals.at("f") );
    auto g = expr( locals.at("g") );
    auto r_1 = expr( locals.at("r_1") );
    auto r_2 = un - r_1*p_exact;
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
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<Dim,Dynamic>>( runtimeGeometryOrder ) );
    quickstart::checkMeshOrder( mesh );
    toc("mesh",true);

    // Match at least the degree of the HDG mass terms. On curved faces,
    // using a lower rule for representable boundary data can spoil exactness.
    // This is a polynomial estimate; nonpolynomial data may need rhs_quad.
    auto integrateRhs = [&]( auto const& range, auto const& expression )
    {
        int const q = quadratureOrder == quad_order_from_expression
                          ? std::max( static_cast<int>( expression.polynomialOrder() ), 2*order ) * mesh->order()
                          : quadratureOrder;
        return integrate( _range=range, _expr=expression, _quad=q );
    };

    // ****** Hybrid-mixed formulation ******
    // We treat Vh, Wh, and Mh separately
    tic();

    auto Vh = Pdhv<Dynamic>( mesh, runtimeOrder );
    auto Wh = Pdh<Dynamic>( mesh, runtimeOrder );
    auto face_mesh = createSubmesh( _mesh=mesh, _range=faces(mesh ), _update=0 );
    auto Mh = Pdh<Dynamic>( face_mesh, runtimeOrder );

    toc("spaces",true);
    auto P0dh = Pdh<Dynamic>( mesh, RuntimeOrder{0} );
    auto Xh = Pdh<Dynamic>( face_mesh, RuntimeOrder{0} );
    auto uf = Xh->element(cst(1.));

    cout << "Exact potential if applicable: " << p_exact_str << "\n"
         << "Exact flux if applicable: " << u_exact_str << "\n";

    cout << "#elts: " << mesh->numGlobalElements() << std::endl
         << "#faces: " << mesh->numGlobalFaces() << std::endl
         << "#facesMh: " << face_mesh->numGlobalElements() << std::endl
         << "Vh<" << order << "> : " << Vh->nDof() << std::endl
         << "Wh<" << order << "> : " << Wh->nDof() << std::endl
         << "Mh<" << order << "> : " << Mh->nDof() << std::endl;
    cout << mesh->numGlobalElements()  << " " << mesh->numGlobalFaces() << " "
         << Vh->nDof() << " " << Wh->nDof() << " " << Mh->nDof() << std::endl;

    int status_cg = 0;
    if ( boption( "solvecg" ) == true )
    {
        Feel::cout << "-- CG<" << order+1 << "> starts ----------------------------------------------------------\n";
        auto cgXh = Pch<Dynamic>( mesh, runtimePostOrder );
        Feel::cout << "cgXh<" << order+1 << "> : " << cgXh->nDof() << std::endl;
        auto u = cgLaplacian( _space=cgXh, _data=std::tuple{k,f,p_exact,un,r_1,r_2} );
#if defined(FEELPP_HAS_SYMPY)
        if ( u )
            status_cg = check( checker( _name= "L2/H1 convergence cG",
                                        _solution_key="p",
                                        _gradient_key="grad_p",
                                        _inputs=locals
                                       ), *u );
#endif
        Feel::cout << "-- CG<" << order+1 << "> done ----------------------------------------------------------\n";
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
    rhs(1_c) += integrateRhs( elements(mesh), f*id(w) );

    rhs(2_c) += integrateRhs( markedfaces(mesh,"Neumann"), id(l)*un  );
    rhs(2_c) += integrateRhs( markedfaces(mesh,"Dirichlet"), id(l)*p_exact );
    rhs(2_c) += integrateRhs( markedfaces(mesh, "Robin"), id(l)*r_2);

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

    toc("a(0,0)",Environment::logVerbosityLevel()>0);

    tic();
    a(0_c,1_c) += integrate(_range=elements(mesh),_expr=-(idt(p)*div(v)));
    toc("a(0,1)",Environment::logVerbosityLevel()>0);

    tic();
    a(0_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=( idt(phat)*(leftface(normal(v))+
                                               rightface(normal(v)))) );
    a(0_c,2_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=idt(phat)*(normal(v)));
    toc("a(0,2)",Environment::logVerbosityLevel()>0);

    //
    // Second row a(1_c,:)
    //
    tic();
    a(1_c,0_c) += integrate(_range=elements(mesh),_expr=(id(w)*divt(u)));
    toc("a(1,0)",Environment::logVerbosityLevel()>0);

    tic();
    a(1_c,1_c) += integrate(_range=internalfaces(mesh),
                            _expr=tau_constant *
                            ( leftfacet( idt(p))*leftface(id(w)) +
                              rightfacet( idt(p))*rightface(id(w) )));
    a(1_c,1_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=(tau_constant * id(w)*idt(p)));
    toc("a(1,1)",Environment::logVerbosityLevel()>0);

    tic();
    a(1_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=-tau_constant * idt(phat) *
                            ( leftface( id(w) )+
                              rightface( id(w) )));
    a(1_c,2_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=-tau_constant * idt(phat) * id(w) );
    toc("a(1,2)",Environment::logVerbosityLevel()>0);

    //
    // Third row a(2_c,:)
    //

    tic();
    a(2_c,0_c) += integrate(_range=internalfaces(mesh),
                            _expr=( id(l)*(leftfacet(normalt(u))+rightfacet(normalt(u))))
                            //_expr=( cst(2.)*(leftfacet(trans(idt(u))*N())+rightfacet(trans(idt(u))*N())) ),
                            );
    toc("a(2,0).1",Environment::logVerbosityLevel()>0);

    tic();
    // BC
    a(2_c,0_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=( id(l)*(normalt(u))));
    toc("a(2,0).3",Environment::logVerbosityLevel()>0);

    tic();
    a(2_c,1_c) += integrate(_range=internalfaces(mesh),
                            _expr=tau_constant * id(l) * ( leftfacet( idt(p) )+
                                                           rightfacet( idt(p) )));

    a(2_c,1_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=tau_constant * id(l) * ( idt(p) ) );
    toc("a(2,1)",Environment::logVerbosityLevel()>0);

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

    toc("a(2,2)",Environment::logVerbosityLevel()>0);


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
    auto Whp = Pdh<Dynamic>( mesh, runtimePostOrder );
    auto pps = product( Whp );
    auto PP = pps.element();
    auto ppp = PP(0_c);
    toc("postproceSsing.space",Environment::logVerbosityLevel()>0);
    tic();
    tic();
    auto b = blockform2( pps, solve::strategy::local, backend() );
    b( 0_c, 0_c ) = integrate( _range=elements(mesh), _expr=inner(gradt(ppp),grad(ppp)));
    toc("postprocessing.assembly.a",Environment::logVerbosityLevel()>0);
    tic();
    auto ell = blockform1( pps, solve::strategy::local, backend() );
    ell(0_c) = integrate( _range=elements(mesh), _expr=-lambda*grad(ppp)*idv(up));
    toc("postprocessing.assembly.l",Environment::logVerbosityLevel()>0);
    toc("postprocessing.assembly",Environment::logVerbosityLevel()>0);

    tic();
    tic();
    b.solve( _solution=PP, _rhs=ell, _name="sc.post", _local=true);
    toc("postprocessing.solve.local",Environment::logVerbosityLevel()>0);
    ppp=PP(0_c);
    tic();
    quickstart::matchCellMean( ppp, pp, P0dh );
    toc("postprocessing.solve.correction",Environment::logVerbosityLevel()>0);
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
    e->add( "potentialpp", ppp );
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
        return hdg_laplacian<FEELPP_DIM>();
    }
    catch( ... )
    {
        handleExceptions();
    }
    return 1;
}
