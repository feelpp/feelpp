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
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelmesh/complement.hpp>
namespace Feel {


inline
po::options_description
makeOptions()
{
    po::options_description hdgoptions( "HDG options" );
    hdgoptions.add_options()
        ( "k", po::value<std::string>()->default_value( "-1" ), "diffusion coefficient" )
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
        ;
    return hdgoptions;
}
 
inline
AboutData
makeAbout()
{
    AboutData about( "qs_hdg_laplacian_ibc" ,
                     "qs_hdg_laplacian_ibc" ,
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

    tic();
    auto mesh = loadMesh( new Mesh<Simplex<Dim>> );
    int nbIbc = nelements(markedfaces(mesh,"Ibc"),true) >= 1 ? 1 : 0;
    int nbIbcOde = nelements(markedfaces(mesh,"IbcOde"),true) >= 1 ? 2 : 0;
    std::map<std::string,std::pair<int,std::string>> ibcs;
    if ( nbIbc )
    {
        ibcs["Ibc"]= std::pair{0,"Ibc"};
    }
    if ( nbIbcOde )
    {
        ibcs["Ibc"]= std::pair{0,"IbcOde"};
        ibcs["Ode"]= std::pair{1,"IbcOde"};
    }
    toc("mesh",true);

    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;
    
    auto K = expr(soption("k"));
    auto lambda = cst(1.)/K;
    
#if defined(FEELPP_HAS_SYMPY)

    std::map<std::string,std::string> locals{{"dim",std::to_string(Dim)},{"k",soption("k")},{"p",soption("solution.sympy.p")},{"grad_p",""}, {"u",""}, {"un",""}, {"f",""}};
    Feel::pyexprFromFile( Environment::expand(soption("pyexpr.filename")), locals );

    for( auto d: locals )
        cout << d.first << ":" << d.second << std::endl;

    std::string p_exact_str = locals.at("p");
    std::string u_exact_str = locals.at("u");
    auto p_exact = expr( p_exact_str );
    auto u_exact = expr<Dim,1>( u_exact_str );
    auto un_exact = expr( locals.at("un") );
    auto f_exact = expr( locals.at("f") );
    std::map<std::string,double> ibc_exact_map;
    for( auto const& [ibc_type, ibc_data ] : ibcs )
        if ( ibc_data.second == "IbcOde" && ibc_type == "Ibc" )
            ibc_exact_map[ibc_type] = 0;
        else
            ibc_exact_map[ibc_type] = integrate(markedfaces(mesh,ibc_data.second), un_exact).evaluate()(0,0);

#else
    std::string p_exact_str = soption("solution.p");
    std::string u_exact_str = soption("solution.u");
    auto p_exact = expr(p_exact_str);
    auto u_exact = expr<Dim,1>(u_exact_str);
    auto un_exact = trans(u_exact)*N();
    auto f_exact = expr( soption( "functions.f") );
    auto ibc_exact = expr( soption("functions.i") );
#endif

    // ****** Hybrid-mixed formulation ******
    // We treat Vh, Wh, and Mh separately
    tic();

    auto Vh = Pdhv<OrderP>( mesh, true );
    auto Wh = Pdh<OrderP>( mesh, true );
    auto select_faces = [mesh]( auto const& ewrap ) {
        auto const& e = unwrap_ref( ewrap );
        if ( e.hasMarker() && ( e.marker().value() == mesh->markerName( "Ibc" )  ||
                                e.marker().value() == mesh->markerName( "IbcOde" ) ) )
            return true;
        return false; };
    auto complement_integral_bdy = complement(faces(mesh),select_faces);
    auto complement_integral_bdy_boundary = complement(boundaryfaces(mesh),select_faces);

    auto face_mesh = createSubmesh( _mesh=mesh, _range=complement_integral_bdy, _update=0 );
    // auto face_mesh = createSubmesh( _mesh=mesh, _range=faces(mesh ), _update=0 );
    auto Mh = Pdh<OrderP>( face_mesh,true );
    auto ibc_mesh = createSubmesh( _mesh=mesh, _range=markedfaces(mesh, {"Ibc","IbcOde"}), _update=0 );
    auto Ch = Pch<0>( ibc_mesh, true );

    toc("spaces",true);
    auto P0dh = Pdh<0>(mesh);
    auto Xh = Pdh<0>(face_mesh);
    auto uf = Xh->element(cst(1.));
    auto cgXh = Pch<OrderP>(mesh);

    cout << "Exact potential if applicable: " << p_exact_str << "\n"
         << "Exact flux if applicable: " << u_exact_str << "\n";

    cout << "#elts: " << mesh->numGlobalElements() << std::endl
         << "#faces: " << mesh->numGlobalFaces() << std::endl
         << "#facesMh: " << face_mesh->numGlobalElements() << std::endl
         << "Vh<" << OrderP << "> : " << Vh->nDof() << std::endl
         << "Wh<" << OrderP << "> : " << Wh->nDof() << std::endl
         << "Mh<" << OrderP << "> : " << Mh->nDof() << std::endl;
    if( nbIbc > 0 || nbIbcOde)
        cout << "Ch<0> : " << Ch->nDof() << std::endl;
    cout << mesh->numGlobalElements()  << " " << mesh->numGlobalFaces() << " "
         << Vh->nDof() << " " << Wh->nDof() << " " << Mh->nDof() << " "
         << cgXh->nDof() << std::endl;

    if ( boption( "solvecg" ) == true )
    {
        auto u = cgXh->element();
        tic();
        tic();
        auto l = form1( _test=cgXh );
        l = integrate( _range=elements(mesh), _expr=-f_exact*id(u) );
        l += integrate(_range=markedfaces(mesh,"Neumann"),
                       _expr=id(u)*un_exact );
        toc("cg.assembly.l",FLAGS_v>0);
        tic();
        auto a = form2(_trial=cgXh, _test=cgXh );
        a = integrate( _range=elements(mesh), _expr=gradt(u)*trans(grad(u)) );
        a += on(_range=markedfaces(mesh,"Dirichlet"),
                _rhs=l,
                _element=u,
                _expr=p_exact);
        toc("cg.assembly.a",FLAGS_v>0);
        toc("cg.assembly",FLAGS_v>0);
        tic();
        a.solve( _rhs=l,_solution=u );
        toc("cg.solve",FLAGS_v>0);
        auto norms = [=]( std::string const& solution ) ->std::map<std::string,double>
            {
                tic();
                double l2 = normL2(_range=elements(mesh), _expr=idv(u)-expr(solution) );
                toc("L2 error norm");
                tic();
                double h1 = normH1(_range=elements(mesh), _expr=idv(u)-expr(solution), _grad_expr=gradv(u)-grad<2>(expr(solution)) );
                toc("H1 error norm");
                tic();
                double semih1 = normL2(_range=elements(mesh), _expr=gradv(u)-grad<2>(expr(solution)) );
                toc("semi H1 error norm");
                return { { "L2", l2 }, {  "H1", h1 }, {"semih1",semih1} };
            };

        int status = checker(p_exact_str).runOnce( norms, rate::hp( mesh->hMax(), cgXh->fe()->order() ) );

    }
    auto u = Vh->element( "u" );
    auto v = Vh->element( "v" );
    auto p = Wh->element( "p" );
    auto q = Wh->element( "q" );
    auto w = Wh->element( "w" );
    auto phat = Mh->element( "phat" );
    auto l = Mh->element( "lambda" );
    auto mu = Ch->element( "mu" );
    auto nu = Ch->element( "nu" );

    tic();
    auto ibcSpaces = dynProductPtr( nbIbc+nbIbcOde, Ch);
    //ibcSpaces->setProperties( {"Ibc","Ode" } );
    std::vector<std::string> props;
    if ( nbIbc )
        props.push_back( "Ibc" );
    if ( nbIbcOde )
    {
        props.push_back( "Ibc" );
        props.push_back( "Ode" );
    }
    ibcSpaces->setProperties( props );
    auto ps = product2( ibcSpaces, Vh, Wh, Mh );
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
    for( auto const& [ ibc_type, ibc_data ]  : ibcs )
    {
        auto const& [ibc_space_index,ibc_marker] = ibc_data;
        Feel::cout << "rhs assembly " << ibc_type << " -- block:" << ibc_space_index << " -- marker: " << ibc_marker << " -- value:" << ibc_exact_map[ibc_type] << std::endl;
        double meas = integrate( _range=markedfaces(mesh, ibc_marker ),
                                 _expr=cst(1.)).evaluate()(0,0);
        rhs(3_c,ibc_space_index) += integrate( _range=markedfaces(mesh, ibc_marker ),
                                               _expr=ibc_exact_map[ibc_type]*id(nu)/meas);
    }

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
                            _expr=( idt(phat)*(leftface(normal(v))+
                                               rightface(normal(v)))) );

    a(0_c,2_c) += integrate(_range=complement_integral_bdy_boundary,
                            _expr=idt(phat)*(normal(v)));
    toc("a(0,2)",FLAGS_v>0);

    tic();
    for( auto const& [ibc_type,ibc_data]  : ibcs )
    {
        auto const& [ibc_space_index,ibc_marker] = ibc_data;
        if ( ibc_type == "Ibc" )
            a(0_c,3_c,0,ibc_space_index) += integrate(_range=markedfaces(mesh,ibc_marker),
                                                      _expr=idt(mu)*(trans(id(u))*N()) );
    }
    toc("a(0,3)",FLAGS_v>0);

    //
    // Second row a(1_c,:)
    //
    tic();
    a(1_c,0_c) += integrate(_range=elements(mesh),_expr=-(id(w)*divt(u)));
    toc("a(1,0)",FLAGS_v>0);

    tic();
    a(1_c,1_c) += integrate(_range=internalfaces(mesh),
                            _expr=-tau_constant *
                            ( leftfacet( idt(p))*leftface(id(w)) +
                              rightfacet( idt(p))*rightface(id(w) )));
    a(1_c,1_c) += integrate(_range=boundaryfaces(mesh),
                            _expr=-(tau_constant * id(w)*idt(p)));
    toc("a(1,1)",FLAGS_v>0);

    tic();
    a(1_c,2_c) += integrate(_range=internalfaces(mesh),
                            _expr=tau_constant * idt(phat) *
                            ( leftface( id(w) )+
                              rightface( id(w) )));

    a(1_c,2_c) += integrate(_range=complement_integral_bdy_boundary,
                            _expr=tau_constant * idt(phat) * id(w) );
    toc("a(1,2)",FLAGS_v>0);

    tic();
    for( auto const& [ibc_type,ibc_data]  : ibcs )
    {
        auto const& [ibc_space_index,ibc_marker] = ibc_data;
        if ( ibc_type == "Ibc" )
            a(1_c,3_c,0,ibc_space_index) += integrate( _range=markedfaces(mesh,ibc_marker),
                                                       _expr=tau_constant*idt(mu)*id(w) );
        
    }
    toc("a(1,3)", FLAGS_v>0);

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
                            _expr=-(1-0.5*boption("sc.condense"))*tau_constant * idt(phat) * id(l) );
    a(2_c,2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=-tau_constant * idt(phat) * id(l)  );
    a(2_c,2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                            _expr=idt(phat) * id(l) );
    toc("a(2,2)",FLAGS_v>0);

    //
    // Fourth row a(3_c,:)
    //
    
    for( auto const& [ibc_type,ibc_data]  : ibcs )
    {
        auto const& [ibc_space_index,ibc_marker] = ibc_data;
        Feel::cout << "row 3_c assembly " << ibc_type << " -- block:" << ibc_space_index << " -- marker: " << ibc_marker << std::endl;
        tic();
        if ( ibc_type == "Ibc" )
            a(3_c,0_c,ibc_space_index,0) += integrate( _range=markedfaces(mesh,ibc_marker),
                                                       _expr=(trans(idt(u))*N())*id(nu) );
        toc("a(3,0)",FLAGS_v>0);

        tic();
        if ( ibc_type == "Ibc" )
            a(3_c,1_c,ibc_space_index,0) += integrate( _range=markedfaces(mesh,ibc_marker),
                                                       _expr=tau_constant*idt(p)*id(nu) );
        toc("a(3,1)",FLAGS_v>0);

        if ( ibc_type == "Ode" )
            a(3_c,3_c,ibc_space_index-1,ibc_space_index) += integrate( _range=markedfaces(mesh,ibc_marker),
                                                                       _expr=-id(nu)*idt(mu) );
        
        tic();
        double c = 1.;
        if ( ibc_type == "Ibc" )
            c = -tau_constant.expression().value();
        a(3_c,3_c,ibc_space_index,ibc_space_index) += integrate( _range=markedfaces(mesh,ibc_marker),
                                                                 _expr=c*id(nu)*idt(mu) );
        toc("a(3,3)",FLAGS_v>0);
    }
    toc("matrices",true);
        
    
    
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
    toc("postprocessing.space",FLAGS_v>0);
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
    toc("postprocessing.solve",FLAGS_v>0);
    toc("postprocessing",true);

    

    tic();
    v.on( _range=elements(mesh), _expr=u_exact );
    q.on( _range=elements(mesh), _expr=p_exact );
    auto e = exporter( _mesh=mesh );
    e->setMesh( mesh );
    e->addRegions();
    e->add( "flux", U(0_c) );
    e->add( "potential", U(1_c) );
    e->add( "potentialpp", PP(0_c) );
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
    auto norms_ppp = [=]( std::string const& solution ) ->std::map<std::string,double>
        {
            tic();
            double l2err_p = 1e+30;
            l2err_p = normL2( elements(mesh), expr(solution) - idv(ppp) );
            toc("L2 error norm ppp");
            return { { "L2", l2err_p }};
        };
#if 0
    int status = checker().runOnce( {norms_p,norms_u},
                                    { rate::hp( mesh->hMax(), Vh->fe()->order() ), rate::hp( mesh->hMax(), Vh->fe()->order() ) } );
#else
    int status1 = checker("L2 Convergence potential",p_exact_str).runOnce( norms_p, rate::hp( mesh->hMax(), Vh->fe()->order() ) );
    int status2 = checker("L2 Convergence flux",u_exact_str).runOnce( norms_u, rate::hp( mesh->hMax(), Vh->fe()->order() ) );
    int status3 = checker("L2 Convergence post-processed potential",p_exact_str).runOnce( norms_ppp, rate::hp( mesh->hMax(), Vh->fe()->order()+1 ) );
#endif

    
    // end::check[]

    return status1 || status2 || status3;
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
        return !hdg_laplacian<FEELPP_DIM,1>();
    if ( ioption( "order" ) == 2 )
        return !hdg_laplacian<FEELPP_DIM,2>();
#if 0

    if ( ioption( "order" ) == 3 )
        return !hdg_laplacian<FEELPP_DIM,3>();
    if ( ioption( "order" ) == 4 )
        return !hdg_laplacian<FEELPP_DIM,4>();
#endif
    return 1;
}
