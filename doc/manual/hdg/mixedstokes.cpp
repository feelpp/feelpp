#ifndef _MIXEDSTOKES_HPP
#define _MIXEDSTOKES_HPP

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/complement.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <feel/feelts/bdf.hpp>
#include <boost/algorithm/string.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>

using namespace Feel;
using namespace Feel::vf;

int main(int argc, char** argv)
{
    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=about(_name="mixed-stokes",
                                        _author="Romain Hild",
                                        _email=""),
                           _desc=feel_options()
                           );

    using convex_type = Simplex<FEELPP_DIM>;
    using mesh_type = Mesh<convex_type>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using face_convex_type = Simplex<FEELPP_DIM-1>;
    using face_mesh_type = Mesh<face_convex_type>;
    using face_mesh_ptrtype = std::shared_ptr<face_mesh_type>;

    auto mesh = loadMesh( new Mesh<Simplex<FEELPP_DIM> >);
    auto face_mesh = createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0);
    auto Wh = Pdhm<FEELPP_ORDER>( mesh );
    auto Vh = Pdhv<FEELPP_ORDER>( mesh );
    auto Ph = Pdh<FEELPP_ORDER>( mesh );
    auto Mh = Pdhv<FEELPP_ORDER>( face_mesh );
    auto Ch = Pdh<0>( face_mesh );

    // auto Ps = product(Wh,Vh,Ph,Mh,Ch);
    auto Ps = product(Wh,Vh,Ph,Mh);
    auto bbf = blockform2(Ps);
    auto blf = blockform1(Ps);

    auto L = Wh->element();
    auto G = Wh->element();
    auto u = Vh->element();
    auto v = Vh->element();
    auto p = Ph->element();
    auto q = Ph->element();
    auto uhat = Mh->element();
    auto mu = Mh->element();
    auto rho = Ch->element();
    auto phi = Ch->element();

    auto I = Id<FEELPP_DIM>();
    auto nu = doption("parameters.nu");
    auto tau = doption("parameters.tau");
    auto f = expr<FEELPP_DIM,1>(soption("functions.f"));
    auto u_expr = expr<FEELPP_DIM,1>(soption("functions.u"));
    auto p_expr = expr(soption("functions.p"));


    // (L,G)
    bbf(0_c,0_c) = integrate(elements(mesh), inner(idt(L),id(G)));
    // (u,div(G))
    bbf(1_c,0_c) += integrate(elements(mesh), inner(idt(u),div(G)));
    // <uhat,Gn>
    bbf(3_c,0_c) += integrate(internalfaces(mesh),
                              trans(idt(uhat))*leftface(id(G)*N())
                              + trans(idt(uhat))*rightface(id(G)*N()) );
    bbf(3_c,0_c) += integrate(boundaryfaces(mesh),
                              -trans(idt(uhat))*(id(G)*N()) );

    // (nu L, grad(v))
    bbf(0_c,1_c) += integrate(elements(mesh), nu*inner(idt(L),grad(v)));
    // (pI,grad(v))
    bbf(2_c,1_c) += integrate(elements(mesh), inner(idt(p)*I,grad(v)));

    // <-nu Ln,v>
    bbf(0_c,1_c) += integrate(internalfaces(mesh),
                              -nu*leftfacet(trans(idt(L)*N()))*leftface(id(v))
                              -nu*rightfacet(trans(idt(L)*N()))*rightface(id(v)) );
    bbf(0_c,1_c) += integrate(boundaryfaces(mesh),
                              -nu*trans(idt(L)*N())*id(v) );

    // <tau (u ox n)n,v>
    bbf(1_c,1_c) += integrate(internalfaces(mesh),
                              tau*trans(leftfacet(idt(u)*trans(N())*N()))*leftface(id(v))
                              + tau*trans(rightfacet(idt(u)*trans(N())*N()))*rightface(id(v)) );
    bbf(1_c,1_c) += integrate(boundaryfaces(mesh),
                              tau*trans(idt(u)*trans(N())*N())*id(v) );

    // <pIn,v>
    bbf(2_c,1_c) += integrate(internalfaces(mesh),
                             leftfacet(trans(idt(p)*I*N()))*leftface(id(v))
                             + rightfacet(trans(idt(p)*I*N()))*rightface(id(v)) );
    bbf(2_c,1_c) += integrate(boundaryfaces(mesh),
                              trans(idt(p)*I*N())*id(v) );

    // <-tau (uhat ox n)n,v>
    bbf(3_c,1_c) += integrate(internalfaces(mesh),
                              - tau*trans(leftfacet(idt(uhat)*trans(N())*N()))*leftface(id(v))
                              - tau*trans(rightfacet(idt(uhat)*trans(N())*N()))*rightface(id(v)) );
    bbf(3_c,1_c) += integrate(boundaryfaces(mesh),
                              - tau*trans(idt(uhat)*trans(N())*N())*id(v) );

    // (-u,grad(q))
    bbf(1_c,2_c) += integrate(elements(mesh), -trans(idt(u))*trans(grad(q)) );

    // <uhat.n,q> lack <uhat.n,-qbar>
    bbf(3_c,2_c) += integrate(internalfaces(mesh),
                              trans(idt(uhat))*leftface(N()*id(q))
                              + trans(idt(uhat))*rightface(N()*id(q)) );
    bbf(3_c,2_c) += integrate(boundaryfaces(mesh),
                              trans(idt(uhat))*N()*id(q) );

    // <-nu Ln, mu>
    bbf(0_c,3_c) += integrate(internalfaces(mesh),
                              -nu*leftfacet(trans(idt(L)*N()))*leftface(id(mu))
                              -nu*rightfacet(trans(idt(L)*N()))*rightface(id(mu)) );
    bbf(0_c,3_c) += integrate(boundaryfaces(mesh),
                              -nu*trans(idt(L)*N())*id(mu) );
    // <tau (u ox n)n, mu>
    bbf(1_c,3_c) += integrate(internalfaces(mesh),
                              tau*trans(leftfacet(idt(u)*trans(N())*N()))*leftface(id(mu))
                              + tau*trans(rightfacet(idt(u)*trans(N())*N()))*rightface(id(mu)) );
    bbf(1_c,3_c) += integrate(boundaryfaces(mesh),
                              tau*trans(idt(u)*trans(N())*N())*id(mu) );

    // <pIn, mu>
    bbf(2_c,3_c) += integrate(internalfaces(mesh),
                             leftfacet(trans(idt(p)*I*N()))*leftface(id(mu))
                             + rightfacet(trans(idt(p)*I*N()))*rightface(id(mu)) );
    bbf(2_c,3_c) += integrate(boundaryfaces(mesh),
                              trans(idt(p)*I*N())*id(mu) );

    // <-tau (uhat ox n)n,mu>
    bbf(3_c,3_c) += integrate(internalfaces(mesh),
                              - tau*trans(leftfacet(idt(uhat)*trans(N())*N()))*leftface(id(mu))
                              - tau*trans(rightfacet(idt(uhat)*trans(N())*N()))*rightface(id(mu)) );
    bbf(3_c,3_c) += integrate(boundaryfaces(mesh),
                              - tau*trans(idt(uhat)*trans(N())*N())*id(mu) );

    // // <uhat.n,phi>
    // bbf(3_c,4_c) += integrate(internalfaces(mesh),
    //                           leftfacet(trans(idt(uhat))*N())*id(phi)
    //                           + rightfacet(trans(idt(uhat))*N())*id(phi) );
    // bbf(3_c,4_c) += integrate(boundaryfaces(mesh),
    //                           trans(idt(uhat))*N()*id(phi) );

    // (p,1) ?? 0_c,1_c,2_c,3_c,4_c ??
    // bbf(2_c,0_c) += integrate(elements(mesh), idt(p) );
    // bbf(2_c,1_c) += integrate(elements(mesh), idt(p) );
    bbf(2_c,2_c) += integrate(elements(mesh), idt(p) );
    // bbf(2_c,3_c) += integrate(faces(mesh), leftfacet(idt(p)) + rightfacet(idt(p)) );
    // bbf(2_c,4_c) += integrate(faces(mesh), leftfacet(idt(p)) + rightfacet(idt(p)) );


    // blf(1_c) = integrate(elements(mesh), inner(f,id(v)) );
    blf(1_c) = integrate(elements(mesh), trans(grad(u_expr)*u_expr)*id(v) );
    blf(3_c) += integrate(boundaryfaces(mesh), trans((-nu*grad(u_expr)+p_expr*I)*N())*id(mu) );
    // bbf(1_c,0_c).on(boundaryfaces(mesh), u_expr);
    // bbf(1_c,1_c).on(boundaryfaces(mesh), u_expr);
    // bbf(1_c,2_c).on(boundaryfaces(mesh), u_expr);
    // bbf(1_c,3_c).on(boundaryfaces(mesh), u_expr);
    // bbf(1_c,4_c).on(boundaryfaces(mesh), u_expr);

    auto U = Ps.element();
    bbf.solve(_rhs=blf, _solution=U);

    auto u_ex = Vh->element();
    u_ex.on(elements(mesh), u_expr);
    auto p_ex = Ph->element();
    p_ex.on(elements(mesh), p_expr);

    auto e = exporter(mesh);
    // e->add("L", U(0_c));
    e->add("u", U(1_c));
    e->add("p", U(2_c));
    e->add("u_ex", u_ex);
    e->add("p_ex", p_ex);
    e->save();

    return 0;
}

#endif
