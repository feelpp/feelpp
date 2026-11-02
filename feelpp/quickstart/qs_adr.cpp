/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2014-2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include "discretization.hpp"
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf_eval.hpp>

int main(int argc, char**argv )
{
    // tag::env[]
    using namespace Feel;
    using Feel::cout;
    try
    {
        po::options_description adroptions( "Advection Diffusion Reaction options" );
        adroptions.add_options()
            ( "order", po::value<int>()->default_value( 2 ), "approximation order" )
            ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" )
            ( "stabilisation", po::value<bool>()->default_value( false ), "enable/disable stabilisation disable=0, enable=1" )
            ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
            ;

        Environment env( _argc=argc, _argv=argv,
                        _desc=adroptions,
                        _about=about(_name="qs_adr",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));
        // end::env[]

        // tag::mesh_space[]
        tic();
        auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,Dynamic>>( quickstart::geometryOrder() ));
        quickstart::checkMeshOrder( mesh );
        toc("loadMesh");

        tic();
        auto Vh = Pch<Dynamic>( mesh, quickstart::checkedOrder( "order" ) );
        auto u = Vh->element("u");
        auto gamma = expr(soption(_name="functions.gamma"));
        auto beta = expr<FEELPP_DIM,1>(soption(_name="functions.beta"));
        auto sigma = expr<FEELPP_DIM,FEELPP_DIM>(soption(_name="functions.sigma"));
        auto f = expr( soption(_name="functions.f"), "f" );
        auto r_1 = expr( soption(_name="functions.a"), "a" ); // Robin left hand side expression
        auto r_2 = expr( soption(_name="functions.b"), "b" ); // Robin right hand side expression
        auto n = expr( soption(_name="functions.c"), "c" ); // Neumann expression
        auto g = expr( soption(_name="functions.g"), "g" );
        auto v = Vh->element( g, "g" );
        toc("Vh");
        // end::mesh_space[]

        // tag::forms[]
        tic();
        auto l = form1( _test=Vh );
        l = integrate(_range=elements(mesh),
                    _expr=f*id(v));
        l+=integrate(_range=markedfaces(mesh,"Robin"), _expr=r_2*id(v));
        l+=integrate(_range=markedfaces(mesh,"Neumann"), _expr=n*id(v));
        toc("l");

        tic();
        auto a = form2( _trial=Vh, _test=Vh);
        a = integrate(_range=elements(mesh),
                    _expr=grad(v)*(sigma*trans(gradt(u))) + (gradt(u)*beta)*id(v) + gamma*idt(u)*id(v) );

        // tag::stab[]
        if ( boption("stabilisation") )
        {
            auto stable = expr( soption(_name="delta",_prefix="functions") );
            auto epsilon=(trace(sigma)/FEELPP_DIM);
            auto delta = stable*constant(1.0)/(1.0/h() + epsilon/(h()*h()));
            auto  L = -epsilon*laplacian(v)+ grad(v)*beta + gamma*id(v);
            auto  Lt = -epsilon*laplaciant(u)+ gradt(u)*beta + gamma*idt(u);

            a+= integrate(_range=elements(mesh), _expr=delta*L*Lt );
            l+= integrate(_range=elements(mesh), _expr=delta*f*L );
        }
        // end::stab[]
        
        a+=integrate(_range=markedfaces(mesh,"Robin"), _expr=r_1*idt(u)*id(v));
        a+=on(_range=markedfaces(mesh,"Dirichlet"), _rhs=l, _element=u, _expr=g );
        //! if no markers Robin Neumann or Dirichlet are present in the mesh then
        //! impose Dirichlet boundary conditions over the entire boundary
        if ( !mesh->hasAnyMarker({"Robin", "Neumann","Dirichlet"}) )
            a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g );
        toc("a");

        tic();
        //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
        if ( !boption( "no-solve" ) )
            a.solve(_rhs=l,_solution=u);
        toc("a.solve");

        // end::forms[]
        cout << "||u-u_h||_L2=" << normL2(_range=elements(mesh), _expr=idv(u)-g) << std::endl;

        // tag::export[]
        tic();
        auto e = exporter( _mesh=mesh );
        e->addRegions();
        e->add( "uh", u );
        e->add( "u", v );
        e->save();
        toc("Exporter");
        return 0;
        // end::export[]
    }
    catch(...)
    {
        handleExceptions();
    }
    return 1;
}
// end::global[]
