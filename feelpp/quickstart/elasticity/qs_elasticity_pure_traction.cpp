/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>
#include "nullspace-rigidbody.hpp"

int main(int argc, char**argv )
{
    using namespace Feel;
    try
    {
        po::options_description laplacianoptions( "Elasticity options" );
        laplacianoptions.add_options()
            ( "domain", po::value<std::string>()->default_value( "" ), "name of the subdomain where rigid body are applied, if empty then it is the full domain" )
            ( "E", po::value<double>()->default_value( 1.0e6 ), "Young modulus" )
            ( "nu", po::value<double>()->default_value( 0.3 ), "Poisson ratio" )
            ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
            ( "weakdir", po::value<bool>()->default_value( false ), "use weak dirichlet" )
            ( "gamma", po::value<double>()->default_value( 100 ), "penalisation term" )
            ( "nullspace", po::value<bool>()->default_value( false ), "add null space" )
            ;

        Environment env( _argc=argc, _argv=argv,
                    _desc=laplacianoptions,
                    _about=about(_name="qs_elasticity",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

        tic();
        auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
        toc("loadMesh");

        tic();
        auto Vh = Pchv<1>( mesh );
        auto V0hv = Pchv<0>( mesh );
    #if FEELPP_DIM == 2
        auto V0h = Pch<0>( mesh );
    #else
        auto V0h = Pchv<0>( mesh );
    #endif
        auto Xh = product( Vh, V0hv, V0h );
        toc("Vh");
        

        auto u = Vh->element("u");
        auto v = Vh->element("v");
        auto mv = V0hv->element("mv");
        auto nv = V0hv->element("nv");
        auto m = V0h->element("m");
        auto n = V0h->element("n");
        auto nu = doption(_name="nu");
        auto E = doption(_name="E");
        auto lambda = E*nu/( (1+nu)*(1-2*nu) );
        auto mu = E/(2*(1+nu));
        auto deft = sym(gradt(u));
        auto def = sym(grad(u));
        auto Id = eye<FEELPP_DIM,FEELPP_DIM>();
        auto sigmat = lambda*trace(deft)*Id + 2*mu*deft;
        auto sigma = lambda*trace(def)*Id + 2*mu*def;
        auto f = expr<FEELPP_DIM,1>( soption(_name="functions.f"), "f" );
        auto g = expr<FEELPP_DIM,1>( soption(_name="functions.g"), "g" );
        auto tau = expr<FEELPP_DIM,1>( soption(_name="functions.tau"), "tau" );
    #if FEELPP_DIM == 2
        auto omega = expr( soption(_name="functions.omega"), "omega" );
    #else
        auto omega = expr<3,1>( soption(_name="functions.omega"), "omega" );
    #endif
        auto rigid_body_domain_name = soption(_name="domain");
        
        tic();
        auto l = blockform1( Xh, solve::strategy::monolithic, backend() );
        l(0_c) = integrate(_range=elements(mesh),
                        _expr=inner(f,id(v)));
        l(0_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                            _expr=inner(g,id(v)));

        auto rigid_body_domain = (rigid_body_domain_name.empty())?elements(mesh):markedelements(mesh, rigid_body_domain_name );
        std::cout << "surface rigid domain=" << integrate( _range=rigid_body_domain, _expr=cst(1.0) ).evaluate() << std::endl; 
        // translation
        l(1_c) += integrate(_range=elements(mesh),
                            _expr=inner(tau,id(mv)));
        // rotation
        l(2_c) += integrate(_range=rigid_body_domain,
                            _expr=inner(omega,id(m)));
        toc("l");

        tic();
        auto a = blockform2( Xh, solve::strategy::monolithic, backend() );
        a(0_c, 0_c) = integrate(_range=elements(mesh),
                                _expr=inner( sigmat, grad(v) ) );
        a(0_c,1_c) += integrate(_range=elements(mesh),
                                _expr=trans(idt(mv))*id(v));
        a(1_c,0_c) += integrate(_range=elements(mesh),
                                _expr=trans(id(mv))*idt(v));
    #if FEELPP_DIM == 2
        a(0_c,2_c) += integrate(_range=rigid_body_domain,
                                _expr=trans(idt(m))*curlx(v));
        a(2_c,0_c) += integrate(_range=rigid_body_domain,
                                _expr=trans(id(m))*curlxt(v));
    #else
        a(0_c,2_c) += integrate(_range=rigid_body_domain,
                                _expr=trans(idt(m))*curl(v));
        a(2_c,0_c) += integrate(_range=rigid_body_domain,
                                _expr=trans(id(m))*curlt(v));
    #endif
        toc("a");

        //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
        if ( !boption( "no-solve" ) )
        {
            tic();
            auto b = backend();
            auto rigidBodyModes = std::make_shared<NullSpace<double>>( b, qsNullSpace( Vh ) );
            b->attachNearNullSpace( rigidBodyModes );
            if ( boption(_name="nullspace") )
                b->attachNullSpace( rigidBodyModes );

            auto U=Xh.element();
            a.solve(_rhs=l,_solution=U);
            u=U(0_c);
            //m=U(1_c);
            toc("a.solve");
            
        }

        tic();
        auto e = exporter( _mesh=mesh );
        e->addRegions();
        e->add( "u", u );
        e->save();
        toc("Exporter");
        return 0;
    }
    catch( ... )
    {
        handleExceptions();
    }
    return 1;
}
