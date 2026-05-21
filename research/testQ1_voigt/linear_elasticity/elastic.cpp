/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>
// #include "nullspace-rigidbody.hpp"


int main(int argc, char**argv )
{
    using namespace Feel;
    try
    {
        po::options_description laplacianoptions( "Elasticity options" );
        laplacianoptions.add_options()
            ( "E", po::value<double>()->default_value( 1.0e6 ), "Young modulus" )
            ( "nu", po::value<double>()->default_value( 0.3 ), "Poisson ratio" )
            ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
            ( "weakdir", po::value<bool>()->default_value( false ), "use weak dirichlet" )
            ( "gamma", po::value<double>()->default_value( 100 ), "penalisation term" )
            ( "moment_x", po::value<bool>()->default_value( false ), "Moment x test" )

            // ( "nullspace", po::value<bool>()->default_value( false ), "add null space" )
            ;

        Environment env( _argc=argc, _argv=argv,
                    _desc=laplacianoptions,
                    _about=about(_name="qs_elasticity",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));


        tic();
        // auto mesh = loadMesh(_mesh=new Mesh<Simplex<3,1>>);
        auto mesh_file = Environment::expand( soption(_name = "gmsh.filename") );
        auto mesh = loadMesh(_mesh = new Mesh<Hypercube<3,1>>(), _filename = mesh_file, _scale = 1, _straighten = false);

        toc("loadMesh");

        tic();
        auto Vh = Pchv<1>( mesh );
        toc("Vh");

        auto u = Vh->element("u");
        auto v = Vh->element("v");
        auto nu = doption(_name="nu");
        auto E = doption(_name="E");
        auto lambda = E*nu/( (1+nu)*(1-2*nu) );
        auto mu = E/(2*(1+nu));
        auto deft = sym(gradt(u));
        auto def = sym(grad(u));
        auto Id = eye<3,3>();
        auto sigmat = lambda*trace(deft)*Id + 2*mu*deft;
        auto sigma = lambda*trace(def)*Id + 2*mu*def;
        auto f = expr<3,1>( soption(_name="functions.f"), "f" );
        auto g = expr<3,1>( soption(_name="functions.g"), "g" );

        tic();
        auto l = form1( _test=Vh );
        // possible d'appliquer juste une force sur un point ?

        if( boption(_name = "moment_x" ) ) {
            auto force = vec( cst(0.), -6*Pz(), 6*(Py() - 0.5) );
            l = integrate(_range= markedfaces(mesh, "ForceApply"), _expr = inner( force, id(v) ));
        }
        else
            l = integrate(_range= markedfaces(mesh, "ForceApply"), _expr = inner( f, id(v) ));


        toc("l");

        tic();
        auto a = form2( _trial=Vh, _test=Vh);
        a = integrate(_range = elements(mesh), _expr = inner( sigmat, grad(v) ));

        // Appliquer sur des points ? avec markedpoints 
        a+=on(_range=markedfaces(mesh,"Dirichlet"), _rhs=l, _element=u, _expr=g );  
        toc("a");


        tic();
        a.solve(_rhs = l, _solution = u, _rebuild=true);
        toc("a.solve");


        tic();
        auto e = exporter( _mesh = mesh );
        e->addRegions();
        e->add( "u", u );
        e->save();
        toc("Exporter");

        l.vector().printMatlab("form1.m");
        a.matrix().printMatlab("form2.m");
        u.printMatlab("solution.m");   // print de la solution noeud par noeud

        return 0;
    }
    catch(...)
    {
        handleExceptions();
    }
    return 1;
}