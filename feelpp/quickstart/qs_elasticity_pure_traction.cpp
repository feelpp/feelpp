/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>
#include "nullspace-rigidbody.hpp"

namespace Feel
{
enum class ConstraintImplType
{
    LagrangeMultiplier,
    Nodal
};
template <typename SpaceType, typename ExprTransT, typename ExprRotT>
auto rigidBodyMotion( std::shared_ptr<SpaceType> const& space,
                        ExprTransT const& tau,
                        ExprRotT const& omega,
                        std::string const& rigid_body_domain_name="" )
{
    auto mesh = space->mesh();
    auto v = space->element();
    auto w = space->element();
    auto rigid_body_domain = (rigid_body_domain_name.empty())?elements(mesh):markedelements(mesh, rigid_body_domain_name );
    std::cout << "surface rigid domain=" << integrate( _range=rigid_body_domain, _expr=cst(1.0) ).evaluate() << std::endl;
    double meas = integrate( _range=elements(mesh), _expr=cst(1.0) ).evaluate()(0,0);
    auto center_of_mass = integrate( _range=elements(mesh), _expr=P()/meas ).evaluate();
    std::cout << fmt::format( "meas={}, center of mass={}", meas, center_of_mass ) << std::endl;
    auto R = rigidRotationMatrix( omega.evaluate() );
    std::cout << "R=" << R << std::endl;
    v.on(_range=elements(mesh), _expr=constant(R)*(P()-constant(center_of_mass))+constant(center_of_mass)-P());
    w.on(_range=elements(mesh), _expr=tau );
    // return translation, rotation
    return std::tuple{w,v};
}
template<ConstraintImplType c_impl, typename SpaceType>
auto elasticity( std::shared_ptr<SpaceType>& Vh )
{
    auto mesh = Vh->mesh();
    toc("Vh");


    auto u = Vh->element("u");
    auto v = Vh->element("v");
    auto w = Vh->element("w");

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
    auto rigid_body_domain = (rigid_body_domain_name.empty())?elements(mesh):markedelements(mesh, rigid_body_domain_name );
    auto getSpace = []( auto const& space )
    {
        auto mesh = space->mesh();
        if constexpr ( c_impl == ConstraintImplType::LagrangeMultiplier )
        {
            auto V0hv = Pchv<0>( mesh );
            #if FEELPP_DIM == 2
            auto V0h = Pch<0>( mesh );
            #else
            auto V0h = Pchv<0>( mesh );
            #endif
            return product( space, V0hv, V0h );
        }  
        if constexpr ( c_impl == ConstraintImplType::Nodal )
        {
            return product( space );
        }
        else
        {
            return product( space );
        }
    };

    auto Xh = getSpace( Vh );
#if 1
    tic();
    auto l = blockform1( Xh, solve::strategy::monolithic, backend() );
    l(0_c) = integrate( _range=elements(mesh),
                        _expr=inner(f,id(v)));
    l(0_c) += integrate( _range=boundaryfaces(mesh),
                            _expr=inner(g,id(v)));

    toc("l");
    tic();
    auto a = blockform2( Xh, solve::strategy::monolithic, backend() );
    a(0_c, 0_c) = integrate(_range=elements(mesh),
                            _expr=inner( sigmat, grad(v) ) );
    toc("a");

    if constexpr ( c_impl == ConstraintImplType::LagrangeMultiplier )
    {
        auto mv = Xh[1]->element();
        auto m = Xh[2]->element();
        // translation
        l(1_c) += integrate(_range=elements(mesh),
                            _expr=inner(0.*tau,id(mv)));
        // rotation
#if FEELPP_DIM == 2
        l(2_c) += integrate(_range=rigid_body_domain,
                            _expr=cst(0.)*id(m));
#else
        l(2_c) += integrate(_range=rigid_body_domain,
                            _expr=inner(0*oneX(),id(m)));
#endif

        tic();
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
    }
    else if constexpr ( c_impl == ConstraintImplType::Nodal )
    {
        auto [trans,rot] = rigidBodyMotion( mesh, omega, tau, rigid_body_domain_name );

    }
    else
    {
        //static_assert( false_v<ConstraintImplType>, "ConstraintImplType not supported" );
    }

    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
    {
        tic();
        auto U=Xh.element();
        a.solve(_rhs=l,_solution=U);
        u=U(0_c);
        //m=U(1_c);
        toc("a.solve");

        auto rbtrans=integrate(_range=elements(mesh),_expr=idv(u)).evaluate();
        auto rbrot = integrate(_range=elements(mesh),_expr=curlv(u)).evaluate();
        std::cout << fmt::format( "rbtrans={}, rbrot={}", rbtrans, rbrot ) << std::endl;

    }
#endif
    return u;
}
}
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
                    _about=about(_name="qs_elasticity_pure_traction",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

        tic();
        auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
        toc("loadMesh");
        auto tau = expr<FEELPP_DIM,1>( soption(_name="functions.tau"), "tau" );
    #if FEELPP_DIM == 2
        auto omega = expr( soption(_name="functions.omega"), "omega" );
    #else
        auto omega = expr<3,1>( soption(_name="functions.omega"), "omega" );
    #endif
        auto rigid_body_domain_name = soption(_name="domain");

        auto Vh = Pchv<1>( mesh );
        auto [trans,rot] = rigidBodyMotion( Vh, tau, omega, rigid_body_domain_name );

        //auto [ u, v, w ] = elasticity<ConstraintImplType::Nodal>( mesh );
        auto ulag = elasticity<ConstraintImplType::LagrangeMultiplier>( Vh );

        tic();
        auto e = exporter( _mesh=mesh );
        e->addRegions();
        e->add( "ulag", ulag );
        //e->add( "ulag", u+trans+rot );
        e->add( "rotation", rot );
        e->add( "translation", trans );
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
