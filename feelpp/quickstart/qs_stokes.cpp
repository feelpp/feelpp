/*
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2014-2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/productfunctionspaces.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feeldiscr/check.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelvf/vf_eval.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelvf/blockforms.hpp>
#include "discretization.hpp"

///[stokes]
template<typename VelocitySpacePtrType, typename PressureSpacePtrType>
int
stokes(VelocitySpacePtrType velocitySpace, PressureSpacePtrType pressureSpace)
{
    using namespace Feel;
    using namespace boost::hana::literals;

    // tag::mesh_space[]
    tic();
    auto mesh = velocitySpace->mesh();
    auto Vh = productFunctionSpacesPtr( velocitySpace, pressureSpace );
    auto U = Vh->element("U");
    auto& u = U.template element<0>();
    auto& p = U.template element<1>();
    auto v = U.template element<0>();
    auto q = U.template element<1>();
    auto mu = doption(_name="mu");
    auto f = vf::expr<FEELPP_DIM,1>( soption(_name="functions.f") );
    auto thechecker = checker(_name="qs_stokes",_solution_key=soption("checker.solution"));
    auto solution = vf::expr<FEELPP_DIM,1>( thechecker.check()? thechecker.solution() : soption(_name="functions.g") );
    auto g = solution;
    toc("Vh");
    // end::mesh_space[]

    // tag::forms[]
    tic();
    auto l = blockform1( Vh->productSpace(), solve::strategy::monolithic, backend() );
    l( 0_c ) += Feel::integrate(_range=elements(mesh),
                        _expr=vf::inner(f,vf::id(v)));
    toc("l");

    tic();
    auto a = blockform2( Vh->productSpace(), solve::strategy::monolithic, backend() );
    auto deft = vf::sym(vf::gradt(u));
    tic();
    a( 0_c, 0_c ) += Feel::integrate(_range=elements(mesh),
                        _expr=vf::inner( 2*mu*deft, vf::grad(v) ) );
    a( 0_c, 1_c ) += Feel::integrate(_range=elements(mesh),
                         _expr=-vf::idt(p)*vf::div(v) );
    a( 1_c, 0_c ) += Feel::integrate(_range=elements(mesh),
                         _expr=vf::id(q)*vf::divt(u) );
    toc("a");

    if ( mesh->hasAnyMarker({"inlet","Dirichlet"}) )
        a.row( 0_c ) += Feel::on(_range=markedfaces(mesh,{"inlet","Dirichlet"}), _rhs=l( 0_c ), _element=u, _expr=g );
    if ( mesh->hasAnyMarker({"wall","letters"}) )
        a.row( 0_c ) += Feel::on(_range=markedfaces(mesh,{"wall","letters"}), _rhs=l( 0_c ), _element=u, _expr=vf::zero<FEELPP_DIM,1>() );
    toc("a");

    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
    {
        auto const status = a.solve(_rhs=l,_solution=U.productElement());
        if ( !status.isConverged() )
            return 1;
    }
    toc("a.solve");
    // end::forms[]

    // tag::export[]
    tic();
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "uh", u );
    e->add( "ph", p );
    if ( thechecker.check() )
    {
        e->add( "u", solution );
    }
    e->save();
    toc("Exporter");
    // end::export[]

    return check( thechecker, u );
}
///[stokes]
int main(int argc, char**argv )
{
    ///[stokes-env]
    using namespace Feel;
    try
    {
	    po::options_description laplacianoptions( "Stokes options" );
	    laplacianoptions.add_options()
            ( "mu", po::value<double>()->default_value( 1.0 ), "viscosity" )
            ( "space", po::value<std::string>()->default_value( "P2P1" ), "space type: P2P1(default), P1P1, P1P0" )
            ( "order", po::value<int>(), "pressure order (default 1, or 0 for P1P0); velocity order is one higher except for equal-order P1P1" )
            ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" )
            ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
	    	;

	    Environment env( _argc=argc, _argv=argv,
                       _desc=laplacianoptions,
                       _about=about(_name="qs_stokes",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));
        ///[stokes-env]

        tic();
        ///[stokes-mesh]
        auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,Dynamic>>( quickstart::geometryOrder() ));
        quickstart::checkMeshOrder( mesh );
        ///[stokes-mesh]
        toc("loadMesh");

        ///[stokes-space]
        auto const space = soption( "space" );
        if ( space != "P2P1" && space != "P1P1" && space != "P1P0" )
            throw std::invalid_argument( "space must be P2P1, P1P1 or P1P0" );
        bool const discontinuousPressure = space == "P1P0";
        auto const pressureOrder = Environment::vm().count( "order" )
                                       ? quickstart::checkedOrder( "order", discontinuousPressure ? 0 : 1 )
                                       : RuntimeOrder{ discontinuousPressure ? 0 : 1 };
        auto const velocityOrder = RuntimeOrder::checked( pressureOrder.value + ( space == "P1P1" ? 0 : 1 ) );
        auto velocitySpace = Pchv<Dynamic>( mesh, velocityOrder );
        int status;
        if ( discontinuousPressure )
            status = stokes( velocitySpace, Pdh<Dynamic>( mesh, pressureOrder ) );
        else
            status = stokes( velocitySpace, Pch<Dynamic>( mesh, pressureOrder ) );
        ///[stokes-space]
        return status;
    }
    catch( ... ) 
    {
        handleExceptions();
    }
    return 1;
}
