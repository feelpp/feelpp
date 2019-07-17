/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 22 Feb 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_FLUIDMECHANICS_HPP
#define FEELPP_FLUIDMECHANICS_HPP 1

#include <functional>
#include <tuple>
#include <feel/feel.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelts/newmark.hpp>
#include <feel/feelmodels/modelproperties.hpp>

namespace Feel
{
enum FluidMechanicsModel
{
    STOKES = 0,
    NAVIER_STOKES
};
enum FluidMechanicsOptions
{
    FM_STEADY=0x1,
    FM_UNSTEADY=0x2,
    FM_LINEARIZED=0x3,
    FM_FULLY_IMPLICIT=0x4
};


po::options_description
fluidmechanics_options(std::string prefix)
{
    po::options_description fluidmechanicsoptions( "FluidMechanics problem options" );
    fluidmechanicsoptions.add_options()
        ( prefixvm( prefix, "filename").c_str(), Feel::po::value<std::string>()->default_value( "" ), "filename that describes the model" )
        ( prefixvm( prefix, "model").c_str(), Feel::po::value<std::string>()->default_value( "Stokes" ), "Stokes, NavierStokes" )
        ( prefixvm( prefix, "verbose").c_str(), Feel::po::value<bool>()->default_value( true ), "verbose" )
        ( prefixvm( prefix, "postprocess.type").c_str(), Feel::po::value<std::string>()->default_value( "" ),
          "evaluate a field : point, integral" )
        ( prefixvm( prefix, "postprocess.filename").c_str(), Feel::po::value<std::string>()->default_value( "" ),
          "filename with a list of points for pointwise evaluation" )
        ( prefixvm( prefix, "postprocess.markers").c_str(), Feel::po::value<std::vector<std::string>>(),
          "list of markers to evaluate pointwise and/or integral(average/flux)" )
        ( prefixvm( prefix, "postprocess.evaluation.velocity").c_str(), Feel::po::value<bool>()->default_value( 0 ),
          "pointwise evaluation for velocity" )
        ( prefixvm( prefix, "postprocess.evaluation.flowrate").c_str(), Feel::po::value<bool>()->default_value( 0 ),
          "pointwise evaluation for flowrate" )
        ( prefixvm( prefix, "postprocess.evaluation.pressure").c_str(), Feel::po::value<bool>()->default_value( 0 ),
          "pointwise evaluation for pressure" )
        ;
    return fluidmechanicsoptions.add( backend_options(prefix) ).add( ts_options(prefix) );;
}

/**
 * @brief a fluid mechanics class parametrized by an approximation space and a model type
 */
template<typename FluidSpaceType>
class FluidMechanics 
{
    
public:
    
    using self_type = FluidMechanics<FluidSpaceType>;
    using fluid_space_type = decay_type<FluidSpaceType>;
    using fluid_space_ptrtype = std::shared_ptr<fluid_space_type>;
    using fluid_type = typename fluid_space_type::element_type;

    // mesh
    using mesh_type = typename fluid_space_type::mesh_type;
    using mesh_ptrtype = typename fluid_space_type::mesh_ptrtype;

    // velocity
    using velocity_space_type = typename fluid_space_type::template sub_functionspace_type<0>;
    using velocity_space_ptrtype = std::shared_ptr<velocity_space_type>;
    using velocity_type = typename velocity_space_type::element_type;
    using velocity_view_type = typename fluid_type::template sub_element_type<0>;

    // pressure
    using pressure_space_type = typename fluid_space_type::template sub_functionspace_type<1>;
    using pressure_space_ptrtype = std::shared_ptr<pressure_space_type>;
    using pressure_type = typename pressure_space_type::element_type;
    using pressure_view_type = typename fluid_type::template sub_element_type<1>;

#if defined( FEELPP_FLUIDMECHANICS_BOUSSINESQ )
    // temperature
    using temperature_space_type = typename fluid_space_type::template sub_functionspace_type<2>;
    using temperature_space_ptrtype = std::shared_ptr<temperature_space_type>;
    using temperature_type = typename temperature_space_type::element_type;
    using temperature_view_type = typename fluid_type::template sub_element_type<2>;
#endif

    // properties : rho, mu, k, Cp
    using property_space_ptrtype = Pdh_ptrtype<mesh_type,0>;
    using properties = typename Pdh_type<mesh_type,0>::element_type;

    // exporter
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_type>>;

    // time stepping
    using ts_ptrtype = std::shared_ptr<Bdf<fluid_space_type>>;
    
    FluidMechanics() = delete;
    FluidMechanics( std::string name, fluid_space_ptrtype Xh, int _Options = FM_LINEARIZED );

    static constexpr int dim = mesh_type::nDim;

    // @return the time stepping strategy
    ts_ptrtype ts() { return t; }

    void setMaterialProperties()
        {
            for( auto const& mat : props.materials() )
            {
                if ( Environment::isMasterRank() )
                    std::cout << "FluidMechanics: set material " << material(mat) << "(marker = " << marker(mat) << "\n";
                P0Rho.on( _range=markedelements(mesh,marker(mat)), _expr=cst(material(mat).rho()) );
                P0Mu.on( _range=markedelements(mesh,marker(mat)), _expr=cst(material(mat).mu() ) );
            }
        }
    template < typename ExprT >
    void updateRho( Expr<ExprT> const& e )
        {
            P0Rho.on( _range=elements(mesh), _expr=e );
        }
    template < typename ExprT >
    void updateMu(Expr<ExprT> const& e )
        {
            P0Mu.on(_range=elements(mesh), _expr=e );
        }
#if defined( FEELPP_FLUIDMECHANICS_BOUSSINESQ )
    // heat diffusion
    template < typename ExprT >
    void updateK(Expr<ExprT> const& e )
        {
            P0K.on(_range=elements(mesh), _expr=e );
        }
#endif
    void init();

    using SolveData = SolverNonLinear<double>::SolveData;
    SolveData solve();

    fluid_type& solution()  { return U; }
    fluid_type const& solution() const  { return U; }
    velocity_view_type const& velocity() const { return U.template element<0>(); }
    pressure_view_type const& pressure() const { return U.template element<1>(); }
#if defined( FEELPP_FLUIDMECHANICS_BOUSSINESQ )
    temperature_view_type const& temperature() const { return U.template element<2>(); }
#endif
    
    void exportResults();
    
    void updateResidual(const vector_ptrtype& X, vector_ptrtype& R);
    void updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J);

private:
    void initModel( mpl::int_<STOKES> );
    void initModel( mpl::int_<NAVIER_STOKES> );
    SolveData solve( mpl::int_<STOKES> );
    SolveData solve( mpl::int_<NAVIER_STOKES> );
    void addDirichlet( form2_type<fluid_space_type,fluid_space_type>& bf,
                       form1_type<fluid_space_type>& lf );
private:
    ModelProperties props;
    int Options;
    bool verbose;
    fluid_space_ptrtype Xh;
    property_space_ptrtype P0h;
    mesh_ptrtype mesh;
    fluid_type U,V;
    velocity_view_type u,v;
    pressure_view_type p,q;
#if defined( FEELPP_FLUIDMECHANICS_BOUSSINESQ )
    temperature_view_type T,w;
#endif
    properties P0Rho, P0Mu;
#if defined( FEELPP_FLUIDMECHANICS_BOUSSINESQ )
    properties P0K, P0Cp;
#endif
    map_vector_field<dim,1,2> dirichlet_conditions;
    map_scalar_field<2> dx_dirichlet_conditions, dy_dirichlet_conditions, dz_dirichlet_conditions;
    map_matrix_field<dim,dim,2> neumann_conditions;
    vector_field_expression<dim,1,2> gravityForce;

    ts_ptrtype t;
    
    backend_ptrtype M_backend;
    vector_ptrtype Res;
    sparse_matrix_ptrtype Jac;

    form1_type<fluid_space_type> l,lt;
    form2_type<fluid_space_type,fluid_space_type> a,at;
    exporter_ptrtype e,et;
};

template<typename FluidSpaceType>
FluidMechanics<FluidSpaceType>::FluidMechanics( std::string n, fluid_space_ptrtype Wh, int _Options )
    :
    props( Environment::expand( soption( _name=prefixvm(n,"filename")) ) ),
    Options( _Options ),
    verbose( boption(_name=prefixvm(n,"verbose")) ),
    Xh( Wh ),
    P0h( Pdh<0>( Xh->mesh() ) ),
    mesh( Xh->mesh() ),
    U( Xh->element() ),
    V( Xh->element() ),
    u( U.template element<0>() ),
    v( V.template element<0>() ),
    p( U.template element<1>() ),
    q( V.template element<1>() ),
    P0Rho( P0h->element() ),
    P0Mu( P0h->element() ),
    t( bdf( _space=Xh, _name=n,_rank_proc_in_files_name=true ) ),
    M_backend( backend( _name=n ) ),
    l(form1( _test=Xh)),
    lt(form1( _test=Xh)),
    a(form2( _trial=Xh, _test=Xh)),
    at(form2( _trial=Xh, _test=Xh)),
    e ( exporter( _mesh=mesh, _prefix=props.name(), _geo="static" ) ),
    et ( exporter( _mesh=mesh, _prefix=props.name()+"+t", _geo="static" ) )
{
    tic();
    tic();
    dirichlet_conditions = props.boundaryConditions().template getVectorFields<dim> ( "velocity", "Dirichlet" );
    dx_dirichlet_conditions = props.boundaryConditions().getScalarFields ( "velocity_x", "Dirichlet" );
    dy_dirichlet_conditions = props.boundaryConditions().getScalarFields ( "velocity_y", "Dirichlet" );
    dz_dirichlet_conditions = props.boundaryConditions().getScalarFields ( "velocity_z", "Dirichlet" );
    
    neumann_conditions = props.boundaryConditions().template getMatrixFields<dim> ( "velocity", "Neumann" );
    toc("FM boundary conditions",verbose||FLAGS_v>0);
    tic();
    this->setMaterialProperties();
    toc("FM material properties",verbose||FLAGS_v>0);
    

    if ( props.model() == "Stokes" )
        this->initModel( mpl::int_<STOKES>() );
    else if ( props.model() == "Navier_Stokes" )
        this->initModel( mpl::int_<NAVIER_STOKES>() );
    
    toc("FluidMechanics constructor", verbose || FLAGS_v > 0 );
    
    
}
template<typename FluidSpaceType>
void
FluidMechanics<FluidSpaceType>::initModel( mpl::int_<STOKES> )
{
    auto deft = sym(gradt( u ));
    auto def = grad( v );
    
    //l = integrate( _range=elements(mesh), _expr=idv(P0Rho)*trans(gravityForce)*id(v) );

    a = integrate( _range=elements( mesh ), _expr=idv(P0Mu)*inner( deft,def ) );
    a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) - divt( u )*id( q ) );

    if ( Options & FM_UNSTEADY )
    {
        LOG(INFO) << "add mass time term";
        a +=integrate( _range=elements( mesh ), _expr=idv(P0Rho)*trans(idt(u))*id(v)*t->polyDerivCoefficient(0) );
    }

    neumann_conditions.setParameterValues( props.parameters().toParameterValues() );
    for( auto const& d : neumann_conditions )
    {
        a += integrate( _range=markedfaces( mesh, marker(d) ), _expr=trans(expression(d)*N())*id(v) );
    }
    
    
}

template<typename FluidSpaceType>
void
FluidMechanics<FluidSpaceType>::initModel( mpl::int_<NAVIER_STOKES> )
{
    initModel( mpl::int_<STOKES>() );
    at.zero();
    at += a;
    lt.zero();
    lt += l;
}

template<typename FluidSpaceType>
void
FluidMechanics<FluidSpaceType>::addDirichlet( form2_type<fluid_space_type, fluid_space_type>& bf,
                                                              form1_type<fluid_space_type>& lf)  
{
    tic();
    dirichlet_conditions.setParameterValues( props.parameters().toParameterValues() );
    dirichlet_conditions.setParameterValues( { {"t",t->time()}}  );
    for( auto const& d : dirichlet_conditions )
    {
        LOG(INFO) << " - dirichlet condition on " << marker(d) << " expr=" << expression(d) << "\n";
        bf += on( _range=markedfaces( mesh, marker(d)), _rhs=lf, _element=u,  _expr=expression(d) );
    }
    toc("FluidMechanics::addDirichlet", verbose || FLAGS_v > 0);
}
template<typename FluidSpaceType>
void
FluidMechanics<FluidSpaceType>::init()
{
    tic();
    // start or restart
    if ( !t->isRestart() )
    {
        t->start();
        //svk.setInitialGuess();
    }
    else
    {
        double ti = t->restart();
        u = t->prior();
        if ( e->doExport() )
            e->restart(ti);
    }
    toc("FluidMechanics::init", verbose || FLAGS_v > 0);
}

template<typename FluidSpaceType>
void
FluidMechanics<FluidSpaceType>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    u = *X;
    

}
template<typename FluidSpaceType>
void
FluidMechanics<FluidSpaceType>::updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    u = *X;


}

template<typename FluidSpaceType>
typename FluidMechanics<FluidSpaceType>::SolveData
FluidMechanics<FluidSpaceType>::solve()
{
    tic();
    SolveData d;
    if ( props.model() == "Stokes" )
        d = solve( mpl::int_<STOKES>() );
    else if ( props.model() == "Navier_Stokes" )
        d = solve( mpl::int_<NAVIER_STOKES>() );
        
    toc("FluidMechanics::solve", verbose || FLAGS_v > 0 );
    return d;
}
template<typename FluidSpaceType>
typename FluidMechanics<FluidSpaceType>::SolveData
FluidMechanics<FluidSpaceType>::solve( mpl::int_<NAVIER_STOKES>)
{
    tic();
    LOG(INFO) << "Solving for NS...\n";
    if ( Options & FM_LINEARIZED )
    {
        LOG(INFO) << " - linearized...\n";
        auto bdf_poly = t->polyDeriv();
        auto rhsu =  bdf_poly.template element<0>();
        lt = integrate( _range=elements(mesh), _expr=idv(P0Rho)*(trans(idv(rhsu))*id(u) ) );
        lt += l;
        
        auto extrap = t->poly();
        auto extrapu = extrap.template element<0>();
        at = integrate( _range=elements(mesh), _expr=idv(P0Rho)*trans(gradt(u)*idv(extrapu))*id(v) );
        at += a;
        this->addDirichlet( at, lt );
    }
    toc("FluidMechanics::solve<NavierStokes> assembly", verbose || FLAGS_v > 0 );
    tic();
    auto d = at.solve( _rhs=lt, _solution=U );
    toc("FluidMechanics::solve<NavierStokes> solve", verbose || FLAGS_v > 0 );
    return d;
}
template<typename FluidSpaceType>
typename FluidMechanics<FluidSpaceType>::SolveData
FluidMechanics<FluidSpaceType>::solve( mpl::int_<STOKES>)
{
    tic();
    this->addDirichlet( a, l );
    toc("FluidMechanics::solve<Stokes> dirichlet", verbose || FLAGS_v > 0 );
    tic();
    auto d = a.solve( _rhs=l, _solution=U );
    toc("FluidMechanics::solve<Stokes> solve", verbose || FLAGS_v > 0 );
    return d;
}
template<typename FluidSpaceType>
void
FluidMechanics<FluidSpaceType>::exportResults()
{
    tic();

    if ( t->isSteady() )
    {
        for ( auto const& o : props.postProcess()["Fields"] )
        {
            if ( o == "velocity" )
                e->add( "velocity", u );
            if ( o == "pressure" )
                e->add( "pressure", p );
        }
#if defined( FEELPP_FLUIDMECHANICS_BOUSSINESQ )
        e->add( "temperature", T );
#endif
        e->save();
    }
    else
    {
        
        //e->step(t->time())->add( "rho", P0Rho );
        //e->step(t->time())->add( "mu", P0Mu );
        for ( auto const& o : props.postProcess()["Force"] )
        {
            auto F = integrate( _range=markedfaces(mesh,o),
                                _expr=-idv(p)*N()+2*idv(P0Mu)*sym(gradv(u))*N() ).evaluate();
            et->step(t->time())->addScalar( "Fx", F(0,0) );
            et->step(t->time())->addScalar( "Fy", F(1,0) );
            if ( is_3d<mesh_type>::value )
                et->step(t->time())->addScalar( "Fz", F(2,0) );
        }
        for ( auto const& o : props.postProcess()["Fields"] )
        {
            if ( o == "velocity" )
                et->step(t->time())->add( "velocity", u );
            if ( o == "pressure" )
                et->step(t->time())->add( "pressure", p );
            if ( o == "vorticity" )
            {
                q.on(_range=elements(mesh), _expr=curlzv(v) );
                et->step(t->time())->add( "vorticity", q );
            }
            if ( o == "cellRe" )
            {
                q.on(_range=elements(mesh), _expr=idv(P0Rho)*norm2(idv(u))*h()/(2*idv(P0Mu)) );
                et->step(t->time())->add( "cellRe", q );
            }
        }
#if defined( FEELPP_FLUIDMECHANICS_BOUSSINESQ )
        et->step(t->time())->add( "temperature", T );
#endif
        et->save();
    }

    toc("FluidMechanics::exportResults", verbose);
}
}

#endif
