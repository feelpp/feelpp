/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
            Daniele Prada <daniele.prada85@gmail.com>
 Date: 16 March 2016
 
 Copyright (C) 2015-present Feel++ Consortium
 
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
#ifndef FEELPP_HDG_DARCY_HPP
#define FEELPP_HDG_DARCY_HPP 1

#include <functional>
#include <tuple>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelts/newmark.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelalg/vectorblock.hpp>

namespace Feel
{

namespace HDG
{

po::options_description
darcy_options(std::string prefix)
{
    po::options_description darcyoptions( "Options for solving Darcy-like problems with HDG" );
    darcyoptions.add_options()
        ( prefixvm( prefix, "filename").c_str(), Feel::po::value<std::string>()->default_value( "" ), "json file describing model properties" )
       ( prefixvm( prefix, "tau-cst").c_str(), Feel::po::value<double>()->default_value( 1. ), "stabilization constant" )
        /*
         The order of the stabilization function greatly affects the order of convergence of
         the HDG method (see Cockburn, Dong, Guzman, Restelli, Sacco,
         "A hybridizable discontinuous Galerkin method for steady-state
         convection-diffusion-reaction problems").
         */
       ( prefixvm( prefix, "tau-order").c_str(), Feel::po::value<int>()->default_value( 0 ), "order of the stabilization function" )
        ( prefixvm( prefix, "verbose").c_str(), Feel::po::value<bool>()->default_value( true ), "verbose" )
        ;
    return darcyoptions.add( backend_options(prefix) ).add( ts_options(prefix) );
}

template<typename PotSpaceType>
class Darcy
{

public:
    using self_type = Darcy<PotSpaceType>;
    using potential_space_type = typename mpl::if_<is_shared_ptr<PotSpaceType>,
                                                  mpl::identity<typename PotSpaceType::element_type>,
                                                  mpl::identity<PotSpaceType>>::type::type;
    using potential_space_ptrtype = std::shared_ptr<potential_space_type>;
    using potential_type = typename potential_space_type::element_type;
    using potential_ptrtype = typename potential_space_type::element_ptrtype;

    static constexpr int order = potential_space_type::nOrder;

    // mesh
    using mesh_type = typename potential_space_type::mesh_type;
    using mesh_ptrtype = typename potential_space_type::mesh_ptrtype;

    // Flux
    using flux_space_type = Pdhv_type<mesh_type, order>;
    using flux_space_ptrtype = Pdhv_ptrtype<mesh_type, order>;
    using flux_type = typename flux_space_type::element_type;
    using flux_ptrtype = typename flux_space_type::element_ptrtype;

    /*
     Lagrange multiplier
     */

    static constexpr int dim = mesh_type::nDim;
    static constexpr int convexOrder = mesh_type::shape_type::nOrder;

    // mesh for the Lagrange multiplier
    typedef mpl::if_< mpl::bool_< mesh_type::shape_type::is_simplex >, Simplex<dim-1, convexOrder, dim> , Hypercube<dim-1, convexOrder, dim> > face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef std::shared_ptr<face_mesh_type> face_mesh_ptrtype;

    using multiplier_space_type = Pdh_type<face_mesh_type, order>;
    using multiplier_space_ptrtype = Pdh_ptrtype<face_mesh_type, order>;
    using multiplier_type = typename multiplier_space_type::element_type;
    using multiplier_ptrtype = typename multiplier_space_type::element_ptrtype;

    // exporter
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_type>>;

    Darcy() = delete;
    Darcy( std::string name, potential_space_ptrtype Wh_ptr_t );

    void setStabilizationConstant( double t ) { tau_constant = t };
    void setStabilizationOrder( int  o ) { tau_order = o };

    void setPermeabilityCoefficient( scalar_field_expression<2> e ) { permeability = e; }
    void setSource( scalar_field_expression<2> e ) { source = e; }

    void solve();

    flux_type const& flux() const { return *up; }
    potential_type const& potential() const { return *pp; }
    multiplier_type const& multiplier() const { return *phatp; }

    void exportResults();

private:
    void initModel();
    void addDirichlet();
    void addNeumann();

private:
    std::string name;
    ModelProperties props;

    std::string model;
    bool verbose;

    // meshes
    mesh_ptrtype mesh;
    face_mesh_ptrtype face_mesh;

    // FE spaces
    flux_space_ptrtype Vh;
    potential_space_ptrtype Wh;
    multiplier_space_ptrtype Mh;

    // FE elements
    flux_ptrtype up;
    potential_ptrtype pp;
    multiplier_ptrtype phatp;

    flux_type u, v;
    potential_type p, w;
    multiplier_type phat, l;

    // Big matrix and vectors associated to the three fields formulation
    sparse_matrix_ptrtype A;
    vector_ptrtype F, U;

    // Auxiliary block vector for storing the solution
    BlocksBaseVector<double> hdg_sol(3);

    // Stabilization function
    double tau_constant;
    int tau_order;

    // Boundary conditions
    map_scalar_field<2> dirichlet_conditions;
    map_vector_field<dim,1,2> neumann_conditions;

    scalar_field_expression<2> permeability;
    scalar_field_expression<2> source;

    backend_ptrtype M_backend;

    exporter_ptrtype e;
};

template<typename PotSpaceType>
Darcy<PotSpaceType>::Darcy( std::string n, potential_space_ptrtype Wh_ptr_t )
    :
    name( n ),
    props( Environment::expand( soption( _name=prefixvm(name,"filename")) ) ),
    model( soption( _name=prefixvm(name,"model")) ),
    verbose( boption(_name=prefixvm(name,"verbose")) ),
    Wh( Wh_ptr_t ),
    mesh( Wh->mesh() ),
    Vh( Pdhv<order>( mesh, true ) ),
    face_mesh( createSubmesh( _mesh=mesh, _range=faces(mesh), _update=0 ) ),
    Mh( Pdh<order>( face_mesh, true ) ),
    up( Vh->elementPtr() ),
    pp( Wh->elementPtr() ),
    phatp( Mh->elementPtr() ),
    u( Vh->element() ),
    v( Vh->element() ),
    p( Wh->element() ),
    w( Wh->element() ),
    phat( Mh->element() ),
    l( Mh->element() ),
    tau_constant( doption("tau-cst") ),
    tau_order( ioption("tau-order") ),
    M_backend( backend( _name=n ) ),
    e( exporter( _mesh=mesh, _prefix=props.shortName(), _geo="static" ) )
{
    tic();

    //    permeability = expr( soption(prefixvm(name, "permeability")) );
    //    source       = expr( soption(prefixvm(name, "source")) );
    permeability = props.functions()["permeability"];
    source       = props.functions()["source"];

    permeability.setParameterValues( props.parameters().toParameterValues() );
    source.setParameterValues( props.parameters().toParameterValues() );

    dirichlet_conditions = props.boundaryConditions().getScalarFields ( "potential", "Dirichlet" );
    neumann_conditions = props.boundaryConditions().template getVectorFields<dim> ( "flux", "Neumann" );

    this -> initModel();

    toc("Darcy constructor", verbose || FLAGS_v > 0);
}

template<typename PotSpaceType>
void
Darcy<PotSpaceType>::initModel()
{
    // build the big matrix associated to bilinear form over Vh x Wh x Mh
    BlocksBaseGraphCSR hdg_graph(3,3);
    hdg_graph(0,0) = stencil( _test=Vh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,0) = stencil( _test=Wh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,0) = stencil( _test=Mh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(0,1) = stencil( _test=Vh,_trial=Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,1) = stencil( _test=Wh,_trial=Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,1) = stencil( _test=Mh,_trial=Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(0,2) = stencil( _test=Vh,_trial=Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,2) = stencil( _test=Wh,_trial=Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,2) = stencil( _test=Mh,_trial=Mh, _diag_is_nonzero=false, _close=false)->graph();

    A = backend()->newBlockMatrix(_block=hdg_graph);

    BlocksBaseVector<double> hdg_vec(3);
    hdg_vec(0,0) = backend()->newVector( Vh );
    hdg_vec(1,0) = backend()->newVector( Wh );
    hdg_vec(2,0) = backend()->newVector( Mh );
    F = backend()->newBlockVector(_block=hdg_vec, _copy_values=false);

    hdg_sol(0,0) = up;
    hdg_sol(1,0) = pp;
    hdg_sol(2,0) = phatp;
    U = backend()->newBlockVector(_block=hdg_sol, _copy_values=false);

    /*
     * Assemble A and F
     */

    auto lambda = cst(1.)/permeability;

    // Building the RHS

    auto rhs2 = form1( _test=Wh, _vector=F,
                       _rowstart=Vh->nLocalDofWithGhost() );

    rhs2 += integrate(_range=elements(mesh),
                      _expr=source*id(w));

    auto a11 = form2( _trial=Vh, _test=Vh,_matrix=A );
    a11 += integrate(_range=elements(mesh),_expr=(trans(lambda*idt(u))*id(v)) );

    auto a12 = form2( _trial=Wh, _test=Vh,_matrix=A,
                      _rowstart=0, _colstart=Vh->nLocalDofWithGhost() );
    a12 += integrate(_range=elements(mesh),_expr=-(idt(p)*div(v)));

    auto a13 = form2( _trial=Mh, _test=Vh,_matrix=A,
                      _rowstart=0, _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());
    a13 += integrate(_range=internalfaces(mesh),
                     _expr=( idt(phat)*leftface(trans(id(v))*N())+
                             idt(phat)*rightface(trans(id(v))*N())) );
    a13 += integrate(_range=boundaryfaces(mesh),
                     _expr=idt(phat)*trans(id(v))*N());

    auto a21 = form2( _trial=Vh, _test=Wh,_matrix=A,
                      _rowstart=Vh->nLocalDofWithGhost(), _colstart=0);
    a21 += integrate(_range=elements(mesh),_expr=(-grad(w)*idt(u)));
    a21 += integrate(_range=internalfaces(mesh),
                     _expr=( leftface(id(w))*leftfacet(trans(idt(u))*N()) ) );
    a21 += integrate(_range=internalfaces(mesh),
                     _expr=(rightface(id(w))*rightfacet(trans(idt(u))*N())) );
    a21 += integrate(_range=boundaryfaces(mesh),
                     _expr=(id(w)*trans(idt(u))*N()));

    auto a22 = form2( _trial=Wh, _test=Wh,_matrix=A,
                      _rowstart=Vh->nLocalDofWithGhost(), _colstart=Vh->nLocalDofWithGhost() );
    a22 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant *
                     ( leftfacet( pow(h(),tau_order)*idt(p))*leftface(id(w)) +
                       rightfacet( pow(h(),tau_order)*idt(p))*rightface(id(w) )));
    a22 += integrate(_range=boundaryfaces(mesh),
                     _expr=(tau_constant * pow(h(),tau_order)*id(w)*idt(p)));

    auto a23 = form2( _trial=Mh, _test=Wh,_matrix=A,
                      _rowstart=Vh->nLocalDofWithGhost(), _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());
    a23 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * idt(phat) *
                     ( leftface( pow(h(),tau_order)*id(w) )+
                       rightface( pow(h(),tau_order)*id(w) )));
    a23 += integrate(_range=boundaryfaces(mesh),
                     _expr=-tau_constant * idt(phat) * pow(h(),tau_order)*id(w) );

    auto a31 = form2( _trial=Vh, _test=Mh,_matrix=A,
                      _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(), _colstart=0);
    a31 += integrate(_range=internalfaces(mesh),
                     _expr=( id(l)*(leftfacet(trans(idt(u))*N())+
                                    rightfacet(trans(idt(u))*N())) ) );

    auto a32 = form2( _trial=Wh, _test=Mh,_matrix=A,
                      _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(), _colstart=Vh->nLocalDofWithGhost());
    a32 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant * id(l) * ( leftfacet( pow(h(),tau_order)*idt(p) )+
                                                    rightfacet( pow(h(),tau_order)*idt(p) )));

    auto a33 = form2(_trial=Mh, _test=Mh,_matrix=A,
                     _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(), _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());
    a33 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * idt(phat) * id(l) * ( leftface( pow(h(),tau_order) )+
                                                                 rightface( pow(h(),tau_order) )));
}

template<typename PotSpaceType>
void
Darcy<PotSpaceType>::addDirichlet()
{
    tic();

    auto rhs3 = form1( _test=Mh, _vector=F,
                       _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());
    auto a33 = form2(_trial=Mh, _test=Mh,_matrix=A,
                     _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(),
                     _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());

    dirichlet_conditions.setParameterValues( props.parameters().toParameterValues() );
    for ( auto const& d : dirichlet_conditions )
    {
        LOG(INFO) << " - dirichlet condition on " << name(d) << " expr=" << expression(d) << "\n";
        rhs3 += integrate(_range=markedfaces(mesh, markers(d)),
                          _expr=id(l)*expression(d));
        a33 += integrate(_range=markedfaces(mesh, markers(d)),
                         _expr=idt(phat) * id(l) );

    }

    toc("Darcy::addDirichlet", verbose || FLAGS_v > 0);
}

template<typename PotSpaceType>
void
Darcy<PotSpaceType>::addNeumann()
{
    tic();

    auto rhs3 = form1( _test=Mh, _vector=F,
                       _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());
    auto a31 = form2( _trial=Vh, _test=Mh, _matrix=A,
                      _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(),
                      _colstart=0);
    auto a32 = form2( _trial=Wh, _test=Mh, _matrix=A,
                      _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(),
                      _colstart=Vh->nLocalDofWithGhost());
    auto a33 = form2(_trial=Mh, _test=Mh,_matrix=A,
                     _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(),
                     _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());

    neumann_conditions.setParameterValues( props.parameters().toParameterValues() );
    for ( auto const& d : neumann_conditions )
    {
        LOG(INFO) << " - Neumann condition on " << name(d) << " expr=" << expression(d) << "\n";
        rhs3 += integrate(_range=markedfaces(mesh, markers(d)),
                          _expr=id(l)*expression(d));
        a31 += integrate(_range=markedfaces(mesh, markers(d)),
                         _expr=( id(l)*(trans(idt(u))*N()) ));
        a32 += integrate(_range=markedfaces(mesh, markers(d)),
                         _expr=tau_constant * id(l) * ( pow(h(),tau_order)*idt(p) ) );
        a33 += integrate(_range=markedfaces(mesh, markers(d)),
                         _expr=-tau_constant * idt(phat) * id(l) * ( pow(h(),tau_order) ) );
    }

    toc("Darcy::addNeumann", verbose || FLAGS_v > 0);
}

template<typename PotSpaceType>
void
Darcy<PotSpaceType>::solve()
{
    tic();
    this -> addDirichlet();
    this -> addNeumann();
    backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );
    hdg_sol.localize(U);
    toc("Darcy::solve", verbose || FLAGS_v > 0 );
}

template<typename PotSpaceType>
void
Darcy<PotSpaceType>::exportResults()
{
    tic();

    for ( auto const& o : props.postProcess()["Fields"] )
    {
        if ( o == "flux" )
            e->add( "flux", *up );
        if ( o == "potential" )
            e->add( "potential", *pp );
        if ( o == "multiplier" )
            e->add( "multiplier", *phatp)
    }
    e->save();

    toc("Darcy::exportResults", verbose);
}

} // HDG

} // Feel

#endif
