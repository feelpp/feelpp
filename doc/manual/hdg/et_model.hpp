/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel++ library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Daniele Prada <daniele.prada85@gmail.com>
       Date: 2016-02-10

  Copyright (C) 2016 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
//#include <feel/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelts/ts.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelmesh/complement.hpp>

#include <boost/optional.hpp>
#include <boost/algorithm/string.hpp>

namespace Feel {

inline
po::options_description
makeETHDGOptions()
{
    po::options_description testhdivoptions( "Electro-Thermal Model options" );
    testhdivoptions.add_options()
        ( "Tw", po::value<double>()->default_value( 1 ), "water temperature" )
        ( "tau_constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "tau_order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( "picard.itol", po::value<double>()->default_value( 1e-4 ), "tolerance" )
        ( "picard.itmax", po::value<int>()->default_value( 10 ), "iterations max" )
        // ("Kp", po::value<double>()->default_value(100.), "PID proportional coefficient")
        // ("Ki", po::value<double>()->default_value(0.), "PID integral coefficient")
        // ("Kd", po::value<double>()->default_value(0.), "PID derivative coefficient")
        ("model_json", po::value<std::string>()->default_value("model.json"), "json file for the model")
        ;
    return testhdivoptions;
}

inline
po::options_description
makeETHDGLibOptions()
{
    po::options_description libOptions( "Lib options" );
    libOptions.add( backend_options( "potential" ) );
    libOptions.add( backend_options( "temperature" ) );
    return libOptions.add( feel_options() );
}

template<int Dim, int OrderP>
class ElectroThermal
{
public:

    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;
    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,1,Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef boost::shared_ptr<face_mesh_type> face_mesh_ptrtype;

    // Ch
    using Ch_t = Pch_type<mesh_type,0>;
    using Ch_ptr_t = Pch_ptrtype<mesh_type,0>;
    using Ch_element_t = typename Ch_t::element_type;
    using Ch_element_ptr_t = typename Ch_t::element_ptrtype;
    // Xh
    using Xh_t = Pch_type<mesh_type,OrderP>;
    using Xh_ptr_t = Pch_ptrtype<mesh_type,OrderP>;
    using Xh_element_t = typename Xh_t::element_type;
    using Xh_element_ptr_t = typename Xh_t::element_ptrtype;
    // Wh
    using Wh_t = Pdh_type<mesh_type,OrderP>;
    using Wh_ptr_t = Pdh_ptrtype<mesh_type,OrderP>;
    using Wh_element_t = typename Wh_t::element_type;
    using Wh_element_ptr_t = typename Wh_t::element_ptrtype;
    // Vh
    using Vh_t = Pdhv_type<mesh_type,OrderP>;
    using Vh_ptr_t = Pdhv_ptrtype<mesh_type,OrderP>;
    using Vh_element_t = typename Vh_t::element_type;
    using Vh_element_ptr_t = typename Vh_t::element_ptrtype;
    // Mh
    using Mh_t = Pdh_type<face_mesh_type,OrderP>;
    using Mh_ptr_t = Pdh_ptrtype<face_mesh_type,OrderP>;
    using Mh_element_t = typename Mh_t::element_type;
    using Mh_element_ptr_t = typename Mh_t::element_ptrtype;

    //! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

private:
    model_prop_ptrtype M_modelProperties;

    Vh_ptr_t M_Vh; // current
    Wh_ptr_t M_Wh; // potential
    Mh_ptr_t M_Mh;
    Xh_ptr_t M_Xh; // temperature
    Ch_ptr_t M_Ch; // constant

    backend_ptrtype M_pBackend;
    backend_ptrtype M_tBackend;
    sparse_matrix_ptrtype M_A;
    sparse_matrix_ptrtype M_A_cst;
    vector_ptrtype M_F;
    BlocksBaseVector<double> M_hdg_sol;
    vector_ptrtype M_U;

    Vh_element_ptr_t M_up;
    Wh_element_ptr_t M_pp;
    Xh_element_ptr_t M_Tp;

    int M_tau_order;

    bool M_integralCondition;
    bool M_isPicard;

    void initGraphs();
    void initGraphsWithIntegralCond();
    void assembleACst();
    void assembleA( int iter = -1 );
    void assembleF( int iter = -1 );
    void solveT( int iter = -1 );

public:
    void init();
    void solve();
    void exportResults();
};

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::init()
{
    M_modelProperties = std::make_shared<model_prop_type>( Environment::expand( soption("model_json") ) );
    if ( boost::icontains(M_modelProperties->model(), "integral" ) )
        M_integralCondition = true;
    else
        M_integralCondition = false;
    if ( boost::icontains(M_modelProperties->model(),"hdg-picard") )
        M_isPicard = true;
    else
        M_isPicard = false;

    cout << "Model : " << M_modelProperties->model()
         << " using ";
    if ( M_integralCondition )
        cout << "integral condition on the current and ";
    if ( M_isPicard )
        cout << "Picard algorithm" << std::endl;
    else
        cout << "linear case" << std::endl;

    int proc_rank = Environment::worldComm().globalRank();

    M_tau_order = ioption("tau_order");

    tic();
    auto mesh = loadMesh( new mesh_type);
    auto complement_integral_bdy = complement(faces(mesh),[&mesh]( auto const& e ) { return e.marker().value() == mesh->markerName( "V1" ); });
    auto face_mesh = createSubmesh( mesh, complement_integral_bdy, EXTRACTION_KEEP_MESH_RELATION, 0 );
    auto face_mesh_bottom = createSubmesh( mesh, markedfaces(mesh,"V1"), EXTRACTION_KEEP_MESH_RELATION, 0 );

    toc("mesh",true);

    // ****** Hybrid-mixed formulation ******
    // We treat Vh, Wh, and Mh separately
    tic();

    M_Vh = Pdhv<OrderP>( mesh, true );
    M_Wh = Pdh<OrderP>( mesh, true );
    M_Mh = Pdh<OrderP>( face_mesh, true );
    M_Xh = Pch<OrderP>( mesh );
    M_Ch = Pch<0>( mesh );

    toc("spaces",true);

    cout << "Vh<" << OrderP << "> : " << M_Vh->nDof() << std::endl
         << "Wh<" << OrderP << "> : " << M_Wh->nDof() << std::endl
         << "Mh<" << OrderP << "> : " << M_Mh->nDof() << std::endl
         << "Xh<" << OrderP << "> : " << M_Xh->nDof() << std::endl;
    if ( M_integralCondition )
        cout << "Ch<" << 0 << "> : " << M_Ch->nDof() << std::endl;

    M_pBackend = backend( _name="potential", _rebuild=true);
    M_tBackend = backend( _name="temperature", _rebuild=true);

    M_up = M_Vh->elementPtr( "u" );
    M_pp = M_Wh->elementPtr( "p" );
    M_Tp = M_Xh->elementPtr( "T" );

    tic();
    if ( M_integralCondition )
        this->initGraphsWithIntegralCond();
    else
        this->initGraphs();
    toc("matrices",true);
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::solve()
{
    assembleACst();
    if ( boost::icontains(M_modelProperties->model(),"hdg-picard") )
    {
        int itmax = ioption( "picard.itmax" );
        double tol = doption( "picard.itol" );

        auto mesh = M_Vh->mesh();

        auto po = M_Wh->element();
        auto To = M_Xh->element();

        // picard loop
        cout << "  #iteration incrp incrT current" << std::endl;
        int it = 0;
        double incrp, incrT;
        do
        {
            assembleA(it);
            assembleF(it);
            M_pBackend->solve( _matrix=M_A, _rhs=M_F, _solution=M_U );
            M_hdg_sol.localize(M_U);

            solveT(it);

            incrp = normL2( _range=elements(mesh), _expr=idv(*M_pp)-idv(po) );
            incrT = normL2( _range=elements(mesh), _expr=idv(*M_Tp)-idv(To) );
            po = *M_pp;
            To=*M_Tp;

            // compute current
            double I = integrate(_range=markedfaces(mesh,"bottom"),
                                 _expr=inner(idv(*M_up),N())).evaluate()(0,0);
            cout << "  picard #" << it << " " << incrp << " " << incrT << " " << I << std::endl;
        } while ( ( incrp > tol || incrT > tol ) && ( ++it < itmax ) );

    }
    else {
        assembleA();
        assembleF();
        M_pBackend->solve( _matrix=M_A, _rhs=M_F, _solution=M_U );
        M_hdg_sol.localize(M_U);
        solveT();
    }
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::initGraphs()
{
    auto phatp = M_Mh->elementPtr( "phat" );

    BlocksBaseGraphCSR hdg_graph(3,3);
    hdg_graph(0,0) = stencil( _test=M_Vh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,0) = stencil( _test=M_Wh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,0) = stencil( _test=M_Mh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();

    hdg_graph(0,1) = stencil( _test=M_Vh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,1) = stencil( _test=M_Wh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,1) = stencil( _test=M_Mh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();

    hdg_graph(0,2) = stencil( _test=M_Vh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,2) = stencil( _test=M_Wh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,2) = stencil( _test=M_Mh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();

    M_A = M_pBackend->newBlockMatrix(_block=hdg_graph);
    M_A_cst = M_pBackend->newBlockMatrix(_block=hdg_graph);

    BlocksBaseVector<double> hdg_vec(3);
    hdg_vec(0,0) = M_pBackend->newVector( M_Vh );
    hdg_vec(1,0) = M_pBackend->newVector( M_Wh );
    hdg_vec(2,0) = M_pBackend->newVector( M_Mh );
    M_F = M_pBackend->newBlockVector(_block=hdg_vec, _copy_values=false);

    M_hdg_sol = BlocksBaseVector<double>(3);
    M_hdg_sol(0,0) = M_up;
    M_hdg_sol(1,0) = M_pp;
    M_hdg_sol(2,0) = phatp;
    M_U = M_pBackend->newBlockVector(_block=M_hdg_sol, _copy_values=false);
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::initGraphsWithIntegralCond()
{
    auto phatp = M_Mh->elementPtr( "phat" );
    auto mup = M_Ch->elementPtr( "c1" );

    BlocksBaseGraphCSR hdg_graph(4,4);
    hdg_graph(0,0) = stencil( _test=M_Vh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,0) = stencil( _test=M_Wh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,0) = stencil( _test=M_Mh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,0) = stencil( _test=M_Ch,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();

    hdg_graph(0,1) = stencil( _test=M_Vh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,1) = stencil( _test=M_Wh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,1) = stencil( _test=M_Mh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,1) = stencil( _test=M_Ch,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();

    hdg_graph(0,2) = stencil( _test=M_Vh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,2) = stencil( _test=M_Wh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,2) = stencil( _test=M_Mh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,2) = stencil( _test=M_Ch,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();

    hdg_graph(0,3) = stencil( _test=M_Vh,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,3) = stencil( _test=M_Wh,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,3) = stencil( _test=M_Mh,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,3) = stencil( _test=M_Ch,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();

    M_A = M_pBackend->newBlockMatrix(_block=hdg_graph);

    BlocksBaseVector<double> hdg_vec(4);
    hdg_vec(0,0) = M_pBackend->newVector( M_Vh );
    hdg_vec(1,0) = M_pBackend->newVector( M_Wh );
    hdg_vec(2,0) = M_pBackend->newVector( M_Mh );
    hdg_vec(3,0) = M_pBackend->newVector( M_Ch );
    M_F = M_pBackend->newBlockVector(_block=hdg_vec, _copy_values=false);

    M_hdg_sol = BlocksBaseVector<double>(4);
    M_hdg_sol(0,0) = M_up;
    M_hdg_sol(1,0) = M_pp;
    M_hdg_sol(2,0) = phatp;
    M_hdg_sol(3,0) = mup;
    M_U = M_pBackend->newBlockVector(_block=M_hdg_sol, _copy_values=false);
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::assembleA( int iter )
{
    M_A->zero();
    M_A->addMatrix(1,M_A_cst);

    auto mesh = M_Vh->mesh();

    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );

    auto a11 = form2( _trial=M_Vh, _test=M_Vh,_matrix=M_A );

    for( auto const& pairMat : M_modelProperties->materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        if ( M_isPicard && iter > 0)
        {
            auto alpha = material.getDouble("alpha");
            auto T0 = material.getDouble("T0");
            auto sigma0 = material.getDouble("sigma0");
            auto sigma = material.getScalar("sigma", "T", idv(M_Tp));
            sigma.setParameterValues({{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});

            // (sigma^-1 j, v)
            a11 += integrate(_range=markedelements(mesh,marker), _expr=(trans(idt(u))*id(v))/sigma );
        }
        else {
            // (sigma^-1 j, v)
            a11 += integrate(_range=markedelements(mesh,marker), _expr=(trans(idt(u))*id(v))/material.getScalar("sigma0") );
        }
    }
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::assembleACst()
{
    M_A_cst->zero();
    auto mesh = M_Vh->mesh();

    // stabilisation parameter
    auto tau_constant = cst(doption("tau_constant"));
    auto Tw = M_modelProperties->parameters()["Tw"].value();

    auto T = M_Xh->element( "T" );
    auto q = M_Xh->element( "q" );

    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    auto p = M_Wh->element( "p" );
    auto w = M_Wh->element( "w" );
    auto nu = M_Ch->element( "nu" );
    auto phat = M_Mh->element( "phat" );
    auto l = M_Mh->element( "lambda" );


    // Building the LHS
    auto a12 = form2( _trial=M_Wh, _test=M_Vh,_matrix=M_A_cst,
                      _rowstart=0,
                      _colstart=1 );
    auto a13 = form2( _trial=M_Mh, _test=M_Vh,_matrix=M_A_cst,
                      _rowstart=0,
                      _colstart=2);
    auto a21 = form2( _trial=M_Vh, _test=M_Wh,_matrix=M_A_cst,
                      _rowstart=1,
                      _colstart=0);
    auto a22 = form2( _trial=M_Wh, _test=M_Wh,_matrix=M_A_cst,
                      _rowstart=1,
                      _colstart=1 );
    auto a23 = form2( _trial=M_Mh, _test=M_Wh,_matrix=M_A_cst,
                      _rowstart=1,
                      _colstart=2);
    auto a31 = form2( _trial=M_Vh, _test=M_Mh,_matrix=M_A_cst,
                      _rowstart=2,
                      _colstart=0);
    auto a32 = form2( _trial=M_Wh, _test=M_Mh,_matrix=M_A_cst,
                      _rowstart=2,
                      _colstart=1);
    auto a33 = form2(_trial=M_Mh, _test=M_Mh,_matrix=M_A_cst,
                     _rowstart=2,
                     _colstart=2);
    auto a55 = form2(_trial=M_Xh, _test=M_Xh,_matrix=M_A_cst,
                     _rowstart=M_integralCondition?4:3,
                     _colstart=M_integralCondition?4:3);

    // -(p,div(v))
    a12 += integrate(_range=elements(mesh),_expr=-(idt(p)*div(v)));

    // <phat,v.n>_Gamma
    a13 += integrate(_range=internalfaces(mesh),
                     _expr=( idt(phat)*leftface(trans(id(v))*N())+
                             idt(phat)*rightface(trans(id(v))*N())) );
    a13 += integrate(_range=boundaryfaces(mesh),
                     _expr=idt(phat)*trans(id(v))*N());

    // -(j, grad(w))
    a21 += integrate(_range=elements(mesh),_expr=(-grad(w)*idt(u)));
    // <j.n,w>_Gamma
    a21 += integrate(_range=internalfaces(mesh),
                     _expr=( leftface(id(w))*leftfacet(trans(idt(u))*N()) ) );
    a21 += integrate(_range=internalfaces(mesh),
                     _expr=(rightface(id(w))*rightfacet(trans(idt(u))*N())) );
    a21 += integrate(_range=boundaryfaces(mesh),
                     _expr=(id(w)*trans(idt(u))*N()));

    // <tau p, w>_Gamma
    a22 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant *
                     ( leftfacet( pow(h(),M_tau_order)*idt(p))*leftface(id(w)) +
                       rightfacet( pow(h(),M_tau_order)*idt(p))*rightface(id(w) )));
    a22 += integrate(_range=boundaryfaces(mesh),
                     _expr=(tau_constant * pow(h(),M_tau_order)*id(w)*idt(p)));

    // <-tau phat, w>_Gamma
    a23 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * idt(phat) *
                     ( leftface( pow(h(),M_tau_order)*id(w) )+
                       rightface( pow(h(),M_tau_order)*id(w) )));
    a23 += integrate(_range=boundaryfaces(mesh),
                     _expr=-tau_constant * idt(phat) * pow(h(),M_tau_order)*id(w) );


    // <j.n,mu>_Omega/Gamma
    a31 += integrate(_range=internalfaces(mesh),
                     _expr=( id(l)*(leftfacet(trans(idt(u))*N())+
                                    rightfacet(trans(idt(u))*N())) ) );

    // <tau p, mu>_Omega/Gamma
    a32 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant * id(l) * ( leftfacet( pow(h(),M_tau_order)*idt(p) )+
                                                    rightfacet( pow(h(),M_tau_order)*idt(p) )));

    // <-tau phat, mu>_Omega/Gamma
    a33 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * idt(phat) * id(l) * ( leftface( pow(h(),M_tau_order) )+
                                                                 rightface( pow(h(),M_tau_order) )));


    // BC

    auto itField = M_modelProperties->boundaryConditions().find( "potential");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                // <phat, mu>_{Gamma_in,Gamma_out}
                a33 += integrate(_range=markedfaces(mesh,marker),
                                 _expr=idt(phat) * id(l) );
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                // <j.n,mu>_Gamma_E
                a31 += integrate(_range=markedfaces(mesh,marker),
                                 _expr=( id(l)*(trans(idt(u))*N()) ));

                // <tau p, mu>_Gamma_E
                a32 += integrate(_range=markedfaces(mesh,marker),
                                 _expr=tau_constant * id(l) * ( pow(h(),M_tau_order)*idt(p) ) );

                // <-tau phat, mu>_Gamma_E
                a33 += integrate(_range=markedfaces(mesh,marker),
                                 _expr=-tau_constant * idt(phat) * id(l) * ( pow(h(),M_tau_order) ) );
            }
        }
    }

    if ( M_integralCondition )
    {
        auto a14 = form2(_trial=M_Ch, _test=M_Vh,_matrix=M_A_cst,
                         _rowstart=0,
                         _colstart=3);
        auto a34 = form2(_trial=M_Ch, _test=M_Mh,_matrix=M_A_cst,
                         _rowstart=2,
                         _colstart=3);
        auto a41 = form2(_trial=M_Vh, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=0);
        auto a42 = form2(_trial=M_Wh, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=1);
        auto a43 = form2(_trial=M_Mh, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=2);
        auto a24 = form2(_trial=M_Ch, _test=M_Wh,_matrix=M_A_cst,
                         _rowstart=1,
                         _colstart=3);
        auto a44 = form2(_trial=M_Ch, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=3);

        itField = M_modelProperties->boundaryConditions().find( "current");
        // only if model contains integral
        if ( itField != M_modelProperties->boundaryConditions().end() )
        {
            auto mapField = (*itField).second;
            auto itType = mapField.find( "Integral" );
            if ( itType != mapField.end() )
            {
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();

                    // <lambda, v.n>_Gamma_out
                    a14 += integrate( _range=markedfaces(mesh,marker), _expr=trans(id(u))*N()*idt(nu) );

                    // <lambda, tau w>_Gamma_out
                    a24 += integrate( _range=markedfaces(mesh,marker),
                                      _expr=tau_constant * ( pow(h(),M_tau_order)*id(w) ) * idt(nu) );

                    // <lambda, -tau mu>_Gamma_out
                    a34 += integrate(_range=markedfaces(mesh,marker),
                                     _expr=-tau_constant * idt(nu) * id(l) * ( pow(h(),M_tau_order) ));

                    // <j.n, m>_Gamma_out
                    a41 += integrate( _range=markedfaces(mesh,marker), _expr=trans(idt(u))*N()*id(nu) );

                    // <tau p, m>_Gamma_out
                    a42 += integrate( _range=markedfaces(mesh,marker), _expr=tau_constant *
                                      ( pow(h(),M_tau_order)*idt(p) ) * id(nu) );

                    // <tau phat, m>_Gamma_out
                    a43 += integrate(_range=markedfaces(mesh,marker),
                                     _expr=-tau_constant * id(nu) * idt(phat) * ( pow(h(),M_tau_order) ));
                    a44 += integrate( _range=markedfaces(mesh,marker), _expr=-pow(h(),M_tau_order)*id(nu)*idt(nu) );
                }
            }
        }
    }
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::assembleF( int iter )
{
    M_F->zero();
    auto mesh = M_Vh->mesh();
    auto nu = M_Ch->element( "nu" );
    auto l = M_Mh->element( "lambda" );

    // Building the RHS

    auto rhs3 = form1( _test=M_Mh, _vector=M_F,
                       _rowstart=2);

    auto itField = M_modelProperties->boundaryConditions().find( "potential");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());
                // <V, mu>_Gamma_out
                rhs3 += integrate(_range=markedfaces(mesh,marker),
                                  _expr=id(l)*g);
            }
        }
    }

    if ( M_integralCondition )         // only if model contains integral
    {
        auto rhs4 = form1( _test=M_Ch, _vector=M_F,
                           _rowstart=3);
        itField = M_modelProperties->boundaryConditions().find( "current");
        if ( itField != M_modelProperties->boundaryConditions().end() )
        {
            auto mapField = (*itField).second;
            auto itType = mapField.find( "Integral" );
            if ( itType != mapField.end() )
            {
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();
                    auto g = expr(exAtMarker.expression());
                    // <I_target,m>_Gamma_out
                    rhs4 += integrate(_range=markedfaces(mesh,marker),
                                      _expr=g*id(nu));
                }
            }
        }
    }
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::solveT( int iter )
{
    auto mesh = M_Vh->mesh();

    auto Tw = M_modelProperties->parameters()["Tw"].value();

    auto T = M_Xh->element( "T" );
    auto q = M_Xh->element( "q" );

    auto ma = M_tBackend->newMatrix( _test=M_Xh, _trial=M_Xh );
    auto vf = M_tBackend->newVector( M_Xh );
    auto a = form2( _test=M_Xh, _trial=M_Xh, _matrix=ma );
    auto f = form1( _test=M_Xh, _vector=vf );

    for( auto const& pairMat : M_modelProperties->materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        if ( M_isPicard && iter > 0)
        {
            auto alpha = material.getDouble("alpha");
            auto T0 = material.getDouble("T0");
            auto sigma0 = material.getDouble("sigma0");
            auto sigma = material.getScalar("sigma", "T", idv(M_Tp));
            sigma.setParameterValues({{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
            auto k0 = material.getDouble("k0");
            auto k = material.getScalar("k", "T", idv(M_Tp));
            k.setParameterValues({{"k0",k0},{"T0",T0},{"alpha",alpha}});
            // (k grad(T), grad(q))_Omega
            a += integrate(_range=markedelements(mesh,marker),
                           _expr=k*gradt(T)*trans(grad(q)) );
            // (j^2/sigma,q)_Omega
            f += integrate(_range=markedelements(mesh, marker),
                           _expr=inner(idv(*M_up))/sigma * id(q));
        }
        else {
            // (k grad(T), grad(q))
            a += integrate(_range=markedelements(mesh,marker),
                           _expr=material.getScalar("k0")*gradt(T)*trans(grad(q)) );
            // (j^2/sigma,q)_Omega
            f += integrate(_range=markedelements(mesh, marker),
                           _expr=inner(idv(*M_up))/material.getScalar("sigma0") * id(q));
        }
    }

    auto itField = M_modelProperties->boundaryConditions().find( "temperature");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g1 = expr(exAtMarker.expression1());
                // <hT, q>_Gamma_C
                a += integrate(_range=markedfaces(mesh,marker),
                               _expr=g1*idt(T)*id(q));
                std::string gst = exAtMarker.expression2();
                auto g2 = expr(gst, {{"Tw",Tw}});
                // <hT_w,q>_Gamma_C
                f += integrate(_range=markedfaces(mesh,marker),
                               _expr=g2*id(q));
            }
        }
    }

    M_tBackend->solve( _matrix=ma, _rhs=vf, _solution=M_Tp);
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::exportResults()
{
    auto e = exporter( M_Vh->mesh());
    e->add("potential", *M_pp);
    e->add("current", *M_up);
    e->add("temperature", *M_Tp);
    e->save();
}

} // Feel
