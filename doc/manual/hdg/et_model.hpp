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

#include <boost/optional.hpp>
#include <boost/algorithm/string.hpp>

namespace Feel {



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

    sparse_matrix_ptrtype M_A;
    vector_ptrtype M_F;
    BlocksBaseVector<double> M_hdg_sol;
    vector_ptrtype M_U;

    Vh_element_ptr_t M_up;
    Wh_element_ptr_t M_pp;
    Xh_element_ptr_t M_Tp;

    boost::optional<double> M_V;

    int M_tau_order;

    bool M_integralCondition;
    bool M_isPicard;

    void initGraphs();
    void initGraphsWithIntegralCond();
    void assembleA( int iter = -1 );
    void assembleF( int iter = -1 );

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

    M_V = boost::none;

    M_tau_order = ioption("tau_order");

    tic();
    auto mesh = loadMesh( new mesh_type);
    auto face_mesh = createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );
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
         << "Xh<" << OrderP << "> : " << M_Xh->nDof() << std::endl
         << "Ch<" << 0 << "> : " << M_Ch->nDof() << std::endl;

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
            backend(_rebuild=true)->solve( _matrix=M_A, _rhs=M_F, _solution=M_U );
            M_hdg_sol.localize(M_U);
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
        backend(_rebuild=true)->solve( _matrix=M_A, _rhs=M_F, _solution=M_U );
        M_hdg_sol.localize(M_U);
    }
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::initGraphs()
{
    auto phatp = M_Mh->elementPtr( "phat" );

    BlocksBaseGraphCSR hdg_graph(4,4);
    hdg_graph(0,0) = stencil( _test=M_Vh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,0) = stencil( _test=M_Wh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,0) = stencil( _test=M_Mh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,0) = stencil( _test=M_Xh,_trial=M_Vh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,1) = stencil( _test=M_Vh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,1) = stencil( _test=M_Wh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,1) = stencil( _test=M_Mh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,1) = stencil( _test=M_Xh,_trial=M_Wh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,2) = stencil( _test=M_Vh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,2) = stencil( _test=M_Wh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,2) = stencil( _test=M_Mh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,2) = stencil( _test=M_Xh,_trial=M_Mh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,3) = stencil( _test=M_Vh,_trial=M_Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(1,3) = stencil( _test=M_Wh,_trial=M_Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(2,3) = stencil( _test=M_Mh,_trial=M_Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(3,3) = stencil( _test=M_Xh,_trial=M_Xh, _diag_is_nonzero=false, _close=false)->graph();

    M_A = backend()->newBlockMatrix(_block=hdg_graph);

    BlocksBaseVector<double> hdg_vec(4);
    hdg_vec(0,0) = backend()->newVector( M_Vh );
    hdg_vec(1,0) = backend()->newVector( M_Wh );
    hdg_vec(2,0) = backend()->newVector( M_Mh );
    hdg_vec(3,0) = backend()->newVector( M_Xh );
    M_F = backend()->newBlockVector(_block=hdg_vec, _copy_values=false);

    M_hdg_sol = BlocksBaseVector<double>(4);
    M_hdg_sol(0,0) = M_up;
    M_hdg_sol(1,0) = M_pp;
    M_hdg_sol(2,0) = phatp;
    M_hdg_sol(3,0) = M_Tp;
    M_U = backend()->newBlockVector(_block=M_hdg_sol, _copy_values=false);
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::initGraphsWithIntegralCond()
{
    auto phatp = M_Mh->elementPtr( "phat" );
    auto mup = M_Ch->elementPtr( "c1" );

    BlocksBaseGraphCSR hdg_graph(5,5);
    hdg_graph(0,0) = stencil( _test=M_Vh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,0) = stencil( _test=M_Wh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,0) = stencil( _test=M_Mh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,0) = stencil( _test=M_Ch,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(4,0) = stencil( _test=M_Xh,_trial=M_Vh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,1) = stencil( _test=M_Vh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,1) = stencil( _test=M_Wh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,1) = stencil( _test=M_Mh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,1) = stencil( _test=M_Ch,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(4,1) = stencil( _test=M_Xh,_trial=M_Wh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,2) = stencil( _test=M_Vh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,2) = stencil( _test=M_Wh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,2) = stencil( _test=M_Mh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,2) = stencil( _test=M_Ch,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(4,2) = stencil( _test=M_Xh,_trial=M_Mh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,3) = stencil( _test=M_Vh,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,3) = stencil( _test=M_Wh,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,3) = stencil( _test=M_Mh,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,3) = stencil( _test=M_Ch,_trial=M_Ch, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(4,3) = stencil( _test=M_Xh,_trial=M_Ch, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,4) = stencil( _test=M_Vh,_trial=M_Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(1,4) = stencil( _test=M_Wh,_trial=M_Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(2,4) = stencil( _test=M_Mh,_trial=M_Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(3,4) = stencil( _test=M_Ch,_trial=M_Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(4,4) = stencil( _test=M_Xh,_trial=M_Xh, _diag_is_nonzero=false, _close=false)->graph();

    M_A = backend()->newBlockMatrix(_block=hdg_graph);

    BlocksBaseVector<double> hdg_vec(5);
    hdg_vec(0,0) = backend()->newVector( M_Vh );
    hdg_vec(1,0) = backend()->newVector( M_Wh );
    hdg_vec(2,0) = backend()->newVector( M_Mh );
    hdg_vec(3,0) = backend()->newVector( M_Ch );
    hdg_vec(4,0) = backend()->newVector( M_Xh );
    M_F = backend()->newBlockVector(_block=hdg_vec, _copy_values=false);

    M_hdg_sol = BlocksBaseVector<double>(5);
    M_hdg_sol(0,0) = M_up;
    M_hdg_sol(1,0) = M_pp;
    M_hdg_sol(2,0) = phatp;
    M_hdg_sol(3,0) = mup;
    M_hdg_sol(4,0) = M_Tp;
    M_U = backend()->newBlockVector(_block=M_hdg_sol, _copy_values=false);
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::assembleA( int iter )
{
    M_A->zero();
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
    size_type VhDof = M_Vh->nLocalDofWithGhost();
    size_type WhDof = M_Wh->nLocalDofWithGhost();
    size_type MhDof = M_Mh->nLocalDofWithGhost();
    size_type ChDof = M_integralCondition ? M_Ch->nLocalDofWithGhost() : 0;

    auto a11 = form2( _trial=M_Vh, _test=M_Vh,_matrix=M_A );
    auto a12 = form2( _trial=M_Wh, _test=M_Vh,_matrix=M_A,
                      _rowstart=0,
                      _colstart=VhDof );
    auto a13 = form2( _trial=M_Mh, _test=M_Vh,_matrix=M_A,
                      _rowstart=0,
                      _colstart=VhDof+WhDof);
    auto a21 = form2( _trial=M_Vh, _test=M_Wh,_matrix=M_A,
                      _rowstart=VhDof,
                      _colstart=0);
    auto a22 = form2( _trial=M_Wh, _test=M_Wh,_matrix=M_A,
                      _rowstart=VhDof,
                      _colstart=VhDof );
    auto a23 = form2( _trial=M_Mh, _test=M_Wh,_matrix=M_A,
                      _rowstart=VhDof,
                      _colstart=VhDof+WhDof);
    auto a31 = form2( _trial=M_Vh, _test=M_Mh,_matrix=M_A,
                      _rowstart=VhDof+WhDof,
                      _colstart=0);
    auto a32 = form2( _trial=M_Wh, _test=M_Mh,_matrix=M_A,
                      _rowstart=VhDof+WhDof,
                      _colstart=VhDof);
    auto a33 = form2(_trial=M_Mh, _test=M_Mh,_matrix=M_A,
                     _rowstart=VhDof+WhDof,
                     _colstart=VhDof+WhDof);
    auto a55 = form2(_trial=M_Xh, _test=M_Xh,_matrix=M_A,
                     _rowstart=VhDof+WhDof+MhDof+ChDof,
                     _colstart=VhDof+WhDof+MhDof+ChDof);

    a12 += integrate(_range=elements(mesh),_expr=-(idt(p)*div(v)));

    a13 += integrate(_range=internalfaces(mesh),
                     _expr=( idt(phat)*leftface(trans(id(v))*N())+
                             idt(phat)*rightface(trans(id(v))*N())) );
    a13 += integrate(_range=boundaryfaces(mesh),
                     _expr=idt(phat)*trans(id(v))*N());


    a21 += integrate(_range=elements(mesh),_expr=(-grad(w)*idt(u)));
    a21 += integrate(_range=internalfaces(mesh),
                     _expr=( leftface(id(w))*leftfacet(trans(idt(u))*N()) ) );
    a21 += integrate(_range=internalfaces(mesh),
                     _expr=(rightface(id(w))*rightfacet(trans(idt(u))*N())) );
    a21 += integrate(_range=boundaryfaces(mesh),
                     _expr=(id(w)*trans(idt(u))*N()));

    a22 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant *
                     ( leftfacet( pow(h(),M_tau_order)*idt(p))*leftface(id(w)) +
                       rightfacet( pow(h(),M_tau_order)*idt(p))*rightface(id(w) )));
    a22 += integrate(_range=boundaryfaces(mesh),
                     _expr=(tau_constant * pow(h(),M_tau_order)*id(w)*idt(p)));

    a23 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * idt(phat) *
                     ( leftface( pow(h(),M_tau_order)*id(w) )+
                       rightface( pow(h(),M_tau_order)*id(w) )));
    a23 += integrate(_range=boundaryfaces(mesh),
                     _expr=-tau_constant * idt(phat) * pow(h(),M_tau_order)*id(w) );


    a31 += integrate(_range=internalfaces(mesh),
                     _expr=( id(l)*(leftfacet(trans(idt(u))*N())+
                                    rightfacet(trans(idt(u))*N())) ) );

    a32 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant * id(l) * ( leftfacet( pow(h(),M_tau_order)*idt(p) )+
                                                    rightfacet( pow(h(),M_tau_order)*idt(p) )));

    a33 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * idt(phat) * id(l) * ( leftface( pow(h(),M_tau_order) )+
                                                                 rightface( pow(h(),M_tau_order) )));

    for( auto const& pairMat : M_modelProperties->materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;

        if ( M_isPicard && iter > 0)
        {
            auto sigma0 = material.getDouble("sigma0");
            auto alpha = material.getDouble("alpha");
            auto sigma = material.getScalar("sigma", "T", idv(M_Tp));
            sigma.setParameterValues({{"sigma0",sigma0},{"alpha",alpha},{"Tw",Tw}});
            auto k0 = material.getDouble("k0");
            auto k = material.getScalar("k", "T", idv(M_Tp));
            k.setParameterValues({{"k0",k0},{"Tw",Tw},{"alpha",alpha}});

            a11 += integrate(_range=markedelements(mesh,marker), _expr=(trans(idt(u))*id(v))/sigma );
            a55 += integrate(_range=markedelements(mesh,marker), _expr=k*gradt(T)*trans(grad(T)) );
        }
        else {
            a11 += integrate(_range=markedelements(mesh,marker), _expr=(trans(idt(u))*id(v))/material.getScalar("sigma0") );
            a55 += integrate(_range=markedelements(mesh,marker), _expr=material.getScalar("k0")*gradt(T)*trans(grad(T)) );
        }
    }


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
                a31 += integrate(_range=markedfaces(mesh,marker),
                                 _expr=( id(l)*(trans(idt(u))*N()) ));

                a32 += integrate(_range=markedfaces(mesh,marker),
                                 _expr=tau_constant * id(l) * ( pow(h(),M_tau_order)*idt(p) ) );

                a33 += integrate(_range=markedfaces(mesh,marker),
                                 _expr=-tau_constant * idt(phat) * id(l) * ( pow(h(),M_tau_order) ) );
            }
        }
    }

    if ( M_integralCondition )
    {
        auto a14 = form2(_trial=M_Ch, _test=M_Vh,_matrix=M_A,
                         _rowstart=0,
                         _colstart=VhDof+WhDof+MhDof);
        auto a34 = form2(_trial=M_Ch, _test=M_Mh,_matrix=M_A,
                         _rowstart=VhDof+WhDof,
                         _colstart=VhDof+WhDof+MhDof);
        auto a41 = form2(_trial=M_Vh, _test=M_Ch,_matrix=M_A,
                         _rowstart=VhDof+WhDof+MhDof,
                         _colstart=0);
        auto a42 = form2(_trial=M_Wh, _test=M_Ch,_matrix=M_A,
                         _rowstart=VhDof+WhDof+MhDof,
                         _colstart=VhDof);
        auto a43 = form2(_trial=M_Mh, _test=M_Ch,_matrix=M_A,
                         _rowstart=VhDof+WhDof+MhDof,
                         _colstart=VhDof+WhDof);
        auto a24 = form2(_trial=M_Ch, _test=M_Wh,_matrix=M_A,
                         _rowstart=VhDof,
                         _colstart=VhDof+WhDof+MhDof);
        itField = M_modelProperties->boundaryConditions().find( "current");
        // only if model contains integral
        if ( itField != M_modelProperties->boundaryConditions().end() )
        {
            auto mapField = (*itField).second;
            auto itType = mapField.find( "Dirichlet" );
            if ( itType != mapField.end() )
            {
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();

                    a14 += integrate( _range=markedfaces(mesh,marker), _expr=trans(id(u))*N()*idt(nu) );

                    a24 += integrate( _range=markedfaces(mesh,marker),
                                      _expr=tau_constant * ( pow(h(),M_tau_order)*id(w) ) * idt(nu) );

                    a34 += integrate(_range=markedfaces(mesh,marker),
                                     _expr=-tau_constant * idt(nu) * id(l) * ( pow(h(),M_tau_order) ));

                    a41 += integrate( _range=markedfaces(mesh,marker), _expr=trans(idt(u))*N()*id(nu) );

                    a42 += integrate( _range=markedfaces(mesh,marker), _expr=tau_constant *
                                      ( pow(h(),M_tau_order)*idt(p) ) * id(nu) );

                    a43 += integrate(_range=markedfaces(mesh,marker),
                                     _expr=-tau_constant * id(nu) * idt(phat) * ( pow(h(),M_tau_order) ));

                }
            }
        }
    }

    itField = M_modelProperties->boundaryConditions().find( "temperature");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression1());
                a55 += integrate(_range=markedfaces(mesh,marker), _expr=g*idt(T)*id(q));
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
    auto q = M_Xh->element( "q" );

    auto Tw = M_modelProperties->parameters()["Tw"].value();

    // Building the RHS

    size_type VhDof = M_Vh->nLocalDofWithGhost();
    size_type WhDof = M_Wh->nLocalDofWithGhost();
    size_type MhDof = M_Mh->nLocalDofWithGhost();
    size_type ChDof = M_integralCondition ? M_Ch->nLocalDofWithGhost() : 0;

    auto rhs3 = form1( _test=M_Mh, _vector=M_F,
                       _rowstart=VhDof+WhDof);
    auto rhs4 = form1( _test=M_Xh, _vector=M_F,
                       _rowstart=VhDof+WhDof+MhDof+ChDof);

    for( auto const& pairMat : M_modelProperties->materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        if ( M_isPicard && iter > 0 )
        {
            auto sigma0 = material.getDouble("sigma0");
            auto alpha = material.getDouble("alpha");
            auto sigma = material.getScalar("sigma", "T", idv(M_Tp));
            sigma.setParameterValues({{"sigma0",sigma0},{"alpha",alpha},{"Tw",Tw}});

            rhs4 += integrate(_range=markedelements(mesh, marker),_expr=inner(idv(*M_up))/sigma * id(q));
        }
        else
        {
            rhs4 += integrate(_range=markedelements(mesh, marker),_expr=inner(idv(*M_up))/material.getScalar("sigma0") * id(q));
        }
    }

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
                if ( M_V )
                {
                    double V = *M_V;
                    rhs3 += integrate(_range=markedfaces(mesh,marker),
                                      _expr=id(l)*V);
                }
                else {
                    auto g = expr(exAtMarker.expression());
                    rhs3 += integrate(_range=markedfaces(mesh,marker),
                                      _expr=id(l)*g);
                }
            }
        }
    }

    if ( M_integralCondition )         // only if model contains integral
    {
        auto rhs5 = form1( _test=M_Ch, _vector=M_F,
                           _rowstart=VhDof+WhDof+MhDof);
        itField = M_modelProperties->boundaryConditions().find( "current");
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
                    rhs5 += integrate(_range=markedfaces(mesh,marker),
                                      _expr=g*id(nu));
                }
            }
        }
    }

    itField = M_modelProperties->boundaryConditions().find( "temperature");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                std::string gst = exAtMarker.expression2();
                auto g = expr(gst, {{"Tw",Tw}});
                rhs4 += integrate(_range=markedfaces(mesh,marker),
                                  _expr=g*id(q));
            }
        }
    }
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
