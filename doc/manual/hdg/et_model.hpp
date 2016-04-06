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

    double M_V;
    int M_tau_order;

    void assembleA();
    void assembleF();

public:
    void init();
    void solve();
};

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::init()
{
    M_modelProperties = std::make_shared<model_prop_type>( Environment::expand( soption("model_json") ) );

    int proc_rank = Environment::worldComm().globalRank();

    M_V=doption("V");
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
    //auto face_mesh_bottom = createSubmesh( mesh, markedfaces(mesh,"bottom"), EXTRACTION_KEEP_MESH_RELATION, 0 );
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
    auto phatp = M_Mh->elementPtr( "phat" );
    M_Tp = M_Xh->elementPtr( "T" );
    auto mup = M_Ch->elementPtr( "c1" );

    tic();

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
    toc("matrices",true);
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::solve()
{
    assembleA();
    assembleF();
    backend(_rebuild=true)->solve( _matrix=M_A, _rhs=M_F, _solution=M_U );
    M_hdg_sol.localize(M_U);
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::assembleA()
{
    M_A->zero();
    auto mesh = M_Vh->mesh();

    // stabilisation parameter
    // auto k = expr(soption("k"));
    // auto sigma = expr(soption("sigma"));
    auto tau_constant = cst(doption("tau_constant"));

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

    auto a11 = form2( _trial=M_Vh, _test=M_Vh,_matrix=M_A );
    auto a12 = form2( _trial=M_Wh, _test=M_Vh,_matrix=M_A,
                      _rowstart=0,
                      _colstart=M_Vh->nLocalDofWithGhost() );
    auto a13 = form2( _trial=M_Mh, _test=M_Vh,_matrix=M_A,
                      _rowstart=0,
                      _colstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost());
    auto a14 = form2(_trial=M_Ch, _test=M_Vh,_matrix=M_A,
                     _rowstart=0,
                     _colstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost()+M_Mh->nLocalDofWithGhost());
    auto a21 = form2( _trial=M_Vh, _test=M_Wh,_matrix=M_A,
                      _rowstart=M_Vh->nLocalDofWithGhost(),
                      _colstart=0);
    auto a22 = form2( _trial=M_Wh, _test=M_Wh,_matrix=M_A,
                      _rowstart=M_Vh->nLocalDofWithGhost(),
                      _colstart=M_Vh->nLocalDofWithGhost() );
    auto a23 = form2( _trial=M_Mh, _test=M_Wh,_matrix=M_A,
                      _rowstart=M_Vh->nLocalDofWithGhost(),
                      _colstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost());
    auto a24 = form2(_trial=M_Ch, _test=M_Wh,_matrix=M_A,
                     _rowstart=M_Vh->nLocalDofWithGhost(),
                     _colstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost()+M_Mh->nLocalDofWithGhost());
    auto a31 = form2( _trial=M_Vh, _test=M_Mh,_matrix=M_A,
                      _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost(),
                      _colstart=0);
    auto a32 = form2( _trial=M_Wh, _test=M_Mh,_matrix=M_A,
                      _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost(),
                      _colstart=M_Vh->nLocalDofWithGhost());
    auto a33 = form2(_trial=M_Mh, _test=M_Mh,_matrix=M_A,
                     _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost(),
                     _colstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost());
    auto a34 = form2(_trial=M_Ch, _test=M_Mh,_matrix=M_A,
                     _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost(),
                     _colstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost()+M_Mh->nLocalDofWithGhost());
    auto a41 = form2(_trial=M_Vh, _test=M_Ch,_matrix=M_A,
                     _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost()+M_Mh->nLocalDofWithGhost(),
                     _colstart=0);
    auto a42 = form2(_trial=M_Wh, _test=M_Ch,_matrix=M_A,
                     _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost()+M_Mh->nLocalDofWithGhost(),
                     _colstart=M_Vh->nLocalDofWithGhost());
    auto a43 = form2(_trial=M_Mh, _test=M_Ch,_matrix=M_A,
                     _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost()+M_Mh->nLocalDofWithGhost(),
                     _colstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost());
    auto a55 = form2(_trial=M_Xh, _test=M_Xh,_matrix=M_A,
                     _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost()+M_Mh->nLocalDofWithGhost()+M_Ch->nLocalDofWithGhost(),
                     _colstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost()+M_Mh->nLocalDofWithGhost()+M_Ch->nLocalDofWithGhost());

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
        a11 += integrate(_range=markedelements(mesh,marker),_expr=(trans(idt(u))*id(v))/material.getScalar("sigma") );
        a55 += integrate(_range=markedelements(mesh,marker), _expr=material.getScalar("k")*gradt(T)*trans(grad(T)));
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
                std::string gst = exAtMarker.expression1();
                auto mats = M_modelProperties->materials();
                auto mat = mats.material( marker );
                auto h = boost::lexical_cast<double>(mat.getString( "h" ));
                auto g = expr(gst, {{"h",h}});
                a55 += integrate(_range=markedfaces(mesh,marker), _expr=g*idt(T)*id(q));
            }
        }
    }
}

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::assembleF()
{
    M_F->zero();
    auto mesh = M_Vh->mesh();
    auto nu = M_Ch->element( "nu" );
    auto l = M_Mh->element( "lambda" );
    auto q = M_Xh->element( "q" );

    // Building the RHS

    auto rhs3 = form1( _test=M_Mh, _vector=M_F,
                       _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost());
    auto rhs4 = form1( _test=M_Xh, _vector=M_F,
                       _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost()+M_Mh->nLocalDofWithGhost()+M_Ch->nLocalDofWithGhost());
    auto rhs5 = form1( _test=M_Ch, _vector=M_F,
                       _rowstart=M_Vh->nLocalDofWithGhost()+M_Wh->nLocalDofWithGhost()+M_Mh->nLocalDofWithGhost());

    for( auto const& pairMat : M_modelProperties->materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        rhs4 += integrate(_range=markedelements(mesh, marker),_expr=inner(idv(*M_up))/material.getScalar("sigma") * id(q));
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
                auto g = expr(exAtMarker.expression());
                rhs3 += integrate(_range=markedfaces(mesh,marker),
                                  _expr=id(l)*g);
            }
        }
    }

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

    auto Tw = doption("Tw");
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
                auto mats = M_modelProperties->materials();
                auto mat = mats.material( marker );
                auto h = boost::lexical_cast<double>(mat.getString( "h" ));
                auto g = expr(gst, {{"h",h}, {"Tw",Tw}});
                rhs4 += integrate(_range=markedfaces(mesh,marker),
                                  _expr=g*id(q));
            }
        }
    }
}


} // Feel
