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
#include <feel/feel.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>

namespace Feel {
namespace FeelModels{


template<int Dim, int OrderP>
class ElectroThermal
    :
public Application
{
    typedef Application super;

public:

    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename std::shared_ptr<backend_type> backend_ptrtype ;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,1,Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef std::shared_ptr<face_mesh_type> face_mesh_ptrtype;

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


    /**
     * Constructor
     */
    ElectroThermal()
        :
        super()

    {
        this->changeRepository( boost::format( "ElectroThermal/%1%/h_%2%/" )
                                % this->about().appName()
                                % doption("gmsh.hsize")
                                );
    }

    template<typename MatrixType, typename VectorType, typename Vht, typename Wht, typename Mht, typename Cht, typename Xht,
             typename Vhet, typename Whet, typename Xhet>
    void
    assemble_A_and_F( MatrixType A,
                      VectorType F,
                      Vht Vh,
                      Wht Wh,
                      Mht Mh,
                      Cht Ch,
                      Xht Xh,
                      Vhet,
                      Whet,
                      Xhet );

    void run();
    double V;
};

template<int Dim, int OrderP>
void
ElectroThermal<Dim, OrderP>::run()
{
    int proc_rank = Environment::worldComm().globalRank();

    bool check_errors = boption("check_errors");

    V=doption("V");
    tic();
    auto mesh = loadMesh( new mesh_type );
    toc("mesh",true);

    // ****** Hybrid-mixed formulation ******
    // We treat Vh, Wh, and Mh separately
    tic();

    auto Vh = Pdhv<OrderP>( mesh, true );
    auto Wh = Pdh<OrderP>( mesh, true );
    auto face_mesh = createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );
    //auto face_mesh_bottom = createSubmesh( mesh, markedfaces(mesh,"bottom"), EXTRACTION_KEEP_MESH_RELATION, 0 );
    auto Mh = Pdh<OrderP>( face_mesh, true );
    auto Xh = Pch<OrderP>( mesh );
    auto Ch = Pch<0>( mesh );

    toc("spaces",true);

    cout << "Vh<" << OrderP << "> : " << Vh->nDof() << std::endl
         << "Wh<" << OrderP << "> : " << Wh->nDof() << std::endl
         << "Mh<" << OrderP << "> : " << Mh->nDof() << std::endl
         << "Xh<" << OrderP << "> : " << Xh->nDof() << std::endl
         << "Ch<" << 0 << "> : " << Ch->nDof() << std::endl;

    auto up = Vh->elementPtr( "u" );
    auto pp = Wh->elementPtr( "p" );
    auto phatp = Mh->elementPtr( "phat" );
    auto Tp = Xh->elementPtr( "T" );
    auto mup = Ch->elementPtr( "c1" );

    auto u = Vh->element( "u" );
    auto v = Vh->element( "v" );
    auto p = Wh->element( "p" );
    auto q = Xh->element( "q" );
    auto T = Xh->element( "T" );
    auto w = Wh->element( "w" );
    auto phat = Mh->element( "phat" );
    auto l = Mh->element( "lambda" );
    auto nu = Ch->element( "nu" );

    // Number of dofs associated with each space U and P
    auto nDofu = u.functionSpace()->nDof();
    auto nDofp = p.functionSpace()->nDof();
    auto nDofphat = phat.functionSpace()->nDof();
    auto nDofT = T.functionSpace()->nDof();

    tic();
#if 1
    BlocksBaseGraphCSR hdg_graph(5,5);
    hdg_graph(0,0) = stencil( _test=Vh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,0) = stencil( _test=Wh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,0) = stencil( _test=Mh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,0) = stencil( _test=Ch,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(4,0) = stencil( _test=Xh,_trial=Vh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,1) = stencil( _test=Vh,_trial=Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,1) = stencil( _test=Wh,_trial=Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,1) = stencil( _test=Mh,_trial=Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,1) = stencil( _test=Ch,_trial=Wh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(4,1) = stencil( _test=Xh,_trial=Wh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,2) = stencil( _test=Vh,_trial=Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,2) = stencil( _test=Wh,_trial=Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,2) = stencil( _test=Mh,_trial=Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,2) = stencil( _test=Ch,_trial=Mh, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(4,2) = stencil( _test=Xh,_trial=Mh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,3) = stencil( _test=Vh,_trial=Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(1,3) = stencil( _test=Wh,_trial=Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(2,3) = stencil( _test=Mh,_trial=Ch, _diag_is_nonzero=false, _close=false)->graph();
    hdg_graph(3,3) = stencil( _test=Ch,_trial=Ch, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(4,3) = stencil( _test=Xh,_trial=Ch, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    hdg_graph(0,4) = stencil( _test=Vh,_trial=Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(1,4) = stencil( _test=Wh,_trial=Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(2,4) = stencil( _test=Mh,_trial=Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(3,4) = stencil( _test=Ch,_trial=Xh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    hdg_graph(4,4) = stencil( _test=Xh,_trial=Xh, _diag_is_nonzero=false, _close=false)->graph();

    auto A = backend()->newBlockMatrix(_block=hdg_graph);

    #else
    // build the big matrix associated to bilinear form over Vh x Wh x Mh
    auto A = backend()->newBlockMatrix(_block=csrGraphBlocks(Xh,Wh,Mh,Ch,Xh));
#endif
    BlocksBaseVector<double> hdg_vec(5);
    hdg_vec(0,0) = backend()->newVector( Vh );
    hdg_vec(1,0) = backend()->newVector( Wh );
    hdg_vec(2,0) = backend()->newVector( Mh );
    hdg_vec(3,0) = backend()->newVector( Ch );
    hdg_vec(4,0) = backend()->newVector( Xh );
    auto F = backend()->newBlockVector(_block=hdg_vec, _copy_values=false);

    BlocksBaseVector<double> hdg_sol(5);
    hdg_sol(0,0) = up;
    hdg_sol(1,0) = pp;
    hdg_sol(2,0) = phatp;
    hdg_sol(3,0) = mup;
    hdg_sol(4,0) = Tp;
    auto U = backend()->newBlockVector(_block=hdg_sol, _copy_values=false);
    toc("matrices",true);

    int itmax = ioption( "picard.itmax" );
    double tol = doption( "picard.itol" );
    auto bdfT = bdf(_space=Xh);

    auto po = Wh->element();
    auto To = Xh->element();
    cout << "#iteration time current" << std::endl;
    auto e = exporter(_mesh=mesh);
    double Kp  = doption( _name="Kp" );
    double Ki  = doption( _name="Ki" );
    double Kd  = doption( _name="Kd" );
    double I_target=doption("I_target");
    double I_measured = 0;
    double I_error  = I_target-I_measured;
    double I_error_last = I_error ;
    double I_In = 0.5*I_error ;

#if 1
    assemble_A_and_F( A, F, Vh, Wh, Mh, Ch, Xh, up, pp, Tp );
    backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );
    hdg_sol.localize(U);
    double I =  integrate(_range=markedfaces(mesh,"bottom"),
                          _expr=inner(idv(*up),N())).evaluate()(0,0);
    cout << "I=" << I << std::endl;

    // ****** Compute error ******

    auto sigma = expr(soption("sigma"));
    auto p_exact = expr(soption("p_exact"));
    auto gradp_exact = grad<Dim>(p_exact);
    auto u_exact = -sigma*trans(gradp_exact);

    if (check_errors)
    {
        bool has_dirichlet = nelements(markedfaces(mesh,{"top"}), true) >= 1;
        auto l2err_u = normL2( _range=elements(mesh), _expr=u_exact - idv(*up) );

        double l2err_p = 1e+30;
        if ( has_dirichlet )
        {
            l2err_p = normL2( _range=elements(mesh), _expr=p_exact - idv(*pp) );
        }
        else
        {
            auto mean_p_exact = mean( elements(mesh), p_exact )(0,0);
            auto mean_p = mean( elements(mesh), idv(*pp) )(0,0);
            l2err_p = normL2( elements(mesh),
                              (p_exact - cst(mean_p_exact)) - (idv(*pp) - cst(mean_p)) );
        }
        cout << "||u_exact - u||_L2 = " << l2err_u << std::endl;
        cout << "||p_exact - p||_L2 = " << l2err_p << std::endl;
    }

    e->add( "V", *pp );
    e->add( "J", *up );
    e->add( "T", *Tp );
    e->save();
#else
    // time loop for control
    for ( bdfT->start(T) ;
          (!bdfT->isFinished() );
          bdfT->next(T) )
    {
        // picard loop
        cout << "  #iteration incrp incrT current" << std::endl;
        int it = 0;
        double incrp, incrT;
        do
        {
            assemble_A_and_F( A, F, Vh, Wh, Mh, Ch, Xh, up, pp, Tp );
            tic();
            backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );
            hdg_sol.localize(U);
            toc("solve",false);

            incrp = normL2( _range=elements(mesh), _expr=idv(*pp)-idv(po) );
            incrT = normL2( _range=elements(mesh), _expr=idv(*Tp)-idv(To) );

            po = *pp;
            To=*Tp;
            // compute current
            double I = integrate(_range=markedfaces(mesh,"bottom"),
                                 _expr=inner(idv(*up),N())).evaluate()(0,0);
            cout << "  picard #" << it << " " << incrp << " " << incrT << " " << I << std::endl;
        } while ( ( incrp > tol || incrT > tol ) && ( ++it < itmax ) );
        T = *Tp;
        I_measured = integrate(_range=markedfaces(mesh,"bottom"),
                               _expr=inner(idv(*up),N())).evaluate()(0,0);


        //update error
        I_error_last = I_error ;
        I_error = I_target-I_measured;
        I_In += 0.5*I_error ;
        cout << bdfT->iteration() << " " << bdfT->time() << "s  "  << V << " "
             << I_measured << " " << I_error << " " << I_error_last << " "
             << I_target << std::endl;
        V = Kp*I_error + Ki*I_In + Kd*(I_error- I_error_last)/0.5 ;


        e->step(bdfT->time())->add( "V", *pp );
        e->step(bdfT->time())->add( "J", *up );
        e->step(bdfT->time())->add( "T", *Tp );
        e->save();
    }
#endif
}

template<int Dim, int OrderP>
template<typename MatrixType, typename VectorType, typename Vht, typename Wht, typename Mht, typename Cht, typename Xht,
         typename Vhet, typename Whet, typename Xhet>
void
ElectroThermal<Dim, OrderP>::assemble_A_and_F( MatrixType A,
                                               VectorType F,
                                               Vht Vh,
                                               Wht Wh,
                                               Mht Mh,
                                               Cht Ch,
                                               Xht Xh,
                                               Vhet Jp,
                                               Whet Vp,
                                               Xhet Tp
                                               )
{
    A->zero();
    F->zero();
    // stabilisation parameter
    int M_tau_order = ioption("tau_order");

    auto k = expr(soption("k"));
    auto sigma = expr(soption("sigma"));
    auto tau_constant = cst(doption("tau_constant"));
    auto mesh = Vh->mesh();

    auto T = Xh->element( "T" );
    auto q = Xh->element( "q" );

    auto u = Vh->element( "u" );
    auto v = Vh->element( "v" );
    auto p = Wh->element( "p" );
    auto w = Wh->element( "w" );
    auto nu = Ch->element( "nu" );
    auto phat = Mh->element( "phat" );
    auto l = Mh->element( "lambda" );

    bool check_errors = boption("check_errors");

    auto p_exact = expr(soption("p_exact"));
    auto gradp_exact = grad<Dim>(p_exact);
    auto u_exact = -sigma*trans(gradp_exact);
    auto f = -sigma*laplacian(p_exact);

    // Building the RHS

    auto rhs2 = form1( _test=Wh, _vector=F,
                       _rowstart=Vh->nLocalDofWithGhost() );

    rhs2 += integrate(_range=elements(mesh), _expr=f*id(w));

    auto rhs3 = form1( _test=Mh, _vector=F,
                       _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());
    if (check_errors)
    {
        rhs3 += integrate(_range=markedfaces(mesh,{"top"}),
                          _expr=id(l)*p_exact);
        rhs3 += integrate(_range=markedfaces(mesh,"R"),
                          _expr=-id(l)*(sigma*gradp_exact*N()) );
    }
    else
        rhs3 += integrate(_range=markedfaces(mesh,"top"),
                          _expr=id(l)*V);

    auto rhs4 = form1( _test=Xh, _vector=F,
                       _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost()+Mh->nLocalDofWithGhost()+Ch->nLocalDofWithGhost());
    rhs4 += integrate(_range=elements(mesh),_expr=inner(idv(*Jp))/sigma * id(q));
    rhs4 += integrate(_range=markedfaces(mesh,"R"),
                      _expr=doption("h")*doption("Tw")*id(q));

    auto rhs5 = form1( _test=Ch, _vector=F,
                       _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost()+Mh->nLocalDofWithGhost());
    if (check_errors)
        rhs5 += integrate(_range=markedfaces(mesh,"bottom"),
                          _expr=(-sigma*gradp_exact*N())*id(nu));
    else
        rhs5 += integrate(_range=markedfaces(mesh,"bottom"),
                          _expr=doption("I_target")*id(nu));

    auto a11 = form2( _trial=Vh, _test=Vh,_matrix=A );
    a11 += integrate(_range=elements(mesh),_expr=(trans(idt(u))*id(v))/sigma );


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

    auto a14 = form2(_trial=Ch, _test=Vh,_matrix=A,
                     _rowstart=0,
                     _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost()+Mh->nLocalDofWithGhost());
    a14 += integrate( _range=markedfaces(mesh,"bottom"), _expr=trans(id(u))*N()*idt(nu) );

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
                     ( leftfacet( pow(h(),M_tau_order)*idt(p))*leftface(id(w)) +
                       rightfacet( pow(h(),M_tau_order)*idt(p))*rightface(id(w) )));
    a22 += integrate(_range=boundaryfaces(mesh),
                     _expr=(tau_constant * pow(h(),M_tau_order)*id(w)*idt(p)));

    auto a23 = form2( _trial=Mh, _test=Wh,_matrix=A,
                      _rowstart=Vh->nLocalDofWithGhost(), _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());
    a23 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * idt(phat) *
                     ( leftface( pow(h(),M_tau_order)*id(w) )+
                       rightface( pow(h(),M_tau_order)*id(w) )));

    a23 += integrate(_range=boundaryfaces(mesh),
                     _expr=-tau_constant * idt(phat) * pow(h(),M_tau_order)*id(w) );

    auto a24 = form2(_trial=Ch, _test=Wh,_matrix=A,
                     _rowstart=Vh->nLocalDofWithGhost(),
                     _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost()+Mh->nLocalDofWithGhost());
    a24 += integrate( _range=markedfaces(mesh,"bottom"), _expr=tau_constant *
                      ( pow(h(),M_tau_order)*id(w) ) * idt(nu) );

    auto a31 = form2( _trial=Vh, _test=Mh,_matrix=A,
                      _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(), _colstart=0);
    a31 += integrate(_range=internalfaces(mesh),
                     _expr=( id(l)*(leftfacet(trans(idt(u))*N())+
                                    rightfacet(trans(idt(u))*N())) ) );

    // BC
    a31 += integrate(_range=markedfaces(mesh,"R"),
                     _expr=( id(l)*(trans(idt(u))*N()) ));

    auto a32 = form2( _trial=Wh, _test=Mh,_matrix=A,
                      _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(), _colstart=Vh->nLocalDofWithGhost());
    a32 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant * id(l) * ( leftfacet( pow(h(),M_tau_order)*idt(p) )+
                                                    rightfacet( pow(h(),M_tau_order)*idt(p) )));

    a32 += integrate(_range=markedfaces(mesh,"R"),
                     _expr=tau_constant * id(l) * ( pow(h(),M_tau_order)*idt(p) ) );


    auto a33 = form2(_trial=Mh, _test=Mh,_matrix=A,
                     _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(), _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());

    a33 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * idt(phat) * id(l) * ( leftface( pow(h(),M_tau_order) )+
                                                                 rightface( pow(h(),M_tau_order) )));

    a33 += integrate(_range=markedfaces(mesh,"R"),
                     _expr=-tau_constant * idt(phat) * id(l) * ( pow(h(),M_tau_order) ) );
    a33 += integrate(_range=markedfaces(mesh,{"top"}),
                     _expr=idt(phat) * id(l) );

    auto a34 = form2(_trial=Ch, _test=Mh,_matrix=A,
                     _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost(), _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost()+Mh->nLocalDofWithGhost());

    a34 += integrate(_range=markedfaces(mesh,"bottom"),
                     _expr=-tau_constant * idt(nu) * id(l) * ( pow(h(),M_tau_order) ));

    auto a41 = form2(_trial=Vh, _test=Ch,_matrix=A,
                     _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost()+Mh->nLocalDofWithGhost(),
                     _colstart=0);
    a41 += integrate( _range=markedfaces(mesh,"bottom"), _expr=trans(idt(u))*N()*id(nu) );

    auto a42 = form2(_trial=Wh, _test=Ch,_matrix=A,
                     _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost()+Mh->nLocalDofWithGhost(),
                     _colstart=Vh->nLocalDofWithGhost());
    a42 += integrate( _range=markedfaces(mesh,"bottom"), _expr=tau_constant *
                      ( pow(h(),M_tau_order)*idt(p) ) * id(nu) );

    auto a43 = form2(_trial=Mh, _test=Ch,_matrix=A,
                     _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost()+Mh->nLocalDofWithGhost(),
                     _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost());

    a43 += integrate(_range=markedfaces(mesh,"bottom"),
                     _expr=-tau_constant * id(nu) * idt(phat) * ( pow(h(),M_tau_order) ));

    auto a55 = form2(_trial=Xh, _test=Xh,_matrix=A,
                     _rowstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost()+Mh->nLocalDofWithGhost()+Ch->nLocalDofWithGhost(),
                     _colstart=Vh->nLocalDofWithGhost()+Wh->nLocalDofWithGhost()+Mh->nLocalDofWithGhost()+Ch->nLocalDofWithGhost());
    a55 += integrate(_range=elements(mesh), _expr=k*gradt(T)*trans(grad(T)));
    a55 += integrate(_range=markedfaces(mesh,"R"), _expr=doption("h")*idt(T)*id(q));
}

} // FeelModels
} // Feel
