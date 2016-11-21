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

namespace Feel {


inline
po::options_description
makeOptions()
{
    po::options_description testhdivoptions( "test h_div options" );
    testhdivoptions.add_options()
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ( "k", po::value<std::string>()->default_value( "-1" ), "k" )
        ( "p_exact", po::value<std::string>()->default_value( "(1/(2*Pi*Pi))*sin(Pi*x)*sin(Pi*y):x:y" ), "p exact" )
        // ( "p_exacts", po::value<std::vector<std::string> >(), "p exact to test" )
        ( "d1", po::value<double>()->default_value( 0.5 ), "d1 (stabilization term)" )
        ( "d2", po::value<double>()->default_value( 0.5 ), "d2 (stabilization term)" )
        ( "d3", po::value<double>()->default_value( 0.5 ), "d3 (stabilization term)" )
        ( "nb_refine", po::value<int>()->default_value( 4 ), "nb_refine" )
        ( "use_hypercube", po::value<bool>()->default_value( true ), "use hypercube or a given geometry" )
        // begin{dp}
        /*
         Stabilization function for hybridized methods.
         */
        ( "tau_constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "tau_order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        // end{dp}
        ;
    return testhdivoptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "hdg" ,
                     "hdg" ,
                     "0.1",
                     "convergence test for Hdg problem (using Hdiv conforming elts)",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Daniele Prada", "developer", "daniele.prada85@gmail.com", "" );
    return about;

}

template<int Dim, int OrderP>
class Hdg
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
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;

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

    using Vh_t =  Pdhv_type<mesh_type,OrderP>;
    using Vh_ptr_t =  Pdhv_ptrtype<mesh_type,OrderP>;
    using Wh_t =  Pdh_type<mesh_type,OrderP>;
    using Wh_ptr_t =  Pdh_ptrtype<mesh_type,OrderP>;
    using Mh_t =  Pdh_type<face_mesh_type,OrderP>;
    using Mh_ptr_t =  Pdh_ptrtype<face_mesh_type,OrderP>;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    Hdg()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        exporter( Exporter<mesh_type>::New( Environment::about().appName() ) ),
        M_d1( doption("d1") ),
        M_d2( doption("d2") ),
        M_d3( doption("d3") ),
        // stabilisation parameter
        M_tau_constant( doption("tau_constant") ),
        M_tau_order( ioption("tau_order") )
    {
        this->changeRepository( boost::format( "benchmark_hdg/%1%/h_%2%/" )
                                % this->about().appName()
                                % doption("hsize")
                                );
    }

    /**
     * run the application
     */
    void convergence();

    template<typename MatrixType, typename VectorType, typename VhType, typename WhType, typename MhType,
             typename ExprP, typename ExprGradP>
    void
    assemble_A_and_F( MatrixType A,
                      VectorType F,
                      VhType Vh,
                      WhType Wh,
                      MhType Mh,
                      ExprP p_exact,
                      ExprGradP gradp_exact );

private:
    //! linear algebra backend
    backend_ptrtype M_backend;
    //! mesh characteristic size
    double meshSize;
    //! exporter factory
    export_ptrtype exporter;

    double M_d1, M_d2, M_d3;

    // begin{dp}
    double M_tau_constant;
    int    M_tau_order;
    // end{dp}
}; //Hdg

template<int Dim, int OrderP>
void
Hdg<Dim, OrderP>::convergence()
{
    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;

    auto K = expr(soption("k"));
    auto lambda = cst(1.)/K;

    // Exact solutions
    auto p_exact = expr(soption("p_exact"));
    auto gradp_exact = grad<Dim>(p_exact);
    auto u_exact = -K*trans(gradp_exact);
    auto f = -K*laplacian(p_exact);

    cout << "k : " << K /*<< "\tlambda : " << lambda*/ << std::endl;
    cout << "p : " << p_exact << std::endl;
    cout << "gradp : " << gradp_exact << std::endl;
    //cout << "u : " << u_exact << std::endl;
    cout << "f : " << laplacian(p_exact) << std::endl;


    // Coeff for stabilization terms
    // auto d1 = cst(M_d1);
    // auto d2 = cst(M_d2);
    // auto d3 = cst(M_d3);
    auto tau_constant = cst(M_tau_constant);

    std::ofstream out;
    MasterStream cvg( out );
    cvg.open("data.csv");
    cvg << "hsize" << "\t" << "nDofu" << "\t" << "l2err_u"  << "\t" << "nDofp"  << "\t" << "l2err_p" << "\n";

    double current_hsize = meshSize;
    for(int i=0; i<ioption("nb_refine"); i++)
    {
        mesh_ptrtype mesh;
        std::string mesh_name;
        if( boption("use_hypercube"))
            mesh_name = (boost::format( "%1%-%2%D" ) % "hypercube" % Dim ).str();
        else
        {
            Feel::fs::path mypath(soption( _name="gmsh.filename" ));
            mesh_name = mypath.stem().string();
        }

        tic();
        // be careful to delete the mesh if you rerun the benchmarks
        if( !fs::exists(mesh_name+".msh") )
        {
            //std::cout << "createGMSHmesh" << std::endl;
            LOG(INFO) << "[Hdg] Mesh has been created \n";
            if( boption("use_hypercube"))
                mesh = createGMSHMesh( _mesh=new mesh_type,
                                       _desc=domain( _name = mesh_name ,
                                                     _shape = "hypercube",
                                                     _usenames = true,
                                                     _dim = Dim,
                                                     _h = meshSize,
                                                     _xmin=0,_xmax=2,
                                                     _ymin=0,_ymax=2,
                                                     _zmin=0,_zmax=2 ) );
            else
                mesh = loadMesh( new mesh_type);
        }
        else
        {
            //std::cout << "loadGMSHmesh" << std::endl;
            LOG(INFO) << "[Hdg] Mesh has been loaded (refine level = " << i << ") \n";
            mesh = loadGMSHMesh( _mesh=new mesh_type,
                                 _filename=mesh_name+".msh",
                                 _refine=i,
                                 _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
        }
        toc("mesh",true);

        // ****** Hybrid-mixed formulation ******
        // We treat Vh, Wh, and Mh separately
        tic();

        Vh_ptr_t Vh = Pdhv<OrderP>( mesh, true );
        Wh_ptr_t Wh = Pdh<OrderP>( mesh, true );
        auto face_mesh = createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );
        Mh_ptr_t Mh = Pdh<OrderP>( face_mesh,true );

        toc("spaces",true);


        size_type nFaceInParallelMesh = nelements(faces(mesh),true) - nelements(interprocessfaces(mesh),true)/2;
        CHECK( nelements(elements(face_mesh),true) == nFaceInParallelMesh  ) << "something wrong with face mesh" << nelements(elements(face_mesh),true) << " " << nFaceInParallelMesh;
        auto Xh = Pdh<0>(face_mesh);
        auto uf = Xh->element(cst(1.));
        CHECK( uf.size() == nFaceInParallelMesh ) << "check faces failed " << uf.size() << " " << nFaceInParallelMesh;

        cout << "Vh<" << OrderP << "> : " << Vh->nDof() << std::endl
             << "Wh<" << OrderP << "> : " << Wh->nDof() << std::endl
             << "Mh<" << OrderP << "> : " << Mh->nDof() << std::endl;

        auto up = Vh->elementPtr( "u" );
        auto pp = Wh->elementPtr( "p" );
        auto phatp = Mh->elementPtr( "phat" );

        auto u = Vh->element( "u" );
        auto v = Vh->element( "v" );
        auto p = Wh->element( "p" );
        auto q = Wh->element( "q" );
        auto w = Wh->element( "w" );
        auto phat = Mh->element( "phat" );
        auto l = Mh->element( "lambda" );

        // Number of dofs associated with each space U and P
        auto nDofu = u.functionSpace()->nDof();
        auto nDofp = p.functionSpace()->nDof();
        auto nDofphat = phat.functionSpace()->nDof();

        tic();
        // build csr graph blocks from set of function spaces Vh Wh and Mh
#if 0
        auto A = backend()->newBlockMatrix(_block=csrGraphBlocks(Vh,Wh,Mh));
#else
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
        auto A = backend()->newBlockMatrix(_block=hdg_graph);
#endif


        // build vector blocks from sub-vector of size given by set of function spaces
        BlocksBaseVector<double> hdg_vec(3);
        hdg_vec(0,0) = backend()->newVector( Vh );
        hdg_vec(1,0) = backend()->newVector( Wh );
        hdg_vec(2,0) = backend()->newVector( Mh );
        auto F = backend()->newBlockVector(_block=hdg_vec, _copy_values=false);

        // build vector blocks from concatenation of subvectors up pp and phatp
        auto hdg_sol = vectorBlocks(up,pp,phatp);
        auto U = backend()->newBlockVector(_block=hdg_sol, _copy_values=false);

        assemble_A_and_F( A, F, Vh, Wh, Mh, p_exact, gradp_exact );
        toc("matrices",true);

        tic();
        backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );

        hdg_sol.localize(U);
        toc("solve",true);

        cout << "[Hdg] solve done" << std::endl;

        // ****** Compute error ******


        tic();

        bool has_dirichlet = nelements(markedfaces(mesh,"Dirichlet"),true) >= 1;


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
        toc("error");

        cout << "[" << i << "]||u_exact - u||_L2 = " << l2err_u << std::endl;
        cout << "[" << i << "]||p_exact - p||_L2 = " << l2err_p << std::endl;
        cvg << current_hsize
            << "\t" << Vh->nDof() << "\t" << l2err_u
            << "\t" << Wh->nDof() << "\t" << l2err_p << std::endl;


        tic();
        std::string exportName =  ( boost::format( "%1%-refine-%2%" ) % this->about().appName() % i ).str();
        std::string uName = ( boost::format( "velocity-refine-%1%" ) % i ).str();
        std::string u_exName = ( boost::format( "velocity-ex-refine-%1%" ) % i ).str();
        std::string pName = ( boost::format( "potential-refine-%1%" ) % i ).str();
        std::string p_exName = ( boost::format( "potential-ex-refine-%1%" ) % i ).str();

        v.on( _range=elements(mesh), _expr=u_exact );
        q.on( _range=elements(mesh), _expr=p_exact );
        export_ptrtype exporter_cvg( export_type::New( exportName ) );

        exporter_cvg->step( i )->setMesh( mesh );
        exporter_cvg->step( i )->add( uName, *up );
        exporter_cvg->step( i )->add( pName, *pp );
        exporter_cvg->step( i )->add( u_exName, v );
        exporter_cvg->step( i )->add( p_exName, q );
        exporter_cvg->save();

        current_hsize /= 2.0;
        toc("export");

    }

    cvg.close();
}

template<int Dim, int OrderP>
template<typename MatrixType, typename VectorType, typename VhType, typename WhType, typename MhType,
         typename ExprP, typename ExprGradP>
void
Hdg<Dim, OrderP>::assemble_A_and_F( MatrixType A,
                                    VectorType F,
                                    VhType Vh,
                                    WhType Wh,
                                    MhType Mh,
                                    ExprP p_exact,
                                    ExprGradP gradp_exact )
{
    auto K = expr(soption("k"));
    auto lambda = cst(1.)/K;
    auto tau_constant = cst(M_tau_constant);
    auto f = -K*laplacian(p_exact);
    auto mesh = Vh->mesh();

    auto u = Vh->element( "u" );
    auto v = Vh->element( "v" );
    auto p = Wh->element( "p" );
    auto q = Wh->element( "q" );
    auto w = Wh->element( "w" );
    auto phat = Mh->element( "phat" );
    auto l = Mh->element( "lambda" );

    // Building the RHS
    //
    // This is only a part of the RHS - how to build the whole RHS? Is it right to
    // imagine we moved it to the left? SKIPPING boundary conditions for the moment.
    // How to identify Dirichlet/Neumann boundaries?
    auto rhs2 = form1( _test=Wh, _vector=F,
                       _rowstart=1 );

    rhs2 += integrate(_range=elements(mesh),
                      _expr=f*id(w));

    cout << "rhs2 works fine" << std::endl;

    // begin dp: changed signs (to move terms to the left)
    auto rhs3 = form1( _test=Mh, _vector=F,
                       _rowstart=2);
    rhs3 += integrate(_range=markedfaces(mesh,"Neumann"),
                      _expr=-id(l)*K*gradp_exact*N());
    rhs3 += integrate(_range=markedfaces(mesh,"Dirichlet"),
                      _expr=id(l)*p_exact);
    // end dp

    cout << "rhs3 works fine" << std::endl;

    auto a11 = form2( _trial=Vh, _test=Vh,_matrix=A );
    a11 += integrate(_range=elements(mesh),_expr=(trans(lambda*idt(u))*id(v)) );

    cout << "a11 works fine" << std::endl;

    auto a12 = form2( _trial=Wh, _test=Vh,_matrix=A,
                      _rowstart=0, _colstart=1 );
    a12 += integrate(_range=elements(mesh),_expr=-(idt(p)*div(v)));

    cout << "a12 works fine" << std::endl;

    // begin dp: added extended pattern, multiplied by 0.5 when integrating over internalfaces
    auto a13 = form2( _trial=Mh, _test=Vh,_matrix=A,
                      _rowstart=0, _colstart=2);

    a13 += integrate(_range=internalfaces(mesh),
                     _expr=( idt(phat)*leftface(trans(id(v))*N())+
                             idt(phat)*rightface(trans(id(v))*N())) );
    a13 += integrate(_range=boundaryfaces(mesh),
                     _expr=idt(phat)*trans(id(v))*N());
    // end dp

    cout << "a13 works fine" << std::endl;

    // begin dp: added extended pattern
    auto a21 = form2( _trial=Vh, _test=Wh,_matrix=A,
                      _rowstart=1, _colstart=0);

    // end dp
    a21 += integrate(_range=elements(mesh),_expr=(-grad(w)*idt(u)));
    cout << " . a211 ok" << std::endl;
    a21 += integrate(_range=internalfaces(mesh),
                     _expr=( leftface(id(w))*leftfacet(trans(idt(u))*N()) ) );
    cout << " . a212l ok" << std::endl;
    a21 += integrate(_range=internalfaces(mesh),
                     _expr=(rightface(id(w))*rightfacet(trans(idt(u))*N())) );
    cout << " . a212r ok" << std::endl;
    a21 += integrate(_range=boundaryfaces(mesh),
                     _expr=(id(w)*trans(idt(u))*N()));
    cout << " . a213 ok" << std::endl;
    cout << "a21 works fine" << std::endl;

    // begin dp: added extended pattern
    auto a22 = form2( _trial=Wh, _test=Wh,_matrix=A,
                      _rowstart=1, _colstart=1 );

    // end dp
    a22 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant *
                     ( leftfacet( pow(h(),M_tau_order)*idt(p))*leftface(id(w)) +
                       rightfacet( pow(h(),M_tau_order)*idt(p))*rightface(id(w) )));
    a22 += integrate(_range=boundaryfaces(mesh),
                     _expr=(tau_constant * pow(h(),M_tau_order)*id(w)*idt(p)));

    cout << "a22 works fine" << std::endl;

    // begin dp: added extended pattern, multiplied by 0.5
    auto a23 = form2( _trial=Mh, _test=Wh,_matrix=A,
                      _rowstart=1, _colstart=2);
    a23 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * idt(phat) *
                     ( leftface( pow(h(),M_tau_order)*id(w) )+
                       rightface( pow(h(),M_tau_order)*id(w) )));
    // end dp
    a23 += integrate(_range=boundaryfaces(mesh),
                     _expr=-tau_constant * idt(phat) * pow(h(),M_tau_order)*id(w) );

    cout << "a23 works fine" << std::endl;

    // begin dp: added extended pattern, multiplied by 0.5
    auto a31 = form2( _trial=Vh, _test=Mh,_matrix=A,
                      _rowstart=2, _colstart=0);
    a31 += integrate(_range=internalfaces(mesh),
                     _expr=( id(l)*(leftfacet(trans(idt(u))*N())+
                                    rightfacet(trans(idt(u))*N())) ) );
    // end dp

    // BC
    a31 += integrate(_range=markedfaces(mesh,"Neumann"),
                     _expr=( id(l)*(trans(idt(u))*N()) ));

    cout << "a31 works fine" << std::endl;

    // begin dp: added extended pattern, mulitplied by 0.5
    auto a32 = form2( _trial=Wh, _test=Mh,_matrix=A,
                      _rowstart=2, _colstart=1);
    a32 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant * id(l) * ( leftfacet( pow(h(),M_tau_order)*idt(p) )+
                                                    rightfacet( pow(h(),M_tau_order)*idt(p) )));
    // end do
    a32 += integrate(_range=markedfaces(mesh,"Neumann"),
                     _expr=tau_constant * id(l) * ( pow(h(),M_tau_order)*idt(p) ) );

    cout << "a32 works fine" << std::endl;

    auto a33 = form2(_trial=Mh, _test=Mh,_matrix=A,
                     _rowstart=2, _colstart=2);
    // begin dp: mulitplied by 0.25
    a33 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * idt(phat) * id(l) * ( leftface( pow(h(),M_tau_order) )+
                                                                 rightface( pow(h(),M_tau_order) )));
    // end dp
    a33 += integrate(_range=markedfaces(mesh,"Neumann"),
                     _expr=-tau_constant * idt(phat) * id(l) * ( pow(h(),M_tau_order) ) );
    a33 += integrate(_range=markedfaces(mesh,"Dirichlet"),
                     _expr=idt(phat) * id(l) );

    cout << "a33 works fine" << std::endl;

}


} // Feel
