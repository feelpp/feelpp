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
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelmesh/complement.hpp>

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
        ( "nb_refine", po::value<int>()->default_value( 4 ), "nb_refine" )
        ( "use_hypercube", po::value<bool>()->default_value( false ), "use hypercube or a given geometry" )
        // begin{dp}
        /*
         Stabilization function for hybridized methods.
         */
        ( "tau_constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "tau_order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( "nb_ibc", po::value<int>()->default_value( 1 ), "number of integral boundary condition" )
        // end{dp}
        ;
    return testhdivoptions.add( Feel::feel_options() ).add( backend_options("sc")).add( backend_options("hdg"));
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

    using Vh_t =  Pdhv_type<mesh_type,OrderP>;
    using Vh_ptr_t =  Pdhv_ptrtype<mesh_type,OrderP>;
    using Wh_t =  Pdh_type<mesh_type,OrderP>;
    using Wh_ptr_t =  Pdh_ptrtype<mesh_type,OrderP>;
    using Mh_t =  Pdh_type<face_mesh_type,OrderP>;
    using Mh_ptr_t =  Pdh_ptrtype<face_mesh_type,OrderP>;

    using Ch_t = Pch_type<face_mesh_type,0>;
    using Ch_ptr_t = Pch_ptrtype<face_mesh_type,0>;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef std::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    Hdg()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        exporter( Exporter<mesh_type>::New( Environment::about().appName() ) ),
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

private:
    //! linear algebra backend
    backend_ptrtype M_backend;
    //! mesh characteristic size
    double meshSize;
    //! exporter factory
    export_ptrtype exporter;

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
    auto tau_constant = cst(M_tau_constant);

    std::ofstream out;
    MasterStream cvg( out );
    cvg.open("data.csv");
    cvg << "hsize" << "\t" << "nDofu" << "\t" << "l2err_u"  << "\t" << "nDofp"  << "\t" << "l2err_p" << "\n";

    double current_hsize = meshSize;
    for(int i=0; i<ioption("nb_refine"); i++)
    {
        backend(_rebuild=true);
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

        auto complement_integral_bdy = complement(faces(mesh),[&mesh]( auto const& e ) {
                for( int i = 0; i < ioption("nb_ibc"); ++i )
                {
                    auto ind = i == 0 ? "" : std::to_string(i+1);
                    std::string marker = boost::str(boost::format("Ibc%1%") % ind );
                    if ( e.hasMarker() && e.marker().value() == mesh->markerName( marker ) )
                        return true;
                }
                return false;
            });
        auto gammaMinusIntegral = complement(boundaryfaces(mesh), [&mesh]( auto const& e ) {
                for( int i = 0; i < ioption("nb_ibc"); ++i )
                {
                    auto ind = i == 0 ? "" : std::to_string(i+1);
                    std::string marker = boost::str(boost::format("Ibc%1%") % ind );
                    if ( e.hasMarker() && e.marker().value() == mesh->markerName( marker ) )
                        return true;
                }
                return false;
            });
        auto face_mesh = createSubmesh( mesh, complement_integral_bdy, EXTRACTION_KEEP_MESH_RELATION, 0 );
        std::vector<std::string> ibc_markers(ioption("nb_ibc"));
        for( int i = 0; i < ioption("nb_ibc"); ++i )
        {
            auto ind = i == 0 ? "" : std::to_string(i+1);
            std::string marker = boost::str(boost::format("Ibc%1%") % ind );
            ibc_markers.push_back( marker );
        }
        auto ibc_mesh = createSubmesh( mesh, markedfaces(mesh,ibc_markers), EXTRACTION_KEEP_MESH_RELATION, 0 );
        toc("mesh",true);

        // ****** Hybrid-mixed formulation ******
        // We treat Vh, Wh, and Mh separately
        tic();

        Vh_ptr_t Vh = Pdhv<OrderP>( mesh, true );
        Wh_ptr_t Wh = Pdh<OrderP>( mesh, true );
        Mh_ptr_t Mh = Pdh<OrderP>( face_mesh,true );

        Ch_ptr_t Ch = Pch<0>(ibc_mesh, true );

        toc("spaces",true);

        LOG(INFO)<< "number elemnts:" << mesh->numGlobalElements() << std::endl;
        LOG(INFO) << "number local elemnts:" << mesh->numElements() << std::endl;
        LOG(INFO) << "number local active elements:" << nelements(elements(mesh),false) << std::endl;
        size_type nFaceInParallelMesh = nelements(faces(mesh),true) - nelements(interprocessfaces(mesh),true)/2;
        //CHECK( nelements(elements(face_mesh),true) == nFaceInParallelMesh  ) << "something wrong with face mesh" << nelements(elements(face_mesh),true) << " " << nFaceInParallelMesh;
        auto Xh = Pdh<0>(face_mesh);
        auto uf = Xh->element(cst(1.));
        //CHECK( uf.size() == nFaceInParallelMesh ) << "check faces failed " << uf.size() << " " << nFaceInParallelMesh;

        cout << "Vh<" << OrderP << "> : " << Vh->nDof() << std::endl
             << "Wh<" << OrderP << "> : " << Wh->nDof() << std::endl
             << "Mh<" << OrderP << "> : " << Mh->nDof() << std::endl
             << "Ch<" << 0 << "> : " << Ch->nDof() << std::endl;

        auto u = Vh->element( "u" );
        auto v = Vh->element( "v" );
        auto p = Wh->element( "p" );
        auto q = Wh->element( "q" );
        auto w = Wh->element( "w" );
        auto phat = Mh->element( "phat" );
        auto l = Mh->element( "lambda" );
        auto nu = Ch->element( "nu" );
        auto uI = Ch->element( "uI" );

        // Number of dofs associated with each space U and P
        auto nDofu = u.functionSpace()->nDof();
        auto nDofp = p.functionSpace()->nDof();
        auto nDofphat = phat.functionSpace()->nDof();

        tic();
        auto ibcSpaces = std::make_shared<ProductSpace<Ch_ptr_t,true> >( ioption("nb_ibc"), Ch);
        auto ps = product2( ibcSpaces, Vh, Wh, Mh );

        solve::strategy strategy = boption("sc.condense")?solve::strategy::static_condensation:solve::strategy::monolithic;

        auto a = blockform2( ps, strategy, backend() );
        auto rhs = blockform1( ps, strategy, backend() );

        auto K = expr(soption("k"));
        auto lambda = cst(1.)/K;
        auto tau_constant = cst(M_tau_constant);
        auto f = -K*laplacian(p_exact);

        auto sc_param = boption("sc.condense") ? 0.5 : 1.0;

        // Building the RHS
        rhs(1_c) += integrate(_range=elements(mesh),
                              _expr=-f*id(w));

        cout << "rhs2 works fine" << std::endl;

        rhs(2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                              _expr=-id(l)*K*gradp_exact*N());
        rhs(2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                              _expr=id(l)*p_exact);

        cout << "rhs3 works fine" << std::endl;

        for( int i = 0; i < ioption("nb_ibc"); ++i )
        {
            std::string marker = boost::str(boost::format("Ibc%1%") % (i == 0 ? "" : std::to_string(i+1)) );
            double intjn = integrate(_range=markedfaces(mesh,"Ibc"),
                                     _expr=inner(u_exact,N())).evaluate()(0,0);
            double meas = integrate( _range=markedfaces(mesh, "Ibc"),
                                     _expr=cst(1.0)).evaluate()(0,0);
            rhs(3_c, i) += integrate( _range=markedfaces(mesh, "Ibc"),
                                     _expr=intjn*id(nu)/meas );

            cout << "rhs(4[" << i << "]) works fine on " << marker << " for g=" << intjn << std::endl;
        }

        //
        // First row a(0_c,:)
        //
        a(0_c,0_c) += integrate(_range=elements(mesh),_expr=(trans(lambda*idt(u))*id(v)) );

        cout << "a11 works fine" << std::endl;

        a(0_c,1_c) += integrate(_range=elements(mesh),_expr=-(idt(p)*div(v)));

        cout << "a12 works fine" << std::endl;

        a(0_c,2_c) += integrate(_range=internalfaces(mesh),
                                _expr=( idt(phat)*leftface(trans(id(v))*N())+
                                        idt(phat)*rightface(trans(id(v))*N())) );

        a(0_c,2_c) += integrate(_range=gammaMinusIntegral,
                                _expr=idt(phat)*trans(id(v))*N());

        cout << "a13 works fine" << std::endl;

        //
        // Second row a(1_c,:)
        //
#if 0
        a(1_c,0_c) += integrate(_range=elements(mesh),_expr=(-grad(w)*idt(u)));
        a(1_c,0_c) += integrate(_range=internalfaces(mesh),
                                _expr=( leftface(id(w))*leftfacet(trans(idt(u))*N()) ) );
        a(1_c,0_c) += integrate(_range=internalfaces(mesh),
                                _expr=(rightface(id(w))*rightfacet(trans(idt(u))*N())) );
        a(1_c,0_c) += integrate(_range=boundaryfaces(mesh),
                                _expr=(id(w)*trans(idt(u))*N()));
#else
        a(1_c,0_c) += integrate(_range=elements(mesh),_expr=-(id(w)*divt(u)));
#endif
        cout << "a21 works fine" << std::endl;

        a(1_c,1_c) += integrate(_range=internalfaces(mesh),
                                _expr=-tau_constant *
                                ( leftfacet( pow(h(),M_tau_order)*idt(p))*leftface(id(w)) +
                                  rightfacet( pow(h(),M_tau_order)*idt(p))*rightface(id(w) )));
        a(1_c,1_c) += integrate(_range=boundaryfaces(mesh),
                                _expr=-(tau_constant * pow(h(),M_tau_order)*id(w)*idt(p)));

        cout << "a22 works fine" << std::endl;

        a(1_c,2_c) += integrate(_range=internalfaces(mesh),
                                _expr=tau_constant * idt(phat) *
                                ( leftface( pow(h(),M_tau_order)*id(w) )+
                                  rightface( pow(h(),M_tau_order)*id(w) )));

        a(1_c,2_c) += integrate(_range=gammaMinusIntegral,
                                _expr=tau_constant * idt(phat) * pow(h(),M_tau_order)*id(w) );

        cout << "a23 works fine" << std::endl;

        //
        // Third row a(2_c,:)
        //
        a(2_c,0_c) += integrate(_range=internalfaces(mesh),
                                _expr=( id(l)*(leftfacet(trans(idt(u))*N())+
                                               rightfacet(trans(idt(u))*N())) ) );
        // BC
        a(2_c,0_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                                _expr=( id(l)*(trans(idt(u))*N()) ));

        cout << "a31 works fine" << std::endl;

        a(2_c,1_c) += integrate(_range=internalfaces(mesh),
                                _expr=tau_constant * id(l) * ( leftfacet( pow(h(),M_tau_order)*idt(p) )+
                                                               rightfacet( pow(h(),M_tau_order)*idt(p) )));
        a(2_c,1_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                                _expr=tau_constant * id(l) * ( pow(h(),M_tau_order)*idt(p) ) );

        cout << "a32 works fine" << std::endl;

        a(2_c,2_c) += integrate(_range=internalfaces(mesh),
                                _expr=-sc_param*tau_constant * idt(phat) * id(l) * ( leftface( pow(h(),M_tau_order) )+
                                                                            rightface( pow(h(),M_tau_order) )));
        a(2_c,2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                                _expr=-tau_constant * idt(phat) * id(l) * ( pow(h(),M_tau_order) ) );
        a(2_c,2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                                _expr=idt(phat) * id(l) );

        cout << "a33 works fine" << std::endl;

        for( int i = 0; i < ioption("nb_ibc"); ++i )
        {
            std::string marker = boost::str(boost::format("Ibc%1%") % (i == 0 ? "" : std::to_string(i+1)) );
            cout << "applying ibc condition on a for marker " << marker << std::endl;
            // <lambda, v.n>_Gamma_I
            a( 0_c, 3_c, 0, i ) += integrate( _range=markedfaces(mesh,"Ibc"),
                                              _expr= idt(uI) * (trans(id(v))*N()) );

            cout << "a(1,4[" << i << "]) works fine" << std::endl;

            // <lambda, tau w>_Gamma_I
            a( 1_c, 3_c, 1, i ) += integrate( _range=markedfaces(mesh,"Ibc"),
                                              _expr=tau_constant * idt(uI) * id(w) * ( pow(h(),M_tau_order)) );

            cout << "a(2,4[" << i << "]) works fine" << std::endl;

            // <j.n, m>_Gamma_I
            a( 3_c, 0_c, i, 0 ) += integrate( _range=markedfaces(mesh,"Ibc"),
                                              _expr=(trans(idt(u))*N()) * id(nu) );

            cout << "a(4[" << i << "],1) works fine" << std::endl;

            // <tau p, m>_Gamma_I
            a( 3_c, 1_c, i, 1 ) += integrate( _range=markedfaces(mesh,"Ibc"),
                                              _expr=tau_constant *idt(p)  * id(nu)* ( pow(h(),M_tau_order)) );

            cout << "a(4[" << i << "],2) works fine" << std::endl;

            // -<lambda2, m>_Gamma_I
            a( 3_c, 3_c, i, i ) += integrate( _range=markedfaces(mesh,"Ibc"),
                                              _expr=-tau_constant * id(nu) *idt(uI)* ( pow(h(),M_tau_order)) );

            cout << "a(4[" << i << "],4[" << i << "]) works fine" << std::endl;

        }
        toc("matrices",true);

        tic();
        auto U=ps.element();

        a.solve( _solution=U, _rhs=rhs, _condense=boption("sc.condense"), _name="hdg");
        toc("solve",true);

        cout << "[Hdg] solve done" << std::endl;

        // ****** Compute error ******
        auto up = U(0_c);
        auto pp = U(1_c);
        tic();

        bool has_dirichlet = nelements(markedfaces(mesh,"Dirichlet"),true) >= 1;


        auto l2err_u = normL2( _range=elements(mesh), _expr=u_exact - idv(up) );
        double l2err_p = 1e+30;
        if ( has_dirichlet )
        {
            l2err_p = normL2( _range=elements(mesh), _expr=p_exact - idv(pp) );
        }
        else
        {
            auto mean_p_exact = mean( elements(mesh), p_exact )(0,0);
            auto mean_p = mean( elements(mesh), idv(pp) )(0,0);
            l2err_p = normL2( elements(mesh),
                              (p_exact - cst(mean_p_exact)) - (idv(pp) - cst(mean_p)) );
        }
        toc("error");

        cout << "[" << i << "]||u_exact - u||_L2 = " << l2err_u << std::endl;
        cout << "[" << i << "]||p_exact - p||_L2 = " << l2err_p << std::endl;
        double intjnex = integrate(_range=markedfaces(mesh,"Ibc"),
                                   _expr=inner(u_exact,N())).evaluate()(0,0);
        double intjn = integrate(_range=markedfaces(mesh,"Ibc"),
                                 _expr=inner(idv(up),N())).evaluate()(0,0);
        cout << "[" << i << "]int u_exact.n  = " << intjnex << std::endl;
        cout << "[" << i << "]int u.n  = " << intjn << std::endl;
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
        exporter_cvg->step( i )->add( uName, U(0_c) );
        exporter_cvg->step( i )->add( pName, U(1_c) );
        exporter_cvg->step( i )->add( u_exName, v );
        exporter_cvg->step( i )->add( p_exName, q );
        exporter_cvg->save();

        current_hsize /= 2.0;
        toc("export");

    }

    cvg.close();
}



} // Feel
