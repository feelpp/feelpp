/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel++ library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Daniele Prada <daniele.prada85@gmail.com>
			 Lorenzo Sala <sala@unistra.fr>
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

namespace Feel {


inline
po::options_description
makeOptions()
{
    po::options_description testhdivoptions( "test h_div options" );
    testhdivoptions.add_options()
        ( "convtest", po::value<bool>()->default_value( 1 ), "1: convergence test, 0:otherwise" )
        ( "hsize", po::value<double>()->default_value( 0.8 ), "mesh size" )
        ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ( "lambda", po::value<std::string>()->default_value( "1" ), "lambda" )
        ( "mu", po::value<std::string>()->default_value( "1" ), "mu" )
        ( "v_testfun", po::value<std::string>()->default_value( "{1,0,0,1}:x:y" ), "value for the test function" )
        ( "w_testfun", po::value<std::string>()->default_value( "{1,0,0,1}:x:y" ), "value for the test function" )
        ( "u_exact", po::value<std::string>()->default_value( "{(1/(2*Pi*Pi))*sin(Pi*x)*cos(Pi*y),(1/(2*Pi*Pi))*cos(Pi*x)*sin(Pi*y)}:x:y" ), "u exact" )
        ( "f", po::value<std::string>()->default_value( "{-3.0*sin(pi*x)*cos(pi*y),-3.0*sin(pi*y)*cos(pi*x)}:x:y"  ), "divergence of the stress tensor")
        ( "load", po::value<std::string>()->default_value( "{0,0,0,0}:x:y"  ), "load")
        // ( "u_exacts", po::value<std::vector<std::string> >(), "u exact to test" )
        ( "hface", po::value<int>()->default_value( 0 ), "hface" )
        ( "nb_refine", po::value<int>()->default_value( 4 ), "nb_refine" )
        ( "use_hypercube", po::value<bool>()->default_value( true ), "use hypercube or a given geometry" )
        ( "tau_constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "tau_order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ;
    return testhdivoptions.add( Feel::feel_options() ).add( backend_options("sc"));
}

inline
AboutData
makeAbout()
{
    AboutData about( "hdg" ,
                     "hdg" ,
                     "0.1",
                     "convergence test for an HDG method for linear elasticity",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Daniele Prada", "developer", "daniele.prada85@gmail.com", "" );
    about.addAuthor( "Lorenzo Sala", "developer", "sala@unistra.fr", "" );

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

    using Vh_t =  Pdhms_type<mesh_type,OrderP>;
    using Vh_ptr_t =  Pdhms_ptrtype<mesh_type,OrderP>;
    using Wh_t =  Pdhv_type<mesh_type,OrderP>;
    using Wh_ptr_t =  Pdhv_ptrtype<mesh_type,OrderP>;
    using Mh_t =  Pdhv_type<face_mesh_type,OrderP>;
    using Mh_ptr_t =  Pdhv_ptrtype<face_mesh_type,OrderP>;
    using M0h_t =  Pdh_type<face_mesh_type,0>;
    using M0h_ptr_t =  Pdh_ptrtype<face_mesh_type,0>;

    using product_space_type = ProductSpaces<Vh_ptr_t,Wh_ptr_t,Mh_ptr_t>;
    using product_space_ptrtype = boost::shared_ptr<product_space_type>;

    using blockform2_type = BlockBilinearForm<ProductSpaces< Vh_ptr_t, Wh_ptr_t, Mh_ptr_t > &>;
    using blockform2_ptrtype = boost::shared_ptr<blockform2_type>;
    using blockform1_type = BlockLinearForm<ProductSpaces< Vh_ptr_t, Wh_ptr_t, Mh_ptr_t > &>;
    using blockform1_ptrtype = boost::shared_ptr<blockform1_type>;
 
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
        // stabilisation parameter
        M_tau_constant( doption("tau_constant") ),
        M_tau_order( ioption("tau_order") )
    {
        BOOST_ASSERT(OrderP >= 1);
        if(Environment::isMasterRank())
            std::cout << "[TestHDiv]\n";
        this->changeRepository( boost::format( "benchmark_hdg/%1%/h_%2%/" )
                                % this->about().appName()
                                % doption("hsize")
                                );
    }

    /**
     * run the application
     */
    void convergence();

    template<typename MatrixType, typename VectorType, typename VhType, typename WhType, typename MhType,typename ExprU, typename ExprSigma, typename ExprDivSigma>
    void
    assemble_A_and_F( MatrixType A,
                      VectorType F,
                      VhType Vh,
                      WhType Wh,
                      MhType Mh,
                      ExprU u_exact,
                      ExprSigma sigma_exact,
                      ExprDivSigma div_sigma_exact );

private:
    //! linear algebra backend
    backend_ptrtype M_backend;
    //! mesh characteristic size
    double meshSize;
    //! exporter factory
    export_ptrtype exporter;

    double M_tau_constant;
    int    M_tau_order;
}; //Hdg

template<int Dim, int OrderP>
void
Hdg<Dim, OrderP>::convergence()
{
    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;

    auto lambda = expr(soption("lambda"));
    auto mu     = expr(soption("mu"));

    // Exact solutions
    auto u_exact = expr<Dim,1>(soption("u_exact"));
    auto gradu_exact = grad(u_exact);
    auto eps_exact   = cst(0.5) * ( gradu_exact + trans(gradu_exact) );
    auto sigma_exact = lambda * trace(eps_exact) * eye<Dim>() + cst(2.) * mu * eps_exact;
    auto f = expr<Dim,1>(soption("f"));

    cout << "lambda : " << lambda      << std::endl;
    cout << "mu     : " << mu          << std::endl;
    cout << "u      : " << u_exact     << std::endl;
    //    cout << "sigma  : " << sigma_exact << std::endl;
    cout << "f      : " << f           << std::endl;


    // Coeff for stabilization terms
    auto tau_constant = cst(M_tau_constant);

    std::ofstream out;
    MasterStream cvg( out );
    cvg.open("data.csv");
    cvg << "hsize" << "\t" << "nDofsigma" << "\t" << "l2err_sigma"  << "\t" << "nDofu"  << "\t" << "l2err_u" << "\n";

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
            {
                mesh = createGMSHMesh( _mesh=new mesh_type,
                                       _desc=domain( _name = mesh_name ,
                                                     _shape = "hypercube",
                                                     _usenames = true,
                                                     _dim = Dim,
                                                     _h = meshSize,
                                                     _xmin=0,_xmax=2,
                                                     _ymin=0,_ymax=2,
                                                     _zmin=0,_zmax=2 ) );
                mesh->addMarkerName( "clamped",( Dim==2 )?1:19, (Dim==2)?1:2);
                mesh->addMarkerName( "tip",( Dim==2)?3:27, (Dim==2)?1:2);
            }
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

        Vh_ptr_t Vh = Pdhms<OrderP>( mesh, true );
        Wh_ptr_t Wh = Pdhv<OrderP>( mesh, true );
        auto face_mesh = createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );
        Mh_ptr_t Mh = Pdhv<OrderP>( face_mesh,true );

        toc("spaces",true);


        size_type nFaceInParallelMesh = nelements(faces(mesh),true) - nelements(interprocessfaces(mesh),true)/2;
        CHECK( nelements(elements(face_mesh),true) == nFaceInParallelMesh  ) << "something wrong with face mesh" << nelements(elements(face_mesh),true) << " " << nFaceInParallelMesh;
        auto Xh = Pdh<0>(face_mesh);
        auto uf = Xh->element(cst(1.));
        CHECK( uf.size() == nFaceInParallelMesh ) << "check faces failed " << uf.size() << " " << nFaceInParallelMesh;

        cout << "Vh<" << OrderP   << "> : " << Vh->nDof() << std::endl
             << "Wh<" << OrderP+1 << "> : " << Wh->nDof() << std::endl
             << "Mh<" << OrderP   << "> : " << Mh->nDof() << std::endl;
		/*
      	auto sigmap = Vh->elementPtr( "sigma" );
        auto up     = Wh->elementPtr( "u" );
        auto uhatp  = Mh->elementPtr( "uhat" );
		*/

        auto sigma = Vh->element( "sigma" );
        auto v     = Vh->element( "v" );
        auto u     = Wh->element( "u" );
        auto w     = Wh->element( "w" );
        auto uhat  = Mh->element( "uhat" );
        auto m     = Mh->element( "m" );

        // Number of dofs associated with each space
        auto nDofsigma = sigma.functionSpace()->nDof();
        auto nDofu     = u.functionSpace()->nDof();
        auto nDofuhat  = uhat.functionSpace()->nDof();


    	auto lambda = expr(soption("lambda"));
    	auto mu     = expr(soption("mu"));
    	auto c1     = cst(0.5)/mu;
   	 	auto c2     = -lambda/(cst(2.) * mu * (cst(Dim)*lambda + cst(2.)*mu));



		tic();
        auto ps = product( Vh, Wh, Mh );
		auto a = blockform2( ps );
		auto rhs = blockform1( ps );

	/*	
        // build the big matrix associated to bilinear form over Vh x Wh x Mh
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
        BlocksBaseVector<double> hdg_vec(3);
        hdg_vec(0,0) = backend()->newVector( Vh );
        hdg_vec(1,0) = backend()->newVector( Wh );
        hdg_vec(2,0) = backend()->newVector( Mh );
        auto F = backend()->newBlockVector(_block=hdg_vec, _copy_values=false);
		
        auto hdg_sol = vectorBlocks(sigmap,up,uhatp);
        auto U = backend()->newBlockVector(_block=hdg_sol, _copy_values=false);
		
    	assemble_A_and_F( A, F, Vh, Wh, Mh, u_exact, sigma_exact, f );
	*/	



    // Building the RHS	
    M0h_ptr_t M0h = Pdh<0>( face_mesh,true );
    auto H     = M0h->element( "H" );
    if ( ioption("hface" ) == 0 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hMax()) );
    else if ( ioption("hface" ) == 1 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hMin()) );
    else if ( ioption("hface" ) == 2 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hAverage()) );
    else
        H.on( _range=elements(face_mesh), _expr=h() );
	
	rhs(1_c) += integrate(_range=elements(mesh),
                      _expr=trans(f)*id(w));

    cout << "rhs2 works fine" << std::endl;

    // in convergence test Neumann condition is given from the displacement and
    // constitutive law
    if ( boption( "convtest" ) )
        rhs(2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                          _expr=trans(id(m))*(sigma_exact*N()));
    else
    {
        auto load = expr<2,2>( soption("load") );
        rhs(2_c) += integrate(_range=markedfaces(mesh,"Tip"),
                          _expr=trans(id(m))*(load*N()));
    }
	
    rhs(2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                      _expr=trans(id(m))*u_exact);


    cout << "rhs3 works fine" << std::endl;

	// Building the matrix
    a( 0_c, 0_c ) +=  integrate(_range=elements(mesh),_expr=(c1*inner(idt(sigma),id(v))) );
    a( 0_c, 0_c ) += integrate(_range=elements(mesh),_expr=(c2*trace(idt(sigma))*trace(id(v))) );

	a( 0_c, 1_c ) += integrate(_range=elements(mesh),_expr=(trans(idt(u))*div(v)));
 
    a( 0_c, 2_c) += integrate(_range=internalfaces(mesh),
                                _expr=-( trans(idt(uhat))*leftface(id(v)*N())+
                                         trans(idt(uhat))*rightface(id(v)*N())) );
    a( 0_c, 2_c) += integrate(_range=boundaryfaces(mesh),
                                _expr=-trans(idt(uhat))*(id(v)*N()));
 
    a( 1_c, 0_c) += integrate(_range=elements(mesh),
                                _expr=(trans(id(w))*divt(sigma)));
	// begin dp: here we need to put the projection of u on the faces
    a( 1_c, 1_c) += integrate(_range=internalfaces(mesh),_expr=-tau_constant *
    	         ( leftfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*leftface(id(w)) +
        	       rightfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*rightface(id(w) )));
 
    a( 1_c, 1_c) += integrate(_range=boundaryfaces(mesh),
    			_expr=-(tau_constant * pow(idv(H),M_tau_order)*trans(idt(u))*id(w)));
 
    a( 1_c, 2_c) += integrate(_range=internalfaces(mesh), _expr=tau_constant *
             ( leftfacet(trans(idt(uhat)))*leftface( pow(idv(H),M_tau_order)*id(w))+
               rightfacet(trans(idt(uhat)))*rightface( pow(idv(H),M_tau_order)*id(w) )));
 
    a( 1_c, 2_c) += integrate(_range=boundaryfaces(mesh),
             _expr=tau_constant * trans(idt(uhat)) * pow(idv(H),M_tau_order)*id(w) );
 
 
    a( 2_c, 0_c) += integrate(_range=internalfaces(mesh),
             _expr=( trans(id(m))*(leftfacet(idt(sigma)*N())+
                     rightfacet(idt(sigma)*N())) ) );
	a( 2_c, 1_c) += integrate(_range=internalfaces(mesh),
                _expr=-tau_constant * trans(id(m)) * (leftfacet( pow(idv(H),M_tau_order)*idt(u) )+
                    rightfacet( pow(idv(H),M_tau_order)*idt(u) )));

    a( 2_c, 2_c) += integrate(_range=internalfaces(mesh),
	    _expr=cst(0.5)*tau_constant * trans(idt(uhat)) * id(m) * ( leftface( pow(idv(H),M_tau_order) )+
                    rightface( pow(idv(H),M_tau_order) )));
    
   	a( 2_c, 2_c) += integrate(_range=markedfaces(mesh,{"Dirichlet"}),
                                _expr=trans(idt(uhat)) * id(m) );
    
    
	a( 2_c, 0_c) += integrate(_range=markedfaces(mesh,{"Neumann","Tip"}),
                        		_expr=( trans(id(m))*(idt(sigma)*N()) ));

    a( 2_c, 1_c) += integrate(_range=markedfaces(mesh,{"Neumann","Tip"}),
 		                        _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tau_order)*idt(u) ) );

    a( 2_c, 2_c) += integrate(_range=markedfaces(mesh,{"Neumann","Tip"}),
                        _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tau_order) ) );





		 
	toc("matrices",true);

	// a.close();
	// rhs.close();

    tic();
    // backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );
    // hdg_sol.localize(U);	
	auto U = ps.element();
	a.solve( _solution=U, _rhs=rhs);
    toc("solve",true);
	cout << "[Hdg] solve done" << std::endl;
		
	auto sigmap = U(0_c);
	auto up = U(1_c);
 	auto uhatp = U(2_c);

	Feel::cout << "sigma: \t" << sigmap << std::endl;
	Feel::cout << "u: \t" << up << std::endl;
	Feel::cout << "uhat: \t" << uhatp << std::endl;
			
    // ****** Compute error ******
    tic();
        bool has_dirichlet = nelements(markedfaces(mesh,"Dirichlet"),true) >= 1;
        BOOST_ASSERT(has_dirichlet);

        /*
        How does Feel++ handle BC on single components? Say, Dirichlet on u_x and
        Neumann on dot(sigma*n, e_y)?
        */

    auto l2err_sigma = normL2( _range=elements(mesh), _expr=sigma_exact - idv(sigmap) );
    auto l2err_u     = normL2( _range=elements(mesh), _expr=u_exact - idv(up) );

    toc("error");

    cout << "[" << i << "]||sigma_exact - sigma||_L2 = " << l2err_sigma << std::endl;
    cout << "[" << i << "]||u_exact - u||_L2 = " << l2err_u << std::endl;
        cvg << current_hsize
            << "\t" << Vh->nDof() << "\t" << l2err_sigma
            << "\t" << Wh->nDof() << "\t" << l2err_u << std::endl;


        tic();
        std::string exportName =  ( boost::format( "%1%-refine-%2%" ) % this->about().appName() % i ).str();
        std::string sigmaName = ( boost::format( "stress-refine-%1%" ) % i ).str();
        std::string sigma_exName = ( boost::format( "stress-ex-refine-%1%" ) % i ).str();
        std::string uName = ( boost::format( "displacement-refine-%1%" ) % i ).str();
        std::string u_exName = ( boost::format( "displacement-ex-refine-%1%" ) % i ).str();

        v.on( _range=elements(mesh), _expr=sigma_exact );
        w.on( _range=elements(mesh), _expr=u_exact );
        export_ptrtype exporter_cvg( export_type::New( exportName ) );

        exporter_cvg->step( i )->setMesh( mesh );
        exporter_cvg->step( i )->add( sigmaName, sigmap );
        exporter_cvg->step( i )->add( uName, up );
        exporter_cvg->step( i )->add( sigma_exName, v );
        exporter_cvg->step( i )->add( u_exName, w );
        exporter_cvg->save();

        current_hsize /= 2.0;
        toc("export");

    }

    cvg.close();
}

template<int Dim, int OrderP>
template<typename MatrixType, typename VectorType, typename VhType, typename WhType, typename MhType,typename ExprU, typename ExprSigma, typename ExprDivSigma>
void
Hdg<Dim, OrderP>::assemble_A_and_F( MatrixType A,
                                    VectorType F,
                                    VhType Vh,
                                    WhType Wh,
                                    MhType Mh,
                                    ExprU u_exact,
                                    ExprSigma sigma_exact,
                                    ExprDivSigma div_sigma_exact )
{

    auto lambda = expr(soption("lambda"));
    auto mu     = expr(soption("mu"));
    auto c1     = cst(0.5)/mu;
    auto c2     = -lambda/(cst(2.) * mu * (cst(Dim)*lambda + cst(2.)*mu));
    auto tau_constant = cst(M_tau_constant);
    auto mesh = Vh->mesh();
    auto face_mesh = Mh->mesh();	//createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );
    M0h_ptr_t M0h = Pdh<0>( face_mesh,true );

    auto sigma = Vh->element( "sigma" );
    auto v     = Vh->element( "v" );
    auto u     = Wh->element( "u" );
    auto w     = Wh->element( "w" );
    auto uhat  = Mh->element( "uhat" );
    auto m     = Mh->element( "m" );
	
    auto H     = M0h->element( "H" );
    if ( ioption("hface" ) == 0 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hMax()) );
    else if ( ioption("hface" ) == 1 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hMin()) );
    else if ( ioption("hface" ) == 2 )
        H.on( _range=elements(face_mesh), _expr=cst(mesh->hAverage()) );
    else
        H.on( _range=elements(face_mesh), _expr=h() );
	

    // Building the RHS

     // Building the RHS
	auto rhs2 = form1( _test=Wh, _vector=F,
                      _rowstart=1 );

   	rhs2 += integrate(_range=elements(mesh),
                     _expr=trans(div_sigma_exact)*id(w));

   cout << "rhs2 works fine" << std::endl;

   auto rhs3 = form1( _test=Mh, _vector=F,
                      _rowstart=2);

   // in convergence test Neumann condition is given from the displacement and
   // constitutive law
   if ( boption( "convtest" ) )
       rhs3 += integrate(_range=markedfaces(mesh,"Neumann"),
                         _expr=trans(id(m))*(sigma_exact*N()));
   else
   {
       auto load = expr<2,2>( soption("load") );
       rhs3 += integrate(_range=markedfaces(mesh,"Tip"),
                         _expr=trans(id(m))*(load*N()));
   }
   rhs3 += integrate(_range=markedfaces(mesh,"Dirichlet"),
                     _expr=trans(id(m))*u_exact);

   cout << "rhs3 works fine" << std::endl;

   auto a11 = form2( _trial=Vh, _test=Vh,_matrix=A );
   auto a11_b1 = form2( _trial=Vh, _test=Vh );
   auto a11_b2 = form2( _trial=Vh, _test=Vh );

   a11 += integrate(_range=elements(mesh),_expr=(c1*inner(idt(sigma),id(v))) );
   a11 += integrate(_range=elements(mesh),_expr=(c2*trace(idt(sigma))*trace(id(v))) );

   a11_b1 = integrate(_range=elements(mesh),_expr=(c1*inner(idt(sigma),id(v))) );
   a11_b2 = integrate(_range=elements(mesh),_expr=(c2*trace(idt(sigma))*trace(id(v))) );

   cout << "a11 works fine" << std::endl;

    auto a12 = form2( _trial=Wh, _test=Vh,_matrix=A,
                      _rowstart=0, _colstart=1 );
    a12 += integrate(_range=elements(mesh),_expr=(trans(idt(u))*div(v)));

    auto a12_b = form2( _trial=Wh, _test=Vh );
    a12_b = integrate(_range=elements(mesh),_expr=(trans(idt(u))*div(v)));

    cout << "a12 works fine" << std::endl;

    auto a13 = form2( _trial=Mh, _test=Vh,_matrix=A,
                      _rowstart=0, _colstart=2);

    a13 += integrate(_range=internalfaces(mesh),
                     _expr=-( trans(idt(uhat))*leftface(id(v)*N())+
                              trans(idt(uhat))*rightface(id(v)*N())) );
    a13 += integrate(_range=boundaryfaces(mesh),
                     _expr=-trans(idt(uhat))*(id(v)*N()));

    auto a21 = form2( _trial=Vh, _test=Wh,_matrix=A,
                      _rowstart=1, _colstart=0);


    a21 += integrate(_range=elements(mesh),_expr=(trans(id(w))*divt(sigma)));

    // begin dp: here we need to put the projection of u on the faces
    auto a22 = form2( _trial=Wh, _test=Wh,_matrix=A,
                      _rowstart=1, _colstart=1);

    a22 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant *
                     ( leftfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*leftface(id(w)) +
                       rightfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*rightface(id(w) )));
    a22 += integrate(_range=boundaryfaces(mesh),
                     _expr=-(tau_constant * pow(idv(H),M_tau_order)*trans(idt(u))*id(w)));


    auto a23 = form2( _trial=Mh, _test=Wh,_matrix=A,
                      _rowstart=1, _colstart=2);
    a23 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant *
                     ( leftfacet(trans(idt(uhat)))*leftface( pow(idv(H),M_tau_order)*id(w))+
                       rightfacet(trans(idt(uhat)))*rightface( pow(idv(H),M_tau_order)*id(w) )));

    a23 += integrate(_range=boundaryfaces(mesh),
                     _expr=tau_constant * trans(idt(uhat)) * pow(idv(H),M_tau_order)*id(w) );

    auto a31 = form2( _trial=Vh, _test=Mh,_matrix=A,
                      _rowstart=2, _colstart=0);
    a31 += integrate(_range=internalfaces(mesh),
                     _expr=( trans(id(m))*(leftfacet(idt(sigma)*N())+
                                    rightfacet(idt(sigma)*N())) ) );

    // BC
    a31 += integrate(_range=markedfaces(mesh,{"Neumann","Tip"} ),
                     _expr=( trans(id(m))*(idt(sigma)*N()) ));

    auto a32 = form2( _trial=Wh, _test=Mh,_matrix=A,
                      _rowstart=2, _colstart=1);
    a32 += integrate(_range=internalfaces(mesh),
                     _expr=-tau_constant * trans(id(m)) * (leftfacet( pow(idv(H),M_tau_order)*idt(u) )+
                                                    rightfacet( pow(idv(H),M_tau_order)*idt(u) )));

    a32 += integrate(_range=markedfaces(mesh,{"Neumann","Tip"}),
                     _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tau_order)*idt(u) ) );

    auto a33 = form2(_trial=Mh, _test=Mh,_matrix=A,
                     _rowstart=2, _colstart=2);

    a33 += integrate(_range=internalfaces(mesh),
                     _expr=tau_constant * trans(idt(uhat)) * id(m) * ( leftface( pow(idv(H),M_tau_order) )+
                                                                 rightface( pow(idv(H),M_tau_order) )));

    a33 += integrate(_range=markedfaces(mesh,{"Neumann","Tip"}),
                     _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tau_order) ) );
    a33 += integrate(_range=markedfaces(mesh,{"Dirichlet"}),
                     _expr=trans(idt(uhat)) * id(m) );


}


} // Feel
