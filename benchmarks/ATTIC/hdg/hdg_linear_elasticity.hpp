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
        ( "hsize", po::value<double>()->default_value( 0.8 ), "mesh size" )
        ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ( "lambda", po::value<std::string>()->default_value( "1" ), "lambda" )
        ( "mu", po::value<std::string>()->default_value( "1" ), "mu" )
        ( "u_exact", po::value<std::string>()->default_value("empty"), "u exact" )
		( "f", po::value<std::string>()->default_value("empty"), "divergence of the stress tensor")
        ( "hface", po::value<int>()->default_value( 0 ), "hface" )
        ( "use_hypercube", po::value<bool>()->default_value( false ), "use hypercube or a given geometry" )
        ( "tau_constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "tau_order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
		( "use-sc", po::value<bool>()->default_value(true), "use static condensation")
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


template<int Dim, int OrderP, int OrderG=1>
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
    typedef Simplex<Dim,OrderG> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,OrderG,Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef std::shared_ptr<face_mesh_type> face_mesh_ptrtype;

	static const uint16_type expr_order = OrderP+4;

    using Vh_t =  Pdhms_type<mesh_type,OrderP>;
    using Vh_ptr_t =  Pdhms_ptrtype<mesh_type,OrderP>;
    using Wh_t =  Pdhv_type<mesh_type,OrderP>;
    using Wh_ptr_t =  Pdhv_ptrtype<mesh_type,OrderP>;
    using Mh_t =  Pdhv_type<face_mesh_type,OrderP>;
    using Mh_ptr_t =  Pdhv_ptrtype<face_mesh_type,OrderP>;
    using M0h_t =  Pdh_type<face_mesh_type,0>;
    using M0h_ptr_t =  Pdh_ptrtype<face_mesh_type,0>;
    using Ch_ptr_t = Pchv_ptrtype<face_mesh_type,0>;


    using product_space_type = ProductSpaces<Vh_ptr_t,Wh_ptr_t,Mh_ptr_t>;
    using product_space_ptrtype = std::shared_ptr<product_space_type>;

    using blockform2_type = BlockBilinearForm<ProductSpaces< Vh_ptr_t, Wh_ptr_t, Mh_ptr_t > &>;
    using blockform2_ptrtype = std::shared_ptr<blockform2_type>;
    using blockform1_type = BlockLinearForm<ProductSpaces< Vh_ptr_t, Wh_ptr_t, Mh_ptr_t > &>;
    using blockform1_ptrtype = std::shared_ptr<blockform1_type>;

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
    void run();


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

template<int Dim, int OrderP, int OrderG>
void
Hdg<Dim, OrderP, OrderG>::run()
{
    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;

	double sc_param = boption("use-sc") ? 0.5 : 1.0;
    
    auto lambda = expr(soption("lambda"));
    auto mu     = expr(soption("mu"));

    // Exact solutions
    auto u_exact = expr<Dim,1,expr_order>(soption("u_exact"));
	if (soption("u_exact") == "empty")
	{
		if (Dim == 2)
			u_exact = expr<Dim,1,expr_order>("{(1/(2*Pi*Pi))*sin(Pi*x)*cos(Pi*y),(1/(2*Pi*Pi))*cos(Pi*x)*sin(Pi*y)}:x:y");
		else
			u_exact = expr<Dim,1,expr_order>("{ cos(Pi*x)*cos(Pi*y)*cos(Pi*z), cos(Pi*y)*sin(Pi*x)*sin(Pi*z), cos(Pi*x)*cos(Pi*z)*sin(Pi*y)  }:x:y:z");
	}
    auto gradu_exact = grad<Dim,expr_order>(u_exact);
    auto eps_exact   = cst(0.5) * ( gradu_exact + trans(gradu_exact) );
    auto sigma_exact = lambda * trace(eps_exact) * eye<Dim>() + cst(2.) * mu * eps_exact;
    auto f = expr<Dim,1,expr_order>(soption("f"));
	if (soption("f") == "empty")
	{
		if (Dim == 2)
			f = expr<Dim,1,expr_order>("{-3.0*sin(pi*x)*cos(pi*y),-3.0*sin(pi*y)*cos(pi*x)}:x:y");
		else
			f = expr<Dim,1,expr_order>("{ 2*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z) - 2*pi^2*cos(pi*x)*sin(pi*y)*sin(pi*z) - 5*pi^2*cos(pi*x)*cos(pi*y)*cos(pi*z), 2*pi^2*cos(pi*z)*sin(pi*x)*sin(pi*y) - 5*pi^2*cos(pi*y)*sin(pi*x)*sin(pi*z) - 2*pi^2*cos(pi*x)*cos(pi*y)*sin(pi*z), 2*pi^2*cos(pi*y)*sin(pi*x)*sin(pi*z) - 5*pi^2*cos(pi*x)*cos(pi*z)*sin(pi*y) - 2*pi^2*cos(pi*z)*sin(pi*x)*sin(pi*y)  }:x:y:z");
	}

    cout << "lambda : " << lambda      << std::endl;
    cout << "mu     : " << mu          << std::endl;
    cout << "u      : " << u_exact     << std::endl;
    //    cout << "sigma  : " << sigma_exact << std::endl;
    cout << "f      : " << f           << std::endl;


    // Coeff for stabilization terms
    auto tau_constant = cst(M_tau_constant);


	int status = 0;
    double current_hsize = meshSize;

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
            mesh->addMarkerName( "Dirichlet",( Dim==2 )?1:19, (Dim==2)?1:2);
            mesh->addMarkerName( "Neumann",( Dim==2)?3:27, (Dim==2)?1:2);
        }
        else
            mesh = loadMesh( new mesh_type);
    }
	else
    {
        //std::cout << "loadGMSHmesh" << std::endl;
        LOG(INFO) << "[Hdg] Mesh has been loaded \n";
        mesh = loadGMSHMesh( _mesh=new mesh_type,
                             _filename=mesh_name+".msh",
                             _refine=0,
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


	// auto lambda = expr(soption("lambda"));
	// auto mu     = expr(soption("mu"));
	auto c1     = cst(0.5)/mu;
 	auto c2     = -lambda/(cst(2.) * mu * (cst(Dim)*lambda + cst(2.)*mu));



	tic();
    solve::strategy s = boption("use-sc")?solve::strategy::static_condensation:solve::strategy::monolithic;
	




#if 1	
	auto ps = product( Vh, Wh, Mh );
	auto a = blockform2( ps , s, backend() );
	auto rhs = blockform1( ps , s, backend() );
#else
    Ch_ptr_t Ch = Pchv<0>(face_mesh, true);
	auto ibcSpaces = std::make_shared<ProductSpace<Ch_ptr_t,true> >(0, Ch);
	auto ps = product2( ibcSpaces, Vh, Wh, Mh );

	auto A = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(ps), backend() );
	auto F = makeSharedVectorCondensed<value_type>(s, blockVector(ps), backend(), false);
    auto a = blockform2(ps, A); 
    auto rhs = blockform1(ps, F);
#endif

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

    rhs(1_c) += integrate(_range=elements(mesh), _expr=trans(f)*id(w));
    // rhs(1_c) += integrate(_range=elements(mesh), _expr= divv(sigma_exact)*id(w));

    cout << "rhs2 works fine" << std::endl;

    // in convergence test Neumann condition is given from the displacement and
    // constitutive law
    rhs(2_c) += integrate(_range=markedfaces(mesh,"Neumann"), _expr=trans(id(m))*(sigma_exact*N()));

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
                              _expr=sc_param*tau_constant * trans(idt(uhat)) * id(m) * ( leftface( pow(idv(H),M_tau_order) )+
                                                                                          rightface( pow(idv(H),M_tau_order) )));

    a( 2_c, 2_c) += integrate(_range=markedfaces(mesh,"Dirichlet"),
                              _expr=trans(idt(uhat)) * id(m) );


    a( 2_c, 0_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                              _expr=( trans(id(m))*(idt(sigma)*N()) ));

    a( 2_c, 1_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                              _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tau_order)*idt(u) ) );

    a( 2_c, 2_c) += integrate(_range=markedfaces(mesh,"Neumann"),
                              _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tau_order) ) );

    toc("matrices",true);

    tic();
    auto U = ps.element();
    auto Ue = ps.element();
    a.solve( _solution=U, _rhs=rhs, _rebuild=true, _condense=boption("use-sc"));

    toc("solve",true);
    cout << "[Hdg] solve done" << std::endl;

    auto sigmap = U(0_c);
    auto up = U(1_c);
    auto uhatp = U(2_c);

    Ue(0_c).on( _range=elements(mesh), _expr=sigma_exact );
    Ue(1_c).on( _range=elements(mesh), _expr=u_exact );

#if 0
        Feel::cout << "sigma exact: \t" << Ue(0_c) << std::endl;
        Feel::cout << "sigma: \t" << sigmap << std::endl;
        Feel::cout << "u exact: \t" << Ue(1_c) << std::endl;
        Feel::cout << "u: \t" << up << std::endl;
        Feel::cout << "uhat: \t" << uhatp << std::endl;
#endif


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
	auto h1err_u 	 = normH1(_range=elements(mesh), _expr=idv(up)- u_exact, _grad_expr=gradv(up)-grad(u_exact) );
   	toc("error");

	if ( !checker().check() )
	{
    	cout << "||sigma_exact - sigma||_L2 = " << l2err_sigma << std::endl;
    	cout << "||u_exact - u||_L2 = " << l2err_u << std::endl;
    	cout << "||u_exact - u||_H1 = " << h1err_u << std::endl;
	}

	// CHECKER
	if ( checker().check() )
	{

		// compute l2 and h1 norm of u-u_h where u=solution
		auto norms = [=]( std::string const& u_exact ) ->std::map<std::string,double>
		{
			tic();
			double l2 = l2err_u; 
			toc("L2 error norm");
			
			tic();
			double h1 = h1err_u; 
			toc("H1 error norm");
			
			return { { "L2", l2 } , {  "H1", h1 } };
		};
		
		status = checker().runOnce( norms, rate::hp( mesh->hMax(), Wh->fe()->order() ) );
		
	}

    tic();
    std::string exportName =  ( boost::format( "%1%" ) % this->about().appName() ).str();
    std::string sigmaName = "stress";
    std::string sigma_exName = "stress-ex";
    std::string uName = "displacement";
    std::string u_exName = "displacement-ex";

    v.on( _range=elements(mesh), _expr=sigma_exact , _quad=_Q<expr_order>());
    w.on( _range=elements(mesh), _expr=u_exact , _quad=_Q<expr_order>());
    export_ptrtype exporter_cvg( export_type::New( exportName ) );

    exporter_cvg->setMesh( mesh );
    exporter_cvg->add( sigmaName, sigmap );
    exporter_cvg->add( uName, up );
    exporter_cvg->add( sigma_exName, v );
    exporter_cvg->add( u_exName, w );
    exporter_cvg->save();

    toc("export");

	// return !status;
}



} // Feel
