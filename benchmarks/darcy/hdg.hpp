/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Author(s): Cecile Daversin  <cecile.daversin@lncmi.cnrs.fr>
       Date: 2011-12-07

  Copyright (C) 2011 UJF
  Copyright (C) 2011 CNRS

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
/**
   \file hdg.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2014-03-05
 */
#include <feel/feel.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

using namespace Feel;

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
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

// begin{dp}. OrderU dropped
// end{dp}
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
    // begin{dp}. The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,1> convex_type_multiplier;
    typedef Mesh<convex_type_multiplier> mesh_type_multiplier;
    typedef boost::shared_ptr<mesh_type_multiplier> mesh_ptrtype_multiplier;
    // end{dp}

    //! the basis type of our approximation space
    //    typedef bases<RaviartThomas<OrderU> > RT_basis_type; //RT vectorial space
    typedef bases<Lagrange<OrderP,Vectorial> > lagrange_basis_v_type; //Lagrange vectorial space
    typedef bases<Lagrange<OrderP,Scalar> > lagrange_basis_s_type; //Lagrange scalar space
    typedef bases<Lagrange<OrderP,Scalar> > lagrange_multiplier_basis_type; //PK scalar space
    // begin{dp}. How to define this global basis?
    // The problem is that P^k on tetrahedra is not the same as P^k on a faces.
    // Should we treat Vh, Wh, Mh separately?
    typedef bases<Lagrange<OrderP,Vectorial>, Lagrange<OrderP,Scalar>, Lagrange<OrderP,Scalar> > basis_type; //For Hdg : (u,p,phat) (\in [L2]^n x L2)
    // end{dp}

    //! the approximation function space type
    // begin{dp}. Do FunctionSpace objects support two meshes (Elements/Faces)?
    // If not, we definitely need to work with Vh, Wh, Mh separately (or Vh x Wh, Mh, at least)
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    // end{dp}
    typedef FunctionSpace<mesh_type, lagrange_basis_s_type> lagrange_space_s_type;
    typedef FunctionSpace<mesh_type, lagrange_basis_v_type> lagrange_space_v_type;
    typedef FunctionSpace<mesh_type, RT_basis_type> RT_space_type;
    // begin{dp}. To deal with Mh separately
    typedef FunctionSpace<mesh_ptrtype_multiplier, lagrange_multiplier_basis_type> lagrange_space_multiplier_type;
    // end{dp}

    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<lagrange_space_s_type> lagrange_space_s_ptrtype;
    typedef boost::shared_ptr<lagrange_space_v_type> lagrange_space_v_ptrtype;
    typedef boost::shared_ptr<RT_space_type> RT_space_ptrtype;
    // begin{dp}. To deal with Mh separately
    typedef boost::shared_ptr<lagrange_space_multiplier_type> lagrange_space_multiplier_ptrtype;
    // end{dp}

    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;

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
        // begin{dp}
        M_tau_constant( doption("tau_constant") ),
        M_tau_order( ioption("tau_order") )
        // end{dp}
    {
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

    auto k = expr(soption("k"));
    auto lambda = cst(1.)/k;

    // Exact solutions
    auto p_exact = expr(soption("p_exact"));
    auto gradp_exact = grad<Dim>(p_exact);
    auto u_exact = -k*trans(gradp_exact);
    auto f = -k*laplacian(p_exact);

    if( Environment::isMasterRank()){
        std::cout << "k : " << k /*<< "\tlambda : " << lambda*/ << std::endl;
        std::cout << "p : " << p_exact << std::endl;
        std::cout << "gradp : " << gradp_exact << std::endl;
        // std::cout << "u : " << u_exact << std::endl;
        // std::cout << "divu : " << f << std::endl;
    }


    // Coeff for stabilization terms
    // auto d1 = cst(M_d1);
    // auto d2 = cst(M_d2);
    // auto d3 = cst(M_d3);
    auto tau_constant = cst(M_tau_constant);

    std::ofstream cvg_p, cvg_u, cvg_divu, cvg_projL2, cvg_projHDIV, cvg_inter;
    if( proc_rank == 0 )
        {
            // Open
            cvg_p.open( "convergence_p.dat", std::ios::out | std::ios::trunc);
            cvg_u.open( "convergence_u.dat", std::ios::out | std::ios::trunc);
            cvg_divu.open( "convergence_divu.dat", std::ios::out | std::ios::trunc);
            cvg_projL2.open( "convergence_projL2.dat", std::ios::out | std::ios::trunc);
            cvg_projHDIV.open( "convergence_projHDIV.dat", std::ios::out | std::ios::trunc);
            cvg_inter.open( "convergence_inter.dat", std::ios::out | std::ios::trunc);

            // Head
            cvg_u << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\n";
            cvg_p << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\t" << "h1err" << "\n";
            cvg_divu << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\n";
            cvg_projL2 << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\n";
            cvg_projHDIV << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\n";
            cvg_inter << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\n";
        }

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
                        mesh = loadMesh( new mesh_type );
                }
            else
                {
                    //std::cout << "loadGMSHmesh" << std::endl;
                    LOG(INFO) << "[Hdg] Mesh has been loaded (refine level = " << i << ") \n";
                    mesh = loadGMSHMesh( _mesh=new mesh_type,
                                         _filename=mesh_name+".msh",
                                         _refine=i );
                }
            toc("mesh");

            // ****** Hybrid-mixed formulation ******
            // We try to treat Vh, Wh, and Mh separately
            tic();
            // Stuff from darcy.hpp:
            //
            // space_ptrtype Yh = space_type::New( mesh );
            // RT_space_ptrtype RTh = RT_space_type::New( mesh ); //Dh<OrderU>( mesh );
            // lagrange_space_s_ptrtype Xh = lagrange_space_s_type::New( mesh );
            // lagrange_space_v_ptrtype Xhvec = lagrange_space_v_type::New( mesh );

            // begin{dp}.
            lagrange_space_v_ptrtype Vh = lagrange_space_v_type::New( mesh );
            lagrange_space_s_ptrtype Wh = lagrange_space_s_type::New( mesh );
            // >>>>>>>>>>> CRITICAL QUESTION: is faces(mesh) a mesh in R^{n-1}????? <<<<<<<<<<<
            lagrange_space_multiplier_ptrtype Mh = lagrange_space_multiplier_type::New( faces(mesh) );
            // end{dp}

            toc("spaces");
            if( Environment::isMasterRank() )
                // begin{dp}
                std::cout << "Vh<" << OrderP << "> : " << Vh->nDof() << std::endl
                          << "Wh<" << OrderP << "> : " << Wh->nDof() << std::endl
                          << "Mh<" << OrderP << "> : " << Mh->nDof() << std::endl;
                // end{dp}

            // auto U_rt = Yh->element( "(u,p)" ); //trial
            // auto V_rt = Yh->element( "(v,q)" ); //test

            // auto u_rt = U_rt.template element<0>( "u" ); // velocity field
            // auto v_rt = V_rt.template element<0>( "v" );
            // auto p_rt = U_rt.template element<1>( "p" ); // potential field
            // auto q_rt = V_rt.template element<1>( "q" );

            // begin{dp}
            auto u_rt = Vh->element( "u" );
            auto v_rt = Vh->element( "v" );
            auto p_rt = Wh->element( "p" );
            auto q_rt = Wh->element( "q" );
            auto phat_rt = Mh->element( "phat" );
            auto lambda_rt = Mh->element( "lambda" );
            // end{dp}

            // Number of dofs associated with each space U and P
            auto nDofu = u_rt.functionSpace()->nDof();
            auto nDofp = p_rt.functionSpace()->nDof();
            // begin{dp}. This should be different from nDofp
            auto nDofphat = phat_rt.functionSpace()->nDof();
            // end{dp}

            tic();
            // begin{dp}. Question: Is there a way to create the cartesian product
            // Vh x Wh x M at this point? Or do we need to create it from
            // the beginning of the program (but remember initial remarks)?
            //
            // Stuff from darcy.hpp:
            //
            // auto F_rt = M_backend->newVector( Yh );
            // auto hdgRT_rhs = form1( _test=Yh, _vector=F_rt );
            // // fq
            // hdgRT_rhs += integrate( _range=elements(mesh), _expr = -f*id(q_rt) );
            // // GLS(Hdiv) stabilization terms
            // hdgRT_rhs += integrate( _range=elements(mesh)
            //                           , _expr = d2*lambda*f*div(v_rt) );
            //
            // end{dp}


            a11 = form2( _trial=Vh, _test=Vh );
            a11 += integrate(_range=elements(mesh),_expr=(trans(K*idt(u))*id(v)) );

            // begin{dp}
            // Watch out the negative sign!
            a12 = form2( _trial=Wh, _test=Vh );
            a12 += integrate(_range=elements(mesh),_expr=-(idt(p)*div(v)));

            a13 = form2( _trial=Mh, _test=Vh );
            a13 += integrate(_range=internalfaces(mesh),
                            _expr=( idt(phat)*leftface(trans(id(v))*N())+
                                    idt(phat)*rightface(trans(id(v))*N())) );
            a13 += integrate(_range=boundaryfaces(mesh),
                             _expr=(idt(phat)*trans(id(v))*N()));

            a21 = form2( _trial=Vh, _test=Wh );
            a21 += integrate(_range=elements(mesh),_expr=(-grad(w)*idt(u)));
            a21 += integrate(_range=internalfaces(mesh),
                             _expr=( leftface(id(w)*trans(idt(u))*N())+
                                     rightface(id(w)*trans(idt(u))*N())) );
            a21 += integrate(_range=boundaryfaces(mesh),
                             _expr=(id(w)*trans(idt(u))*N()));

            a22 = form2( _trial=Wh, _test=Wh );
            a22 += integrate(_range=internalfaces(mesh),
                            _expr=tau_constant * ( leftface( pow(h(),M_tau_order)*idt(p)*id(w) )+
                                                   rightface( pow(h(),M_tau_order)*idt(p)*id(w) )));
            a22 += integrate(_range=boundaryfaces(mesh),
                             _expr=(tau_constant * pow(h(),M_tau_order)*id(w)*idt(p)));

            a23 = form2( _trial=Mh, _test=Wh );
            a23 += integrate(_range=internalfaces(mesh),
                            _expr=-tau_constant * idt(phat) * ( leftface( pow(h(),M_tau_order)*id(w) )+
                                                   rightface( pow(h(),M_tau_order)*id(w) )));
            a23 += integrate(_range=boundaryfaces(mesh),
                             _expr=-tau_constant * idt(phat) * pow(h(),M_tau_order)*id(w) );

            a31 = form2( _trial=Vh, _test=Mh );
            a31 += integrate(_range=internalfaces(mesh),
                             _expr=( id(l)*leftface(trans(idt(u))*N())+

            a32 = form2( _trial=Wh, _test=Mh );
            a32 += integrate(_range=internalfaces(mesh),
                             _expr=tau_constant * id(l) * ( leftface( pow(h(),M_tau_order)*idt(p) )+
                                                   rightface( pow(h(),M_tau_order)*idt(p) )));

            a33 = form2(_trial=Mh, _test=Mh);
            a23 += integrate(_range=internalfaces(mesh),
                             _expr=-tau_constant * idt(phat) * id(l) * ( leftface( pow(h(),M_tau_order) )+
                                                   rightface( pow(h(),M_tau_order) )));


            // Building the RHS
            //
            // This is only a part of the RHS - how to build the whole RHS? Is it right to
            // imagine we moved it to the left? SKIPPING boundary conditions for the moment.
            // How to identify Dirichlet/Neumann boundaries?
            auto rhs2 = form1( _test=Wh );
            rhs2 += integrate(_range=elements(mesh),
                              _expr=-f*id(w));

            // auto M_rt = M_backend->newMatrix( Yh, Yh );
            // auto hdgRT = form2( _test=Yh, _trial=Yh, _matrix=M_rt);
            // // Lambda u v
            // hdgRT += integrate( _range=elements(mesh),
            //                       _expr = lambda*trans(idt(u_rt))*id(v_rt) );
            // // p div(v)
            // hdgRT += integrate( _range=elements(mesh),
            //                       _expr = -idt(p_rt)*div(v_rt) );
            // // div(u) q
            // hdgRT += integrate( _range=elements(mesh),
            //                       _expr = -divt(u_rt)*id(q_rt) );
            // // GLS(Hdiv) stabilization terms
            // hdgRT += integrate( _range=elements(mesh),
            //                       _expr = d1*k*(trans(lambda*idt(u_rt)+trans(gradt(p_rt)))*(lambda*id(v_rt)+trans(grad(q_rt)))) );
            // hdgRT += integrate( _range=elements(mesh),
            //                       _expr = d2*( lambda*divt(u_rt)*div(v_rt) ));
            // // only homogenous k
            // hdgRT += integrate( _range=elements(mesh),
            //                       _expr = d3*( lambda*trans(curlt(u_rt))*curl(v_rt) ));

            // hdgRT += on( _range=boundaryfaces(mesh), _rhs=F_rt,  _element=u_rt,  _expr=u_exact );

            // end{dprada}
            toc("matrices");

            tic();
            // begin{dp}. Solve problem: how to build the global entities?
            //            backend(_rebuild=true)->solve( _matrix=M_rt, _solution=U_rt, _rhs=F_rt );
            // end{dp}
            toc("solve");
            if( Environment::isMasterRank())
                std::cout << "[Hdg] solve done" << std::endl;

            // ****** Compute error ******
            tic();

            auto mean_p_exact = mean( elements(mesh), p_exact )(0,0);
            auto mean_p = mean( elements(mesh), idv(p_rt) )(0,0);

            auto l2err_u = normL2( _range=elements(mesh), _expr=u_exact - idv(u_rt) );
            auto l2err_p = normL2( elements(mesh),
                                   (p_exact - cst(mean_p_exact)) - (idv(p_rt) - cst(mean_p)) ); // moyenne nulle
            auto l2err_divu = normL2( _range=elements(mesh), _expr=f - divv(u_rt) );
            auto h1err_p = normH1( elements(mesh),
                                   (p_exact - cst(mean_p_exact)) - (idv(p_rt) - cst(mean_p)),
                                   _grad_expr=trans(gradp_exact) - trans(gradv(p_rt)) );
            auto divu_h = vf::project(Wh, elements(mesh), divv(u_rt) ); 

            // ****** Projection Operators (L2 - HDIV) : check proj( div u ) = proj( f ) ******
            // L2 projection
            auto l2_v = opProjection( _domainSpace=Vh, _imageSpace=Vh, _type=L2 ); //l2 vectorial proj
            auto l2_s = opProjection( _domainSpace=Wh, _imageSpace=Wh, _type=L2 ); //l2 scalar proj
            auto E_l2 = l2_v->project( _expr= u_exact );
            auto l2error_L2 = normL2( _range=elements(mesh), _expr=divv(E_l2) - f );

            // auto hdiv = opProjection( _domainSpace=RTh, _imageSpace=RTh, _type=HDIV ); //hdiv proj (RT elts)
            // auto E_hdiv = hdiv->project( _expr= (u_exact), _div_expr=f );
            // auto l2error_HDIV = normL2( _range=elements(mesh), _expr=divv(E_hdiv) - f );

            // interpolant
            v_rt.on( elements(mesh), u_exact);
            auto l2error_inter = normL2( elements(mesh), idv(v_rt) - u_exact);

            toc("error");

            if( proc_rank == 0 )
                {
                    std::cout << "[" << i << "]||u_exact - u||_L2 = " << l2err_u << std::endl;
                    std::cout << "[" << i << "]||p_exact - p||_L2 = " << l2err_p << std::endl;
                    std::cout << "[" << i << "]||f - divu||_L2 = " << l2err_divu << std::endl;
                    std::cout << "[" << i << "]||p_exact - p||_H1 = " << h1err_p << std::endl;
                    std::cout << "[" << i << "]||I_h(p_exact) - p_exact||_L2 = " << l2error_inter << std::endl;

                    cvg_u << current_hsize << "\t" << nDofu << "\t" << l2err_u << "\n";
                    cvg_p << current_hsize << "\t" << nDofp << "\t" << l2err_p << "\t" << h1err_p << "\n";
                    cvg_divu << current_hsize << "\t" << nDofu << "\t" << l2err_divu << "\n";
                    cvg_projL2 << current_hsize << "\t" << Vh->nDof() << "\t" << l2error_L2 << "\n";
                    // cvg_projHDIV << current_hsize << "\t" << RTh->nDof() << "\t" << l2error_HDIV << "\n";
                    // cvg_inter << current_hsize << "\t" << RTh->nDof() << "\t" << l2error_inter << "\n";
                }

            // ****** Export results ******
            tic();
            std::string exportName =  ( boost::format( "%1%-refine-%2%" ) % this->about().appName() % i ).str();
            std::string uName = ( boost::format( "velocity-refine-%1%" ) % i ).str();
            std::string u_exName = ( boost::format( "velocity-ex-refine-%1%" ) % i ).str();
            std::string pName = ( boost::format( "potential-refine-%1%" ) % i ).str();
            std::string p_exName = ( boost::format( "potential-ex-refine-%1%" ) % i ).str();

            auto u_ex_proj = vf::project(Xhvec, elements(mesh), u_exact );
            auto p_ex_proj = vf::project(p_rt.functionSpace(), elements(mesh), p_exact );

            export_ptrtype exporter_cvg( export_type::New( exportName ) );

            exporter_cvg->step( i )->setMesh( mesh );
            exporter_cvg->step( i )->add( uName, u_rt );
            exporter_cvg->step( i )->add( pName, p_rt );
            exporter_cvg->step( i )->add( u_exName, u_ex_proj );
            exporter_cvg->step( i )->add( p_exName, p_ex_proj );
            exporter_cvg->save();

            current_hsize /= 2.0;
            toc("export");
        }

    if( proc_rank == 0 )
        {
            cvg_u.close();
            cvg_p.close();
            cvg_divu.close();
            cvg_projL2.close();
            cvg_projHDIV.close();
            cvg_inter.close();
        }
}
