/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
   \file darcy.hpp
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
        ( "k1", po::value<double>()->default_value( 0.0 ), "k1 : k1=0 corresponds to homogeneous case" )
        ( "k2", po::value<double>()->default_value( 1.0 ), "k2" )
        ( "k_2D", po::value<std::string>()->default_value( "k1*(x-2)*x*(y-2)*y+k2:x:y:k1:k2" ), "k" )
        ( "u_exact_2D", po::value<std::string>()->default_value( "{-1/(2*Pi)*cos( Pi*x )*sin( Pi*y ), -1/(2*Pi)*sin( Pi*x )*cos( Pi*y )}:x:y" ), "u exact" )
        ( "p_exact_2D", po::value<std::string>()->default_value( "(1/(2*Pi*Pi))*sin( Pi*x )*sin( Pi*y ):x:y" ), "p exact" )
        ( "k_3D", po::value<std::string>()->default_value( "k1*(x-2)*x*(y-2)*y*(z-2)*z+k2:x:y:z:k1:k2" ), "k" )
        ( "u_exact_3D", po::value<std::string>()->default_value( "{-1/(2*Pi)*cos( Pi*x )*sin( Pi*y )*sin( Pi*z ),-1/(2*Pi)*sin( Pi*x )*cos( Pi*y )*sin( Pi*z ),-1/(2*Pi)*sin( Pi*x )*sin( Pi*y )*cos( Pi*z )}:x:y:z" ), "u exact" )
        ( "p_exact_3D", po::value<std::string>()->default_value( "(1/(2*Pi*Pi))*sin( Pi*x )*sin( Pi*y )*sin( Pi*z ):x:y:z" ), "p exact" )
        ( "d1", po::value<double>()->default_value( 0.0 ), "d1 (stabilization term)" )
        ( "d2", po::value<double>()->default_value( 0.0 ), "d2 (stabilization term)" )
        ( "nb_refine", po::value<int>()->default_value( 5 ), "nb_refine" )
        ;
    return testhdivoptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "darcy" ,
                     "darcy" ,
                     "0.1",
                     "convergence test for Darcy problem (using Hdiv conforming elts)",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

template<int Dim, int OrderU, int OrderP>
class Darcy
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

    //! the basis type of our approximation space
    typedef bases<RaviartThomas<OrderU> > RT_basis_type; //RT vectorial space
    typedef bases<Lagrange<OrderP,Vectorial> > lagrange_basis_v_type; //Lagrange vectorial space
    typedef bases<Lagrange<OrderP,Scalar> > lagrange_basis_s_type; //Lagrange scalar space
    typedef bases<Lagrange<0,Scalar> > lagrange_multiplier_basis_type; //P0 scalar space
    typedef bases< RaviartThomas<OrderU>, Lagrange<OrderP,Scalar> > basis_type; //For Darcy : (u,p) (\in H_div x L2)

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef FunctionSpace<mesh_type, lagrange_basis_s_type> lagrange_space_s_type;
    typedef FunctionSpace<mesh_type, lagrange_basis_v_type> lagrange_space_v_type;
    typedef FunctionSpace<mesh_type, RT_basis_type> RT_space_type;

    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<lagrange_space_s_type> lagrange_space_s_ptrtype;
    typedef boost::shared_ptr<lagrange_space_v_type> lagrange_space_v_ptrtype;
    typedef boost::shared_ptr<RT_space_type> RT_space_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    Darcy()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        exporter( Exporter<mesh_type>::New( this->vm() ) ),
        M_k1( doption("k1") ),
        M_k2( doption("k2") ),
        M_d1( doption("d1") ),
        M_d2( doption("d2") )
    {
        std::cout << "[TestHDiv]\n";
        this->changeRepository( boost::format( "/benchmark_darcy/%1%/h_%2%/%3%-%4%/" )
                                % this->about().appName()
                                % doption("hsize")
                                % M_k1
                                % M_k2
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

    double M_k1, M_k2, M_d1, M_d2;

}; //Darcy

template<int Dim, int OrderU, int OrderP>
void
Darcy<Dim, OrderU, OrderP>::convergence()
{
    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;

    auto k1 = cst(M_k1);
    auto k2 = cst(M_k2);

    auto k = expr(soption((boost::format( "k_%1%D" ) % Dim ).str()));
    auto lambda = cst(1.)/k;
    std::cout << "k : " << k << std::endl;

    // Exact solutions
    auto u_exact = expr<Dim,1>(soption((boost::format( "u_exact_%1%D" ) % Dim ).str()));
    auto divu_exact = div(u_exact);
    auto gradu_exact = grad(u_exact);
    // auto graddivu_exact = grad<Dim>(div(u_exact));
    auto f = divu_exact;
    std::cout << "u : " << u_exact << std::endl;
    std::cout << "divu : " << divu_exact << std::endl;
    std::cout << "gradu : " << gradu_exact << std::endl;

    auto p_exact = expr(soption((boost::format( "p_exact_%1%D" ) % Dim ).str()));
    auto gradp_exact = grad<Dim>(p_exact);
    std::cout << "p : " << p_exact << std::endl;
    std::cout << "gradp : " << gradp_exact << std::endl;


    // Coeff for stabilization terms
    auto d1 = cst(M_d1);
    auto d2 = cst(M_d2);

    std::ofstream cvg_p, cvg_u, cvg_divu, cvg_projL2, cvg_projHDIV;
    if( proc_rank == 0 )
        {
            // Open
            cvg_p.open( "convergence_p.dat", std::ios::out | std::ios::trunc);
            cvg_u.open( "convergence_u.dat", std::ios::out | std::ios::trunc);
            cvg_divu.open( "convergence_divu.dat", std::ios::out | std::ios::trunc);
            cvg_projL2.open( "convergence_projL2.dat", std::ios::out | std::ios::trunc);
            cvg_projHDIV.open( "convergence_projHDIV.dat", std::ios::out | std::ios::trunc);

            // Head
            cvg_u << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\t" << "h1err" << "\n";
            cvg_p << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\t" << "h1err" << "\n";
            cvg_divu << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\t" << "h1err" << "\n";
            cvg_projL2 << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\n";
            cvg_projHDIV << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\n";
        }

    double current_hsize = meshSize;
    for(int i=0; i<ioption("nb_refine"); i++)
        {
            mesh_ptrtype mesh;
            std::string mesh_name = (boost::format( "%1%-%2%D" ) % "hypercube" % Dim ).str();

            if( !fs::exists(mesh_name+".msh") )
                {
                    //std::cout << "createGMSHmesh" << std::endl;
                    LOG(INFO) << "[Darcy] Mesh has been created \n";
                    mesh = createGMSHMesh( _mesh=new mesh_type,
                                           _desc=domain( _name = mesh_name ,
                                                         _shape = "hypercube",
                                                         _usenames = true,
                                                         _dim = Dim,
                                                         _h = meshSize,
                                                         _xmin=0,_xmax=2,
                                                         _ymin=0,_ymax=2,
                                                         _zmin=0,_zmax=2 ) );
                }
            else
                {
                    //std::cout << "loadGMSHmesh" << std::endl;
                    LOG(INFO) << "[Darcy] Mesh has been loaded (refine level = " << i << ") \n";
                    mesh = loadGMSHMesh( _mesh=new mesh_type,
                                         _filename=mesh_name+".msh",
                                         _refine=i );
                }

            // ****** Dual-mixed formulation - Hdiv ******
            space_ptrtype Yh = space_type::New( mesh );
            RT_space_ptrtype RTh = RT_space_type::New( mesh ); //Dh<OrderU>( mesh );
            lagrange_space_s_ptrtype Xh = lagrange_space_s_type::New( mesh );
            lagrange_space_v_ptrtype Xhvec = lagrange_space_v_type::New( mesh );

            auto U_rt = Yh->element( "(u,p)" ); //trial
            auto V_rt = Yh->element( "(v,q)" ); //test

            auto u_rt = U_rt.template element<0>( "u" ); // velocity field
            auto v_rt = V_rt.template element<0>( "v" );
            auto p_rt = U_rt.template element<1>( "p" ); // potential field
            auto q_rt = V_rt.template element<1>( "q" );

            // Number of dofs associated with each space U and P
            auto nDofu = u_rt.functionSpace()->nDof();
            auto nDofp = p_rt.functionSpace()->nDof();

            auto F_rt = M_backend->newVector( Yh );
            auto darcyRT_rhs = form1( _test=Yh, _vector=F_rt );
            // fq
            darcyRT_rhs += integrate( _range=elements(mesh), _expr = f*id(q_rt) );
            // GLS(Hdiv) stabilization terms
            darcyRT_rhs += integrate( _range=elements(mesh), _expr = -d2*lambda*f*div(v_rt) );

            auto M_rt = M_backend->newMatrix( Yh, Yh );
            auto darcyRT = form2( _test=Yh, _trial=Yh, _matrix=M_rt);
            // Lambda u v
            darcyRT += integrate( _range=elements(mesh), _expr = -lambda*trans(idt(u_rt))*id(v_rt) );
            // p div(v)
            darcyRT += integrate( _range=elements(mesh), _expr = idt(p_rt)*div(v_rt) );
            // div(u) q
            darcyRT += integrate( _range=elements(mesh), _expr = divt(u_rt)*id(q_rt) );
            // GLS(Hdiv) stabilization terms
            darcyRT += integrate( _range=elements(mesh), _expr = d1*k*(trans(lambda*idt(u_rt)+trans(gradt(p_rt)))*(lambda*id(v_rt)+trans(grad(q_rt)))) );
            darcyRT += integrate( _range=elements(mesh), _expr = -d2*( lambda*divt(u_rt)*div(v_rt) ));

            // darcyRT += on( _range=boundaryfaces(mesh), _rhs=F_rt,  _element=u_rt,  _expr=u_exact );
            darcyRT += on( _range=boundaryfaces(mesh), _rhs=F_rt,  _element=p_rt,  _expr=p_exact );

            // Solve problem
            backend(_rebuild=true)->solve( _matrix=M_rt, _solution=U_rt, _rhs=F_rt );

            std::cout << "[Darcy] RT solve done" << std::endl;

            // ****** Compute error ******
            auto l2err_u = normL2( _range=elements(mesh), _expr=u_exact - idv(u_rt) );
            auto l2err_p = normL2( _range=elements(mesh), _expr=p_exact - idv(p_rt) );
            auto l2err_divu = normL2( _range=elements(mesh), _expr=divu_exact - divv(u_rt) );

            auto h1err_u = normH1( _range=elements(mesh), _expr=u_exact - idv(u_rt), _grad_expr=gradu_exact - gradv(u_rt) );
            auto h1err_p = normH1( _range=elements(mesh), _expr=p_exact - idv(p_rt), _grad_expr=trans(gradp_exact) - trans(gradv(p_rt)) );
            auto divu_rt = vf::project(Xh, elements(mesh), divv(u_rt) ); //TODO : Raviart-Thomas interpolant !
            // auto h1err_divu = normH1( _range=elements(mesh), _expr=divu_exact - divv(u_rt), _grad_expr=graddivu_exact - trans(gradv(divu_rt)) );

            // ****** Projection Operators (L2 - HDIV) : check proj( div u ) = proj( f ) ******
            // L2 projection
            auto l2_v = opProjection( _domainSpace=Xhvec, _imageSpace=Xhvec, _type=L2 ); //l2 vectorial proj
            auto l2_s = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=L2 ); //l2 scalar proj
            auto E_l2 = l2_v->project( _expr= u_exact );
            auto l2error_L2 = normL2( _range=elements(mesh), _expr=divv(E_l2) - f );

            auto hdiv = opProjection( _domainSpace=RTh, _imageSpace=RTh, _type=HDIV ); //hdiv proj (RT elts)
            auto E_hdiv = hdiv->project( _expr= (u_exact), _div_expr=f );
            auto l2error_HDIV = normL2( _range=elements(mesh), _expr=divv(E_hdiv) - f );

            if( proc_rank == 0 )
                {
                    std::cout << "[" << i << "]||u_exact - u||_L2 = " << l2err_u << std::endl;
                    std::cout << "[" << i << "]||p_exact - p||_L2 = " << l2err_p << std::endl;
                    std::cout << "[" << i << "]||divu_exact - divu||_L2 = " << l2err_divu << std::endl;

                    std::cout << "[" << i << "]||u_exact - u||_H1 = " << h1err_u << std::endl;
                    std::cout << "[" << i << "]||p_exact - p||_H1 = " << h1err_p << std::endl;
                    // std::cout << "[" << i << "]||divu_exact - divu||_H1 = " << h1err_divu << std::endl;

                    cvg_u << current_hsize << "\t" << nDofu << "\t" << l2err_u << "\t" << h1err_u << "\n";
                    cvg_p << current_hsize << "\t" << nDofp << "\t" << l2err_p << "\t" << h1err_p << "\n";
                    // cvg_divu << current_hsize << "\t" << nDofu << "\t" << l2err_divu << "\t" << h1err_divu << "\n";
                    cvg_projL2 << current_hsize << "\t" << Xhvec->nDof() << "\t" << l2error_L2 << "\n";
                    cvg_projHDIV << current_hsize << "\t" << RTh->nDof() << "\t" << l2error_HDIV << "\n";
                }

            // ****** Export results ******
            std::string exportName =  ( boost::format( "%1%-refine-%2%" ) % this->about().appName() % i ).str();
            std::string uName = ( boost::format( "velocity-refine-%1%" ) % i ).str();
            std::string u_exName = ( boost::format( "velocity-ex-refine-%1%" ) % i ).str();
            std::string pName = ( boost::format( "potential-refine-%1%" ) % i ).str();
            std::string p_exName = ( boost::format( "potential-ex-refine-%1%" ) % i ).str();

            auto u_ex_proj = vf::project(Xhvec, elements(mesh), u_exact );
            auto p_ex_proj = vf::project(p_rt.functionSpace(), elements(mesh), p_exact );

            export_ptrtype exporter_cvg( export_type::New( this->vm(), exportName ) );

            exporter_cvg->step( 0 )->setMesh( mesh );
            exporter_cvg->step( 0 )->add( uName, u_rt );
            exporter_cvg->step( 0 )->add( pName, p_rt );
            exporter_cvg->step( 0 )->add( u_exName, u_ex_proj );
            exporter_cvg->step( 0 )->add( p_exName, p_ex_proj );
            exporter_cvg->save();

            current_hsize /= 2.0;
        }

    if( proc_rank == 0 )
        {
            cvg_u.close();
            cvg_p.close();
            cvg_divu.close();
            cvg_projL2.close();
            cvg_projHDIV.close();
        }
}
