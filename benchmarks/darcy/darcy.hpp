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
        ( "d1", po::value<double>()->default_value( 0.0 ), "d1 (stabilization term)" )
        ( "d2", po::value<double>()->default_value( 0.0 ), "d2 (stabilization term)" )
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
    void convergence(int nb_refine);

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
Darcy<Dim, OrderU, OrderP>::convergence(int nb_refine)
{
    int proc_rank = Environment::worldComm().globalRank();
    auto Pi = M_PI;

    auto k1 = cst(M_k1);
    auto k2 = cst(M_k2);

#if defined(FEELPP_2D)
    auto k = k1*(Px()-2)*Px()*(Py()-2)*Py() + k2;
    auto lambda = cst(1.)/k;

    auto f = k*sin( Pi*Px() )*sin( Pi*Py() )
        -(k1/Pi)*( (Px()-1)*(Py()-2)*Py()*cos( Pi*Px() )*sin( Pi*Py() ) + (Px()-2)*Px()*(Py()-1)*sin( Pi*Px() )*cos( Pi*Py() ) );

    // Exact solution
    auto dkdx = 2*k1*(Px()-1)*Py()*(Py()-2);
    auto dkdy = 2*k1*(Py()-1)*Px()*(Px()-2);
    auto coeff_uex = -k/(2*Pi); //c
    auto dcdx = -k1*(Px()-1)*(Py()-2)*Py()/Pi;
    auto dcdy = -k1*(Py()-1)*(Px()-2)*Px()/Pi;

    auto u_exact = coeff_uex*vec(cos( Pi*Px() )*sin( Pi*Py() ), sin( Pi*Px() )*cos( Pi*Py() ));
    auto p_exact = (1/(2*Pi*Pi))*sin( Pi*Px() )*sin( Pi*Py() );
    auto divu_exact = f;

    auto gradu_exact = mat<2,2>( dcdx*cos(Pi*Px())*sin(Pi*Py()) + (k/2)*sin(Pi*Px())*sin(Pi*Py()),
                                 dcdy*cos(Pi*Px())*sin(Pi*Py()) - (k/2)*cos(Pi*Px())*cos(Pi*Py()),
                                 dcdx*sin(Pi*Px())*cos(Pi*Py()) - (k/2)*cos(Pi*Px())*cos(Pi*Py()),
                                 dcdy*sin(Pi*Px())*cos(Pi*Py()) + (k/2)*sin(Pi*Px())*sin(Pi*Py()) );
    auto gradp_exact = -lambda*u_exact;
    auto graddivu_exact = vec( dkdx*sin(Pi*Px())*sin(Pi*Py()) + k*Pi*cos(Pi*Px())*sin(Pi*Py())
                               -(k1/Pi)*( Py()*(Py()-2)*sin(Pi*Py())*( cos(Pi*Px())-Pi*(Px()-1)*sin(Pi*Px()) )
                                          + (Py()-1)*cos(Pi*Py())*( 2*(Px()-1)*sin(Pi*Px()) + Pi*Px()*(Px()-2)*cos(Pi*Px()) )),
                               dkdy*sin(Pi*Px())*sin(Pi*Py()) + k*Pi*sin(Pi*Px())*cos(Pi*Py())
                               -(k1/Pi)*( (Px()-1)*cos(Pi*Px())*( 2*(Py()-1)*sin(Pi*Py()) + Pi*Py()*(Py()-2)*cos(Pi*Py()) )
                                          + Px()*(Px()-2)*sin(Pi*Px())*( cos(Pi*Py()) - Pi*(Py()-1)*sin(Pi*Py()) )) );

#elif defined(FEELPP_3D)
    auto k = k1*(Px()-2)*Px()*(Py()-2)*Py()*(Pz()-2)*Pz() + k2;
    auto lambda = cst(1.)/k;

    auto f = (3/2)*k*sin( Pi*Px() )*sin( Pi*Py() )*sin( Pi*Pz() )
        -(k1/Pi)*( (Px()-1)*(Py()-2)*Py()*(Pz()-2)*Pz()*cos( Pi*Px() )*sin( Pi*Py() )*sin( Pi*Pz() )
                   + (Px()-2)*Px()*(Py()-1)*(Pz()-2)*Pz()*sin( Pi*Px() )*cos( Pi*Py() )*sin( Pi*Pz() )
                   + (Px()-2)*Px()*(Py()-2)*Py()*(Pz()-1)*sin( Pi*Px() )*sin( Pi*Py() )*cos( Pi*Pz() ) );

    // Exact solution
    auto dkdx = 2*k1*(Px()-1)*Py()*(Py()-2)*Pz()*(Pz()-2);
    auto dkdy = 2*k1*(Py()-1)*Px()*(Px()-2)*Pz()*(Pz()-2);
    auto dkdz = 2*k1*(Pz()-1)*Px()*(Px()-2)*Py()*(Py()-2);
    auto coeff_uex = -k/(2*Pi); //c
    auto dcdx = -k1*(Px()-1)*Py()*(Py()-2)*Pz()*(Pz()-2)/Pi;
    auto dcdy = -k1*(Py()-1)*Px()*(Px()-2)*Pz()*(Pz()-2)/Pi;
    auto dcdz = -k1*(Pz()-1)*Px()*(Px()-2)*Py()*(Py()-2)/Pi;

    auto u_exact = coeff_uex*vec(cos( Pi*Px() )*sin( Pi*Py() )*sin( Pi*Pz() ),
                                 sin( Pi*Px() )*cos( Pi*Py() )*sin( Pi*Pz() ),
                                 sin( Pi*Px() )*sin( Pi*Py() )*cos( Pi*Pz() ));
    auto p_exact = (1/(2*Pi*Pi))*sin( Pi*Px() )*sin( Pi*Py() )*sin( Pi*Pz() );
    auto divu_exact = f;

    auto gradu_exact = mat<3,3>( dcdx*cos(Pi*Px())*sin(Pi*Py())*sin(Pi*Pz()) + (k/2)*sin(Pi*Px())*sin(Pi*Py())*sin(Pi*Pz()),
                                 dcdy*cos(Pi*Px())*sin(Pi*Py())*sin(Pi*Pz()) - (k/2)*cos(Pi*Px())*cos(Pi*Py())*sin(Pi*Pz()),
                                 dcdz*cos(Pi*Px())*sin(Pi*Py())*sin(Pi*Pz()) - (k/2)*cos(Pi*Px())*sin(Pi*Py())*cos(Pi*Pz()),
                                 dcdx*sin(Pi*Px())*cos(Pi*Py())*sin(Pi*Pz()) - (k/2)*cos(Pi*Px())*cos(Pi*Py())*sin(Pi*Pz()),
                                 dcdy*sin(Pi*Px())*cos(Pi*Py())*sin(Pi*Pz()) + (k/2)*sin(Pi*Px())*sin(Pi*Py())*sin(Pi*Pz()),
                                 dcdz*sin(Pi*Px())*cos(Pi*Py())*sin(Pi*Pz()) - (k/2)*sin(Pi*Px())*cos(Pi*Py())*cos(Pi*Pz()),
                                 dcdx*sin(Pi*Px())*sin(Pi*Py())*cos(Pi*Pz()) - (k/2)*cos(Pi*Px())*sin(Pi*Py())*cos(Pi*Pz()),
                                 dcdy*sin(Pi*Px())*sin(Pi*Py())*cos(Pi*Pz()) - (k/2)*sin(Pi*Px())*cos(Pi*Py())*cos(Pi*Pz()),
                                 dcdz*sin(Pi*Px())*sin(Pi*Py())*cos(Pi*Pz()) + (k/2)*sin(Pi*Px())*sin(Pi*Py())*sin(Pi*Pz()) );
    auto gradp_exact = -lambda*u_exact;
    auto graddivu_exact = vec( (3/2)*( dkdx*sin(Pi*Px())*sin(Pi*Py())*sin(Pi*Pz()) + k*Pi*cos(Pi*Px())*sin(Pi*Py())*sin(Pi*Pz()) )
                               -(k1/Pi)*( Py()*(Py()-2)*(Pz()-2)*Pz()*sin(Pi*Py())*sin(Pi*Pz())*( cos(Pi*Px())-Pi*(Px()-1)*sin(Pi*Px()) )
                                          + (Py()-1)*(Pz()-2)*Pz()*cos(Pi*Py())*sin(Pi*Pz())*( 2*(Px()-1)*sin(Pi*Px()) + Pi*Px()*(Px()-2)*cos(Pi*Px()) )
                                          + (Py()-2)*Py()*(Pz()-1)*sin(Pi*Py())*cos(Pi*Pz())*( 2*(Px()-1)*sin(Pi*Px()) + Pi*Px()*(Px()-2)*cos(Pi*Px()) )),
                               (3/2)*( dkdy*sin(Pi*Px())*sin(Pi*Py())*sin(Pi*Pz()) + k*Pi*sin(Pi*Px())*cos(Pi*Py())*sin(Pi*Pz()) )
                               -(k1/Pi)*( (Px()-1)*(Pz()-2)*Pz()*cos(Pi*Px())*sin(Pi*Pz())*( 2*(Py()-1)*sin(Pi*Py()) + Pi*Py()*(Py()-2)*cos(Pi*Py()) )
                                          + Px()*(Px()-2)*(Pz()-2)*Pz()*sin(Pi*Px())*sin(Pi*Pz())*( cos(Pi*Py()) - Pi*(Py()-1)*sin(Pi*Py()) )
                                          + Px()*(Px()-2)*(Pz()-1)*sin(Pi*Px())*cos(Pi*Pz())*( 2*(Py()-1)*sin(Pi*Py()) + Pi*Py()*(Py()-2)*cos(Pi*Py()) )),
                               (3/2)*( dkdz*sin(Pi*Px())*sin(Pi*Py())*sin(Pi*Pz()) + k*Pi*sin(Pi*Px())*sin(Pi*Py())*cos(Pi*Pz()) )
                               -(k1/Pi)*( (Px()-1)*(Py()-2)*Py()*cos(Pi*Px())*sin(Pi*Py())*( 2*(Pz()-1)*sin(Pi*Pz()) + Pi*Pz()*(Pz()-2)*cos(Pi*Pz()) )
                                          + Px()*(Px()-2)*(Py()-1)*sin(Pi*Px())*cos(Pi*Py())*( 2*(Pz()-1)*sin(Pi*Pz()) + Pi*Pz()*(Pz()-2)*cos(Pi*Pz()) )
                                          + Px()*(Px()-2)*Py()*(Py()-2)*sin(Pi*Px())*sin(Pi*Py())*( cos(Pi*Pz()) - Pi*(Pz()-1)*sin(Pi*Pz()) )) );
#endif

    auto K = k*eye<2>();
    auto Lambda = lambda*eye<2>();

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
    for(int i=0; i<nb_refine; i++)
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
                                                         //_dim = 2,
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

            auto u_rt = U_rt.template element<0>( "u" ); //velocity field
            auto v_rt = V_rt.template element<0>( "v" ); // potential field
            auto p_rt = U_rt.template element<1>( "p" );
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

            // Solve problem
            backend(_rebuild=true)->solve( _matrix=M_rt, _solution=U_rt, _rhs=F_rt );

            std::cout << "[Darcy] RT solve done" << std::endl;

            // ****** Compute error ******
            auto l2err_u = normL2( _range=elements(mesh), _expr=u_exact - idv(u_rt) );
            auto l2err_p = normL2( _range=elements(mesh), _expr=p_exact - idv(p_rt) );
            auto l2err_divu = normL2( _range=elements(mesh), _expr=divu_exact - divv(u_rt) );

            auto h1err_u = normH1( _range=elements(mesh), _expr=u_exact - idv(u_rt), _grad_expr=gradu_exact - gradv(u_rt) );
            auto h1err_p = normH1( _range=elements(mesh), _expr=p_exact - idv(p_rt), _grad_expr=gradp_exact - trans(gradv(p_rt)) );
            auto divu_rt = vf::project(Xh, elements(mesh), divv(u_rt) ); //TODO : Raviart-Thomas interpolant !
            auto h1err_divu = normH1( _range=elements(mesh), _expr=divu_exact - divv(u_rt), _grad_expr=graddivu_exact - trans(gradv(divu_rt)) );

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
                    std::cout << "[" << i << "]||divu_exact - divu||_H1 = " << h1err_divu << std::endl;

                    cvg_u << current_hsize << "\t" << nDofu << "\t" << l2err_u << "\t" << h1err_u << "\n";
                    cvg_p << current_hsize << "\t" << nDofp << "\t" << l2err_p << "\t" << h1err_p << "\n";
                    cvg_divu << current_hsize << "\t" << nDofu << "\t" << l2err_divu << "\t" << h1err_divu << "\n";
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
