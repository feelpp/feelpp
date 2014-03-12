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
    typedef Simplex<2,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    typedef bases<RaviartThomas<OrderU> > RT_basis_type; //RT vectorial space
    typedef bases<Lagrange<OrderP,Vectorial> > lagrange_basis_v_type; //Lagrange vectorial space
    typedef bases<Lagrange<OrderP,Scalar> > lagrange_basis_s_type; //Lagrange scalar space
    typedef bases<Lagrange<0,Scalar> > lagrange_multiplier_basis_type; //P0 scalar space
    typedef bases< RaviartThomas<0>, Lagrange<1,Scalar> > basis_type; //For Darcy : (u,p) (\in H_div x L2)

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
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        std::cout << "[TestHDiv]\n";
        this->changeRepository( boost::format( "/benchmark_darcy/%1%/h_%2%/" )
                                % this->about().appName()
                                % this->vm()["hsize"].template as<double>()
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

}; //Darcy


template<int Dim, int OrderU, int OrderP>
void
Darcy<Dim, OrderU, OrderP>::convergence(int nb_refine)
{
    auto Pi = M_PI;

    auto alpha = cst(1.); //temporarly
    auto k = alpha;
    auto K = k*eye<2>();
    auto lambda = cst(1.)/alpha;
    auto Lambda = lambda*eye<2>();
    auto f = sin( Pi*Px() )*sin( Pi*Py() );

    // Coeff for stabilization terms
    auto d1 = cst(1.0/2.0);
    auto d2 = cst(1.0/2.0);

    std::ofstream cvg_p("convergence_p.dat", ios::out | ios::trunc);
    std::ofstream cvg_u("convergence_u.dat", ios::out | ios::trunc);
    std::ofstream cvg_divu("convergence_divu.dat", ios::out | ios::trunc);
    std::ofstream cvg_projL2("convergence_projL2.dat", ios::out | ios::trunc);
    std::ofstream cvg_projHDIV("convergence_projHDIV.dat", ios::out | ios::trunc);

    cvg_u << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\t" << "h1err" << "\n";
    cvg_p << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\t" << "h1err" << "\n";
    cvg_divu << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\t" << "h1err" << "\n";
    cvg_projL2 << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\n";
    cvg_projHDIV << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\n";

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
                                                         _dim = 2,
                                                         _h = meshSize,
                                                         _xmin=0,_xmax=2,
                                                         _ymin=0,_ymax=2 ) );
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
            lagrange_space_s_ptrtype Xh = lagrange_space_s_type::New( mesh );
            lagrange_space_v_ptrtype Xhvec = lagrange_space_v_type::New( mesh );

            auto U_rt = Yh->element( "(u,p)" ); //trial
            auto V_rt = Yh->element( "(v,q)" ); //test

            auto u_rt = U_rt.element<0>( "u" ); //velocity field
            auto v_rt = V_rt.element<0>( "v" ); // potential field
            auto p_rt = U_rt.element<1>( "p" );
            auto q_rt = V_rt.element<1>( "q" );

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
            auto u_exact = -(1/(2*Pi))*vec(cos( Pi*Px() )*sin( Pi*Py() ), sin( Pi*Px() )*cos( Pi*Py() ));
            auto p_exact = (1/(2*Pi*Pi))*sin( Pi*Px() )*sin( Pi*Py() );
            auto divu_exact = f;

            auto gradu_exact = (1/2)*mat<2,2>(sin(Pi*Px())*sin(Pi*Py()),-cos(Pi*Px())*cos(Pi*Py()),-cos(Pi*Px())*cos(Pi*Py()),sin(Pi*Px())*sin(Pi*Py()));
            auto gradp_exact = -1*u_exact;
            auto graddivu_exact = Pi*vec( cos(Pi*Px())*sin(Pi*Py()), sin(Pi*Px())*cos(Pi*Py()) );

            auto l2err_u = normL2( _range=elements(mesh), _expr=u_exact - idv(u_rt) );
            auto l2err_p = normL2( _range=elements(mesh), _expr=p_exact - idv(p_rt) );
            auto l2err_divu = normL2( _range=elements(mesh), _expr=divu_exact - divv(u_rt) );

            auto h1err_u = normH1( _range=elements(mesh), _expr=u_exact - idv(u_rt), _grad_expr=gradu_exact - gradv(u_rt) );
            auto h1err_p = normH1( _range=elements(mesh), _expr=p_exact - idv(p_rt), _grad_expr=gradp_exact - trans(gradv(p_rt)) );
            auto divu_rt = vf::project(Xh, elements(mesh), divv(u_rt) ); //TODO : Raviart-Thomas interpolant !
            auto h1err_divu = normH1( _range=elements(mesh), _expr=divu_exact - divv(u_rt), _grad_expr=graddivu_exact - trans(gradv(divu_rt)) );

            std::cout << "[" << i << "]||u_exact - u||_L2 = " << l2err_u << std::endl;
            std::cout << "[" << i << "]||p_exact - p||_L2 = " << l2err_p << std::endl;
            std::cout << "[" << i << "]||divu_exact - divu||_L2 = " << l2err_divu << std::endl;

            std::cout << "[" << i << "]||u_exact - u||_H1 = " << h1err_u << std::endl;
            std::cout << "[" << i << "]||p_exact - p||_H1 = " << h1err_p << std::endl;
            std::cout << "[" << i << "]||divu_exact - divu||_H1 = " << h1err_divu << std::endl;

            // double current_hsize = (i!=0)?meshSize/(2*i):meshSize;
            cvg_u << current_hsize << "\t" << nDofu << "\t" << l2err_u << "\t" << h1err_u << "\n";
            cvg_p << current_hsize << "\t" << nDofp << "\t" << l2err_p << "\t" << h1err_p << "\n";
            cvg_divu << current_hsize << "\t" << nDofu << "\t" << l2err_divu << "\t" << h1err_divu << "\n";

            // ****** Projection Operators (L2 - HDIV) : check proj( div u ) = proj( f ) ******
            auto RTh = Dh<0>( mesh );

            // L2 projection
            auto l2_v = opProjection( _domainSpace=Xhvec, _imageSpace=Xhvec, _type=L2 ); //l2 vectorial proj
            auto l2_s = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=L2 ); //l2 scalar proj
            auto E_l2 = l2_v->project( _expr= trans(u_exact) );
            auto l2error_L2 = normL2( _range=elements(mesh), _expr=divv(E_l2) - f );

            auto hdiv = opProjection( _domainSpace=RTh, _imageSpace=RTh, _type=HDIV ); //hdiv proj (RT elts)
            auto E_hdiv = hdiv->project( _expr= trans(u_exact), _div_expr=f );
            auto l2error_HDIV = normL2( _range=elements(mesh), _expr=divv(E_hdiv) - f );

            cvg_projL2 << current_hsize << "\t" << Xhvec->nDof() << "\t" << l2error_L2 << "\n";
            cvg_projHDIV << current_hsize << "\t" << RTh->nDof() << "\t" << l2error_HDIV << "\n";

            // ****** Export results ******
            auto u_ex_proj = vf::project(Xhvec, elements(mesh), u_exact );
            auto p_ex_proj = vf::project(p_rt.functionSpace(), elements(mesh), p_exact );

            export_ptrtype exporter_cvg( export_type::New( this->vm(),
                                                           ( boost::format( "%1%-refine-%2%" )
                                                             % this->about().appName()
                                                             % i ).str() ) );

            exporter_cvg->step( 0 )->setMesh( mesh );
            exporter_cvg->step( 0 )->add( "velocity", u_rt );
            exporter_cvg->step( 0 )->add( "potential", p_rt );
            exporter_cvg->step( 0 )->add( "velocity_ex", u_ex_proj );
            exporter_cvg->step( 0 )->add( "potential_ex", p_ex_proj );
            exporter_cvg->save();

            current_hsize /= 2.0;
        }

    cvg_u.close();
    cvg_p.close();
    cvg_divu.close();
    cvg_projL2.close();
    cvg_projHDIV.close();
}
