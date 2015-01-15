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
   \file curl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2014-03-05
 */

#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description testhcurloptions( "test h_curl options" );
    testhcurloptions.add_options()
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ( "penaldir", po::value<double>()->default_value( 50 ), "penalization term)" )
        ;
    return testhcurloptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "curl-curl" ,
                     "curl-curl" ,
                     "0.1",
                     "convergence test for curl-curl formulation (using Hcurl conforming elts)",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

template<int Dim, int OrderU, int OrderP>
class CurlFormulation
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
    typedef bases<Nedelec<OrderU,NedelecKind::NED1> > N_basis_type;
    typedef bases<Lagrange<OrderP,Vectorial> > lagrange_basis_v_type; //Lagrange vectorial space
    typedef bases<Lagrange<OrderP,Scalar> > lagrange_basis_s_type; //Lagrange scalar space

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, lagrange_basis_s_type> lagrange_space_s_type;
    typedef FunctionSpace<mesh_type, lagrange_basis_v_type> lagrange_space_v_type;
    typedef FunctionSpace<mesh_type, N_basis_type> N_space_type;

    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<lagrange_space_s_type> lagrange_space_s_ptrtype;
    typedef boost::shared_ptr<lagrange_space_v_type> lagrange_space_v_ptrtype;
    typedef boost::shared_ptr<N_space_type> N_space_ptrtype;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    CurlFormulation()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        this->changeRepository( boost::format( "/benchmark_curl/%1%/h_%2%/" )
                                % this->about().appName()
                                % doption("hsize")
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

}; //CurlFormulation

template<int Dim, int OrderU, int OrderP>
void
CurlFormulation<Dim, OrderU, OrderP>::convergence(int nb_refine)
{
    int proc_rank = Environment::worldComm().globalRank();

#if defined(FEELPP_2D)

    //rhs
    auto f = vec( 3-Py()*Py(), 3-Px()*Px() ); //f = curl(curl(u_exact)) + u_exact
    auto u_exact = vec( 1-Py()*Py(), 1-Px()*Px() );
    auto gradu_exact = mat<Dim,Dim>( cst(0.), -2*Py(), cst(0.), -2*Px() );

#elif defined(FEELPP_3D)

    std::cout << "3D convergence study not yet implemented!" << std::endl;
    auto f = vec( cst(0.), cst(0.), cst(0.) );

    auto u_exact = vec( cst(0.), cst(0.), cst(0.) );
    auto gradu_exact = mat<Dim,Dim>( cst(0.), cst(0.), cst(0.),
                                     cst(0.), cst(0.), cst(0.),
                                     cst(0.), cst(0.), cst(0.) );

#endif

    std::ofstream cvg_u, cvg_projL2;
    if( proc_rank == 0 )
        {
            // Open
            cvg_u.open( "convergence_u.dat", std::ios::out | std::ios::trunc);
            cvg_projL2.open( "convergence_projL2.dat", std::ios::out | std::ios::trunc);

            // Head
            cvg_u << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\t" << "h1err" << "\n";
            cvg_projL2 << "hsize" << "\t" << "nDof" << "\t" << "l2err" << "\n";
        }

    double current_hsize = meshSize;
    for(int i=0; i<nb_refine; i++)
        {
            mesh_ptrtype mesh;
            std::string mesh_name = (boost::format( "%1%-%2%D" ) % "hypercube" % Dim ).str();

            if( !fs::exists(mesh_name+".msh") )
                {
                    LOG(INFO) << "[CurlFormulation] Mesh has been created \n";
                    mesh = createGMSHMesh( _mesh=new mesh_type,
                                           _desc=domain( _name = mesh_name ,
                                                         _shape = "hypercube",
                                                         //_usenames = true,
                                                         _dim = Dim,
                                                         _h = meshSize,
                                                         _xmin=-1,_xmax=1,
                                                         _ymin=-1,_ymax=1,
                                                         _zmin=-1,_zmax=1 ) );

                    LOG(INFO) << "[CurlFormulation] Mesh has been created \n";
                }
            else
                {
                    //std::cout << "loadGMSHmesh" << std::endl;
                    LOG(INFO) << "[CurlFormulation] Mesh has been loaded (refine level = " << i << ") \n";
                    mesh = loadGMSHMesh( _mesh=new mesh_type,
                                         _filename=mesh_name+".msh",
                                         _refine=i );
                }

            // ****** curl-curl formulation - Hcurl ******
            // auto Vh = Pchv<OrderP>( mesh );
            // auto Xh = Ned1h<OrderU>( mesh );
            auto Vh = Pchv<1>( mesh );
            auto Xh = Ned1h<0>( mesh );
            auto u = Xh->element();
            auto phi = Xh->element();
            //auto penaldir = doption("penaldir");
            auto penaldir = 50;

            auto v = Vh->element( u_exact );

            //variationnal formulation : curl(curl(u)) + u = f
            auto l = form1( _test=Xh );
            l = integrate( _range=elements(mesh), _expr=trans(f)*id(phi) );
            //Dirichlet bc on weak form
            l += integrate(boundaryfaces(mesh), - curlx(phi)*cross(idv(v),N())
                           + penaldir*trans( cross(idv(v),N()) )*cross(id(phi),N())/hFace() );

            auto a = form2( _test=Xh, _trial=Xh );
            a = integrate(elements(mesh), curlxt(u)*curlx(phi) + trans(idt(u))*id(phi) );

            //Dirichlet bc on weak form
            a += integrate(boundaryfaces(mesh), -curlxt(u)*cross(id(phi),N())
                           - curlx(phi)*cross(idt(u),N())
                           + penaldir*trans( cross(idt(u),N()) )*cross(id(phi),N())/hFace() );

            a.solve( _solution=u, _rhs=l, _rebuild=true);

            // ****** Compute error ******
            auto l2err_u = normL2( _range=elements(mesh), _expr=u_exact - idv(u) );
            auto h1err_u = normH1( _range=elements(mesh), _expr=u_exact - idv(u), _grad_expr=gradu_exact - gradv(u) );

            if( proc_rank == 0 )
                {
                    std::cout << "[" << i << "]||u_exact - u||_L2 = " << l2err_u << std::endl;
                    std::cout << "[" << i << "]||u_exact - u||_H1 = " << h1err_u << std::endl;

                    cvg_u << current_hsize << "\t" << Xh->nDof() << "\t" << l2err_u << "\t" << h1err_u << "\n";
                }

            // ****** Export results ******
            std::string exportName =  ( boost::format( "%1%-refine-%2%" ) % this->about().appName() % i ).str();
            std::string uName = ( boost::format( "sol-refine-%1%" ) % i ).str();
            std::string u_exName = ( boost::format( "sol-ex-refine-%1%" ) % i ).str();

            //auto u_ex_proj = vf::project(Vh, elements(mesh), u_exact );

            export_ptrtype exporter_cvg( export_type::New( this->vm(), exportName ) );

            exporter_cvg->step( 0 )->setMesh( mesh );
            exporter_cvg->step( 0 )->add( uName, u );
            exporter_cvg->step( 0 )->add( u_exName, v );
            exporter_cvg->save();

            current_hsize /= 2.0;
        }

    if( proc_rank == 0 )
        {
            cvg_u.close();
            cvg_projL2.close();
        }
}
