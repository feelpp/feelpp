/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-01-09

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2009-2013 Feel++ Consortium

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
   \file steadyns.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-02-05
 */
#include <feel/feel.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feel.hpp>

#include <boost/noncopyable.hpp>
#include <boost/signals2.hpp>
#include <boost/format.hpp>

#include <feel/feelpoly/im.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>




inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description  steadynsoptions( "Navier Stokes problem options" );
    steadynsoptions.add_options()
        ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
        ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
        ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
        ( "rho", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
        ( "hsize", Feel::po::value<double>()->default_value( 0.01 ), "first h value to start convergence" )
        ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
        ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
        ( "beta", Feel::po::value<double>()->default_value( 0.0 ), "convection coefficient" )
        ( "2D.u_exact_x", Feel::po::value<std::string>()->default_value( "" ), "velocity first component" )
        ( "2D.u_exact_y", Feel::po::value<std::string>()->default_value( "" ), "velocity second component" )
        ( "2D.u_exact_z", Feel::po::value<std::string>()->default_value( "" ), "velocity third component" )
        ( "2D.p_exact", Feel::po::value<std::string>()->default_value( "" ), "pressure" )
        ( "export-matlab", "export matrix and vectors in matlab" )
        ( "testcase", Feel::po::value<std::string>()->default_value( "kovasnay" ), "test case" )
        ;
    return  steadynsoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "Steady_Ns" ,
                           "Steady_Ns" ,
                           "0.1",
                           "Stokes equation on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2009-2012 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "" );
    return about;

}


namespace Feel
{
using namespace Feel::vf;

/**
 * \class Stokes class
 * \brief solves the stokes equations
 *
 */

class Steady_Ns
    :
public Application
{
    typedef Application super;
public:

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;
    typedef Simplex<2> convex_type;

    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/

    typedef Lagrange<2, Vectorial> basis_u_type;
    typedef Lagrange<1, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;
#if 0
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
#elif  1
    typedef bases<basis_u_type,basis_p_type> basis_type;
#endif

    /*space*/
    //# marker2 #
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //# endmarker2 #

    /* functions */
    //# marker3 #
    typedef typename space_type::element_type element_type;
    //# endmarker3 #

    /* BDF */
    typedef Bdf<space_type >  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;


    /* export */
    typedef Exporter<mesh_type> export_type;

    FEELPP_DONT_INLINE
    Steady_Ns( );

    // init mesh and space
    FEELPP_DONT_INLINE
    void init();

    /**
     * run the convergence test
     */
    FEELPP_DONT_INLINE
    void run();

    //Export
    FEELPP_DONT_INLINE void exportResults( element_type& U, double t );
    FEELPP_DONT_INLINE void exportResults( element_type& U, int i );

private:

    backend_ptrtype M_backend;
    double meshSize;

    double mu;
    double rho;

    double penalbc;

    mesh_ptrtype mesh;
    space_ptrtype Vh;
    sparse_matrix_ptrtype M,D;
    vector_ptrtype F;

    /*Timer*/
    std::map<std::string,std::pair<boost::timer,double> > timers;
    std::ofstream timings;

    /*BDF */
    bdf_ptrtype M_bdf;

    boost::shared_ptr<export_type> exporter;
}; // Steady_Ns

Steady_Ns::Steady_Ns( )
    :
    super( ),
    M_backend( backend_type::build( soption("backend") ) ),
    meshSize( this->vm()["hsize"].as<double>() ),
    mu( this->vm()["mu"].as<value_type>() ),
    rho( this->vm()["rho"].as<value_type>() ),
    penalbc( this->vm()["bccoeff"].as<value_type>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{

}

void Steady_Ns::exportResults( element_type& U, double t )
{
    exporter->step( t )->setMesh ( mesh );
    exporter->step( t )->add( "u", U. element<0>() );
    exporter->step( t )->add( "p", U. element<1>() );
    exporter->save();
}
void Steady_Ns::exportResults( element_type& U, int i )
{
    exporter->step( i )->setMesh( U.functionSpace()->mesh());
    exporter->step( i )->add( "u", U. element<0>() );
    exporter->step( i )->add( "p", U. element<1>() );
    exporter->save();
}



void Steady_Ns::init()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    mesh = loadMesh( _mesh=new mesh_type );
    if ( Environment::isMasterRank() )
    {
    std::cout << "number of elements of 2D: " << mesh->numElements() << "\n";
    }
    Vh = space_type::New( mesh );
}

void Steady_Ns::run()
{
    this->init();

    auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            auto U = Vh->element( "(u,p)" );
            auto V = Vh->element( "(v,q)" );
            auto u = U.element<0>( "u" );
            auto v = V.element<0>( "u" );
            auto p = U.element<1>( "p" );
            auto q = V.element<1>( "p" );
#if 0
            auto lambda = U.element<2>();
            auto nu = V.element<2>();
#endif
            auto mu=0.035;
            auto rho=1056;

            auto deft = sym(gradt( u ));
            auto def = sym(grad( v ));

            if (!J) J = backend()->newMatrix( Vh, Vh );
            auto a = form2( _test=Vh, _trial=Vh, _matrix=J );

            a = integrate( elements( mesh ), 2*mu*inner(deft,def) );
            a += integrate( elements( mesh ), - id(q)*divt(u) -idt(p)*div(v) );
            // Convective terms
#if 1       //Beta=u
            a += integrate( elements( mesh ), rho*trans(id(v))*gradv(u)*idt(u));
            a += integrate( elements( mesh ), rho*trans(id(v))*gradt(u)*idv(u));
#elif 0     //Beta=u_exact
            std::string u1_str = option(_name="2D.u_exact_x").as<std::string>();
            std::string u2_str = option(_name="2D.u_exact_y").as<std::string>();
            auto vars=symbols<2>();
            auto u1 = parse( u1_str, vars );
            auto u2 = parse( u2_str, vars );
            matrix u_exact_g = matrix(2,1);
            u_exact_g = u1,u2 ;
            auto u_exact = expr<2,1,7>( u_exact_g, vars, "u_exact" );

            a += integrate( elements( mesh ), trans(id(v))*gradt(u)*u_exact);
#endif

#if 0       //Use lagrange multiplier
            a += integrate(elements(mesh), id(q)*idt(lambda)+idt(p)*id(nu));
#endif

#if 0       //Weak Dirichlet conditions on all the boundary faces
            a += integrate( boundaryfaces( mesh ),-trans( -idt(p)*N()+mu*gradt(u)*N() )*id( v ));
            a += integrate( boundaryfaces( mesh ),-trans( -id(p)*N()+mu*grad(u)*N() )*idt( u ));
            a += integrate( boundaryfaces( mesh ), +penalbc*inner( idt( u ),id( v ) )/hFace() );
#endif

#if 0     //Neumann BC on inlet and outlet and Dirichlet BC on the wall
            a += integrate( markedfaces(mesh,"wall1"),-trans( -idt(p)*N()+mu*gradt(u)*N() )*id( v ));
            a += integrate( markedfaces(mesh,"wall2"),-trans( -idt(p)*N()+mu*gradt(u)*N() )*id( v ));
            std::cout << "Neumann " << "\n";
            a += integrate( markedfaces(mesh,"wall1"),-trans( -id(p)*N()+mu*grad(u)*N() )*idt( u ));
            a += integrate( markedfaces(mesh,"wall2"),-trans( -id(p)*N()+mu*grad(u)*N() )*idt( u ));

            a += integrate( markedfaces(mesh,"wall1"), +penalbc*trans( idt( u ) )*id( v )/hFace() );
            a += integrate( markedfaces(mesh,"wall2"), +penalbc*trans( idt( u ) )*id( v )/hFace() );
#endif

#if 1     //(For Passerni test case)Dirichlet Boundary condition on the inlet and wall and stress free on the outlet
            a += integrate( markedfaces(mesh,"wall"),-trans( -idt(p)*N()+2*mu*deft*N() )*id( v ));
            a += integrate( markedfaces(mesh,"wall"),-trans( -id(p)*N()+2*mu*def*N() )*idt( u ));
            a += integrate( markedfaces(mesh,"wall"), +penalbc*trans( idt( u ) )*id( v )/hFace() );

            a += integrate( markedfaces(mesh,"inlet"),-trans( -idt(p)*N()+2*mu*deft*N() )*id( v ));
            a += integrate( markedfaces(mesh,"inlet"),-trans( -id(p)*N()+2*mu*def*N() )*idt( u ));
            a += integrate( markedfaces(mesh,"inlet"), +penalbc*trans( idt( u ) )*id( v )/hFace() );

#endif


            //Time
            a+= integrate(elements(mesh), trans(rho*idt(u))*id(v)*M_bdf->polyDerivCoefficient(0));


        };

    auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            auto U = Vh->element( "(u,p)" );
            auto V = Vh->element( "(v,q)" );
            auto u = U.element<0>( "u" );
            auto v = V.element<0>( "u" );
            auto p = U.element<1>( "p" );
            auto q = V.element<1>( "p" );
#if 0
            auto lambda = U.element<2>();
            auto nu = V.element<2>();
#endif



            /////////////////////////////// Use Ginac //////////////////////////////////////
            std::string u1_str = option(_name="2D.u_exact_x").as<std::string>();
            std::string u2_str = option(_name="2D.u_exact_y").as<std::string>();
            std::string p_str = option(_name="2D.p_exact").as<std::string>();

             LOG(INFO) << "ux = " << u1_str<<"\n";
             LOG(INFO) << "uy = " << u2_str<<"\n";
             LOG(INFO) << "p = " << p_str<<"\n";
            auto vars=symbols<2>();
            auto u1 = parse( u1_str, vars );
            auto u2 = parse( u2_str, vars );
            matrix u_exact_g = matrix(2,1);
            u_exact_g = u1,u2 ;
            auto p_exact_g = parse( p_str, vars );

            LOG(INFO) << "u_exact = " << u_exact_g;
            LOG(INFO) << "p_exact = " << p_exact_g;

            auto u_exact = expr<2,1,7>( u_exact_g, vars, "u_exact" );
            auto p_exact = expr<7>( p_exact_g, vars, "p_exact" );
            auto beta=u_exact;

            auto gradu_exact_g = grad( u_exact_g, vars );
            auto divu_exact_g = div( u_exact_g, vars );
            LOG(INFO) << "gradu_exact_g = " << gradu_exact_g;
            LOG(INFO) << "divu_exact_g = " << divu_exact_g;
            auto gradu_exact = expr<2,2,7>( gradu_exact_g, vars, "gradu_exact" );
            auto divu_exact = expr<1,1,7>( divu_exact_g, vars, "divu_exact" );
            auto convection=gradu_exact*beta;

            auto f_g = rho*gradu_exact_g*u_exact_g -2*mu*laplacian( u_exact_g, vars ) + grad( p_exact_g, vars ).transpose();
            auto f = expr<2,1,7>( f_g, vars, "f" );
            LOG(INFO) << "f = " << f_g << "\n";
            ///////////////////////////////////////////////////////////////////////////

            U=*X;

            auto r = form1( _test=Vh, _vector=R );
            r = integrate( elements( mesh ),-inner( f,id( v ) ) );
            r += integrate( elements( mesh ), trace(trans(2*mu*gradv( u ))*grad( v )) );
            r +=  integrate( elements( mesh ),-idv(p)*div(v) - id(q)*divv(u));
            // convective terms
#if 1       //Beta=u
            r += integrate( elements( mesh ), trans(rho*gradv( u )*idv(u))*id(v));
#elif 0     //Beta=u_exact
            r += integrate( elements( mesh ), trans(gradv( u )*u_exact)*id(v));
#endif


#if 0       //Lagrange multiplier
            r += integrate ( elements( mesh ), +id( q )*idv( lambda )+(idv( p )-p_exact)*id( nu ) );
#endif

            auto SigmaNv = ( -idv( p )*N() + 2*mu*sym(gradv( u ))*N() );
            auto SigmaN = ( -id( q )*N() + 2*mu*sym(grad( v ))*N() );

#if 0       //Weak Dirichlet BC on all the boundary faces
            r +=integrate ( boundaryfaces(mesh), - trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) - u_exact ) + penalbc*trans( idv( u ) - u_exact )*id( v )/hFace() );
#endif

#if 0     //Neumann BC on the inlet and outlet faces and Dirichlet BC on the wall
            r +=integrate ( markedfaces(mesh,"inlet"),  -trans( -P_inlet*N() )*id( v ) );
            r +=integrate ( markedfaces(mesh,"outlet"), -trans( -P_outlet*N())*id( v ) );
            r +=integrate ( markedfaces(mesh,"wall1"),  -trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) - u_exact ) + penalbc*trans( idv( u ) - u_exact )*id( v )/hFace() );
            r +=integrate ( markedfaces(mesh,"wall2"),  -trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) - u_exact ) + penalbc*trans( idv( u ) - u_exact )*id( v )/hFace() );
#endif

#if 1    //(For Passerni test case)Dirichlet Boundary condition on the inlet and wall and stress free on the outlet
            r +=integrate ( markedfaces(mesh,"wall"),  -trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) - u_exact ) + penalbc*trans( idv( u ) - u_exact )*id( v )/hFace() );
            r +=integrate ( markedfaces(mesh,"inlet"),  -trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) - u_exact ) + penalbc*trans( idv( u ) - u_exact )*id( v )/hFace() );
#endif



            //Time
            auto bdf_poly = M_bdf->polyDeriv();
            auto bdfu_poly = bdf_poly.element<0>();
            r += integrate( elements( mesh ), -rho*trans(idv(bdfu_poly))*id(v));
            r += integrate( elements( mesh ), rho*trans(idv(u))*id(v)*M_bdf->polyDerivCoefficient(0));

        };


    auto U = Vh->element( "(u,p)" );
    auto V = Vh->element( "(v,q)" );
    auto u = U.element<0>( "u" );
    auto v = V.element<0>( "u" );
    auto p = U.element<1>( "p" );
    auto q = V.element<1>( "p" );
#if 0
    auto lambda = U.element<2>();
    auto nu = V.element<2>();
#endif


    /////////////////////////////// Use Ginac //////////////////////////////////////
    std::string u1_str = option(_name="2D.u_exact_x").as<std::string>();
    std::string u2_str = option(_name="2D.u_exact_y").as<std::string>();
    std::string p_str = option(_name="2D.p_exact").as<std::string>();
    auto vars=symbols<2>();
    auto u1 = parse( u1_str, vars );
    auto u2 = parse( u2_str, vars );
    matrix u_exact_g = matrix(2,1);
    u_exact_g = u1,u2 ;
    auto p_exact_g = parse( p_str, vars );
    auto u_exact = expr<2,1,7>( u_exact_g, vars, "u_exact" );
    auto p_exact = expr<7>( p_exact_g, vars, "p_exact" );
    auto beta=u_exact;

    auto gradu_exact_g = grad( u_exact_g, vars );
    auto divu_exact_g = div( u_exact_g, vars );
    auto gradu_exact = expr<2,2,7>( gradu_exact_g, vars, "gradu_exact" );
    auto divu_exact = expr<1,1,7>( divu_exact_g, vars, "divu_exact" );
    auto convection=gradu_exact*beta;

    //auto f_g = gradu_exact_g*u_exact_g -mu*laplacian( u_exact_g, vars ) + grad( p_exact_g, vars ).transpose();
    //auto f = expr<2,1,7>( f_g, vars, "f" );
    auto f=vec(cst(0.),cst(0.));
    ///////////////////////////////////////////////////////////////////////////////

    u=vf::project(Vh->template functionSpace<0>(), elements(mesh), zero<2,1>());
    p=vf::project(Vh->template functionSpace<1>(), elements(mesh), constant(0.0));

    //u=vf::project(Vh->functionSpace<0>(), elements(mesh),u_exact);
    //p=vf::project(Vh->functionSpace<1>(), elements(mesh),p_exact);

    // -- INITIALIZATION -- //
    backend()->nlSolver()->residual =Residual;
    backend()->nlSolver()->jacobian =Jacobian;

    M_bdf=bdf(_space=Vh);
    M_bdf->initialize(U);
    M_bdf->start();
    // Temporal loop
    for( M_bdf->start() ; M_bdf->isFinished()==false ; M_bdf->next() )
        {
            std::cout << "Time : " << M_bdf->time() << "\n";
            LOG( INFO ) << "Time : " << M_bdf->time() << "\n";

            timers["solve"].first.restart();
            backend()->nlSolve( _solution = U );
            this->exportResults(U, M_bdf->time() );
            M_bdf->shiftRight(U);
            timers["solve"].second=timers["solve"].first.elapsed();
        }
    U.save(_path=".");
    timers["all"].second=timers["all"].first.elapsed();




    // Solving
    //backend()->nlSolve( _solution=U );

#if 0
    u=vf::project(Vh->functionSpace<0>(), elements(mesh), u_exact);
    p=vf::project(Vh->functionSpace<1>(), elements(mesh), p_exact);

    backend()->nlSolve( _solution=U );
#endif


    double u_errorL2 = integrate( elements( u.mesh() ), trans( idv( u )-u_exact )*( idv( u )-u_exact ) ).evaluate()( 0, 0 );
    std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";

    double p_errorL2 = integrate( elements( u.mesh() ), ( idv( p ) - p_exact )*( idv( p )-p_exact )).evaluate()( 0, 0 );
    std::cout << "||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";

    v = vf::project( u.functionSpace(), elements( u.mesh() ), u_exact );
    q = vf::project( p.functionSpace(), elements( p.mesh() ), p_exact );
      if ( exporter->doExport() )
          {
              exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
              exporter->step( 0 )->addRegions();
              exporter->step( 0 )->add( "u", U.element<0>() );
              exporter->step( 0 )->add( "p", U.element<1>() );
              exporter->save();
                }
};
}


int main( int argc, char** argv )
{

    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="steady_ns",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );
    Feel::Steady_Ns Steady_Ns;
    Steady_Ns.run();
}
