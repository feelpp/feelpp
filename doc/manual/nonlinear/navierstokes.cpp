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
   \file navierstokes.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-02-05
 */
#include <feel/feel.hpp>
#include <feel/feelalg/backend.hpp>
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description  navierstokesoptions( "Navier Stokes problem options" );
     navierstokesoptions.add_options()
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.01 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return  navierstokesoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "navier_stoke" ,
                           "navier_stokes" ,
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

class Navier_Stokes
    :
public Application
{
    typedef Application super;
public:

#ifndef FEELPP_BC
#define  FEELPP_BC Neumann
#endif
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
#elif  (FEELPP_BC == Neumann)
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

    /* export */
    typedef Exporter<mesh_type> export_type;

    FEELPP_DONT_INLINE
    Navier_Stokes( );

    // init mesh and space
    FEELPP_DONT_INLINE
    void init();

    /**
     * run the convergence test
     */
    FEELPP_DONT_INLINE
    void run();

private:
    //FEELPP_DONT_INLINE
    //void exportResults( element_type& u, element_type& v );

private:

    backend_ptrtype M_backend;
    double meshSize;

    double mu;
    double penalbc;

    mesh_ptrtype mesh;
    space_ptrtype Vh;
    sparse_matrix_ptrtype M,D;
    vector_ptrtype F;

    boost::shared_ptr<export_type> exporter;
}; // Navier_Stokes

Navier_Stokes::Navier_Stokes( )
    :
    super( ),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( this->vm()["hsize"].template as<double>() ),
    mu( this->vm()["mu"].template as<value_type>() ),
    penalbc( this->vm()["bccoeff"].template as<value_type>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{

}

void Navier_Stokes::init()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    mesh = createGMSHMesh( _mesh=new mesh_type,
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                            _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % convex_type().dimension() % 1 ).str() ,
                                          _shape="hypercube",
                                          _dim=convex_type().dimension(),
                                          _h=meshSize ) );
    /*
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name="kovaznay",
                                         _usenames=false,
                                         _shape="hypercube",
                                         _h=meshSize,
                                         _xmin=-0.5, _xmax=1,
                                         _ymin=-0.5, _ymax=1.5 ) );*/
     //mesh = unitSquare();
    mesh->addMarkerName( "inlet", 1, 1 );
    mesh->addMarkerName( "outlet", 3, 1 );
    mesh->addMarkerName( "wall1", 2, 1 );
    mesh->addMarkerName( "wall2", 4, 1 );
    std::cout << "number of elements of 2D: " << mesh->numElements() << "\n";

    Vh = space_type::New( mesh );
}

void Navier_Stokes::run()
{
    this->init();

    auto U = Vh->element( "(u,p)" );
    auto V = Vh->element( "(v,q)" );
    auto u = U.template element<0>( "u" );
    auto v = V.template element<0>( "u" );
    auto p = U.template element<1>( "p" );
    auto q = V.template element<1>( "p" );
#if 0
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
#endif


    auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            auto U = Vh->element( "(u,p)" );
            auto V = Vh->element( "(v,q)" );
            auto u = U.template element<0>( "u" );
            auto v = V.template element<0>( "u" );
            auto p = U.template element<1>( "p" );
            auto q = V.template element<1>( "p" );
#if 0
            auto lambda = U.template element<2>();
            auto nu = V.template element<2>();
#endif
            auto mu=1;

            if (!J) J = backend()->newMatrix( Vh, Vh );
            std::cout << "coucou1 " << "\n";
            auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
            a = integrate( elements( mesh ), mu*trace(gradt( u )*trans(grad( v )) ) );
            std::cout << "coucou2 " << "\n";
            a += integrate( elements( mesh ), - id(q)*divt(u) -idt(p)*div(v) );
            std::cout << "coucou3 " << "\n";
            // Convective terms
            //a += integrate( elements( mesh ), trans(id(v))*gradv(u)*idt(u));
            //a += integrate( elements( mesh ), trans(id(v))*gradt(u)*idv(u));
            std::cout << "coucou5 " << "\n";

#if 0
            a += integrate(elements(mesh), id(q)*idt(lambda)+idt(p)*id(nu));
#endif

            //a += integrate(elements(mesh), idt(p)*id(nu));

#if 0
            std::cout << "Dirichlet " << "\n";
            //Weak Dirichlet conditions
            a += integrate( boundaryfaces( mesh ),-trans( -idt(p)*N()+mu*gradt(u)*N() )*id( v ));
            a += integrate( boundaryfaces( mesh ),-trans( -id(p)*N()+mu*grad(u)*N() )*idt( u ));
            a += integrate( boundaryfaces( mesh ), +penalbc*trans(idt( u ))*id( v )/hFace() );
#elif 1
            //Neumann BC
            a += integrate( markedfaces(mesh,"wall1"),-trans( -idt(p)*N()+mu*gradt(u)*N() )*id( v ));
            a += integrate( markedfaces(mesh,"wall2"),-trans( -idt(p)*N()+mu*gradt(u)*N() )*id( v ));
            std::cout << "Neumann " << "\n";
            a += integrate( markedfaces(mesh,"wall1"),-trans( -id(p)*N()+mu*grad(u)*N() )*idt( u ));
            a += integrate( markedfaces(mesh,"wall2"),-trans( -id(p)*N()+mu*grad(u)*N() )*idt( u ));

            a += integrate( markedfaces(mesh,"wall1"), +penalbc*trans( idt( u ) )*id( v )/hFace() );
            a += integrate( markedfaces(mesh,"wall2"), +penalbc*trans( idt( u ) )*id( v )/hFace() );
#endif




        };

    auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            std::cout << "coucou6 " << "\n";
            auto U = Vh->element( "(u,p)" );
            auto V = Vh->element( "(v,q)" );
            auto u = U.template element<0>( "u" );
            auto v = V.template element<0>( "u" );
            auto p = U.template element<1>( "p" );
            auto q = V.template element<1>( "p" );
#if 0
            auto lambda = U.template element<2>();
            auto nu = V.template element<2>();
#endif

            //Solution Kovasnay
            /*auto mu=0.035;
            double lambdaa = 1./( 2.*mu ) - math::sqrt( 1./( 4.*mu*mu ) + 4.*pi*pi );
            auto u1 = 1. - exp( lambdaa * Px() ) * cos( 2.*pi*Py() );
            auto u2 = ( lambdaa/( 2.*pi ) ) * exp( lambdaa * Px() ) * sin( 2.*pi*Py() );
            auto u_exact = vec( u1,u2 ) ;
            auto p_exact = -0.5*exp( 2.*lambdaa*Px() ) ;
            auto f1 = ( -mu*( -lambdaa*lambdaa*exp( lambdaa*Px() )*cos( 2.0*pi*Py() )+4.0*exp( lambdaa*Px() )*cos( 2.0*pi*Py() )*pi*pi )-lambdaa*exp( 2.0*lambdaa*Px() ) );
            auto f2 = ( -mu*( lambdaa*lambdaa*lambdaa*exp( lambdaa*Px() )*sin( 2.0*pi*Py() )/pi/2.0-2.0*lambdaa*exp( lambdaa*Px() )*sin( 2.0*pi*Py() )*pi ) );
            auto f = vec( f1,f2 ); //+ convection;*/

            //Solution Poiseuille
            auto mu=1;
            auto rayon=1;
            auto P_inlet = 1.;
            auto P_outlet = 0.;
            auto L=1;

            auto u_exact=vec(Py()*(1-(Py()))/(2*L),cst(0.) );
            auto p_exact=(P_outlet-P_inlet)*Px() + P_inlet;
            auto f=vec(cst(0.) , cst(0.) );

            U=*X;

            auto r = form1( _test=Vh, _vector=R );
            r = integrate( elements( mesh ),-inner( f,id( v ) ) );
            r += integrate( elements( mesh ), trace(trans(mu*gradv( u ))*grad( v )) );
            r +=  integrate( elements( mesh ),-idv(p)*div(v) - id(q)*divv(u));
            // convective terms
            //r += integrate( elements( mesh ), trans(gradv( u )*idv(u))*id(v));


#if 0
            r += integrate ( elements( mesh ), +id( q )*idv( lambda )+idv( p )*id( nu ) );
#endif

#if 0
            //Weak Dirichlet
            auto SigmaNv = ( -idv( p )*N() + mu*gradv( u )*N() );
            auto SigmaN = ( -id( q )*N() + mu*grad( v )*N() );
            r +=integrate ( boundaryfaces(mesh), - trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) - u_exact ) + penalbc*trans( idv( u ) - u_exact )*id( v )/hFace() );
#elif 1
            //Neumann
            auto SigmaNv = ( -idv( p )*N() + mu*gradv( u )*N() );
            auto SigmaN = ( -id( q )*N() + mu*grad( v )*N() );
            r +=integrate ( markedfaces(mesh,"inlet"),  -trans( -P_inlet*N() )*id( v ) );
            r +=integrate ( markedfaces(mesh,"outlet"), -trans( -P_outlet*N())*id( v ) );
            r +=integrate ( markedfaces(mesh,"wall1"),  -trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) - u_exact ) + penalbc*trans( idv( u ) - u_exact )*id( v )/hFace() );
            r +=integrate ( markedfaces(mesh,"wall2"),  -trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) - u_exact ) + penalbc*trans( idv( u ) - u_exact )*id( v )/hFace() );
#endif

        };

    //Poiseuille
    auto rayon=1;
    auto P_inlet = 1.;
    auto P_outlet = 0.;
    auto L=1;
    auto u_exact=vec(Py()*(1-Py())/(2*L),cst(0.) );
    auto p_exact=(P_outlet-P_inlet)*Px() + P_inlet;
    auto f=vec(cst(0.) , cst(0.) );

    //Solution Kovasnay
    /*auto mu=0.035;
    double lambdaa = 1./( 2.*mu ) - math::sqrt( 1./( 4.*mu*mu ) + 4.*pi*pi );
    auto u1 = 1. - exp( lambdaa * Px() ) * cos( 2.*pi*Py() );
    auto u2 = ( lambdaa/( 2.*pi ) ) * exp( lambdaa * Px() ) * sin( 2.*pi*Py() );
    auto u_exact = vec( u1,u2 ) ;
    auto p_exact = -0.5*exp( 2.*lambdaa*Px() );
    auto f1 = ( -mu*( -lambdaa*lambdaa*exp( lambdaa*Px() )*cos( 2.0*pi*Py() )+4.0*exp( lambdaa*Px() )*cos( 2.0*pi*Py() )*pi*pi )-lambdaa*exp( 2.0*lambdaa*Px() ) );
    auto f2 = ( -mu*( lambdaa*lambdaa*lambdaa*exp( lambdaa*Px() )*sin( 2.0*pi*Py() )/pi/2.0-2.0*lambdaa*exp( lambdaa*Px() )*sin( 2.0*pi*Py() )*pi ) );
    auto f = vec( f1,f2 ); //+ convection;*/

    u=vf::project(Vh->functionSpace<0>(), elements(mesh), zero<2,1>());
    p=vf::project(Vh->functionSpace<1>(), elements(mesh), constant(0.0));

    backend()->nlSolver()->residual = Residual;
    backend()->nlSolver()->jacobian = Jacobian;
    backend()->nlSolve( _solution=U );

#if 0
    u=vf::project(Vh->functionSpace<0>(), elements(mesh), u_exact);
    p=vf::project(Vh->functionSpace<1>(), elements(mesh), p_exact);

    backend()->nlSolve( _solution=U );
#endif


    vector_ptrtype R = backend()->newVector( Vh );
    vector_ptrtype X = backend()->newVector( Vh );
    *X = U;
    Residual( X,R );
    std::cout << "residual -= " << R->dot( X ) << "\n";

    //sparse_matrix_ptrtype J = backend()->newMatrix( Vh, Vh );
    //Jacobian( X,J );
    //backend()->solve( _matrix=J, _solution=U, _rhs=R );





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
                exporter->step( 0 )->add( "u", U.template element<0>() );
                exporter->step( 0 )->add( "p", U.template element<1>() );
                exporter->step( 0 )->add( "u_exact", v);
                exporter->step( 0 )->add( "p_exact", q );
                exporter->save();
            }
        //this->exportResults( U, V );
};
}


int main( int argc, char** argv )
{

    using namespace Feel;
#if 0
    auto appli_name="steady_ns_dirichlet";
#elif 1
    auto appli_name="steady_ns_neumann";
#endif
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name= appli_name,
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );
    Feel::Navier_Stokes Navier_Stokes;
    Navier_Stokes.run();
}
