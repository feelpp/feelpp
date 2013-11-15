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
   \file drivencavity.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-02-05
 */
#include <feel/feel.hpp>
#include <feel/feelalg/backend.hpp>
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description  drivencavityoptions( "Driven Cavity problem options" );
     drivencavityoptions.add_options()
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return  drivencavityoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "driven_cavity" ,
                           "driven_cavity" ,
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

class Driven_Cavity
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
    //#if defined( FEELPP_USE_LM )
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
    //#else
    //typedef bases<basis_u_type,basis_p_type> basis_type;
    //#endif


    typedef bases<basis_u_type> basis_type_U;
    typedef FunctionSpace<mesh_type, basis_type_U> space_type_U;
    typedef boost::shared_ptr<space_type_U> space_ptrtype_U;

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
    Driven_Cavity( );

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
    space_ptrtype_U Uh;
    sparse_matrix_ptrtype M,D;
    vector_ptrtype F;

    boost::shared_ptr<export_type> exporter;
}; // Driven_Cavity

Driven_Cavity::Driven_Cavity( )
    :
    super( ),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( this->vm()["hsize"].template as<double>() ),
    mu( this->vm()["mu"].template as<value_type>() ),
    penalbc( this->vm()["bccoeff"].template as<value_type>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{

}

void Driven_Cavity::init()
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
     //mesh = unitSquare();
    mesh->addMarkerName( "inlet", 1, 1 );
    mesh->addMarkerName( "outlet", 3, 1 );
    mesh->addMarkerName( "wall1", 2, 1 );
    mesh->addMarkerName( "wall2", 4, 1 );
    std::cout << "number of elements of 2D: " << mesh->numElements() << "\n";

    Vh = space_type::New( mesh );
}

void Driven_Cavity::run()
{
    this->init();

    auto U = Vh->element( "(u,p)" );
    auto V = Vh->element( "(v,q)" );
    auto u = U.template element<0>( "u" );
    auto v = V.template element<0>( "u" );
    auto p = U.template element<1>( "p" );
    auto q = V.template element<1>( "p" );
    // #if defined( FEELPP_USE_LM )
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
    //#endif


    auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            auto U = Vh->element( "(u,p)" );
            auto V = Vh->element( "(v,q)" );
            auto u = U.template element<0>( "u" );
            auto v = V.template element<0>( "u" );
            auto p = U.template element<1>( "p" );
            auto q = V.template element<1>( "p" );
            //#if defined( FEELPP_USE_LM )
            auto lambda = U.template element<2>();
            auto nu = V.template element<2>();
            auto mu=0.001;
            //#endif

            if (!J) J = backend()->newMatrix( Vh, Vh );
            std::cout << "coucou1 " << "\n";
            auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
            a += integrate( elements( mesh ), mu*inner(gradt( u ),grad( v )) );
            std::cout << "coucou2 " << "\n";
            a += integrate( elements( mesh ), id(q)*divt(u) -idt(p)*div(v) );
            std::cout << "coucou3 " << "\n";
            // Convective terms
            a += integrate( elements( mesh ), trans(id(v))*gradv(u)*idt(u));
            a += integrate( elements( mesh ), trans(id(v))*gradt(u)*idv(u));
            std::cout << "coucou5 " << "\n";

            //#if defined( FEELPP_USE_LM )
            a += integrate(elements(mesh), id(q)*idt(lambda)+idt(p)*id(nu));
            //#elif
            //a += integrate(elements(mesh), idt(p)*id(nu));

            //Weak Dirichlet conditions
            a += integrate( boundaryfaces( mesh ),-trans( -idt(p)*N()+mu*gradt(u)*N() )*id( v ));//
            a += integrate( boundaryfaces( mesh ),-trans( -id(p)*N()+mu*grad(u)*N() )*idt( u ));//
            a += integrate( boundaryfaces( mesh ), +penalbc*inner( idt( u ),id( v ) )/hFace() );

          

        };

    auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            std::cout << "coucou6 " << "\n";
            auto U = Vh->element( "(u,p)" );
            auto V = Vh->element( "(v,q)" );
            auto u = U.template element<0>( "u" );
            //auto u_exact = U.template element<0>( "u_exact" );
            auto v = V.template element<0>( "u" );
            auto p = U.template element<1>( "p" );
            auto q = V.template element<1>( "p" );
            //#if defined( FEELPP_USE_LM )
            auto lambda = U.template element<2>();
            auto nu = V.template element<2>();
            //#endif
            auto mu=0.001;
            std::cout << "coucou7 " << "\n";


            auto uex=vec(cst(1.) , cst(0.) );
            auto u_exact=vf::project(Vh->functionSpace<0>(), markedfaces(mesh, 4), uex );


            U=*X;
            std::cout << "coucou8 " << "\n";
            auto r = form1( _test=Vh, _vector=R );
            std::cout << "coucou9 " << "\n";
            //r += integrate( elements( mesh ),-inner( f,id( v ) ) );
            std::cout << "coucou10 " << "\n";
            r += integrate( elements( mesh ), trans(gradv( u )*idv(u))*id(v));
            std::cout << "coucou11 " << "\n";
            r += integrate( elements( mesh ), inner(mu*gradv( u ),grad( v )) );
            std::cout << "coucou12 " << "\n";
            r +=  integrate( elements( mesh ),-idv(p)*div(v) + id(q)*divv(u));
            std::cout << "coucou13 " << "\n";

           
            //#if defined( FEELPP_USE_LM )
            r += integrate ( elements( mesh ), +id( q )*idv( lambda )+idv( p )*id( nu ) );
            //#endif


            //Weak Dirichlet
            auto SigmaNv = ( -idv( p )*N() + mu*gradv( u )*N() );
            auto SigmaN = ( -id( q )*N() + mu*grad( v )*N() );
            r +=integrate ( boundaryfaces(mesh), - trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) -idv(u_exact) ) + penalbc*trans( idv( u ) - idv(u_exact) )*id( v )/hFace() );


        };
    std::cout << "coucou14 " << "\n";
    u=vf::project(Vh->functionSpace<0>(), elements(mesh), zero<2,1>());
    //std::cout << "coucou15 " << "\n";
    p=vf::project(Vh->functionSpace<1>(), elements(mesh), constant(0.0));
    //std::cout << "coucou16 " << "\n";

    backend()->nlSolver()->residual = Residual;
    //std::cout << "coucou17 " << "\n";
    backend()->nlSolver()->jacobian = Jacobian;
    //std::cout << "coucou18 " << "\n";
    backend()->nlSolve( _solution=U );
    //std::cout << "coucou19 " << "\n";

        if ( exporter->doExport() )
            {
                exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
                exporter->step( 0 )->addRegions();
                auto v = U.functionSpace()->template functionSpace<0> ()->element();
                v = U.template element<0>();
                exporter->step( 0 )->add( "u", U.template element<0>() );
                exporter->step( 0 )->add( "p", U.template element<1>() );
                exporter->save();
            }
        //this->exportResults( U, V );
};
}


int main( int argc, char** argv )
{

    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="driven_cavity",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );
    Feel::Driven_Cavity Driven_Cavity;
    Driven_Cavity.run();
}
