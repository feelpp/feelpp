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
#include<feel/feeldiscr/projector.hpp>
#include<feel/feelvf/operations.hpp>
#include <boost/noncopyable.hpp>
#include <boost/signals2.hpp>
#include <boost/format.hpp>

#include <feel/feelpoly/im.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelpde/preconditionerbtpcd.hpp>


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
    typedef bases<basis_u_type,basis_p_type> basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;


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

}; // Steady_Ns

Steady_Ns::Steady_Ns( )
:
super( ),
//M_backend( backend_type::build( this->vm() ) ),
mu( this->vm()["mu"].as<value_type>() ),
rho( this->vm()["rho"].as<value_type>() ),
penalbc( this->vm()["bccoeff"].as<value_type>() )
{

}


void Steady_Ns::init()
{
    double meshSize = option(_name="gmsh.hsize").as<double>();
#if 0
    GeoTool::Rectangle R( meshSize,"myRectangle",GeoTool::Node(0,0),GeoTool::Node(5,1));
    R.setMarker(_type="line",_name="inlet",_marker4=true);
    R.setMarker(_type="line",_name="outlet",_marker2=true);
    R.setMarker(_type="line",_name="wall",_marker1=true,_marker3=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new Mesh<Simplex<2>>,_name="qs_stokes");
#else
    mesh = loadMesh(_mesh=new mesh_type);
#endif


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
        //std::cout<< " In Jacobian function \n";
        auto U = Vh->element( "(u,p)" );
        auto V = Vh->element( "(v,q)" );
        auto u = U.element<0>( "u" );
        auto v = V.element<0>( "u" );
        auto p = U.element<1>( "p" );
        auto q = V.element<1>( "p" );

        auto deft = sym(gradt( u ));
        auto def = sym(grad( v ));

        //std::cout<< "        Start assembling the bilinear form a  \n";
        if (!J) J = backend()->newMatrix( Vh, Vh );
        auto a = form2( _test=Vh, _trial=Vh, _matrix=J );

        a = integrate( elements( mesh ), 2*mu*inner(deft,def) );
        a += integrate( elements( mesh ), - id(q)*divt(u) -idt(p)*div(v) );

        //std::cout<< "        Assembling convectiv terms \n";
        // Convective terms
        a += integrate( elements( mesh ), rho*trans(id(v))*gradv(u)*idt(u));
        a += integrate( elements( mesh ), rho*trans(id(v))*gradt(u)*idv(u));

        //std::cout<< "        Setting boundary condition on the wall \n";
        //(For Passerni test case)Dirichlet Boundary condition on the inlet and wall and stress free on the outlet
        a += integrate( markedfaces(mesh,"wall"),-trans( -idt(p)*N()+2*mu*deft*N() )*id( v ));
        a += integrate( markedfaces(mesh,"wall"),-trans( -id(p)*N()+2*mu*def*N() )*idt( u ));
        a += integrate( markedfaces(mesh,"wall"), +penalbc*trans( idt( u ) )*id( v )/hFace() );

        //std::cout<< "        Setting boundary condition on the inlet \n";
        a += integrate( markedfaces(mesh,"inlet"),-trans( -idt(p)*N()+2*mu*deft*N() )*id( v ));
        a += integrate( markedfaces(mesh,"inlet"),-trans( -id(p)*N()+2*mu*def*N() )*idt( u ));
        a += integrate( markedfaces(mesh,"inlet"), +penalbc*trans( idt( u ) )*id( v )/hFace() );
        //std::cout<< "        DONE \n";


    };

    auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
    {
        //std::cout<< " In residual function \n";
        auto U = Vh->element( "(u,p)" );
        auto V = Vh->element( "(v,q)" );
        auto u = U.element<0>( "u" );
        auto v = V.element<0>( "u" );
        auto p = U.element<1>( "p" );
        auto q = V.element<1>( "p" );

        /////////////////////////////// Use Ginac //////////////////////////////////////
        auto u_exact = expr<2,1>( soption(_name="functions.g") );
        auto f = vec(cst(0.), cst(0.));

        U=*X;
        //std::cout<< "        Start assembling r \n";
        auto r = form1( _test=Vh, _vector=R );
        r = integrate( elements( mesh ),-inner( f,id( v ) ) );
        r += integrate( elements( mesh ), trace(trans(2*mu*gradv( u ))*grad( v )) );
        r +=  integrate( elements( mesh ),-idv(p)*div(v) - id(q)*divv(u));
        //std::cout<< "        Assembling convectiv terms \n";
        // convective terms
        r += integrate( elements( mesh ), trans(rho*gradv( u )*idv(u))*id(v));

        auto SigmaNv = ( -idv( p )*N() + 2*mu*sym(gradv( u ))*N() );
        auto SigmaN = ( -id( q )*N() + 2*mu*sym(grad( v ))*N() );

        //std::cout<< "        Setting boundary condition on the wall \n";
        //Dirichlet Boundary condition on the inlet and wall and stress free on the outlet
        r +=integrate ( markedfaces(mesh,"wall"),  -trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) - vec(cst(0.),cst(0.)) ) + penalbc*trans( idv( u ) -vec(cst(0.),cst(0.)) )*id( v )/hFace() );
        //std::cout<< "        Setting boundary condition on the inlet \n";
        r +=integrate ( markedfaces(mesh,"inlet"),  -trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) - u_exact ) + penalbc*trans( idv( u ) - u_exact )*id( v )/hFace() );
        //std::cout<< "        DONE! \n";


    };


    auto U = Vh->element( "(u,p)" );
    auto V = Vh->element( "(v,q)" );
    auto u = U.element<0>( "u" );
    auto v = V.element<0>( "u" );
    auto p = U.element<1>( "p" );
    auto q = V.element<1>( "p" );

    auto g = expr<2,1>( soption(_name="functions.g") );
    //std::cout<< " Iniliatising p \n";
    p=vf::project(Vh->functionSpace<1>(), elements(mesh), cst(0.));
    //std::cout<< " Iniliatising u \n";
    u=vf::project(Vh->functionSpace<0>(), elements(mesh), vec(cst(0.),cst(0.)));
    

    // -- INITIALIZATION -- //
    //std::cout<< " Iniliatising the residual \n";
    backend()->nlSolver()->residual =Residual;
    //std::cout<< " Iniliatising the Jacobian \n";
    backend()->nlSolver()->jacobian =Jacobian;

    // Solving
    if ( boption("btpcd") )
    {
        std::map<std::string,std::set<flag_type>> bcs;
        bcs["Dirichlet"].insert(mesh->markerName("inlet"));
        bcs["Dirichlet"].insert(mesh->markerName("wall"));
        bcs["Neumann"].insert(mesh->markerName("outlet"));
        auto a_btpcd = btpcd( _space=Vh, _bc=bcs );
        a_btpcd->update( zero<2,1>(), g );
        backend()->nlSolve(_solution=U,_backend=backend(),_prec=a_btpcd );

    }
    else
        backend()->nlSolve( _solution=U );


    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->add( "v", v );
    e->add( "p", p );
    e->save();

};
}


int main( int argc, char** argv )
{
    
    using namespace Feel;
    po::options_description  steadynsoptions( "Navier Stokes problem options" );
    steadynsoptions.add_options()
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "rho", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ;

    Environment env( _argc=argc, _argv=argv,
                    _desc=steadynsoptions,
                    _about=about(_name="steady_ns",
                                 _author="Christophe Prud'homme",
                                 _email="christophe.prudhomme@feelpp.org") );
    Feel::Steady_Ns Steady_Ns;
    Steady_Ns.run();
}
