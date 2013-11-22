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
         ( "continuation", Feel::po::value<bool>()->default_value( 0 ), "use continuation algorithm: =0 (no) =1(yes)" )
         ( "Re", Feel::po::value<double>()->default_value( 100 ), "Reynolds number" )
         ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
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

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Ranine Tarabay", "developer", "ranine.a.tarabay@gmail.com", "" );
    return about;

}


namespace Feel
{
template<int Dim>
class DrivenCavity : public Application
{
    typedef Application super;
public:


    typedef double value_type;

    typedef Simplex<Dim> convex_type;
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

    DrivenCavity( );

    // init mesh and space
    FEELPP_DONT_INLINE void init();

    FEELPP_DONT_INLINE void run();

    FEELPP_DONT_INLINE void exportResults( element_type const& U );

    FEELPP_DONT_INLINE void Jacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J);
    FEELPP_DONT_INLINE void Residual(const vector_ptrtype& X, vector_ptrtype& R);

private:

    double meshSize;

    double Re;
    double penalbc;

    mesh_ptrtype mesh;
    space_ptrtype Vh;

    boost::shared_ptr<export_type> exporter;
}; // DrivenCavity

template<int Dim>
DrivenCavity<Dim>::DrivenCavity( )
    :
    super( ),
    Re( this->vm()["Re"].template as<value_type>() ),
    penalbc( this->vm()["bccoeff"].template as<value_type>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{

}

template<int Dim>
void DrivenCavity<Dim>::init()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }
    mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    std::cout << "number of elements of 2D: " << mesh->numElements() << "\n";

    Vh = space_type::New( mesh );
}

template<int Dim>
void
DrivenCavity<Dim>::Jacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J)
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

    //#endif

    if (!J) J = backend(_name="newtonns")->newMatrix( Vh, Vh );
    auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
    a = integrate( elements( mesh ), inner(gradt( u ),grad( v ))/Re );
    a += integrate( elements( mesh ), id(q)*divt(u) -idt(p)*div(v) );
    // Convective terms
    a += integrate( elements( mesh ), trans(id(v))*gradv(u)*idt(u));
    a += integrate( elements( mesh ), trans(id(v))*gradt(u)*idv(u));

    //#if defined( FEELPP_USE_LM )
    a += integrate(elements(mesh), id(q)*idt(lambda)+idt(p)*id(nu));
    //#elif
    //a += integrate(elements(mesh), idt(p)*id(nu));

    //Weak Dirichlet conditions
    a += integrate( boundaryfaces( mesh ),-trans( -idt(p)*N()+gradt(u)*N()/Re )*id( v ));//
    a += integrate( boundaryfaces( mesh ),-trans( -id(p)*N()+grad(u)*N()/Re )*idt( u ));//
    a += integrate( boundaryfaces( mesh ), +penalbc*inner( idt( u ),id( v ) )/hFace() );


}

template<int Dim>
void
DrivenCavity<Dim>::Residual(const vector_ptrtype& X, vector_ptrtype& R)
{
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

    auto uex=unitX();
    auto u_exact=vf::project(Vh->template functionSpace<0>(), markedfaces(mesh, "wall2"), uex );


    U=*X;
    auto r = form1( _test=Vh, _vector=R );
    //r += integrate( elements( mesh ),-inner( f,id( v ) ) );
    r = integrate( elements( mesh ), trans(gradv( u )*idv(u))*id(v));//convective term
    r += integrate( elements( mesh ), inner(gradv( u ),grad( v ))/Re );
    r +=  integrate( elements( mesh ),-idv(p)*div(v) + id(q)*divv(u));
    //#if defined( FEELPP_USE_LM )
    r += integrate ( elements( mesh ), +id( q )*idv( lambda )+idv( p )*id( nu ) );
    //#endif


    //Weak Dirichlet
    auto SigmaNv = ( -idv( p )*N() + gradv( u )*N()/Re );
    auto SigmaN = ( -id( q )*N() + grad( v )*N()/Re );
    r +=integrate ( boundaryfaces(mesh), - trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) -idv(u_exact) ) + penalbc*trans( idv( u ) - idv(u_exact) )*id( v )/hFace() );


}

template<int Dim>
void DrivenCavity<Dim>::exportResults( element_type const& U )
{
    auto uex=unitX();
    auto u_exact=vf::project(Vh->template functionSpace<0>(), markedfaces(mesh, "wall2"), uex );

    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->addRegions();
        exporter->step( 0 )->add( "u", U.template element<0>() );
        exporter->step( 0 )->add( "p", U.template element<1>() );
        exporter->step( 0 )->add( "uex", u_exact);
        exporter->save();
    }

}
template<int Dim>
void DrivenCavity<Dim>::run()
{
    this->init();

    auto U = Vh->element( "(u,p)" );
    auto V = Vh->element( "(v,q)" );
    auto u = U.template element<0>( "u" );
    auto v = V.template element<0>( "u" );
    auto p = U.template element<1>( "p" );
    auto q = V.template element<1>( "p" );
    //#if defined( FEELPP_USE_LM )
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
    //#endif



    u=vf::project(Vh->template functionSpace<0>(), elements(mesh), zero<Dim,1>());
    p=vf::project(Vh->template functionSpace<1>(), elements(mesh), constant(0.0));

    std::cout << "Initializing residual " << "\n";
    backend(_name="newtonns")->nlSolver()->residual = boost::bind( &DrivenCavity::Residual,
                                                                   boost::ref( *this ), _1, _2 );
    std::cout << "Initializing the Jacobian matrix" << "\n";
    backend(_name="newtonns")->nlSolver()->jacobian = boost::bind( &DrivenCavity::Jacobian,
                                                                   boost::ref( *this ), _1, _2 );

    if ( option(_name="continuation" ).as<bool>() )
    {
        double ReTarget = Re;
        int N = std::ceil( std::log( Re ) );
        for( int i  = 0; i <= N; ++i )
        {
            Re = std::exp( std::log(1)+i*(std::log(ReTarget)-std::log(1))/double(N));
            std::cout << "Start solving for Reynolds = " << Re << "\n";
            backend(_name="newtonns")->nlSolve( _solution=U );
            std::cout << "Done!" << "\n";
        }
    }
    else
    {
        std::cout << "Start solving for Reynolds = " << Re << "\n";
        backend(_name="newtonns")->nlSolve( _solution=U );
        std::cout << "Done!" << "\n";
    }

    this->exportResults( U );
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
    Feel::DrivenCavity<2> DrivenCavity;
    DrivenCavity.run();
}
