/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-04

  Copyright (C) 2009 Christophe Prud'homme
  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)

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
   \file stokes.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#include <feel/feel.hpp>

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description stokesoptions( "Stokes options" );
    stokesoptions.add_options()
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ;
    return stokesoptions.add( Feel::feel_options() ) ;
}




namespace Feel
{
/**
 * \class Stokes class
 * \brief solves the stokes equations
 *
 */
class Stokes
    :
public Simget
{
    typedef Simget super;
public:


    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*mesh*/
    typedef Simplex<2> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    //# marker1 #,
    typedef Lagrange<2, Vectorial> basis_u_type;
    typedef Lagrange<1, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;
    // use lagrange multipliers to ensure zero mean pressure
#if defined( FEELPP_USE_LM )
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
#else
    typedef bases<basis_u_type,basis_p_type> basis_type;
#endif
    //# endmarker1 #

    /*space*/
    //# marker2 #
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //# endmarker2 #

    /* functions */
    //# marker3 #
    typedef space_type::element_type element_type;
    //# endmarker3 #

    /* export */
    typedef Exporter<mesh_type> export_type;

    FEELPP_DONT_INLINE
    Stokes();

    // init mesh and space
    FEELPP_DONT_INLINE
    void init();

    /**
     * run the convergence test
     */
    FEELPP_DONT_INLINE
    void run();

private:


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename ExprUExact, typename ExprPExact>
    FEELPP_DONT_INLINE
    void exportResults( ExprUExact uexact, ExprPExact pexact,
                        element_type& u, element_type& v );

private:

    backend_ptrtype M_backend;
    double meshSize;

    double mu;
    double penalbc;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
}; // Stokes


Stokes::Stokes()
    :
    super(),
    M_backend( backend_type::build( soption("backend") ) ),
    meshSize( this->vm()["hsize"].as<double>() ),
    mu( this->vm()["mu"].as<value_type>() ),
    penalbc( this->vm()["bccoeff"].as<value_type>() )

{

}

void
Stokes::init()
{
    Environment::changeRepository( boost::format( "doc/manual/tutorial/%1%/%2%/P%3%P%4%/h_%5%/" )
                                   % this->about().appName()
                                   % convex_type::name()
                                   % basis_u_type::nOrder % basis_p_type::nOrder
                                   % this->vm()["hsize"].as<double>() );


    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % convex_type().dimension() % 1 ).str() ,
                                         _shape="hypercube",
                                         _dim=convex_type().dimension(),
                                         _h=meshSize ) );


    Xh = space_type::New( mesh );

}
void
Stokes::run()
{
    mpi::timer chrono;
    this->init();
    LOG(INFO) << "chrono init: " << chrono.elapsed() << "\n";

    auto U = Xh->element( "(u,p)" );
    auto V = Xh->element( "(u,q)" );
    auto u = U.element<0>( "u" );
    auto v = V.element<0>( "u" );
    auto p = U.element<1>( "p" );
    auto q = V.element<1>( "p" );
#if defined( FEELPP_USE_LM )
    auto lambda = U.element<2>();
    auto nu = V.element<2>();
#endif
    //# endmarker4 #

    LOG(INFO) << "[dof]         number of dof: " << Xh->nDof() << "\n";
    LOG(INFO) << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    LOG(INFO) << "[dof]      number of dof(U): " << Xh->functionSpace<0>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(U): " << Xh->functionSpace<0>()->nLocalDof()  << "\n";
    LOG(INFO) << "[dof]      number of dof(P): " << Xh->functionSpace<1>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(P): " << Xh->functionSpace<1>()->nLocalDof()  << "\n";

    LOG(INFO) << "Data Summary:\n";
    LOG(INFO) << "   hsize = " << meshSize << "\n";
    LOG(INFO) << "  export = " << this->vm().count( "export" ) << "\n";
    LOG(INFO) << "      mu = " << mu << "\n";
    LOG(INFO) << " bccoeff = " << penalbc << "\n";

    symbol x("x"),y("y");
    auto u1 = -256*y*(y-1)*(2*y-1)*x*x*(x-1)*(x-1);
    auto u2 = 256*x*(x-1)*(2*x-1)*y*y*(y-1)*(y-1);
    matrix u_exact_g = matrix(2,1);
    u_exact_g = u1,u2;
    auto p_exact_g = (x-0.5)*(y-0.5 );

    auto u_exact = expr<2,1,2>( u_exact_g, {x,y} );
    auto p_exact = expr( p_exact_g, {x,y} );
	auto f_g = -mu*laplacian( u_exact_g, {x,y} ) + grad( p_exact_g, {x,y} ).transpose();
    LOG(INFO) << "f = " << f_g << "\n";

    //# marker5 #
    auto deft = gradt( u );
    auto def = grad( v );
    //# endmarker5 #

    //# marker6 #
    // total stress tensor (trial)
    auto SigmaNt = -idt( p )*N()+mu*deft*N();

    // total stress tensor (test)
    auto SigmaN = -id( p )*N()+mu*def*N();
    //# endmarker6 #

    // f is such that f = \Delta u_exact + \nabla p_exact
    auto f = expr<2,1,2>( f_g, {x,y} );

    auto F = M_backend->newVector( Xh );
    auto D =  M_backend->newMatrix( Xh, Xh );

    chrono.restart();
    // right hand side
    auto stokes_rhs = form1( _test=Xh, _vector=F );
    stokes_rhs += integrate( elements( mesh ),inner( f,id( v ) ) );
    stokes_rhs += integrate( boundaryfaces( mesh ), inner( u_exact,-SigmaN+penalbc*id( v )/hFace() ) );

    LOG(INFO) << "chrono lhs: " << chrono.elapsed() << "\n";
    LOG(INFO) << "[stokes] vector local assembly done\n";

    /*
     * Construction of the left hand side
     */
    //# marker7 #
    auto stokes = form2( _test=Xh, _trial=Xh, _matrix=D );

    stokes += integrate( elements( mesh ), mu*inner( deft,def ) );
    LOG(INFO) << "chrono mu*inner(deft,def): " << chrono.elapsed() << "\n";
    chrono.restart();
    stokes +=integrate( elements( mesh ), - div( v )*idt( p ) + divt( u )*id( q ) );
    LOG(INFO) << "chrono (u,p): " << chrono.elapsed() << "\n";
    chrono.restart();
#if defined( FEELPP_USE_LM )
    stokes +=integrate( elements( mesh ), id( q )*idt( lambda ) + idt( p )*id( nu ) );
    LOG(INFO) << "chrono (lambda,p): " << chrono.elapsed() << "\n";
    chrono.restart();
#endif

    stokes +=integrate( boundaryfaces( mesh ), -inner( SigmaNt,id( v ) ) );
    stokes +=integrate( boundaryfaces( mesh ), -inner( SigmaN,idt( u ) ) );
    stokes +=integrate( boundaryfaces( mesh ), +penalbc*inner( idt( u ),id( v ) )/hFace() );
    LOG(INFO) << "chrono bc: " << chrono.elapsed() << "\n";
    chrono.restart();
    //# endmarker7 #


    chrono.restart();
    M_backend->solve( _matrix=D, _solution=U, _rhs=F );
    LOG(INFO) << "chrono solver: " << chrono.elapsed() << "\n";

    this->exportResults( u_exact, p_exact, U, V );
    LOG(INFO) << "chrono export: " << chrono.elapsed() << "\n";

} // Stokes::run



template<typename ExprUExact, typename ExprPExact>
void
Stokes::exportResults( ExprUExact u_exact, ExprPExact p_exact,
                       element_type& U, element_type& V )
{
    auto u = U.element<0>();
    auto p = U.element<1>();

    auto v = V.element<0>();
    auto q = V.element<1>();
#if defined( FEELPP_USE_LM )
    auto lambda = U.element<2>();
    auto nu = V.element<2>();
    LOG(INFO) << "value of the Lagrange multiplier lambda= " << lambda( 0 ) << "\n";

#endif

    double u_errorL2 = normL2( _range=elements( u.mesh() ), _expr=( idv( u )-u_exact ) );
    LOG(INFO) << "||u_error||_2 = " << u_errorL2 << "\n";;

//! [mean]
    double mean_p = mean( _range=elements( u.mesh() ), _expr=idv( p ) )(0,0);
    LOG(INFO) << "[stokes] mean(p)=" << mean_p << "\n";
//! [mean]

//! [norml2]
    double p_errorL2 = normL2( _range=elements( u.mesh() ), _expr=( idv( p )-mean_p - p_exact ) );
    LOG(INFO) << "||p_error||_2 = " << p_errorL2 << "\n";;
//! [norml2]

    double mean_div_u = mean( _range=elements( u.mesh() ), _expr=divv( u ) )(0,0);
    LOG(INFO) << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = normL2( _range=elements( u.mesh() ), _expr=divv( u ) );
    LOG(INFO) << "[stokes] ||div(u)||_2=" << div_u_error_L2 << "\n";

    v = vf::project( u.functionSpace(), elements( u.mesh() ), u_exact );
    q = vf::project( p.functionSpace(), elements( p.mesh() ), p_exact );

    boost::shared_ptr<export_type> exporter( export_type::New() );
    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->addRegions();
        auto v = U.functionSpace()->functionSpace<0> ()->element();
        v = U.element<0>();

#if defined( FEELPP_USE_LM )
        exporter->step( 0 )->add( std::vector<std::string>({"u","p","l"}), U );
        exporter->step( 0 )->add( std::vector<std::string>({"u_exact","p_exact","l_exact"}), V );
#else
        exporter->step( 0 )->add( std::vector<std::string>({"u","p"}), U );
        exporter->step( 0 )->add( std::vector<std::string>({"u_exact","p_exact"}), V );
#endif

        exporter->save();
    }

} // Stokes::export
} // Feel

int
main( int argc, char** argv )
{

    using namespace Feel;

#if defined( FEELPP_USE_LM )
    std::string name = "stokes_lm";
#else
    std::string name = "stokes";
#endif
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name=name.c_str(),
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    // SOME BAD ELEMENTS
    // P1/P0 : locking
    //typedef Feel::Stokes<Simplex<2>, Lagrange<1, Vectorial>,Lagrange<0, Scalar,Discontinuous> > stokes_type;
    // P1/P1 : spurious modes
    //typedef Feel::Stokes<Simplex<2>, Lagrange<1, Vectorial>,Lagrange<1, Scalar> > stokes_type;

    // SOME GOOD ELEMENTS
    // P2/P1
    // CR0/P0
    //typedef Feel::Stokes<Simplex<2>, CrouzeixRaviart<1, Vectorial>,Lagrange<0, Scalar,Discontinuous> > stokes_type;


    /* define and run application */
    Stokes stokes;
    stokes.run();
}





