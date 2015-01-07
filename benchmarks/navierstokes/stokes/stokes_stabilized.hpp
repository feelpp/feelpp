/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2006-06-15

   Copyright (C) 2006 EPFL

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
   \file stokes_stabilized.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2006-06-15
*/

#ifndef __stokesStabilized
#define __stokesStabilized 1

#include <feel/options.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>

#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <feel/feelalg/backend.hpp>

#include <feel/feelvf/vf.hpp>



// tag of the lid
const int LID = 1;
std::string createDomain( double );

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description stokesStabilizedoptions( "StokesStabilized options" );
    stokesStabilizedoptions.add_options()
    ( "nu", Feel::po::value<double>()->default_value( 0.01 ), "viscosity value" )
    ( "hsize", Feel::po::value<double>()->default_value( 1 ), "first h value to start convergence" )
    ( "penalisation", Feel::po::value<double>()->default_value( 1.0 ), "penalisation parameter" )
    ( "export", "export results(ensight, data file(1D)" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "bc" )
    ( "b", Feel::po::value<double>()->default_value( 0 ), "b" )
    ( "matrix", Feel::po::value<bool>()->default_value( 0 ), "print matrix" )
    ( "errors", Feel::po::value<bool>()->default_value( 1 ), "print errors" )
    ( "timings", Feel::po::value<bool>()->default_value( 1 ), "print timings" )
    ( "graphics", Feel::po::value<bool>()->default_value( 0 ), "export graphics" )

    ;
    Feel::po::options_description solveroptions( "algebraic solver options" );
    solveroptions.add_options()
    ( "tolerance", Feel::po::value<double>()->default_value( 1.e-15 ), "tolerance of the iterative solvers" )
    ( "maxiter", Feel::po::value<int>()->default_value( 200 ), "set maximum number of iterations" )
    ( "residual", Feel::po::value<std::string>()->default_value( "" ), "test of residual for GMRES" )
    ( "drop_tolerance", Feel::po::value<double>()->default_value( 1.e-5 ), "drop tolerance for the preconditioner" )
    ( "fillin", Feel::po::value<int>()->default_value( 1 ), "fillin for the preconditioner" )
    ;
    return stokesStabilizedoptions.add( solveroptions ).add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "stokesStabilized" ,
                           "stokesStabilized" ,
                           "0.1",
                           "Stabilized Stokes (interior penalty) solver",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2007 Universite Joseph Fourier" );

    about.addAuthor( "Goncalo Pena", "developer", "goncalo.pena@epfl.ch", "" );
    return about;

}


namespace Feel
{
template< int Order >
class StokesStabilized
    :
public Application
{
    typedef Application super;

public:

    // -- TYPEDEFS --
    static const uint16_type uOrder = Order;
    static const uint16_type pOrder = uOrder;

    static const uint16_type Dim = 2;

    typedef double value_type;

    typedef Simplex<Dim, 1,Dim> entity_type;

    /*matrix*/
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::vector_type vector_type;

    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;


    /*mesh*/
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef fusion::vector<fem::Lagrange<Dim, uOrder, Vectorial, Continuous, double, Simplex, PointSetWarpBlend>,
            fem::Lagrange<Dim, pOrder, Scalar, Continuous, double, Simplex, PointSetWarpBlend>
            > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_0_type;
    typedef typename element_type::template sub_element<1>::type element_1_type;

    typedef typename space_type::template sub_functionspace<0>::type velocity_space_ptrtype;
    typedef typename space_type::template sub_functionspace<1>::type pressure_space_ptrtype;

    typedef typename velocity_space_ptrtype::element_type velocity_space_type;
    typedef typename pressure_space_ptrtype::element_type pressure_space_type;

    typedef typename velocity_space_type::element_type velocity_element_type;
    typedef typename pressure_space_type::element_type pressure_element_type;

    /*quadrature*/
    template<int i>
    struct MyIm
    {
        typedef IM<Dim, i, value_type, Simplex> type;
    };

    /* export */
    typedef Exporter<mesh_type> export_type;

    StokesStabilized( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        meshSize( doption("hsize") ),
        nu( doption("nu") ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
    {
        VLOG(1) << "[StokesStabilized] hsize = " << meshSize << "\n";
        VLOG(1) << "[StokesStabilized] nu = " << nu << "\n";
        VLOG(1) << "[StokesStabilized] export = " << this->vm().count( "export" ) << "\n";

    }


    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize );

    /**
     * run the convergence test
     */
    void run( );

private:

    void exportResults( double time, element_type& u );

private:

    backend_ptrtype M_backend;

    double meshSize;
    double nu;

    boost::shared_ptr<export_type> exporter;

}; // StokesStabilized



template<int Order>
typename StokesStabilized<Order>::mesh_ptrtype
StokesStabilized<Order>::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );

    GmshHypercubeDomain<2,1,2,Simplex> td;
    td.setCharacteristicLength( meshSize );

    td.setY( std::make_pair( -0.5, 1.0 ) );
    td.setX( std::make_pair( -0.5, 1.5 ) );

    std::string meshName = td.generate( "square" );

    ImporterGmsh<mesh_type> import( meshName );

    mesh->accept( import );

    return mesh;
} // StokesStabilized::createMesh


template<int Order>
void
StokesStabilized<Order>::run( )
{
    this->changeRepository( boost::format( "benchmarks/%1%/%2%D/P%3%P%4%/h_%5%/" )
                            % this->about().appName()
                            % Dim
                            % uOrder % pOrder
                            % doption("hsize")
                          );
    /*
     * logs will be in <feel repo>/<app name>/<entity>/P<p>/h_<h>
     */
    this->setLogs();

    using namespace Feel::vf;

    mesh_ptrtype mesh = createMesh( meshSize );

    space_ptrtype Xh = space_type::New( mesh );

    element_type U( Xh, "U" );
    element_type V( Xh, "V" );

    element_0_type u( U.template element<0>() );
    element_0_type v( U.template element<0>() );
    element_1_type p( U.template element<1>() );
    element_1_type q( U.template element<1>() );

    U.container() = ublas::scalar_vector<double>( U.size(), 0.0 );
    /*
     * a quadrature rule for numerical integration
     */

    std::cout << "Number of dofs: " << Xh->nDof() << "\n";

    double pi = 4.0 * math::atan( double( 1.0 ) );
    const double penalisation = doption("penalisation");
    const double tolerance = doption("tolerance");
    const int maxiter = ioption("maxiter");
    const int fillin = ioption("fillin");
    const int bctype = ioption("bctype");
    const double b = doption("b");
    const bool printMatrix = boption("matrix");

    const double drop_tolerance = doption("drop_tolerance");
    const std::string residual = soption("residual");
    double lambda = 1./( 2.*nu ) - math::sqrt( 1./( 4.*nu*nu ) + 4.*pi*pi );

    const bool printErrors = boption("errors");
    const bool printTimings = boption("timings");
    const bool exportGraphics = boption("graphics");

    AUTO( u1, val( 1. - exp( lambda * Px() ) * cos( 2.*pi*Py() ) ) );
    AUTO( u2, val( ( lambda/( 2.*pi ) ) * exp( lambda * Px() ) * sin( 2.*pi*Py() ) ) );

    AUTO( u_exact, u1*oneX() + u2*oneY() );

    AUTO( du_dx, val( -lambda*exp( lambda * Px() )*cos( 2.*pi*Py() ) ) );
    AUTO( du_dy, val( 2*pi*exp( lambda * Px() )*sin( 2.*pi*Py() ) ) );
    AUTO( dv_dx, val( ( lambda*lambda/( 2*pi ) )*exp( lambda * Px() )*sin( 2.*pi*Py() ) ) );
    AUTO( dv_dy, val( lambda*exp( lambda * Px() )*cos( 2.*pi*Py() ) ) );

    AUTO( grad_exact, ( mat<2,2>( du_dx, du_dy, dv_dx, dv_dy ) ) );

    AUTO( beta, b*( oneX() + oneY() ) );

    AUTO( convection, grad_exact*beta );

    AUTO( p_exact, val( ( 1-exp( 2.*lambda*Px() ) )/2.0 ) );

    AUTO( f1, val( exp( lambda * Px() )*( ( lambda*lambda - 4.*pi*pi )*nu*cos( 2.*pi*Py() ) - lambda*exp( lambda * Px() ) ) ) );
    AUTO( f2, val( exp( lambda * Px() )*nu*( lambda/( 2.*pi ) )*sin( 2.*pi*Py() )*( -lambda*lambda +4*pi*pi ) ) );

    AUTO( f, f1*oneX() + f2*oneY() + convection );

    size_type pattern = Pattern::COUPLED|Pattern::EXTENDED;

    mpi::timer timer;

    timer.restart();
    sparse_matrix_ptrtype Dcst( M_backend->newMatrix( Xh,Xh ) );


    form2( Xh, Xh, Dcst, _init=true, _pattern = pattern ) = integrate( elements( mesh ), typename MyIm<2*uOrder-1>::type(),
            ( nu*( trace( trans( gradt( u ) )*grad( v ) ) ) + trans( gradt( u )*beta )*id( v )
              - div( v ) * idt( p )
              + divt( u ) * id( q ) )
                                                                     );


    double aux = timer.elapsed();

    double p_term = double( uOrder );
    p_term = math::pow( p_term, 7./2. );
    double p_term2 = math::pow( double( uOrder ), 1./2. );

    timer.restart();

    if ( b == 0 )
        form2( Xh, Xh, Dcst ) += integrate( internalfaces( mesh ), typename MyIm<2*( pOrder-1 )>::type(),
                                            constant( penalisation )*hFace()*hFace()
                                            *( trans( jumpt( gradt( p ) ) )*jump( grad( q ) ) )
                                            /p_term
                                          );

    else
        form2( Xh, Xh, Dcst ) += integrate( internalfaces( mesh ), typename MyIm<2*( pOrder-1 )>::type(),
                                            constant( penalisation )*hFace()*hFace()*vf::min( 1,b*hFace()/( nu*p_term2 ) )
                                            *( trans( jumpt( gradt( p ) ) )*jump( grad( q ) ) )
                                            /( p_term*b )
                                          );

    if ( printTimings == 1 )
        std::cout << "Stabilisation term assembly: " << timer.elapsed() << "\n";

    timer.restart();
    form2( Xh, Xh, Dcst ) += integrate( boundaryfaces( mesh ), typename MyIm<uOrder+pOrder>::type(),
                                        trans( idt( p )*N() )*id( v )
                                      );

    form2( Xh, Xh, Dcst ) += integrate( boundaryfaces( mesh ), typename MyIm<2*uOrder-1>::type(),
                                        - nu*trans( gradt( u )*N() )*id( v )
                                      );
    aux += timer.elapsed();

    if ( printTimings == 1 )
        std::cout << "Stokes matrix assembly: " << aux << "\n";

    timer.restart();
    Dcst->close();

    if ( printTimings == 1 )
        std::cout << "Stokes matrix global assembly: " << timer.elapsed() << "\n";


    if ( printMatrix )
        Dcst->printMatlab( "Dcst.m" );

    timer.restart();
    vector_ptrtype Fcst( M_backend->newVector( Xh ) );
    form1( Xh, Fcst, _init=true ) = integrate( elements( mesh ), typename MyIm<20>::type(),
                                    trans( f )*id( v )
                                             );

    Fcst->close();

    if ( printTimings == 1 )
        std::cout << "RHS assembly: " << timer.elapsed() << "\n";

    timer.restart();
    form2( Xh, Xh, Dcst ) += on( boundaryfaces( mesh ), U.template element<0>(), Fcst, u_exact );

    if ( printTimings == 1 )
        std::cout << "Boundary conditions: " << timer.elapsed() << "\n";


    if ( printTimings == 0 )
    {

        backend_ptrtype ns( backend_type::build( this->vm() ) );
#if 0
        ns.set_residual( residual );
        ns.set_verbose( "all" );
        ns.set_tol( tolerance );
        ns.set_maxiter( maxiter );
        ns.set_fillin( fillin );
        ns.set_drop( drop_tolerance );


        PreconditionerIfpack ns_prec( ns.get_options(), "ILU" );
        ns_prec.buildPreconditioner( *Dcst );

        op_mat_ptrtype nsOp = new op_mat_type( Dcst, ns.get_options(), "NS", ns_prec.getPrec() );
        vector_ptrtype X( M_backend->newVector( Xh ) );
        X->zero();
        X->close();

        std::cout << "Solve linear system...\n";
        nsOp->ApplyInverse( Fcst->vec(), X->vec() );

        backend_type::Epetra2Ublas( *X, U );

#else
        vector_ptrtype X( M_backend->newVector( Xh ) );
        X->zero();
        X->close();
        ns->solve( Dcst, Dcst, X, Fcst );
        U = *X;

#endif


    }

    double area = 0;
    double p_l1_error = 0;
    double mean_value = 0;
    double p_error = 0;
    double H1_u_error = 0;
    double divergence = 0;


    if ( printTimings == 0 )
    {
        area = integrate( elements( mesh ), typename MyIm<1>::type(),
                          constant( 1.0 ) ).evaluate()( 0,0 );

        p_l1_error = integrate( elements( mesh ), typename MyIm<8*pOrder>::type(),
                                idv( U.template element<1>() ) - p_exact ).evaluate()( 0,0 );

        mean_value = p_l1_error / area;

        p_error = integrate( elements( mesh ), typename MyIm<8*pOrder>::type(),
                             ( idv( U.template element<1>() ) - p_exact - constant( mean_value ) )^2
                           ).evaluate()( 0,0 );
    }

    timer.restart();
    double L2_u_error = integrate( elements( mesh ), typename MyIm<8*uOrder>::type(),
                                   trans( idv( U.template element<0>() )-u_exact )
                                   *( idv( U.template element<0>() )-u_exact )
                                 ).evaluate()( 0,0 );

    if ( printTimings )
        std::cout << "Timing L2 error for the velocity: " << timer.elapsed() << "\n";

    if ( printTimings == 0 )
    {
        H1_u_error = integrate( elements( mesh ), typename MyIm<8*uOrder>::type(),
                                trace( ( gradv( U.template element<0>() )-grad_exact )
                                       *trans( gradv( U.template element<0>() )-grad_exact ) )
                              ).evaluate()( 0,0 );


        divergence = integrate( elements( mesh ), typename MyIm<2*( uOrder-1 )>::type(),
                                trace( gradv( U.template element<0>() ) )
                                *trace( gradv( U.template element<0>() ) )
                              ).evaluate()( 0,0 );

    }

    if ( exportGraphics )
    {
        this->exportResults( 0, U );

        V.template element<0>() = project( Xh->template functionSpace<0>(), elements( mesh ), u_exact );
        V.template element<1>() = project( Xh->template functionSpace<1>(), elements( mesh ), p_exact );


        this->exportResults( 1, V );
    }

    if ( printErrors )
    {
        std::cout << "L2 error for the pressure: " << math::sqrt( p_error ) << "\n";
        std::cout << "L2 error for the velocity: " << math::sqrt( L2_u_error ) << "\n";
        std::cout << "H1 error for the velocity: " << math::sqrt( H1_u_error + L2_u_error ) << "\n";
        std::cout << "divergence for the velocity: " << math::sqrt( divergence ) << "\n";
    }


} // StokesStabilized::run


template<int Order>
void
StokesStabilized<Order>::exportResults( double time,
                                        element_type& U )
{
    exporter->step( time )->setMesh( U.functionSpace()->mesh() );
    exporter->step( time )->add( "u", U.template element<0>().comp( X ) );
    exporter->step( time )->add( "v", U.template element<0>().comp( Y ) );
    exporter->step( time )->add( "pressure", U.template element<1>() );
    exporter->save();
} // StokesStabilized::export

// instantiation
template<int Order> const uint16_type StokesStabilized<Order>::uOrder;
template<int Order> const uint16_type StokesStabilized<Order>::pOrder;
template<int Order> const uint16_type StokesStabilized<Order>::Dim;




} // Feel

#endif // __stokesStabilized
