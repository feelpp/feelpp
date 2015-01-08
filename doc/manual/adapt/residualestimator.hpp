/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Stephane Veys <stephane.veys@gmail.com>
       Date: 2010-08-05

  Copyright (C) 2010 Universite Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@gmail.com>
   \date 2010-08-05
 */
#include <feel/feel.hpp>
#include <feel/feeldiscr/elementdiv.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>

using namespace Feel;
using namespace boost::numeric::ublas;

/**
 * \file residualestimator.html
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline
AboutData
makeAbout()
{
    AboutData about( "residualestimator" ,
                     "residualestimator" ,
                     "0.2",
                     "nD(n=1,2,3) Residual Estimator on Laplacian equation",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "StÃÂÃÂÃÂÃÂÃÂÃÂÃÂÃÂ©phane Veys", "developer", "stephane.veys@gmail.com", "" );
    return about;

}
/**
 * Laplacian Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f \f$ on \f$\Omega\f$ and \f$ u= g \f$ on \f$\Gamma \f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 * \tparam Order the approximation order
 */
template<int Dim, int Order = 1>
class ResidualEstimator
    :
public Simget
{
    typedef Simget super;
public:
    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef boost::shared_ptr<backend_type> backend_ptrtype;


    //! sparse matrix type associated with backend
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    //! sparse matrix type associated with backend (shared_ptr<> type)
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    //! vector type associated with backend
    typedef typename backend_type::vector_type vector_type;
    //! vector type associated with backend (shared_ptr<> type)
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order 1
    typedef Simplex<Dim> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! function space that holds piecewise constant (\f$P_0\f$) functions (e.g. to store material properties or partitioning
    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar,Discontinuous> > > p0_space_type;
    typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
    //! an element type of the \f$P_0\f$ discontinuous function space
    typedef typename p0_space_type::element_type p0_element_type;


    //! the basis type of our approximation space
    typedef bases<Lagrange<1,Scalar> > p1_basis_type;
    typedef FunctionSpace<mesh_type, p1_basis_type> p1_space_type;
    typedef boost::shared_ptr<p1_space_type> p1_space_ptrtype;
    typedef typename p1_space_type::element_type p1_element_type;


    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type,1> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    ResidualEstimator( AboutData const& about )
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        M_backendP1( backend_type::build( soption("backend") ) ),
        meshSize( 0.1 ),
        exporter( export_type::New( "gmsh", this->about().appName() ) ),
        order( 1 ),
        dim( 1 ),
        shape( "hypercube" ),
        fn( 1 ),
        alpha( 3 ),
        beta( 10 ),
        weakdir( 1 ),
        error_type( 1 ),
        tol( 1e-2 ),
        penaldir( 50 )

    {
    }
    ResidualEstimator()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        M_backendP1( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        exporter( export_type::New( "gmsh", this->about().appName() ) ),
        order( ioption("order") ),
        dim( ioption("dim") ),
        shape( soption("shape") ),
        fn( ioption("fn") ),
        alpha( doption("alpha") ),
        beta( doption("beta") ),
        weakdir( ioption("weakdir") ),
        error_type( ioption("adapt-error-type") ),
        tol( doption("adapt-tolerance") ),
        penaldir( doption("penaldir") )

    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

    // this function will move to the mesh library in the mesh class
    BOOST_PARAMETER_CONST_MEMBER_FUNCTION(
        ( mesh_ptrtype ), // return type
        adapt,    // 2. function name

        tag,           // 3. namespace of tag types

        ( required
          ( h, * ) ) // 4. one required parameter, and

        ( optional
          ( maxit,           *( boost::is_integral<mpl::_> ), 10 )
          ( hmin,            *( boost::is_arithmetic<mpl::_> ), 1e-2 )
          ( hmax,            *( boost::is_arithmetic<mpl::_> ), 2 )
          ( model,           *, "" )
          ( statistics,      *( boost::is_integral<mpl::_> ), 0 )
          ( update,          *( boost::is_integral<mpl::_> ), MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER )
          ( collapseOnBoundary, *( boost::is_integral<mpl::_> ), true )
          ( collapseOnBoundaryTolerance, *( boost::is_arithmetic<mpl::_> ), 1e-6 ) ) // 5. optional
    )
    {

    }

private:

    //! linear algebra backend
    backend_ptrtype M_backend;
    backend_ptrtype M_backendP1;

    //! mesh characteristic size
    double meshSize;

    //! dimension
    int dim;

    //! order of finite element approximation
    int order;

    //! shape of the domain
    std::string shape;

    //! function id
    int fn;
    //! parameter
    double alpha;
    double beta;
    bool weakdir;
    double penaldir;
    int error_type;
    double tol;


    //! mesh
    mesh_ptrtype mesh;

    //! exporter
    export_ptrtype exporter;

    //! Piecewise constant functions space
    p0_space_ptrtype P0h;
    p1_space_ptrtype P1h;

    //double estimatorH1, estimatorL2, estimator;
    p1_element_type  h_new;
    std::string msh_name;
    bool first_time;

    int tag_Neumann;
    int tag_Dirichlet;
}; // ResidualEstimator

template<int Dim, int Order>
void
ResidualEstimator<Dim,Order>::run()
{
    if ( dim && dim != Dim ) return ;

    if ( order && order != Order ) return ;

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute ResidualEstimator<" << Dim << ">\n";
    std::vector<double> X;
    X.push_back( meshSize );

    if ( shape == "hypercube" )
        X.push_back( 1 );

    else if ( shape == "ellipsoid" )
        X.push_back( 2 );

    else // default is simplex
        X.push_back( 0 );

    X.push_back( fn );
    X.push_back( alpha );
    X.push_back( beta );
    X.push_back( weakdir );
    X.push_back( penaldir );
    first_time=true;
    X.push_back( first_time );
    X.push_back( error_type );
    X.push_back( tol );
    std::vector<double> Y( 4 );
    run( X.data(), X.size(), Y.data(), Y.size() );

}
template<int Dim, int Order>
void
ResidualEstimator<Dim,Order>::run( const double* X, unsigned long P, double* Y, unsigned long n )
{
    /*
     * set parameters
     */
    meshSize = X[0];

    if ( X[1] == 0 ) shape = "simplex";

    if ( X[1] == 1 ) shape = "hypercube";

    if ( X[1] == 2 ) shape = "ellipsoid";

    fn = X[2];
    alpha = X[3];
    beta = X[4];
    weakdir = X[5];
    penaldir = X[6];
    first_time = X[7];
    error_type = X[8];
    tol = X[9];

    double estimatorH1, estimatorL2, estimator;

    if ( first_time )
    {
        if ( !this->vm().count( "nochdir" ) )
            Environment::changeRepository( boost::format( "doc/tutorial/%1%/%2%-%3%/P%4%/h_%5%/" )
                                           % this->about().appName()
                                           % shape
                                           % Dim
                                           % Order
                                           % meshSize );

        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _parametricnodes=1,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                             _usenames=true,
                                             _shape=shape,
                                             _dim=Dim,
                                             _h=X[0] ),
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
        tag_Neumann = mesh->markerName( "Neumann" );
        tag_Dirichlet = mesh->markerName( "Dirichlet" );

    }//end if(first_time)


    /**
     * The function space and some associated elements(functions) are then defined
     * \snippet residualestimator.hpp toto1
     */
    // [toto1]
    P0h = p0_space_type::New( mesh );
    P1h = p1_space_type::New( mesh );
    space_ptrtype Xh = space_type::New( mesh );
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    // [toto1]



    /** define \f$g\f$ the expression of the exact solution and
     * \f$f\f$ the expression of the right hand side such that \f$g\f$
     * is the exact solution
     * \snippet residualestimator.hpp toto2
     */
    // [toto2]
    //# marker1 #
    value_type pi = M_PI;
    //! deduce from expression the type of g (thanks to keyword 'auto')
    auto fn1 = ( 1-Px()*Px() )*( 1-Py()*Py() )*( 1-Pz()*Pz() )*exp( beta*Px() );
    auto fn2 = sin( alpha*pi*Px() )*cos( alpha*pi*Py() )*cos( alpha*pi*Pz() )*exp( beta*Px() );
    auto g = chi( fn==1 )*fn1 + chi( fn==2 )*fn2;
    auto grad_g =
        trans( chi( fn==1 )*( ( -2*Px()*( 1-Py()*Py() )*( 1-Pz()*Pz() )*exp( beta*Px() )+beta*fn1 )*unitX()+
                              -2*Py()*( 1-Px()*Px() )*( 1-Pz()*Pz() )*exp( beta*Px() )*unitY()+
                              -2*Pz()*( 1-Px()*Px() )*( 1-Py()*Py() )*exp( beta*Px() )*unitZ() )+
               chi( fn==2 )*( +( alpha*pi*cos( alpha*pi*Px() )*cos( alpha*pi*Py() )*cos( alpha*pi*Pz() )*exp( beta*Px() )+beta*fn2 )*unitX()+
                              -alpha*pi*sin( alpha*pi*Px() )*sin( alpha*pi*Py() )*cos( alpha*pi*Pz() )*exp( beta*Px() )*unitY()+
                              -alpha*pi*sin( alpha*pi*Px() )*cos( alpha*pi*Py() )*sin( alpha*pi*Pz() )*exp( beta*Px() )*unitZ() ) );
    //! deduce from expression the type of laplacian  (thanks to keyword 'auto')
    auto minus_laplacian_g =
        ( chi( fn == 1 )*( 2*( 1-Py()*Py() )*( 1-Pz()*Pz() )*exp( beta*Px() ) + 4*beta*Px()*( 1-Py()*Py() )*( 1-Pz()*Pz() )*exp( beta*Px() ) - beta*beta *fn1 +
                           2*( 1-Px()*Px() )*( 1-Pz()*Pz() )*exp( beta*Px() )*chi( Dim >= 2 ) +
                           2*( 1-Px()*Px() )*( 1-Py()*Py() )*exp( beta*Px() )*chi( Dim == 3 ) )  +
          chi( fn == 2 )*  (
              exp( beta*Px() )*( Dim*alpha*alpha*pi*pi*sin( alpha*pi*Px() )-beta*beta*sin( alpha*pi*Px() )-2*beta*alpha*pi*cos( alpha*pi*Px() ) )*
              ( cos( alpha*pi*Py() )*chi( Dim>=2 ) + chi( Dim==1 ) ) * ( cos( alpha*pi*Pz() )*chi( Dim==3 ) + chi( Dim<=2 ) )
          )

        );

    //# endmarker1 #
    // [toto2]


    using namespace Feel::vf;

    /**
     * Construction of the right hand side. F is the vector that holds
     * the algebraic representation of the right habd side of the
     * problem
     * \snippet residualestimator.hpp toto3
     */
    // [toto3]
    //# marker2 #
    vector_ptrtype F( M_backend->newVector( Xh ) );
    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements( mesh ), minus_laplacian_g*id( v ) )+
        integrate( markedfaces( mesh, tag_Neumann ),
                   grad_g*vf::N()*id( v ) );

    //# endmarker2 #
    if ( this->comm().size() != 1 || weakdir )
    {
        //# marker41 #
        form1( _test=Xh, _vector=F ) +=
            integrate( markedfaces( mesh,tag_Dirichlet ), g*( -grad( v )*vf::N()+penaldir*id( v )/hFace() ) );
        //# endmarker41 #
    }

    F->close();

    // [toto3]

    /**
     * create the matrix that will hold the algebraic representation
     * of the left hand side
     * \snippet residualestimator.hpp toto4
     */
    //# marker3 #
    // [toto4]
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );
    // [toto4]

    /**
      * assemble \f$\int_\Omega \nu \nabla u \cdot \nabla v\f$
     * \snippet residualestimator.hpp toto5
    */
    // [toto5]
    form2( _test=Xh, _trial=Xh, _matrix=D, _init=true ) =
        integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) );
    // [toto5]
    //# endmarker3 #

    if ( this->comm().size() != 1 || weakdir )
    {
        /** weak dirichlet conditions treatment for the boundaries marked 1 and 3
         * -# assemble \f$\int_{\partial \Omega} -\nabla u \cdot \mathbf{n} v\f$
         * -# assemble \f$\int_{\partial \Omega} -\nabla v \cdot \mathbf{n} u\f$
         * -# assemble \f$\int_{\partial \Omega} \frac{\gamma}{h} u v\f$
     * \snippet residualestimator.hpp toto6
         */
        // [toto6]
        //# marker10 #
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            integrate( markedfaces( mesh,tag_Dirichlet ),
                       -( gradt( u )*vf::N() )*id( v )
                       -( grad( v )*vf::N() )*idt( u )
                       +penaldir*id( v )*idt( u )/hFace() );
        D->close();
        //# endmarker10 #
        // [toto6]
    }

    else
    {
        /** strong(algebraic) dirichlet conditions treatment for the boundaries marked 1 and 3
         * -# first close the matrix (the matrix must be closed first before any manipulation )
         * -# modify the matrix by cancelling out the rows and columns of D that are associated with the Dirichlet dof
     * \snippet residualestimator.hpp toto7
         */
        // [toto7]
        //# marker5 #
        D->close();
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            on( markedfaces( mesh, tag_Dirichlet ), u, F, g );
        //# endmarker5 #
        // [toto7]

    }


    /**
      * solve the system
     * \snippet residualestimator.hpp toto8
     */
        // [toto8]
    //! solve \f$ D u = F \f$
    backend_type::build( soption("backend") )->solve( _matrix=D, _solution=u, _rhs=F );
        // [toto8]

    /**
   * compute the \f$L_2\f$ norm of the error
     * \snippet residualestimator.hpp toto9
     */
     // [toto9]
    //# marker7 #
    double L2exact = math::sqrt( integrate( elements( mesh ),g*g ).evaluate()( 0,0 ) );
    double L2error2 =integrate( elements( mesh ),( idv( u )-g )*( idv( u )-g ) ).evaluate()( 0,0 );
    double L2error = math::sqrt( L2error2 );
    // [toto9]
    double semiH1error2 = integrate( elements( mesh ),( gradv( u )-grad_g )*( gradv( u )-grad_g ) ).evaluate()( 0,0 );
    double H1error = math::sqrt( L2error2+semiH1error2 );
    double H1exact = math::sqrt( integrate( elements( mesh ),g*g+grad_g*trans( grad_g ) ).evaluate()( 0,0 ) );


    auto RealL2ErrorP0 = integrate( elements( mesh ), ( idv( u )-g )*( idv( u )-g ), _Q<10>() ).broken( P0h );
    auto RealSemiH1ErrorP0 = integrate( elements( mesh ), ( gradv( u )-grad_g )* trans( gradv( u )-grad_g ), _Q<10>() ).broken( P0h );
    auto H1RealErrorP0 = RealL2ErrorP0;
    H1RealErrorP0 += RealSemiH1ErrorP0;
    H1RealErrorP0 = H1RealErrorP0.sqrt();
    /*******************residual estimator**********************/
    //the source terme is given by : minus_laplacian_g


    auto term1 = vf::h()*( minus_laplacian_g+trace( hessv( u ) ) );
    auto term2 = jumpv( gradv( u ) );
    auto term3 = gradv( u )*vf::N()-grad_g*N();


    auto rho = integrate( elements( mesh ), term1*term1, _Q<10>() ).broken( P0h ).sqrt();


    rho += integrate( internalfaces( mesh ),0.25*vf::h()*term2*term2 ).broken( P0h ).sqrt();

    rho += integrate( markedfaces( mesh,tag_Neumann ),
                      vf::h()*term3*term3 ).broken( P0h ).sqrt();


    auto h=vf::project( P0h, elements( mesh ), vf::h() );
    auto npen=vf::project( P0h, elements( mesh ), vf::nPEN() );
    auto H1estimator = rho;

    //if we use real error for adaptation then H1errorP1 will be the projection of H1RealErrorP0 on P1 space
    //else H1errorP1 will be the projection of H1estimator on P1 space
    p1_element_type H1errorP1;

    if ( error_type==2 )
    {
        H1errorP1 = div( vf::sum( P1h, idv( H1RealErrorP0 )*meas() ), vf::sum( P1h, meas() ) );
    }

    else if ( error_type==1 )
    {
        H1errorP1 = div( vf::sum( P1h, idv( H1estimator )*meas() ), vf::sum( P1h, meas() ) );
    }

    else
    {
        std::cout<<"Problem with parameter adapt-error-type, please choice between 1 and 2"<<std::endl;
        return;
    }

    //auto new_hsize=vf::project( P1h, elements(mesh), vf::pow(vf::pow(vf::h(),Order)*(1e-4)/idv(H1estimatorP1),1./Order));
    //auto new_hsize=vf::project( P1h, elements(mesh), vf::h()*(1e-2)/idv(H1estimatorP1) );

    estimatorH1=math::sqrt( H1estimator.pow( 2 ).sum() );
    estimatorL2=math::sqrt( element_product( H1estimator,h ).pow( 2 ).sum() );


    Y[0] = L2error/L2exact;
    Y[1] = H1error/H1exact;
    Y[2] = estimatorL2/L2exact;
    Y[3] = estimatorH1/H1exact;

    h_new = P1h->element();

    h_new = vf::project( P1h, elements( mesh ),
                         vf::max( vf::pow(
                                      vf::pow( vf::h(),Order )*( tol )/idv( H1errorP1 ),
                                      1./Order ),
                                  doption("adapt-hmin") ) );
    /**********************end of residual estimaor*************/


    //! save the results
    //! project the exact solution
    element_type exact_solution( Xh, "exact_solution" );
    exact_solution = vf::project( Xh, elements( mesh ), g );
    element_type u_minus_exact( Xh, "u-g" );
    u_minus_exact = vf::project( Xh, elements( mesh ), idv( u )-g );


    if ( exporter->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( 0 )->setMesh( mesh );
        exporter->step( 0 )->add( "unknown", u );
        exporter->step( 0 )->add( "exact solution", exact_solution );
        exporter->step( 0 )->add( "u - exact solution", u_minus_exact ) ;
        exporter->step( 0 )->add( "nPEN" , npen );
        exporter->step( 0 )->add( "H1 error estimator P0" , H1estimator );
        exporter->step( 0 )->add( "H1 Real error P0" , H1RealErrorP0 );
        exporter->step( 0 )->add( "H1 error P1 (used for determination of new hsize)" , H1errorP1 );
        exporter->step( 0 )->add( "new hsize" , h_new );

        exporter->save();
        LOG(INFO) << "exportResults done\n";
    }


    LOG(INFO)<< " real L2 error : "<<Y[0]<<"\n";
    LOG(INFO)<< " estimated L2 error "<<Y[2]<<"\n";
    LOG(INFO)<< " real H1 error : "<<Y[1]<<"\n";
    LOG(INFO)<< " estimated H1 error "<<Y[3]<<"\n";


    std::ostringstream geostr;

    if ( boption("gmshmodel") )
    {
        if ( boption("gmshgeo") )
            geostr << shape << "-" << Dim << ".geo";

        else
            geostr << shape << "-" << Dim << ".msh";
    }

    mesh  = adapt( _h=h_new,
                   _model=geostr.str(),
                   _hmin=doption("adapt-hmin"),
                   _hmax=doption("adapt-hmax") );
} // ResidualEstimator::run







