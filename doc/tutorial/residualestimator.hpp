/* -*- mode: c++; coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
             Stéphane Veys <stephane.veys@gmail.com>
       Date: 2010-08-05

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file residualestimator.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \author Stéphane Veys <stephane.veys@gmail.com>
   \date 2010-08-05
 */
/** include predefined life command line options */
#include <life/options.hpp>

/** include linear algebra backend */
#include <life/lifealg/backend.hpp>

/** include function space class */
#include <life/lifediscr/functionspace.hpp>

/** include helper function to define \f$P_0\f$ functions associated with regions  */
#include <life/lifediscr/region.hpp>

/** include integration methods */
#include <life/lifepoly/im.hpp>

/** include gmsh mesh importer */
#include <life/lifefilters/gmsh.hpp>

/** include exporter factory class */
#include <life/lifefilters/exporter.hpp>

/** include  polynomialset header */
#include <life/lifepoly/polynomialset.hpp>

/** include  the header for the variational formulation language (vf) aka FEEL++ */
#include <life/lifevf/vf.hpp>

#if defined (HAVE_MADLIB_H)
#include <MAdLib.h>
#endif // HAVE_MADLIB_H



/** use Life namespace */
using namespace Life;
using namespace Life::vf;
using namespace boost::numeric::ublas;

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Life::Application subclass.
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
                     Life::AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    about.addAuthor("Stéphane Veys", "developer", "stephane.veys@gmail.com", "");
    return about;

}
/**
 * \class ResidualEstimator
 *
 * Laplacian Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
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
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    ResidualEstimator( AboutData const& about )
        :
        super( about ),
        M_backend( backend_type::build() ),
        M_backendP1( backend_type::build() ),
        meshSize( 0.1 ),
        exporter( Exporter<mesh_type>::New( "gmsh", this->about().appName() ) ),
        order( 1 ),
        dim( 1 ),
        shape( "hypercube" ),
        fn( 1 ),
        alpha( 3 ),
        beta( 10 ),
        weakdir( 1 ),
        tol( 1e-2 ),
        penaldir( 50 )

    {
    }
    ResidualEstimator( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        M_backend( backend_type::build( this->vm() ) ),
        M_backendP1( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( "gmsh", this->about().appName() ) ),
        order( this->vm()["order"].template as<int>() ),
        dim( this->vm()["dim"].template as<int>() ),
        shape( this->vm()["shape"].template as<std::string>() ),
        fn( this->vm()["fn"].template as<int>() ),
        alpha( this->vm()["alpha"].template as<double>() ),
        beta( this->vm()["beta"].template as<double>() ),
        weakdir( this->vm()["weakdir"].template as<int>() ),
        tol( this->vm()["tol"].template as<double>() ),
        penaldir( this->vm()["penaldir"].template as<double>() )

    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N);

    mesh_ptrtype adapt_mesh( p1_element_type& sizefield );
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
    double tol;

    //! mesh
    mesh_ptrtype mesh;

    //! exporter
    export_ptrtype exporter;

    //! Piecewise constant functions space
    p0_space_ptrtype P0h;
    p1_space_ptrtype P1h;

    double estimatorH1, estimatorL2, estimator;
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
    std::vector<double> Y( 4 );
    first_time=true;
    run( X.data(), X.size(), Y.data(), Y.size() );

    mesh  = adapt_mesh( h_new );
    first_time=false;
    run( X.data(), X.size(), Y.data(), Y.size() );


}
template<int Dim, int Order>
void
ResidualEstimator<Dim,Order>::run( const double* X, unsigned long P, double* Y, unsigned long n)
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





    if(first_time){
      if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "doc/tutorial/%1%/%2%-%3%/P%4%/h_%5%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );

        mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=(boost::format( "%1%-%2%" ) % shape % Dim).str() ,
                                                      _usenames=true,
                                                      _shape=shape,
                                                      _dim=Dim,
                                                      _h=X[0] ),
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
        tag_Neumann = mesh->markerName("Neumann");
        tag_Dirichlet = mesh->markerName("Dirichlet");

    }//end if(first_time)
    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    P0h = p0_space_type::New( mesh );
    P1h = p1_space_type::New( mesh );
    space_ptrtype Xh = space_type::New( mesh );
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    /** \endcode */

    /** define \f$g\f$ the expression of the exact solution and
     * \f$f\f$ the expression of the right hand side such that \f$g\f$
     * is the exact solution
     */
    /** \code */
    //# marker1 #
    value_type pi = M_PI;
    //! deduce from expression the type of g (thanks to keyword 'auto')
    auto fn1 = (1-Px()*Px())*(1-Py()*Py())*(1-Pz()*Pz())*exp(beta*Px());
    auto fn2 = sin(alpha*pi*Px())*cos(alpha*pi*Py())*cos(alpha*pi*Pz())*exp(beta*Px());
    auto g = chi(fn==1)*fn1 + chi(fn==2)*fn2;
    auto grad_g =
        trans(chi(fn==1)*((-2*Px()*(1-Py()*Py())*(1-Pz()*Pz())*exp(beta*Px())+beta*fn1)*unitX()+
                          -2*Py()*(1-Px()*Px())*(1-Pz()*Pz())*exp(beta*Px())*unitY()+
                          -2*Pz()*(1-Px()*Px())*(1-Py()*Py())*exp(beta*Px())*unitZ())+
              chi(fn==2)*(+(alpha*pi*cos(alpha*pi*Px())*cos(alpha*pi*Py())*cos(alpha*pi*Pz())*exp(beta*Px())+beta*fn2)*unitX()+
                          -alpha*pi*sin(alpha*pi*Px())*sin(alpha*pi*Py())*cos(alpha*pi*Pz())*exp(beta*Px())*unitY()+
                          -alpha*pi*sin(alpha*pi*Px())*cos(alpha*pi*Py())*sin(alpha*pi*Pz())*exp(beta*Px())*unitZ()));
    //! deduce from expression the type of laplacian  (thanks to keyword 'auto')
    auto minus_laplacian_g =
        (chi( fn == 1 )*( 2*(1-Py()*Py())*(1-Pz()*Pz())*exp(beta*Px()) + 4*beta*Px()*(1-Py()*Py())*(1-Pz()*Pz())*exp(beta*Px()) - beta*beta *fn1 +
                          2*(1-Px()*Px())*(1-Pz()*Pz())*exp(beta*Px())*chi(Dim >= 2) +
                          2*(1-Px()*Px())*(1-Py()*Py())*exp(beta*Px())*chi(Dim == 3) )  +
         chi( fn == 2 )*  (
             exp(beta*Px())*(Dim*alpha*alpha*pi*pi*sin(alpha*pi*Px())-beta*beta*sin(alpha*pi*Px())-2*beta*alpha*pi*cos(alpha*pi*Px()))*
             ( cos(alpha*pi*Py())*chi(Dim>=2) + chi(Dim==1)) * ( cos(alpha*pi*Pz())*chi(Dim==3) + chi(Dim<=2) )
             )

	 );

    //# endmarker1 #
    /** \endcode */




    using namespace Life::vf;

    /**
     * Construction of the right hand side. F is the vector that holds
     * the algebraic representation of the right habd side of the
     * problem
     */
    /** \code */
    //# marker2 #
    vector_ptrtype F( M_backend->newVector( Xh ) );
    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements(mesh), minus_laplacian_g*id(v) )+
        integrate( markedfaces( mesh, tag_Neumann ),
                   grad_g*vf::N()*id(v) );
    //# endmarker2 #
    if ( this->comm().size() != 1 || weakdir )
    {
        //# marker41 #
        form1( _test=Xh, _vector=F ) +=
            integrate( markedfaces(mesh,tag_Dirichlet), g*(-grad(v)*vf::N()+penaldir*id(v)/hFace()) );
        //# endmarker41 #
    }
    F->close();

    /** \endcode */

    /**
     * create the matrix that will hold the algebraic representation
     * of the left hand side
     */
    //# marker3 #
    /** \code */
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );
    /** \endcode */

    //! assemble $\int_\Omega \nu \nabla u \cdot \nabla v$
    /** \code */
    form2( Xh, Xh, D, _init=true ) =
        integrate( elements(mesh), gradt(u)*trans(grad(v)) );
    /** \endcode */
    //# endmarker3 #

    if ( this->comm().size() != 1 || weakdir )
    {
        /** weak dirichlet conditions treatment for the boundaries marked 1 and 3
         * -# assemble \f$\int_{\partial \Omega} -\nabla u \cdot \mathbf{n} v\f$
         * -# assemble \f$\int_{\partial \Omega} -\nabla v \cdot \mathbf{n} u\f$
         * -# assemble \f$\int_{\partial \Omega} \frac{\gamma}{h} u v\f$
         */
        /** \code */
        //# marker10 #
        form2( Xh, Xh, D ) +=
            integrate( markedfaces(mesh,tag_Dirichlet),
                       -(gradt(u)*vf::N())*id(v)
                       -(grad(v)*vf::N())*idt(u)
                       +penaldir*id(v)*idt(u)/hFace());
        D->close();
        //# endmarker10 #
        /** \endcode */
    }
    else
    {
        /** strong(algebraic) dirichlet conditions treatment for the boundaries marked 1 and 3
         * -# first close the matrix (the matrix must be closed first before any manipulation )
         * -# modify the matrix by cancelling out the rows and columns of D that are associated with the Dirichlet dof
         */
        /** \code */
        //# marker5 #
        D->close();
        form2( Xh, Xh, D ) +=
            on( markedfaces(mesh, tag_Dirichlet), u, F, g );
        //# endmarker5 #
        /** \endcode */

    }
    /** \endcode */


    //! solve the system
    /** \code */
    //! solve \f$ D u = F \f$
    backend_type::build()->solve( _matrix=D, _solution=u, _rhs=F );
    /** \endcode */

    //! compute the \f$L_2$ norm of the error
    /** \code */
    //# marker7 #
    double L2exact = math::sqrt(integrate(elements(mesh),g*g ).evaluate()(0,0));
    double L2error2 =integrate(elements(mesh),(idv(u)-g)*(idv(u)-g) ).evaluate()(0,0);
    double L2error =   math::sqrt( L2error2 );
    double semiH1error2 = integrate(elements(mesh),(gradv(u)-grad_g)*(gradv(u)-grad_g)).evaluate()(0,0);
    double H1error = math::sqrt(L2error2+semiH1error2);
    double H1exact = math::sqrt(integrate(elements(mesh),g*g+grad_g*trans(grad_g) ).evaluate()(0,0));


    /*******************residual estimator**********************/
    //the source terme is given by : minus_laplacian_g

    auto term1 = vf::h()*(minus_laplacian_g+trace(hessv(u)));
    auto term2 = jumpv(gradv(u));
    auto term3 = gradv(u)*vf::N()-grad_g*N();

    auto rho = integrate(elements(mesh), term1*term1 ).broken(P0h).sqrt();

    rho += integrate(internalfaces(mesh),0.25*vf::h()*term2*term2).broken( P0h ).sqrt();

    rho += integrate( markedfaces(mesh,tag_Neumann),
                      vf::h()*term3*term3 ).broken( P0h ).sqrt();

    auto h=vf::project(P0h, elements(mesh), vf::h() );
    auto npen=vf::project(P0h, elements(mesh), vf::nPEN() );
    auto H1estimator = element_product( rho, npen.sqrt() );

    auto H1estimatorP1 = element_div( vf::sum( P1h, idv(H1estimator)*meas()), vf::sum( P1h, meas()) );
    //auto new_hsize=vf::project( P1h, elements(mesh), vf::pow(vf::pow(vf::h(),Order)*(1e-4)/idv(H1estimatorP1),1./Order));
    //auto new_hsize=vf::project( P1h, elements(mesh), vf::h()*(1e-2)/idv(H1estimatorP1) );

    double estimatorH1=math::sqrt(H1estimator.pow(2).sum());
    double estimatorL2=math::sqrt(element_product(H1estimator,h).pow(2).sum());


    Y[0] = L2error/L2exact;
    Y[1] = H1error/H1exact;
    Y[2] = estimatorL2/L2exact;
    Y[3] = estimatorH1/H1exact;

    h_new = P1h->element();
    //value_type tol = this->vm()["tol"].as<double>();
    h_new = vf::project( P1h, elements(mesh),
                         vf::max(vf::pow(
                                     vf::pow( vf::h(),Order)*(tol)/idv(H1estimatorP1),
                                     1./Order),
                                 this->vm()["hmin"].template as<double>() ) );
    /**********************end of residual estimaor*************/


    //! save the results
    /** \code */
    //! project the exact solution
    element_type exact_solution( Xh, "exact_solution" );
    exact_solution = vf::project( Xh, elements(mesh), g );
    element_type u_minus_exact( Xh, "u-g" );
    u_minus_exact = vf::project( Xh, elements(mesh), idv(u)-g );


    if ( exporter->doExport() )
    {
        Log() << "exportResults starts\n";

        exporter->step(0)->setMesh( mesh );
        exporter->step(0)->add( "unknown", u );
        exporter->step(0)->add( "exact solution", exact_solution);
        exporter->step(0)->add( "u - exact solution", u_minus_exact) ;
        exporter->step(0)->add( "nPEN" , npen );
        exporter->step(0)->add( "H1 error estimator P0" , H1estimator);
        exporter->step(0)->add( "H1 error estimator P1" , H1estimatorP1);
        exporter->step(0)->add( "new hsize" , h_new);

        exporter->save();
        Log() << "exportResults done\n";
    }
    /** \endcode */

    Log()<< " real L2 error : "<<Y[0]<<"\n";
    Log()<< " estimated L2 error "<<Y[2]<<"\n";
    Log()<< " real H1 error : "<<Y[1]<<"\n";
    Log()<< " estimated H1 error "<<Y[3]<<"\n";

} // ResidualEstimator::run




template<int Dim, int Order>
typename ResidualEstimator<Dim,Order>::mesh_ptrtype
ResidualEstimator<Dim,Order>::adapt_mesh( p1_element_type& sizefield )
{
#if defined (HAVE_MADLIB_H)
    saveGMSHMesh( _filename="inputmesh.msh",
                  _mesh=sizefield.mesh() );
    std::ostringstream mshstr;
    mshstr << shape << "-" << Dim << ".msh";
    //msh_name = mshstr.str();
    MAd::pGModel model = 0;
    if ( this->vm()["gmshmodel"].template as<bool>() == true )
    {
        GM_create(&model,"theModel");
        if ( this->vm()["gmshgeo"].template as<bool>() == true )
        {
            std::ostringstream geostr;
            geostr << shape << "-" << Dim <<  ".geo";
            GM_readFromGEO( model, geostr.str().c_str() );
        }
        else
        {
            GM_readFromMSH(model, mshstr.str().c_str());
        }
    }
    MAd::pMesh amesh = MAd::M_new(model);
    MAd::M_load(amesh,mshstr.str().c_str());

    MAd::PWLSField * sizeField = new MAd::PWLSField(amesh);
    sizeField->setCurrentSize();
    auto _elit = P1h->mesh()->beginElement();
    auto _elen = P1h->mesh()->endElement();
    for( ; _elit != _elen; ++_elit )
    {
        for( int l = 0; l < _elit->numPoints; ++l )
        {
            int dof = P1h->dof()->localToGlobal( _elit->id(), l, 0 ).get<0>();
            int pid = _elit->point( l ).id()+1;
            sizeField->setSize( pid , sizefield(dof) );
        }
    }

    MAd::MeshAdapter* ma = new MAd::MeshAdapter(amesh,sizeField);

    std::cout << "Statistics before optimization: \n";
    ma->printStatistics(std::cout);
    ma->writePos("meanRatioBefore.pos",MAd::OD_MEANRATIO);

    // Optimize
    // ---------
    std::cout << "Optimizing mesh...\n\n";
    ma->run();

    // Outputs final mesh
    // -------------------
    std::cout << "Statistics after optimization: \n";
    ma->printStatistics(std::cout);
    ma->writePos("meanRatioAfter.pos",MAd::OD_MEANRATIO);

    MAd::M_writeMsh (amesh, "result.msh", 2, NULL);

    return loadGMSHMesh( _mesh=new mesh_type,
                         _filename="result.msh",
                         _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
    //# endmarker7 #
    /** \endcode */
#endif // HAVE_MADLIB_H
}// ResidualEstimator::adapt_mesh



