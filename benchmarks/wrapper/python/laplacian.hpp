/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


  This file is part of the Feel library

  Author(s):
  Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Perrimond Benoit <Benoit.Perrimond@bvra.e.ujf-grenoble.fr>
  Vincent Chabannes <vincent.chabannes@gmail.com>

  Date: 2008-02-07

  Copyright (C) 2008-2010 Université Joseph Fourier (Grenoble I)

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
   \file laplacian.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Perrimond Benoit <Benoit.Perrimond@bvra.e.ujf-grenoble.fr>
   \author Vincent Chabannes <vincent.chabannes@gmail.com>
   \date 2008-02-07
 */
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>


#include <feel/feelvf/vf.hpp>
#include <fstream>
#include <sstream>

namespace Feel
{

inline
po::options_description
makeOptions()
{
    po::options_description laplacianoptions( "Laplacian options" );
    laplacianoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ( "nu", po::value<double>()->default_value( 1 ), "coef diffusion" )
    ( "beta", po::value<double>()->default_value( 1 ), "coef reaction " )
    ( "gammabc", po::value<double>()->default_value( 20 ), "weak Dirichlet penalisation parameter " )

    ( "weak", "use weak dirichlet conditions" )
    ;
    return laplacianoptions.add( Feel::feel_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "laplacian2" ,
                     "laplacian2" ,
                     "0.2",
                     "nD(n=1,2,3) Laplacian on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008 Universit� Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Benoit Perrimond", "developer", "Benoit.Perrimond@bvra.e.ujf-grenoble.fr", "" );
    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@gmail.com", "" );
    return about;

}


/**
 * Laplacian Solver using discontinous approximation spaces
 *
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 */
template<int Dim, int Order, int RDim = Dim, template<uint16_type,uint16_type,uint16_type> class Entity=Simplex>
class Laplacian
    :
public Application
{
    typedef Application super;
public:


    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim, 1,RDim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous> p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;
    typedef fusion::vector<Lagrange<Order+1, Scalar> > exact_basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef FunctionSpace<mesh_type, exact_basis_type, value_type> exact_space_type;
    typedef boost::shared_ptr<exact_space_type> exact_space_ptrtype;
    typedef typename exact_space_type::element_type exact_element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    Laplacian( AboutData const& ad )
        :
        super( ad ),
        backend(),
        meshSize( 0.1 ),

        M_use_weak_dirichlet( true ),
        M_gammabc( 20 )

        //exporter( Exporter<mesh_type>::New( this->about().appName() ) )
    {
        this->setLogs();

        if ( M_use_weak_dirichlet )
            LOG(INFO)  << "use weak Dirichlet BC\n";

        if ( exporter && exporter->doExport() )
            LOG(INFO)  << "export results to ensight format\n";
    }

    Laplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        backend( backend_type::build( this->vm() ) ),
        meshSize( doption("hsize") ),

        M_use_weak_dirichlet( this->vm().count( "weak" ) ),
        M_gammabc( doption("gammabc") ),

        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
    {
        if ( M_use_weak_dirichlet )
            LOG(INFO)  << "use weak Dirichlet BC\n";

        if ( exporter->doExport() )
            LOG(INFO)  << "export results to ensight format\n";
    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize );

    /**
     * run the convergence test
     */
    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );
private:

    /**
     * solve the system D u = F
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type& v );

private:

    backend_ptrtype backend;

    double meshSize;

    bool M_use_weak_dirichlet;
    double M_gammabc;

    export_ptrtype exporter;

}; // Laplacian

template<int Dim, int Order, int RDim, template<uint16_type,uint16_type,uint16_type> class Entity>
typename Laplacian<Dim,Order,RDim,Entity>::mesh_ptrtype
Laplacian<Dim,Order,RDim,Entity>::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );

    GmshHypercubeDomain td( entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,entity_type::is_hypercube );
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( -1, 1 ) );

    if ( ( Dim==1 ) && ( RDim==2 ) )  td.setY( std::make_pair( 1, 1 ) );

    if ( Dim>=2 )  td.setY( std::make_pair( -1, 1 ) );

    if ( ( Dim==2 ) && ( RDim==3 ) )  td.setZ( std::make_pair( 1, 1 ) );

    if ( Dim==3 )  td.setZ( std::make_pair( -1, 1 ) );

    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    return mesh;
} // Laplacian::createMesh

template<int Dim, int Order, int RDim, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, Order, RDim, Entity>::run()
{
    std::vector<double> X( 1 );
    std::vector<double> Y( 2 );
    X[0] = doption("hsize");
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template<int Dim, int Order, int RDim, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, Order, RDim, Entity>::run( const double* X, unsigned long P,
        double* Y, unsigned long YN )
{
    boost::timer t1;

    value_type nu = 1;
    value_type beta = 1;

    meshSize = X[0];

    using namespace Feel::vf;

    /*
     * allocate the backend: this must be placed here as the matrices, vectors
     * and preconditioners allocated witha previous backend may be invalidated
     * by a change of h for example in the case of pythin scripting
     */
    //backend = backend_ptrtype( backend_type::build( BACKEND_PETSC ) );
    backend = backend_ptrtype( backend_type::build( soption("backend") ) );

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createMesh( meshSize );
    LOG(INFO)  << "mesh created in " << t1.elapsed() << "s\n";
    t1.restart();

    /*
     * The function space and some associate elements are then defined
     */
    space_ptrtype Xh = space_type::New( mesh );
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    LOG(INFO)  << "[functionspace] Number of dof " << Xh->nLocalDof() << "\n";
    LOG(INFO)  << "function space and elements created in " << t1.elapsed() << "s\n";
    t1.restart();

    exact_space_ptrtype Eh = exact_space_type::New( mesh );
    exact_element_type fproj( Eh, "f" );
    exact_element_type gproj( Eh, "g" );
    LOG(INFO)  << "[functionspace] Number of dof " << Eh->nLocalDof() << "\n";
    LOG(INFO)  << "function space and elements created in " << t1.elapsed() << "s\n";
    t1.restart();




    value_type pi = M_PI;
    AUTO( g, sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() ) );
    AUTO( f, ( pi*pi*Dim*nu+beta )*g );
    AUTO( zf, 0*Px()+0*Py() );

    int tag1,tag2;

    if ( ( Dim == 1 ) || ( Dim == 2 ) )
    {
        tag1 = 1;
        tag2 = 3;
    }

    else if ( Dim == 3 )
    {
        tag1 = 15;
        tag2 = 23;
    }

    fproj = vf::project( Eh, elements( mesh ), f );
    gproj = vf::project( Eh, elements( mesh ), g );
    v = vf::project( Xh, elements( mesh ), g );

    // Construction of the right hand side

    vector_ptrtype F( backend->newVector( Xh ) );


    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements( mesh ), idv( fproj )*id( v ) );

    if ( M_use_weak_dirichlet )
    {
        form1( Xh, _vector=F ) +=
            integrate( markedfaces( mesh,tag1 ),
                       zf*( -nu*grad( v )*N()+M_gammabc*id( v )/hFace() ) );
        form1( Xh, _vector=F ) +=
            integrate( markedfaces( mesh,tag2 ),
                       zf*( -nu*grad( v )*N()+M_gammabc*id( v )/hFace() ) );

    }

    F->close();
    LOG(INFO)  << "F assembled in " << t1.elapsed() << "s\n";
    t1.restart();

    //Construction of the left hand side

    sparse_matrix_ptrtype D( backend->newMatrix( Xh, Xh ) );


    form2( Xh, Xh, _matrix=D, _init=true );
    LOG(INFO)  << "D initialized in " << t1.elapsed() << "s\n";
    t1.restart();

    form2( Xh, Xh, _matrix=D ) +=
        integrate( elements( mesh ),
                   nu*( gradt( u )*trans( grad( v ) ) )
                   + beta*( idt( u )*id( v ) ) );
    LOG(INFO)  << "D stiffness+mass assembled in " << t1.elapsed() << "s\n";
    t1.restart();

    if ( M_use_weak_dirichlet )
    {

        form2( Xh, Xh, _matrix=D ) += integrate( markedfaces( mesh,tag1 ),
                                         ( - nu*trans( id( v ) )*( gradt( u )*N() )
                                           - nu*trans( idt( u ) )*( grad( v )*N() )
                                           + M_gammabc*trans( idt( u ) )*id( v )/hFace() ) );
        form2( Xh, Xh, _matrix=D ) += integrate( markedfaces( mesh,tag2 ),
                                         ( - nu*trans( id( v ) )*( gradt( u )*N() )
                                           - nu*trans( idt( u ) )*( grad( v )*N() )
                                           + M_gammabc*trans( idt( u ) )*id( v )/hFace() ) );
        LOG(INFO)  << "D weak bc assembled in " << t1.elapsed() << "s\n";
        t1.restart();

    }

    D->close();

    LOG(INFO)  << "D assembled in " << t1.elapsed() << "s\n";


    if ( ! M_use_weak_dirichlet )
    {
        t1.restart();
        form2( Xh, Xh, _matrix=D ) +=
            on( markedfaces( mesh, tag1 ), u, F, g )+
            on( markedfaces( mesh, tag2 ), u, F, g );
        LOG(INFO)  << "Strong Dirichlet assembled in " << t1.elapsed() << "s on faces " << tag1 << " and " << tag2 << " \n";
    }

    t1.restart();

    this->solve( D, u, F );

    LOG(INFO)  << "solve in " << t1.elapsed() << "s\n";
    t1.restart();

    double L2error2 =integrate( elements( mesh ),
                                ( idv( u )-idv( gproj ) )*( idv( u )-idv( gproj ) ) ).evaluate()( 0, 0 );
    double L2error =   math::sqrt( L2error2 );

    LOG(INFO)  << "||error||_L2=" << L2error << "\n";
    LOG(INFO)  << "L2 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();


    double semiH1error2 =integrate( elements( mesh ),
                                    ( gradv( u )-gradv( gproj ) )*trans( gradv( u )-gradv( gproj ) ) ).evaluate()( 0, 0 ) ;

    LOG(INFO)  << "semi H1 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();

    double H1error =   math::sqrt( semiH1error2+L2error2 );


    LOG(INFO)  << "||error||_H1=" << H1error << "\n";
    LOG(INFO)  << "H1 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();

    this->exportResults( u, v );
    Y[0] = L2error;
    Y[1] = H1error;

} // Laplacian::run

template<int Dim, int Order, int RDim, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, Order, RDim, Entity>::solve( sparse_matrix_ptrtype& D,
        element_type& u,
        vector_ptrtype& F )
{
    vector_ptrtype U( backend->newVector( u.functionSpace() ) );
    backend->solve( D, D, U, F );
    u = *U;
} // Laplacian::solve


template<int Dim, int Order, int RDim, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, Order, RDim, Entity>::exportResults( element_type& U, element_type& v )
{
    if ( exporter && exporter->doExport() )
    {
        LOG(INFO)  << "exportResults starts\n";

        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->add( "u", U );
        exporter->step( 0 )->add( "exact", v );

        exporter->save();
    }
} // Laplacian::export

}
