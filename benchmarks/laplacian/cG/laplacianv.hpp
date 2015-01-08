/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


  This file is part of the Feel library

  Author(s):
  Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Perrimond Benoit <Benoit.Perrimond@bvra.e.ujf-grenoble.fr>
  Vincent Chabannes <vincent.chabannes@gmail.com>

  Date: 2008-02-07

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file laplacianv.cpp
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
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>


#include <feel/feelvf/vf.hpp>
#include <fstream>
#include <sstream>

#include <feel/feelcore/applicationxml.hpp>
#include <feel/feelcore/xmlparser.hpp>

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description laplacianvoptions( "Laplacian Vectorial options" );
    laplacianvoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ( "nu", po::value<double>()->default_value( 1 ), "coef diffusion" )
    ( "beta", po::value<double>()->default_value( 1 ), "coef reaction " )

    ( "weak", "fully weak formulation (including Dirichlet conditions)" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return laplacianvoptions.add( Feel::feel_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "laplacianv" ,
                     "laplacianv" ,
                     "0.2",
                     "nD(n=1,2,3) Laplacian Vectorial on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Benoit Perrimond", "developer", "Benoit.Perrimond@bvra.e.ujf-grenoble.fr", "" );
    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@gmail.com", "" );
    return about;

}


/**
 * LaplacianV Solver using discontinous approximation spaces
 *
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 */
template<int Dim, int Order, int RDim>
class LaplacianV
    :
public ApplicationXML
{
    typedef ApplicationXML super;
public:

    // -- TYPEDEFS --
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Simplex<Dim, 1,RDim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Vectorial> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    LaplacianV( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        backend( backend_type::build( this->vm() ) ),
        meshSize( doption("hsize") ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
    {

        Parameter h;

        if ( Dim == 1 )
            h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.01:0.09:0.2" );

        else if ( Dim == 2 )
            h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.04:0.09:0.2" );

        else if ( Order < 5 )
            h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.06:0.09:0.2" );

        else
            h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.1:0.15:0.2" );

        this->
        addParameter( Parameter( _name="dim",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Dim  ).c_str() ) )
        .addParameter( Parameter( _name="order",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Order  ).c_str() ) )
        .addParameter( Parameter( _name="beta",_type=CONT_ATTR,_latex="\\beta",_values="0.01:1:10" ) )
        .addParameter( Parameter( _name="nu",_type=CONT_ATTR,_latex="\\nu",_values="0.01:1:10" ) )
        .addParameter( h );

        std::vector<Parameter> depend;
        std::vector<std::string> funcs;
        depend.push_back( h );
        std::ostringstream oss;
        oss << "h**" << boost::lexical_cast<std::string>( Order + 1  ) ;
        funcs.push_back( oss.str() );
        oss.str( "" );
        std::vector<std::string> funcs2;
        oss << "h**" << boost::lexical_cast<std::string>( Order ) ;
        funcs2.push_back( oss.str() );

        this->
        addOutput( Output( _name="norm_L2",_latex="\\left\\| . \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs ) )
        .addOutput( Output( _name="norm_H1",_latex="\\left\\| . \\right\\|_{H^1}",_dependencies=depend,_funcs=funcs2 ) );

    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize );

    /**
     * run the convergence test
     */
    void run();

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

    export_ptrtype exporter;

}; // LaplacianV

template<int Dim, int Order, int RDim> const uint16_type LaplacianV<Dim,Order, RDim>::imOrder;

template<int Dim, int Order, int RDim>
typename LaplacianV<Dim,Order, RDim>::mesh_ptrtype
LaplacianV<Dim,Order, RDim>::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );

    GmshHypercubeDomain td( entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,entity_type::is_hypercube );
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( -1, 1 ) );

    if ( Dim>=2 )  td.setY( std::make_pair( -1, 1 ) );

    if ( Dim==3 )  td.setZ( std::make_pair( -1, 1 ) );

    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    return mesh;
} // LaplacianV::createMesh


template<int Dim, int Order, int RDim>
void
LaplacianV<Dim, Order, RDim>::run()
{
    this->addParameterValue( Dim )
    .addParameterValue( Order )
    .addParameterValue( doption("beta") )
    .addParameterValue( doption("nu") )
    .addParameterValue( doption("hsize") );

    if ( this->preProcessing() == RUN_EXIT ) return;

    //    int maxIter = 10.0/meshSize;
    using namespace Feel::vf;

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createMesh( meshSize );

    /*
     * The function space and some associate elements are then defined
     */
    space_ptrtype Xh = space_type::New( mesh );
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";
    value_type nu = doption("nu");
    value_type beta = doption("beta");


    value_type pi = M_PI;
    AUTO( g, sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() )*one() );
    AUTO( f, ( pi*pi*Dim+1 )*g );
    AUTO( zf, 0*one() );

    boost::timer t1;

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

    // Construction of the right hand side

    vector_ptrtype F( backend->newVector( Xh ) );

    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements( mesh ),
                   trans( f )*id( v ) );
    double M_gammabc = 20;

    if ( this->vm().count( "weak" ) )
    {
        form1( Xh, F ) +=
            integrate( markedfaces( mesh,tag1 ),
                       trans( zf )*( -nu*grad( v )*N()+M_gammabc*id( v )/hFace() ) );
        form1( Xh, F ) +=
            integrate( markedfaces( mesh,tag2 ),
                       trans( zf )*( -nu*grad( v )*N()+M_gammabc*id( v )/hFace() ) );

    }

    F->close();
    LOG(INFO) << "F assembled in " << t1.elapsed() << "s\n";
    t1.restart();

    //Construction of the left hand side

    sparse_matrix_ptrtype D( backend->newMatrix( Xh, Xh ) );

    form2( Xh, Xh, D, _init=true );
    LOG(INFO) << "D initialized in " << t1.elapsed() << "s\n";
    t1.restart();

    form2( Xh, Xh, D ) +=
        integrate( elements( mesh ),
                   nu*( trace( gradt( u )*trans( grad( v ) ) ) )
                   + beta*( trans( idt( u ) )*id( v ) ) );
    LOG(INFO) << "D stiffness+mass assembled in " << t1.elapsed() << "s\n";
    t1.restart();

    if ( this->vm().count( "weak" ) )
    {
        form2( Xh, Xh, D ) += integrate( markedfaces( mesh,tag1 ),
                                         ( - nu*trans( id( v ) )*( gradt( u )*N() )
                                           - nu*trans( idt( u ) )*( grad( v )*N() )
                                           + M_gammabc*trans( idt( u ) )*id( v )/hFace() ) );
        form2( Xh, Xh, D ) += integrate( markedfaces( mesh,tag2 ),
                                         ( - nu*trans( id( v ) )*( gradt( u )*N() )
                                           - nu*trans( idt( u ) )*( grad( v )*N() )
                                           + M_gammabc*trans( idt( u ) )*id( v )/hFace() ) );
        LOG(INFO) << "D weak bc assembled in " << t1.elapsed() << "s\n";
        t1.restart();
    }

    D->close();



    if ( !this->vm().count( "weak" ) )
    {
        form2( Xh, Xh, D ) +=
            on( markedfaces( mesh, tag1 ), u, F, zf )+
            on( markedfaces( mesh, tag2 ), u, F, zf );
        LOG(INFO) << "D strong bc assembled in " << t1.elapsed() << "s\n";
        t1.restart();
    }


    if ( this->vm().count( "export-matlab" ) )
    {
        D->printMatlab( "D.m" );
        F->printMatlab( "F.m" );
    }

    this->solve( D, u, F );

    LOG(INFO) << "solve in " << t1.elapsed() << "s\n";
    t1.restart();

    double L2error2 =integrate( elements( mesh ),
                                trans( idv( u )-g )*( idv( u )-g ) ).evaluate()( 0, 0 );
    double L2error =   math::sqrt( L2error2 );

    LOG(INFO) << "||error||_L2=" << L2error << "\n";
    LOG(INFO) << "L2 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();


    v = vf::project( Xh, elements( mesh ), g );
    double semiH1error2 =integrate( elements( mesh ),
                                    trace( ( gradv( u )-gradv( v ) )*trans( gradv( u )-gradv( v ) ) ) ).evaluate()( 0, 0 ) ;

    LOG(INFO) << "semi H1 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();

    double H1error =   math::sqrt( semiH1error2+L2error2 );
    LOG(INFO) << "||error||_H1=" << H1error << "\n";


    LOG(INFO) << "H1 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();

    this->exportResults( u, v );


    this->addOutputValue( L2error ).addOutputValue( H1error );
    this->postProcessing();

} // LaplacianV::run

template<int Dim, int Order, int RDim>
void
LaplacianV<Dim, Order, RDim>::solve( sparse_matrix_ptrtype& D,
                                     element_type& u,
                                     vector_ptrtype& F )
{


    vector_ptrtype U( backend->newVector( u.functionSpace() ) );
    backend->solve( D, D, U, F );
    u = *U;
} // LaplacianV::solve


template<int Dim, int Order, int RDim>
void
LaplacianV<Dim, Order, RDim>::exportResults( element_type& U, element_type& E )
{
    LOG(INFO) << "exportResults starts\n";

    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );

        exporter->step( 0 )->add( "pid",
                                  regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );
        exporter->step( 0 )->add( "u", U );
        exporter->step( 0 )->add( "exact", E );

        exporter->save();
    }
} // LaplacianV::export

