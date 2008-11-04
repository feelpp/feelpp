/* -*- mode: c++ -*-


  This file is part of the Life library

  Author(s):
  Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
  Perrimond Benoit <Benoit.Perrimond@bvra.e.ujf-grenoble.fr>
  Vincent Chabannes <vincent.chabannes@gmail.com>

  Date: 2008-02-07

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file laplacianv.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \author Perrimond Benoit <Benoit.Perrimond@bvra.e.ujf-grenoble.fr>
   \author Vincent Chabannes <vincent.chabannes@gmail.com>
   \date 2008-02-07
 */
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include <life/options.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/region.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>


#include <life/lifevf/vf.hpp>
#include <fstream>

using namespace Life;

inline
po::options_description
makeOptions()
{
    po::options_description laplacianvoptions("Laplacian Vectorial options");
    laplacianvoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.5 ), "mesh size")
        ("nu", po::value<double>()->default_value( 1 ), "coef diffusion")
        ("beta", po::value<double>()->default_value( 1 ), "coef reaction " )
        ;
    return laplacianvoptions.add( Life::life_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "laplacianv" ,
                     "laplacianv" ,
                     "0.2",
                     "nD(n=1,2,3) Laplacian Vectorial on simplices or simplex products",
                     Life::AboutData::License_GPL,
                     "Copyright (c) 2008 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    about.addAuthor("Benoit Perrimond", "developer", "Benoit.Perrimond@bvra.e.ujf-grenoble.fr", "");
    about.addAuthor("Vincent Chabannes", "developer", "vincent.chabannes@gmail.com", "");
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
    public Application
{
    typedef Application super;
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
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<fem::Lagrange<Dim, 0, Scalar, Discontinuous> > > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<fem::Lagrange<Dim, Order, Vectorial, Continuous, double, Simplex> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /*quadrature*/
    typedef IM<Dim, imOrder, value_type, Simplex> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef typename export_type::timeset_type timeset_type;

    LaplacianV( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm()["exporter"].template as<std::string>() )->setOptions( this->vm() ) ),
        timeSet( new timeset_type( "laplacianv" ) )
    {
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "laplacianv" );
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
    typename export_type::timeset_ptrtype timeSet;

}; // LaplacianV

template<int Dim, int Order, int RDim> const uint16_type LaplacianV<Dim,Order, RDim>::imOrder;

template<int Dim, int Order, int RDim>
typename LaplacianV<Dim,Order, RDim>::mesh_ptrtype
LaplacianV<Dim,Order, RDim>::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );

    GmshTensorizedDomain<entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,Simplex> td;
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( -1, 1 ) );
    if(Dim>=2)  td.setY( std::make_pair( -1, 1 ) );
    if(Dim==3)  td.setZ( std::make_pair( -1, 1 ) );
    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    return mesh;
} // LaplacianV::createMesh


template<int Dim, int Order, int RDim>
void
LaplacianV<Dim, Order, RDim>::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    //    int maxIter = 10.0/meshSize;
    using namespace Life::vf;

    this->changeRepository( boost::format( "%1%/nu_%2%/beta_%3%/%4%/P%5%/h_%6%/" )
                            % this->about().appName()
                            % this->vm()["nu"].template as<double>()
                            % this->vm()["beta"].template as<double>()
                            % entity_type::name()
                            % Order
                            % this->vm()["hsize"].template as<double>()
                            );



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

    Log() << "Number of dof " << Xh->nLocalDof() << "\n";
    value_type nu = this->vm()["nu"].template as<double>();
    value_type beta = this->vm()["beta"].template as<double>();


    value_type pi = M_PI;
    AUTO( g, sin(pi*Px())*cos(pi*Py())*cos(pi*Pz())*one() );
    AUTO( f, (pi*pi*Dim+1)*g );
    AUTO( zf, 0*one() );

    boost::timer t1;

    // Construction of the right hand side

    vector_ptrtype F( backend->newVector( Xh ) );

    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements(mesh), im_type(),
                   trans(f)*id(v) );

    F->close();
    Log() << "F assembled in " << t1.elapsed() << "s\n";
    t1.restart();

    //Construction of the left hand side

    sparse_matrix_ptrtype D( backend->newMatrix( Xh, Xh ) );

    form2( Xh, Xh, D, _init=true ) =
        integrate( elements(mesh), im_type(),
                   nu*(trace(gradt(u)*trans(grad(v))))
                   + beta*(trans(idt(u))*id(v)) );

    D->close();

    Log() << "D assembled in " << t1.elapsed() << "s\n";
    t1.restart();
    if ( Dim == 1 )
        {
            form2( Xh, Xh, D ) +=
                on( markedfaces(mesh, 1), u, F, zf )+
                on( markedfaces(mesh, 2), u, F, zf );
        }
    else if ( Dim == 2 )
        {
            form2( Xh, Xh, D ) +=
                on( markedfaces(mesh, 1), u, F, zf )+
                on( markedfaces(mesh, 3), u, F, zf );

        }
    else if ( Dim == 3 )
        {
            form2( Xh, Xh, D ) +=
                on( markedfaces(mesh, 15), u, F, zf )+
                on( markedfaces(mesh, 23), u, F, zf );

        }
    Log() << "D+Dirichlet assembled in " << t1.elapsed() << "s\n";
    t1.restart();

    this->solve( D, u, F );

    Log() << "solve in " << t1.elapsed() << "s\n";
    t1.restart();

    double L2error2 =integrate( elements(mesh), im_type(),
                                trans(idv(u)-g)*(idv(u)-g) ).evaluate()( 0, 0 );
    double L2error =   math::sqrt( L2error2 );

    Log() << "||error||_L2=" << L2error << "\n";
    Log() << "L2 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();


    v = project( Xh, elements(mesh), g );
    double semiH1error2 =integrate( elements(mesh), IM<Dim, 2*Order-2, value_type, Simplex>(),
                                    trace((gradv(u)-gradv(v))*trans(gradv(u)-gradv(v))) ).evaluate()( 0, 0 ) ;

    Log() << "semi H1 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();

    double H1error =   math::sqrt( semiH1error2+L2error2 );
    Log() << "||error||_H1=" << H1error << "\n";


    Log() << "H1 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();

    this->exportResults( u, v );


    this->changeRepository( boost::format( "%1%/TestConv/%2%D/P%3%/nu_%4%/beta_%5%/")
                            % this->about().appName()
                            % Dim
                            % Order
                            % this->vm()["nu"].template as<double>()
                            % this->vm()["beta"].template as<double>()
                            );



    std::ofstream fichier("resultTest.txt",std::ios_base::app);
    fichier <<this->vm()["hsize"].template as<double>()
            <<" "<<H1error<<" "<<L2error<<"\n";

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
    Log() << "exportResults starts\n";

    typename timeset_type::step_ptrtype timeStep = timeSet->step( 0 );
    timeStep->setMesh( U.functionSpace()->mesh() );

    timeStep->add( "pid",
                   regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );
    timeStep->add( "u", U );
    timeStep->add( "exact", E );

    exporter->save();
} // LaplacianV::export

