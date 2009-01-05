/* -*- mode: c++ coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-07

  Copyright (C) 2008 Universite Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file laplacian.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-07
 */
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


using namespace Life;
using namespace Life::vf;

inline
po::options_description
makeOptions()
{
    po::options_description laplacianoptions("Laplacian options");
    laplacianoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.5 ), "mesh size")
        ("nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient")
        ("weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ("penaldir", Life::po::value<double>()->default_value( 10 ),
         "penalisation parameter for the weak boundary Dirichlet formulation")
        ;
    return laplacianoptions.add( Life::life_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "laplacian" ,
                     "laplacian" ,
                     "0.2",
                     "nD(n=1,2,3) Laplacian on simplices or simplex products",
                     Life::AboutData::License_GPL,
                     "Copyright (c) 2008 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


/**
 * Laplacian Solver using discontinous approximation spaces
 *
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 */
template<int Dim>
class Laplacian
    :
    public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type Order = 2;
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
    typedef Simplex<Dim, 1,Dim> entity_type;
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /* store partitioning */
    typedef FunctionSpace<mesh_type, fusion::vector<fem::Lagrange<Dim, 0, Scalar, Discontinuous> > > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<fem::Lagrange<Dim, Order, Scalar, Continuous, double, Simplex> > basis_type;

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

    Laplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm()["exporter"].template as<std::string>() )->setOptions( this->vm() ) ),
        timeSet( new timeset_type( "laplacian" ) )
    {
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "laplacian" );
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
    void exportResults( element_type& u );

private:

    backend_ptrtype M_backend;

    double meshSize;

    export_ptrtype exporter;
    typename export_type::timeset_ptrtype timeSet;

}; // Laplacian

template<int Dim> const uint16_type Laplacian<Dim>::Order;
template<int Dim> const uint16_type Laplacian<Dim>::imOrder;

template<int Dim>
typename Laplacian<Dim>::mesh_ptrtype
Laplacian<Dim>::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );

    GmshTensorizedDomain<entity_type::nDim,
                         entity_type::nOrder,
                         entity_type::nDim,
                         Simplex> td;
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( -1, 1 ) );
    td.setY( std::make_pair( -1, 1 ) );
    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    return mesh;
} // Laplacian::createMesh


template<int Dim>
void
Laplacian<Dim>::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    /*
     * we change to the directory where the results and logs will be
     * stored
     */
    this->changeRepository( boost::format( "doc/tutorial/%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
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


    value_type pi = M_PI;
    AUTO( g, sin(pi*Px())*cos(pi*Py())*cos(pi*Pz()) );
    AUTO( f, pi*pi*Dim*g );

    bool weakdir = this->vm()["weakdir"].template as<int>();
    value_type penaldir = this->vm()["penaldir"].template as<double>();


    /*
     * Construction of the right hand side
     */
    vector_ptrtype F( M_backend->newVector( Xh ) );

    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements(mesh), im_type(),
                   f*id(v) );
    if ( Application::nProcess () != 1 || weakdir )
        {
            form1( _test=Xh, _vector=F ) +=
                integrate( markedfaces(mesh,1), im_type(),
                           g*(-grad(v)*N()+penaldir*id(v)/hFace()) ) +
                integrate( markedfaces(mesh,3), im_type(),
                           g*(-grad(v)*N()+penaldir*id(v)/hFace()) );
        }
    F->close();

    /*
     * Construction of the left hand side
     */
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );

    value_type nu = this->vm()["nu"].template as<double>();


    form2( Xh, Xh, D, _init=true ) =
        integrate( elements(mesh), im_type(),
                   nu*gradt(u)*trans(grad(v)) );
    if ( Application::nProcess () != 1 || weakdir )
        {
            form2( Xh, Xh, D ) +=
                integrate( markedfaces(mesh,1), im_type(),
                           -(gradt(u)*N())*id(v)
                           -(grad(v)*N())*idt(u)
                           +penaldir*id(v)*idt(u)/hFace()) +
                integrate( markedfaces(mesh,3), im_type(),
                           -(gradt(u)*N())*id(v)
                           -(grad(v)*N())*idt(u)
                           +penaldir*id(v)*idt(u)/hFace());
            D->close();
        }
    else
        {
            D->close();
            form2( Xh, Xh, D ) +=
                on( markedfaces(mesh, 1), u, F, g )+
                on( markedfaces(mesh, 3), u, F, g );

        }

    this->solve( D, u, F );

    double L2error2 =integrate(elements(mesh), im_type(),
                               (idv(u)-g)*(idv(u)-g) ).evaluate()(0,0);
    double L2error =   math::sqrt( L2error2 );
    Log() << "||error||_L2=" << L2error << "\n";

    this->exportResults( u );
} // Laplacian::run

template<int Dim>
void
Laplacian<Dim>::solve( sparse_matrix_ptrtype& D,
                       element_type& u,
                       vector_ptrtype& F )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    M_backend->solve( D, D, U, F );
    u = *U;
} // Laplacian::solve


template<int Dim>
void
Laplacian<Dim>::exportResults( element_type& U )
{
    Log() << "exportResults starts\n";

    typename timeset_type::step_ptrtype timeStep = timeSet->step( 0 );
    timeStep->setMesh( U.functionSpace()->mesh() );

    timeStep->add( "pid",
                   regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );
    timeStep->add( "u", U );

    exporter->save();
} // Laplacian::export




int
main( int argc, char** argv )
{
    Laplacian<2> laplacian( argc, argv, makeAbout(), makeOptions() );

    laplacian.run();
}





