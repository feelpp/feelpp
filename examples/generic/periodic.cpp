/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-10-03

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
   \file periodic.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-10-03
 */
#include <life/options.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/region.hpp>
#include <life/lifediscr/operatorlinear.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>


#include <life/lifevf/vf.hpp>




inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description periodicoptions("Periodic Laplacian options");
    periodicoptions.add_options()
        ("hsize", Life::po::value<double>()->default_value( 0.1 ), "mesh size in domain")

        ("penalbc", Life::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions")

        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return periodicoptions.add( Life::life_options() );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "periodic" ,
                           "periodic" ,
                           "0.1",
                           "nD(n=1,2,3) Periodic Laplacian on simplices or simplex products",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2008 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}

std::pair<std::string,std::string> createRing( int Dim, double h, double rmin, double rmax );

namespace Life
{
using namespace vf;
/**
 * Fat boundary method for the laplacian
 *
 */
template<int Dim, int Order>
class PeriodicLaplacian
    :
    public Application
{
    typedef Application super;
public:

#define Entity Simplex

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
    typedef Entity<Dim,1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type, Periodic<2,4,value_type> > functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;

    /*quadrature*/
    template<int IMORDER> struct MyIM : public IM<Dim, IMORDER, value_type, Entity> {};

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /** constructor */
    PeriodicLaplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od );

    /** mesh generation */
    mesh_ptrtype createMesh();

    /**
     * run the convergence test
     */
    void run();

private:



    /**
     * solve the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type&, element_type& );

private:

    backend_ptrtype M_backend;

    double h;
    double penalisation_bc;

    mesh_ptrtype mesh;

    functionspace_ptrtype Xh;

    export_ptrtype exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;

}; // Periodic

template<int Dim, int Order>
PeriodicLaplacian<Dim,Order>::PeriodicLaplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),

    // Data
    h( this->vm()["hsize"].template as<double>() ),
    penalisation_bc( this->vm()["penalbc"].template as<value_type>() ),

    // spaces
    Xh(),

    // export
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),

    //
    timers()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }



    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % h
                            );

    Log() << "create mesh\n";
    mesh = createMesh();

    Log() << "create space\n";
    node_type trans(2);
    trans[0]=0;
    trans[1]=2;
    Xh = functionspace_type::New( mesh, MESH_COMPONENTS_DEFAULTS, Periodic<2,4,value_type>( trans ) );

    Log() << "print space info\n";
    Xh->printInfo();

    Log() << "Constructor done\n";
}
template<int Dim, int Order>
typename PeriodicLaplacian<Dim,Order>::mesh_ptrtype
PeriodicLaplacian<Dim,Order>::createMesh()
{
    timers["mesh"].first.restart();
    mesh_ptrtype mesh( new mesh_type );
    //mesh->setRenumber( false );

    GmshTensorizedDomain<entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,Entity> td;
    td.setVersion( 2 );
    td.setCharacteristicLength( h );
    td.setX( std::make_pair( -1., 1. ) );
    td.setY( std::make_pair( -1., 1. ) );
    std::string fname = td.generate( "square" );

    ImporterGmsh<mesh_type> import( fname );
    import.setVersion( "2.0" );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    Log() << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;
} // PeriodicLaplacian::createMesh

template<int Dim, int Order>
void
PeriodicLaplacian<Dim, Order>::run()
{
    //    int maxIter = 10.0/meshSize;
    using namespace Life::vf;

    timers["init"].first.restart();

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

#if 0
    AUTO( g, Px()*Py()+2*Px()+1 );
    AUTO( grad_g, vec(
                      Py()+2,
                      Px()
                      )
          );
    AUTO( f, 0 );
#else
    AUTO( g, sin(M_PI*Px())*cos(M_PI*Py()));
    AUTO( grad_g, vec(
                      +M_PI*cos(M_PI*Px())*cos(M_PI*Py()),
                      -M_PI*sin(M_PI*Px())*sin(M_PI*Py())
                      )
          );
    AUTO( f, 2*M_PI*M_PI*g );
#endif

    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    form2( Xh, Xh, M, _init=true ) = integrate( elements( mesh ), MyIM<2*(Order-1)>(), gradt(u)*trans(grad(v)) );
    form2( Xh, Xh, M ) += integrate( markedfaces( mesh, 1 ), MyIM<2*(Order-1)>(),
                                     -gradt(u)*N()*id(v)
                                     -grad(v)*N()*idt(u)
                                     + penalisation_bc*id(u)*idt(v)/hFace() );
    form2( Xh, Xh, M ) += integrate( markedfaces( mesh, 3 ), MyIM<2*(Order-1)>(),
                                     -gradt(u)*N()*id(v)
                                     -grad(v)*N()*idt(u)
                                     + penalisation_bc*id(u)*idt(v)/hFace() );

    M->close();

    double area = integrate( elements(mesh), MyIM<0>(), constant(1.0) ).evaluate()( 0, 0);
    double mean = integrate( elements(mesh), MyIM<5>(), g ).evaluate()( 0, 0)/area;
    Log() << "int g  = " << mean << "\n";
    vector_ptrtype F( M_backend->newVector( Xh ) );
    form1( Xh, F, _init=true ) = ( integrate( elements( mesh ), MyIM<Order+5>(), f*id(v) )
                                   //+integrate( boundaryfaces( mesh ), MyIM<Order+5>(), (trans(grad_g)*N())*id(v) )

                                   );
    F->close();

    if ( this->vm().count( "export-matlab" ) )
        {
            M->printMatlab( "M.m" );
            F->printMatlab( "F.m" );
        }

    this->solve( M, u, F );

    Log() << "area   = " << area << "\n";
    Log() << "int g  = " << integrate( elements(mesh), MyIM<5>(), g ).evaluate()( 0, 0)/area << "\n";
    Log() << "int u  = " << integrate( elements(mesh), MyIM<Order>(), idv(u) ).evaluate()( 0, 0)/area << "\n";
    Log() << "error  = " << math::sqrt( integrate( elements(mesh), MyIM<10>(), (idv(u)-g)*(idv(u)-g) ).evaluate()( 0, 0) ) << "\n";

    v = vf::project( Xh, elements(mesh), g );

    element_type e( Xh, "e" );
    e = vf::project( Xh, elements(mesh), idv(u)-g );


    exportResults( u, v, e );


} // PeriodicLaplacian::run

template<int Dim, int Order>
void
PeriodicLaplacian<Dim, Order>::solve( sparse_matrix_ptrtype& D,
                                element_type& u,
                                vector_ptrtype& F )
{
    timers["solver"].first.restart();


    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    M_backend->solve( D, D, U, F );
    u = *U;

    timers["solver"].second = timers["solver"].first.elapsed();
    Log() << "[timer] solve: " << timers["solver"].second << "\n";
} // PeriodicLaplacian::solve


template<int Dim, int Order>
void
PeriodicLaplacian<Dim, Order>::exportResults( element_type& U, element_type& V, element_type& E )

{
    timers["export"].first.restart();

    Log() << "exportResults starts\n";

    exporter->step(1.)->setMesh( U.functionSpace()->mesh() );
    exporter->step(1.)->add( "u", U );
    exporter->step(1.)->add( "exact", V );
    exporter->step(1.)->add( "error", E );

    exporter->save();
    timers["export"].second = timers["export"].first.elapsed();
    Log() << "[timer] exportResults(): " << timers["export"].second << "\n";
} // PeriodicLaplacian::export
} // Life




int
main( int argc, char** argv )
{
    using namespace Life;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 4;

    typedef Life::PeriodicLaplacian<nDim, nOrder> laplacian_periodic_type;

    /* define and run application */
    laplacian_periodic_type laplacian_periodic( argc, argv, makeAbout(), makeOptions() );

    laplacian_periodic.run();
}









