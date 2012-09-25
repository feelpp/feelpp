/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-10-03

  Copyright (C) 2008-2012 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-10-03
 */
#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description periodicoptions( "Periodic Laplacian options" );
    periodicoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "mesh size in domain" )

    ( "penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions" )

    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return periodicoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "periodic" ,
                           "periodic" ,
                           "0.1",
                           "nD(n=1,2,3) Periodic Laplacian on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008-2012 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

std::pair<std::string,std::string> createRing( int Dim, double h, double rmin, double rmax );

namespace Feel
{
using namespace vf;
/**
 * Fat boundary method for the laplacian
 *
 */
template<int Dim, int Order>
class PeriodicLaplacian
    :
public Simget
{
    typedef Simget super;
public:

#define Entity Simplex

    // -- TYPEDEFS --
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef bases<Lagrange<Order, Scalar> > basis_type;
    typedef bases<Lagrange<Order, Vectorial> > basis_vec_type;

    typedef bases<Lagrange<Order, Scalar>, Lagrange<Order, Scalar> > basis_composite_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, Periodic<> > functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef FunctionSpace<mesh_type, basis_vec_type, Periodic<> > functionspace_vec_type;
    typedef boost::shared_ptr<functionspace_vec_type> functionspace_vec_ptrtype;

    //typedef FunctionSpace<mesh_type, basis_composite_type, periodic<Periodic<>, Periodic<> > > functionspace_composite_type;
    typedef FunctionSpace<mesh_type, basis_composite_type, Periodic<> > functionspace_composite_type;
    typedef boost::shared_ptr<functionspace_composite_type> functionspace_composite_ptrtype;


    typedef typename functionspace_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /** constructor */
    PeriodicLaplacian();

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type&, element_type&, element_type& );

private:

    backend_ptrtype M_backend;

    double h;
    double penalisation_bc;

    mesh_ptrtype mesh;

    functionspace_ptrtype Xh;
    functionspace_composite_ptrtype Xhc;
    functionspace_vec_ptrtype Xh_vec;

    export_ptrtype exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;

}; // Periodic

template<int Dim, int Order>
PeriodicLaplacian<Dim,Order>::PeriodicLaplacian()
    :
    super(),
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
    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % h
                          );

    LOG(INFO) << "create mesh\n";
    const std::string shape = "hypercube";
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                         _usenames=false,
                                         _shape=shape,
                                         _dim=Dim,
                                         _h=h,
                                         _xmin=-1, _ymin=-1 ) );


    LOG(INFO) << "create space\n";
    node_type trans( 2 );
    trans[0]=0;
    trans[1]=2;
    Xh = functionspace_type::New( _mesh=mesh, _periodicity=Periodic<>( 2, 4, trans ) );

    Xh_vec = functionspace_vec_type::New( _mesh=mesh, _periodicity=Periodic<>( 2, 4, trans ) );

    //Xhc = functionspace_composite_type::New( _mesh=mesh, _periodicity=fusion::make_vector( Periodic<>(2, 4, trans), Periodic<>(2, 4, trans) ) );
    Xhc = functionspace_composite_type::New( _mesh=mesh, _periodicity=Periodic<>(2, 4, trans) );

    LOG(INFO) << "Xh print space info\n";
    Xh->printInfo();
    LOG(INFO) << "Xh_vec print space info\n";
    Xh_vec->printInfo();

    LOG(INFO) << "Constructor done\n";
}

template<int Dim, int Order>
void
PeriodicLaplacian<Dim, Order>::run()
{
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
    auto g = sin( M_PI*Px() )*cos( M_PI*Py() );
    auto grad_g = vec( +M_PI*cos( M_PI*Px() )*cos( M_PI*Py() ),
                       -M_PI*sin( M_PI*Px() )*sin( M_PI*Py() ) );
    auto f = 2*M_PI*M_PI*g;
#endif
    LOG(INFO) << "Number of Periodic elements: " << std::distance( Xh->dof()->beginPeriodicElements(),
                                                                   Xh->dof()->endPeriodicElements() ) << "\n";
    auto M = M_backend->newMatrix( Xh, Xh );

    form2( Xh, Xh, M ) = integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) );
    form2( Xh, Xh, M ) += integrate( _range=markedfaces( mesh, 1 ),
                                     _expr=-gradt( u )*N()*id( v )
                                     -grad( v )*N()*idt( u )
                                     + penalisation_bc*id( u )*idt( v )/hFace() );

    form2( Xh, Xh, M ) += integrate( _range=markedfaces( mesh, 3 ),
                                     _expr=-gradt( u )*N()*id( v )
                                     -grad( v )*N()*idt( u )
                                     + penalisation_bc*id( u )*idt( v )/hFace() );

    auto F = M_backend->newVector( Xh );
    form1( Xh, F ) = ( integrate( _range=elements( mesh ), _expr=f*id( v ) )
                       //+integrate( boundaryfaces( mesh ), _Q<Order+5>(), (trans(grad_g)*N())*id(v) )
        );

    if ( this->vm().count( "export-matlab" ) )
    {
        LOG(INFO) << "Exporting lhs and rhs to matlab...\n";
        M->printMatlab( "M.m" );
        F->printMatlab( "F.m" );
    }

    backend()->solve( _matrix=M, _solution=u, _rhs=F );

    LOG(INFO) << "mean(g)  = " << mean( _range=elements( mesh ), _expr=g ) << "\n";
    LOG(INFO) << "mean(u)  = " << mean( _range=elements( mesh ), _expr=idv( u ) ) << "\n";
    LOG(INFO) << "error  = " << normL2( _range=elements( mesh ), _expr=( idv( u )-g ) ) << "\n";
    auto bdy1 = mean( markedfaces( mesh,2 ), idv( u ) );
    auto bdy2 = mean( markedfaces( mesh,4 ), idv( u ) );
    LOG(INFO) << "error mean periodic  boundary 1 - 2  = " << bdy1-bdy2 << "\n";

    v = vf::project( Xh, elements( mesh ), g );

    element_type e( Xh, "e" );
    e = vf::project( Xh, elements( mesh ), idv( u )-g );

    element_type j( Xh, "Py" );
    j = vf::project( Xh, elements( mesh ), Py() );


    exportResults( u, v, e, j );


} // PeriodicLaplacian::run

template<int Dim, int Order>
void
PeriodicLaplacian<Dim, Order>::exportResults( element_type& U, element_type& V, element_type& E, element_type& F )

{
    timers["export"].first.restart();

    LOG(INFO) << "exportResults starts\n";

    exporter->step( 1. )->setMesh( U.functionSpace()->mesh() );
    exporter->step( 1. )->add( "u", U );
    exporter->step( 1. )->add( "exact", V );
    exporter->step( 1. )->add( "error", E );
    exporter->step( 1. )->add( "f", F );

    exporter->save();
    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] exportResults(): " << timers["export"].second << "\n";
} // PeriodicLaplacian::export
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    Application app;
    app.add( new PeriodicLaplacian<2,4>() );
    app.run();
}









