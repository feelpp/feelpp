/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-22

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
   \file resistance.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-22
 */
#include <feel/feel.hpp>




inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description resistanceoptions( "Resistance Laplacian options" );
    resistanceoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "mesh size in domain" )
    ( "penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions" )

    ( "T0", Feel::po::value<double>()->default_value( 300 ), "Temperature imposed at the left wall" )
    ( "k1", Feel::po::value<double>()->default_value( 0.2 ), "conductivity of material 1" )
    ( "k2", Feel::po::value<double>()->default_value( 2 ), "conductivity of material 2" )
    ( "conductance", Feel::po::value<double>()->default_value( 100 ), "Conductance between the domain 1 and 2(temperature discontinuity)" )
    ( "Q", Feel::po::value<double>()->default_value( 1000 ), "Heat flux" )


    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return resistanceoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "resistance" ,
                           "resistance" ,
                           "0.1",
                           "nD(n=1,2,3) Resistance Laplacian on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008-2010 Université Joseph Fourier" );

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
class ResistanceLaplacian
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
    typedef Entity<Dim, 1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Entity<1, 1,Dim> line_entity_type;
    typedef Mesh<line_entity_type> line_mesh_type;
    typedef boost::shared_ptr<line_mesh_type> line_mesh_ptrtype;

    typedef DiscontinuousInterfaces<fusion::vector<mpl::vector<mpl::int_<4>, mpl::int_<6>, mpl::int_<7> > > > discontinuity_type;
    typedef bases<Lagrange<Order, Scalar, discontinuity_type> > basis_type;
    typedef bases<Lagrange<Order-1, Vectorial> > vectorial_basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef FunctionSpace<mesh_type, vectorial_basis_type> vectorial_functionspace_type;
    typedef boost::shared_ptr<vectorial_functionspace_type> vectorial_functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar, Discontinuous> > > p0_space_type;
    typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
    typedef typename p0_space_type::element_type p0_element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /** constructor */
    ResistanceLaplacian();

    /** mesh generation */
    mesh_ptrtype createMesh();

    line_mesh_ptrtype createLine();

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
    void exportResults( p0_element_type& k, element_type& u, element_type&, element_type& );

private:

    backend_ptrtype M_backend;

    double h;
    double penalisation_bc;

    mesh_ptrtype mesh;
    line_mesh_ptrtype line_mesh;

    functionspace_ptrtype Xh;
    vectorial_functionspace_ptrtype Yh;

    export_ptrtype exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;

}; // Resistance

template<int Dim, int Order>
ResistanceLaplacian<Dim,Order>::ResistanceLaplacian()
    :
    super(),
    M_backend( backend_type::build( soption("backend") ) ),

    // Data
    h( doption("hsize") ),
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

    LOG(INFO) << "create mesh\n";
    mesh = createMesh();
    line_mesh = createLine();

    LOG(INFO) << "create space\n";
    node_type trans( 2 );
    trans[0]=0;
    trans[1]=2;
    Xh = functionspace_type::New( _mesh=mesh );
    Yh = vectorial_functionspace_type::New( _mesh=mesh );

    LOG(INFO) << "print space info\n";
    Xh->printInfo();

    LOG(INFO) << "Constructor done\n";
}
template<int Dim, int Order>
typename ResistanceLaplacian<Dim,Order>::mesh_ptrtype
ResistanceLaplacian<Dim,Order>::createMesh()
{
    timers["mesh"].first.restart();
    mesh_ptrtype mesh( new mesh_type );
    //mesh->setRenumber( false );

    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << 2 << ";\n"
         << "a=" << -1 << ";\n"
         << "b=" << 1 << ";\n"
         << "c=" << -1 << ";\n"
         << "d=" << 1 << ";\n"
         << "h=" << 0.1 << ";\n"
         << "Point(1) = {a,c,0.0,h};\n"
         << "Point(2) = {b,c,0.0,h};\n"
         << "Point(3) = {b,d,0.0,h};\n"
         << "Point(4) = {a,d,0.0,h};\n"
         << "Point(5) = {0,c,0.0,h};\n"
         << "Point(6) = {0,d,0.0,h};\n"
         << "Point(7) = {a,0,0.0,h};\n"
         << "Point(8) = {b,0,0.0,h};\n"
         << "Point(9) = {0,0,0.0,h};\n"
         << "Line(1) = {1,5};\n"
         << "Line(2) = {5,2};\n"
         << "Line(3) = {2,8};\n"
         << "Line(4) = {8,3};\n"
         << "Line(5) = {3,6};\n"
         << "Line(6) = {6,4};\n"
         << "Line(7) = {4,7};\n"
         << "Line(8) = {7,1};\n"
         << "/* discontinuity (vertical line) */\n"
         << "Line(9) = {5, 9};\n"
         << "Line(10) = {9, 6};\n"
         << "/* horizontal line through square */\n"
         << "Line(11) = {7, 9};\n"
         << "Line(12) = {9, 8};\n"
         << "\n"
         << "Line Loop(19) = {3, -12, -9, 2};\n"
         << "Plane Surface(20) = {19};\n"
         << "Line Loop(21) = {4, 5, -10, 12};\n"
         << "Plane Surface(22) = {21};\n"
         << "Line Loop(23) = {6, 7, 11, 10};\n"
         << "Plane Surface(24) = {23};\n"
         << "Line Loop(25) = {8, 1, 9, -11};\n"
         << "Plane Surface(26) = {25};\n"
         << "\n"
         << "Physical Line(\"Tflux\") = {3, 4};\n"
         << "Physical Line(\"Tfixed\") = {8, 7};\n"
         << "Physical Line(\"Tinsulated\") = {1, 2, 6, 5};\n"
         << "Physical Line(\"Tdiscontinuity\") = {10, 9};\n"
         << "Physical Line(\"Tline\") = {11, 12};\n"
         << "\n"
         << "Physical Surface(\"k1\") = {20, 22};\n"
         << "Physical Surface(\"k2\") = {24, 26};\n";

    Gmsh gmsh;
    std::string fname = gmsh.generate( "square", ostr.str()  ).template get<0>();

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    LOG(INFO) << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;
} // ResistanceLaplacian::createMesh

template<int Dim, int Order>
typename ResistanceLaplacian<Dim,Order>::line_mesh_ptrtype
ResistanceLaplacian<Dim,Order>::createLine()
{
    timers["mesh"].first.restart();
    line_mesh_ptrtype mesh( new line_mesh_type );
    //mesh->setRenumber( false );

    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << 2 << ";\n"
         << "a=" << -1 << ";\n"
         << "b=" << 1 << ";\n"
         << "h=" << h << ";\n"
         << "Point(1) = {a,0,0.0,h};\n"
         << "Point(2) = {0,0,0.0,h};\n"
         << "Point(3) = {b,0,0.0,h};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Physical Line(\"line\") = {1,2};\n";


    Gmsh gmsh;
    std::string fname = gmsh.generate( "line", ostr.str()  ).template get<0>();

    ImporterGmsh<line_mesh_type> import( fname );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    LOG(INFO) << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;
} // ResistanceLaplacian::createMesh

template<int Dim, int Order>
void
ResistanceLaplacian<Dim, Order>::run()
{
    //    int maxIter = 10.0/meshSize;
    using namespace Feel::vf;

    timers["init"].first.restart();

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    double T0= doption("T0");
    double k1= doption("k1");
    double k2= doption("k2");
    double c= doption("conductance");
    double Q= doption("Q");
    form2( Xh, Xh, _matrix=M, _init=true ) = ( integrate( markedelements( mesh,  mesh->markerName( "k1" ) ),
                                       k1*gradt( u )*trans( grad( v ) ) )+
                                       integrate( markedelements( mesh,  mesh->markerName( "k2" ) ),
                                               k2*gradt( u )*trans( grad( v ) ) ) );

    form2( Xh, Xh, _matrix=M ) += integrate( markedfaces( mesh, mesh->markerName( "Tfixed" ) ),
                                     -k1*gradt( u )*N()*id( v )
                                     -k1*grad( v )*N()*idt( u )
                                     + penalisation_bc*id( u )*idt( v )/hFace() );
    auto N21 = vec( constant( -1. ),constant( 0. ) );
    form2( Xh, Xh, _matrix=M ) += integrate( markedfaces( mesh, mesh->markerName( "Tdiscontinuity" ) ),
                                     c*( trans( jump( id( v ) ) )*N21 )*( trans( jumpt( idt( u ) ) )*N21 ) );

    M->close();

    vector_ptrtype F( M_backend->newVector( Xh ) );
    form1( Xh, _vector=F, _init=true ) = ( integrate( markedfaces( mesh, mesh->markerName( "Tflux" ) ), Q*id( v ) )+
                                   integrate( markedfaces( mesh, mesh->markerName( "Tfixed" ) ),
                                           T0*( -k1*grad( v )*N()+ penalisation_bc*id( v )/hFace() ) ) );

    F->close();

    if ( this->vm().count( "export-matlab" ) )
    {
        M->printMatlab( "M.m" );
        F->printMatlab( "F.m" );
    }

    this->solve( M, u, F );

    double meas = integrate( markedfaces( mesh, mesh->markerName( "Tdiscontinuity" ) ),constant( 1.0 ) ).evaluate()( 0,0 ) ;
    double mean_jump = integrate( markedfaces( mesh, mesh->markerName( "Tdiscontinuity" ) ),
                                  trans( jumpv( idv( u ) ) )*N21 ).evaluate()( 0,0 );
    std::cout <<  "int ([[T]]) = " << mean_jump << "\n";
    LOG(INFO) <<  "int ([[T]]) = " << mean_jump << "\n";
    std::cout <<  "mean([[T]]) = " << mean_jump/meas << "\n";
    LOG(INFO) <<  "mean([[T]]) = " << mean_jump/meas << "\n";

    p0_space_ptrtype P0h = p0_space_type::New( mesh );
    p0_element_type k( P0h, "k" );

    k = vf::project( P0h, elements( mesh ),
                     ( emarker()==mesh->markerName( "k1" ) )*k1 +
                     ( emarker()==mesh->markerName( "k2" ) )*k2 );
    std::cout << "flux = " << integrate( markedfaces( mesh, mesh->markerName( "Tdiscontinuity" ) ),
                                         leftfacev( idv( k )*gradv( u )*N() )+
                                         rightfacev( idv( k )*gradv( u )*N() ) ).evaluate()( 0,0 ) << "\n";


    exportResults( k, u, u, u );


} // ResistanceLaplacian::run

template<int Dim, int Order>
void
ResistanceLaplacian<Dim, Order>::solve( sparse_matrix_ptrtype& D,
                                        element_type& u,
                                        vector_ptrtype& F )
{
    timers["solver"].first.restart();


    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    M_backend->solve( D, D, U, F );
    u = *U;

    timers["solver"].second = timers["solver"].first.elapsed();
    LOG(INFO) << "[timer] solve: " << timers["solver"].second << "\n";
} // ResistanceLaplacian::solve


template<int Dim, int Order>
void
ResistanceLaplacian<Dim, Order>::exportResults( p0_element_type& k, element_type& U, element_type& V, element_type& E )

{
    timers["export"].first.restart();

    LOG(INFO) << "exportResults starts\n";

    exporter->step( 1. )->setMesh( U.functionSpace()->mesh() );
    exporter->step( 1. )->add( "k", k );
    exporter->step( 1. )->add( "u", U );
    //exporter->step(1.)->add( "exact", V );


    typename vectorial_functionspace_type::element_type g( Yh, "grad u" );
    g = vf::project( Yh, elements( Yh->mesh() ), trans( gradv( U ) ) );
    exporter->step( 1. )->add( "grad(u)", g );

    exporter->save();
    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] exportResults(): " << timers["export"].second << "\n";

    std::ofstream os( "profile.dat" );
    os.precision( 4 );
    os.width( 10 );
    os.setf( std::ios::right );

    for ( size_type i = 0; i < line_mesh->numPoints(); ++i )
    {
        if ( i != 1 )
            os  << line_mesh->point( i ).node()[0] << " " << U( line_mesh->point( i ).node() ) << std::endl;
    }

    // this is the last point of the 1D mesh ie y = 0.13m
    os  << line_mesh->point( 1 ).node()[0] << " " << U( line_mesh->point( 1 ).node() ) << std::endl;

} // ResistanceLaplacian::export
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );
    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 1;

    typedef Feel::ResistanceLaplacian<nDim, nOrder> laplacian_resistance_type;

    /* define and run application */
    laplacian_resistance_type laplacian_resistance;

    laplacian_resistance.run();
}
