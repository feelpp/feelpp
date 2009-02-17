/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-01-22

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file resistance.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-01-22
 */
#include <fstream>

#include <life/options.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/region.hpp>
#include <life/lifediscr/operatorlinear.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>


#include <life/lifevf/vf.hpp>




inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description resistanceoptions("Resistance Laplacian options");
    resistanceoptions.add_options()
        ("hsize", Life::po::value<double>()->default_value( 0.1 ), "mesh size in domain")
        ("penalbc", Life::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions")

        ("T0", Life::po::value<double>()->default_value( 300 ), "Temperature imposed at the left wall")
        ("k1", Life::po::value<double>()->default_value( 0.2 ), "conductivity of material 1")
        ("k2", Life::po::value<double>()->default_value( 2 ), "conductivity of material 2")
        ("c", Life::po::value<double>()->default_value( 100 ), "Conductance between the domain 1 and 2(temperature discontinuity)")
        ("Q", Life::po::value<double>()->default_value( 1000 ), "Heat flux")


        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return resistanceoptions.add( Life::life_options() );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "resistance" ,
                           "resistance" ,
                           "0.1",
                           "nD(n=1,2,3) Resistance Laplacian on simplices or simplex products",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2008-2009 Université Joseph Fourier");

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
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Entity<1, 1,Dim> line_entity_type;
    typedef Mesh<GeoEntity<line_entity_type> > line_mesh_type;
    typedef boost::shared_ptr<line_mesh_type> line_mesh_ptrtype;

    typedef bases<fem::Lagrange<Dim, Order, Scalar, Continuous, double, Entity> > basis_type;
    typedef bases<fem::Lagrange<Dim, Order-1, Vectorial, Discontinuous, double, Entity> > vectorial_basis_type;
    typedef DiscontinuousInterfaces<fusion::vector<mpl::vector<mpl::int_<4>, mpl::int_<5>, mpl::int_<6> > > >  discontinuity_type;
    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, discontinuity_type> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef FunctionSpace<mesh_type, vectorial_basis_type> vectorial_functionspace_type;
    typedef boost::shared_ptr<vectorial_functionspace_type> vectorial_functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;

    typedef FunctionSpace<mesh_type, fusion::vector<fem::Lagrange<Dim, 0, Scalar, Discontinuous> > > p0_space_type;
    typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
    typedef typename p0_space_type::element_type p0_element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef typename export_type::timeset_type timeset_type;

    /** constructor */
    ResistanceLaplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od );

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
    typename export_type::timeset_ptrtype timeSet;

    std::map<std::string,std::pair<boost::timer,double> > timers;

}; // Resistance

template<int Dim, int Order>
ResistanceLaplacian<Dim,Order>::ResistanceLaplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),

    // Data
    h( this->vm()["hsize"].template as<double>() ),
    penalisation_bc( this->vm()["penalbc"].template as<value_type>() ),

    // spaces
    Xh(),

    // export
    exporter( Exporter<mesh_type>::New( this->vm()["exporter"].template as<std::string>() )->setOptions( this->vm() ) ),
    timeSet( new timeset_type( "resistance" ) ),

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

    timeSet->setTimeIncrement( 1.0 );
    exporter->addTimeSet( timeSet );
    exporter->setPrefix( "resistance" );

    Log() << "create mesh\n";
    mesh = createMesh();
    line_mesh = createLine();

    Log() << "create space\n";
    node_type trans(2);
    trans[0]=0;
    trans[1]=2;
    Xh = functionspace_type::New( mesh, MESH_COMPONENTS_DEFAULTS );
    Yh = vectorial_functionspace_type::New( mesh, MESH_COMPONENTS_DEFAULTS );

    Log() << "print space info\n";
    Xh->printInfo();

    Log() << "Constructor done\n";
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
         << "h=" << h << ";\n"
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
         << "Line(9) = {7,0};\n"
         << "Line(9) = {0,8};\n"
         << "Line(10) = {5,0};\n"
         << "Line(11) = {0,6};\n"
         << "Line Loop(8) = {1,2,7,6};\n"
         << "Line Loop(9) = {3,4,5,-7};\n"
         << "Plane Surface(10) = {8};\n"
         << "Plane Surface(11) = {9};\n"
         << "Physical Line(\"insulated\") = {2,3,5,6};\n"
         << "Physical Line(\"Tfixed\") = {1};\n"
         << "Physical Line(\"flux\") = {4};\n"
         << "Physical Line(\"discontinuity\") = {7};\n"
         << "Physical Surface(\"k1\") = {10};\n"
         << "Physical Surface(\"k2\") = {11};\n";


    Gmsh gmsh;
    std::string fname = gmsh.generate( "square", ostr.str()  );

    ImporterGmsh<mesh_type> import( fname );
    import.setVersion( "2.0" );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    Log() << "[timer] createMesh(): " << timers["mesh"].second << "\n";
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
    std::string fname = gmsh.generate( "line", ostr.str()  );

    ImporterGmsh<line_mesh_type> import( fname );
    import.setVersion( "2.0" );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    Log() << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;
} // ResistanceLaplacian::createMesh

template<int Dim, int Order>
void
ResistanceLaplacian<Dim, Order>::run()
{
    //    int maxIter = 10.0/meshSize;
    using namespace Life::vf;

    timers["init"].first.restart();

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    double T0= this->vm()["T0"].template as<double>();
    double k1= this->vm()["k1"].template as<double>();
    double k2= this->vm()["k2"].template as<double>();
    double c= this->vm()["c"].template as<double>();
    double Q= this->vm()["Q"].template as<double>();
    form2( Xh, Xh, M, _init=true ) = ( integrate( markedelements( mesh,  mesh->markerName( "k1" ) ), _Q<2*(Order-1)>(),
                                                  k1*gradt(u)*trans(grad(v)) )+
                                       integrate( markedelements( mesh,  mesh->markerName( "k2" ) ), _Q<2*(Order-1)>(),
                                                  k2*gradt(u)*trans(grad(v)) ) );

    form2( Xh, Xh, M ) += integrate( markedfaces( mesh, mesh->markerName( "Tfixed" ) ), _Q<2*(Order-1)>(),
                                     -k1*gradt(u)*N()*id(v)
                                     -k1*grad(v)*N()*idt(u)
                                     + penalisation_bc*id(u)*idt(v)/hFace() );
    AUTO(N21,vec(constant(-1.),constant(0.)) );
    form2( Xh, Xh, M ) += integrate( markedfaces( mesh, mesh->markerName( "discontinuity" ) ), _Q<2*(Order)>(),
                                     c*(trans(jump(id(v)))*N21)*(trans(jumpt( idt( u ) ))*N21) );

    M->close();

    vector_ptrtype F( M_backend->newVector( Xh ) );
    form1( Xh, F, _init=true ) = ( integrate( markedfaces( mesh, mesh->markerName( "flux" ) ), _Q<Order+5>(), Q*id(v) )+
                                   integrate( markedfaces( mesh, mesh->markerName( "Tfixed" ) ), _Q<Order+5>(),
                                              T0*(-k1*grad(v)*N()+ penalisation_bc*id(v)/hFace() ) ) );

    F->close();

    if ( this->vm().count( "export-matlab" ) )
        {
            M->printMatlab( "M.m" );
            F->printMatlab( "F.m" );
        }

    this->solve( M, u, F );

    double meas = integrate( markedfaces( mesh, mesh->markerName( "discontinuity" ) ), _Q<0>(),constant(1.0)).evaluate()(0,0) ;
    double mean_jump = integrate( markedfaces( mesh, mesh->markerName( "discontinuity" ) ), _Q<Order>(),
                                  trans(jumpv(idv(u)))*N21).evaluate()(0,0);
    std::cout <<  "int ([[T]]) = " << mean_jump << "\n";
    Log() <<  "int ([[T]]) = " << mean_jump << "\n";
    std::cout <<  "mean([[T]]) = " << mean_jump/meas << "\n";
    Log() <<  "mean([[T]]) = " << mean_jump/meas << "\n";

    p0_space_ptrtype P0h = p0_space_type::New( mesh );
    p0_element_type k( P0h, "k" );

    k = vf::project( P0h, elements( mesh ),
                     (emarker()==mesh->markerName( "k1" ))*k1 +
                     (emarker()==mesh->markerName( "k2" ))*k2 );
    std::cout << "flux = " << integrate( markedfaces( mesh, mesh->markerName( "discontinuity" ) ), _Q<Order>(),
                                         jumpv(idv(k)*gradv(u))).evaluate()(0,0) << "\n";


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
    Log() << "[timer] solve: " << timers["solver"].second << "\n";
} // ResistanceLaplacian::solve


template<int Dim, int Order>
void
ResistanceLaplacian<Dim, Order>::exportResults( p0_element_type& k, element_type& U, element_type& V, element_type& E )

{
    timers["export"].first.restart();

    Log() << "exportResults starts\n";
    typename timeset_type::step_ptrtype timeStep = timeSet->step( 0.0 );
    timeStep->setMesh( U.functionSpace()->mesh() );
    timeStep->add( "k", k );
    timeStep->add( "u", U );
    //timeStep->add( "exact", V );


    typename vectorial_functionspace_type::element_type g( Yh, "grad u" );
    g = vf::project( Yh, elements( Yh->mesh() ), trans(gradv(U)) );
    timeStep->add( "grad(u)", g );

    exporter->save();
    timers["export"].second = timers["export"].first.elapsed();
    Log() << "[timer] exportResults(): " << timers["export"].second << "\n";

    std::ofstream os( "profile.dat" );
    os.precision( 4 );
    os.width(10 );
    os.setf( std::ios::right );
    for( size_type i = 0; i < line_mesh->numPoints(); ++i )
        {
            if ( i != 1 )
                os  << line_mesh->point( i ).node()[0] << " " << U( line_mesh->point( i ).node() ) << std::endl;
        }
    // this is the last point of the 1D mesh ie y = 0.13m
    os  << line_mesh->point( 1 ).node()[0] << " " << U( line_mesh->point( 1 ).node() ) << std::endl;

} // ResistanceLaplacian::export
} // Life




int
main( int argc, char** argv )
{
    using namespace Life;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 1;

    typedef Life::ResistanceLaplacian<nDim, nOrder> laplacian_resistance_type;

    /* define and run application */
    laplacian_resistance_type laplacian_resistance( argc, argv, makeAbout(), makeOptions() );

    laplacian_resistance.run();
}










