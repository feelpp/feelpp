/* -*- mode: c++: coding: utf-8 -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-16

  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)

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
   \file test_importergmsh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-06-16
 */
#include <sstream>

#include <boost/timer.hpp>
// Boost.Test
// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE mesh filter testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include <feel/feelcore/feel.hpp>


using boost::unit_test::test_suite;

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;


template<int Dim, template <uint16_type,uint16_type,uint16_type> class Entity = Simplex>
void
checkCreateGmshMesh( std::string const& shape, std::string const& convex = "Simplex" )
{
    typedef Mesh<Entity<Dim,1,Dim> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh;
    // simplex
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name=( boost::format( "%1%-%2%-%3%" )  % shape % convex % Dim ).str() ,
                                         _convex=convex,
                                         _usenames=true,
                                         _addmidpoint=false,
                                         _shape=shape,
                                         _h=0.5 ) );

    auto neumann = markedfaces( mesh, "Neumann" );
    BOOST_CHECK_NE( std::distance( neumann.template get<1>(), neumann.template get<2>() ), 0 );
    auto dirichlet = markedfaces( mesh, "Dirichlet" );
    BOOST_CHECK_NE( std::distance( dirichlet.template get<1>(), dirichlet.template get<2>() ), 0 );
    BOOST_CHECK_EQUAL( std::distance( neumann.template get<1>(), neumann.template get<2>() )+
                       std::distance( dirichlet.template get<1>(), dirichlet.template get<2>() ),
                       std::distance( mesh->beginFaceOnBoundary(), mesh->endFaceOnBoundary() ) );

}
BOOST_AUTO_TEST_SUITE( gmshsuite )
Feel::Environment env( boost::unit_test::framework::master_test_suite().argc,
                       boost::unit_test::framework::master_test_suite().argv );

typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( gmshsimplex, T, dim_types )
{
    checkCreateGmshMesh<T::value>( "simplex" );
}
BOOST_AUTO_TEST_CASE_TEMPLATE( gmshhypercube_simplex, T, dim_types )
{
    checkCreateGmshMesh<T::value>( "hypercube" );
}
BOOST_AUTO_TEST_CASE_TEMPLATE( gmshhypercube_hypercube, T, dim_types )
{
    checkCreateGmshMesh<T::value, Hypercube>( "hypercube", "Hypercube" );
}
BOOST_AUTO_TEST_CASE_TEMPLATE( gmshellipsoid, T, dim_types )
{
    checkCreateGmshMesh<T::value>( "ellipsoid" );
}

BOOST_AUTO_TEST_CASE( gmshgeo )
{
    typedef Mesh<Simplex<2,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh;
    // simplex
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename="feel.geo",
                                      _dim=2,
                                      _order=1,
                                      _h=0.2 ) );
    auto letters = markedfaces( mesh, "letters" );
    BOOST_CHECK_NE( std::distance( letters.template get<1>(), letters.template get<2>() ), 0 );
    auto wall = markedfaces( mesh, "wall" );
    BOOST_CHECK_NE( std::distance( wall.template get<1>(), wall.template get<2>() ), 0 );
    auto inlet = markedfaces( mesh, "inlet" );
    BOOST_CHECK_NE( std::distance( inlet.template get<1>(), inlet.template get<2>() ), 0 );
    auto outlet = markedfaces( mesh, "outlet" );
    BOOST_CHECK_NE( std::distance( outlet.template get<1>(), outlet.template get<2>() ), 0 );
    auto markedelts = markedelements( mesh, "feel" );
    BOOST_CHECK_EQUAL( std::distance( mesh->beginElement(), mesh->endElement() ),
                       std::distance( markedelts.template get<1>(), markedelts.template get<2>() ) );
}

BOOST_AUTO_TEST_CASE( gmshpartgeo )
{
    typedef Mesh<Simplex<2,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh;
    // simplex
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename="feel.geo",
                                      _dim=2,
                                      _order=1,
                                      _h=0.2 ),
                           _partitions=2 );
}

#if defined( FEELPP_HAS_TBB )
template<typename elt_iterator>
class tbb_check_mesh
{
public:
    tbb_check_mesh()
        :
        count( 0 )
    {}
    tbb_check_mesh( tbb_check_mesh& o, tbb::split )
        :
        count( o.count )
    {}
    void operator() ( const tbb::blocked_range<elt_iterator >& r )
    {
        for ( auto _elt = r.begin(); _elt != r.end(); ++_elt, ++count )
        {

        }

    }
    void join( tbb_check_mesh& other )
    {
        count += other.count;
    }
    double count;
};
BOOST_AUTO_TEST_CASE( gmshgeo_tbb )
{
    typedef Mesh<Simplex<2,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh;
    // simplex
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename="feel.geo",
                                      _dim=2,
                                      _order=1,
                                      _h=0.2 ) );
    std::vector<boost::reference_wrapper<const mesh_type::element_type> > v;
    BOOST_FOREACH( const mesh_type::element_type& i,mesh->elements() )
    v.push_back( boost::cref( i ) );
    tbb::blocked_range<decltype( v.begin() )> r( v.begin(), v.end() );
    BOOST_TEST_MESSAGE( "range size=" << r.size() << "\n" );
    tbb_check_mesh<decltype( v.begin() )> counter;
    tbb::tick_count parallel_t0 = tbb::tick_count::now();
    tbb::parallel_reduce( r, counter );
    tbb::tick_count parallel_t1 = tbb::tick_count::now();

    BOOST_CHECK_EQUAL( counter.count, mesh->numElements() );

    tbb::tick_count serial_t0 = tbb::tick_count::now();
    int count = 0;

    for ( auto _elt = mesh->beginElement(); _elt != mesh->endElement(); ++_elt, ++count )
    {

    }

    tbb::tick_count serial_t1 = tbb::tick_count::now();
    BOOST_CHECK_EQUAL( count, mesh->numElements() );

    BOOST_TEST_MESSAGE( "Serial version ran in " << ( serial_t1 - serial_t0 ).seconds() << " seconds" << "\n"
                        << "Parallel version ran in " <<  ( parallel_t1 - parallel_t0 ).seconds() << " seconds" << "\n"
                        << "Resulting in a speedup of " << ( serial_t1 - serial_t0 ).seconds() / ( parallel_t1 - parallel_t0 ).seconds() << "\n" );

}
#endif // FEELPP_HAS_TBB

BOOST_AUTO_TEST_CASE_TEMPLATE( gmshimportexport, T, dim_types )
{
    BOOST_TEST_MESSAGE( "[gmshimportexport] for dimension " << T::value << "\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh,meshimp;
    // simplex
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name=( boost::format( "simplex-%1%" )  % T::value ).str() ,
                                         _usenames=true,
                                         _addmidpoint=false,
                                         _shape="simplex",
                                         _dim=T::value,
                                         _h=0.5 ) );

    std::ostringstream estr;
    estr << "gmshexp-" << T::value;
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    export_ptrtype exporter( Exporter<mesh_type>::New( "gmsh", estr.str() ) );
    exporter->step( 0 )->setMesh( mesh );
    exporter->save();
    std::ostringstream fstr;
    fstr << "gmshexp-" << T::value << "-1_0.msh";
    meshimp = loadGMSHMesh( _mesh=new mesh_type,
                            _filename=fstr.str(),
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );

    auto elements1 = elements( mesh );
    auto elements2 = elements( meshimp );
    BOOST_CHECK_EQUAL( std::distance( elements1.template get<1>(), elements1.template get<2>() ),
                       std::distance( elements2.template get<1>(), elements2.template get<2>() ) );

    auto neumann1 = markedfaces( mesh, "Neumann" );
    auto neumann2 = markedfaces( meshimp, "Neumann" );
    BOOST_CHECK_EQUAL( std::distance( neumann1.template get<1>(), neumann1.template get<2>() ),
                       std::distance( neumann2.template get<1>(), neumann2.template get<2>() ) );
    auto dirichlet1 = markedfaces( mesh, "Dirichlet" );
    auto dirichlet2 = markedfaces( meshimp, "Dirichlet" );
    BOOST_CHECK_EQUAL( std::distance( dirichlet1.template get<1>(), dirichlet1.template get<2>() ),
                       std::distance( dirichlet2.template get<1>(), dirichlet2.template get<2>() ) );
    BOOST_WARN_EQUAL( std::distance( mesh->beginFaceOnBoundary(), mesh->endFaceOnBoundary() ),
                      std::distance( meshimp->beginFaceOnBoundary(), meshimp->endFaceOnBoundary() ) );
    BOOST_CHECK_EQUAL( std::distance( mesh->beginElement(), mesh->endElement() ),
                       std::distance( meshimp->beginElement(), meshimp->endElement() ) );
    BOOST_TEST_MESSAGE( "[gmshimportexport] for dimension " << T::value << " done.\n" );
}

/*
int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}
*/

BOOST_AUTO_TEST_CASE_TEMPLATE( meditimport, T, dim_types )
{
    using namespace Feel::vf;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh,meshimp;
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=mshconvert( "Cylref.mesh" ),
                           _physical_are_elementary_regions=true,
                           _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );

    std::ostringstream estr;
    estr << "gmshexp-" << T::value;
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    export_ptrtype exporter( Exporter<mesh_type>::New( "gmsh", estr.str() ) );
    exporter->step( 0 )->setMesh( mesh );
    exporter->save();
    std::ostringstream fstr;
    fstr << "gmshexp-" << T::value << "-1_0.msh";
    meshimp = loadGMSHMesh( _mesh=new mesh_type,
                            _filename=fstr.str(),
                            _physical_are_elementary_regions=true,
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );

    auto neumann1 = markedfaces( mesh, "Neumann" );
    auto neumann2 = markedfaces( meshimp, mesh->markerName( "Neumann" ) );
    BOOST_CHECK_EQUAL( std::distance( neumann1.template get<1>(), neumann1.template get<2>() ),
                       std::distance( neumann2.template get<1>(), neumann2.template get<2>() ) );
    auto dirichlet1 = markedfaces( mesh, "Dirichlet" );
    auto dirichlet2 = markedfaces( meshimp, mesh->markerName( "Dirichlet" ) );
    BOOST_CHECK_EQUAL( std::distance( dirichlet1.template get<1>(), dirichlet1.template get<2>() ),
                       std::distance( dirichlet2.template get<1>(), dirichlet2.template get<2>() ) );
    BOOST_CHECK_EQUAL( std::distance( mesh->beginFaceOnBoundary(), mesh->endFaceOnBoundary() ),
                       std::distance( meshimp->beginFaceOnBoundary(), meshimp->endFaceOnBoundary() ) );
    BOOST_CHECK_EQUAL( std::distance( mesh->beginElement(), mesh->endElement() ),
                       std::distance( meshimp->beginElement(), meshimp->endElement() ) );
    BOOST_CHECK_EQUAL( integrate( elements( mesh ), cst( 1. ) ).evaluate(),
                       integrate( elements( meshimp ), cst( 1. ) ).evaluate() );
    BOOST_CHECK_EQUAL( integrate( boundaryfaces( mesh ), cst( 1. ) ).evaluate() ,
                       integrate( boundaryfaces( meshimp ), cst( 1. ) ).evaluate() );
}

BOOST_AUTO_TEST_SUITE_END()
