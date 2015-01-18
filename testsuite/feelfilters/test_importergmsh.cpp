/* -*- mode: c++: coding: utf-8 -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-06-16
 */
#include <sstream>

#include <boost/timer.hpp>
// Boost.Test
// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE mesh filter testsuite
#include <testsuite/testsuite.hpp>
#include <boost/mpl/list.hpp>

#include <feel/feelcore/feel.hpp>


using boost::unit_test::test_suite;

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelfilters/straightenmesh_impl.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/geo.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/convert2msh.hpp>
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
    BOOST_TEST_MESSAGE("Checking meshes for shape : " << shape << " using convex " << convex << " in " << Dim << "D..." );

    BOOST_CHECK_NE( nelements(markedfaces(mesh, "Neumann"),true), 0 );
    BOOST_CHECK_NE( nelements(markedfaces(mesh, "Dirichlet" ),true), 0 );
    BOOST_CHECK_EQUAL( nelements(markedfaces(mesh, "Dirichlet"),false)+nelements(markedfaces(mesh, "Neumann"),false),
                       nelements(boundaryfaces(mesh),false) );
}

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( gmshsuite )


typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2> > dim_types;

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

    BOOST_CHECK_NE(nelements(markedfaces( mesh, "letters" ),true), 0 );
    BOOST_CHECK_NE(nelements(markedfaces( mesh, "wall" ),true), 0 );
    BOOST_CHECK_NE(nelements(markedfaces( mesh, "inlet" ),true), 0 );
    BOOST_CHECK_NE(nelements(markedfaces( mesh, "outlet" ),true), 0 );
    BOOST_CHECK_NE(nelements(markedelements( mesh, "feel" ),true), 0 );
    BOOST_CHECK_EQUAL( std::distance( mesh->beginElementWithProcessId(), mesh->endElementWithProcessId() ),
                       nelements(markedelements( mesh, "feel" ),false) );
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
    if ( T::value==1 && Environment::worldComm().size()>1) return;

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

#if 0
    std::ostringstream estr;
    estr << "gmshexp-" << T::value;
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    export_ptrtype exporter( Exporter<mesh_type>::New( "gmsh", estr.str() ) );
    exporter->step( 0 )->setMesh( mesh );
    exporter->save();
    std::ostringstream fstr;
    fstr << "gmshexp-" << T::value << "-1_0.msh";
#else
    std::ostringstream fstr;
    fstr << "gmshexp-" << T::value << ".msh";
    saveGMSHMesh(_mesh=mesh,_filename=fstr.str() );
#endif


    meshimp = loadGMSHMesh( _mesh=new mesh_type,
                            _filename=fstr.str(),
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );

    BOOST_CHECK_EQUAL( nelements( elements(mesh) ),nelements( elements(meshimp) ) );
    BOOST_CHECK_EQUAL( nelements( markedfaces(mesh,"Neumann") ),nelements( markedfaces(meshimp,"Neumann") ) );
    BOOST_CHECK_EQUAL( nelements( markedfaces(mesh,"Dirichlet") ),nelements( markedfaces(meshimp,"Dirichlet") ) );

    BOOST_WARN_EQUAL( std::distance( mesh->beginFaceOnBoundary(), mesh->endFaceOnBoundary() ),
                      std::distance( meshimp->beginFaceOnBoundary(), meshimp->endFaceOnBoundary() ) );
    BOOST_CHECK_EQUAL( std::distance( mesh->beginElement(), mesh->endElement() ),
                       std::distance( meshimp->beginElement(), meshimp->endElement() ) );

    double r1 = integrate( _range=boundaryfaces( mesh ), _expr=cst( 1. ) ).evaluate()(0,0);
    double r2 = integrate( _range=boundaryfaces( meshimp ), _expr=cst( 1. ) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( std::abs(r1-r2),1e-12 );

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
BOOST_AUTO_TEST_CASE( meditimport )
{
    using namespace Feel::vf;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh,meshimp;
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=convert2msh( "Cylref.mesh" ),
                           _physical_are_elementary_regions=true,
                           _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    mesh->addMarkerName("inlet",1,2);
    mesh->addMarkerName("outlet",2,2);
    mesh->addMarkerName("wall",3,2);

#if 0
    std::ostringstream estr;
    estr << "gmshexp-" << T::value;
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    export_ptrtype exporter( Exporter<mesh_type>::New( "gmsh", estr.str() ) );
    exporter->step( 0 )->setMesh( mesh );
    exporter->save();
    std::ostringstream fstr;
    fstr << "gmshexp-" << T::value << "-1_0.msh";
#else
    std::ostringstream fstr;
    fstr << "gmshexp-" << 3 << ".msh";
    saveGMSHMesh(_mesh=mesh,_filename=fstr.str() );
#endif

    meshimp = loadGMSHMesh( _mesh=new mesh_type,
                            _filename=fstr.str(),
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );

    BOOST_CHECK_EQUAL( nelements( elements(mesh) ),nelements( elements(meshimp) ) );
    BOOST_CHECK_EQUAL( nelements( markedfaces(mesh,"inlet") ),nelements( markedfaces(meshimp,"inlet") ) );
    BOOST_CHECK_EQUAL( nelements( markedfaces(mesh,"outlet") ),nelements( markedfaces(meshimp,"outlet") ) );
    BOOST_CHECK_EQUAL( nelements( markedfaces(mesh,"wall") ),nelements( markedfaces(meshimp,"wall") ) );

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
