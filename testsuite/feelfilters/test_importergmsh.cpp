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
#include <cmath>
#include <array>
#include <type_traits>

// Boost.Test
// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE mesh filter testsuite
#include <feel/feelcore/testsuite.hpp>
#include <boost/test/data/test_case.hpp>

#include <feel/feelcore/feel.hpp>


using boost::unit_test::test_suite;
namespace bdata = boost::unit_test::data;

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelfilters/straightenmesh_impl.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/geo.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/convert2msh.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;


template<int Dim, template <int,int,int> class Entity = Simplex>
void
checkCreateGmshMesh( std::string const& shape, std::string const& convex = "Simplex" )
{
    typedef Mesh<Entity<Dim,1,Dim> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

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

namespace
{
template <typename Callable>
void runForDim( int dim, Callable&& c )
{
    switch ( dim )
    {
    case 1:
        c( std::integral_constant<int, 1>{} );
        break;
    case 2:
        c( std::integral_constant<int, 2>{} );
        break;
    case 3:
        c( std::integral_constant<int, 3>{} );
        break;
    default:
        BOOST_FAIL( "Unsupported dimension " << dim );
    }
}

template <int Dim>
void checkGmshImportExport()
{
    if ( Environment::worldComm().size() > 1 )
        return;

    BOOST_TEST_MESSAGE( "[gmshimportexport] for dimension " << Dim << "\n" );
    using mesh_type = Mesh<Simplex<Dim,1>>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    mesh_ptrtype mesh, meshimp;
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name=( boost::format( "simplex-%1%" )  % Dim ).str() ,
                                         _usenames=true,
                                         _addmidpoint=false,
                                         _shape="simplex",
                                         _dim=Dim,
                                         _h=0.5 ) );

    std::ostringstream fstr;
    fstr << "gmshexp-" << Dim << ".msh";
    saveGMSHMesh( _mesh=mesh, _filename=fstr.str() );

    meshimp = loadGMSHMesh( _mesh=new mesh_type,
                            _filename=fstr.str(),
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );

    BOOST_CHECK_EQUAL( nelements( elements( mesh ) ), nelements( elements( meshimp ) ) );
    BOOST_CHECK_EQUAL( nelements( markedfaces( mesh, "Neumann" ) ), nelements( markedfaces( meshimp, "Neumann" ) ) );
    BOOST_CHECK_EQUAL( nelements( markedfaces( mesh, "Dirichlet" ) ), nelements( markedfaces( meshimp, "Dirichlet" ) ) );
    BOOST_WARN_EQUAL( nelements( boundaryfaces( mesh ) ), nelements( boundaryfaces( meshimp ) ) );
    BOOST_CHECK_EQUAL( std::distance( mesh->beginElement(), mesh->endElement() ),
                       std::distance( meshimp->beginElement(), meshimp->endElement() ) );

    double r1 = integrate( _range=boundaryfaces( mesh ), _expr=cst( 1. ) ).evaluate()( 0,0 );
    double r2 = integrate( _range=boundaryfaces( meshimp ), _expr=cst( 1. ) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( std::abs( r1-r2 ), 1e-12 );

    BOOST_TEST_MESSAGE( "[gmshimportexport] for dimension " << Dim << " done.\n" );
}

using periodic_mesh_type = Mesh<Simplex<2,1>>;
using periodic_mesh_ptrtype = std::shared_ptr<periodic_mesh_type>;

bool gmshDefaultIsV4()
{
    std::string version = FEELPP_GMSH_FORMAT_VERSION;
    return !version.empty() && version.front() == '4';
}

std::string periodicSquareGeoDescription( std::string const& mshVersion, double h = 0.2 )
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << mshVersion << ";\n"
         << "h=" << h << ";\n"
         << "Point(1) = {0,0,0,h};\n"
         << "Point(2) = {1,0,0,h};\n"
         << "Point(3) = {1,1,0,h};\n"
         << "Point(4) = {0,1,0,h};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,4};\n"
         << "Line(4) = {4,1};\n"
         << "Line Loop(5) = {1,2,3,4};\n"
         << "Plane Surface(6) = {5};\n"
         << "Physical Surface(\"Omega\") = {6};\n"
         << "Physical Line(\"Bottom\") = {1};\n"
         << "Physical Line(\"Right\") = {2};\n"
         << "Physical Line(\"Top\") = {3};\n"
         << "Physical Line(\"Left\") = {4};\n"
         << "Periodic Curve {2} = {4} Translate {1,0,0};\n"
         << "Periodic Curve {3} = {1} Translate {0,1,0};\n";
    return ostr.str();
}

periodic_mesh_ptrtype createPeriodicMeshFromGmsh( GMSH_FORMAT format, std::string const& prefix,
                                                  std::string const& mshVersion = FEELPP_GMSH_FORMAT_VERSION )
{
    Gmsh gmsh;
    gmsh.setDimension( 2 );
    gmsh.setOrder( 1 );
    gmsh.setVersion( mshVersion, format );
    gmsh.setPrefix( prefix );

    std::string fname;
    bool generated_or_modified = false;
    boost::tie( fname, generated_or_modified ) = gmsh.generate( prefix, periodicSquareGeoDescription( mshVersion ), true );
    Feel::detail::ignore_unused_variable_warning( generated_or_modified );

    return loadGMSHMesh( _mesh=new periodic_mesh_type,
                         _filename=fname,
                         _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
}

void checkPeriodicMeshData( periodic_mesh_ptrtype const& mesh )
{
    BOOST_REQUIRE( mesh );
    BOOST_CHECK( mesh->isPeriodic() );

    auto const& periodicEntities = mesh->periodicEntities();
    BOOST_CHECK_GE( periodicEntities.size(), static_cast<std::size_t>( 2 ) );

    int matchedEntities = 0;
    for ( auto const& e : periodicEntities )
    {
        if ( e.dim != 1 )
            continue;

        double tx = 0.0;
        double ty = 0.0;
        bool known = false;
        if ( e.slave == 2 && e.master == 4 )
        {
            tx = 1.0;
            ty = 0.0;
            known = true;
        }
        else if ( e.slave == 3 && e.master == 1 )
        {
            tx = 0.0;
            ty = 1.0;
            known = true;
        }

        if ( !known )
            continue;

        BOOST_CHECK( !e.correspondingVertices.empty() );
        int checkedPairs = 0;
        int mappedPairs = 0;
        for ( auto const& [slavePointId, masterPointId] : e.correspondingVertices )
        {
            auto pitSlave = mesh->pointIterator( slavePointId );
            auto pitMaster = mesh->pointIterator( masterPointId );
            BOOST_REQUIRE( pitSlave != mesh->endPoint() );
            BOOST_REQUIRE( pitMaster != mesh->endPoint() );

            auto const& pSlave = pitSlave->second;
            auto const& pMaster = pitMaster->second;
            BOOST_CHECK_SMALL( std::abs( pSlave.node()[0] - pMaster.node()[0] - tx ), 1e-11 );
            BOOST_CHECK_SMALL( std::abs( pSlave.node()[1] - pMaster.node()[1] - ty ), 1e-11 );
            if ( pSlave.id() != pMaster.id() )
            {
                ++checkedPairs;
                if ( pSlave.masterId() == pMaster.id() )
                    ++mappedPairs;
            }
        }
        BOOST_CHECK_GT( checkedPairs, 0 );
        BOOST_CHECK_GT( mappedPairs, 0 );
        ++matchedEntities;
    }
    BOOST_CHECK_EQUAL( matchedEntities, 2 );

    int periodicPointCount = 0;
    for ( auto it = mesh->beginPoint(); it != mesh->endPoint(); ++it )
    {
        if ( it->second.masterId() != it->second.id() )
            ++periodicPointCount;
    }
    BOOST_CHECK_GT( periodicPointCount, 0 );
}
}

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( gmshsuite )

BOOST_DATA_TEST_CASE( gmshsimplex, bdata::make( std::array<int,3>{ { 1,2,3 } } ), dim )
{
    BOOST_TEST_CONTEXT( "shape=simplex dim=" << dim )
    {
        runForDim( dim, []( auto d )
        {
            constexpr int Dim = decltype( d )::value;
            checkCreateGmshMesh<Dim>( "simplex" );
        } );
    }
}
BOOST_DATA_TEST_CASE( gmshhypercube_simplex, bdata::make( std::array<int,3>{ { 1,2,3 } } ), dim )
{
    BOOST_TEST_CONTEXT( "shape=hypercube simplex-convex dim=" << dim )
    {
        runForDim( dim, []( auto d )
        {
            constexpr int Dim = decltype( d )::value;
            checkCreateGmshMesh<Dim>( "hypercube" );
        } );
    }
}
BOOST_DATA_TEST_CASE( gmshhypercube_hypercube, bdata::make( std::array<int,3>{ { 1,2,3 } } ), dim )
{
    BOOST_TEST_CONTEXT( "shape=hypercube hypercube-convex dim=" << dim )
    {
        runForDim( dim, []( auto d )
        {
            constexpr int Dim = decltype( d )::value;
            checkCreateGmshMesh<Dim, Hypercube>( "hypercube", "Hypercube" );
        } );
    }
}
BOOST_DATA_TEST_CASE( gmshellipsoid, bdata::make( std::array<int,3>{ { 1,2,3 } } ), dim )
{
    BOOST_TEST_CONTEXT( "shape=ellipsoid dim=" << dim )
    {
        runForDim( dim, []( auto d )
        {
            constexpr int Dim = decltype( d )::value;
            checkCreateGmshMesh<Dim>( "ellipsoid" );
        } );
    }
}

BOOST_AUTO_TEST_CASE( gmshgeo )
{
    using mesh_type = Mesh<Simplex<2,1>>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    // simplex
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename="feel.geo",
                                      _dim=2,
                                      _order=1,
                                      _h=0.2 ) );

    BOOST_CHECK_NE(nelements(markedfaces( mesh, "letters" ),true), 0 );
    BOOST_CHECK_NE(nelements(markedfaces( mesh, "wall" ),true), 0 );
    BOOST_CHECK_NE(nelements(markedfaces( mesh, "inlet" ),true), 0 );
    BOOST_CHECK_NE(nelements(markedfaces( mesh, "outlet" ),true), 0 );
    BOOST_CHECK_NE(nelements(markedelements( mesh, "feel" ),true), 0 );
    BOOST_CHECK_EQUAL( nelements(elements( mesh ),false),
                       nelements(markedelements( mesh, "feel" ),false) );
}

BOOST_AUTO_TEST_CASE( gmshpartgeo )
{
    using mesh_type = Mesh<Simplex<2,1>>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    // simplex
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename="feel.geo",
                                      _dim=2,
                                      _order=1,
                                      _h=0.2 ),
                           _partitions=2 );
}

BOOST_AUTO_TEST_CASE( gmshperiodic_import_v4_ascii )
{
    if ( Environment::isParallel() || !gmshDefaultIsV4() )
        return;

    auto mesh = createPeriodicMeshFromGmsh( GMSH_FORMAT_ASCII, "periodic-import-v4-ascii" );
    checkPeriodicMeshData( mesh );
}

BOOST_AUTO_TEST_CASE( gmshperiodic_import_v4_binary )
{
    if ( Environment::isParallel() || !gmshDefaultIsV4() )
        return;

    auto mesh = createPeriodicMeshFromGmsh( GMSH_FORMAT_BINARY, "periodic-import-v4-binary" );
    checkPeriodicMeshData( mesh );
}

BOOST_AUTO_TEST_CASE( gmshperiodic_import_v2_ascii )
{
    if ( Environment::isParallel() )
        return;

    auto mesh = createPeriodicMeshFromGmsh( GMSH_FORMAT_ASCII, "periodic-import-v2-ascii", "2" );
    checkPeriodicMeshData( mesh );
}

BOOST_AUTO_TEST_CASE( gmshperiodic_import_v2_binary )
{
    if ( Environment::isParallel() )
        return;

    auto mesh = createPeriodicMeshFromGmsh( GMSH_FORMAT_BINARY, "periodic-import-v2-binary", "2" );
    checkPeriodicMeshData( mesh );
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
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh;
    // simplex
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename="feel.geo",
                                      _dim=2,
                                      _order=1,
                                      _h=0.2 ) );
    std::vector<boost::reference_wrapper<const mesh_type::element_type> > v;
    for( const mesh_type::element_type& i : mesh->elements( ) )
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

BOOST_DATA_TEST_CASE( gmshimportexport, bdata::make( std::array<int,3>{ { 1,2,3 } } ), dim )
{
    BOOST_TEST_CONTEXT( "gmshimportexport dim=" << dim )
    {
        runForDim( dim, []( auto d )
        {
            constexpr int Dim = decltype( d )::value;
            checkGmshImportExport<Dim>();
        } );
    }
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
BOOST_AUTO_TEST_CASE( supportedgmshmesh_import )
{
    using namespace Feel::vf;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    if ( Environment::isParallel() )
      return;
    std::vector<std::string> supported_formats{".mesh"};
#ifdef FEELPP_HAS_GMSH_HAS_MED
    supported_formats.push_back(".med");
#endif
#ifdef FEELPP_HAS_GMSH_HAS_CGNS
    supported_formats.push_back(".cgns");
#endif
    for (auto format: supported_formats)
      {
	std::string filename{ std::string{"tripod"}+format };
	if ( ! fs::exists( fs::path(Environment::findFile(filename) ) ) )
	  throw std::logic_error( std::string( "file not found ") + filename );
	mesh_ptrtype mesh,meshimp;
	mesh = loadMesh( _mesh=new mesh_type,
			 _filename=filename,
			 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
	mesh->addMarkerName("inlet",1,2);
	mesh->addMarkerName("outlet",2,2);
	mesh->addMarkerName("wall",3,2);

#if 0
	std::ostringstream estr;
	estr << "gmshexp-" << T::value;
	typedef Exporter<mesh_type> export_type;
	typedef std::shared_ptr<export_type> export_ptrtype;
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

	BOOST_CHECK_EQUAL( nelements( boundaryfaces( mesh ) ),  nelements( boundaryfaces( meshimp ) ) );
	BOOST_CHECK_EQUAL( std::distance( mesh->beginElement(), mesh->endElement() ),
			   std::distance( meshimp->beginElement(), meshimp->endElement() ) );
	BOOST_CHECK_CLOSE( integrate( _range=elements( mesh ), _expr=cst( 1. ) ).evaluate()(0,0),
                       integrate( _range=elements( meshimp ), _expr=cst( 1. ) ).evaluate()(0,0), 1e-6 );
	BOOST_CHECK_CLOSE( integrate( _range=boundaryfaces( mesh ), _expr=cst( 1. ) ).evaluate()(0,0) ,
                       integrate( _range=boundaryfaces( meshimp ), _expr=cst( 1. ) ).evaluate()(0,0), 1e-6 );

	BOOST_TEST_MESSAGE( "[supportedgmshmesh_import] mesh format: " << format << " done.\n" );
      }
}


BOOST_AUTO_TEST_SUITE_END()
