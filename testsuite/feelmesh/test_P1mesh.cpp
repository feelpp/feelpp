
#define BOOST_TEST_MODULE test_P1mesh
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>



using namespace Feel;
using namespace Feel::vf;

namespace test_P1mesh
{


/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_P1mesh options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ;
    return desc_options.add( Feel::feel_options() );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "Test_P1mesh" ,
                     "Test_P1mesh" ,
                     "0.1",
                     "test P1mesh",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2012 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}

template <uint32_type OrderGeo>
void
test2dP1mesh()
{
    typedef Backend<double> backend_type;
    typedef Mesh<Simplex<2,OrderGeo,2> > mesh_type;

    typedef bases<Lagrange<3,Vectorial,Continuous,PointSetFekete> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;

    typedef Mesh<Simplex<2,1,2> > mesh_P1_type;
    typedef bases<Lagrange<1,Scalar,Continuous,PointSetFekete> > basis_P1_type;
    typedef FunctionSpace<mesh_P1_type, basis_P1_type> space_P1_type;

    WorldComm myWorldComm;
    auto meshSize = doption(_name="hsize");

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 0.6,0 );
    GeoTool::Circle C( meshSize,"OMEGA",x1,x2 );
    C.setMarker( _type="line",_name="Sortie",_markerAll=true );
    C.setMarker( _type="surface",_name="OmegaFluide",_markerAll=true );
    auto mesh = C.createMesh( _mesh = new mesh_type,
                              _name="test2dP1mesh_domain"+mesh_type::shape_type::name(),
                              _partitions=myWorldComm.localSize() );
    auto meshP1 = mesh->createP1mesh();
    auto XhP1 = space_P1_type::New( _mesh=meshP1 );
    auto u = XhP1->element();
    auto v = XhP1->element();

    auto mybackend = backend_type::build( soption( _name="backend" ) );
    auto A = mybackend->newMatrix(_test=XhP1, _trial=XhP1);
    auto F = mybackend->newVector(XhP1);
    form2( _test=XhP1, _trial=XhP1, _matrix=A ) =
        integrate( _range=elements( meshP1 ),
                   _expr=gradt(u)*trans(grad(v)) );
    form1( _test=XhP1, _vector=F ) =
        integrate( _range=elements( meshP1 ),
                   _expr=id(v) );
    A->close();F->close();
    form2( _test=XhP1, _trial=XhP1, _matrix=A ) +=
        on ( _range=boundaryfaces(meshP1),
             _element=u,
             _rhs=F,
             _expr=cst(0.) );
    mybackend->solve(_matrix=A,_solution=u,_rhs=F);

    auto myexporter = Exporter<mesh_P1_type>::New( Environment::vm(), "test2dP1mesh_MyExport" );
    myexporter->step(0)->setMesh( meshP1 );
    myexporter->step(0)->add( "test2dP1mesh_uP1", u );
    myexporter->save();
}



} // namespace test_P1mesh

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_P1mesh::makeAbout(), test_P1mesh::makeOptions() )
/**
 * main code
 */
BOOST_AUTO_TEST_SUITE( interp_P1mesh )

BOOST_AUTO_TEST_CASE( interp_P1mesh )
{
    using namespace Feel::vf;
    using namespace test_P1mesh;

    test_P1mesh::test2dP1mesh<2>();

}

BOOST_AUTO_TEST_SUITE_END()

