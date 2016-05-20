
#define BOOST_TEST_MODULE test_submesh1d
#include <testsuite/testsuite.hpp>
#include <feel/feelvf/vf.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_submesh1d )

BOOST_AUTO_TEST_CASE( test_submesh11d1 )
{
    using namespace Feel;
    
    auto mesh3d = loadMesh(_mesh=new Mesh<Simplex<3>>);
    auto Xh3 = THch<1>( mesh3d );
    auto U3 = Xh3->element();
    auto u3 = U3.element<0>();
    auto p3 = U3.element<1>();
    
    auto mesh1d = createSubmesh(mesh3d, markededges(mesh3d,"centerline"));
    
    size_type nbElements = nelements( elements(mesh1d), true );
    BOOST_CHECK( nbElements >0 );

}
BOOST_AUTO_TEST_SUITE_END()

