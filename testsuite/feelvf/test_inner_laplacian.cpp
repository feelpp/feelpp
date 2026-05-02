#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_inner_laplacian
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( inner_laplacian_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    auto mesh = unitSphere();

    auto Vh = Pchv<2>(mesh);
    auto u = Vh->element();
    auto v = Vh->element();
    auto w = Vh->element();
    auto f = expr<3,1>("{x^2,y^2,z^2}:x:y:z");
    auto lapf = expr<3,1>("{2,2,2}:x:y:z");

    w.on( _range=elements( mesh ), _expr=f );

    auto a = form2( _test=Vh, _trial=Vh);
    a = integrate( _range=elements( mesh ), _expr=inner(laplaciant(u),laplacian(v)), _quad=_Q<4>() );

    auto interpolation_error = normL2( _range=elements( mesh ), _expr=laplacianv(w)-lapf, _quad=_Q<4>() );
    auto assembled_energy = a.matrixPtr()->energy( w, w );
    auto expected_energy = integrate( _range=elements( mesh ), _expr=inner( lapf, lapf ), _quad=_Q<4>() ).evaluate()( 0, 0 );

    BOOST_CHECK_SMALL( interpolation_error, 1e-10 );
    BOOST_CHECK( assembled_energy > 0. );
    BOOST_CHECK_CLOSE( assembled_energy, expected_energy, 1e-9 );
}

BOOST_AUTO_TEST_SUITE_END()
