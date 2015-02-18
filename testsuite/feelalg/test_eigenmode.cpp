#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_eigenmode
#include <testsuite/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelalg/solvereigen.hpp>

using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( eigenmode_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    auto mesh = unitSphere();

    auto Vh = Ned1h<0>(mesh);
    auto u = Vh->element();
    auto v = Vh->element();

    auto a = form2( _test=Vh, _trial=Vh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(u))*curl(v));
    auto matA = a.matrixPtr();
    matA->close();

    auto b = form2( _test=Vh, _trial=Vh);
    b = integrate( elements(mesh), trans(idt( u ))*id( v ) );
    auto matB = b.matrixPtr();
    matB->close();

    auto modes = eigs( _matrixA=matA,
                       _matrixB=matB,
                       _solver=KRYLOVSCHUR,
                       _problem=GHEP,
                       _transform=SINVERT,
                       _spectrum=SMALLEST_MAGNITUDE,
                       _nev=15,
                       _ncv=30
                       );

    if ( Environment::rank() == 0 )
        BOOST_TEST_MESSAGE( "number of modes : " << modes.size() );
    for( auto const& mode: modes )
    {
        auto Av = backend()->newVector( matA->mapRowPtr());
        auto Bv = backend()->newVector( matB->mapRowPtr());
        matA->multVector(boost::get<2>(mode.second), Av);
        matB->multVector(boost::get<2>(mode.second), Bv);
        Bv->scale(boost::get<0>(mode.second));
        auto AvNorm = Av->l2Norm();
        auto BvNorm = Bv->l2Norm();

        if ( Environment::rank() == 0 )
            BOOST_CHECK_CLOSE( AvNorm, BvNorm, 1e-8 );
    }
}

BOOST_AUTO_TEST_SUITE_END()
