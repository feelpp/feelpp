#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_eigenmode
#include <feel/feelcore/testsuite.hpp>

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

    tic();
    auto a = form2( _test=Vh, _trial=Vh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(u))*curl(v));
    toc("assembly a");
    tic();
    auto b = form2( _test=Vh, _trial=Vh);
    b = integrate( _range=elements(mesh), _expr=trans(idt( u ))*id( v ) );
    toc("assembly b");
    tic();
    auto modes = veigs( _formA=a,
                        _formB=b,
                        _solver="krylovschur",
                        _problem="ghep",
                        _transform="shift_invert",
#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR >= 9)
                        _spectrum="target_real",
#else
                        _spectrum="smallest_magnitude",
#endif
                        _nev=15,
                        _ncv=30,
                        _verbose=1
                        );
    toc("eigen solve");
    if ( Environment::rank() == 0 )
        BOOST_TEST_MESSAGE( "number of modes : " << modes.size() );
    tic();
    auto matA = a.matrixPtr();
    auto matB = b.matrixPtr();
    for( auto const& mode: modes )
    {
        auto Av = backend()->newVector( matA->mapRowPtr() );
        auto Bv = backend()->newVector( matB->mapRowPtr() );
        matA->multVector(mode.second, *Av);
        matB->multVector(mode.second, *Bv);
        Av->close();
        Bv->close();
        Bv->scale(mode.first);
        Bv->close();
        auto AvNorm = Av->l2Norm();
        auto BvNorm = Bv->l2Norm();

        if ( Environment::isMasterRank() )
            BOOST_CHECK_CLOSE( AvNorm, BvNorm, 1e-8 );
    }
    toc("post process");
}

BOOST_AUTO_TEST_SUITE_END()
