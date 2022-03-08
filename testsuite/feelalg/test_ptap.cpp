#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_eigenmode
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feeldiscr/ned1h.hpp>

using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( eigenmode_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3> >);

    auto Nh = Ned1h<0>( mesh );
    auto u = Nh->element();
    auto v = Nh->element();

    auto a = form2( _test=Nh, _trial=Nh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(u))*curl(v));
    auto matA = a.matrixPtr();
    matA->close();

    auto cTilde = backend()->newMatrix(_test=Nh, _trial=Nh);
    std::vector<uint32_type> rows( Nh->nLocalDof() );
    std::iota( rows.begin(), rows.end(), 0 );
    std::vector<uint32_type> indexesToKeep( Nh->nLocalDof() );
    std::iota( indexesToKeep.begin(), indexesToKeep.end(), 0 );
    indexesToKeep.erase( indexesToKeep.begin()+2);

    auto C = cTilde->createSubMatrix(rows, indexesToKeep);
    C->close();
    auto aHat = backend()->newMatrix();
    tic();
    backend()->PtAP( matA, C, aHat );
    toc("PtAP");
    aHat->close();
}

BOOST_AUTO_TEST_SUITE_END()
