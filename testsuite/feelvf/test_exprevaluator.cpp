#define BOOST_TEST_MODULE exprevalauator testsuite

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelvf/expressionevaluator.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( exprevaluator_suite )

BOOST_AUTO_TEST_CASE( test_1 )
{
    auto Dmu = ParameterSpace<>::New(1);
    auto muMin = Dmu->element();
    muMin << 1;
    Dmu->setMin(muMin);
    auto muMax = Dmu->element();
    muMax << 2;
    Dmu->setMax(muMax);
    auto mu = Dmu->element();
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3> > );
    auto r = elements(mesh);
    auto Vh = Pch<1>(mesh);
    auto u = Vh->element();
    u.on(_range=elements(mesh), _expr=cst(2.));
    auto ex = vec(mu(0)*cos(cst(0.)),idv(u)*cos(cst(1.)));
    auto ee = ExpressionEvaluatorNonLinear<decltype(r), decltype(ex), decltype(u)>( r, ex, mu, u);
    ee.init(5);
    ee.update(mu);
    double res0 = 0., res1 = 0., res0loc = 0., res1loc = 0.;
    for( auto const& eltWrap : elements(mesh) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        if ( elt.processId() != Environment::rank() )
            continue;

        ee.update( eltWrap );
        for ( uint16_type q = 0; q < ee.nPoints(); ++q )
        {
            res0loc += ee.weight(q)*ee.eval(q, 0);
            res1loc += ee.weight(q)*ee.eval(q, 1);
        }
    }
    mpi::all_reduce(Environment::worldComm().globalComm(), res0loc, res0, std::plus<double>());
    mpi::all_reduce(Environment::worldComm().globalComm(), res1loc, res1, std::plus<double>());
    auto s0 = integrate(_range=elements(mesh), _expr=inner(ex,vec(cst(1.),cst(0.))) ).evaluate()(0,0);
    auto s1 = integrate(_range=elements(mesh), _expr=inner(ex,vec(cst(0.),cst(1.))) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE(res0, s0, 1e-10);
    BOOST_CHECK_CLOSE(res1, s1, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
