#define BOOST_TEST_MODULE exprevalauator testsuite

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelvf/expressionevaluator.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;

BOOST_AUTO_TEST_SUITE( exprevaluator_suite )

BOOST_AUTO_TEST_CASE( test_1 )
{
    auto Dmu = ParameterSpace<>::New(1);
    auto mu = Dmu->element();
    auto mesh = loadMesh( new Mesh<Simplex<3> > );
    auto r = elements(mesh);
    auto Vh = Pch<1>(mesh);
    auto u = Vh->element();
    auto ex = cos(cst(0.));
    auto ee = ExpressionEvaluatorNonLinear<decltype(r), decltype(ex), decltype(u)>( r, ex, mu, u);
}

BOOST_AUTO_TEST_CASE( test_2 )
{
    auto Dmu = ParameterSpace<>::New(1);
    auto mu = Dmu->element();
    auto mesh = loadMesh( new Mesh<Simplex<3> > );
    auto r = elements(mesh);
    auto Vh = Pch<1>(mesh);
    auto u = Vh->element();
    auto ex = vec(cos(cst(0.)),cos(cst(1.)));
    auto ee = ExpressionEvaluatorParam<decltype(r), decltype(ex)>( r, ex, mu);
}

BOOST_AUTO_TEST_CASE( test_3 )
{
    auto Dmu = ParameterSpace<>::New(1);
    auto mu = Dmu->element();
    auto mesh = loadMesh( new Mesh<Simplex<3> > );
    auto r = elements(mesh);
    auto Vh = Pch<1>(mesh);
    auto u = Vh->element();
    auto ex = vec(cos(cst(0.)),cos(cst(1.)));
    auto ee = ExpressionEvaluatorNonLinear<decltype(r), decltype(ex), decltype(u)>( r, ex, mu, u);
}

BOOST_AUTO_TEST_SUITE_END()
