#define BOOST_TEST_MODULE eq testsuite

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/empiricalquadrature.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test EQ Options" );

    options.add( feel_options() )
        .add(eq_options());
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_eq" ,
                     "test_eq" ,
                     "0.1",
                     "EQ test",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2017 Feel++ Consortium" );

    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( eq_suite )

BOOST_AUTO_TEST_CASE( test_1 )
{
    auto Dmu = ParameterSpace<>::New(0);
    Dmu->setDimension(2);
    auto mu_min = Dmu->element();
    mu_min << -4, 4;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 4, 10;
    Dmu->setMax( mu_max );
    Dmu->setParameterName(0, "p1");
    Dmu->setParameterName(1, "p2");

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3> > );
    auto r = elements(mesh);

    auto eq = EmpiricalQuadrature(r,Dmu);
    auto mu = Dmu->element();

    auto g1 = cst(2.)*sin(cst_ref(mu(0)))*cos(cst_ref(mu(1)));
    eq.addExpression(g1,mu);

    auto g2 = cst(3.)*cos(cst_ref(mu(1)))*sin(cst_ref(mu(0)))*P();
    eq.addExpression(g2,mu,1);

    auto Xh = Pch<1>(mesh, r);
    auto u = Xh->element();
    using fn_t = std::function<bool(decltype(mu) const&, decltype(Xh->element())&)>;
    fn_t lambda = [&Xh,&mesh,&r](decltype(mu) const& mu, decltype(u)& v) mutable -> bool {
                      auto lhs = form2(_test=Xh,_trial=Xh);
                      auto rhs = form1(_test=Xh);
                      lhs = integrate(_range=r, _expr=inner(gradt(v),grad(v)) );
                      rhs = integrate(_range=r, _expr=cst(mu(0))*cst(mu(1))*id(v));
                      lhs += on(_rhs=rhs, _element=v, _range=boundaryfaces(mesh), _expr=cst(0.));
                      auto r = lhs.solve(_rhs=rhs, _solution=v, _rebuild=true);
                      return r.isConverged();
                  };
    auto g3 = idv(u)*cst(2.)*inner(P());
    eq.addExpression(g3, mu, u, lambda );

    int err = eq.offline();

    BOOST_CHECK( err == 0 );

    if( !err )
    {
        auto trainset = Dmu->sampling();
        trainset->randomize(10);

        for(auto const& m : *trainset )
        {
            mu = m;
            auto i11 = eq.evaluate(m, 0);
            auto i12 = integrate(_range=r, _expr=g1).evaluate()(0,0);
            BOOST_CHECK_SMALL( i11 - i12 , 1e-8 );

            auto i21 = eq.evaluate(m, 1);
            auto i22 = integrate(_range=r, _expr=g2).evaluate()(1,0);
            BOOST_CHECK_SMALL( i21 - i22 , 1e-8 );

            auto i31 = eq.evaluate(m, 2);
            auto i32 = integrate(_range=r, _expr=g3).evaluate()(0,0);
            BOOST_CHECK_SMALL( i31 - i32 , 1e-8 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
