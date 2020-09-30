#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE geim testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/exporter.hpp>
#include <feel/feelcrb/geim.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test GEIM Options" );

    options.add( feel_options() );
    options.add_options()
        ("nb-sensors", po::value<int>()->default_value(5), "nb of sensors in x direction")
        ("radius", po::value<double>()->default_value(0.2), "radius of sensors")
        ("trainset-size", po::value<int>()->default_value(100), "number of parameter in the trainset")
        ("geim.max-dimension", po::value<int>()->default_value(10), "maximum number of basis")
        ("geim.tolerance", po::value<double>()->default_value(1e-5), "tolerance")
        ;
    options.add(backend_options("geim"));

    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_geim" ,
                     "test_geim" ,
                     "0.1",
                     "GEIM test",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2017 Feel++ Consortium" );

    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( geimsuite )

BOOST_AUTO_TEST_CASE( test_eim1 )
{
    using mesh_type = Mesh<Simplex<2>>;
    using space_type = Pch_type<mesh_type, 1>;
    using geim_type = GEIM<space_type>;
    using vector_ptrtype = typename geim_type::vector_ptrtype;
    using parameterspace_type = typename geim_type::parameterspace_type;
    using parameter_type = typename geim_type::parameter_type;
    using element_type = typename geim_type::element_type;

    int nbSensors = ioption("nb-sensors");
    double r = doption("radius");
    int trainsetSize = ioption("trainset-size");

    auto mesh = unitSquare();
    auto Xh = Pch<1>(mesh);
    auto u = Xh->element();

    std::vector<vector_ptrtype> sigmas;
    for( double x = 0; x < 1; x += 1./nbSensors )
    {
        for( double y = 0; y < 1; y += 1./nbSensors )
        {
            auto f = form1( _test=Xh );
            f += integrate( _range=elements(mesh),
                            _expr=id(u)*exp(-trans(P()-vec(cst(x),cst(y)))*(P()-vec(cst(x),cst(y)))
                                            /(2*r*r)) );
            sigmas.push_back(f.vectorPtr());
        }
    }

    auto Dmu = parameterspace_type::New(2);
    auto muMin = Dmu->element();
    muMin << -2, -2;
    Dmu->setMin(muMin);
    auto muMax = Dmu->element();
    muMax << 2, 2;
    Dmu->setMax(muMax);
    auto Pset = Dmu->sampling();
    Pset->randomize(trainsetSize, true, "toto", false);

    auto solver = [mesh,Xh](parameter_type const& mu) -> element_type {
                      auto u = Xh->element();
                      auto a = form2(_test=Xh, _trial=Xh);
                      auto f = form1(_test=Xh);
                      a = integrate(_range=elements(mesh), _expr=inner(gradt(u),grad(u)) );
                      a += on(_range=boundaryfaces(mesh), _rhs=f, _element=u,
                              _expr=mu(0)*Px()*Px()*Py()+mu(1)*Py()*Py()*Px());
                      a.solve(_rhs=f, _solution=u, _name="geim");
                      return u;
                  };

    geim_type geim = geim_type(sigmas, Pset, solver);
    geim.offline();
    int M = geim.size();
    auto B = geim.matrixB();
    for(int i = 0; i < M; ++i )
    {
        for(int j = 0; j < i; ++j )
            BOOST_CHECK_SMALL(B(j,i), 1e-12);
        BOOST_CHECK_SMALL(B(i,i)-1, 1e-12);
    }
    Feel::cout << B << std::endl;

    auto Pset2 = Dmu->sampling();
    Pset2->randomize(20, true, "", false);

    std::vector<double> errors(Pset2->size());
    for( auto const& mu : *Pset2 )
    {
        auto phi = solver(mu);
        auto I = geim.interpolant(phi);
        auto err = normL2(_range=elements(mesh), _expr=idv(phi)-idv(I));
        errors.push_back(err);
    }
    double max = *max_element(errors.begin(), errors.end());
    Feel::cout << "maximum error = " << max << std::endl;
    BOOST_CHECK_SMALL( max , 1e-8 );
}


BOOST_AUTO_TEST_SUITE_END()
