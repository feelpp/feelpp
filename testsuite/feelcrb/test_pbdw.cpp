#define BOOST_TEST_MODULE pbdw testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/exporter.hpp>
#include <feel/feelcrb/pbdw.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test PBDW Options" );

    options.add( feel_options() );
    options.add_options()
        ("nb-sensors", po::value<int>()->default_value(4), "nb of sensors in x direction")
        ("radius", po::value<double>()->default_value(0.025), "radius of sensors")
        ("trainset-size", po::value<int>()->default_value(10), "number of parameter in the trainset")
        ( "do-ortho", po::value<bool>()->default_value(false), "")
        ( "ortho-tol", po::value<double>()->default_value(1e-2), "")
        ;
    options.add(pbdw_options());

    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_pbdw" ,
                     "test_pbdw" ,
                     "0.1",
                     "PBDW test",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2017 Feel++ Consortium" );

    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( pbdwsuite )

BOOST_AUTO_TEST_CASE( test_pbdw_offline )
{
    using parameterspace_type = ParameterSpaceX;
    using parameter_type = typename parameterspace_type::element_type;
    using mesh_type = Mesh<Simplex<2> >;
    using space_type = Pch_type<mesh_type, 1>;
    using reducedspace_type = ReducedBasisSpace<space_type>;
    using element_type = typename space_type::element_type;
    using element_ptrtype = std::shared_ptr<element_type>;
    using sensorbase_type = SensorBase<space_type>;
    using sensorbase_ptrtype = std::shared_ptr<sensorbase_type>;
    using sensor_type = SensorGaussian<space_type>;
    using node_t = typename sensor_type::node_t;

    int nbSensors = ioption("nb-sensors");
    double r = doption("radius");
    int N = ioption("trainset-size");
    int M = nbSensors*nbSensors;

    auto mesh = loadMesh( _mesh=new mesh_type, _filename="test_geim.geo");
    auto Xh = Pch<1>(mesh);
    auto u = Xh->element();
    auto XR = std::make_shared<reducedspace_type>(Xh);

    auto Dmu = parameterspace_type::New(4);
    auto muMin = Dmu->element();
    muMin << 0.1, -2, -2, -2;
    Dmu->setMin(muMin);
    auto muMax = Dmu->element();
    muMax << 3, 2, 2, 2;
    Dmu->setMax(muMax);
    auto Pset = Dmu->sampling();
    Pset->randomize(N);

    SensorMap<space_type> sigmas;
    for( int i = 1; i < nbSensors+1; ++i )
    {
        for( int j = 1; j < nbSensors+1; ++j )
        {
            node_t n(2);
            n(0) = i*1./(nbSensors+1);
            n(1) = j*1./(nbSensors+1);
            std::string name = "sensors_"+std::to_string((i-1)*nbSensors+(j-1));
            auto s = std::make_shared<sensor_type>(Xh, n, r,name );
            sigmas[name] = s;
        }
    }

    auto solver = [mesh,Xh](parameter_type const& mu) -> element_type {
                      auto u = Xh->element();
                      auto a = form2(_test=Xh, _trial=Xh);
                      auto f = form1(_test=Xh);
                      a = integrate(_range=markedelements(mesh, "Omega1"),
                                    _expr=inner(gradt(u),grad(u)) );
                      a += integrate(_range=markedelements(mesh, "Omega2"),
                                     _expr=mu(0)*inner(gradt(u),grad(u)) );
                      f = integrate(_range=markedfaces(mesh, "left"), _expr=mu(1)*id(u));
                      f += integrate(_range=markedfaces(mesh, "right"), _expr=mu(2)*id(u));
                      a += on(_range=markedfaces(mesh, std::list<std::string>({"bottom","top"})),
                              _rhs=f, _element=u, _expr=cst(mu(3)));
                      a.solve(_rhs=f, _solution=u);
                      return u;
                  };

    for( auto const& mu : *Pset )
    {
        auto u = solver(mu);
        XR->addPrimalBasisElement(u);
    }
    if( boption("do-ortho") )
        XR->orthonormalize("L2", false, doption("ortho-tol"));

    PBDW<decay_type<decltype(XR)>> pbdw("test_pbdw", XR, sigmas);
    pbdw.offline();
    auto m = pbdw.matrix();

    Pset->randomize(20);
    std::vector<double> errors;
    for( auto const& mu : *Pset )
    {
        auto phi = solver(mu);
        auto vn = pbdw.sensors().apply(phi);
        auto I = pbdw.solution(vn);
        errors.push_back( normL2(_range=elements(mesh),_expr=idv(phi)-idv(I)) );
    }
    double max = *max_element(errors.begin(), errors.end());
    BOOST_CHECK_SMALL( max , 1e-5 );
}

BOOST_AUTO_TEST_CASE( test_pbdw_online )
{
    using parameterspace_type = ParameterSpaceX;
    using parameter_type = typename parameterspace_type::element_type;
    using mesh_type = Mesh<Simplex<2> >;
    using space_type = Pch_type<mesh_type, 1>;
    using reducedspace_type = ReducedBasisSpace<space_type>;
    using element_type = typename space_type::element_type;
    using element_ptrtype = std::shared_ptr<element_type>;
    using sensorbase_type = SensorBase<space_type>;
    using sensorbase_ptrtype = std::shared_ptr<sensorbase_type>;
    using sensor_type = SensorGaussian<space_type>;
    using node_t = typename sensor_type::node_t;

    double r = doption("radius");
    auto pbdw = PBDW<reducedspace_type>("test_pbdw", crb::load::fe);
    auto M = pbdw.dimensionM();
    auto nbSensors = std::sqrt(M);
    auto N = pbdw.dimensionN();
    auto Xh = pbdw.functionSpace();
    auto u = Xh->element();
    auto mesh = Xh->mesh();

    auto solver = [mesh,Xh](parameter_type const& mu) -> element_type {
                      auto u = Xh->element();
                      auto a = form2(_test=Xh, _trial=Xh);
                      auto f = form1(_test=Xh);
                      a = integrate(_range=markedelements(mesh, "Omega1"),
                                    _expr=inner(gradt(u),grad(u)) );
                      a += integrate(_range=markedelements(mesh, "Omega2"),
                                     _expr=mu(0)*inner(gradt(u),grad(u)) );
                      f = integrate(_range=markedfaces(mesh, "left"), _expr=mu(1)*id(u));
                      f += integrate(_range=markedfaces(mesh, "right"), _expr=mu(2)*id(u));
                      a += on(_range=markedfaces(mesh, std::list<std::string>({"bottom","top"})),
                              _rhs=f, _element=u, _expr=cst(mu(3)));
                      a.solve(_rhs=f, _solution=u);
                      return u;
                  };

    // outputs
    std::vector<vector_ptrtype> Fs(3);
    auto f1 = form1(_test=Xh);
    double area1 = integrate(_range=markedelements(mesh, "Omega1"), _expr=cst(1.0)).evaluate()(0,0);
    f1 = integrate(_range=markedelements(mesh, "Omega1"), _expr=id(u)/cst(area1) );
    Fs[0] = f1.vectorPtr();
    auto f2 = form1(_test=Xh);
    double area2 = integrate(_range=markedelements(mesh, "Omega2"), _expr=cst(1.0)).evaluate()(0,0);
    f2 = integrate(_range=markedelements(mesh, "Omega2"), _expr=id(u)/cst(area2) );
    Fs[1] = f2.vectorPtr();
    auto f3 = form1(_test=Xh);
    double area3 = integrate(_range=boundaryfaces(mesh), _expr=cst(1.0)).evaluate()(0,0);
    f3 = integrate(_range=boundaryfaces(mesh), _expr=id(u)/cst(area3) );
    Fs[2] = f3.vectorPtr();
    pbdw.setOutputs(Fs);
    auto pbdwOnline = PBDWOnline("test_pbdw");

    // parameters
    auto Dmu = parameterspace_type::New(4);
    auto muMin = Dmu->element();
    muMin << 0.1, -2, -2, -2;
    Dmu->setMin(muMin);
    auto muMax = Dmu->element();
    muMax << 3, 2, 2, 2;
    Dmu->setMax(muMax);
    auto Pset = Dmu->sampling();
    Pset->randomize(20);

    std::vector<double> errors, errorsOutput1, errorsOutput2, errorsOutput3;
    for( auto const& mu : *Pset )
    {
        auto phi = solver(mu);
        auto vn = pbdw.sensors().apply(phi);
        auto I = pbdw.solution(vn);
        errors.push_back( normL2(_range=elements(mesh),_expr=idv(phi)-idv(I)) );

        auto outs1 = pbdwOnline.outputs(vn);
        auto outs2 = pbdwOnline.outputs(vn({1,2,3,4,5,6,8,9}), {1, 2, 3, 4, 5, 6, 8, 9});
        auto outs3 = pbdwOnline.outputsWithout(vn({0,1,2,3,5,6,7,8,9,10,12,13,14,15}),
                                               std::vector<std::string>{"sensors_5","sensors_12"});
        for( int i = 0; i < Fs.size(); ++i )
        {
            errorsOutput1.push_back( std::abs(outs1(i)-inner_product(*Fs[i], phi)) );
            errorsOutput2.push_back( std::abs(outs2(i)-inner_product(*Fs[i], phi)) );
            errorsOutput3.push_back( std::abs(outs3(i)-inner_product(*Fs[i], phi)) );
        }
    }
    double max = *max_element(errors.begin(), errors.end());
    BOOST_CHECK_SMALL( max , 1e-5 );
    double maxOutput1 = *max_element(errorsOutput1.begin(), errorsOutput1.end());
    BOOST_CHECK_SMALL( maxOutput1 , 1e-5 );
    double maxOutput2 = *max_element(errorsOutput2.begin(), errorsOutput2.end());
    BOOST_CHECK_SMALL( maxOutput2 , 1e-4 );
    double maxOutput3 = *max_element(errorsOutput3.begin(), errorsOutput3.end());
    BOOST_CHECK_SMALL( maxOutput3 , 1e-5 );
}

BOOST_AUTO_TEST_SUITE_END()
