#define BOOST_TEST_MODULE geim testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmor/geim.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/loadmesh.hpp>
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
        ("trainset-size", po::value<int>()->default_value(8), "number of parameter in the trainset")
        ;
    options.add(geim_options());

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

BOOST_AUTO_TEST_CASE( test_eim_offline )
{
    using mesh_type = Mesh<Simplex<2>>;
    using space_type = Pch_type<mesh_type, 1>;
    using geim_type = GEIM<space_type>;
    using vector_ptrtype = typename geim_type::vector_ptrtype;
    using parameterspace_type = typename geim_type::parameterspace_type;
    using parameter_type = typename geim_type::parameter_type;
    using element_type = typename geim_type::element_type;
    using linearform_type = FsFunctionalLinear<space_type>;
    using linearform_ptrtype = std::shared_ptr<linearform_type>;

    int nbSensors = ioption("nb-sensors");
    double r = doption("radius");
    int trainsetSize = ioption("trainset-size");

    auto mesh = loadMesh( _mesh=new mesh_type );
    auto Xh = Pch<1>(mesh);
    auto u = Xh->element();

    std::vector<linearform_ptrtype> sigmas;
    for( double x = 0; x < 1; x += 1./nbSensors )
    {
        for( double y = 0; y < 1; y += 1./nbSensors )
        {
            auto f = functionalLinear(_space=Xh);
            *f = integrate( _range=elements(mesh),
                            _expr=id(u)*exp(-trans(P()-vec(cst(x),cst(y)))*(P()-vec(cst(x),cst(y)))
                                            /(2*r*r)) );
            sigmas.push_back(f);
        }
    }

    auto Dmu = parameterspace_type::New(4);
    auto muMin = Dmu->element();
    muMin << 0.1, -2, -2, -2;
    Dmu->setMin(muMin);
    auto muMax = Dmu->element();
    muMax << 3, 2, 2, 2;
    Dmu->setMax(muMax);
    auto Pset = Dmu->sampling();
    Pset->randomize(trainsetSize);

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
                      a.solve(_rhs=f, _solution=u, _name="geim");
                      return u;
                  };

    geim_type geim("test_geim", sigmas, Pset, solver, uuids::nil_uuid(), true);
    geim.offline();

    int M = geim.dimension();
    auto B = geim.matrixB();

    for(int i = 0; i < M; ++i )
    {
        for(int j = 0; j < i; ++j )
            BOOST_CHECK_SMALL(B(j,i), 1e-6);
        BOOST_CHECK_SMALL(B(i,i)-1, 1e-6);
    }
    BOOST_TEST_MESSAGE("B="<< B);

    auto Pset2 = Dmu->sampling();
    Pset2->randomize(20);
    std::vector<double> errors;
    for( auto const& mu : *Pset2 )
    {
        auto phi = solver(mu);
        auto I = geim.interpolant(phi);
        auto err = normL2(_range=elements(mesh), _expr=idv(phi)-idv(I));
        errors.push_back(err);
    }
    double max = *max_element(errors.begin(), errors.end());
    BOOST_TEST_MESSAGE("offline maximum error ="<< max);
    BOOST_CHECK_SMALL( max , 1e-4 );
}

BOOST_AUTO_TEST_CASE( test_eim_online )
{
    using mesh_type = Mesh<Simplex<2>>;
    using space_type = Pch_type<mesh_type, 1>;
    using geim_type = GEIM<space_type>;
    using vectorN_type = typename geim_type::vectorN_type;
    using vector_ptrtype = typename geim_type::vector_ptrtype;
    using parameterspace_type = typename geim_type::parameterspace_type;
    using parameter_type = typename geim_type::parameter_type;
    using element_type = typename geim_type::element_type;

    int nbSensors = ioption("nb-sensors");
    double r = doption("radius");

    geim_type geim("test_geim", crb::load::fe);
    BOOST_TEST_MESSAGE("B="<< geim.matrixB());
    auto Xh = geim.space();
    auto mesh = Xh->mesh();
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

    auto Dmu = parameterspace_type::New(4);
    auto muMin = Dmu->element();
    muMin << 0.1, -2, -2, -2;
    Dmu->setMin(muMin);
    auto muMax = Dmu->element();
    muMax << 3, 2, 2, 2;
    Dmu->setMax(muMax);
    auto Pset = Dmu->sampling();
    Pset->randomize(20);

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
                      a.solve(_rhs=f, _solution=u, _name="geim");
                      return u;
                  };


    auto indices = geim.indices();
    auto M = geim.dimension();
    vectorN_type vn(M);

    std::vector<double> errors;
    for( auto const& mu : *Pset )
    {
        auto phi = solver(mu);
        for( int i = 0; i < M; ++i)
            vn(i) = inner_product(unwrap_ptr(sigmas[indices[i]]),phi);
        auto I = geim.interpolant(vn);
        auto err = normL2(_range=elements(mesh), _expr=idv(phi)-idv(I));
        errors.push_back(err);
    }
    double max = *max_element(errors.begin(), errors.end());
    BOOST_TEST_MESSAGE("online maximum error ="<< max);
    BOOST_CHECK_SMALL( max , 1e-4 );
}
BOOST_AUTO_TEST_SUITE_END()
