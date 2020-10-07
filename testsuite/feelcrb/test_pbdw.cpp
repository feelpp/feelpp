#define BOOST_TEST_MODULE pbdw testsuite

#include <feel/feelcore/testsuite.hpp>
#define USE_BOOST_TEST 1

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
        ;
    options.add(backend_options("pbdw"));
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

bool orthonormalize(std::vector<std::shared_ptr<Pch_element_type<Mesh<Simplex<2>>,1>>>& wn )
{
    int N = wn.size();
    if( N == 0 )
        return false;
    auto Xh = wn[0]->functionSpace();
    auto mesh = Xh->mesh();
    auto mf = form2(_test=Xh, _trial=Xh);
    mf = integrate(elements(mesh), inner(id(wn[0]),idt(wn[0])));
    auto m = mf.matrixPtr();
    for ( size_type i = 0; i < N; ++i )
    {
        auto & wni = unwrap_ptr( wn[i] );
        for ( size_type j = 0; j < i; ++j )
        {
            auto const& wnj = unwrap_ptr( wn[j] );
            double __rij_pr = m->energy( wni, wnj );
            wni.add( -__rij_pr, wnj );
        }
    }

    for ( size_type i = 0; i < N; ++i )
    {
        auto & wni = unwrap_ptr( wn[i] );
        double __rii_pr = math::sqrt( m->energy( wni, wni ) );
        wni.scale( 1./__rii_pr );
    }
    double max = std::numeric_limits<double>::lowest();
    for ( size_type i = 0; i < N; ++i )
    {
        auto & wni = unwrap_ptr( wn[i] );
        for ( size_type j = 0; j < i; ++j )
        {
            auto const& wnj = unwrap_ptr( wn[j] );
            double __rij_pr = m->energy( wni, wnj );
            if( std::abs(__rij_pr) > max )
                max = std::abs(__rij_pr);
        }
    }
    return max < 1e-2;
}

#if USE_BOOST_TEST
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

    std::vector<sensorbase_ptrtype> sigmas;
    for( int i = 1; i < nbSensors+1; ++i )
    {
        for( int j = 1; j < nbSensors+1; ++j )
        {
            node_t n(2);
            n(0) = i*1./(nbSensors+1);
            n(1) = j*1./(nbSensors+1);
            auto s = std::make_shared<sensor_type>(Xh, n, r,
                                                   "sensors_"+std::to_string((i-1)*nbSensors+(j-1)));
            sigmas.push_back(s);
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
                      a.solve(_rhs=f, _solution=u, _name="pbdw");
                      return u;
                  };

    for( auto const& mu : *Pset )
    {
        auto u = solver(mu);
        XR->addPrimalBasisElement(u);
    }
    if( boption("do-ortho") )
    {
        auto wn = XR->primalRB();
        do {} while( !orthonormalize(wn) );
    }

    PBDW pbdw("test_pbdw", XR, sigmas);
    pbdw.offline();
    auto m = pbdw.matrix();

    Pset->randomize(20);
    std::vector<double> errors;
    for( auto const& mu : *Pset )
    {
        auto phi = solver(mu);
        auto vn = vectorN_type(M);
        for(int i = 0; i < M; ++i)
            vn(i) = (*sigmas[i])(phi);
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
    auto mesh = Xh->mesh();

    std::vector<sensorbase_ptrtype> sigmas;
    for( int i = 1; i < nbSensors+1; ++i )
    {
        for( int j = 1; j < nbSensors+1; ++j )
        {
            node_t n(2);
            n(0) = i*1./(nbSensors+1);
            n(1) = j*1./(nbSensors+1);
            auto s = std::make_shared<sensor_type>(Xh, n, r,
                                                   "sensors_"+std::to_string((i-1)*nbSensors+(j-1)));
            sigmas.push_back(s);
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
                      a.solve(_rhs=f, _solution=u, _name="pbdw");
                      return u;
                  };

    auto Dmu = parameterspace_type::New(4);
    auto muMin = Dmu->element();
    muMin << 0.1, -2, -2, -2;
    Dmu->setMin(muMin);
    auto muMax = Dmu->element();
    muMax << 3, 2, 2, 2;
    Dmu->setMax(muMax);
    auto Pset = Dmu->sampling();
    Pset->randomize(20);

    std::vector<double> errors;
    for( auto const& mu : *Pset )
    {
        auto phi = solver(mu);
        auto vn = vectorN_type(M);
        for(int i = 0; i < M; ++i)
            vn(i) = (*sigmas[i])(phi);
        auto I = pbdw.solution(vn);
        errors.push_back( normL2(_range=elements(mesh),_expr=idv(phi)-idv(I)) );
    }
    double max = *max_element(errors.begin(), errors.end());
    BOOST_CHECK_SMALL( max , 1e-5 );
}

BOOST_AUTO_TEST_SUITE_END()

#else
int main(int argc, char** argv)
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

    Environment env(_argc=argc, _argv=argv,
                    _about=makeAbout(), _desc=makeOptions());

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

    std::vector<sensorbase_ptrtype> sigmas;
    for( int i = 1; i < nbSensors+1; ++i )
    {
        for( int j = 1; j < nbSensors+1; ++j )
        {
            node_t n(2);
            n(0) = i*1./(nbSensors+1);
            n(1) = j*1./(nbSensors+1);
            auto s = std::make_shared<sensor_type>(Xh, n, r,
                                                   "sensors_"+std::to_string((i-1)*nbSensors+(j-1)));
            sigmas.push_back(s);
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
                      a.solve(_rhs=f, _solution=u, _name="pbdw");
                      return u;
                  };

    for( auto const& mu : *Pset )
    {
        auto u = solver(mu);
        XR->addPrimalBasisElement(u);
    }
    if( boption("do-ortho") )
    {
        auto wn = XR->primalRB();
        do {} while( !orthonormalize(wn) );
    }

    auto e1 = exporter(_mesh=mesh,_name="test_pbdw_rb");
    int i = 0;
    for( auto const& b : XR->primalRB() )
        e1->add("b_"+std::to_string(i++), *b);
    for( auto const& s : sigmas )
        e1->add(s->name(), Xh->element(s->containerPtr()));
    e1->save();

    PBDW pbdw("test_pbdw", XR, sigmas);
    pbdw.offline();
    auto m = pbdw.matrix();

    Pset->randomize(4);
    std::vector<double> errors;
    auto e = exporter(_mesh=mesh,_name="test_pbdw");
    i = 0;
    for( auto const& mu : *Pset )
    {
        auto phi = solver(mu);
        auto vn = vectorN_type(M);
        for(int i = 0; i < M; ++i)
            vn(i) = (*sigmas[i])(phi);
        auto I = pbdw.solution(vn);
        e->step(i)->add("phi", phi);
        e->step(i)->add("I",I);
        i++;
        errors.push_back( normL2(_range=elements(mesh),_expr=idv(phi)-idv(I)) );
    }
    e->save();
    double max = *max_element(errors.begin(), errors.end());
    Feel::cout << "maximum error = " << max << std::endl;
}
#endif
