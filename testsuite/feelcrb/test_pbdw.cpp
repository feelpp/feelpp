#define BOOST_TEST_MODULE pbdw testsuite

//#include <feel/feelcore/testsuite.hpp>
#define USE_BOOST_TEST 0

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
        ("nb-sensors", po::value<int>()->default_value(5), "nb of sensors in x direction")
        ("radius", po::value<double>()->default_value(0.2), "radius of sensors")
        ("trainset-size", po::value<int>()->default_value(100), "number of parameter in the trainset")
        ("do-ortho", po::value<bool>()->default_value(false), "")
        ( "geim.dimension-max", po::value<int>()->default_value(10), "maximum number of basis" )
        ( "geim.tolerance", po::value<double>()->default_value(1e-10), "tolerance" )
        ( "geim.rebuild-database", po::value<bool>()->default_value(true), "rebuild the database" )
        ( "geim.db.load", po::value<int>()->default_value(2), "=0 use db.filename, =1 use last DB created =2 use last DB modified =3 use db.id =4 create new db" )
        ( "geim.db.filename", po::value<std::string>()->default_value(""), "path to the db when db.load or db.update = 0" )
        ( "geim.db.id", po::value<std::string>()->default_value(""), "id of the db when db.load or db.update = 3" )
        ;
    options.add(backend_options("pbdw"));

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

void orthonormalize(std::vector<std::shared_ptr<Pch_element_type<Mesh<Simplex<2>>,1>>>& wn )
{
    int N = wn.size();
    for ( size_type i = 0; i < N; ++i )
    {
        auto & wni = unwrap_ptr( wn[i] );
        for ( size_type j = 0; j < i; ++j )
        {
            auto const& wnj = unwrap_ptr( wn[j] );
            double __rij_pr = inner_product( wni, wnj );
            wni.add( -__rij_pr, wnj );
        }
    }

    for ( size_type i = 0; i < N; ++i )
    {
        auto & wni = unwrap_ptr( wn[i] );
        double __rii_pr = math::sqrt( inner_product( wni, wni ) );
        Feel::cout << "norm " << i << " = " << __rii_pr << std::endl;
        wni.scale( 1./__rii_pr );
    }
}

#if USE_BOOST_TEST
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( pbdwsuite )

BOOST_AUTO_TEST_CASE( test_pbdw_offline )
{
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

    std::vector<vector_ptrtype> sigmas;
    for( double x = 0; x < 1; x += 1./nbSensors )
    {
        for( double y = 0; y < 1; y += 1./nbSensors )
        {
            auto f = form1( _test=Xh );
            auto ex = exp(-trans(P()-vec(cst(x),cst(y)))*(P()-vec(cst(x),cst(y)))/(2*r*r));
            auto n = integrate(_range=elements(mesh),_expr=ex).evaluate()(0,0);
            Feel::cout << "sigma("<<x<<","<<y<<") = " << n << std::endl;
            f += integrate( _range=elements(mesh),
                            _expr=id(u)*ex );
            sigmas.push_back(f.vectorPtr());
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

    auto e1 = exporter(_mesh=mesh,_name="test_geim_rb");
    int i = 0;
    for( auto const& mu : *Pset )
    {
        auto u = solver(mu);
        e1->step(i++)->add("b",u);
        XR->addPrimalBasisElement(u);
    }
    e1->save();
    if( boption("do-ortho") )
    {
        auto wn = XR->primalRB();
        orthonormalize(wn);
    }

    PBDW pbdw("test_pbdw", XR, sigmas);
    pbdw.offline();
    auto m = pbdw.matrix();
    Feel::cout << "matrix =\n" << m << std::endl;

    Pset->randomize(20);
    std::vector<double> errors;
    auto e = exporter(_mesh=mesh,_name="test_pbdw");
    i = 0;
    for( auto const& mu : *Pset )
    {
        auto phi = solver(mu);
        auto vn = vectorN_type(M);
        for(int i = 0; i < M; ++i)
            vn(i) = inner_product(unwrap_ptr(sigmas[i]),phi);
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
