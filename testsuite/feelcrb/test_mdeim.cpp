#define BOOST_TEST_MODULE mdeim testsuite

#include <testsuite/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelcrb/deim.hpp>
#include <feel/feelcrb/parameterspace.hpp>

using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test MDEIM Options" );

    options.add( feel_options() )
        .add(deimOptions());
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_mdeim" ,
                     "test_mdeim" ,
                     "0.1",
                     "MDEIM test",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2017 Feel++ Consortium" );

    return about;
}

class MDeimTest :
    public boost::enable_shared_from_this<MDeimTest>
{
public :
    typedef MDeimTest self_type;

    typedef double value_type;
    typedef Mesh<Simplex<2> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    typedef ParameterSpace<2> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;

    typedef MDEIM<self_type> mdeim_type;
    typedef boost::shared_ptr<mdeim_type> mdeim_ptrtype;

    MDeimTest() :
        Dmu( parameterspace_type::New(2) )
    {
        auto mesh = loadMesh( _mesh=new mesh_type, _filename="test_deim.geo");
        auto Xh = Pch<2>( mesh );
        auto Yh = Pch<3>( mesh );
        auto u = Xh->element();
        auto v = Yh->element();
        M = backend()->newMatrix( _test=Xh, _trial=Yh );
        M1 = backend()->newMatrix( _test=Xh, _trial=Yh );
        M2 = backend()->newMatrix( _test=Xh, _trial=Yh );
        M3 = backend()->newMatrix( _test=Xh, _trial=Yh );
        M4 = backend()->newMatrix( _test=Xh, _trial=Yh );

        auto f1 = form2( _test=Xh, _trial=Yh, _matrix=M1 );
        auto f2 = form2( _test=Xh, _trial=Yh, _matrix=M2 );
        auto f3 = form2( _test=Xh, _trial=Yh, _matrix=M3 );
        auto f4 = form2( _test=Xh, _trial=Yh, _matrix=M4 );

        f1 = integrate( markedelements(mesh,"Omega1"), id(u)*idt(v) );
        f2 = integrate( markedelements(mesh,"Omega2"), grad(u)*trans(gradt(v)) );
        f3 = integrate( markedfaces(mesh,"Gamma2"), grad(u)*N()*idt(v) );
        f4 = integrate( markedfaces(mesh,"Gamma1"), id(u)*idt(v) );
    }

    void run()
    {
        auto mu_min = Dmu->element();
        mu_min << 0, -10;
        auto mu_max = Dmu->element();
        mu_max << 5, 15;
        Dmu->setMin( mu_min );
        Dmu->setMax( mu_max );

        auto Pset = Dmu->sampling();
        std::vector<size_type> Ne(2);
        Ne[0] = 10;
        Ne[1] = 10;
        Pset->equidistributeProduct( Ne , true , "mdeim_test_sampling" );

        mdeim_ptrtype M_deim ( new mdeim_type( Dmu, Pset ) );
        M_deim->assemble = boost::bind( &MDeimTest::assemble, boost::ref(*this), _1  );

        M_deim->run();
        int m = M_deim->size();

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "Number of mode in DEIM : " << m );
        BOOST_CHECK( m==4 );


        Ne[0] = 100;
        Ne[1] = 100;
        Pset->equidistributeProduct( Ne , true , "mdeim_test_sampling" );

        auto base = M_deim->q();
        for ( auto const& mu : *Pset )
        {
            auto coeff = M_deim->beta(mu);

            assemble(mu);

            for ( int i=0; i<m; i++ )
                M->addMatrix( -coeff(i), base[i] );

            double norm = M->linftyNorm();
            BOOST_CHECK_SMALL( norm, 1e-9 );
        }

    }

    sparse_matrix_ptrtype assemble( parameter_type mu)
    {
        M->zero();
        M->addMatrix( mu[0], M1 );
        M->addMatrix( mu[1], M2 );
        M->addMatrix( mu[1]*mu[0] + mu[0]*mu[0], M3 );
        M->addMatrix( mu[1]*mu[1] + 37*mu[0], M4 );
        return M;
    }

private :
    sparse_matrix_ptrtype M, M1, M2, M3, M4;
    parameterspace_ptrtype Dmu;

};

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( mdeim_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    boost::shared_ptr<MDeimTest> m( new MDeimTest );
    m->run();
}

BOOST_AUTO_TEST_SUITE_END()
