#define BOOST_TEST_MODULE mdeim testsuite

#include <testsuite.hpp>

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
        .add(deimOptions())
        .add(crbSEROptions());
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

    typedef Pch_type<mesh_type,2> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    MDeimTest() :
        M_backend(backend()),
        Dmu( parameterspace_type::New(2) )
    {
        auto mesh = loadMesh( _mesh=new mesh_type, _filename="test_deim.geo");
        Xh = Pch<2>( mesh );

        setFunctionSpaces(Xh);
    }

    void setModelOnlineDeim( std::string name )
    {
        M_backend = backend( _name=name, _kind="eigen_dense", _worldcomm=Environment::worldCommSeq() );
    }

    void setFunctionSpaces( space_ptrtype Rh )
    {
        Xh = Rh;
         M = M_backend->newMatrix(Xh,Xh);
    }

    uuids::uuid uuid() const { return boost::uuids::nil_uuid(); }
    parameterspace_ptrtype parameterSpace() { return Dmu;}

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

        Environment::setOptionValue("deim.rebuild-db", true );
        auto M_deim = mdeim( _model=this->shared_from_this(),
                             _sampling=Pset );

        M_deim->run();
        int m = M_deim->size();

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "Number of mode in DEIM : " << m );
        BOOST_CHECK( m==4 );

        Pset->randomize( 100 , true , "mdeim_test_sampling" );

        auto base = M_deim->q();
        for ( auto const& mu : *Pset )
        {
            auto coeff = M_deim->beta(mu);

            auto M = assembleForMDEIM(mu);

            for ( int i=0; i<m; i++ )
                M->addMatrix( -coeff(i), base[i] );

            double norm = M->linftyNorm();
            BOOST_CHECK_SMALL( norm, 1e-9 );
        }

    }

    sparse_matrix_ptrtype assembleForMDEIM( parameter_type mu)
    {
        auto mesh = Xh->mesh();
        auto u = Xh->element();
        auto v = Xh->element();

        auto f = form2( _test=Xh, _trial=Xh, _backend=M_backend, _matrix=M );
        f = integrate( markedelements(mesh,"Omega1"), mu[0]*id(u)*idt(v) );
        f += integrate( markedelements(mesh,"Omega2"), mu[1]*grad(u)*trans(gradt(v)) );
        f += integrate( markedfaces(mesh,"Gamma2"), (mu[1]*mu[0] + mu[0]*mu[0])*grad(u)*N()*idt(v) );
        f += integrate( markedfaces(mesh,"Gamma1"), (mu[1]*mu[1])*id(u)*idt(v) );

        return M;
    }

    // These 3 functions are only needed for compilation
    sparse_matrix_ptrtype assembleForMDEIMnl( parameter_type mu, element_type const& u)
    {
        return M;
    }
    space_ptrtype functionSpace()
    {
        return Xh;
    }
    element_type solve( parameter_type const& mu )
    {
        return Xh->element();
    }
    std::string modelName()
    {
        return "test_mdeim";
    }



private :
    space_ptrtype Xh;
    backend_ptrtype M_backend;
    parameterspace_ptrtype Dmu;
    sparse_matrix_ptrtype M;

};

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( mdeim_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    boost::shared_ptr<MDeimTest> m( new MDeimTest );
    m->run();
}

BOOST_AUTO_TEST_SUITE_END()
