#define BOOST_TEST_MODULE deim testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelmor/deim.hpp>
#include <feel/feelmor/parameterspace.hpp>

using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test DEIM Options" );

    options.add( feel_options() )
        .add(deimOptions())
        .add(crbSEROptions());
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_deim" ,
                     "test_deim" ,
                     "0.1",
                     "DEIM test",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2017 Feel++ Consortium" );

    return about;
}

template <int POrder,bool IsVectorial=false,bool IsMatricial=false>
class DeimTest :
    public std::enable_shared_from_this<DeimTest<POrder,IsVectorial,IsMatricial>>
{
public :
    typedef DeimTest<POrder,IsVectorial,IsMatricial> self_type;
    static const bool is_vect = IsVectorial;
    static const bool is_mat = IsMatricial;
    static const bool by_block = false;

    typedef double value_type;
    typedef Mesh<Simplex<2> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Pch_type<mesh_type,POrder> scalarspace_type;
    typedef Pchv_type<mesh_type,POrder> vectorialspace_type;
    typedef typename mpl::if_< mpl::bool_<is_vect>,
                               vectorialspace_type,
                               scalarspace_type>::type space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef Backend<double> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename mpl::if_<mpl::bool_<is_mat>,
                              sparse_matrix_type,
                              vector_type >::type tensor_type;
    typedef std::shared_ptr<tensor_type> tensor_ptrtype;

    typedef ParameterSpace<2> parameterspace_type;
    typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef std::shared_ptr<sampling_type> sampling_ptrtype;
    typedef parameterspace_type::element_type parameter_type;

    typedef DEIMBase<parameterspace_type,space_type,tensor_type> deim_type;
    typedef std::shared_ptr<deim_type> deim_ptrtype;


    DeimTest( std::string const& prefix = "") :
        M_backend( backend() ),
        Dmu( parameterspace_type::New(2) ),
        M_uuid( Environment::randomUUID( true ) ),
        M_prefix(prefix)
    {
        auto mesh = loadMesh( _mesh=new mesh_type, _filename="test_deim.geo");
        Xh = space_type::New( mesh );
        setFunctionSpaces(Xh);
    }

    void setModelOnlineDeim( std::string name )
    {
        M_backend = backend( _name=name, _kind="eigen_dense", _worldcomm=Environment::worldCommSeqPtr() );
    }

    void setFunctionSpaces( space_ptrtype Rh )
    {
        Xh = Rh;
        if ( is_mat )
            M = M_backend->newMatrix(_test=Xh,_trial=Xh);
        else
            V = M_backend->newVector(Xh);
    }

    uuids::uuid uuid() const { return M_uuid; }
    std::string prefix() const { return M_prefix; }
    parameterspace_ptrtype parameterSpace() { return Dmu;}
    void setOnlineModel() {}

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
        Pset->equidistributeProduct( Ne , true , "deim_test_sampling" );

        Environment::setOptionValue("deim.rebuild-database", true );
        initDeim( Pset );

        M_deim->run();
        int m = M_deim->size();
        int real_m = is_mat? 4:3;

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "Number of mode in DEIM : " << m );
        BOOST_CHECK( m==real_m );

        Pset = Dmu->sampling();
        Pset->randomize( 100 , true , "deim_test_sampling" );

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "Compare expansion with monolithic" );
        auto base = M_deim->q();
        for ( auto const& mu : *Pset )
        {
            auto coeff = M_deim->beta(mu);
            betas.push_back( coeff );

            double norm=0;
            if ( is_mat )
            {
                assembleForMDEIM(mu,0);
                for ( int i=0; i<m; i++ )
                    add( -coeff(i), base[i] );
                norm = M->linftyNorm();
            }
            else
            {
                assembleForDEIM(mu,0);
                for ( int i=0; i<m; i++ )
                    add( -coeff(i), base[i] );
                norm = V->linftyNorm();
            }

            BOOST_CHECK_SMALL( norm, 1e-9 );
        }


        // We rebuild a new DEIM object with same uuid so he will reload the db
        Environment::setOptionValue("deim.rebuild-database", false );
        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "Rebuild and check" );
        initDeim( Pset );

        // Here we check if the betas are well computed without rebuilding the basis tensors
        int i = 0;
        for ( auto const& mu : *Pset )
        {
            auto coeff = M_deim->beta(mu);
            for ( int k=0; k<coeff.size(); k++)
                BOOST_CHECK_CLOSE( coeff(k), betas[i](k) , 1e-9 );
            i++;
        }

        // Now we rebuild the basis tensors from the data in the db
        // and we check if the new basis tensors are the equal to the previous.
        // We also recheck the coeff after the rebuild
        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "Reassemble and check" );
        M_deim->run();
        auto base2 = M_deim->q();
        for ( int m=0; m<base.size(); m++ )
        {
            add(  -1, base2[m], base[m] );
            double norm = base[m]->linftyNorm();
            BOOST_CHECK_SMALL( norm, 1e-9 );
        }
        i = 0;
        for ( auto const& mu : *Pset )
        {
            auto coeff = M_deim->beta(mu);
            for ( int k=0; k<coeff.size(); k++)
                BOOST_CHECK_CLOSE( coeff(k), betas[i](k) , 1e-9 );
            i++;
        }
    }


    vector_ptrtype assembleForDEIM( parameter_type const& mu, int const& tag )
    {
        return assembleForDEIM( mu, mpl::bool_<is_vect>() );
    }
    vector_ptrtype assembleForDEIM( parameter_type const& mu, mpl::bool_<false> )
    {
        auto mesh = Xh->mesh();
        auto u = Xh->element();

        auto f = form1( _test=Xh, _trial=Xh, _backend=M_backend, _vector=V );
        f = integrate( _range=markedelements(mesh,"Omega1"), _expr=math::cos(mu[0])*id(u) );
        f += integrate( _range=markedelements(mesh,"Omega2"), _expr=math::sin(mu[1])*id(u) );
        f += integrate( _range=markedfaces(mesh,"Gamma1"), _expr=math::exp(mu[0])*mu[1]*id(u) );

        return V;
    }
    vector_ptrtype assembleForDEIM( parameter_type const& mu, mpl::bool_<true> )
    {
        auto mesh = Xh->mesh();
        auto u = Xh->element();

        auto f = form1( _test=Xh, _trial=Xh, _backend=M_backend, _vector=V );
        f = integrate( _range=markedelements(mesh,"Omega1"), _expr=math::cos(mu[0])*trans(id(u))*oneX() );
        f += integrate( _range=markedelements(mesh,"Omega2"), _expr=math::sin(mu[1])*trans(id(u))*oneY() );
        f += integrate( _range=markedfaces(mesh,"Gamma1"),_expr=math::exp( mu[0])*mu[1]*trans(id(u))*N() );

        return V;
    }

    sparse_matrix_ptrtype assembleForMDEIM( parameter_type const& mu, int const& tag )
    {
        auto mesh = Xh->mesh();
        auto u = Xh->element();
        auto v = Xh->element();

        auto f = form2( _test=Xh, _trial=Xh, _backend=M_backend, _matrix=M );
        f = integrate( _range=markedelements(mesh,"Omega1"), _expr=mu[0]*inner(id(u),idt(v)) );
        f += integrate( _range=markedelements(mesh,"Omega2"), _expr=mu[1]*inner(grad(u),gradt(v)) );
        f += integrate( _range=markedfaces(mesh,"Gamma2"), _expr=(mu[1]*mu[0] + mu[0]*mu[0])*inner(grad(u)*N(),idt(v)) );
        f += integrate( _range=markedfaces(mesh,"Gamma1"), _expr=(mu[1]*mu[1])*inner(id(u),idt(v)) );

        return M;
    }

    // These functions are only needed for compilation
    vector_ptrtype assembleForDEIMnl( parameter_type const& mu, element_type const& u, int const& tag )
    {
        return V;
    }
    sparse_matrix_ptrtype assembleForMDEIMnl( parameter_type const& mu, element_type const& u, int const& tag )
    {
        return M;
    }
    space_ptrtype functionSpace()
    {
        return Xh;
    }
    virtual std::pair<element_type,bool> safeSolve( parameter_type const& mu )
    {
        auto u = Xh->element();
        return std::make_pair( u, true );
    }
    std::string modelName()
    {
        return "test_deim";
    }
    void initOnlineModel()
    {}
    typename space_type::mesh_support_vector_type
        functionspaceMeshSupport( mesh_ptrtype const& mesh ) const
    {
        return typename space_type::mesh_support_vector_type();
    }
 private :
    void initDeim( sampling_ptrtype Pset )
    {
        return initDeim( Pset, mpl::bool_<is_mat>() );
    }
    void initDeim( sampling_ptrtype Pset, mpl::bool_<true> )
    {
        M_deim = mdeim( _model=this->shared_from_this(),
                        _sampling=Pset );
    }
    void initDeim( sampling_ptrtype Pset, mpl::bool_<false> )
    {
        M_deim = deim( _model=this->shared_from_this(),
                       _sampling=Pset );
    }

    //! evaluate V= V+a*vec
    void add( double const& a, vector_ptrtype vec, vector_ptrtype v=nullptr )
    {
        if ( !v )
            V->add( a, vec );
        else
            v->add( a, vec );
    }
    //! evaluate M= M+a*mat
    void add(  double const& a, sparse_matrix_ptrtype mat, sparse_matrix_ptrtype m=nullptr )
    {
        if ( !m )
            M->addMatrix( a, mat );
        else
            m->addMatrix( a, mat );
    }
private :
    space_ptrtype Xh;
    backend_ptrtype M_backend;
    parameterspace_ptrtype Dmu;
    vector_ptrtype V;
    sparse_matrix_ptrtype M;

    uuids::uuid M_uuid;
    deim_ptrtype M_deim;
    std::vector<vectorN_type> betas;

    std::string M_prefix;
};

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( deim_suite )

BOOST_AUTO_TEST_CASE( test_3 )
{
    std::shared_ptr<DeimTest<3>> m( new DeimTest<3> );
    m->run();
}

BOOST_AUTO_TEST_CASE( test_2v )
{
    std::shared_ptr<DeimTest<2,true>> m( new DeimTest<2,true> );
    m->run();
}

BOOST_AUTO_TEST_CASE( test_m1 )
{
    std::shared_ptr<DeimTest<1,false,true>> m( new DeimTest<1,false,true> );
    m->run();
}

BOOST_AUTO_TEST_CASE( test_m1v )
{
    std::shared_ptr<DeimTest<1,true,true>> m( new DeimTest<1,true,true> );
    m->run();
}

BOOST_AUTO_TEST_SUITE_END()
