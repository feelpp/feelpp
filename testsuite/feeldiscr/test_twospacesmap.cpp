#define BOOST_TEST_MODULE test_twospaces_map

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/thch.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/twospacesmap.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_twospacesmap )

BOOST_AUTO_TEST_CASE( onespace )
{
    typedef Mesh<Simplex<2,1>> mesh_type;
    typedef Pchv_type<mesh_type,2> space_type;

    auto mesh = unitCircle();
    auto Xh = Pchv<2>( mesh );
    auto u = Xh->element();

    u.on( _range=elements(mesh), _expr=20*abs(sin(2*Pi*Px())*exp(Px()*Py()+Py()*Py()))*oneX()
          + 20*abs(cos(2*Pi*Px())*exp(Px()*Py()+Py()*Py()))*oneY() );

    int size = u.size();
    std::vector<int> randoms;
    std::set<int> element_ids;

    std::srand (std::time(NULL));
    if ( Environment::isMasterRank() )
        for ( int i=0; i<10; i++ )
            randoms.push_back( std::rand() % size );
    boost::mpi::broadcast( Environment::worldComm(), randoms, Environment::worldComm().masterRank() );

    for ( auto index : randoms )
    {
        auto searchGpDof = Xh->dof()->searchGlobalProcessDof( index );
        if ( boost::get<0>( searchGpDof ) )
        {
            size_type gpdof = boost::get<1>( searchGpDof );
            for ( auto const& dof : Xh->dof()->globalDof( gpdof ) )
            {
                size_type eltId = dof.second.elementId();
                if ( Xh->mesh()->element( eltId ).isGhostCell() )
                    continue;
                element_ids.insert( eltId );

            }
        }
    }

    auto newmesh = createSubmesh( _mesh=mesh, _range=idelements(mesh,element_ids.begin(), element_ids.end()) );
    saveGMSHMesh(_mesh=newmesh, _filename="submesh.msh" );
    Environment::worldComm().barrier();

    auto seqmesh = loadMesh( _mesh=new mesh_type,
                             _filename="submesh.msh", _worldcomm= Environment::worldCommSeqPtr() );

    auto Rh = space_type::New(_mesh=seqmesh );
    auto my_map = TwoSpacesMap<space_type>( Rh, Xh );

    auto ur = Rh->element();
    ur.on( _range=elements(seqmesh), _expr=20*abs(sin(2*Pi*Px())*exp(Px()*Py()+Py()*Py()))*oneX()
           + 20*abs(cos(2*Pi*Px())*exp(Px()*Py()+Py()*Py())
                    )*oneY() );

    for ( auto index : randoms )
    {
        double value = 0;
        int proc_number = Xh->map().procOnGlobalCluster(index);

        if ( Environment::worldComm().globalRank()==proc_number )
        {
            auto searchGpDof = Xh->dof()->searchGlobalProcessDof( index );
            CHECK( boost::get<0>( searchGpDof ) ) << "GPDof not found\n";
            size_type gpdof = boost::get<1>( searchGpDof );
            value = u( gpdof );
        }
        boost::mpi::broadcast( Environment::worldComm(), value, proc_number );

        int index_r = my_map.clusterToSequential( index );

        double value_r = ur( index_r );
        BOOST_CHECK_CLOSE( value, value_r, 1e-9 );
    }

    my_map.project( ur, u );

    double diff = normL2( _range=elements(newmesh), _expr=idv(u)-idv(ur) );
    BOOST_CHECK_SMALL( diff, 1e-9 );
}


BOOST_AUTO_TEST_CASE( composite )
{
    typedef Mesh<Simplex<2,1>> mesh_type;
    typedef THch_type<1,mesh_type> space_type;

    auto mesh = unitCircle();
    auto Xh = space_type::New( mesh );
    auto Xh0 = Xh->functionSpace<0>();
    auto Xh1 = Xh->functionSpace<1>();

    auto U = Xh->element();
    auto u = U.template element<0>();
    auto p = U.template element<1>();
    auto expr_u = 20*sin(2*Pi*Px())*exp(Px()*Py()+Py()*Py())*oneX() + 20*cos(2*Pi*Px())*exp(Px()*Py()+Py()*Py())*oneY();
    auto expr_p = cos(Pi*Px()*pow(Py(),3))*exp( Px()-Py()*(Px() - Py() ) );
    u.on( _range=elements(mesh), _expr=expr_u );
    p.on( _range=elements(mesh), _expr=expr_p  );

    int size = U.size();
    std::vector<int> randoms;
    std::set<int> element_ids;

    // build a random set of index in the range of the vector
    std::srand (std::time(NULL));
    if ( Environment::isMasterRank() )
        for ( int i=0; i<10; i++ )
            randoms.push_back( std::rand() % size );
    boost::mpi::broadcast( Environment::worldComm(), randoms, Environment::worldComm().masterRank() );

    for ( auto index : randoms )
    {
        auto searchGpDof = Xh->dof()->searchGlobalProcessDof( index );
        if ( boost::get<0>( searchGpDof ) )
        {
            size_type gpdof = boost::get<1>( searchGpDof );
            int space = Xh->dof()->databaseIndexFromContainerId( gpdof );
            gpdof = Xh->dof()->containerIdToDofId( space, gpdof );

            if ( space==0 )
            {
                for ( auto const& dof : Xh0->dof()->globalDof( gpdof ) )
                {
                    size_type eltId = dof.second.elementId();
                    if ( Xh0->mesh()->element( eltId ).isGhostCell() )
                        continue;
                    element_ids.insert( eltId );
                }
            }
            else if ( space==1 )
            {
                for ( auto const& dof : Xh1->dof()->globalDof( gpdof ) )
                {
                    size_type eltId = dof.second.elementId();
                    if ( Xh1->mesh()->element( eltId ).isGhostCell() )
                        continue;
                    element_ids.insert( eltId );
                }
            }
        }
    }

    auto newmesh = createSubmesh( _mesh=mesh, _range=idelements(mesh,element_ids.begin(), element_ids.end()) );
    saveGMSHMesh(_mesh=newmesh, _filename="submesh.msh" );
    Environment::worldComm().barrier();

    auto seqmesh = loadMesh( _mesh=new mesh_type,
                             _filename="submesh.msh", _worldcomm= Environment::worldCommSeqPtr() );

    auto Rh = space_type::New(_mesh=seqmesh );

    auto my_map = TwoSpacesMap<space_type>( Rh, Xh );

    auto Ur = Rh->element();
    auto ur = Ur.template element<0>();
    auto pr = Ur.template element<1>();
    ur.on( _range=elements(seqmesh), _expr=expr_u );
    pr.on( _range=elements(seqmesh), _expr=expr_p );

    for ( auto index : randoms )
    {
        double value = 0;
        int proc_number = Xh->map().procOnGlobalCluster(index);

        if ( Environment::worldComm().globalRank()==proc_number )
        {
            auto searchGpDof = Xh->dof()->searchGlobalProcessDof( index );
            CHECK( boost::get<0>( searchGpDof ) ) << "GPDof not found\n";
            size_type gpdof = boost::get<1>( searchGpDof );
            value = U( gpdof );
        }
        boost::mpi::broadcast( Environment::worldComm(), value, proc_number );

        int index_r = my_map.clusterToSequential( index );
        CHECK( index_r>=0 )<<"index "<< index <<" was not found in p_to_s map\n";
        double value_r = Ur( index_r );

        BOOST_CHECK_CLOSE( value, value_r, 1e-9 );
    }
    my_map.project( Ur, U );

    double diff = normL2( _range=elements(newmesh), _expr=idv(u)-idv(ur) );
    BOOST_CHECK_SMALL( diff, 1e-9 );

    diff = normL2( _range=elements(newmesh), _expr=idv(p)-idv(pr) );
    BOOST_CHECK_SMALL( diff, 1e-9 );

    auto Rh0 = Rh->functionSpace<0>();
    auto Rh1 = Rh->functionSpace<1>();

    auto v = Xh0->element();
    auto q = Xh1->element();
    auto vr = Rh0->element();
    auto qr = Rh1->element();

    v.on( _range=elements(mesh), _expr=expr_u );
    q.on( _range=elements(mesh), _expr=expr_p );

    my_map.template project<0>( vr, v );
    my_map.template project<1>( qr, q );

    diff = normL2( _range=elements(newmesh), _expr=idv(v)-idv(vr) );
    BOOST_CHECK_SMALL( diff, 1e-9 );

    diff = normL2( _range=elements(newmesh), _expr=idv(q)-idv(qr) );
    BOOST_CHECK_SMALL( diff, 1e-9 );

}




BOOST_AUTO_TEST_SUITE_END()
