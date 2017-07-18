#define BOOST_TEST_MODULE test_reducedglobaldof

#include <testsuite.hpp>

#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;

#define POUT std::cout<<"["<<Environment::worldComm().globalRank()<<"] "

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_reducedglobaldof )

BOOST_AUTO_TEST_CASE( globaldof_recovery )
{
    auto mesh = unitCircle();
    auto Xh = Pch<1>( mesh );
    auto u = Xh->element();

    auto expr = 20*abs(sin(2*Pi*Px())*exp(Px()*Py()+Py()*Py()));

    auto f = form1( Xh );
    f = integrate( elements(mesh), expr*id(u) );
    auto V = f.vectorPtr();

    // get the maximum entry of the vector and corresponding index (dof)
    auto absV = V->clone();
    *absV = *V;
    absV->abs();
    int index = 0;
    double max = absV->maxWithIndex( &index );
    int proc_n = V->map().procOnGlobalCluster(index);

    // get all the list of all the elements containing the dof
    std::set<int> element_ids;
    std::vector<std::pair<int,int>> local_id_dof;
    if ( Environment::worldComm().globalRank()==proc_n )
        for ( auto const& dof : Xh->dof()->globalDof(index-V->map().firstDofGlobalCluster()) )
        {
            element_ids.insert( dof.second.elementId() );
            local_id_dof.push_back( std::make_pair(dof.second.elementId(),dof.second.localDof()) );
        }

    // create a submesh and subspace on these elements
    auto newmesh = createSubmesh( mesh, idelements(mesh,element_ids.begin(), element_ids.end()) );
    auto Rh = Pch<1>(newmesh);
    auto v = Rh->element();

    auto fR = form1(Rh);
    fR = integrate( elements(newmesh), expr*id(v) );
    auto VR = fR.vectorPtr();

    // check if the indexR in the reduced space correspond to the good entry
    double valueR=0;
    for ( auto const& l : local_id_dof )
    {
        auto sub_id = newmesh->meshToSubMesh( l.first );

        auto indexR = Rh->dof()->localToGlobalId( sub_id, l.second );
        proc_n = VR->map().procOnGlobalCluster(indexR);
        if ( Environment::worldComm().globalRank()==proc_n )
            valueR = math::abs(VR->operator()( indexR - VR->map().firstDofGlobalCluster() ));

        BOOST_CHECK_CLOSE( max, valueR, 1e-9 );
    }

    boost::mpi::broadcast( Environment::worldComm(), valueR, proc_n );
    BOOST_CHECK_CLOSE( max, valueR, 1e-9 );
}

BOOST_AUTO_TEST_CASE( parallel_to_seq )
{
    typedef Mesh<Simplex<2,1>> mesh_type;
    typedef Pch_type<mesh_type,1> space_type;

    auto expr = 20*abs(sin(2*Pi*Px())*exp(Px()*Py()+Py()*Py()));
    auto mesh = unitCircle();
    auto Xh = space_type::New( _mesh=mesh );
    auto u = Xh->element();

    auto f = form1( Xh );
    f = integrate( elements(mesh), expr*id(u) );
    auto V = f.vectorPtr();

    // get the maximum entry of the vector and corresponding index (dof)
    auto absV = V->clone();
    *absV = *V;
    absV->abs();
    int index = 0;
    double max = absV->maxWithIndex( &index );

    int proc_n = V->map().procOnGlobalCluster(index);

    std::set<int> element_ids;
    std::set<int> points_id;
    int elt_id=-1;
    int l_dof=-1;

    if ( Environment::worldComm().globalRank()==proc_n )
    {
        for ( auto const& dof : Xh->dof()->globalDof(index-V->map().firstDofGlobalCluster()) )
        {
            element_ids.insert( dof.second.elementId() );
            if ( elt_id==-1)
            {
                elt_id=dof.second.elementId();
                l_dof=dof.second.localDof();
            }
        }

        auto elt = mesh->element(elt_id);
        for ( int p=0; p<elt.nPoints(); p++ )
            points_id.insert( elt.point(p).id()+1 );
    }
    boost::mpi::broadcast( Environment::worldComm(), points_id, proc_n );
    boost::mpi::broadcast( Environment::worldComm(), l_dof, proc_n );

    // build submesh and reload it in sequential mode
    auto newmesh = createSubmesh( mesh, idelements(mesh,element_ids.begin(), element_ids.end()) );
    saveGMSHMesh(_mesh=newmesh, _filename="submesh.msh" );
    auto seqmesh = loadMesh( _mesh=new mesh_type,
                             _filename="submesh.msh", _worldcomm= Environment::worldCommSeq() );


    // build a map to find elt id from points ids
    std::map<std::set<int>,int> elts_map;
    for ( auto const& eltWrap : elements(seqmesh) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        std::set<int> pts_id;
        for ( int p=0; p<elt.nPoints(); p++ )
            pts_id.insert( elt.point(p).id() );
        elts_map[pts_id] = elt.id();
    }

    // find the corresponding element in the sequential submesh
    auto map_it = elts_map.find( points_id );
    int elt_idR=-1;
    if ( map_it!=elts_map.end() )
        elt_idR = map_it->second;
    else
        POUT<<"elt id not found in map\n";

    auto Rh = space_type::New(seqmesh,
                              _worldscomm=std::vector<WorldComm>(1,Environment::worldCommSeq()) );
    auto v = Rh->element();
    auto fR = form1(Rh);
    fR = integrate( elements(seqmesh), expr*id(v) );
    auto VR = fR.vectorPtr();
    auto indexR = Rh->dof()->localToGlobalId( elt_idR, l_dof );

    double valueR = VR->operator()( indexR );
    BOOST_CHECK_CLOSE( max, valueR, 1e-9 );
}



BOOST_AUTO_TEST_SUITE_END()
