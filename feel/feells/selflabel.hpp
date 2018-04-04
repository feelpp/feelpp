#ifndef SELFLABEL_HPP
#define SELFLABEL_HPP 1

#include <boost/serialization/set.hpp>

namespace Feel {

template< typename SpaceType, typename SpaceP0Type >
class SelfLabel
{

public :

    typedef SelfLabel<SpaceType, SpaceP0Type > self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef SpaceType space_type;
    typedef boost::shared_ptr< space_type > space_ptrtype;
    typedef SpaceP0Type space_P0_type;
    typedef boost::shared_ptr<space_P0_type> space_P0_ptrtype;

    typedef typename space_type::element_ptrtype element_ptrtype;
    typedef typename space_P0_type::element_ptrtype elementP0_ptrtype;

    typedef typename space_type::mesh_ptrtype mesh_ptrtype;

    const int markerfluid = 0;
    const int eltIsInPhiNeg = 1;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    SelfLabel( space_ptrtype const& space, space_P0_ptrtype const& space_P0 )
        :
            M_Xh( space ),
            M_Xh0( space_P0 ),
            M_label( space->elementPtr() ),
            M_labelP0( space_P0->elementPtr() )
    {}

    static self_ptrtype New( space_ptrtype const& space, space_P0_ptrtype const& space_P0 )
    {
        self_ptrtype sl( new self_type(space, space_P0) );
        return sl;
    }

    //--------------------------------------------------------------------//
    void setLabel( element_ptrtype const& l ) { *M_label = *l; }
    element_ptrtype const& getLabel() const { return M_label; }

    elementP0_ptrtype const& getP0Label() const { return M_labelP0; }

    void updateLabel( element_ptrtype const& newphi );

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private :
    void generateSubmesh();
    void propagateLabel( int labelValue, elementP0_ptrtype labelOnSubMesh );

    void pushLabelOnMesh( elementP0_ptrtype labelOnSubMesh );

    elementP0_ptrtype makeMarkerElementsWithLabelOnSubmesh();

    void clean()
    {
        M_phi.reset();
        M_XhP0Submesh.reset();
        M_submesh.reset();
    }

    //--------------------------------------------------------------------//
    space_ptrtype M_Xh;
    space_P0_ptrtype M_Xh0;
    space_P0_ptrtype M_XhP0Submesh;

    mesh_ptrtype M_submesh;

    element_ptrtype M_label;
    elementP0_ptrtype M_labelP0;

    element_ptrtype M_phi;
};


template<typename SpaceType, typename SpaceP0Type>
void
SelfLabel<SpaceType, SpaceP0Type>::generateSubmesh()
{
    // return sub-mesh associated to elements having at least one dof not in fluid
    auto mark = M_Xh0->elementPtr();

    for( auto const &it_elt : elements( M_Xh->mesh() ) )
    {
        auto const& elt = boost::unwrap_ref( it_elt );
        for( int j=0; j<space_type::fe_type::nDof; ++j )
        {
            if ( M_phi->localToGlobal(elt.id(), j, 0)<0 )
            {
                mark->assign(elt.id(), 0, 0, eltIsInPhiNeg);
                break;
            }
        }
    }

    M_Xh->mesh()->updateMarker2( *mark );

    M_submesh = createSubmesh( M_Xh->mesh(), marked2elements(M_Xh->mesh(), eltIsInPhiNeg) );
} // getsubmesh


template<typename SpaceType, typename SpaceP0Type>
typename SelfLabel<SpaceType, SpaceP0Type>::elementP0_ptrtype
SelfLabel<SpaceType, SpaceP0Type>::makeMarkerElementsWithLabelOnSubmesh()
{
    // return marker on the M_submesh, where elements having at least one dof with a label by this label

    auto mark = M_XhP0Submesh->elementPtr();
    *mark = vf::project( 
            _space=M_XhP0Submesh, 
            _range=elements(M_submesh), 
            _expr=cst( markerfluid ) 
            );

    for (auto const& it_elt_submesh : elements(M_submesh) )
    {
        auto const& elt_submesh = boost::unwrap_ref( it_elt_submesh );
        const size_type idOnMesh = M_submesh->subMeshToMesh( elt_submesh.id() );

        for (int j=0; j<space_type::fe_type::nDof; ++j)
        {

            if ( std::abs( M_label->localToGlobal( idOnMesh, j, 0 ) - markerfluid ) > 1e-6 )
            {
                // label is not fluid
                mark->assign( elt_submesh.id(), 0, 0, M_label->localToGlobal( idOnMesh, j, 0 ) );
                break;
            }
        }
    }

    return mark;
}


template<typename SpaceType, typename SpaceP0Type>
void
SelfLabel<SpaceType, SpaceP0Type>::updateLabel(element_ptrtype const& newphi)
{
    M_phi = newphi;

    this->generateSubmesh();
    M_XhP0Submesh = space_P0_type::New(M_submesh);

    auto labelOnSubMesh = this->makeMarkerElementsWithLabelOnSubmesh();
    M_submesh->updateMarker2( *labelOnSubMesh );

    const int nbLabels = boost::math::iround(labelOnSubMesh->max());

    // to do in multi-thread ?
    for (int i=1; i<=nbLabels; ++i)
        this->propagateLabel( i, labelOnSubMesh );

    // to enable only in debug maybe
    //this->checkEltsAreLabeled();

    // label is updated here
    this->pushLabelOnMesh( labelOnSubMesh );

    this->clean();
} // updateLabel



template<typename SpaceType, typename SpaceP0Type>
void
SelfLabel<SpaceType, SpaceP0Type>::propagateLabel( int labelValue, elementP0_ptrtype labelOnSubMesh )
{
    std::unordered_set< size_type > eltsToVisit;
    // communicate to neighbors which ghost has to be visited
    int neighborSubdomains = M_submesh->neighborSubdomains().size();
    int nbRequest = 2*neighborSubdomains;
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;

    for ( auto const& it_elt : marked2elements(M_submesh, labelValue) )
    {
        auto const& elt = boost::unwrap_ref( it_elt );
        eltsToVisit.insert( elt.id() );
    }

     //auto exp = exporter(_mesh=M_submesh, _name="propagation");
     //int iter=0;
     //exp->step(iter++)->add("labelOnSubMesh", *labelOnSubMesh);
     //exp->save();

    bool eltsToVisitAreEmptyOnAllProc = false;

    while( !eltsToVisitAreEmptyOnAllProc )
    {
        // contains the dof ghost which have been marked
        std::map<rank_type, std::set<size_type> > dataToRecv;
        std::map<rank_type, std::set<size_type> > dataToSend;

        while( ! eltsToVisit.empty() )
        {
            auto elt_id = eltsToVisit.begin();
            auto const & elt = M_submesh->element( *elt_id );

             //exp->step(iter++)->add("labelOnSubMesh", *labelOnSubMesh);
             //exp->save();

            for (uint16_type face_id = 0; face_id<elt.nTopologicalFaces(); ++face_id)
            {
                size_type elt_neigh_id = elt.neighbor( face_id );
                if ( elt_neigh_id == invalid_size_type_value )
                    continue;
                auto const & elt_neigh = M_submesh->element( elt_neigh_id );

                // pid of the subdomain owning the element
                const rank_type pid = elt_neigh.processId();

                // element is a ghost, need to tell it's owner to visit it
                if (pid != M_submesh->worldComm().localRank() )
                {
                    size_type idOnOwner = M_submesh->element( elt_neigh_id ).idInOthersPartitions(pid) ;
                    dataToSend[pid].insert( idOnOwner );
                }
                else
                {
                    // check if the neighbor elt is already labelled
                    if ( std::abs(int(labelOnSubMesh->localToGlobal(elt_neigh_id, 0, 0))-labelValue) > 1e-6 )
                        eltsToVisit.insert(elt_neigh_id);
                }
            }

            labelOnSubMesh->assign(*elt_id, 0,0, labelValue);
            eltsToVisit.erase( elt_id );
        }

        // communicate to neighbors which ghost has to be visited
        cptRequest=0;
        for ( rank_type neighborRank : M_submesh->neighborSubdomains() )
        {
            reqs[cptRequest++] = M_submesh->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank] );
            reqs[cptRequest++] = M_submesh->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank] );
        }

        mpi::wait_all(reqs, reqs + nbRequest);

        // add the elements added by others
        for (auto const& id_set_pair : dataToRecv )
        {
            for (auto const& elt_candidate_id : id_set_pair.second)
            {
                if ( std::abs(int(labelOnSubMesh->localToGlobal(elt_candidate_id, 0, 0))-labelValue) > 1e-6 )
                    eltsToVisit.insert(elt_candidate_id);
            }
        }

        bool eltstovisit_isempty = eltsToVisit.empty();
        // loop until all the proc have a null value of neighbors to add
        eltsToVisitAreEmptyOnAllProc = mpi::all_reduce(M_submesh->worldComm(),
                eltstovisit_isempty,
                std::logical_and<bool>() );
    }

     //exp->step(iter++)->add("labelOnSubMesh", *labelOnSubMesh);
     //exp->save();

    delete [] reqs;

}// propagateLabel



template<typename SpaceType, typename SpaceP0Type>
void
SelfLabel<SpaceType, SpaceP0Type>::pushLabelOnMesh( elementP0_ptrtype labelOnSubMesh )
{
    // take a marker P0 of labels on a M_submesh
    // create a Pn (levelset order) label field (L0)
    // all the dof phi<0 on a elt labeled i have value i
    // all dof phi>0 have value markerfluid

    M_labelP0->setConstant( markerfluid );
    M_label->setConstant( markerfluid );

    for (auto const& it_elt_submesh : elements(M_submesh) )
    {
        auto const& elt_submesh = boost::unwrap_ref( it_elt_submesh );
        const size_type idOnMesh = M_submesh->subMeshToMesh( elt_submesh.id() );

        const int labelElt = boost::math::iround( labelOnSubMesh->localToGlobal( elt_submesh.id(), 0, 0 ) );

        M_labelP0->assign( idOnMesh, 0, 0, labelElt );

        for (int j=0; j<space_type::fe_type::nDof; ++j)
        {
            M_label->assign( idOnMesh, j, 0,
                    M_phi->localToGlobal(idOnMesh, j, 0) < 0 ? labelElt : markerfluid );
        }
    }
} // pushLabelOnMesh


template<typename SpaceType, typename SpaceP0Type>
boost::shared_ptr< SelfLabel<SpaceType, SpaceP0Type> > 
selfLabel( 
        boost::shared_ptr<SpaceType> const& space, 
        boost::shared_ptr<SpaceP0Type> const& spaceP0 
        )
{
    auto sl = SelfLabel<SpaceType, SpaceP0Type>::New(space, spaceP0);
    return sl;
}


} // namespace Feel




#endif
