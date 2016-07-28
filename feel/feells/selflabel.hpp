#ifndef SELFLABEL_HPP
#define SELFLABEL_HPP 1

#include <boost/serialization/set.hpp>

namespace Feel
{

template< typename space_type, typename space_P0_type >
class SelfLabel
{

public :

    typedef SelfLabel< space_type, space_P0_type > self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef boost::shared_ptr< space_type > space_ptrtype;
    typedef boost::shared_ptr<space_P0_type> space_P0_ptrtype;

    typedef typename space_type::element_ptrtype element_ptrtype;
    typedef typename space_P0_type::element_ptrtype elementP0_ptrtype;

    typedef typename space_type::mesh_ptrtype mesh_ptrtype;

    const int markerfluid = 0;
    const int eltIsInPhiNeg = 1;

    SelfLabel( space_ptrtype space, space_P0_ptrtype space_P0 )
        :
        Xh( space ),
        Xh0(space_P0)
        {
        }

    static self_ptrtype New( space_ptrtype space, space_P0_ptrtype space_P0 )
        {
            self_ptrtype sl( new self_type(space, space_P0) );
            return sl;
        }

    void setLabel( element_ptrtype l )
        { this->label = l; }

    void updateLabel( element_ptrtype newphi );

    element_ptrtype getLabel()
        { return label; }

    elementP0_ptrtype getP0Label()
        { return labelP0; }

private :

    space_ptrtype Xh;
    space_P0_ptrtype Xh0;
    space_P0_ptrtype XhP0Submesh;

    mesh_ptrtype submesh;

    element_ptrtype label;
    elementP0_ptrtype labelP0;

    element_ptrtype phi;

    void genereatesubmesh();
    void propagateLabel( int labelValue, elementP0_ptrtype labelOnSubMesh );

    void pushLabelOnMesh( elementP0_ptrtype labelOnSubMesh );

    elementP0_ptrtype makeMarkerElementsWithLabelOnSubmesh();

    void clean()
        {
            phi.reset();
            XhP0Submesh.reset();
            submesh.reset();
        }

}; // SelfLabel


template< typename space_type, typename space_P0_type >
void
SelfLabel<space_type, space_P0_type>::genereatesubmesh()
{
    // return sub-mesh associated to elements having at least one dof not in fluid
    auto mark = Xh0->elementPtr();

    for( auto const &elt : elements( Xh->mesh() ) )
    {
        for( int j=0; j<space_type::fe_type::nDof; ++j)
        {
            if ( phi->localToGlobal(elt.id(), j, 0)<0 )
            {
                mark->assign(elt.id(), 0, 0, eltIsInPhiNeg);
                break;
            }
        }
    }

    Xh->mesh()->updateMarker2( *mark );

    submesh = createSubmesh( Xh->mesh(), marked2elements(Xh->mesh(), eltIsInPhiNeg) );
} // getsubmesh


template< typename space_type, typename space_P0_type >
typename SelfLabel<space_type, space_P0_type>::elementP0_ptrtype
SelfLabel<space_type, space_P0_type>::makeMarkerElementsWithLabelOnSubmesh()
{
    // return marker on the submesh, where elements having at least one dof with a label by this label

    auto mark = XhP0Submesh->elementPtr();
    *mark = vf::project( _space=XhP0Submesh, _range=elements(submesh), _expr=cst( markerfluid ) );

    for (auto const& elt_submesh : elements(submesh) )
    {
        const size_type idOnMesh = submesh->subMeshToMesh( elt_submesh.id() );

        for (int j=0; j<space_type::fe_type::nDof; ++j)
        {

            if ( std::abs( label->localToGlobal( idOnMesh, j, 0 ) - markerfluid ) > 1e-6 )
            {
                // label is not fluid
                mark->assign( elt_submesh.id(), 0, 0, label->localToGlobal( idOnMesh, j, 0 ) );
                break;
            }
        }
    }

    return mark;
}


template< typename space_type, typename space_P0_type >
void
SelfLabel<space_type, space_P0_type>::updateLabel(element_ptrtype newphi)
{
    this->phi = newphi;

    this->genereatesubmesh();
    XhP0Submesh = space_P0_type::New(submesh);

    auto labelOnSubMesh = this->makeMarkerElementsWithLabelOnSubmesh();
    submesh->updateMarker2( *labelOnSubMesh );

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



template< typename space_type, typename space_P0_type >
void
SelfLabel<space_type, space_P0_type>::propagateLabel( int labelValue, elementP0_ptrtype labelOnSubMesh )
{
    std::unordered_set< size_type > eltsToVisit;
    // communicate to neighbors which ghost has to be visited
    int neighborSubdomains = submesh->neighborSubdomains().size();
    int nbRequest = 2*neighborSubdomains;
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;

    for ( auto const& elt : marked2elements(submesh, labelValue) )
    {
        eltsToVisit.insert( elt.id() );
    }

    // auto exp = exporter(_mesh=submesh, _name="propagation");
    // int iter=0;

    bool eltsToVisitAreEmptyOnAllProc = false;

    while( !eltsToVisitAreEmptyOnAllProc )
    {
        // contains the dof ghost which have been marked
        std::map<rank_type, std::set<size_type> > dataToRecv;
        std::map<rank_type, std::set<size_type> > dataToSend;

        while( ! eltsToVisit.empty() )
        {
            auto elt_id = eltsToVisit.begin();
            auto const & elt = submesh->element( *elt_id );

            // exp->step(iter++)->add("labelOnSubMesh", *labelOnSubMesh);
            // exp->save();

            for (uint16_type face_id = 0; face_id<elt.nTopologicalFaces(); ++face_id)
            {
                auto const & elt_neigh = elt.neighbor( face_id );
                size_type elt_neigh_id = elt_neigh.first;

                if ( elt_neigh_id == invalid_size_type_value )
                    continue;

                // pid of the subdomain owning the element
                const rank_type pid = elt_neigh.second;

                // element is a ghost, need to tell it's owner to visit it
                if (pid != submesh->worldComm().localRank() )
                {
                    size_type idOnOwner = submesh->element( elt_neigh_id, pid ).idInOthersPartitions(pid) ;
                    dataToSend[pid].insert( idOnOwner );
                }
                else
                {
                    // check if the neighbor elt is already labelized
                    if ( std::abs(int(labelOnSubMesh->localToGlobal(elt_neigh_id, 0, 0))-labelValue) > 1e-6 )
                        eltsToVisit.insert(elt_neigh_id);
                }
            }

            labelOnSubMesh->assign(*elt_id, 0,0, labelValue);
            eltsToVisit.erase( elt_id );
        }

        // communicate to neighbors which ghost has to be visited
        cptRequest=0;
        for ( rank_type neighborRank : submesh->neighborSubdomains() )
        {
            reqs[cptRequest++] = submesh->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank] );
            reqs[cptRequest++] = submesh->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank] );
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
        eltsToVisitAreEmptyOnAllProc = mpi::all_reduce(submesh->worldComm(),
                                                       eltstovisit_isempty,
                                                       std::logical_and<bool>() );
    }

    delete [] reqs;

}// propagateLabel



template< typename space_type, typename space_P0_type >
void
SelfLabel<space_type, space_P0_type>::pushLabelOnMesh( elementP0_ptrtype labelOnSubMesh )
{
    // take a marker P0 of labels on a submesh
    // create a Pn (levelset order) label field (L0)
    // all the dof phi<0 on a elt labeled i have value i
    // all dof phi>0 have value markerfluid

    labelP0 = Xh0->elementPtr();
    labelP0->setConstant( markerfluid );

    for (auto const& elt_submesh : elements(submesh) )
    {
        const size_type idOnMesh = submesh->subMeshToMesh( elt_submesh.id() );

        const int labelElt = boost::math::iround( labelOnSubMesh->localToGlobal( elt_submesh.id(), 0, 0 ) );

        labelP0->assign( idOnMesh, 0, 0, labelElt );

        for (int j=0; j<space_type::fe_type::nDof; ++j)
        {
            label->assign( idOnMesh, j, 0,
                           phi->localToGlobal(idOnMesh, j, 0) < 0 ? labelElt : markerfluid );
        }
    }
} // pushLabelOnMesh


template< typename space_type, typename space_P0_type >
boost::shared_ptr< SelfLabel<space_type, space_P0_type> >
selfLabel( boost::shared_ptr<space_type> space, boost::shared_ptr<space_P0_type> spaceP0 )
{
    auto sl = SelfLabel<space_type, space_P0_type>::New(space, spaceP0);
    return sl;
}


} // namespace Feel




#endif
