#ifndef SELFLABEL_HPP
#define SELFLABEL_HPP 1

namespace Feel
{

template< typename space_type, typename space_P0_type >
class SelfLabel
{

public :

    typedef boost::shared_ptr< space_type > space_ptrtype;
    typedef boost::shared_ptr<space_P0_type> space_P0_ptrtype;

    const int markerfluid = 0;
    const int eltIsInPhiNeg = 1;

    SelfLabel( space_ptrtype space, space_P0_ptrtype space_P0 )
        :
        Xh( space ),
        Xh0(space_P0)
        {
        }

    void updateLabel( space_ptrtype newphi );

    space_type::element_ptrtype getLabel()
        { return label; }

    space_P0_type::element_ptrtype getP0Label()
        { return labelP0; }

    void markElementsWithLabel(space_type::mesh_ptrtype msh );

    void clean()
        {
            // submesh, subspace Xh0, pointeur sur phi ?
        }

private :

    space_ptrtype Xh;
    space_P0_ptrtype Xh0;

    space_type::mesh_ptrtype submesh;

    space_type::element_ptrtype label;
    space_P0_type::element_ptrtype labelP0;

    space_type::element_ptrtype phi;

    space_type::mesh_ptrtype getsubmesh();

}; // SelfLabel


template< typename space_type >
SelfLabel<space_type>::space_type::mesh_ptrtype
SelfLabel<space_type>::getsubmesh()
{
    // return sub-mesh associated to elements having at least one dof not in fluid
    auto mark = Xh0->elementPtr();

    for( auto const &elt : elements( Xh->mesh() ) )
    {
        for( int j=0; j<space_type::fe_type::nDof; ++j)
        {
            if ( phi->localToGlobal(elt.id(), j, 0)<0 )
            {
                mark->assign(elt.id(), j, 0, eltIsInPhiNeg);
                break;
            }
        }
    }

    Xh->mesh()->updateMarker2( *mark );

    auto submesh = createSubmesh( Xh->mesh(), marked2elements(Xh->mesh(), eltIsInPhiNeg) );
    return submesh;
} // getsubmesh


template< typename space_type >
SelfLabel<space_type>::elementP0_ptrtype
SelfLabel<space_type>::makeMarkerElementsWithLabel( space_type::mesh_ptrtype msh )
{
    // return marker on the elements having at least one dof with a label by this label

    auto mark = XhP0Submesh->elementPtr();
    *mark = vf::project( _space=XhP0Submesh, _range=elements(msh), _expr=cst( markerfluid ) );

    for (auto const& elt : elements(msh) )
    {
        for (int j=0; j<space_type::fe_type::nDof; ++j)
        {
            if ( std::abs( label->localToGlobal( elt.id(), j, 0 ) - markerfluid ) > 1e-6 )
            {
                // label is not fluid
                mark.assign( elt.id(), j, 0, label->localToGlobal( elt.id(), j, 0 ) );
                break;
            }
        }
    }

    return mark;
}


template< typename space_type >
void
SelfLabel<space_type>::updateLabel(space_ptrtype newphi)
{
    this->phi = newphi;

    auto submesh = this->getsubmesh();
    // make subspace P0 XhP0Submesh

    auto labelOnSubMesh = this->makeMarkerElementsWithLabel( submesh );
    submesh->updateMarker2( *labelOnSubMesh );

    const int nbLabels = boost::iround(labelOnSubMesh->max());

    // to do in multi-thread ?
    for (int i=1; i<=nbLabels; ++i)
        this->propagateLabel( i, labelOnSubMesh );

    // checker que tous les elements du sous maillage sont marqués
    // enable only in debug maybe
    this->checkEltsAreLabeled();

    // label is updated here
    this->pushLabelOnMesh( labelOnSubMesh );

} // updateLabel



template< typename space_type >
void
SelfLabel<space_type>::propagateLabel( int labelValue, elementP0_ptrtype labelOnSubMesh )
{
    std::unordered_set< size_type > eltsToVisit;

    for ( auto const& elt : marked2elements(mesh, labelValue) )
    {
        eltsToVisit.insert( elt.id() );
    }

    while( ! eltsToVisit.empty() )
    {
        auto elt = eltsToVisit.begin();

        for (uint16_type face_id = 0; face_id<elt.nTopologicalFaces(); ++face_id)
        {
            auto const & elt_neigh = elt.neighbor( face_id );
            size_type elt_neigh_id = elt_neigh.first;

            // pid auquel appartient l'element
            rank_type pid = elt_neigh.second;

            // pid courant : mesh->worldComm().localRank();

            // check if the neighbor elt is already labelized
            if ( std::abs(labelOnSubMesh->localToGlobal(elt.id(), 0, 0)-labelValue) > 1e-6 )
                eltsToVisit.insert(elt_neigh_id);
        }

        labelOnSubMesh->assign(elt.id(), 0,0, labelValue);
        eltsToVisit.erase( elt );
    }

}// propagateLabel



template< typename space_type >
void
SelfLabel<space_type>::checkEltsAreLabeled()
{
    for (auto const & elt : elements(submesh))
    {
        CHECK( vf::abs(elt.marker2().value()-markerfluid)>1e-6 )
            <<"element " << elt.id() << " has still a marker 'fluid'. It is supposed to have the value of one of the domain label\n";
    }
}// checkEltsAreLabeled


template< typename space_type >
void
SelfLabel<space_type>::pushLabelOnMesh( elementP0_ptrtype labelOnSubMesh )
{
    // take a marker P0 of labels on a submesh
    // create a Pn (levelset order) label field (L0)
    // all the dof phi<0 on a elt labeled i have value i
    // all dof phi>0 have value markerfluid

    labelP0 = Xh0->elementPtr();
    labelP0->setConstant( markerfluid );
    labelP0->on(_range=marked2element(XhP0->mesh(), eltIsInPhiNeg), _expr=idv(labelOnSubMesh) );

    for (auto const& elt: elements(Xh->mesh()))
    {
        const int labelElt = boost::iround(labelP0->localToGlobal(elt.id(), 0, 0));
        for (int j=0; j<space_type::fe_type::nDof; ++j)
        {
            label->assign( elt.id(), j, 0,
                           phi->localToGlobal(elt.id(), j, 0) < 0 ? labelElt : markerfluid );
        }
    }

} // pushLabelOnMesh



} // namespace Feel




#endif
