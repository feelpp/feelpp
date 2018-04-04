#ifndef _DISTLABELS_FMS_HPP
#define _DISTLABELS_FMS_HPP 1

//#define FM_EXPORT 1

#include <feel/feells/fastmarchingbase.hpp>
#include <feel/feells/selflabel.hpp>

namespace Feel {

template<typename FunctionSpaceType, typename PeriodicityType =NoPeriodicity>
class LabelDistanceFMS : public FastMarchingBase<FunctionSpaceType, PeriodicityType, double>
{
public:
    typedef FastMarchingBase<FunctionSpaceType, PeriodicityType, double> super_type;
    typedef LabelDistanceFMS<FunctionSpaceType, PeriodicityType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef typename super_type::functionspace_type functionspace_type;
    typedef typename super_type::functionspace_ptrtype functionspace_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::element_ptrtype element_ptrtype;
    typedef typename super_type::value_type value_type;
    static const uint16_type Dim = super_type::Dim;

    typedef typename super_type::periodicity_type periodicity_type;

    //--------------------------------------------------------------------//
    LabelDistanceFMS( 
            functionspace_ptrtype const& space,
            periodicity_type periodicity =NoPeriodicity() );

    static self_ptrtype New(
            functionspace_ptrtype const& space,
            periodicity_type periodicity =NoPeriodicity() );

    //--------------------------------------------------------------------//
    void setSelfLabel( element_ptrtype const& label ) { M_selfLabel = label; }
    //--------------------------------------------------------------------//
    element_ptrtype const& getNearestNeighbourLabel() const { return M_NNLabel; }
    element_ptrtype const& getNearestNeighbourDistance() const { return M_NNDistance; }
    element_ptrtype const& getNextNearestNeighbourLabel() const { return M_nextNNLabel; }
    element_ptrtype const& getNextNearestNeighbourDistance() const { return M_nextNNDistance; }

private:
    void initMarch( element_type const& phi, bool useMarker2AsMarkerDone );
    void processDof( size_type idOnProc, value_type val, value_type const& opt_data );
    void updateHeap( size_type idDone );

    element_ptrtype getDistance() const { return super_type::getDistance(); }

    //--------------------------------------------------------------------//
    value_type fmsNNDistRec( 
            std::vector<size_type> & ids,
            size_type idClose,
            value_type phiOld ) const;
    value_type fmsNextNNDistRec( 
            std::vector<size_type> & ids,
            size_type idClose,
            value_type phiOld ) const;
    //--------------------------------------------------------------------//
    element_ptrtype M_selfLabel;
    element_ptrtype M_label;
    element_ptrtype M_labelDist;
    element_ptrtype M_NNLabel;
    element_ptrtype M_NNDistance;
    element_ptrtype M_nextNNLabel;
    element_ptrtype M_nextNNDistance;

    std::map<rank_type, std::vector<boost::tuple<size_type, value_type, value_type>>> M_dataToSend, M_dataToRecv;

#if defined( FM_EXPORT )
    boost::shared_ptr<Exporter<typename super_type::mesh_type>> M_ex;
    int M_count_iteration;
#endif
};

#define LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS \
    template<typename FunctionSpaceType, typename PeriodicityType> \
    /**/ 
#define LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE \
    LabelDistanceFMS<FunctionSpaceType, PeriodicityType> \
    /**/ 

LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS
LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::LabelDistanceFMS(
        functionspace_ptrtype const& space,
        periodicity_type periodicity )
:
    super_type( space, periodicity ),
    M_label( space->elementPtr() ),
    M_labelDist( space->elementPtr() ),
    M_NNLabel( space->elementPtr() ),
    M_NNDistance( space->elementPtr() ),
    M_nextNNLabel( space->elementPtr() ),
    M_nextNNDistance( space->elementPtr() )
{
#if defined( FM_EXPORT )
    M_ex = exporter(_mesh=this->functionSpace()->mesh(), _name="fastmarchin");
    M_count_iteration = 0;
#endif
}

LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS
typename LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::self_ptrtype
LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::New(
        functionspace_ptrtype const& space,
        periodicity_type periodicity )
{
    self_ptrtype fm( new self_type(space, periodicity) );
    return fm;
}

LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS
void
LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::initMarch( 
        element_type const& phi,
        bool useMarker2AsMarkerDone )
{
    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    *(this->M_distance) = vf::project(
            this->functionSpace(), 
            elements(this->functionSpace()->mesh()), 
            idv(phi) );

    *M_NNDistance = vf::project(
            this->functionSpace(), 
            elements(this->functionSpace()->mesh()), 
            idv(phi) );

    M_label->setConstant( invalid_uint16_type_value );
    //*M_label = vf::project(
            //this->functionSpace(), 
            //elements(this->functionSpace()->mesh()), 
            //vf::cst(invalid_uint16_type_value ) );
    M_labelDist->setConstant( 1e8 );
    M_NNLabel->setConstant( invalid_uint16_type_value );
    M_nextNNLabel->setConstant( invalid_uint16_type_value );

    //*(this->M_status) = vf::project( 
            //this->functionSpace(), 
            //elements(this->functionSpace()->mesh()), 
            //vf::cst(super_type::FAR) );
    this->M_status->setConstant( super_type::FAR );

    std::set<size_type> doneIds;

    auto rangeElements = this->functionSpace()->mesh()->elementsWithProcessId();
    auto it = std::get<0>( rangeElements );
    auto en = std::get<1>( rangeElements );

    for ( ; it!=en ; ++it )
    {
        auto const& elt = boost::unwrap_ref( *it );

        std::vector<size_type> idsPhiPlus;
        std::vector<size_type> idsPhiMinus;
        uint16_type L0PhiMinus = invalid_uint16_type_value;
        uint16_type L0PhiPlus = invalid_uint16_type_value;
        for ( uint16_type j = 0; j < ndofv; ++j )
        {
            size_type index = this->functionSpace()->dof()->localToGlobal(elt, j, 0).index();
            if ( phi[index] < 0.0 )
            {
                idsPhiMinus.push_back( index );
                L0PhiMinus = (*M_selfLabel)[index];
            }
            else if ( phi[index] > 0.0 )
            {
                idsPhiPlus.push_back( index );
                L0PhiPlus = (*M_selfLabel)[index];
            }
        }
        //if the element has nodes with positive and negative values
        //mark as done
        if ( idsPhiPlus.size() != ndofv && idsPhiMinus.size() != ndofv )
        {
            for ( uint16_type j = 0; j < idsPhiPlus.size(); ++j )
            {
                doneIds.insert( idsPhiPlus[j] );
                (*this->M_status)[idsPhiPlus[j]] = super_type::DONE;
                (*M_NNLabel)[idsPhiPlus[j]] = L0PhiMinus;
                (*M_label)[idsPhiPlus[j]] = L0PhiMinus;
            }
            for ( uint16_type j = 0; j < idsPhiMinus.size(); ++j )
            {
                doneIds.insert( idsPhiMinus[j] );
                (*this->M_status)[idsPhiMinus[j]] = super_type::DONE;
                (*M_NNLabel)[idsPhiMinus[j]] = L0PhiPlus;
                (*M_label)[idsPhiMinus[j]] = L0PhiPlus;
            }
        }
    }

    // communicate the DONE list between all the proc
    this->reduceDonePoints( doneIds );
    if (Environment::worldComm().size() > 1)
    {
        sync( *M_NNLabel, "min" );
        sync( *M_label, "min" );
    }

    // initialize close distances in heap and mark close points in (*M_status) array
    for ( auto dit = doneIds.begin(); dit != doneIds.end(); ++dit )
       this->updateHeap( *dit );

#if defined( FM_EXPORT )
    M_count_iteration++;
    M_ex->step(M_count_iteration)->addRegions( "pid", "pid" );
    M_ex->step(M_count_iteration)->add("distance", *(this->M_distance));
    M_ex->step(M_count_iteration)->add("label", *M_label);
    M_ex->step(M_count_iteration)->add("labelDist", *M_labelDist);
    M_ex->step(M_count_iteration)->add("status", *(this->M_status));
    M_ex->step(M_count_iteration)->add("NNlabel", *M_NNLabel);
    M_ex->step(M_count_iteration)->add("NNdistance", *M_NNDistance);
    M_ex->step(M_count_iteration)->add("nextNNlabel", *M_nextNNLabel);
    M_ex->step(M_count_iteration)->add("nextNNdistance", *M_nextNNDistance);
    M_ex->save();
#endif
}

LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS
void
LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::processDof( size_type idOnProc, value_type val, value_type const& opt_data )
{
#if defined( FM_EXPORT )
    std::cout << "Processing dof " << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster( idOnProc )
        << " with value " << val << " and label " << opt_data
        << " on proc " << this->functionSpace()->worldComm().localRank() << "\n";
#endif

    if( (*M_NNLabel)[idOnProc] == invalid_uint16_type_value )
    {
        (*M_NNDistance)[idOnProc] = val;
        //(*M_NNLabel)[idOnProc] = (*M_label)[idOnProc];
        (*M_NNLabel)[idOnProc] = opt_data;
        (*M_labelDist)[idOnProc] = 1e8;
    }
    else
    {
        (*M_nextNNDistance)[idOnProc] = val;
        //(*M_nextNNLabel)[idOnProc] = (*M_label)[idOnProc];
        (*M_nextNNLabel)[idOnProc] = opt_data;
        this->setDofStatus( idOnProc, super_type::DONE );
    }

    this->setDofDistance(idOnProc, val);
    (*M_label)[idOnProc] = opt_data;
}

LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS
void
LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::updateHeap( size_type idDone )
{
    // nbr : neighbours of the node just done
    std::set<size_type> const & nbrs = this->neighbors().find(idDone)->second;

    std::vector<size_type> ids( 1, idDone );
    for ( auto n0it = nbrs.begin(); n0it != nbrs.end(); ++n0it )
    {
        if( (*M_label)[idDone] != (*M_NNLabel)[*n0it] 
                && (*M_label)[idDone] != (*M_nextNNLabel)[*n0it]
                && (*M_label)[idDone] != (*M_selfLabel)[*n0it] )
        {
            //if (this->getDofStatus(*n0it) == super_type::CLOSE )
            if (this->getDofStatus(*n0it) == super_type::FAR )
                this->setDofStatus(*n0it, super_type::CLOSE);

            bool hasNNLabel = ( (*M_NNLabel)[*n0it] != invalid_uint16_type_value );
            bool hasNextNNLabel = ( (*M_nextNNLabel)[*n0it] != invalid_uint16_type_value );
            if( !hasNNLabel || !hasNextNNLabel )
            {
                /* to give a reference, compute phi with only one DONE neighbours
                   it is sure that *n0it is a DONE neighbours */
                // one neighbor
                ids.push_back(*n0it);
                value_type phiNew = this->fmsDistN( ids, *(this->M_distance) );
                ids.pop_back();

                bool updateLabel = false;
                if( std::abs(phiNew) < std::abs((*M_labelDist)[*n0it]) )
                {
                    (*M_labelDist)[*n0it] = phiNew;
                    (*M_label)[*n0it] = (*M_label)[idDone];
                    updateLabel = true;
                }

                /*compute all the phi possible with all the neighbors around 
                 * and returns the smallest one*/
                if( !hasNNLabel )
                    phiNew = this->fmsNNDistRec( ids, *n0it, phiNew );
                else if( !hasNextNNLabel )
                    phiNew = this->fmsNextNNDistRec( ids, *n0it, phiNew );

                value_type label = (*M_label)[idDone];

#if defined( FM_EXPORT )
                std::cout << "\tUpdating dof " << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster( *n0it )
                    << " with value " << phiNew << " and label " << label
                    << " on proc " << this->functionSpace()->worldComm().localRank() << "\n";
#endif

                //this->heap().change( std::make_pair( phiNew, *n0it ), std::vector<value_type>(1, label) );
                this->heap().change( std::make_pair( phiNew, *n0it ) );
                if(updateLabel)
                {
                    this->heap().dataAtIndex( *n0it ) = label;
                }
            } // if CLOSE
        }
    } // loop over neighbor 0

#if defined( FM_EXPORT )
    M_count_iteration++;
    M_ex->step(M_count_iteration)->addRegions( "pid", "pid" );
    M_ex->step(M_count_iteration)->add("distance", *(this->M_distance));
    M_ex->step(M_count_iteration)->add("label", *M_label);
    M_ex->step(M_count_iteration)->add("labelDist", *M_labelDist);
    M_ex->step(M_count_iteration)->add("status", *(this->M_status));
    M_ex->step(M_count_iteration)->add("NNlabel", *M_NNLabel);
    M_ex->step(M_count_iteration)->add("NNdistance", *M_NNDistance);
    M_ex->step(M_count_iteration)->add("nextNNlabel", *M_nextNNLabel);
    M_ex->step(M_count_iteration)->add("nextNNdistance", *M_nextNNDistance);
    M_ex->save();
#endif
}

LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS
typename LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::value_type
LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::fmsNNDistRec( 
        std::vector<size_type> & ids,
        size_type idClose,
        value_type phiOld ) const
{
    /*search for all neighbors with L1 done
      recalculate phi with all the neighbors using NNDistance
      returns the smallest phi  */

    // only allows for getting at maximum 2 nodes in 2d and 3 nodes in 3d
    if ( ids.size() >= Dim )
        return phiOld;

    value_type phiNew(phiOld);

    std::set<size_type> const & nbrs = this->M_neighbors.find(idClose)->second;
    for (auto nit = nbrs.begin(); nit != nbrs.end(); ++nit )
    {
        // if the next neighbors is not a neighbor with L1 done
        if( (*M_NNLabel)[*nit] == invalid_uint16_type_value )
            continue;

        // if this node has already been accepted, stop
        bool unique = true;
        for ( auto idsit=ids.begin(); idsit != ids.end(); ++idsit )
            unique &= ( *idsit != *nit );
        if ( !unique ) // points must be unique
            continue;

        if ( this->M_neighbors.find(ids[0])->second.find(*nit) ==
                this->M_neighbors.find(ids[0])->second.end() )
            continue;

        // one neighbor more
        ids.push_back(*nit);

        ids.push_back(idClose);
        value_type phiCand = this->fmsDistN( ids, *M_NNDistance );
        ids.pop_back();
        if ( phiCand != 0.0 )
            phiNew = super_type::closerOne( phiCand, phiNew );

        phiNew = fmsNNDistRec( ids, idClose, phiNew );

        ids.pop_back();
    }
    return phiNew;
}

LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS
typename LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::value_type
LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::fmsNextNNDistRec( 
        std::vector<size_type> & ids,
        size_type idClose,
        value_type phiOld ) const
{
    /*search for all neighbors with L2 done
      recalculate phi with all the neighbors using nextNNDistance
      returns the smallest phi  */

    // only allows for getting at maximum 2 nodes in 2d and 3 nodes in 3d
    if ( ids.size() >= Dim )
        return phiOld;

    value_type phiNew(phiOld);

    std::set<size_type> const & nbrs = this->M_neighbors.find(idClose)->second;
    for (auto nit = nbrs.begin(); nit != nbrs.end(); ++nit )
    {
        // if the next neighbors is not a neighbor with L2 done
        if( (*M_nextNNLabel)[*nit] == invalid_uint16_type_value )
            continue;

        // if this node has already been accepted, stop
        bool unique = true;
        for ( auto idsit=ids.begin(); idsit != ids.end(); ++idsit )
            unique &= ( *idsit != *nit );
        if ( !unique ) // points must be unique
            continue;

        if ( this->M_neighbors.find(ids[0])->second.find(*nit) ==
                this->M_neighbors.find(ids[0])->second.end() )
            continue;

        // one neighbor more
        ids.push_back(*nit);

        ids.push_back(idClose);
        value_type phiCand = this->fmsDistN( ids, *M_nextNNDistance );
        ids.pop_back();
        if ( phiCand != 0.0 )
            phiNew = super_type::closerOne( phiCand, phiNew );

        phiNew = fmsNextNNDistRec( ids, idClose, phiNew );

        ids.pop_back();
    }
    return phiNew;
}

#undef LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE
#undef LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS

} // namespace Feel

#endif
