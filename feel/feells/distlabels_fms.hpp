#ifndef _DISTLABELS_FMS_HPP
#define _DISTLABELS_FMS_HPP 1

#define FM_EXPORT 1

#include <feel/feells/fastmarchingbase.hpp>
#include <feel/feells/selflabel.hpp>

namespace Feel {

template<typename FunctionSpaceType, typename PeriodicityType =NoPeriodicity>
class LabelDistanceFMS : public FastMarchingBase<FunctionSpaceType, PeriodicityType>
{
public:
    typedef FastMarchingBase<FunctionSpaceType, PeriodicityType> super_type;
    typedef LabelDistanceFMS<FunctionSpaceType, PeriodicityType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef typename super_type::functionspace_type functionspace_type;
    typedef typename super_type::functionspace_ptrtype functionspace_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::element_ptrtype element_ptrtype;
    typedef typename super_type::value_type value_type;

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
    element_ptrtype getNearestNeighbourLabel() const { return M_NNLabel; }
    element_ptrtype getNearestNeighbourDistance() const { return this->getDistance(); }
    element_ptrtype getNextNearestNeighbourLabel() const { return M_nextNNLabel; }
    element_ptrtype getNextNearestNeighbourDistance() const { return M_nextNNDistance; }

private:
    void initMarch( element_type const& phi, bool useMarker2AsMarkerDone );
    void processDof( size_type idOnProc, value_type val );
    void updateHeap( size_type idDone );

    element_ptrtype getDistance() const { return super_type::getDistance(); }

    //--------------------------------------------------------------------//
    element_ptrtype M_selfLabel;
    element_ptrtype M_label;
    element_ptrtype M_NNLabel;
    element_ptrtype M_NNDistance;
    element_ptrtype M_nextNNLabel;
    element_ptrtype M_nextNNDistance;

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
    //auto doneIds = super_type::initMarch( phi, useMarker2AsMarkerDone );

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
    M_NNLabel->setConstant( invalid_uint16_type_value );
    M_nextNNLabel->setConstant( invalid_uint16_type_value );

    //*(this->M_status) = vf::project( 
            //this->functionSpace(), 
            //elements(this->functionSpace()->mesh()), 
            //vf::cst(super_type::FAR) );
    this->M_status->setConstant( super_type::FAR );

    std::set<size_type> doneIds;

    auto it = this->functionSpace()->mesh()->beginElementWithProcessId();
    auto en = this->functionSpace()->mesh()->endElementWithProcessId();

    for ( ; it!=en ; ++it )
    {
        std::vector<size_type> idsPhiPlus;
        std::vector<size_type> idsPhiMinus;
        uint16_type L0PhiMinus = invalid_uint16_type_value;
        uint16_type L0PhiPlus = invalid_uint16_type_value;
        for ( uint16_type j = 0; j < ndofv; ++j )
        {
            size_type index = this->functionSpace()->dof()->localToGlobal(*it, j, 0).index();
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

    // initialize close distances in heap and mark close points in (*M_status) array
    for ( auto dit = doneIds.begin(); dit != doneIds.end(); ++dit )
        this->updateHeap( *dit );

#if defined( FM_EXPORT )
    M_count_iteration++;
    M_ex->step(M_count_iteration)->add("label", *M_label);
    M_ex->step(M_count_iteration)->add("status", *(this->M_status));
    M_ex->step(M_count_iteration)->add("NNlabel", *M_NNLabel);
    M_ex->step(M_count_iteration)->add("NNdistance", *M_NNDistance);
    M_ex->save();
#endif
}

LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS
void
LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE::processDof( size_type idOnProc, value_type val )
{
    if( (*M_NNLabel)[idOnProc] == invalid_uint16_type_value )
    {
        (*M_NNDistance)[idOnProc] = val;
        (*M_NNLabel)[idOnProc] = (*M_label)[idOnProc];
    }
    else
    {
        (*M_nextNNDistance)[idOnProc] = val;
        (*M_nextNNLabel)[idOnProc] = (*M_label)[idOnProc];
        this->setDofStatus( idOnProc, super_type::DONE );
    }

    this->setDofDistance(idOnProc, val);
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
                && (*M_label)[idDone] != (*M_nextNNLabel)[*n0it] )
        {
            if (this->getDofStatus(*n0it) == super_type::FAR )
               this->setDofStatus(*n0it, super_type::CLOSE);
            if (this->getDofStatus(*n0it) == super_type::CLOSE )
            {
                /* to give a reference, compute phi with only one DONE neighbours
                   it is sure that *n0it is a DONE neighbours */
                // one neighbor
                ids.push_back(*n0it);
                value_type phiNew = this->fmsDistN( ids );
                ids.pop_back();

                /*compute all the phi possible with all the neighbors around and returns the smallest one*/
                phiNew = this->fmsDistRec( ids, *n0it, phiNew );

                this->heap().change( std::make_pair( phiNew, *n0it ) );

                // update M_label
                (*M_label)[*n0it] = (*M_NNLabel)[idDone];
            } // if CLOSE
        }
    } // loop over neighbor 0

#if defined( FM_EXPORT )
    M_count_iteration++;
    M_ex->step(M_count_iteration)->add("label", *M_label);
    M_ex->step(M_count_iteration)->add("status", *(this->M_status));
    M_ex->step(M_count_iteration)->add("NNlabel", *M_NNLabel);
    M_ex->step(M_count_iteration)->add("NNdistance", *M_NNDistance);
    M_ex->save();
#endif
}

#undef LABELDISTANCEFMS_CLASS_TEMPLATE_TYPE
#undef LABELDISTANCEFMS_CLASS_TEMPLATE_DECLARATIONS

//template<typename FunctionSpaceType, typename PeriodicityType = NoPerio

} // namespace Feel

#endif
