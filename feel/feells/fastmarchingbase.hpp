#ifndef _FASTMARCHINGBASE_HPP
#define _FASTMARCHINGBASE_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>

#include <feel/feells/fmsheap.hpp>
#include <feel/feells/fmspoint.hpp>
#include <feel/feells/lstypes.hpp>
#include <boost/serialization/set.hpp>

namespace Feel {

template<typename FunctionSpaceType, typename PeriodicityType = NoPeriodicity>
class FastMarchingBase
{
public:
    static_assert( FunctionSpaceType::fe_type::nOrder == 1, "FunctionSpaceType needs to be a finite element space of order 1");
    static_assert( ! FunctionSpaceType::is_periodic , "Space for fast marching must be non periodic, but periodicity can be given as second template argument");

    //--------------------------------------------------------------------//
    // Typedefs
    typedef FastMarchingBase<FunctionSpaceType, PeriodicityType> self_type;

    enum status_type {FAR=0, CLOSE=1, DONE=2};

    typedef Backend<double> backend_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef FunctionSpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef typename functionspace_type::value_type value_type;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type geoelement_type;
    static const uint16_type Dim = geoelement_type::nDim;

    typedef PeriodicityType periodicity_type;
    typedef typename periodicity_type::node_type node_type;
    //--------------------------------------------------------------------//
    // Constructor
    FastMarchingBase( functionspace_ptrtype const& space,
            periodicity_type periodicity = NoPeriodicity() );

    virtual ~FastMarchingBase() = default;

    //--------------------------------------------------------------------//
    // Run march
    void run( element_type const& phi, bool useMarker2AsMarkerDone =false );

    //--------------------------------------------------------------------//
    element_ptrtype getDistance() const { return M_distance; }

protected:
    /* The map indexes are the global indexes on the PROC and the set contains 
     * the global indexes on the CLUSTER of its neigbors */
    typedef std::map<size_type, std::set<size_type> > neighbors_type;
    typedef Feel::details::FmsHeap<value_type> heap_type;
    typedef typename heap_type::heap_entry_type heap_entry_type;
    typedef Feel::details::FmsPoint<value_type, Dim> point_type;
    //--------------------------------------------------------------------//

    static value_type closerOne( value_type a, value_type b ) { return (a*a < b*b)? a: b; }

    functionspace_ptrtype const& functionSpace() const { return M_functionspace; }

    void setDofStatus( size_type idOnProc, status_type s ) { (*M_status)[idOnProc] = s; }
    status_type getDofStatus( size_type idOnProc ) const { return static_cast<status_type>( (*M_status)[idOnProc] ); }

    void setDofDistance( size_type idOnProc, value_type v ) { (*M_distance)[idOnProc] = v; }
    value_type getDofDistance( size_type idOnProc ) { return (*M_distance)[idOnProc]; } 

    neighbors_type const& neighbors() const { return M_neighbors; }

    heap_type & heap() { return M_heap; }

    //--------------------------------------------------------------------//
    size_type clusterToProcessor( size_type dof ) const;
    size_type processorToCluster( size_type dof ) const;

    void createPeriodicCorrespondanceTable();

    void reduceDonePoints( std::set<size_type>& doneIds );
    void reduceClosePoints();

    //--------------------------------------------------------------------//

    value_type fmsDistN( std::vector<size_type> const& ids ) const;

    value_type fmsDistRec( std::vector<size_type> & ids,
                           size_type idClose,
                           value_type phiOld ) const;
                    
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    virtual void initMarch( element_type const& phi, bool useMarker2AsMarkerDone );

    virtual void processDof( size_type idOnProc, value_type val ) =0;

    virtual void updateHeap( size_type idDone );

protected:
    //--------------------------------------------------------------------//
    functionspace_ptrtype const& M_functionspace;
    element_ptrtype M_status;
    heap_type M_heap;
    vector_ptrtype M_checkStatus;
    vector_ptrtype M_valueAtClose;
    neighbors_type M_neighbors;
    element_ptrtype M_distance;

    std::map< size_type, size_type> M_ghostClusterToProc;

    periodicity_type M_periodicity;
    boost::bimap< size_type, size_type > M_idTag1_idTag2;
    std::vector<point_type> M_coords;
    node_type M_translation;

    const size_type M_firstDof;
    int M_nbDofTag1;
    int M_nbTotalDone;
};


} // namespace Feel

#include <feel/feells/fastmarchingbase_impl.hpp>

#endif
