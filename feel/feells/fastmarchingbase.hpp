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

template<typename FunctionSpaceType, typename PeriodicityType = NoPeriodicity, typename HeapDataType = Feel::details::none_type>
class FastMarchingBase
{
public:
    static_assert( ! FunctionSpaceType::is_periodic , "Space for fast marching must be non periodic, but periodicity can be given as second template argument");

    //--------------------------------------------------------------------//
    // Typedefs
    typedef FastMarchingBase<FunctionSpaceType, PeriodicityType, HeapDataType> self_type;

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

    typedef HeapDataType heap_data_type;
    //--------------------------------------------------------------------//
    // Constructor
    FastMarchingBase( functionspace_ptrtype const& space,
            periodicity_type periodicity = NoPeriodicity() );

    virtual ~FastMarchingBase() = default;

    //--------------------------------------------------------------------//
    // Run march
    void run( element_type const& phi, bool useMarker2AsMarkerDone =false );
    void run( element_ptrtype const& phi, bool useMarker2AsMarkerDone =false )
    {
        this->run( *phi, useMarker2AsMarkerDone );
    }

    //--------------------------------------------------------------------//
    element_ptrtype getDistance() const { return M_distance; }

protected:
    //--------------------------------------------------------------------//
    template<typename T, typename DataT>
    struct ValueDataType : std::pair<T, DataT>
    {
        typedef T value_type;
        typedef DataT data_type;

        typedef typename std::pair<T, DataT> super_type;
        
        using typename super_type::pair;
        ValueDataType(): super_type() {}
        ValueDataType(T const& v, DataT const& d ): super_type(v,d) {}

        value_type & value() { return this->first; }
        value_type const& value() const { return this->first; }
        data_type data() const { return this->second; }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->first;
            ar & this->second;
        }
    };
    template<typename T>
    struct ValueDataType<T, Feel::details::none_type>
    {
        typedef T value_type;
        typedef Feel::details::none_type data_type;

        ValueDataType() {}
        ValueDataType(T const& v, data_type d = data_type() ): M_value(v) {}

        operator T() { return M_value; }

        value_type & value() { return M_value; }
        value_type const& value() const { return M_value; }
        data_type data() const { return data_type(); }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & M_value;
        }

        value_type M_value;
    };

    typedef ValueDataType<value_type, heap_data_type> value_data_type;
    //--------------------------------------------------------------------//
    /* The map indexes are the global indexes on the PROC and the set contains 
     * the global indexes on the CLUSTER of its neigbors */
    typedef std::map<size_type, std::set<size_type> > neighbors_type;
    typedef Feel::details::FmsHeap<value_type, heap_data_type> heap_type;
    typedef typename heap_type::heap_entry_type heap_entry_type;
    typedef typename heap_type::heap_entry_data_type heap_entry_data_type;
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

    value_type fmsDistN( std::vector<size_type> const& ids, element_type const & __v ) const;

    value_type fmsDistRec( std::vector<size_type> & ids,
                           size_type idClose,
                           value_type phiOld ) const;
                    
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    virtual void initMarch( element_type const& phi, bool useMarker2AsMarkerDone );

    virtual void processDof( size_type idOnProc, value_type val, heap_data_type const& opt_data ) =0;

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
