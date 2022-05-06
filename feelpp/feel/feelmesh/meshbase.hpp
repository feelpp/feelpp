//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 05 Feb 2017
//! @copyright 2005-2006 EPFL
//! @copyright 2011-2017 Feel++ Consortium
//!
#ifndef FEELPP_MESHBASE_HPP
#define FEELPP_MESHBASE_HPP 1

#include <unordered_map>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/context.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelmesh/submeshdata.hpp>
#include <feel/feelmesh/meshsupportbase.hpp>

#if defined(FEELPP_HAS_VTK)
#include <feel/feelcore/disablewarnings.hpp>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <feel/feelcore/reenablewarnings.hpp>
#endif

namespace Feel
{
template <typename IndexT> class SubMeshData;

//!
//! @brief Mesh components enum
//! @ingroup Mesh
//! 
//! Components of a mesh that can be enabled or disabled when calling
//! \c updateForUse()
//!
enum MeshComponents
{
    MESH_UPDATE_EDGES     = ( 1 << 0 ),
    MESH_UPDATE_FACES     = ( 1 << 1 ),
    MESH_CHECK            = ( 1 << 2 ),
    MESH_PARTITION        = ( 1 << 3 ),
    MESH_RENUMBER         = ( 1 << 4 ),
    MESH_ADD_ELEMENTS_INFO = ( 1 << 5 ),
    MESH_PROPAGATE_MARKERS = ( 1 << 6 ),
    MESH_REMOVE_PERIODIC_FACES_FROM_BOUNDARY = ( 1 << 7 ),
    MESH_NO_UPDATE_MEASURES = ( 1 << 8 ),
    MESH_UPDATE_ELEMENTS_ADJACENCY = ( 1 << 9 ),
    MESH_UPDATE_FACES_MINIMAL = ( 1 << 10 ),
    MESH_GEOMAP_NOT_CACHED = ( 1 << 11 ),
    MESH_DO_NOT_UPDATE = (1 << 16 )

};
const uint16_type MESH_ALL_COMPONENTS = MESH_UPDATE_EDGES | MESH_UPDATE_FACES | MESH_CHECK | MESH_PARTITION | MESH_RENUMBER;
const uint16_type MESH_COMPONENTS_DEFAULTS = MESH_RENUMBER | MESH_CHECK;

//!
//! @brief base mesh class
//! @ingroup Mesh
//!
//! @author Christophe Prud'homme
//! @see
//!/
template<typename IndexT = uint32_type>
class FEELPP_EXPORT MeshBase : public CommObject, public JournalWatcher
{
    using super2 = JournalWatcher;
public:


    /** @name Typedefs
     */
    //@{
    using super = CommObject;

    //!
    //! type of indices : default to uint32_type
    //!
    using index_type = IndexT;
    using index_t = index_type;
    using size_type = index_type;

    /**
     * Tuple that contains
     *
     * -# the index of the face
     *
     * -# the processor id the face belongs to
     */
    typedef boost::tuple<index_type, index_type> face_processor_type;


    typedef SubMeshData<index_type> smd_type;
    typedef std::shared_ptr<smd_type> smd_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    MeshBase()
        :
        MeshBase( 3, 3, Environment::worldCommPtr() )
    {}
    explicit MeshBase( worldcomm_ptr_t const& w )
        :
        MeshBase ( 3, 3, w )
    {}
    MeshBase( MeshBase const& ) = default;
    MeshBase( MeshBase && ) = default;
    ~MeshBase() override;

    /**
     * build from a topological dimension, a real dimension and a communicator
     */
    MeshBase( uint16_type topodim, uint16_type realdim,
              worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        :
        MeshBase( "", topodim, realdim, worldComm )
    {}

    /**
     * build from a name, a topological dimension, a real dimension and a communicator
     */
    MeshBase( std::string const& name, uint16_type topodim, uint16_type realdim,
              worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

    //@}

    /** @name Operator overloads
     */
    //@{

    MeshBase& operator=( MeshBase const& m ) = default;
    MeshBase& operator=( MeshBase && m ) = default;

    //@}

    /** @name Accessors
     */
    //@{

    //! return current shared_ptr of type MeshBase
    virtual std::shared_ptr<MeshBase<IndexT>> shared_from_this_meshbase() = 0;

    //! return current shared_ptr of type MeshBase
    virtual std::shared_ptr<const MeshBase<IndexT>> shared_from_this_meshbase() const = 0;

    /**
     * \return \p true if the mesh is ready for use, \p false
     * otherwise
     */
    bool isUpdatedForUse() const
    {
        return M_is_updated;
    }

    //!
    //! set the topological dimension
    //!
    void setTopologicalDimension( int d ) { M_topodim = d; }

    //!
    //! set the real dimension
    //!
    void setRealDimension( int d ) { M_realdim = d; }
    
    /**
     * \return the number of elements
     */
    virtual size_type numGlobalElements() const = 0;
    /**
     * \return the number of faces
     */
    virtual size_type numGlobalFaces() const = 0;
    /**
     * \return the number of edges
     */
    virtual size_type numGlobalEdges() const = 0;
    /**
     * \return the number of points
     */
    virtual size_type numGlobalPoints() const = 0;
    /**
     * \return the number of vertices
     */
    virtual size_type numGlobalVertices() const = 0;
    /**
     * \return the number of elements (in current process)
     */
    virtual size_type numElements() const = 0;

    /**
     * \return the number of faces (in current process)
     */
    virtual size_type numFaces() const = 0;

    /**
     * \return the number of Points (in current process)
     */
    virtual size_type numPoints() const = 0;

    /**
     * \return the number of vertices (in current process)
     */
    size_type numVertices() const
    {
        return M_n_vertices;
    }

    /**
     * Returns the number of partitions.
     */
    rank_type numberOfPartitions() const
    {
        return M_n_parts;
    }

    /**
     * \return \c true if mesh is partitioned, \c false otherwise
     */
    bool isPartitioned() const;

    /**
     * \return an integer(stored  in a \p Context) that encodes the components to be updated by
     * the mesh data structure.
     * \sa Context, MeshComponents
     */
    Context const& components() const
    {
        return M_components;
    }

    /**
     * \return an integer(stored  in a \p Context) that encodes the components to be updated by
     * the mesh data structure.
     * \sa Context, MeshComponents
     */
    Context&       components()
    {
        return M_components;
    }

    /**
     * \return true if the mesh has parametric nodes
     */
    bool isParametric() const
    {
        return M_is_parametric;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //! set the number of partitions 
    void setNumberOfPartitions( rank_type n )
    {
        M_n_parts = n;
    }

    //! set the number of vertices 
    void setNumVertices( size_type n )
    {
        M_n_vertices = n ;
    }

    /**
     * set the components to be updated by \c updateForUse()
     * \sa updateForUse
     */
    void setComponents( size_type components = MESH_ALL_COMPONENTS )
    {
        M_components = components;
    }

    /**
     * set if the mesh is parametric ( e.g. has parametric nodes )
     */
    void setParametric( bool x )
    {
        M_is_parametric = x;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Empty all containers of the mesh. Must be redefined by derived
     * classes.
     */
    virtual void clear();


    /**
     * update all info for this mesh.
     */
    virtual void updateForUse() = 0;

    /**
     * update all info for this mesh according the \p components
     * \sa components(), setComponents()
     */
    virtual void updateForUse( size_type components );

    //! update the mesh when nodes have moved
    virtual void updateForUseAfterMovingNodes( bool upMeasures = true ) = 0;

    /**
     * Call the default partitioner (currently \p metis_partition()).
     */
    virtual void partition ( const rank_type n_parts ) = 0;

    virtual void meshModified() = 0;

    //! set sub mesh data
    void setSubMeshData( smd_ptrtype smd )
        {
            M_smd = smd;
        }

    //! return sub mesh data
    smd_ptrtype subMeshData() const { return M_smd; }

    //! \return true if mesh holds sub mesh data
    bool hasSubMeshData() const { return M_smd.use_count() > 0; }

    //! \return sub mesh
    typename smd_type::mesh_ptrtype subMesh() const
        {
            CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
            return M_smd->mesh;
        }

    //! \return sub mesh
    typename smd_type::mesh_ptrtype parentMesh() const
        {
            CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
            return M_smd->mesh;
        }

    //! true if is a sub mesh
    bool isSubMesh() const { return !M_smd == false; }

    //! \return true if the mesh is related to the mesh \p m
    bool isSubMeshFrom( MeshBase const* m ) const
        {
            DVLOG(4) << "isSubMeshFrom<mesh_ptrtype> called\n";
            if ( !M_smd ) return false;
            bool res= (M_smd->mesh.get() == m);
            DVLOG(4) << "this isSubMeshFrom m: " << res << "\n";
            return res;
        }
    //! \return true if the mesh is related to the mesh \p m
    bool isSubMeshFrom( std::shared_ptr<MeshBase> m ) const
        {
            return isSubMeshFrom( m.get() );
        }

    bool isParentMeshOf( MeshBase const* m ) const
        {
            DVLOG(4) << "isParentMeshOf<mesh_ptrtype> called\n";
            bool res = m->isSubMeshFrom( this );
            if ( res == false ) return res;
            DVLOG(4) << "this isParentMeshOf m: " << res << "\n";
            return res;
        }
    //! \return true if the mesh is related to the mesh \p m
    bool isParentMeshOf( std::shared_ptr<MeshBase> m ) const
        {
            DVLOG(4) << "isParentMeshOf<mesh_ptrtype> called\n";
            bool res = m->isSubMeshFrom( this );
            if ( res == false ) return res;
            DVLOG(4) << "this isParentMeshOf m: " << res << "\n";
            return res;
        }
    bool isSiblingOf( MeshBase const* m ) const
        {
            DVLOG(4) << "isSibling<mesh_ptrtype> called\n";
            if ( !M_smd || !m->hasSubMeshData() ) return false;
            bool res = M_smd->mesh.get() == m->M_smd->mesh.get();
            if ( res == false ) return res;
            DVLOG(4) << "this isSibling m: " << res << "\n";
            return res;
        }
    //! \return true if the mesh is related to the mesh \p m
    bool isSiblingOf( std::shared_ptr<MeshBase> m ) const
        {
            DVLOG(4) << "isSibling<mesh_ptrtype> called\n";
            if ( !M_smd || !m->hasSubMeshData() ) return false;
            bool res = M_smd->mesh.get() == m->M_smd->mesh.get();
            if ( res == false ) return res;
            DVLOG(4) << "this isSibling m: " << res << "\n";
            return res;
        }
#if 0
    template<typename M>
    bool isSubMeshFrom( std::shared_ptr<M> m ) const
        {
            DVLOG(4) << "isSubMeshFrom<M> called\n";
            return false;
        }
#endif
    template<typename M>
    bool isSameMesh( M const* m ) const
        {
            bool same_mesh = ( dynamic_cast<void const*>( this ) == dynamic_cast<void const*>( m ) );
            return same_mesh;
        }

    template<typename M>
    bool isSameMesh( std::shared_ptr<M> m ) const
        {
            bool same_mesh = ( dynamic_cast<void const*>( this ) == dynamic_cast<void*>( m.get() ) );
            return same_mesh;
        }
    template<typename M>
    bool isRelatedTo( M const* m ) const
        {
            bool same_mesh = isSameMesh(m);
            DVLOG(4) << "same_mesh: " << same_mesh << "\n";
            bool is_submesh_from = isSubMeshFrom( m );
            DVLOG(4) << "isSubMeshFrom: " << is_submesh_from << "\n";
            bool is_parentmesh_of = isParentMeshOf( m );
            DVLOG(4) << "is_parentmesh_of: " << is_parentmesh_of << "\n";
            bool is_sibling_of = isSiblingOf( m );
            DVLOG(4) << "is_sibling_of: " << is_sibling_of << "\n";
            return same_mesh || is_submesh_from || is_parentmesh_of || is_sibling_of;
            //return same_mesh || is_submesh_from || is_parentmesh_of;
        }
    template<typename M>
    bool isRelatedTo( std::shared_ptr<M> m ) const
        {
            bool same_mesh = isSameMesh(m);
            DVLOG(4) << "same_mesh: " << same_mesh << "\n";
            bool is_submesh_from = isSubMeshFrom( m );
            DVLOG(4) << "isSubMeshFrom: " << is_submesh_from << "\n";
            bool is_parentmesh_of = isParentMeshOf( m );
            DVLOG(4) << "is_parentmesh_of: " << is_parentmesh_of << "\n";
            bool is_sibling_of = isSiblingOf( m );
            DVLOG(4) << "is_sibling_of: " << is_sibling_of << "\n";
            return same_mesh || is_submesh_from || is_parentmesh_of || is_sibling_of;
            //return same_mesh || is_submesh_from || is_parentmesh_of;
        }

    //! \return id in parent mesh given the id in the sub mesh
    size_type subMeshToMesh( size_type id ) const
        {
            CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
            return M_smd->bm.left.find( id )->second;
        }

    //! \return id in sub mesh given the id in the parent mesh
    size_type meshToSubMesh( size_type id ) const
        {
            CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
            if ( M_smd->bm.right.find( id ) != M_smd->bm.right.end() )
                return M_smd->bm.right.find( id )->second;
            // the submesh element id has not been found, return invalid value
            return invalid_v<size_type>;
        }

    //! \return ids in sub mesh given the ids in the parent mesh
    std::pair<std::vector<size_type>,bool> meshToSubMesh( std::vector<size_type> const& p, bool add_invalid_indices = false ) const
        {
            CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
            std::vector<size_type> sid;//(std::distance(p.first,p.second) );
            bool has_invalid_values = false;
            std::for_each( p.begin(), p.end(), [&]( auto const& id ){
                    if ( M_smd->bm.right.find( id ) != M_smd->bm.right.end() )
                        sid.push_back( M_smd->bm.right.find( id )->second );
                    else
                    {
                        if ( add_invalid_indices )
                            sid.push_back( invalid_v<size_type> );
                        has_invalid_values = true;
                    }
                    // the submesh element id has not been found, return invalid value
                    //return invalid_v<size_type>;
                });
            return std::make_pair(sid,has_invalid_values);
        }

    //! \return id in parent mesh given the id in the sub mesh
    size_type subMeshToMesh( std::shared_ptr<MeshBase> m, size_type id ) const
        {
            if ( this == m.get() )
                return id;
            if ( isRelatedTo( m ) )
            {
                if ( this->isSubMeshFrom( m ) )
                {
                    CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
                    return M_smd->bm.left.find( id )->second;
                }
                else if ( this->isSiblingOf( m ) )
                {
                    size_type id_in_parent =  M_smd->bm.left.find( id )->second;
                    size_type id_in_sibling =  m->meshToSubMesh( id_in_parent );
                    return id_in_sibling;
                }
            }
            return invalid_v<size_type>;
        }

    //! \return id in sub mesh given the id in the parent mesh
    size_type meshToSubMesh( std::shared_ptr<MeshBase> m, size_type id ) const
        {
            if ( this == m.get() )
                return id;
            if ( isRelatedTo( m ) )
            {
                if ( this->isSubMeshFrom( m ) )
                {
                    CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
                    if ( M_smd->bm.right.find( id ) != M_smd->bm.right.end() )
                        return M_smd->bm.right.find( id )->second;

                }
                else if ( this->isSiblingOf( m ) )
                {
                    size_type id_in_parent =  m->subMeshToMesh( id );
                    size_type id_in_sibling =  this->meshToSubMesh( id_in_parent );
                    return id_in_sibling;
                }
                // the submesh element id has not been found, return invalid value
                // will return invalid_v<size_type>
            }
            return invalid_v<size_type>;
        }

    //! store a mesh which has nodes shared with the current mesh
    template<typename MeshType>
    void addMeshWithNodesShared( std::shared_ptr<MeshType> m )
    {
        auto foundMesh = std::find_if( M_meshesWithNodesShared.begin(), M_meshesWithNodesShared.end(),
                                       [&m](auto const& wptr) { return m == wptr.lock(); });
        if ( foundMesh == M_meshesWithNodesShared.end() )
        {
            M_meshesWithNodesShared.push_back( m );
            m->addMeshWithNodesShared( this->shared_from_this_meshbase() );
        }
    }

    //! store a mesh which has nodes shared with the current mesh
    template<typename MeshType>
    void addMeshWithNodesShared( std::shared_ptr<const MeshType> m )
    {
        this->addMeshWithNodesShared( std::const_pointer_cast<MeshType>( m ) );
    }

    //! get all meshes which share the nodes with the current mesh (current mesh is also include)
    std::vector<std::shared_ptr<MeshBase<IndexT>>>
    meshesWithNodesShared() const
    {
        std::vector<std::shared_ptr<MeshBase<IndexT>>> ret;
        auto selfPtr = std::const_pointer_cast< MeshBase<IndexT>>( this->shared_from_this_meshbase() );
        ret.push_back( selfPtr );
        meshesWithNodesSharedImpl( ret );
        return ret;
    }

    //! attach a mesh support
    void attachMeshSupport( MeshSupportBase* ms ) const;

    //! detach a mesh support
    void detachMeshSupport( MeshSupportBase* ms ) const;

    //! return all mesh supports attached
    std::vector<MeshSupportBase*> const& meshSupportsAttached() const { return M_meshSupportsAttached; }

    //!
    //! @return true if the list strings are all  mesh marker
    //!
    bool hasAllMarkers( std::initializer_list<std::string> l )
        {
            bool f = true;
            for (auto n : l )
            {
                if ( !this->hasMarker( n ) )
                    f = false;
            }
            return f;
        }
    //!
    //! @return true if any of the list string contain a mesh marker
    //!
    bool hasAnyMarker( std::initializer_list<std::string> l )
        {
            for (auto n : l )
            {
                if ( this->hasMarker( n ) )
                    return true;
            }
            return false;
        }
    /**
     * @return true if \p marker exists, false otherwise
     */
    bool
    hasMarker( std::string const& marker ) const
        {
            return markerName( marker ) != invalid_v<size_type>;
        }

    /**
     * @return true if \p marker exists and topological dimension of the entity
     * associated is Dim, false otherwise
     */
    bool
    hasElementMarker( std::string const& marker ) const
        {
            return hasMarker( marker ) && ( markerDim( marker ) == M_topodim );
        }

    /**
     * @return true if \p marker exists and topological dimension of the entity
     * associated is Dim-1, false otherwise
     */
    bool
    hasFaceMarker( std::string const& marker ) const
        {
            return hasMarker( marker ) && ( markerDim( marker ) == M_topodim-1 );
        }

    /**
     * @return true if \p marker exists and topological dimension of the entity
     * associated is Dim-2, false otherwise
     */
    bool
    hasEdgeMarker( std::string const& marker ) const
        {
            return (M_topodim == 3) && hasMarker( marker ) &&  ( markerDim( marker ) == M_topodim-2 );
        }

    /**
     * @return true if \p marker exists and topological dimension of the entity
     * associated is 0, false otherwise
     */
    bool
    hasPointMarker( std::string const& marker ) const
        {
            return hasMarker( marker ) &&  ( markerDim( marker ) == 0 );
        }

    /**
     * @return true if all markers are defined in the mesh, false otherwise
     */
    bool hasMarkers( std::initializer_list<std::string> l ) const
        {
            for( auto m : l )
            {
                if ( !hasMarker( m ) ) return false;
            }
            return true;

        }


    /**
     * @return the id associated to the \p marker
     */
    size_type markerName( std::string const& marker ) const
    {
        auto mit = M_markername.find( marker );
        if (  mit != M_markername.end() )
            return mit->second[0];
        return invalid_v<size_type>;
    }
    /**
     * @return the marker name associated to the \p marker id
     */
    std::string markerName( size_type marker ) const
        {
            for( auto const& n : M_markername )
            {
                if (n.second[0] == marker )
                    return n.first;
            }
            return std::string();
        }
    /**
     * @return the marker name associated to the \p marker id
     */
    std::string markerName( size_type marker, size_type dim ) const
        {
            for( auto const& n : M_markername )
            {
                if (n.second[0] == marker && n.second[1] == dim )
                    return n.first;
            }
            return std::string{};
        }

    /**
     * @return the topological dimension associated to the \p marker
     */
    size_type markerDim( std::string const& marker ) const
    {
        auto mit = M_markername.find( marker );
        if (  mit != M_markername.end() )
            return mit->second[1];
        return invalid_v<size_type>;
    }

    /**
     * @return the marker names
     */
    std::map<std::string, std::vector<size_type> > markerNames() const
    {
        return M_markername;
    }
    /**
     * set the marker names 
     */
    void setMarkerNames( std::map<std::string, std::vector<size_type> > const& m )
    {
        M_markername = m;
    }
    /**
     * add a new marker name
     */
    void addMarkerName( std::pair<std::string, std::vector<size_type> > const& marker )
        {
            M_markername.insert( marker );
        }

    /**
     * add a new marker name
     */
    void addMarkerName( std::string __name, int __id ,int __topoDim )
        {
            M_markername[__name] = { static_cast<size_type>(__id), static_cast<size_type>(__topoDim) };
        }

    /// @return the marker id given the marker name \p marker
    flag_type markerId( boost::any const& marker ) const;

    /// @return the set of marker id given the marker name \p marker
    std::set<flag_type> markersId( boost::any const& marker ) const;


#if defined(FEELPP_HAS_VTK)
    struct MappingDataWithVTK
    {
        MappingDataWithVTK() = default;
        MappingDataWithVTK( MappingDataWithVTK const& ) = default;
        MappingDataWithVTK( MappingDataWithVTK && ) = default;
        std::unordered_map<size_type,size_type> mapPointsFeelIdToVTKId;
        std::unordered_map<size_type,size_type> mapElementsFeelIdToVTKId;
    };
    typedef std::pair< vtkSmartPointer<vtkUnstructuredGrid>, std::shared_ptr<MappingDataWithVTK> > vtk_export_type;
    //!
    //! exporter to VTK data structure 
    //!
    virtual vtk_export_type exportVTK( bool exportMarkers=false, std::string const& vtkFieldNameMarkers="markers" ) const = 0;
#endif // FEELPP_HAS_VTK

    //@}



protected:

    /**
     * set to the flag whether the mesh is updated for proper use
     */
    void setUpdatedForUse( bool u )
    {
        M_is_updated = u;
    }

    /**
     * After loading/defining a mesh, we want to have as much locality
     * as possible (elements/faces/nodes to be contiguous). In order
     * to do that the mesh elements/faces/nodes are renumbered. That
     * will be then most helpful when generating the \p Dof table.
     * This procedure should work also with
     * \p comm().size() == 1
     *
     */
    virtual void renumber() = 0;

    /**
     * update the entities of co-dimension 1
     */
    virtual void updateEntitiesCoDimensionOne() = 0;

    /**
     * update the entities of co-dimension 2
     */
    virtual void updateEntitiesCoDimensionTwo() = 0;

    /**
     * check mesh connectivity
     */
    virtual void check() const = 0;

    /**
     * check elements orientation and fix it if needed
     */
    //virtual void checkAndFixPermutation() = 0;


private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & M_components;
            ar & M_is_updated;
            ar & M_is_parametric;
            ar & M_n_vertices;
            ar & M_n_parts;
            ar & M_markername;
        }

    void meshesWithNodesSharedImpl( std::vector<std::shared_ptr<MeshBase<IndexT>>> & ret ) const
    {
        for ( std::weak_ptr<MeshBase<IndexT>> m : M_meshesWithNodesShared )
        {
            if ( m.expired() )
                continue;
            auto otherPtr = m.lock();
            auto foundOther = std::find_if(ret.begin(), ret.end(),
                                           [&otherPtr](auto const& ptr) { return otherPtr == ptr; });
            if ( foundOther == ret.end() )
            {
                ret.push_back( otherPtr );
                otherPtr->meshesWithNodesSharedImpl( ret );
            }
        }
    }

    //! remove a mesh which has nodes shared with the current mesh (typically when the mesh is deleted)
    void removeMeshWithNodesShared( MeshBase<IndexT> * m );

protected:

    /**
     * marker name dictionnary ( std::string -> <int,int> )
     * get<0>() provides the id
     * get<1>() provides the topological dimension
     */
    std::map<std::string, std::vector<size_type> > M_markername;

private:

    /**
     * topological dimension
     */
    uint16_type M_topodim;

    /**
     * real dimension
     */
    uint16_type M_realdim;

    /**
     * encodes components that should be updated
     */
    Context M_components;

    /**
     * \p true if mesh ready to be used, \p false otherwise
     */
    bool M_is_updated;

    /**
     * \p true if the mesh is parametric (e.g. has parametric nodes), \p false otherwise
     */
    bool M_is_parametric;

    /**
     * number of vertices
     */
    size_type M_n_vertices;

    /**
     * The number of partitions the mesh has.  This is set by
     * the partitioners, and may not be changed directly by
     * the user.
     * \note The number of partitions *need not* equal
     * M_comm.size(), consider for example the case
     * where you simply want to partition a mesh on one
     * processor and view the result in GMV.
     */
    rank_type M_n_parts;

    // sub mesh data
    smd_ptrtype M_smd;

    //! meshes with nodes shared
    std::vector<std::weak_ptr<MeshBase<IndexT>>> M_meshesWithNodesShared;

    //! mesh supports attached
    mutable std::vector<MeshSupportBase*> M_meshSupportsAttached;

};
#if !defined(FEELPP_INSTANTIATE_NOEXTERN_MESHBASE)
extern template class MeshBase<uint32_type>;
#endif
}
#endif /* __MeshBase_H */
