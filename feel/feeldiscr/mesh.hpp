/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-05

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2012 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file mesh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-05
 */
#ifndef FEELPP_MESH_HPP
#define FEELPP_MESH_HPP 1

#include <boost/version.hpp>
#include <boost/unordered_map.hpp>

#include <boost/foreach.hpp>
#include <boost/signals2/signal.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/base_object.hpp>
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdivision-by-zero"
#endif
#include <boost/archive/text_oarchive.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/mpi/operations.hpp>

#include <feel/feelcore/context.hpp>

#include <feel/feelcore/functors.hpp>
#include <feel/feelmesh/mesh0d.hpp>
#include <feel/feelmesh/mesh1d.hpp>
#include <feel/feelmesh/mesh2d.hpp>
#include <feel/feelmesh/mesh3d.hpp>
#include <feel/feelmesh/meshutil.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelmesh/enums.hpp>

#include <feel/feelpoly/geomap.hpp>
#include <feel/feelalg/boundingbox.hpp>
#include <feel/feelpoly/geomapinv.hpp>




#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/logical/and.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/list/at.hpp>
#include <boost/preprocessor/list/cat.hpp>
#include <boost/preprocessor/list/for_each_product.hpp>
#include <boost/preprocessor/logical/or.hpp>
#include <boost/preprocessor/tuple/to_list.hpp>
#include <boost/preprocessor/tuple/eat.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/facilities/identity.hpp>

#include <boost/enable_shared_from_this.hpp>

namespace Feel
{
const size_type EXTRACTION_KEEP_POINTS_IDS                = ( 1<<0 );
const size_type EXTRACTION_KEEP_EDGES_IDS                 = ( 1<<1 );
const size_type EXTRACTION_KEEP_FACES_IDS                 = ( 1<<2 );
const size_type EXTRACTION_KEEP_VOLUMES_IDS               = ( 1<<3 );
const size_type EXTRACTION_KEEP_ALL_IDS                   = ( EXTRACTION_KEEP_POINTS_IDS |
                                                              EXTRACTION_KEEP_EDGES_IDS |
                                                              EXTRACTION_KEEP_FACES_IDS |
                                                              EXTRACTION_KEEP_VOLUMES_IDS );
const size_type EXTRACTION_KEEP_MESH_RELATION             = ( 1<<4 );
}
#include <feel/feeldiscr/createsubmesh.hpp>

namespace Feel
{
struct PeriodicEntity
{
    PeriodicEntity( int d, int s, int m )
        :
        dim(d),
        slave(s),
        master(m)
        {}
    int dim;
    int slave;
    int master;
    std::map<int,int> correspondingVertices;
};
struct MeshMarkerName
{
	std::string name;
	std::vector<int> ids;

};

std::vector<MeshMarkerName> markerMap( int Dim );



// partitioner class
template<typename Mesh> class Partitioner;

/**
 * @brief unifying mesh class
 * @details This structure is an aggregation of elements, faces,
 * edges(3D) and points data structures and provides a unified
 * interface with respect to the dimension.
 *
 * @tparam GeoShape Geometric entities type
 * @tparam T numerical type for coordinates
 * @tparam Tag = 0 to discriminate between meshes of the same type
 */
template <typename GeoShape, typename T = double, int Tag = 0>
class Mesh
    :
public mpl::if_<mpl::equal_to<mpl::int_<GeoShape::nDim>,mpl::int_<0> >,
    mpl::identity<Mesh0D<GeoShape > >,
    typename mpl::if_<mpl::equal_to<mpl::int_<GeoShape::nDim>,mpl::int_<1> >,
    mpl::identity<Mesh1D<GeoShape > >,
    typename mpl::if_<mpl::equal_to<mpl::int_<GeoShape::nDim>,mpl::int_<2> >,
    mpl::identity<Mesh2D<GeoShape> >,
    mpl::identity<Mesh3D<GeoShape> > >::type>::type>::type::type,
public boost::addable<Mesh<GeoShape,T,Tag> >,
public boost::enable_shared_from_this< Mesh<GeoShape,T,Tag> >
{
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<GeoShape::nDim>,mpl::int_<0> >,
            mpl::identity<Mesh0D<GeoShape> >,
            typename mpl::if_<mpl::equal_to<mpl::int_<GeoShape::nDim>,mpl::int_<1> >,
            mpl::identity<Mesh1D<GeoShape> >,
            typename mpl::if_<mpl::equal_to<mpl::int_<GeoShape::nDim>,mpl::int_<2> >,
            mpl::identity<Mesh2D<GeoShape> >,
            mpl::identity<Mesh3D<GeoShape> > >::type>::type>::type::type super;
public:


    /** @name Constants
     */
    //@{

    static const uint16_type nDim = GeoShape::nDim;
    static const uint16_type nRealDim = GeoShape::nRealDim;
    static const uint16_type Shape = GeoShape::Shape;
    static const uint16_type nOrder = GeoShape::nOrder;
    static const uint16_type tag = Tag;

    //@}
    /** @name Typedefs
     */
    //@{
    typedef Mesh<GeoShape,T,Tag> type;
    typedef boost::shared_ptr<type> ptrtype;
    
    typedef T value_type;
    typedef GeoShape shape_type;
    typedef typename super::return_type return_type;
    typedef typename node<double>::type node_type;

    typedef typename super::super_elements super_elements;
    typedef typename super::elements_type elements_type;
    typedef typename super::element_type element_type;
    typedef typename super::element_iterator element_iterator;
    typedef typename super::element_const_iterator element_const_iterator;

    typedef typename super::super_faces super_faces;
    typedef typename super::faces_type faces_type;
    typedef typename super::face_type face_type;
    typedef typename super::face_iterator face_iterator;
    typedef typename super::face_const_iterator face_const_iterator;

    typedef typename super::edge_type edge_type;

    typedef typename super::points_type points_type;
    typedef typename super::point_type point_type;
    typedef typename super::point_iterator point_iterator;
    typedef typename super::point_const_iterator point_const_iterator;

    typedef typename element_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;

    typedef typename element_type::gm1_type gm1_type;
    typedef boost::shared_ptr<gm1_type> gm1_ptrtype;

    template<size_type ContextID>
    struct gmc
    {
        typedef typename gm_type::template Context<ContextID, element_type> type;
        typedef boost::shared_ptr<type> ptrtype;
    };

    typedef Mesh<shape_type, T, Tag> self_type;
    typedef self_type mesh_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
    typedef self_ptrtype mesh_ptrtype;

    typedef typename element_type::template reference_convex<T>::type reference_convex_type;

    //typedef Partitioner<self_type> partitioner_type;
    //typedef boost::shared_ptr<partitioner_type> partitioner_ptrtype;

    typedef typename super::face_processor_type face_processor_type;
    typedef typename super::face_processor_type element_edge_type;

    typedef typename mpl::if_<mpl::bool_<GeoShape::is_simplex>,
                              mpl::identity< Mesh< Simplex< GeoShape::nDim,1,GeoShape::nRealDim>, value_type, Tag > >,
                              mpl::identity< Mesh< Hypercube<GeoShape::nDim,1,GeoShape::nRealDim>,value_type, Tag > > >::type::type P1_mesh_type;

    typedef boost::shared_ptr<P1_mesh_type> P1_mesh_ptrtype;

    template<int TheTag>
    struct trace_mesh
    {
        typedef typename mpl::if_<mpl::bool_<GeoShape::is_simplex>,
                                  mpl::identity< Mesh< Simplex< GeoShape::nDim-1,nOrder,GeoShape::nRealDim>, value_type, TheTag > >,
                                  mpl::identity< Mesh< Hypercube<GeoShape::nDim-1,nOrder,GeoShape::nRealDim>,value_type, TheTag > > >::type::type type;
        typedef boost::shared_ptr<type> ptrtype;
        typedef boost::shared_ptr<const type> const_ptrtype;
    };
    typedef typename trace_mesh<Tag>::type trace_mesh_type;
    typedef typename trace_mesh<Tag>::ptrtype trace_mesh_ptrtype;


    template<int TheTag>
    struct trace_trace_mesh
    {
        static const uint16_type nDim = (GeoShape::nDim==1)?GeoShape::nDim-1:GeoShape::nDim-2;
        typedef typename mpl::if_<mpl::bool_<GeoShape::is_simplex>,
                                  mpl::identity< Mesh< Simplex<nDim,nOrder,GeoShape::nRealDim>, value_type, TheTag > >,
                                  mpl::identity< Mesh< Hypercube<nDim,nOrder,GeoShape::nRealDim>,value_type, TheTag > > >::type::type type;

        typedef boost::shared_ptr<type> ptrtype;
        typedef boost::shared_ptr<const type> const_ptrtype;
    };
    typedef typename trace_trace_mesh<Tag>::type trace_trace_mesh_type;
    typedef typename trace_trace_mesh<Tag>::ptrtype trace_trace_mesh_ptrtype;


    //@}

    /**
     * Default mesh constructor
     */
    Mesh( WorldComm const& worldComm = Environment::worldComm() );

    ~Mesh()
        {
            VLOG(1) << "Mesh Destructor";
            this->clear();
        }
    void clear()
        {
            VLOG(1) << "Mesh clear()";
            M_gm.reset();
            M_gm1.reset();
            M_tool_localization.reset();
            super::clear();
        }
    /**
     * @brief allocate a new Mesh
     * @return the Mesh shared pointer
     */
    static mesh_ptrtype New()
        {
            return mesh_ptrtype(new mesh_type);
        }

    self_type& operator+=( self_type const& m );

    /** @name Accessors
     */
    //@{



    /**
     * @brief get the global number of elements
     * @details it requires communication in parallel to
     * retrieve and sum the number of elements in each subdomain.
     * @return the global number of elements
     */
    size_type numGlobalElements() const { return M_numGlobalElements; }
    /**
     * @brief get the global number of faces
     * @details it requires communication in parallel to
     * retrieve and sum the number of faces in each subdomain.
     * @return the global number of faces
     */
    size_type numGlobalFaces() const { return M_numGlobalFaces; }

    /**
     * @brief get the global number of edges
     * @details it requires communication in parallel to
     * retrieve and sum the number of edges in each subdomain.
     * @return the global number of edges
     */
    size_type numGlobalEdges() const { return M_numGlobalEdges; }


    /**
     * @brief get the global number of points
     * @details it requires communication in parallel to
     * retrieve and sum the number of points in each subdomain.
     * @return the global number of points
     */
    size_type numGlobalPoints() const { return M_numGlobalPoints; }

    /**
     * @brief get the global number of vertices
     * @details it requires communication in parallel to
     * retrieve and sum the number of vertices in each subdomain.
     * @return the global number of vertices
     */
    size_type numGlobalVertices() const { return M_numGlobalVertices; }

    /**
     * @brief compute the global number of elements,faces,points and vertices
     * @details it requires communications in parallel to
     * retrieve and sum the contribution of each subdomain.
     */
    void updateNumGlobalElements()
    {
        //int ne = numElements();
        int ne = std::distance( this->beginElementWithProcessId( this->worldComm().rank() ),
                                this->endElementWithProcessId( this->worldComm().rank() ) );
        int nf = std::distance( this->beginFaceWithProcessId( this->worldComm().rank() ),
                                this->endFaceWithProcessId( this->worldComm().rank() ) );
        int ned = 0;/*std::distance( this->beginEdgeWithProcessId( this->worldComm().rank() ),
                      this->endEdgeWithProcessId( this->worldComm().rank() ) );*/
        int np = std::distance( this->beginPointWithProcessId( this->worldComm().rank() ),
                                this->endPointWithProcessId( this->worldComm().rank() ) );


        if ( this->worldComm().localSize() >1 )
        {
#if BOOST_VERSION >= 105500
            std::vector<int> globals{ ne, nf, ned, np, (int)this->numVertices() };
            mpi::all_reduce( this->worldComm(), mpi::inplace(globals.data()), 5, std::plus<int>() );
#else
            std::vector<int> locals{ ne, nf, ned, np, (int)this->numVertices() };
            std::vector<int> globals( 5, 0 );
            mpi::all_reduce( this->worldComm(), locals.data(), 5, globals.data(), std::plus<int>() );
#endif
            M_numGlobalElements = globals[0];
            M_numGlobalFaces = globals[1];
            M_numGlobalEdges = globals[2];
            M_numGlobalPoints = globals[3];
            M_numGlobalVertices = globals[4];
        }
        else
        {
            M_numGlobalElements = ne;
            M_numGlobalFaces = nf;
            M_numGlobalEdges = ned;
            M_numGlobalPoints = np;
            M_numGlobalVertices = this->numVertices();
        }
    }
    /**
     * @return the topological dimension
     */
    uint16_type dimension() const
    {
        return nDim;
    }

    /**
     * @return geometric mapping
     */
    gm_ptrtype const& gm() const
    {
        return M_gm;
    }

    /**
         * @return geometric mapping of order 1
         */
    gm1_ptrtype const& gm1() const
    {
        return M_gm1;
    }

    /**
     * @return the geometric mapping
     */
    gm_ptrtype& gm()
    {
        return M_gm;
    }

    /**
        * @return the geometric mapping of order 1
        */
    gm1_ptrtype& gm1()
    {
        return M_gm1;
    }

    /**
     * @return the reference convex associated with the element of the
     * mesh
     */
    reference_convex_type referenceConvex() const
    {
        return reference_convex_type();
    }

    /**
     * @return the face index of the face \p n in the element \p e
     */
    face_processor_type const& localFaceId( element_type const& e,
                                            size_type const n ) const
    {
        return M_e2f.find(std::make_pair(e.id(),n))->second;
    }

    /**
     * @return the face index of the face \p n in the element \p e
     */
    face_processor_type const& localFaceId( size_type const e,
                                            size_type const n ) const
    {
        return M_e2f.find(std::make_pair(e,n))->second;
    }
#if 0
    /**
     * @return the world comm
     */
    WorldComm const& worldComm() const
    {
        return M_worldComm;
    }

    void setWorldComm( WorldComm const& _worldComm )
    {
        M_worldComm = _worldComm;
    }

    mpi::communicator const& comm() const
    {
        return M_worldComm.localComm();
    }
#endif

    /**
     * \return true if the mesh is substructured, false otherwise
     */
    bool subStructuring() const
        {
            return M_substructuring;
        }
    //@}

    /** @name  Mutators
     */
    //@{

    void setSubStructuring( bool s )
        {
            M_substructuring = s;
        }

    /**
     * set the partitioner to \p partitioner
     */
    void setPartitioner( std::string partitioner )
    {
        //M_part = partitioner_ptrtype( partitioner_type::New( partitioner ) );
    }

    /**
     * @return the face index of the face \p n in the element \p e
     */
    face_processor_type& localFaceId( element_type const& e,
                                      size_type const n )
    {
        return M_e2f[std::make_pair(e.id(),n)];
    }

    /**
     * @return the face index of the face \p n in the element \p e
     */
    face_processor_type& localFaceId( size_type const e,
                                      size_type const n )
    {
        return M_e2f[std::make_pair(e,n)];
    }

    /**
     * @return the id associated to the \p marker
     */
    size_type markerName( std::string const& marker ) const
    {
        auto mit = M_markername.find( marker );
        if (  mit != M_markername.end() )
            return mit->second[0];
        return invalid_size_type_value;
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
     * @return the topological dimension associated to the \p marker
     */
    size_type markerDim( std::string const& marker ) const
    {
        auto mit = M_markername.find( marker );
        if (  mit != M_markername.end() )
            return mit->second[1];
        return invalid_size_type_value;
    }

    /**
     * @return the marker names
     */
    std::map<std::string, std::vector<size_type> > markerNames() const
    {
        return M_markername;
    }

    /**
     * @return a localization tool
     */
    struct Localization;
    boost::shared_ptr<typename self_type::Localization> tool_localization()
    {
        return M_tool_localization;
    }

    /**
     * @brief get the average h
     * @return the average h
     */
    value_type hAverage() const { return M_h_avg; }

    /**
     * @brief get the min h
     * @return the min h
     */
    value_type hMin() const { return M_h_min; }

    /**
     * @brief get the max h
     * @return the max h
     */
    value_type hMax() const { return M_h_max; }

    /**
     * @return the measure of the mesh (sum of the measure of the elements)
     */
    value_type measure( bool parallel = true ) const
    {
        if ( parallel )
            return M_meas;
        return M_local_meas;
    }

    /**
     * @return the measure of the mesh (sum of the measure of the elements)
     */
    value_type measureBoundary() const
    {
        return M_measbdy;
    }

    //@}

    /**
     * set the periodic entities
     */
    void setPeriodicEntities( std::vector<PeriodicEntity> const& e ) { M_periodic_entities = e; }

    /**
     * @return true if the mesh has periodic entities
     */
    bool isPeriodic() const { return M_periodic_entities.empty() == false; }

    /** @name  Methods
     */
    //@{

    /**
     * @return true if \p marker exists, false otherwise
     */
    bool
    hasMarker( std::string marker ) const
        {
            return markerName( marker ) != invalid_size_type_value;
        }

    /**
     * @return true if \p marker exists and topological dimension of the entity
     * associated is Dim-1, false otherwise
     */
    bool
    hasFaceMarker( std::string marker ) const
        {
            return ( markerName( marker ) != invalid_size_type_value ) && ( markerDim( marker ) != nDim-1 );
        }

    /**
     * @return true if \p marker exists and topological dimension of the entity
     * associated is Dim-2, false otherwise
     */
    bool
    hasEdgeMarker( std::string marker ) const
        {
            return ( markerName( marker ) != invalid_size_type_value ) && ( markerDim( marker ) != nDim-2 );
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

    /// @return the marker id given the marker name \p marker
    flag_type markerId( boost::any const& marker );

    /**
     * erase element at position \p position
     *
     * @param position \p position is a valid dereferenceable iterator of the index.
     * @param modify if true, update mesh data structure and in particular faces
     *
     * @return An iterator pointing to the element immediately
     * following the one that was deleted, or \c end() if no such element
     * exists.
     */
    element_iterator eraseElement( element_iterator position, bool modify=true );

    /**
     * @brief compute the trace mesh
     * @details compute the trace mesh associated to the
     * range \p range and tag \p TheTag. The \p range provides
     * the iterators over the boundary faces of the mesh.
     * \@code
     * auto trace_mesh1 = mesh->trace(boundaryfaces(mesh));
     * auto trace_mesh2 = mesh->trace(markedfaces(mesh,"marker"));
     * @endcode
     *
     * @param range iterator range
     * @return the computed trace mesh
     */
    template<typename RangeT, int TheTag>
    typename trace_mesh<TheTag>::ptrtype
    trace( RangeT const& range, mpl::int_<TheTag> ) const;

    /**
     * creates a mesh by iterating over the elements between
     * \p begin_elt and and \p end_elt and adding them the the mesh
     * \p mesh
     *
     * \param mesh new mesh to construct
     * \param begin_elt begin iterator
     * \param begin_elt end iterator
     * \param extraction_policies not in use yet
     *
     * \todo make use of \c extraction_policies
     */
    template<int TheTag=Tag>
    typename trace_mesh<TheTag>::ptrtype
    trace() const
    {
        return trace(boundaryfaces(this->shared_from_this()),mpl::int_<TheTag>());
    }

    template<typename RangeT>
    typename trace_mesh<Tag>::ptrtype
    trace( RangeT const& range ) const;


    template<int TheTag=Tag>
    typename trace_mesh<TheTag>::ptrtype
    wireBasket() const
    {
        return wireBasket(boundaryfaces(this->shared_from_this()),mpl::int_<TheTag>());
    }

    template<typename RangeT, int TheTag>
    typename trace_trace_mesh<TheTag>::ptrtype
    wireBasket( RangeT const& range, mpl::int_<TheTag> ) const;

    template<typename RangeT>
    typename trace_trace_mesh<Tag>::ptrtype
    wireBasket( RangeT const& range ) const;


    template<typename Iterator>
    void createSubmesh( self_type& mesh,
                        Iterator const& begin_elt,
                        Iterator const& end_elt,
                        size_type extraction_policies = EXTRACTION_KEEP_ALL_IDS ) const;

    /**
     * A special case of \p createSubmesh() : it creates a mesh by
     * iterating over all elements having the process id \p pid
     *
     * \param mesh new mesh to construct
     * \param pid process id that will select the elements to add
     */
    void createSubmeshByProcessId( self_type& mesh, uint16_type pid ) const
    {
        this->createSubmesh( mesh,
                             this->beginElementWithProcessId( pid ),
                             this->endElementWithProcessId( pid ) );
    }

    /**
     * Create a P1 mesh from the HO mesh
     */
    P1_mesh_ptrtype createP1mesh() const;

    /**
     * update the Marker2 with a range of elements or faces
     * if elements -> update marker2 for this elements
     * if faces -> update marker2 for this faces
     */
    template<typename IteratorRange>
    void updateMarker2WithRange( IteratorRange const& range,flag_type flag )
    {
        const size_type iDim = boost::tuples::template element<0, IteratorRange>::type::value;
        this->updateMarker2WithRange( range,flag,mpl::int_<iDim>() );
    }

    /**
     * sub method of updateMarker2WithRange : MESH_ELEMENTS
     */
    template<typename IteratorRange>
    void updateMarker2WithRange( IteratorRange const& range,flag_type flag, mpl::int_<MESH_ELEMENTS>/**/ )
    {
        const size_type iDim = boost::tuples::template element<0, IteratorRange>::type::value;
        this->updateMarker2WithRangeElements( range,flag );
    }

    /**
     * sub method of updateMarker2WithRange : MESH_FACES
     */
    template<typename IteratorRange>
    void updateMarker2WithRange( IteratorRange const& range,flag_type flag, mpl::int_<MESH_FACES>/**/ )
    {
        this->updateMarker2WithRangeFaces( range,flag );
    }

    /**
     * update the Marker3 with a range of elements or faces
     * if elements -> update marker3 for this elements
     * if faces -> update marker3 for this faces
     */
    template<typename IteratorRange>
    void updateMarker3WithRange( IteratorRange const& range,flag_type flag )
    {
        const size_type iDim = boost::tuples::template element<0, IteratorRange>::type::value;
        this->updateMarker3WithRange( range,flag,mpl::int_<iDim>() );
    }

    /**
     * sub method of updateMarker3WithRange : MESH_ELEMENTS
     */
    template<typename IteratorRange>
    void updateMarker3WithRange( IteratorRange const& range,flag_type flag, mpl::int_<MESH_ELEMENTS>/**/ )
    {
        this->updateMarker3WithRangeElements( range,flag );
    }

    /**
     * sub method of updateMarker3WithRange : MESH_FACES
     */
    template<typename IteratorRange>
    void updateMarker3WithRange( IteratorRange const& range,flag_type flag, mpl::int_<MESH_FACES> /**/ )
    {
        this->updateMarker3WithRangeFaces( range,flag );
    }

    /**
     * Call the default partitioner (currently \p metis_partition()).
     */
    void partition ( const uint16_type n_parts = 1 );

    /**
     * After loading/defining a mesh, we want to have as much locality
     * as possible (elements/faces/nodes to be contiguous). In order
     * to do that the mesh elements/faces/nodes are renumbered. That
     * will be then most helpful when generating the \p Dof table.
     * This procedure should work also with
     * \p comm().size() == 1
     */
    void renumber()
    {
        renumber( mpl::bool_<( nDim > 1 )>() );
    }
    void renumber( std::vector<size_type> const& node_map, mpl::int_<1> );
    void renumber( std::vector<size_type> const& node_map, mpl::int_<2> );
    void renumber( std::vector<size_type> const& node_map, mpl::int_<3> );

    /**
     * This function only take sense in the 3D modal case with a simplex mesh.
     * In order to construct a global C^0 expansion, we need to assure that
     * two contiguous elements share the same top vertex !
     * This can be achieve by a local renumbering of the vertex.
     */

    void localrenumber();

    /**
     * check elements permutation and fix it if needed
     */
    void checkAndFixPermutation();

    /**
     * send the mesh data structure to processor \p p with  \p tag
     */
    void send( int p, int tag );

    /**
     * receive the mesh data structure to processor \p p with  \p tag
     */
    void recv( int p, int tag );

    /**
     * encode the mesh data structure into a tighter data structure and to avoid
     * pointers in order to serialize it for saving/loading and
     * sending/receiving the mesh
     */
    void encode();

    /**
     * decode the mesh data structure from a tighter data structure and to avoid
     * pointers in order to serialize it for saving/loading and
     * sending/receiving the mesh
     */
    void decode();

    BOOST_PARAMETER_MEMBER_FUNCTION( ( void ),
                                     save,
                                     tag,
                                     ( required
                                       ( name,(std::string) )
                                       ( path,* ) )
                                     ( optional
                                       ( type,( std::string ),std::string( "binary" ) )
                                       ( suffix,( std::string ),std::string( "" ) )
                                       ( sep,( std::string ),std::string( "" ) )
                                         ) )
        {
            Feel::detail::ignore_unused_variable_warning( args );

            if ( !fs::exists( fs::path( path ) ) )
            {
                fs::create_directories( fs::path( path ) );
            }

            std::ostringstream os1;
            os1 << name << sep << suffix << "-" << this->worldComm().globalSize() << "." << this->worldComm().globalRank() << ".fdb";
            fs::path p = fs::path( path ) / os1.str();
            fs::ofstream ofs( p );

            if ( type == "binary" )
            {
                boost::archive::binary_oarchive oa( ofs );
                oa << *this;
            }

            else if ( type == "text" )
            {
                boost::archive::text_oarchive oa( ofs );
                oa << *this;
            }

            else if ( type == "xml" )
            {
                //boost::archive::xml_oarchive oa(ofs);
                //oa << *this;
            }
        }
    BOOST_PARAMETER_MEMBER_FUNCTION(
        ( bool ),
        load,
        tag,
        ( required
          ( name,(std::string) )
          ( path,* ) )
        ( optional
          ( update,(size_type), MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES )
          ( type,( std::string ),std::string( "binary" ) )
          ( suffix,( std::string ),std::string( "" ) )
          ( sep,( std::string ),std::string( "" ) )
            )
        )
        {
            Feel::detail::ignore_unused_variable_warning( args );
            std::ostringstream os1;
            os1 << name << sep << suffix << "-" << this->worldComm().globalSize() << "." << this->worldComm().globalRank() << ".fdb";
            fs::path p = fs::path( path ) / os1.str();
            if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                std::cout << "try loading " << p.native()  << "\n";
            if ( !fs::exists( p ) )
            {
                LOG(INFO) << "[mesh::load] failed loading " << p.native() << "\n";
                std::ostringstream os2;
                os2 << name << sep << suffix << "-" << this->worldComm().globalSize() << "." << this->worldComm().globalRank();
                p = fs::path( path ) / os2.str();
                if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                    std::cout << " now try loading " << p.native()  << "\n";

                if ( !fs::exists( p ) )
                {
                    LOG(INFO) << "[mesh::load] failed loading " << p.native() << "\n";
                    return false;
                }
            }

            if ( !fs::is_regular_file( p ) )
            {
                LOG(INFO) << "[mesh::load] failed loading " << p.native() << "\n";
                return false;
            }

            fs::ifstream ifs( p );

            if ( type == "binary" )
            {
                boost::archive::binary_iarchive ia( ifs );
                ia >> *this;
            }

            else if ( type == "text" )
            {
                boost::archive::text_iarchive ia( ifs );
                ia >> *this;
            }

            else if ( type == "xml" )
            {
                //boost::archive::xml_iarchive ia(ifs);
                //ia >> *this;
            }
            if ( update )
            {
                this->components().reset();
                this->components().set( update );
                this->updateForUse();
            }

            else
            {
                this->components().reset();
            }
            return true;
        }


    FEELPP_DEFINE_VISITABLE();
    //@}

    //private:

    /**
     * Finds all the processors that may contain
     * elements that neighbor my elements.  This list
     * is guaranteed to include all processors that border
     * any of my elements, but may include additional ones as
     * well.  This method computes bounding boxes for the
     * elements on each processor and checks for overlaps.
     */
    //void findNeighboringProcessors();

    /**
     * This function checks if the local numbering of the mesh elements
     * is anticlockwise oriented. For the time being, this function
     * only applies to tetrahedra meshes
     */
    void checkLocalPermutation( mpl::bool_<false> ) const {}
    void checkLocalPermutation( mpl::bool_<true> ) const;


    /**
     * update the mesh data structure before using it
     *
     * by default if the number of processors if > 1, we partition
     * the mesh. A different behaviour is controlled by setting
     * properly \p setComponents(), \p components()
     */
    void updateForUse();

    /**
     * update hAverage, hMin, hMax, measure of the mesh and measure of the boundary mesh
     */
    void updateMeasures();

    void meshModified()
        {
            for( auto fit = this->beginFace(), fen = this->endFace(); fit != fen; ++fit )
            {
                if ( fit->isOnBoundary() && !fit->isConnectedTo0() )
                {
                    std::cout << "erase boundary face...\n";
                    this->eraseFace( fit );
                }
                if ( !fit->isOnBoundary() && fit->isConnectedTo0() && !fit->isConnectedTo1() )
                {
                    std::cout << "found boundary face...\n";
                    this->faces().modify( fit, []( face_type& f ){ f.setOnBoundary( true ); } );
                }
            }
            this->setUpdatedForUse( false );
            this->updateForUse();
        }
private:

    void propagateMarkers( mpl::int_<1> ) {}
    void propagateMarkers( mpl::int_<2> );
    void propagateMarkers( mpl::int_<3> );

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            if ( Archive::is_saving::value )
            {
                DVLOG(2) << "Serializing mesh(saving) ...\n";
                DVLOG(2) << "encoding...\n";
                encode();
                DVLOG(2) << "loading markers...\n";
                ar & M_markername;
                DVLOG(2) << "loading pts...\n";
                ar & M_enc_pts;
                DVLOG(2) << "loading faces...\n";
                ar & M_enc_faces;
                DVLOG(2) << "loading elts...\n";
                ar & M_enc_elts;
            }
            if ( Archive::is_loading::value )
            {
                DVLOG(2) << "Serializing mesh(loading) ...\n";
                DVLOG(2) << "loading markers...\n";
                ar & M_markername;
                DVLOG(2) << "loading pts...\n";
                ar & M_enc_pts;
                DVLOG(2) << "loading faces...\n";
                ar & M_enc_faces;
                DVLOG(2) << "loading elts...\n";
                ar & M_enc_elts;
                decode();

            }

        }
public:
    struct Inverse
            :
        public mpl::if_<mpl::bool_<GeoShape::is_simplex>,
            mpl::identity<GeoMapInverse<nDim,nOrder,nRealDim,T,Simplex> >,
            mpl::identity<GeoMapInverse<nDim,nOrder,nRealDim,T,Hypercube> > >::type::type
    {
        typedef typename mpl::if_<mpl::bool_<GeoShape::is_simplex>,
                mpl::identity<GeoMapInverse<nDim,nOrder,nRealDim,T,Simplex> >,
                mpl::identity<GeoMapInverse<nDim,nOrder,nRealDim,T,Hypercube> > >::type::type super;
        typedef typename super::gic_type gic_type;

        typedef typename mpl::if_<mpl::bool_<GeoShape::is_simplex>,
                mpl::identity<GeoMapInverse<nDim,1,nRealDim,T,Simplex> >,
                mpl::identity<GeoMapInverse<nDim,1,nRealDim,T,Hypercube> > >::type::type super1;
        typedef typename super1::gic_type gic1_type;

        Inverse( boost::shared_ptr<self_type> const& m )
            :
            super(),
            M_mesh ( m )
        {}

        ~Inverse()
        {}

        size_type nPointsInConvex( size_type i ) const
        {
            auto itFind = M_pts_cvx.find(i);
            if ( itFind != M_pts_cvx.end() )
                return M_pts_cvx.find(i)->second.size();
            else
                return 0;
        }
        void pointsInConvex( size_type i, std::vector<boost::tuple<size_type, uint16_type > > &itab ) const
        {
            const size_type nPts = this->nPointsInConvex( i );
            itab.resize( nPts );
            if (nPts == 0 ) return;

            auto it = M_pts_cvx.find(i)->second.begin();
            auto const en = M_pts_cvx.find(i)->second.end();
            for ( size_type j = 0 ; it != en ; ++it )
                itab[j++] = boost::make_tuple( it->first, it->second );
        }

        const boost::unordered_map<size_type,node_type> &referenceCoords( void )
        {
            return M_ref_coords;
        }

        /**
         * distribute the points of the mesh in a kdtree
         *
         * extrapolation = false : Only the points inside the mesh are distributed.
         * extrapolation = true  : Try to project the exterior points.
         * \todo : for extrapolation, verify that all the points have been taken
         *        into account, else test them on the frontiere convexes.
         */
        void distribute( bool extrapolation = false );



    private:
        boost::shared_ptr<self_type> M_mesh;
        boost::unordered_map<size_type, boost::unordered_map<size_type,uint16_type > > M_pts_cvx;
        typedef typename std::map<size_type, uint16_type >::const_iterator map_iterator;
        //typedef typename node<value_type>::type node_type;
        boost::unordered_map<size_type,node_type> M_ref_coords;
        boost::unordered_map<size_type,double> M_dist;
        boost::unordered_map<size_type,size_type> M_cvx_pts;

    }; // Inverse

    struct Localization
    {

        typedef Localization localization_type;
        typedef boost::shared_ptr<localization_type> localization_ptrtype;

        typedef typename matrix_node<typename node_type::value_type>::type matrix_node_type;

        typedef boost::weak_ptr<self_type> mesh_ptrtype;

        typedef KDTree kdtree_type;
        typedef typename boost::shared_ptr<KDTree> kdtree_ptrtype;

        // a node x => a list of id elt which contain the node x
        typedef boost::tuple<node_type, std::list<size_type> > node_elem_type;

        typedef typename std::list<boost::tuple<size_type,node_type> > container_output_type;
        typedef typename std::list<boost::tuple<size_type,node_type> >::iterator container_output_iterator_type;
        //map between element id and list of node described in the reference elt
        //typedef std::map<size_type, std::list<boost::tuple<size_type,node_type> > > container_search_type
        typedef std::map<size_type, container_output_type > container_search_type;
        typedef typename container_search_type::const_iterator container_search_const_iterator_type;
        typedef typename container_search_type::iterator container_search_iterator_type;

        // geomap inverse
        typedef typename mpl::if_<mpl::bool_<GeoShape::is_simplex>,
                                  mpl::identity<GeoMapInverse<nDim,nOrder,nRealDim,T,Simplex> >,
                                  mpl::identity<GeoMapInverse<nDim,nOrder,nRealDim,T,Hypercube> > >::type::type gm_inverse_type;
        typedef typename gm_inverse_type::gic_type gmc_inverse_type;

        typedef typename mpl::if_<mpl::bool_<GeoShape::is_simplex>,
                                  mpl::identity<GeoMapInverse<nDim,1,nRealDim,T,Simplex> >,
                                  mpl::identity<GeoMapInverse<nDim,1,nRealDim,T,Hypercube> > >::type::type gm1_inverse_type;
        typedef typename gm1_inverse_type::gic_type gmc1_inverse_type;

        // reference convex
        typedef typename self_type::gm_type::reference_convex_type ref_convex_type;
        typedef typename self_type::gm1_type::reference_convex_type ref_convex1_type;

        /**--------------------------------------------------------------
         * Constructors
         */
        Localization() :
            M_mesh (),
            M_kd_tree( new kdtree_type() ),
            M_isInit( false ),
            M_isInitBoundaryFaces( false ),
            M_doExtrapolation( boption( _name=(boost::format("mesh%1%d.localisation.use-extrapolation") % nDim).str() ) ),
            M_barycenter(),
            M_barycentersWorld()
        {
            DVLOG(2) << "[Mesh::Localization] create Localization tool\n";
            int optNbNeighbor = ioption( _name=(boost::format("mesh%1%d.localisation.nelt-in-leaf-kdtree") % nDim).str() );
            int usedNbNeighbor = ( optNbNeighbor < 0 )? 2*self_type::element_type::numPoints : optNbNeighbor;
            M_kd_tree->nbNearNeighbor( usedNbNeighbor );

            M_resultAnalysis.clear();
            DVLOG(2) << "[Mesh::Localization] create Localization tool done\n";
        }

        Localization( boost::shared_ptr<self_type> m, bool init_b = true ) :
            M_mesh ( m ),
            M_isInit( init_b ),
            M_isInitBoundaryFaces( false ),
            M_doExtrapolation( boption( _name=(boost::format("mesh%1%d.localisation.use-extrapolation") % nDim).str() ) ),
            M_barycenter(),
            M_barycentersWorld()
        {
            if ( this->isInit() )
                this->init();

            int optNbNeighbor = ioption( _name=(boost::format("mesh%1%d.localisation.nelt-in-leaf-kdtree") % nDim).str() );
            int usedNbNeighbor = ( optNbNeighbor < 0 )? 2*self_type::element_type::numPoints : optNbNeighbor;
            M_kd_tree->nbNearNeighbor( usedNbNeighbor );

            M_resultAnalysis.clear();
        }

        Localization( Localization const & L ) :
            M_mesh( L.M_mesh ),
            M_kd_tree( new kdtree_type( *( L.M_kd_tree ) ) ),
            M_geoGlob_Elts( L.M_geoGlob_Elts ),
            M_isInit( L.M_isInit ),
            M_isInitBoundaryFaces( L.M_isInitBoundaryFaces ),
            M_resultAnalysis( L.M_resultAnalysis ),
            M_doExtrapolation( L.M_doExtrapolation ),
            M_gic( L.M_gic ), M_gic1( L.M_gic1 ),
            M_barycenter( L.M_barycenter ),
            M_barycentersWorld( L.M_barycentersWorld )
        {}

        /*--------------------------------------------------------------
         * Define the mesh whith or not init
         */
        void
        setMesh( boost::shared_ptr<self_type> m,bool b=true )
        {
            M_mesh=m;

            if ( b )
                this->init();

            else M_isInit=b;

            M_resultAnalysis.clear();
        }

        /*--------------------------------------------------------------
         * Define if necessary to use extrapolation
         */
        void
        setExtrapolation( bool b )
        {
            M_doExtrapolation = b;
        }

        /*--------------------------------------------------------------
         * Run the init function if necessary
         */
        void updateForUse()
        {
            if ( !this->isInit() )
                this->init();
        }

        /*--------------------------------------------------------------
         * Run the init function if necessary
         */
        void updateForUseBoundaryFaces()
        {
            if ( !this->isInitBoundaryFaces() )
                this->initBoundaryFaces();
        }

        /*--------------------------------------------------------------
         * Access
         */
        bool isInit() const
        {
            return M_isInit;
        }

        bool isInitBoundaryFaces() const
        {
            return M_isInitBoundaryFaces;
        }

        bool doExtrapolation() const
        {
            return M_doExtrapolation;
        }

        mesh_ptrtype mesh()
        {
            return M_mesh;
        }
        mesh_ptrtype mesh() const
        {
            return M_mesh;
        }

        kdtree_ptrtype kdtree()
        {
            return M_kd_tree;
        }

        kdtree_ptrtype const& kdtree() const
        {
            return M_kd_tree;
        }

        node_type const& barycenter() const
        {
            CHECK( this->isInit() || this->isInitBoundaryFaces()  ) << " localization tool not init \n";
            return M_barycenter;
        }

        void computeBarycenter();

        bool hasComputedBarycentersWorld()
        {
#if BOOST_VERSION >= 105600
            return M_barycentersWorld != boost::none;
#else
            return M_barycentersWorld;
#endif
        }

        std::vector<boost::tuple<bool,node_type> > const& barycentersWorld() const
        {
            CHECK( M_barycentersWorld ) << " you must call computeBarycentersWorld() before barycentersWorld() \n";
            return M_barycentersWorld.get();
        }

        void computeBarycentersWorld();

        container_search_type const & result_analysis() const { return M_resultAnalysis;}

        container_search_iterator_type result_analysis_begin()
        {
            return M_resultAnalysis.begin();
        }
        container_search_iterator_type result_analysis_end()
        {
            return M_resultAnalysis.end();
        }

        /*---------------------------------------------------------------
         * True if the node p is in mesh->element(id)
         */
        boost::tuple<bool,node_type,double> isIn( size_type _id, const node_type & _pt ) const;
        boost::tuple<bool,node_type,double> isIn( size_type _id, const node_type & _pt, const matrix_node_type & setPoints, mpl::int_<1> /**/ ) const;
        boost::tuple<uint16_type,std::vector<bool> > isIn( std::vector<size_type> _ids, const node_type & _pt );

        /*---------------------------------------------------------------
         * Research only one element wich contains the node p
         */
        boost::tuple<bool, size_type,node_type> searchElement(const node_type & p);
        /*---------------------------------------------------------------
         * Research only one element wich contains the node p
         */
        boost::tuple<bool, size_type,node_type> searchElement(const node_type & p,
                                                              const matrix_node_type & setPoints,
                                                              mpl::int_<0> /**/ )
        {
            return searchElement( p );
        }

        /*---------------------------------------------------------------
         * Research only one element wich contains the node p and which this elt have as geometric point contain setPoints
         */
        boost::tuple<bool, size_type,node_type> searchElement( const node_type & p,
                                                               const matrix_node_type & setPoints,
                                                               mpl::int_<1> /**/ );

        /*---------------------------------------------------------------
         * Research all elements wich contains the node p
         */
        boost::tuple<bool, std::list<boost::tuple<size_type,node_type> > > searchElements( const node_type & p );

        /*---------------------------------------------------------------
         * Research the element wich contains the node p, forall p in the
         * matrix_node_type m. The result is save by this object
         */
        boost::tuple<std::vector<bool>, size_type> run_analysis(const matrix_node_type & m,
                                                                const size_type & eltHypothetical);

        /*---------------------------------------------------------------
         * Research the element wich contains the node p, forall p in the
         * matrix_node_type m. The result is save by this object
         */
        boost::tuple<std::vector<bool>, size_type>  run_analysis(const matrix_node_type & m,
                                                                 const size_type & eltHypothetical,
                                                                 const matrix_node_type & setPoints,
                                                                 mpl::int_<0> /**/)
        {
            return run_analysis( m,eltHypothetical );
        }

        /*---------------------------------------------------------------
         * Research the element wich contains the node p, forall p in the
         * matrix_node_type m. The result is save by this object
         */
        boost::tuple<std::vector<bool>, size_type> run_analysis(const matrix_node_type & m,
                                                                const size_type & eltHypothetical,
                                                                const matrix_node_type & setPoints,
                                                                mpl::int_<1> /**/);

        /*---------------------------------------------------------------
         * Reset all data
         */
        void reset()
        {
            M_isInit=false;
            M_isInitBoundaryFaces=false;
            this->init();
        }

        void resetBoundaryFaces()
        {
            M_isInit=false;
            M_isInitBoundaryFaces=false;
            this->initBoundaryFaces();
        }

    private :

        /*---------------------------------------------------------------
         *initializes the kd tree and the map between node and list elements(all elements)
         */
        void init();

        /*---------------------------------------------------------------
         *initializes the kd tree and the map between node and list elements(only on boundary)
         */
        void initBoundaryFaces();

        /*---------------------------------------------------------------
         *search near elt in kd tree and get a sorted list
         */
        void searchInKdTree( const node_type & p,
                             std::list< std::pair<size_type, uint> > & listTri );

        /*---------------------------------------------------------------
         * computed barycenter
         */
        node_type computeBarycenter(mpl::int_<1> /**/) const;
        node_type computeBarycenter(mpl::int_<2> /**/) const;
        node_type computeBarycenter(mpl::int_<3> /**/) const;

    private:

        mesh_ptrtype M_mesh;
        kdtree_ptrtype M_kd_tree;
        //map between node and list elements
        std::map<size_type, node_elem_type > M_geoGlob_Elts;
        bool M_isInit,M_isInitBoundaryFaces;
        container_search_type M_resultAnalysis;
        bool M_doExtrapolation;

        ref_convex_type M_refelem;
        ref_convex1_type M_refelem1;
        mutable boost::shared_ptr<gmc_inverse_type> M_gic;
        mutable boost::shared_ptr<gmc1_inverse_type> M_gic1;

        node_type M_barycenter;
        boost::optional<std::vector<boost::tuple<bool,node_type> > > M_barycentersWorld;

    };


    /** @name  Signals
     */
    //@{

#if !defined( __INTEL_COMPILER )
    /**
     * mesh changed its connectivity
     */
    boost::signals2::signal<void ( MESH_CHANGES )> meshChanged;

    template<typename Observer>
    void addObserver( Observer& obs )
    {
				LOG(INFO) << "Observer attached ! \n";
        meshChanged.connect( obs );
    }
#endif // __INTEL_COMPILER

    void removeFacesFromBoundary( std::initializer_list<uint16_type> markers );

    typename std::set<rank_type>::const_iterator beginNeighborSubdomains() const { return M_neighbor_processors.begin(); }
    typename std::set<rank_type>::const_iterator endNeighborSubdomains() const { return M_neighbor_processors.end(); }
    std::set<rank_type> const& neighborSubdomains() const { return M_neighbor_processors; }
    void addNeighborSubdomain( rank_type p ) { M_neighbor_processors.insert( p ); }

    typename std::set<rank_type>::const_iterator beginFaceNeighborSubdomains() const { return M_face_neighbor_processors.begin(); }
    typename std::set<rank_type>::const_iterator endFaceNeighborSubdomains() const { return M_face_neighbor_processors.end(); }
    std::set<rank_type> const& faceNeighborSubdomains() const { return M_face_neighbor_processors; }
    void addFaceNeighborSubdomain( rank_type p ) { M_face_neighbor_processors.insert( p ); }

    //@}

protected:
    /**
     * Update connectivity of entities of codimension 1
     */
    void updateEntitiesCoDimensionOne();
    void updateEntitiesCoDimensionOne(mpl::bool_<true>);
    void updateEntitiesCoDimensionOne(mpl::bool_<false>);

    /**
     * Update in ghost cells of entities of codimension 1
     */
    void updateEntitiesCoDimensionOneGhostCellByUsingBlockingComm();
    void updateEntitiesCoDimensionOneGhostCellByUsingNonBlockingComm();

    /**
     * check mesh connectivity
     */
    void check() const;


private:

    /**
     * \sa renumber()
     */
    void renumber( mpl::bool_<false> ) {}

    /**
     * \sa renumber()
     */
    void renumber( mpl::bool_<true> );

    /**
     * modify edges on boundary in 3D
     */
    void modifyEdgesOnBoundary( face_iterator& face, mpl::bool_<true> );

    /**
     * modify edges on boundary in 2D or 1D
     */
    void modifyEdgesOnBoundary( face_iterator& face, mpl::bool_<false> );

    /**
     * modify element that may touch the boundary through one of its edge in 1D or 2D
     */
    bool modifyElementOnBoundaryFromEdge( element_iterator& e, mpl::bool_<false> );

    /**
     * modify element that may touch the boundary through one of its edge in 3D
     */
    bool modifyElementOnBoundaryFromEdge( element_iterator& e, mpl::bool_<true> );

    /**
     * update entities on boundary (point, edge, face and element)
     */
    void updateOnBoundary();

    /**
     * fix duplication of point in connection1 with 3d mesh at order 3 and 4
     */
    void fixPointDuplicationInHOMesh( element_iterator iv, face_iterator __fit, mpl::true_ );
    void fixPointDuplicationInHOMesh( element_iterator iv, face_iterator __fit, mpl::false_ );

private:

    //! communicator
    size_type M_numGlobalElements, M_numGlobalFaces, M_numGlobalEdges, M_numGlobalPoints, M_numGlobalVertices;

    bool M_is_gm_cached = false;
    gm_ptrtype M_gm;
    gm1_ptrtype M_gm1;

    value_type M_h_avg;
    value_type M_h_min;
    value_type M_h_max;


    //! measure of the mesh
    value_type M_meas, M_local_meas;

    //! measure of the boundary of the mesh
    value_type M_measbdy, M_local_measbdy;

    //!sub structuring
    bool M_substructuring;

    /**
     * The processors who neighbor the current
     * processor
     */
    //std::vector<uint16_type> M_neighboring_processors;
    std::set<rank_type> M_neighbor_processors;
    std::set<rank_type> M_face_neighbor_processors;

    //partitioner_ptrtype M_part;

    /**
     * Arrays containing the global ids of Faces of each element
     */
    boost::unordered_map<std::pair<int,int>,face_processor_type> M_e2f;

    /**
     * Arrays containing the global ids of edges of each element
     */
    boost::multi_array<element_edge_type,2> M_e2e;

    /**
     * marker name disctionnary ( std::string -> <int,int> )
     * get<0>() provides the id
     * get<1>() provides the topological dimension
     */
    std::map<std::string, std::vector<size_type> > M_markername;

    /**
     * periodic entities
     */
    std::vector<PeriodicEntity> M_periodic_entities;

    /**
     * to encode points coordinates
     */
    std::map<uint64_type, boost::tuple<bool, std::vector<int>, std::vector<double> > > M_enc_pts;

    /**
     * to encode elements
     */
    std::map<uint64_type, std::vector<int> > M_enc_faces;

    /**
     * to encode elements
     */
    std::map<uint64_type, std::vector<int> > M_enc_elts;

    /**
     * tool for localize point in the mesh
     */
    boost::shared_ptr<Localization> M_tool_localization;

};

template<typename Shape, typename T, int Tag>
const uint16_type Mesh<Shape, T, Tag>::nDim;
template<typename Shape, typename T, int Tag>
const uint16_type Mesh<Shape, T, Tag>::nOrder;

template<typename Shape, typename T, int Tag>
template<typename RangeT>
typename Mesh<Shape, T, Tag>::template trace_mesh<Tag>::ptrtype
Mesh<Shape, T, Tag>::trace( RangeT const& range ) const
{
    DVLOG(2) << "[trace] extracting " << range.template get<0>() << " nb elements :"
                  << std::distance(range.template get<1>(),range.template get<2>()) << "\n";
    return Feel::createSubmesh<const mesh_type,RangeT,Tag>( this->shared_from_this(), range );

}

template<typename Shape, typename T, int Tag>
template<typename RangeT>
typename Mesh<Shape, T, Tag>::template trace_trace_mesh<Tag>::ptrtype
Mesh<Shape, T, Tag>::wireBasket( RangeT const& range ) const
{
    DVLOG(2) << "[trace] extracting " << range.template get<0>() << " nb elements :"
                  << std::distance(range.template get<1>(),range.template get<2>()) << "\n";
    return Feel::createSubmesh<const mesh_type,RangeT,Tag>( this->shared_from_this(), range );

}


template<int TheTag> struct Tag : public mpl::int_<TheTag>  {};

template<typename Shape, typename T, int Tag>
template<typename RangeT,int TheTag>
typename Mesh<Shape, T, Tag>::template trace_mesh<TheTag>::ptrtype
Mesh<Shape, T, Tag>::trace( RangeT const& range, mpl::int_<TheTag> ) const
{
    DVLOG(2) << "[trace] extracting " << range.template get<0>() << " nb elements :"
                  << std::distance(range.template get<1>(),range.template get<2>()) << "\n";
    return Feel::createSubmesh<const mesh_type,RangeT,TheTag>( this->shared_from_this(), range );

}

template<typename Shape, typename T, int Tag>
template<typename RangeT,int TheTag>
typename Mesh<Shape, T, Tag>::template trace_trace_mesh<TheTag>::ptrtype
Mesh<Shape, T, Tag>::wireBasket( RangeT const& range, mpl::int_<TheTag> ) const
{
    DVLOG(2) << "[trace] extracting " << range.template get<0>() << " nb elements :"
                  << std::distance(range.template get<1>(),range.template get<2>()) << "\n";
    return Feel::createSubmesh<const mesh_type,RangeT,TheTag>( this->shared_from_this(), range );

}

template<typename Shape, typename T, int Tag>
template<typename Iterator>
void
Mesh<Shape, T, Tag>::createSubmesh( self_type& new_mesh,
                               Iterator const& begin_elt,
                               Iterator const& end_elt,
                               size_type extraction_policies ) const
{
#if 0
    new_mesh = Feel::createSubmesh( this->shared_from_this(), boost::make_tuple(mpl::int_<MESH_ELEMENTS>(),begin_elt, end_elt), extraction_policies );
#else
    Context policies( extraction_policies );

    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] start\n";
    // make sure it is all empty
    new_mesh.clear();
    // inherit all the markers
    new_mesh.M_markername = this->markerNames();
    BOOST_FOREACH( auto marker, new_mesh.M_markername )
    {
        LOG(INFO) << "marker name " << marker.first
                  << " id: " << marker.second[0]
                  << " geoe: " << marker.second[1] << "\n";

    }
    // How the nodes on this mesh will be renumbered to nodes
    // on the new_mesh.
    std::vector<size_type> new_node_numbers ( this->numPoints(), invalid_size_type_value );
    std::vector<size_type> new_vertex ( this->numPoints(), 0 );

    // the number of nodes on the new mesh, will be incremented
    unsigned int n_new_nodes = 0;
    unsigned int n_new_elem  = 0;
    size_type n_new_faces = 0;

    for ( Iterator it = begin_elt; it != end_elt; ++it )
    {

        element_type const& old_elem = *it;

        // copy element so that we can modify it
        element_type new_elem = old_elem;

        // Loop over the nodes on this element.
        for ( unsigned int n=0; n < old_elem.nPoints(); n++ )
        {
            FEELPP_ASSERT ( old_elem.point( n ).id() < new_node_numbers.size() ).error( "invalid point id()" );

            if ( new_node_numbers[old_elem.point( n ).id()] == invalid_size_type_value )
            {
                new_node_numbers[old_elem.point( n ).id()] = n_new_nodes;

                DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] insert point " << old_elem.point( n ) << "\n";

                point_type pt( old_elem.point( n ) );
                pt.setId( n_new_nodes );

                // Add this node to the new mesh
                new_mesh.addPoint ( pt );

                DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] number of  points " << new_mesh.numPoints() << "\n";

                // Increment the new node counter
                n_new_nodes++;

                if ( n < element_type::numVertices )
                {
                    FEELPP_ASSERT( new_vertex[old_elem.point( n ).id()] == 0 ).error( "already seen this point?" );
                    new_vertex[old_elem.point( n ).id()]=1;
                }
            }

            // Define this element's connectivity on the new mesh
            FEELPP_ASSERT ( new_node_numbers[old_elem.point( n ).id()] < new_mesh.numPoints() ).error( "invalid connectivity" );

            DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] adding point old(" << old_elem.point( n ).id()
                          << ") as point new(" << new_node_numbers[old_elem.point( n ).id()]
                          << ") in element " << new_elem.id() << "\n";

            new_elem.setPoint( n, new_mesh.point( new_node_numbers[old_elem.point( n ).id()] ) );

        }

        // set id of element
        new_elem.setId ( n_new_elem );

        // increment the new element counter
        n_new_elem++;


        // Add an equivalent element type to the new_mesh
        new_mesh.addElement( new_elem );


        // Maybe add faces for this element
        for ( unsigned int s=0; s<old_elem.numTopologicalFaces; s++ )
        {
            if ( !old_elem.facePtr( s ) ) continue;

#if 0
            std::cout << "local face id: " << s
                      << " global face id: " << old_elem.face( s ).id() << "\n";
#endif
            // only add face on the boundary: they have some data
            // (boundary ids) which cannot be retrieved otherwise
            //if ( old_elem.neighbor(s) == invalid_size_type_value )
            size_type global_face_id = old_elem.face( s ).id();

            if ( this->hasFace( global_face_id ) )
            {

                //std::cout << "found face " << global_face_id << "\n";
                // get the corresponding face
                face_type const& old_face = old_elem.face( s );
                face_type new_face = old_face;

                // disconnect from elements of old mesh,
                // the connection will be redone in
                // \c updateForUse()
                new_face.disconnect();

                //std::cout << "disconnect face\n";
                // update points info
                for ( uint16_type p = 0; p < new_face.nPoints(); ++p )
                {
#if 0
                    std::cout << "add new point " << new_face.point( p ).id() << " to face \n";
                    std::cout << "add old point " << old_face.point( p ).id() << " to face \n";
                    std::cout << "new point id " << new_node_numbers[old_elem.point( old_elem.fToP( s,p ) ).id()] << "\n";
#endif
                    new_face.setPoint( p, new_mesh.point( new_node_numbers[old_elem.point( old_elem.fToP( s,p ) ).id()] ) );

                }

                new_face.setId( n_new_faces++ );
#if 0
                std::cout << "face id" << new_face.id()
                          << " marker1 : " << new_face.marker()
                          << " old marker1 : " << old_face.marker()
                          << "\n";
#endif
                // add it to the list of faces
                new_mesh.addFace( new_face );
            }
        }
    }

    new_mesh.setNumVertices( std::accumulate( new_vertex.begin(), new_vertex.end(), 0 ) );

    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] update face/edge info if necessary\n";
    // Prepare the new_mesh for use
    new_mesh.components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    new_mesh.updateForUse();

    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] stop\n";
#endif
}


template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::P1_mesh_ptrtype
Mesh<Shape, T, Tag>::createP1mesh() const
{

    P1_mesh_ptrtype new_mesh( new P1_mesh_type( this->worldComm() ) );

    // How the nodes on this mesh will be renumbered to nodes on the new_mesh.
    boost::unordered_map<size_type,size_type> new_node_numbers;
    boost::unordered_map<size_type,int> new_vertex;
    std::vector<size_type> new_element_numbers ( this->numElements(), invalid_size_type_value );

    const int nProc = new_mesh->worldComm().localSize();

    // the number of nodes on the new mesh, will be incremented
    unsigned int n_new_nodes = 0;
    unsigned int n_new_elem  = 0;
    size_type n_new_faces = 0;

    // inherit the table of markersName
    BOOST_FOREACH( auto itMark, this->markerNames() )
    {
        new_mesh->addMarkerName( itMark.first,itMark.second[0],itMark.second[1] );
    }

    // data usefull for parallism
    std::map< int, std::set<boost::tuple<size_type,size_type> > > memoryGhostId;
    //std::set< int > setOfRecvProc;
    std::vector<int> nbMsgToRecv( nProc , 0 );

    auto it = this->beginElement();
    auto const en = this->endElement();
    for ( ; it != en; ++it )
    {
        element_type const& old_elem = *it;

        // create a new element
        typename P1_mesh_type::element_type new_elem;
        // set id of element
        new_elem.setId ( n_new_elem );
        new_element_numbers[old_elem.id()]= n_new_elem;
        // set element markers
        new_elem.setMarker( old_elem.marker().value() );
        new_elem.setMarker2( old_elem.marker2().value() );
        new_elem.setMarker3( old_elem.marker3().value() );
        // partitioning update
        new_elem.setProcessIdInPartition( old_elem.pidInPartition() );
        new_elem.setNumberOfPartitions(old_elem.numberOfPartitions());
        new_elem.setProcessId(old_elem.processId());
        new_elem.setNeighborPartitionIds(old_elem.neighborPartitionIds());

        // Loop over the P1 nodes on this element.
        for ( uint16_type n=0; n < element_type::numVertices; n++ )
            {
                auto const& old_point = old_elem.point( n );

                //if ( !new_node_numbers[old_point.id()] )
                if ( new_node_numbers.find( old_point.id() ) == new_node_numbers.end() )
                    {
                        new_node_numbers[old_point.id()] = n_new_nodes;
                        DVLOG(2) << "[Mesh<Shape,T>::createP1mesh] insert point " << old_point << "\n";
                        typename P1_mesh_type::point_type pt( old_point );
                        pt.setId( n_new_nodes );
                        pt.setProcessId(old_point.processId());
                        pt.clearElementsGhost();

                        // Add this node to the new mesh
                        new_mesh->addPoint( pt );
                        DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] number of  points " << new_mesh->numPoints() << "\n";
                        // Increment the new node counter
                        n_new_nodes++;
                        FEELPP_ASSERT( !new_vertex[old_point.id()] ).error( "already seen this point?" );
                        new_vertex[old_point.id()]=1;
                    }
            // Define this element's connectivity on the new mesh
                //FEELPP_ASSERT ( new_node_numbers[old_elem.point( n ).id()] < new_mesh->numPoints() ).error( "invalid connectivity" );
            DVLOG(2) << "[Mesh<Shape,T>::createP1mesh] adding point old(" << old_point.id()
                          << ") as point new(" << new_node_numbers[old_point.id()]
                          << ") in element " << new_elem.id() << "\n";
            // add point in element
            new_elem.setPoint( n, new_mesh->point( new_node_numbers[old_point.id()] ) );
        } //for ( uint16_type n=0; n < element_type::numVertices; n++ )

        // Add an equivalent element type to the new_mesh
        new_mesh->addElement( new_elem );
        // increment the new element counter
        n_new_elem++;

        // Maybe add faces for this element
        for ( unsigned int s=0; s<old_elem.numTopologicalFaces; s++ )
        {
            if ( !old_elem.facePtr( s ) ) continue;
            // only add face on the boundary: they have some data
            // (boundary ids) which cannot be retrieved otherwise
            const size_type global_face_id = old_elem.face( s ).id();
            if ( this->hasFace( global_face_id ) )
            {
                // get the corresponding face
                face_type const& old_face = old_elem.face( s );
                typename P1_mesh_type::face_type new_face;
                // disconnect from elements of old mesh,
                // the connection will be redone in updateForUse()
                new_face.disconnect();
                // is on boundary
                new_face.setOnBoundary( old_face.isOnBoundary() );
                // set id of face
                new_face.setId( n_new_faces );
                // set face markers
                new_face.setMarker( old_face.marker().value() );
                new_face.setMarker2( old_face.marker2().value() );
                new_face.setMarker2( old_face.marker3().value() );
                // partitioning update
                new_face.setProcessIdInPartition( old_face.pidInPartition() );
                new_face.setNumberOfPartitions(old_face.numberOfPartitions());
                new_face.setProcessId(old_face.processId());
                new_face.clearIdInOthersPartitions();
                new_face.setNeighborPartitionIds(old_face.neighborPartitionIds());
                // update P1 points info
                for ( uint16_type p = 0; p < face_type::numVertices; ++p )
                {
                    //new_face.setPoint( p, new_mesh->point( new_node_numbers[old_elem.point( old_elem.fToP( s,p ) ).id()] ) );
                    new_face.setPoint( p, new_mesh->point( new_node_numbers[ old_face.point(p).id()] ) );
                }
                // add it to the list of faces
                new_mesh->addFace( new_face );
                // increment the new face counter
                ++n_new_faces;
            } // if ( this->hasFace( global_face_id ) )
        } // for ( unsigned int s=0; s<old_elem.numTopologicalFaces; s++ )


        if ( it->isGhostCell() )
        {
            DVLOG(2) << "element " << it->id() << " is a ghost cell\n";
            for (auto it_pid=it->idInOthersPartitions().begin(),en_pid=it->idInOthersPartitions().end() ; it_pid!=en_pid ; ++it_pid)
            {
                DVLOG(2) << " " << it_pid->first << "-" << it_pid->second << "-"<<it->pidInPartition()<<"-"<<new_mesh->worldComm().localRank();
                const int procToSend=it_pid->first;
                DCHECK( procToSend!=it->pidInPartition() ) << "invalid\n";
                memoryGhostId[procToSend].insert(boost::make_tuple(new_elem.id(),it_pid->second));
            }
        }
        else if (it->numberOfNeighborPartitions() /*it->numberOfPartitions()*/ >0 )
        {
#if 0
            setOfRecvProc.insert( it->neighborPartitionIds().begin(), it->neighborPartitionIds().end() );
#else
            auto itneighbor = it->neighborPartitionIds().begin();
            auto const enneighbor = it->neighborPartitionIds().end();
            for ( ; itneighbor!=enneighbor ; ++itneighbor )
                nbMsgToRecv[*itneighbor]++;
#endif
        }
    } // end for it

#if 0
    if ( nProc > 1 )
    {
        std::map< int, std::vector<size_type> > memoryMpiMsg;

        auto itghostproc = memoryGhostId.begin();
        auto const enghostproc = memoryGhostId.end();
        for ( ; itghostproc!=enghostproc ; ++itghostproc )
        {
            const int procToSend = itghostproc->first;
            const int sizeMsgToSend = itghostproc->second.size();
            std::vector<size_type> dataToSend( sizeMsgToSend );
            memoryMpiMsg[procToSend].resize( sizeMsgToSend );

            auto itghostelt = itghostproc->second.begin();
            auto const enghostelt = itghostproc->second.end();
            for ( int k=0 ; itghostelt!=enghostelt ; ++itghostelt,++k )
            {
                dataToSend[k] = itghostelt->template get<1>();
                memoryMpiMsg[procToSend][k] = itghostelt->template get<0>();
            }
            new_mesh->worldComm().localComm().send(procToSend, 0, dataToSend);
        }

        auto itrecvproc = setOfRecvProc.begin();
        auto const enrecvproc = setOfRecvProc.end();
        for ( ; itrecvproc!=enrecvproc ; ++itrecvproc )
        {
            const int procToRecv = *itrecvproc;
            std::vector<size_type> dataToRecv;
            new_mesh->worldComm().localComm().recv(procToRecv, 0, dataToRecv);

            const int nbDataToTreat = dataToRecv.size();
            std::vector<size_type> dataToSend( nbDataToTreat );
            for ( int k=0;k<nbDataToTreat;++k )
            {
                CHECK( dataToRecv[k]!=invalid_size_type_value ) << "invalid id recv \n";
                dataToSend[k] = new_element_numbers[ dataToRecv[k] ];
            }
            new_mesh->worldComm().localComm().send(procToRecv, 1, dataToSend);
        }

        itghostproc = memoryGhostId.begin();
        for ( ; itghostproc!=enghostproc ; ++itghostproc )
        {
            const int procToRecv = itghostproc->first;
            std::vector<size_type> dataToRecv;
            new_mesh->worldComm().localComm().recv(procToRecv, 1, dataToRecv);
            const int nbDataToTreat = dataToRecv.size();
            for ( int k=0;k<nbDataToTreat;++k )
            {
                auto eltToUpdate = new_mesh->elementIterator( memoryMpiMsg[procToRecv][k]/*e.id()*/,  procToRecv );
                new_mesh->elements().modify( eltToUpdate, Feel::detail::updateIdInOthersPartitions( procToRecv, dataToRecv[k]/*idEltAsked*/ ) );
            }
        }

    } // if ( nProc > 1 )
#else
    if ( nProc > 1 )
    {
        std::map< int, std::vector<size_type> > memoryMpiMsg;
        std::vector<int> nbMsgToSend( nProc , 0 );

        auto itghostproc = memoryGhostId.begin();
        auto const enghostproc = memoryGhostId.end();
        for ( ; itghostproc!=enghostproc ; ++itghostproc )
        {
            const int procToSend = itghostproc->first;
            const int sizeMsgToSend = itghostproc->second.size();
            memoryMpiMsg[procToSend].resize( sizeMsgToSend );

            auto itghostelt = itghostproc->second.begin();
            auto const enghostelt = itghostproc->second.end();
            for ( int k=0 ; itghostelt!=enghostelt ; ++itghostelt,++k )
            {
                const size_type idInOtherPart = itghostelt->template get<1>();
                new_mesh->worldComm().localComm().send(procToSend, nbMsgToSend[procToSend], idInOtherPart);
                ++nbMsgToSend[procToSend];
                memoryMpiMsg[procToSend][k] = itghostelt->template get<0>();
            }
            //CHECK( nbMsgToSend[procToSend] == sizeMsgToSend ) << "invalid data to send\n";
        }

#if !defined( NDEBUG )
        // check nbMsgToRecv computation
        std::vector<int> nbMsgToRecv2( nProc , 0 );
        mpi::all_to_all( new_mesh->worldComm().localComm(),
                         nbMsgToSend,
                         nbMsgToRecv2 );
        for ( int proc=0; proc<nProc; ++proc )
            CHECK( nbMsgToRecv[proc]==nbMsgToRecv2[proc] ) << "paritioning data incorect "
                                                           << "myrank " << this->worldComm().localRank() << " proc " << proc
                                                           << " nbMsgToRecv[proc] " << nbMsgToRecv[proc]
                                                           << " nbMsgToRecv2[proc] " << nbMsgToRecv2[proc] << "\n";
#endif

        // recv dof asked and re-send dof in this proc
        for ( int procToRecv=0; procToRecv<nProc; ++procToRecv )
        {
            for ( int cpt=0; cpt<nbMsgToRecv[procToRecv]; ++cpt )
            {
                //recv
                size_type idEltRecv;
                new_mesh->worldComm().localComm().recv( procToRecv, cpt, idEltRecv );

                const size_type idEltAsked = new_element_numbers[ idEltRecv ];
                DCHECK( idEltAsked!=invalid_size_type_value ) << "invalid elt id\n";

                new_mesh->worldComm().localComm().send(procToRecv, cpt, idEltAsked);
            }
        }

        itghostproc = memoryGhostId.begin();
        for ( ; itghostproc!=enghostproc ; ++itghostproc )
        {
            const int procToRecv = itghostproc->first;

            auto itghostelt = itghostproc->second.begin();
            auto const enghostelt = itghostproc->second.end();
            for ( int k=0 ; itghostelt!=enghostelt ; ++itghostelt,++k )
            {
                size_type idEltAsked;
                new_mesh->worldComm().localComm().recv(procToRecv, k, idEltAsked);

                auto eltToUpdate = new_mesh->elementIterator( memoryMpiMsg[procToRecv][k]/*e.id()*/,  procToRecv );
                new_mesh->elements().modify( eltToUpdate, Feel::detail::updateIdInOthersPartitions( procToRecv, idEltAsked ) );
            }
        }



    }  // if ( nProc > 1 )
#endif

    new_mesh->setNumVertices( std::accumulate( new_vertex.begin(), new_vertex.end(), 0,
                                               []( size_type lhs, std::pair<size_type,int> const& rhs )
                                               {
                                                   return lhs+rhs.second;
                                               } ) );


    DVLOG(2) << "[Mesh<Shape,T>::createP1mesh] update face/edge info if necessary\n";
    // Prepare the new_mesh for use
    new_mesh->components().set ( MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    // run intensive job
    new_mesh->updateForUse();

    return new_mesh;
}
namespace detail
{
template<typename T>
struct MeshPoints
{
    template<typename MeshType, typename IteratorType>
    MeshPoints( MeshType* mesh, const WorldComm &, IteratorType it, IteratorType en, const bool outer = false, const bool renumber = false, const bool fill = false, const int startIndex = 1 );

    int translatePointIds(std::vector<int32_t> & ids);
    int translateElementIds(std::vector<int32_t> & ids);

    int globalNumberOfPoints() const { return global_npts; }
    int globalNumberOfElements() const { return global_nelts; }

    std::vector<int> numberOfPoints, numberOfElements;
    int global_nelts{0}, global_npts{0};
    std::vector<int32_t> ids;
    std::map<int32_t, int32_t> new2old;
    std::map<int32_t, int32_t> old2new;
    std::map<int32_t, int32_t> nodemap;
    std::vector<T> coords;
    std::vector<int32_t> elemids;
    std::vector<int32_t> elem;
    size_type offsets_pts, global_offsets_pts;
    size_type offsets_elts, global_offsets_elts;

};

/**
 * Builds information around faces/elements for exporting data
 * @param mesh The mesh from which data is extracted
 * @param it Starting iterator over the faces/elements
 * @param en Ending iterator over the faces/elements
 * @param outer If false, the vertices are place in an x1 y1 z1 ... xn yn zn order, otherwise in the x1 ... xn y1 ... yn z1 ... zn
 * @param renumber If true, the vertices will be renumbered with maps to keep the correspondance between the twoi, otherwise the original ids are kept
 * @param fill It true, the method will generate points coordinates that are 3D, even if the point is specified with 1D or 2D coordinates (filled with 0)
 * @param Specify the startIndex of the renumbered points (typically set to 0 or 1, but no restriction). This is only used when renumber is true, otherwise it is not used.
 */
template<typename T>
template<typename MeshType, typename IteratorType>
MeshPoints<T>::MeshPoints( MeshType* mesh, const WorldComm& worldComm, IteratorType it, IteratorType en, const bool outer, const bool renumber, const bool fill, const int startIndex )
{
    std::set<int> nodeset;
    size_type p = 0;
    auto elt_it = it;

    /* Gather all the vertices of which the elements are made up with into a std::set */
    /* build up correspondance arrays between index in nodeset and previous id */
    for( auto eit = it ; eit != en; ++eit )
    {
        auto const& elt = boost::unwrap_ref( *eit );
        for ( size_type j = 0; j < elt.numLocalVertices; j++ )
        {
            int pid = elt.point( j ).id();
            auto ins = nodeset.insert( pid );
            if ( ins.second )
            {
                if ( renumber )
                { ids.push_back( p + startIndex ); }
                else
                { ids.push_back( pid ); }
                old2new[pid]=ids[p];
                new2old[ids[p]]=pid;
                nodemap[pid] = p;
                ++p;
            }
        }
    }
    CHECK( p == ids.size() ) << "Invalid number of points " << ids.size() << "!=" << p;
    int nv = ids.size();

    coords.resize( 3*nv, 0 );

    auto pit = ids.begin();
    auto pen = ids.end();
    //for( auto i = 0; i < nv; ++i )

    /* put coords of each point into the coords array */
    /* if outer is true, the coords are placed like: x1 x2 ... xn y1 y2 ... yn z1 z2 ... zn */
    /* otherwise, the coords are placed like: x1 y1 z1 x2 y2 z2 ... xn yn zn */
    for( int i = 0; pit != pen; ++pit, ++i )
    {
        //CHECK( *pit > 0 ) << "invalid id " << *pit;
        //LOG(INFO) << "p " << i << "/" << nv << " =" << *pit;
        //int pid = (renumber)?nodemap[*pit]+1:*pit;
        int pid = *pit;

        auto const& p = mesh->point( new2old[*pit] );
        if ( outer )
        { coords[i] = ( T ) p.node()[0]; }
        else
        { coords[3*i] = ( T ) p.node()[0]; }

        if ( MeshType::nRealDim >= 2 )
        {
            if ( outer )
            { coords[nv+i] = ( T ) p.node()[1]; }
            else
            { coords[3*i+1] = ( T ) p.node()[1]; }
        }
        /* Fill 2nd components with 0 if told to do so */
        else
        {
            if(fill)
            {
                if ( outer )
                { coords[nv+i] = (T)0; }
                else
                { coords[3*i+1] = (T)0; }
            }
        }

        if ( MeshType::nRealDim >= 3 )
        {
            if ( outer )
            { coords[2*nv+i] = T( p.node()[2] ); }
            else
            { coords[3*i+2] = T( p.node()[2] ); }
        }
        /* Fill 3nd components with 0 if told to do so */
        else
        {
            if(fill)
            {
                if ( outer )
                { coords[2*nv+i] = (T)0; }
                else
                { coords[3*i+2] = (T)0; }
            }
        }
    }

    /* number of local elements */
    int __ne = std::distance( it, en );

    /* only do this resize if we have at least one elements in the iterator */
    /* otherwise it will segfault */
    if(it != en)
    {
        elem.resize( __ne*boost::unwrap_ref( *it ).numLocalVertices );
        //elem.resize( __ne*mesh->numLocalVertices() );
        elemids.resize( __ne );
    }

    /* build the array containing the id of each vertex for each element */
    elt_it = it;
    size_type e=0;
    for (  ; elt_it != en; ++elt_it, ++e )
    {
        auto const& elt = boost::unwrap_ref( *elt_it );
        elemids[e] = elt.id()+1;
        //std::cout << "LocalV = " << elt.numLocalVertices << std::endl;
        //for ( size_type j = 0; j < mesh->numLocalVertices(); j++ )
        for ( size_type j = 0; j < elt.numLocalVertices; j++ )
        {
            //std::cout << "LocalVId = " << j << " " << e*elt.numLocalVertices+j << std::endl;
            //std::cout << elt.point( j ).id() << std::endl;
            // ensight id start at 1
            elem[e*elt.numLocalVertices+j] = old2new[elt.point( j ).id()];
#if 0
            DCHECK( (elem[e*mesh->numLocalVertices()+j] > 0) && (elem[e*mesh->numLocalVertices()+j] <= nv ) )
                << "Invalid entry : " << elem[e*mesh->numLocalVertices()+j]
                << " at index : " << e*mesh->numLocalVertices()+j
                << " element :  " << e
                << " vertex :  " << j;
#endif
        }
    }
#if 0
    CHECK( e==__ne) << "Invalid number of elements, e= " << e << "  should be " << __ne;
    std::for_each( elem.begin(), elem.end(), [=]( int e )
                   { CHECK( ( e > 0) && e <= __nv ) << "invalid entry e = " << e << " nv = " << nv; } );
#endif

    //size_type offset_pts = ids.size()*sizeof(int)+ coords.size()*sizeof(float);
    //size_type offset_elts = elemids.size()*sizeof(int)+ elem.size()*sizeof(int);
#if 0
    size_type offset_pts = coords.size()*sizeof(float);
    size_type offset_elts = elem.size()*sizeof(int);
#else
    size_type offset_pts = nv;
    size_type offset_elts = __ne;
#endif

    /* gather the number of points and elements fo each process */
    std::vector<int> ost{nv, __ne};
    std::vector<std::vector<int>> ospe;

    mpi::all_gather( worldComm.comm(), ost, ospe );

    /* copy information about number of points/elements
     * per process in a local array */
    for( size_type i = 0; i < ospe.size(); i++)
    {
        numberOfPoints.push_back(ospe[i][0]);
        numberOfElements.push_back(ospe[i][1]);
    }

    /* compute offsets to shift the point and element ids */
    /* regarding to the processor rank */
    offsets_pts = 0;
    global_offsets_pts = 0;
    offsets_elts = 0;
    global_offsets_elts = 0;
    for( size_type i = 0; i < ospe.size(); i++ )
    {
        if ( i < worldComm.localRank() )
        {
            offsets_pts += ospe[i][0];
            offsets_elts += ospe[i][1];
        }
        global_offsets_pts += ospe[i][0];
        global_offsets_elts += ospe[i][1];
    }
    global_npts = global_offsets_pts;
    global_nelts = global_offsets_elts;

    //
    //std::cout << "local offset pts : " << offsets_pts << std::endl;
    //std::cout << "local offset elts : " << offsets_elts << std::endl;
    //std::cout << "global offset pts : " << global_offsets_pts << std::endl;
    //std::cout << "global offset elts : " << global_offsets_elts << std::endl;
    //std::cout << "done with offsets" << std::endl;
}

/**
 * Translate the list of points ids to the new global layout
 * @param ids Array of local point ids to be translated
 */
template<typename T>
int MeshPoints<T>::translatePointIds(std::vector<int32_t> & ptids)
{
    for(int i = 0; i < ptids.size(); i++)
    {
        ptids[i] = offsets_pts + old2new[ ptids[i] ];
    }

    return 0;
}

/**
 * Translate the list of element ids to the new global layout
 * @param ids Array of local point ids to be translated
 */
template<typename T>
int MeshPoints<T>::translateElementIds(std::vector<int32_t> & elids)
{
    for(int i = 0; i < elids.size(); i++)
    {
        elids[i] = offsets_elts + elids[i];
    }

    return 0;
}

}

} // Feel


//#if !defined(FEELPP_INSTANTIATION_MODE)
# include <feel/feeldiscr/meshimpl.hpp>
//#endif //

#endif /* FEELPP_MESH_HPP */
