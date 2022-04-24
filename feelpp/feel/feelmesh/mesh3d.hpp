// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
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
//! @copyright 2005,2006 EPFL
//! @copyright 2007-2010 Université Joseph Fourier (Grenoble I)
//! @copyright 2011-2017 Feel++ Consortium
//!
#ifndef FEELPP_MESH3D_HPP
#define FEELPP_MESH3D_HPP 1

#include <cstdlib>
#include <fstream>
#include <iomanip>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/visitor.hpp>

#include <feel/feelmesh/meshbase.hpp>

#include <feel/feelmesh/geoelement.hpp>

#include <feel/feelmesh/edges.hpp>
#include <feel/feelmesh/elements.hpp>
#include <feel/feelmesh/faces.hpp>
#include <feel/feelmesh/functors.hpp>
#include <feel/feelmesh/points.hpp>

namespace Feel
{
//!
//! @brief 3D mesh class
//! @ingroup Mesh
//!
//! @code
//! // create a 3D mesh made of simplex of order 1
//! Mesh3D<Simplex<3,1> > mesh;
//!
//! // create a 3D mesh made of simplex of order 2
//! Mesh3D<Simplex<3,2> > mesh;
//! @endcode
//!
//!  @author Christophe Prud'homme
//!  @see
//!
template <typename Shape, typename T = double, typename IndexT = uint32_type>
class Mesh3D
    : public VisitableBase<>,
      public MeshBase<IndexT>,
      public Elements<Shape, T>,
      public Points<3, T>,
      public Faces<typename Shape::template shape<2>,
                   typename Elements<Shape, T>::element_type>,
      public Edges<typename Shape::template shape<1>,
                   typename Faces<typename Shape::template shape<2>,
                                  typename Elements<Shape, T>::element_type>::face_type>
{
    // check at compilation time that the shape has indeed dimension 2
    BOOST_STATIC_ASSERT( Shape::nDim == 3 );

  public:
    static const uint16_type nDim = 3;
    static const uint16_type nRealDim = 3;

    /** @name Typedefs
     */
    //@{
    typedef typename VisitableBase<>::return_type return_type;

    typedef VisitableBase<> super_visitable;
    typedef MeshBase<IndexT> super;

    using index_type = typename super::index_type;
    using size_type = typename super::size_type;
    
    typedef Elements<Shape, T> super_elements;
    typedef typename super_elements::elements_type elements_type;
    typedef typename super_elements::element_type element_type;
    typedef typename super_elements::element_iterator element_iterator;
    typedef typename super_elements::element_const_iterator element_const_iterator;
    typedef typename super_elements::update_element_neighbor_type update_element_neighbor_type;
    typedef typename element_type::node_type node_type;

    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;

    typedef Points<3, T> super_points;
    typedef typename super_points::points_type points_type;
    typedef typename super_points::point_type point_type;

    typedef Faces<typename Shape::template shape<2>,
                  typename super_elements::element_type>
        super_faces;
    typedef typename super_faces::faces_type faces_type;
    typedef typename super_faces::face_type face_type;

    typedef typename super_faces::face_iterator face_iterator;
    typedef typename super_faces::face_const_iterator face_const_iterator;

    // typedef typename super_faces::location_faces location_faces;
    typedef typename super_faces::face_reference_wrapper_iterator location_face_iterator;
    typedef typename super_faces::face_reference_wrapper_const_iterator location_face_const_iterator;

    typedef Edges<typename Shape::template shape<1>, face_type> super_edges;
    typedef typename super_edges::edges_type edges_type;
    typedef typename super_edges::edge_type edge_type;
    typedef typename super_edges::edge_iterator edge_iterator;
    typedef typename super_edges::edge_const_iterator edge_const_iterator;

    typedef typename std::pair<size_type, size_type> edge_pair_type;

    typedef Mesh3D<Shape, T> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    static const size_type SHAPE = Shape::Shape;

    typedef typename super::face_processor_type face_processor_type;

    /**
     * Tuple that contains
     *
     * -# the index of the edge
     *
     * -# +1 or -1 depending on the orientation
     */
    typedef boost::tuple<size_type, int> element_edge_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     */
    explicit Mesh3D( std::string const& name = "", worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

    Mesh3D( Mesh3D const& m ) = default;
    Mesh3D( Mesh3D&& m ) = default;

    /**
     * destructor
     */
    ~Mesh3D() override;

    //@}

    /** @name Operator overloads
     */
    //@{

    Mesh3D& operator=( Mesh3D const& m ) = default;
    Mesh3D& operator=( Mesh3D&& m ) = default;

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return \p true if all containers are empty, \p false otherwise
     */
    bool isEmpty() const override
    {
        return ( super_elements::isEmpty() &&
                 super_points::isEmpty() &&
                 super_faces::isEmpty() &&
                 super_edges::isEmpty() &&
                 M_e2e.empty() );
    }

    /**
 * \return the number of elements
 */
    size_type numElements() const override
    {
        return this->elements().size();
    }

    /**
 * \return the number of faces in an element
 */
    size_type numLocalFaces() const
    {
        return super_elements::element_type::numLocalFaces;
    }

    //!
    //! the number of topological faces per element
    //! @return the number of topological faces per element
    //!
    uint16_type numLocalTopologicalFaces() const
    {
        return super_elements::element_type::numTopologicalFaces;
    }

    /**
 * \return the number of edges in an element
 */
    size_type numLocalEdges() const
    {
        return super_elements::element_type::numLocalEdges;
    }

    /**
 * \return the number of vertices in an element
 */
    size_type numLocalVertices() const
    {
        return super_elements::element_type::numLocalVertices;
    }

    /**
 * \return the number of faces
 */
    size_type numFaces() const override
    {
        return this->faces().size();
    }

    /**
 * \return the number of edges
 */
    size_type numEdges() const
    {
        return this->edges().size();
    }

    /**
 * \return the number of points
 */
    size_type numPoints() const override
    {
        return this->points().size();
    }

    /**
 * \return the edge index of the edge \p n in the element \p e
 */
    FEELPP_DEPRECATED
    element_edge_type const& localEdgeId( element_type const& e,
                                          size_type const n ) const
    {
        return M_e2e[e.id()][n];
    }

    /**
 * \return the edge index of the edge \p n in the element \p e
 */
    FEELPP_DEPRECATED
    element_edge_type const& localEdgeId( size_type const e,
                                          size_type const n ) const
    {
        return M_e2e[e][n];
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{
    void setWorldComm( worldcomm_ptr_t const& _worldComm ) override
    {
        MeshBase<IndexT>::setWorldComm( _worldComm );
        this->setWorldCommElements( _worldComm );
        this->setWorldCommFaces( _worldComm );
        this->setWorldCommEdges( _worldComm );
        this->setWorldCommPoints( _worldComm );
    }

    /**
     * clear out all data from the mesh, \p isEmpty() should return
     * \p true after a \p clear()
     */
    void clear() override;

    FEELPP_DEFINE_VISITABLE();

    //@}

  protected:
    /**
     * dummy  implementation
     * \see Mesh
     */
    void renumber() override
    {
        FEELPP_ASSERT( 0 )
            .error( "invalid call" );
    }

    /**
     * update the entities of co-dimension 2
     */
    void updateEntitiesCoDimensionTwo() override;

    /**
     * update permutation of entities of co-dimension 1
     */
    void updateEntitiesCoDimensionOnePermutation();

  private:
    friend class boost::serialization::access;
    template <class Archive>
    FEELPP_NO_EXPORT void serialize( Archive& ar, const unsigned int version )
    {
        ar& boost::serialization::base_object<super>( *this );
        DVLOG( 2 ) << "Serializing points\n";
        ar& boost::serialization::base_object<super_points>( *this );
        DVLOG( 2 ) << "Serializing edges\n";
        ar& boost::serialization::base_object<super_edges>( *this );
        DVLOG( 2 ) << "Serializing faces\n";
        ar& boost::serialization::base_object<super_faces>( *this );
        DVLOG( 2 ) << "Serializing elements\n";
        ar& boost::serialization::base_object<super_elements>( *this );
    }

  private:
    /**
     * Determines the permutation a face given the global indices of the vertices (for tetrahedra)
     */
    FEELPP_NO_EXPORT void determineFacePermutation( uint16_type numZeros, std::vector<size_type> const& def,
                                                    std::vector<size_type> const& cur, std::vector<uint32_type>& diff,
                                                    face_permutation_type& permutation, mpl::bool_<true> );

    /**
     * Determines the permutation a face given the global indices of the vertices (for hexahedra)
     */
    FEELPP_NO_EXPORT void determineFacePermutation( uint16_type numZeros, std::vector<size_type> const& def,
                                                    std::vector<size_type> const& cur, std::vector<uint32_type>& diff,
                                                    face_permutation_type& permutation, mpl::bool_<false> );

    /**
     * Arrays containing the global ids of edges of each element
     */
    boost::multi_array<element_edge_type, 2> M_e2e;
};

template <typename GEOSHAPE, typename T, typename IndexT>
Mesh3D<GEOSHAPE, T, IndexT>::Mesh3D( std::string const& name, worldcomm_ptr_t const& worldComm )
    : super_visitable(),
      super( name, 3, nRealDim, worldComm ),
      super_elements( worldComm ),
      super_points( worldComm ),
      super_faces( worldComm ),
      super_edges( worldComm ),
      M_e2e()
{
}
template <typename GEOSHAPE, typename T, typename IndexT>
Mesh3D<GEOSHAPE, T, IndexT>::~Mesh3D()
{
}
template <typename GEOSHAPE, typename T, typename IndexT>
void Mesh3D<GEOSHAPE, T, IndexT>::clear()
{
    super::clear();
    super_elements::clear();
    super_points::clear();
    super_faces::clear();
    super_edges::clear();

    M_e2e.resize( boost::extents[0][0] );
    CHECK( isEmpty() ) << "all mesh containers should be empty after a clear.";
}

template <typename GEOSHAPE, typename T, typename IndexT>
void Mesh3D<GEOSHAPE, T, IndexT>::determineFacePermutation( uint16_type numZeros, std::vector<size_type> const& def,
                                                    std::vector<size_type> const& cur, std::vector<uint32_type>& diff,
                                                    face_permutation_type& permutation, mpl::bool_<true> )
{
    if ( numZeros == 0 )
    {
        for ( uint16_type i = 0; i < def.size(); ++i )
            diff[i] = def[i] - cur[2 - i];
    }

    std::vector<uint32_type>::iterator _id_it = find( diff.begin(),
                                                      diff.end(),
                                                      uint32_type( 0 ) );

    uint16_type pos = distance( diff.begin(), _id_it );

    if ( numZeros == 0 )
    {
        if ( pos == 0 )
            permutation = face_permutation_type::ROTATION_CLOCKWISE;

        else
            permutation = face_permutation_type::ROTATION_ANTICLOCK;
    }

    else if ( numZeros == 1 )
    {
        if ( pos == 0 )
            permutation = face_permutation_type::REVERSE_HYPOTENUSE;

        else if ( pos == 1 )
            permutation = face_permutation_type::REVERSE_HEIGHT;

        else
            permutation = face_permutation_type::REVERSE_BASE;
    }
}

template <typename GEOSHAPE, typename T, typename IndexT>
void Mesh3D<GEOSHAPE, T, IndexT>::determineFacePermutation( uint16_type numZeros, std::vector<size_type> const& def,
                                                    std::vector<size_type> const& cur, std::vector<uint32_type>& diff,
                                                    face_permutation_type& permutation, mpl::bool_<false> )
{
    std::vector<uint32_type>::iterator _id_it = find( diff.begin(),
                                                      diff.end(),
                                                      uint32_type( 0 ) );

    uint16_type pos = distance( diff.begin(), _id_it );

    if ( numZeros == 2 )
    {
        if ( pos == 0 )
            permutation = face_permutation_type::SECOND_DIAGONAL;

        else
            permutation = face_permutation_type::PRINCIPAL_DIAGONAL;
    }

    else if ( numZeros == 0 )
    {
        if ( cur[0] == def[1] )
            if ( cur[2] == def[3] )
                permutation = face_permutation_type::REVERSE_BASE;

            else
                permutation = face_permutation_type::ROTATION_CLOCKWISE;

        else if ( cur[0] == def[2] )
            if ( cur[2] == def[0] )
                permutation = face_permutation_type::ROTATION_ANTICLOCK;

            else
                permutation = face_permutation_type::REVERSE_HEIGHT;

        else
            permutation = face_permutation_type::ROTATION_TWICE_CLOCKWISE;
    }
}

template <typename GEOSHAPE, typename T, typename IndexT>
void Mesh3D<GEOSHAPE, T, IndexT>::updateEntitiesCoDimensionOnePermutation()
{
    tic();
    std::vector<size_type> _left( face_type::numVertices );
    std::vector<size_type> _right( face_type::numVertices );
    std::vector<uint32_type> _diff( face_type::numVertices );

    //determine permutation for the faces
    for ( face_iterator elt_it = this->beginFace();
          elt_it != this->endFace(); ++elt_it )
    {
        auto const& face = elt_it->second;
        face_permutation_type permutation( face_permutation_type::IDENTITY );

        // if on boundary don't do anything
        if ( face.isOnBoundary() || ( face.pos_second() == invalid_uint16_type_value ) )
            continue;

        bool applyOnFirstConnection = false;
        if ( face.isInterProcessDomain() )
        {
            if ( face.partition1() < face.partition2() )
            {
                if ( face.proc_first() == face.partition2() )
                     applyOnFirstConnection = true;
            }
            else
            {
                if ( face.proc_first() == face.partition1() )
                    applyOnFirstConnection = true;
            }
        }

        for ( uint16_type i = 0; i < face_type::numVertices; ++i )
        {
            _left[i] = face.element0().point( face.element0().fToP( face.pos_first(), i ) ).id();
            _right[i] = face.element1().point( face.element1().fToP( face.pos_second(), i ) ).id();
            _diff[i] = _left[i] - _right[i];
        }

        uint16_type _numZeros = count( _diff.begin(), _diff.end(), uint32_type( 0 ) );

        if ( applyOnFirstConnection )
            determineFacePermutation( _numZeros, _right, _left, _diff,
                                      permutation, mpl::bool_<( SHAPE == SHAPE_TETRA )>() );
        else
            determineFacePermutation( _numZeros, _left, _right, _diff,
                                      permutation, mpl::bool_<( SHAPE == SHAPE_TETRA )>() );

        if ( permutation.value() != face_permutation_type::IDENTITY )
        {
            if ( applyOnFirstConnection )
                this->elementIterator( face.ad_first() )->second.setFacePermutation( face.pos_first(), permutation );
            else
                this->elementIterator( face.ad_second() )->second.setFacePermutation( face.pos_second(), permutation );
        }
    }

#if !defined( NDEBUG )
    if ( this->components().test( MESH_UPDATE_FACES ) )
    {
        element_iterator iv, en;
        boost::tie( iv, en ) = this->elementsRange();
        for ( ; iv != en; ++iv )
        {
            auto const& elt = iv->second;
            for ( size_type j = 0; j < numLocalFaces(); j++ )
            {
                FEELPP_ASSERT( elt.facePtr( j ) )
                ( j )( elt.id() ).warn( "invalid element face check" );
            }
        }
    }
#endif
    toc( "[Mesh3D::updateFaces] element/face permutation", FLAGS_v > 1 );
}

template <typename GEOSHAPE, typename T, typename IndexT>
void Mesh3D<GEOSHAPE, T, IndexT>::updateEntitiesCoDimensionTwo()
{
    rank_type currentPid = this->worldComm().localRank();

    //typedef std::unordered_map<std::set<size_type>, edge_type*, Feel::HashTables::HasherContainers<size_type> > pointstoedge_container_type;
    typedef std::unordered_map<std::vector<size_type>, std::tuple<edge_type*, size_type>, Feel::HashTables::HasherContainers<size_type>> pointstoedge_container_type;
    pointstoedge_container_type _edges;
    typename pointstoedge_container_type::iterator _edgeit;

    size_type next_edge = 0;
    bool edgeinserted = false;
    std::vector<size_type> lids( edge_type::numVertices );

    size_type vid, i1, i2;
    const bool updateComponentAddElements = this->components().test( MESH_ADD_ELEMENTS_INFO );
    tic();
    // First We check if we have already Edges stored
    if ( !this->edges().empty() )
    {
        // dump first the existing edges, to maintain the correct numbering
        // if everything is correct, the numbering structure will reflect
        // the actual edge numbering
        edge_iterator eit = this->beginEdge();
        edge_iterator een = this->endEdge();
        for ( ; eit != een; )
        {
            auto& edge = eit->second;
            i1 = edge.point( 0 ).id();
            i2 = edge.point( 1 ).id();
            //std::set<size_type> s( {i1, i2} );
            if ( i1 < i2 )
            {
                lids[0] = i1;
                lids[1] = i2;
            }
            else
            {
                lids[1] = i1;
                lids[0] = i2;
            }

            boost::tie( _edgeit, edgeinserted ) = _edges.emplace( std::piecewise_construct,
                                                                  std::forward_as_tuple( lids ),
                                                                  std::forward_as_tuple( &edge, invalid_v<size_type> ) );

            if ( edgeinserted )
            {
                next_edge = std::max( next_edge, edge.id() + 1 );
                ++eit;
            }
            else
            {
                eit = this->eraseEdge( eit );
            }
#if 0
            FEELPP_ASSERT( edgeinserted )( i1 )( i2 )(j)( this->edge( j ).id() )( _edgeit->second )( this->edge( j ).hasMarker() )
                ( this->edge( _edgeit->second ).point( 0 ).node() )( this->edge( _edgeit->second ).point( 1 ).node() )
                ( this->edge( j ).point( 0 ).node() )( this->edge( j ).point( 1 ).node() ).error( "Two identical Edges stored in EdgeList" );
            FEELPP_ASSERT( _edgeit->second == this->edge( j ).id() )( _edgeit->second )( this->edge( j ).id() ).error( "Edges in EdgeList have inconsistent id" );
#endif
        }
    }
    toc( "[Mesh3D::updateEdges] adding edges already registered", FLAGS_v > 1 );
    tic();

    if ( true ) //this->edges().empty() )
    {
        // We want that the first edges be those on the boundary, in order to obey the paradigm for
        // a Mesh3D
        auto rangeBoundaryFaces = this->facesOnBoundary();
        auto ifa = std::get<0>( rangeBoundaryFaces );
        auto efa = std::get<1>( rangeBoundaryFaces );
        for ( ; ifa != efa; ++ifa )
        {
            auto& bface = boost::unwrap_ref( *ifa );
            for ( uint16_type j = 0; j < face_type::numEdges; j++ )
            {
                i1 = bface.point( face_type::eToP( j, 0 ) ).id();
                i2 = bface.point( face_type::eToP( j, 1 ) ).id();
                //std::set<size_type> s( { i1,i2 } );
                if ( i1 < i2 )
                {
                    lids[0] = i1;
                    lids[1] = i2;
                }
                else
                {
                    lids[1] = i1;
                    lids[0] = i2;
                }

                boost::tie( _edgeit, edgeinserted ) = _edges.emplace( std::piecewise_construct,
                                                                      std::forward_as_tuple( lids ),
                                                                      std::forward_as_tuple( nullptr, invalid_v<size_type> ) );

                if ( edgeinserted )
                {
                    edge_type edg;
                    edg.setProcessIdInPartition( currentPid );
                    edg.setId( next_edge++ );
                    edg.setOnBoundary( true, 0 );
                    edg.setProcessId( invalid_rank_type_value );
                    for ( uint16_type k = 0; k < 2 + face_type::nbPtsPerEdge; k++ )
                        edg.setPoint( k, const_cast<point_type&>( bface.point( face_type::eToP( j, k ) ) ) );

                    auto res = this->addEdge( std::move(edg) );
                    auto& edgeInserted = res.first->second;
                    std::get<0>( _edgeit->second ) = &edgeInserted;
                }
                else
                {
                    auto edgePtr = std::get<0>( _edgeit->second );
                    if ( !edgePtr->isOnBoundary() )
                        edgePtr->setOnBoundary( true, 0 );
                }
            }
        }
    }
    toc( "[Mesh3D::updateEdges] adding boundaryfaces/edges", FLAGS_v > 1 );
    tic();
    edge_permutation_type reversePermutation( edge_permutation_type::REVERSE_PERMUTATION );

    std::vector<uint16_type> myeToP( edge_type::numVertices * element_type::numEdges );
    for ( uint16_type j = 0; j < element_type::numEdges; j++ )
    {
        for ( uint16_type f = 0; f < edge_type::numVertices; ++f )
            myeToP[j * edge_type::numVertices + f] = element_type::eToP( j, f );
    }
    Eigen::Matrix<size_type, element_type::numVertices, 1> pointIdInElt;

    auto elt_it = this->beginOrderedElement();
    auto elt_en = this->endOrderedElement();
    for ( ; elt_it != elt_en; ++elt_it )
    {
        auto& elt = unwrap_ref( *elt_it );
        rank_type eltPid = elt.processId();
        vid = elt.id();

        for ( uint16_type f = 0; f < element_type::numVertices; ++f )
            pointIdInElt[f] = elt.point( f ).id();

        for ( uint16_type j = 0; j < element_type::numEdges; ++j )
        {
            //auto const& pt0 = elt.point( element_type::eToP( j, 0 ) );
            //auto const& pt1 = elt.point( element_type::eToP( j, 1 ) );
            //i1 = pt0.id();
            //i2 = pt1.id();
            i1 = pointIdInElt[myeToP[j * edge_type::numVertices + 0]];
            i2 = pointIdInElt[myeToP[j * edge_type::numVertices + 1]];
            //std::set<size_type> s( {i1, i2} );
            if ( i1 < i2 )
            {
                lids[0] = i1;
                lids[1] = i2;
            }
            else
            {
                lids[1] = i1;
                lids[0] = i2;
            }

            boost::tie( _edgeit, edgeinserted ) = _edges.emplace( std::piecewise_construct,
                                                                  std::forward_as_tuple( lids ),
                                                                  std::forward_as_tuple( nullptr, i1 ) );

            if ( edgeinserted )
            {
                edge_type edg;
                edg.setProcessIdInPartition( currentPid );
                edg.setProcessId( ( elt.isGhostCell() ) ? invalid_rank_type_value : eltPid );
                edg.setId( next_edge++ );
                edg.setOnBoundary( false );
                // update connected element in edge
                if ( updateComponentAddElements )
                {
                    edg.elements().clear();
                    edg.addElement( vid, j );
                }

                // number of points on the edge is 2 (number of
                // vertices) plus the number of points in the
                // interior of the edge
                auto& pt0 = elt.point( element_type::eToP( j, 0 ) );
                auto& pt1 = elt.point( element_type::eToP( j, 1 ) );
                edg.setPoint( 0, pt0 );
                edg.setPoint( 1, pt1 );
                for ( uint16_type k = 2; k < 2 + element_type::nbPtsPerEdge; k++ )
                    edg.setPoint( k, elt.point( element_type::eToP( j, k ) ) );

                // add edge in mesh container
                auto res = this->addEdge( std::move(edg) );
                auto& edgeInserted = res.first->second;
                // update edge pointer in element
                elt.setEdge( j, boost::cref( edgeInserted ) );
                std::get<0>( _edgeit->second ) = &edgeInserted;

                if ( i1 > i2 )
                    elt.setEdgePermutation( j, reversePermutation );
            }
            else
            {
                auto edgePtr = std::get<0>( _edgeit->second );
                // update edge pointer in element
                elt.setEdge( j, boost::cref( *edgePtr ) );
                // set the process id from element (only active element)
                if ( !elt.isGhostCell() && edgePtr->processId() != eltPid )
                    edgePtr->setProcessId( eltPid );
                if ( updateComponentAddElements || edgePtr->hasMarker() )
                {
                    //DLOG_IF(INFO, eit->marker().isOn()) << "found edge " << eit->id() << " with marker:" << eit->marker() << ", adding element id : " << vid <<  "  local edge id " << j;
                    edgePtr->addElement( vid, j );
                }
#if 0
                // update edge orientation in element
                size_type& ptForPermutation = std::get<1>( _edgeit->second );
                if ( ptForPermutation == invalid_v<size_type> )
                    ptForPermutation = i1;
                else if ( ptForPermutation != i1 )
                    elt.setEdgePermutation( j, reversePermutation );
                // if ( i1 != edgePtr->point( 0 ).id() )
                //     elt.setEdgePermutation( j, reversePermutation );
#endif
                if ( i1 > i2 )
                    elt.setEdgePermutation( j, reversePermutation );
            }

        } // for ( uint16_type j = 0; j < element_type::numEdges; ++j )
    }
    toc( "[Mesh3D::updateEdges] adding element/edges", FLAGS_v > 1 );
    tic();
    // update edge pointers in faces
    face_iterator face_it = this->beginFace();
    face_iterator face_en = this->endFace();
    for ( ; face_it != face_en; ++face_it )
    {
        auto& faceModified = face_it->second;
        if ( !faceModified.isConnectedTo0() )
            continue;
        auto const& elt0 = faceModified.element0();
        uint16_type j = faceModified.pos_first();
        for ( uint16_type e = 0; e < face_type::numEdges; ++e )
        {
            auto const& elt_edge = elt0.edge( elt0.f2e( j, e ) );
            faceModified.setEdge( e, elt_edge );
        }
    }
    toc( "[Mesh3D::updateEdges] updating faces/edges", FLAGS_v > 1 );
#if 0
    edge_iterator e_it = this->beginEdge();
    edge_iterator e_en = this->endEdge();

    for ( ; e_it!=e_en; ++e_it )
    {
        // cleanup the edge data structure :

        if ( e_it->numberOfElements() == 0 )
        {
            // remove all edges that are not connected to any elements
            this->eraseEdge( e_it );
        }

    }

    DVLOG(2) << "[Mesh3D::updateEdges] cleaning up edges : " << ti.elapsed() << "\n";
#endif
}

} // namespace Feel

#endif /* __Mesh3D_H */
