/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-14

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)

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
   \file mesh3d.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-14
 */
#ifndef __Mesh3D_H
#define __Mesh3D_H 1


#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/numeric/ublas/io.hpp>



#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/visitor.hpp>

#include <feel/feelmesh/meshbase.hpp>

#include <feel/feelmesh/geoelement.hpp>

#include <feel/feelmesh/elements.hpp>
#include <feel/feelmesh/faces.hpp>
#include <feel/feelmesh/edges.hpp>
#include <feel/feelmesh/points.hpp>
#include <feel/feelmesh/functors.hpp>

namespace Feel
{
/**
 * \class Mesh3D
 * \brief 3D mesh class
 *
 * \code
 * // create a 3D mesh made of simplex of order 1
 * Mesh3D<Simplex<3,1> > mesh;
 *
 * // create a 3D mesh made of simplex of order 2
 * Mesh3D<Simplex<3,2> > mesh;
 * \endcode
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename Shape>
class Mesh3D
    :
public VisitableBase<>,
public MeshBase,
public Elements<Shape>,
public Points<3>,
public Faces<typename Shape::template shape<2>::type,
             typename Elements<Shape>::element_type>,
public Edges<typename Shape::template shape<1>::type,
             typename Faces<typename Shape::template shape<2>::type,
                            typename Elements<Shape>::element_type>::face_type >
{
    // check at compilation time that the shape has indeed dimension 2
    BOOST_STATIC_ASSERT( Shape::nDim == 3 );

    public:


    static const uint16_type nDim = 3;

    /** @name Typedefs
     */
    //@{
    typedef typename VisitableBase<>::return_type return_type;

    typedef VisitableBase<> super_visitable;
    typedef MeshBase super;

    typedef Elements<Shape> super_elements;
    typedef typename super_elements::elements_type elements_type;
    typedef typename super_elements::element_type element_type;
    typedef typename super_elements::element_iterator element_iterator;
    typedef typename super_elements::element_const_iterator element_const_iterator;
    typedef typename super_elements::update_element_neighbor_type update_element_neighbor_type;
    typedef typename element_type::node_type node_type;

    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;

    typedef Points<3> super_points;
    typedef typename super_points::points_type points_type;
    typedef typename super_points::point_type point_type;

    typedef Faces<typename Shape::template shape<2>::type,
    typename super_elements::element_type> super_faces;
    typedef typename super_faces::faces_type faces_type;
    typedef typename super_faces::face_type face_type;

    typedef typename super_faces::face_iterator face_iterator;
    typedef typename super_faces::face_const_iterator face_const_iterator;

    typedef typename super_faces::location_faces location_faces;
    typedef typename super_faces::location_face_iterator location_face_iterator;
    typedef typename super_faces::location_face_const_iterator location_face_const_iterator;


    typedef Edges<typename Shape::template shape<1>::type,face_type> super_edges;
    typedef typename super_edges::edges_type edges_type;
    typedef typename super_edges::edge_type edge_type;
    typedef typename super_edges::edge_iterator edge_iterator;
    typedef typename super_edges::edge_const_iterator edge_const_iterator;

    typedef typename std::pair<size_type, size_type> edge_pair_type;

    typedef Mesh3D<Shape> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

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
    Mesh3D( WorldComm const& worldComm = Environment::worldComm() );

    /**
     * copy constructor
     */
    Mesh3D( Mesh3D const & m );

    /**
     * destructor
     */
    ~Mesh3D();

    //@}

    /** @name Operator overloads
     */
    //@{

    Mesh3D& operator=( Mesh3D const& m );

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return \p true if all containers are empty, \p false otherwise
     */
    bool isEmpty() const
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
size_type numElements() const
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
size_type numFaces() const
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
size_type numPoints() const
{
    return this->points().size();
}

/**
 * \return the edge index of the edge \p n in the element \p e
 */
element_edge_type const& localEdgeId( element_type const& e,
                                      size_type const n ) const
{
    return M_e2e[e.id()][n];
}

/**
 * \return the edge index of the edge \p n in the element \p e
 */
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
virtual void setWorldComm( WorldComm const& _worldComm )
{
    this->setWorldCommMeshBase( _worldComm );
    this->setWorldCommElements( _worldComm );
    this->setWorldCommFaces( _worldComm );
    this->setWorldCommEdges( _worldComm );
    this->setWorldCommPoints( _worldComm );
}

/**
 * clear out all data from the mesh, \p isEmpty() should return
 * \p true after a \p clear()
 */
virtual void clear();


FEELPP_DEFINE_VISITABLE();

//@}



protected:

/**
 * dummy  implementation
 * \see Mesh
 */
void renumber()
{
    FEELPP_ASSERT( 0 ).error( "invalid call" );
}

/**
 * update the entities of co-dimension 2
 */
void updateEntitiesCoDimensionTwo();

/**
 * update permutation of entities of co-dimension 1
 */
void updateEntitiesCoDimensionOnePermutation();

#if 0
/**
 * update the faces information of the mesh
 */
void updateFaces();

/**
 * update the edges information of the mesh
 */
void updateEdges();

/**
 * check the mesh connectivity
 */
void check() const;

#endif // 0

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & boost::serialization::base_object<super>( *this );
            DVLOG(2) << "Serializing points\n";
            ar & boost::serialization::base_object<super_points>( *this );
            DVLOG(2) << "Serializing edges\n";
            ar & boost::serialization::base_object<super_edges>( *this );
            DVLOG(2) << "Serializing faces\n";
            ar & boost::serialization::base_object<super_faces>( *this );
            DVLOG(2) << "Serializing elements\n";
            ar & boost::serialization::base_object<super_elements>( *this );
        }


private:


/**
 * Determines the permutation a face given the global indices of the vertices (for tetrahedra)
 */
void determineFacePermutation( uint16_type numZeros, std::vector<size_type> const& def,
                               std::vector<size_type> const& cur, std::vector<uint32_type>& diff,
                               face_permutation_type& permutation, mpl::bool_<true> );

/**
* Determines the permutation a face given the global indices of the vertices (for hexahedra)
*/
void determineFacePermutation( uint16_type numZeros, std::vector<size_type> const& def,
                               std::vector<size_type> const& cur, std::vector<uint32_type>& diff,
                               face_permutation_type& permutation, mpl::bool_<false> );

/**
 * Arrays containing the global ids of edges of each element
 */
boost::multi_array<element_edge_type,2> M_e2e;
};

template <typename GEOSHAPE>
Mesh3D<GEOSHAPE>::Mesh3D( WorldComm const& worldComm )
    :
    super_visitable(),
    super( worldComm ),
    super_elements( worldComm ),
    super_points( worldComm ),
    super_faces( worldComm ),
    super_edges( worldComm ),
    M_e2e()
{}

template <typename GEOSHAPE>
Mesh3D<GEOSHAPE>::Mesh3D( Mesh3D const & m )
    :
    super_visitable(),
    super( m ),
    super_elements( m ),
    super_points( m ),
    super_faces( m ),
    super_edges( m ),
    M_e2e( m.M_e2e )
{}

template <typename GEOSHAPE>
Mesh3D<GEOSHAPE>::~Mesh3D()
{}

template <typename GEOSHAPE>
Mesh3D<GEOSHAPE>&
Mesh3D<GEOSHAPE>::operator=( Mesh3D const& m )
{
    if ( this != &m )
    {
        super::operator=( m );
        super_elements::operator=( m );
        super_points::operator=( m );
        super_faces::operator=( m );
        super_edges::operator=( m );

        M_e2e = m.M_e2e;
    }

    return *this;
}

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::clear()
{
    this->elements().clear();
    this->points().clear();
    this->faces().clear();
    this->edges().clear();

    M_e2e.resize( boost::extents[0][0] );
    FEELPP_ASSERT( isEmpty() ).error( "all mesh containers should be empty after a clear." );
}

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::determineFacePermutation( uint16_type numZeros, std::vector<size_type> const& def,
        std::vector<size_type> const& cur, std::vector<uint32_type>& diff,
        face_permutation_type& permutation, mpl::bool_<true> )
{
    if ( numZeros == 0 )
    {
        for ( uint16_type i = 0; i < def.size(); ++i )
            diff[i] = def[i] - cur[2-i];
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

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::determineFacePermutation( uint16_type numZeros, std::vector<size_type> const& def,
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

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::updateEntitiesCoDimensionOnePermutation()
{
    boost::timer ti;
    std::vector<size_type> _left( face_type::numVertices );
    std::vector<size_type> _right( face_type::numVertices );
    std::vector<uint32_type> _diff( face_type::numVertices );

    //determine permutation for the faces
    for ( face_iterator elt_it = this->beginFace();
            elt_it != this->endFace(); ++elt_it )
    {
        face_permutation_type permutation( face_permutation_type::IDENTITY );

        // if on boundary don't do anything
        if ( elt_it->isOnBoundary() || ( elt_it->pos_second() == invalid_uint16_type_value ) )
            continue;

        for ( uint16_type i = 0; i < face_type::numVertices; ++i )
        {
            _left[i] = elt_it->element0().point( elt_it->element0().fToP( elt_it->pos_first(), i ) ).id();

            uint16_type right_p = elt_it->element1().fToP( elt_it->pos_second(), i );
            FEELPP_ASSERT( right_p >= 0 && right_p < elt_it->element1().numLocalPoints )( right_p )( elt_it->element1().numLocalPoints )
            ( elt_it->pos_second() )( i ).error( "invalid point index" );
            _right[i] = elt_it->element1().point( right_p ).id();

            _diff[i] = _left[i] - _right[i];
        }

        uint16_type _numZeros = count( _diff.begin(), _diff.end(), uint32_type( 0 ) );

        determineFacePermutation( _numZeros, _left, _right, _diff,
                                  permutation, mpl::bool_<( SHAPE == SHAPE_TETRA )>() );

        if ( permutation.value() != face_permutation_type::IDENTITY )
            this->elements().modify( this->elementIterator( elt_it->ad_second(), elt_it->proc_second() ),
                                     Feel::detail::UpdateFacePermutation<face_permutation_type>( elt_it->pos_second(),
                                             permutation ) );
    }

    element_iterator iv,  en;
    boost::tie( iv, en ) = this->elementsRange();

    for ( ; iv != en; ++iv )
    {
        for ( size_type j = 0; j < numLocalFaces(); j++ )
        {
            FEELPP_ASSERT( iv->facePtr( j ) )( j )( iv->id() ).warn( "invalid element face check" );
        }
    }

    DVLOG(2) << "[Mesh3D::updateFaces] element/face permutation : " << ti.elapsed() << "\n";
}

template <typename GEOSHAPE>
void
Mesh3D<GEOSHAPE>::updateEntitiesCoDimensionTwo()
{
    boost::timer ti;
    std::map<std::set<int>, size_type > _edges;
    typename std::map<std::set<int>, size_type >::iterator _edgeit;
    int next_edge = 0;
    M_e2e.resize( boost::extents[this->numElements()][this->numLocalEdges()] );

    size_type vid, i1, i2;
    element_type ele;
    face_type bele;

    // First We check if we have already Edges stored
    if ( ! this->edges().empty() )
    {
        // dump first the existing edges, to maintain the correct numbering
        // if everything is correct, the numbering structure will reflect
        // the actual edge numbering

        for ( size_type j = 0; j < this->edges().size(); ++j )
        {
            std::set<int> s;
            i1 = ( this->edge( j ).point( 0 ) ).id();
            i2 = ( this->edge( j ).point( 1 ) ).id();
            s.insert( i1 );
            s.insert( i2 );
            bool edgeinserted = false;

            boost::tie( _edgeit, edgeinserted ) = _edges.insert( std::make_pair( s, next_edge ) );

            if ( edgeinserted )
                ++next_edge;

            FEELPP_ASSERT( edgeinserted )( i1 )( i2 )(j)( this->edge( j ).id() )( _edgeit->second )( this->edge( j ).marker() )
                ( this->edge( _edgeit->second ).point( 0 ).node() )( this->edge( _edgeit->second ).point( 1 ).node() )
                ( this->edge( j ).point( 0 ).node() )( this->edge( j ).point( 1 ).node() ).error( "Two identical Edges stored in EdgeList" );
            FEELPP_ASSERT( _edgeit->second == this->edge( j ).id() )( _edgeit->second )( this->edge( j ).id() ).error( "Edges in EdgeList have inconsistent id" );

        }
    }

    DVLOG(2) << "[Mesh3D::updateEdges] adding edges : " << ti.elapsed() << "\n";
    ti.restart();

    edge_type edg;//(this->worldComm());
    edg.setProcessIdInPartition( this->worldComm().localRank() );

    if ( this->edges().empty() )
    {
        // We want that the first edges be those on the boundary, in order to obey the paradigm for
        // a Mesh3D
        const location_faces& location_index=this->faces().template get<Feel::detail::by_location>();
        location_face_iterator ifa;
        location_face_iterator efa;
        boost::tie( ifa, efa ) = location_index.equal_range( boost::make_tuple( true ) );

        for ( ; ifa!=efa; ++ifa )
        {
            for ( uint16_type j = 0; j < face_type::numEdges; j++ )
            {
                std::set<int> s;
                i1 = bele.eToP( j, 0 );
                i2 = bele.eToP( j, 1 );
                // go to global
                i1 = ( ifa->point( i1 ) ).id();
                i2 = ( ifa->point( i2 ) ).id();
                s.insert( i1 );
                s.insert( i2 );
                bool edgeinserted = false;
                boost::tie( _edgeit, edgeinserted ) = _edges.insert( std::make_pair( s, next_edge ) );

                if ( edgeinserted )
                    ++next_edge;

                // reset the process id (edge not connected to an active elt take this value)
                edg.setProcessId( invalid_rank_type_value );

                if ( edgeinserted )
                {
                    // set edge id
                    edg.setId( _edgeit->second );
                    edg.setOnBoundary( true, 0 );

                    for ( uint16_type k = 0; k < 2 + face_type::nbPtsPerEdge; k++ )
                        edg.setPoint( k, ifa->point( bele.eToP( j, k ) ) );

                    // TODO: should assocate a marker to the edge here ?
                    //edg.addElement( ifa->ad_first() );
                    this->addEdge( edg );
                }
                // set the process id from element (only active element)
                if (!ifa->isGhostCell()) this->edges().modify( this->edgeIterator( _edgeit->second ), Feel::detail::UpdateProcessId(ifa->processId()) );
            }
        }
    }

    DVLOG(2) << "[Mesh3D::updateEdges] adding edges : " << ti.elapsed() << "\n";
    ti.restart();

    std::map<size_type, edge_pair_type> _oriented_edges;
    typedef typename std::map<size_type, edge_pair_type>::iterator oe_iterator;

    for ( element_iterator elt_it = this->beginElement();
            elt_it != this->endElement(); ++elt_it )
    {
        vid = elt_it->id();

        for ( uint16_type j = 0; j < element_type::numEdges; ++j )
        {
            uint16_type j1 = ele.eToP( j, 0 );
            uint16_type j2 = ele.eToP( j, 1 );

            // go to global
            i1 = ( elt_it->point( j1 ) ).id();
            i2 = ( elt_it->point( j2 ) ).id();
            std::set<int> s;
            s.insert( i1 );
            s.insert( i2 );
            bool edgeinserted = false;
            boost::tie( _edgeit, edgeinserted ) = _edges.insert( std::make_pair( s, next_edge ) );

            if ( edgeinserted )
                ++next_edge;

            M_e2e[ vid ][ j] = boost::make_tuple( _edgeit->second, 1 );

            // reset the process id (edge not connected to an active elt take this value)
            edg.setProcessId( invalid_rank_type_value );

            if ( edgeinserted )
            {
                FEELPP_ASSERT( _edgeit->second >= this->numEdges() )( _edgeit->second )( this->numEdges() ).error( "invalid edge index" );
                // set edge id
                edg.setId( _edgeit->second );
                // we have already inserted edges on the boundary so
                // this one _is_ not on the boundary
                edg.setOnBoundary( false );

                if ( this->components().test( MESH_ADD_ELEMENTS_INFO ) )
                    edg.addElement( vid );

                // number of points on the edge is 2 (number of
                // vertices) plus the number of points in the
                // interior of the edge
                for ( uint16_type k = 0; k < 2 + element_type::nbPtsPerEdge; k++ )
                    edg.setPoint( k, elt_it->point( elt_it->eToP( j, k ) ) );

                // TODO: should assocate a marker to the edge here ?
                //edg.addElement( vid );
                this->addEdge( edg );
            }

            else
            {
                if ( this->components().test( MESH_ADD_ELEMENTS_INFO ) )
                    this->edges().modify( this->edgeIterator( _edgeit->second ), [vid] ( edge_type& e ) { e.addElement( vid ); } );
            }

            // set the process id from element (only active element)
            if (!elt_it->isGhostCell()) this->edges().modify( this->edgeIterator( _edgeit->second ), Feel::detail::UpdateProcessId(elt_it->processId()) );

            this->elements().modify( elt_it,
                                     Feel::detail::UpdateEdge<edge_type>( j, boost::cref( this->edge( _edgeit->second ) ) ) );
#if !defined(NDEBUG)
            this->elements().modify( elt_it,
                                     [j]( element_type const& e ) { FEELPP_ASSERT( e.edgePtr( j ) )( e.id() )( j ).error( "invalid edge in element" ); } );
#endif
        }
        for ( uint16_type j = 0; j < element_type::numFaces; ++j )
        {
            auto fit = this->faces().iterator_to( elt_it->face(j));
            for ( uint16_type e = 0; e < face_type::numEdges; ++e )
            {
                auto const& elt_edge = elt_it->edge( elt_it->f2e( j, e ) );
                this->faces().modify( fit,
                                      [e,&elt_edge]( face_type& f ) { f.setEdge(e,elt_edge); } );
            }
        }
    }

    DVLOG(2) << "[Mesh3D::updateEdges] updating element/edges : " << ti.elapsed() << "\n";
    ti.restart();

    for ( element_iterator elt_it = this->beginElement();
            elt_it != this->endElement(); ++elt_it )
    {
        for ( size_type j = 0; j < ( size_type )element_type::numEdges; ++j )
        {
            size_type j1 = ele.eToP( j, 0 );
            size_type j2 = ele.eToP( j, 1 );

            // go to global
            i1 = ( elt_it->point( j1 ) ).id();
            i2 = ( elt_it->point( j2 ) ).id();
            std::set<int> s;
            s.insert( i1 );
            s.insert( i2 );
            bool edgeinserted = false;
            boost::tie( _edgeit, edgeinserted ) = _edges.insert( std::make_pair( s, next_edge ) );

            if ( edgeinserted )
                ++next_edge;

            edge_pair_type _current = std::make_pair( i1, i2 );

            edge_permutation_type permutation( edge_permutation_type::IDENTITY );

            oe_iterator _edge_it = _oriented_edges.find( _edgeit->second );

            if (  _edge_it != _oriented_edges.end() )
            {
                edge_pair_type _default = _edge_it->second;

                FEELPP_ASSERT( _default.first == _current.first ||
                               _default.first == _current.second ).error( "invalid edge index" );

                if ( _default.first != _current.first )
                {
                    permutation = edge_permutation_type::REVERSE_PERMUTATION;
                    this->elements().modify( elt_it,
                                             Feel::detail::UpdateEdgePermutation<edge_permutation_type>( j,
                                                     permutation ) );
                }
            }

            else
            {
                _oriented_edges.insert( std::make_pair( _edgeit->second, _current ) );
            }

        }
    }

    DVLOG(2) << "[Mesh3D::updateEdges] updating edges orientation : " << ti.elapsed() << "\n";
    ti.restart();
#if 0
    edge_iterator e_it = this->beginEdge();
    edge_iterator e_en = this->endEdge();

    for ( ; e_it!=e_en; ++e_it )
    {
        // cleanup the edge data structure :

        if ( e_it->numberOfElements() == 0 )
        {
            // remove all edges that are not connected to any elements
            this->edges().erase( e_it );
        }

    }

    DVLOG(2) << "[Mesh3D::updateEdges] cleaning up edges : " << ti.elapsed() << "\n";
#endif
    ti.restart();
}

} // Feel

#endif /* __Mesh3D_H */
