/*
 This file is part of the Feel library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
 * \file geo0D.hpp
 */
#ifndef _GEO0D_HH_
#define _GEO0D_HH_

#include <boost/operators.hpp>

#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/marker.hpp>
#include <feel/feelmesh/meshbase.hpp>

namespace Feel
{
template <typename IndexT>
class MeshBase;

/**
 *  \defgroup GeoXD Basis Geometrical Entities Geo0D and GeoND.
 *  \ingroup Obsolet_Groups
 *  They are intermediate classes used to build the actual Geometry classes
 *
 *  \warning Geo1D/2D/3D are template classes because some of the info is not
 *  known a priori and I want all vector dimensions determined at compile time
 *  to enhance memory access time.
 */

/*@{*/

/**
 * \class Geo0D
 *
 *  Zero dimensional entity.
 *
 * \ingroup GeoXD
 */
template <uint16_type Dim, typename T = double, typename IndexT = uint32_type>
class Geo0D
    : public boost::equality_comparable<Geo0D<Dim, T>>,
      public boost::less_than_comparable<Geo0D<Dim, T>>,
      public boost::less_than_comparable<Geo0D<Dim, T>, IndexT>,
      public GeoEntity<Simplex<0, 1, Dim>, T, IndexT>
{
    typedef GeoEntity<Simplex<0, 1, Dim>, T, IndexT> super;

  public:
    typedef Geo0D<Dim, T, IndexT> self_type;
    static const uint16_type nDim = Dim;
    static const bool is_simplex = true;
    typedef T value_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;
    using node_type = typename node<T, Dim>::type;
    typedef typename node<T, 2>::type parametric_node_type;
    using index_type = typename super::index_type;
    using size_type = typename super::size_type;
    using meshbase_type = MeshBase<index_type>;

    using marker_type = Marker<flag_type/*uint16_type*/>;
    /**
     * default constructor
     *
     */
    Geo0D() : Geo0D( invalid_v<index_type> ) {}

    /**
     * constructor where I give the id and declare if Geo0D object is on a
     * boundary
     *
     * @param id identifier of the Geo0D
     * @param boundary true if on the boundary, false otherwise
     * @param is_vertex true if the point is a vertex
     *
     * @return
     */
    explicit Geo0D( index_type id, bool boundary = false, bool is_vertex = false )
        :
        Geo0D( id, 0., 0., 0., boundary, is_vertex )
    {}

    /**
     * constructor where I give the id, the point coordinate and I declare
     * if the Geo0D object is on a boundary
     *
     * @param id identifier of the node
     * @param x x-coordinate of the node
     * @param y y-coordinate of the node
     * @param z z-coordinate of the node
     * @param boundary true if on the boundary, false otherwise
     * @param is_vertex true if the point is a vertex
     *
     */
    explicit Geo0D( index_type id, value_type x, value_type y, value_type z, bool boundary = false, bool is_vertex = false );

    /**
     * give the point coordinate
     *
     * @param x x-coordinate of the node
     * @param y y-coordinate of the node
     * @param z z-coordinate of the node
     */
    Geo0D( value_type x, value_type y, value_type z )
        :
        Geo0D( invalid_v<index_type>, x, y, z )
    {}

    /**
     * constructor where I give the id, the point coordinate and I declare
     * if the Geo0D object is on a boundary
     *
     * @param id identifier
     * @param __x node coordinate
     * @param boundary true if on the boundary, false otherwise
     * @param is_vertex true if the point is a vertex
     *
     */

    template <typename TheNodeType>
        Geo0D( index_type id, TheNodeType/*node_type*/ const& __x, bool boundary = false, bool is_vertex = false,  typename std::enable_if<!is_ptr_or_shared_ptr<TheNodeType>::value>::type* = nullptr )
        :
        Geo0D( id, new node_type( __x), boundary, is_vertex, false )
    {}

    /**
     * constructor where I give the id, the point coordinate and I declare
     * if the Geo0D object is on a boundary
     *
     * @param id identifier
     * @param __x node coordinate (shared_ptr)
     * @param boundary true if on the boundary, false otherwise
     * @param is_vertex true if the point is a vertex
     */

    Geo0D( index_type id, const node_type * __x, bool boundary = false, bool is_vertex = false, bool isView = false );

    /**
     * constructor where I give the id, the point coordinate and I declare
     * if the Geo0D object is on a boundary
     *
     * @param id identifier
     * @param __x another Geo0D object
     * @param boundary true if on the boundary, false otherwise
     * @param is_vertex true if the point is a vertex
     * @param isView true if the node is shared
     */
    Geo0D( index_type id, self_type const& __x, bool boundary = false, bool is_vertex = false, bool isView = false )
        :
        Geo0D( id, __x.M_node, boundary, is_vertex, isView )
    {}

    /**
     * the point coordinate
     *
     * @param __x node coordinate
     *
     */
    Geo0D( node_type const& __x )
        :
        Geo0D( invalid_v<index_type>, __x )
    {}

    /**
     * the point coordinate expression
     *
     * @param __expr expression for the coordinates
     *
     */
    template <typename AE>
    Geo0D( ublas::vector_expression<AE> const& __expr )
        :
        Geo0D( invalid_v<index_type>, new node_type(  __expr ) )
    {}

    /**
     * copy/move constructors
     */
    Geo0D( Geo0D const& G )
        :
        super( G ),
        M_isView( G.M_isView ),
        M_node( M_isView? G.M_node : new node_type( G.node() ) ),
        M_master_id( G.M_master_id ),
        M_master_vertex( G.M_master_vertex ),
        M_is_vertex( G.M_is_vertex ),
        M_mesh( G.M_mesh ), //??
        M_uv( G.M_uv )
        {}

    Geo0D( Geo0D&& G )
        :
        super( std::move( G ) ),
        M_isView( std::move( G.M_isView ) ),
        M_node( G.M_node ),
        M_master_id( std::move( G.M_master_id ) ),
        M_master_vertex( std::move( G.M_master_vertex ) ),
        M_is_vertex( std::move( G.M_is_vertex ) ),
        M_mesh( std::move( G.M_mesh ) ), //??
        M_uv( std::move( G.M_uv ) )
        {
            G.M_node = nullptr;
        }

    //! Destructor
    ~Geo0D() override
    {
        if ( !M_isView )
            delete M_node;
    }

     /**
     * assignement operators
     */
    Geo0D& operator=( Geo0D const& G )
    {
        if ( this != &G )
        {
            super::operator=( G );
            M_isView = G.M_isView;
            M_node = M_isView? G.M_node : new node_type( G.node() );
            M_master_id = G.M_master_id;
            M_master_vertex = G.M_master_vertex;
            M_is_vertex = G.M_is_vertex;
            M_mesh = G.M_mesh; //??
            M_uv = G.M_uv;
        }
        return *this;
    }
    Geo0D& operator=( Geo0D&& G )
    {
        if ( this != &G )
        {
            super::operator=( std::move( G ) );
            M_isView = std::move( G.M_isView );
            delete M_node;
            M_node = G.M_node;
            G.M_node = nullptr;
            M_master_id = std::move( G.M_master_id );
            M_master_vertex = std::move( G.M_master_vertex );
            M_is_vertex = std::move( G.M_is_vertex );
            M_mesh = std::move( G.M_mesh ); //??
            M_uv = std::move( G.M_uv );
        }
        return *this;
    }

    template <typename AE>
    Geo0D& operator=( ublas::vector_expression<AE> const& expr )
    {
        M_node->operator=( expr );
        return *this;
    }
    template <typename AE>
    Geo0D& operator+=( ublas::vector_expression<AE> const& expr )
    {
        M_node->operator+=( expr );
        return *this;
    }

    value_type& operator[]( int i )
    {
        return M_node->operator[]( i );
    }

    value_type operator[]( int i ) const
    {
        return M_node->operator[]( i );
    }

    value_type& operator()( int i )
    {
        return this->operator[]( i );
    }

    value_type operator()( int i ) const
    {
        return this->operator[]( i );
    }

    /**
     * the master id is different from the id in the context of periodic
     * boundary conditions on the slave side
     * @return the master id
     */
    index_type masterId() const
    {
        return M_master_id;
    }

    self_type const* masterVertex() const
    {
        return M_master_vertex;
    }

    /**
     * @return true if the entity is periodic, false otherwise
     */
    bool isPeriodic() const { return M_master_id != this->M_id; }

    /**
     * \return \p true if point is a vertex, \p false otherwise
     * \attention DO NOT USE YET, returns always false
     */
    bool isVertex() const
    {
        return M_is_vertex;
    }

    /**
     * set the point as a vertex or not using \p v
     */
    void setAsVertex( bool v )
    {
        M_is_vertex = v;
    }

    /**
     * set the mesh to which this geometric entity belongs to
     */
    void setMesh( meshbase_type const* m )
    {
        M_mesh = m;
    }

    /**
     * \return the mesh to which this geometric entity belongs to
     */
    meshbase_type const* mesh() const
    {
        return M_mesh;
    }

    /**
     * @return the node data structure
     */
    /*Geo0D*/node_type const& node() const
    {
        return *M_node;//*this;
    }

    /**
     * @return the node data structure
     */
    matrix_node_type G() const
    {
        matrix_node_type __G( Dim, 1 );
        ublas::column( __G, 0 ) = *M_node;//this;
        return __G;
    }

    /**
     * @return the node data structure
     */
    matrix_node_type vertices() const
    {
        return this->G();
    }

    /**
     * \return the measure of a point
     */
    value_type measure() const override
    {
        return 0;
    }

    /**
     * \return \p true if the point has parametric coordinates
     */
    bool isParametric() const
    {
        return M_uv.has_value();
    }

    /**
     * \return the geometric dimension of the entity the point belongs to
     */
    int gDim() const
    {
        return std::get<0>( M_uv.value() );
    }

    /**
     * \return the geometric dimension of the entity the point belongs to
     */
    int gTag() const
    {
        return std::get<1>( M_uv.value() );
    }

    /**
     * \return the parametric coordinates
     */
    parametric_node_type const& parametricCoordinates() const
    {
        return std::get<2>( M_uv.value() );
    }

    /**
     * \return the first parametric coordinates
     */
    value_type u() const
    {
        return this->parametricCoordinates()[0];
    }

    /**
     * \return the first parametric coordinates
     */
    value_type v() const
    {
        return this->parametricCoordinates()[1];
    }

    /**
     * set the node coordinates
     *
     * @param __n the node coordinates
     */
    void setNode( node_type const& __n )
    {
        *M_node = __n;
    }

    /**
     * \return true if points are equal, false otherwise
     */
    bool operator==( Geo0D const& geo0d ) const
    {
        return this->id() == geo0d.id();
    }

    bool operator<( Geo0D const& e ) const
    {
        return this->id() < e.id();
    }

    bool operator<( index_type __i ) const
    {
        return this->id() < __i;
    }

    /**
     * show the information about the Geo0D
     *
     * @param verbose true if verbose mode, false otherwise
     * @param c the output stream to use
     *
     * @return the output stream
     */
    std::ostream& showMe( bool verbose = false, std::ostream& c = std::cout ) const;

    void setMasterId( index_type id )
    {
        M_master_id = id;
    }
    void setMasterVertex( self_type const* m )
    {
        M_master_vertex = m;
    }

    /**
     * set the point coordinates. This will typically be called when
     * creating faces (points) in 1D.
     */
    void setPoint( uint16_type const /*i*/, self_type const& p )
    {
        *this = p;
    }
    void setPointCoordG( int i, ublas::vector<double> const& p )
    {
        *this = p;
    }

    /**
     * translate the point by \p trans
     */
    self_type& translate( node_type const& trans )
    {
        *this += trans;
        return *this;
    }

    /**
     * set the tags associated to the points
     * - tags[0] physical region
     * - tags[1] elementary region
     * - tags[2] particular region
     */
    void setTags( std::vector<int> const& tags ) override
    {
        this->setMarker( tags[0] );

        if ( tags.size() > 1 )
            this->setMarker2( tags[1] );

        if ( tags.size() > 2 )
            this->setProcessId( tags[2] );
    }
#if 1
    // TODO REMOVE THIS METHOD
    std::vector<int> tags() const
    {
        std::vector<int> thetags( 3 );
        thetags[0] = ( this->hasMarker( 1 ) ) ? this->marker( 1 ).value() : 0; //M_marker1.value();
        thetags[1] = ( this->hasMarker( 2 ) ) ? this->marker( 2 ).value() : 0; //M_marker2.value();
        thetags[2] = this->processId();
        return thetags;
    }
#endif
    /**
     * set the parametric coordinates of the node (if it is on an point, edge or
     * surface geometric entity)
     */
    void setParametricCoordinates( int gdim, int gtag, parametric_node_type const& x )
    {
        M_uv = std::make_tuple( gdim,gtag,x );
    }

    /**
     * set the parametric coordinates of the node (if it is on an point, edge or
     * surface geometric entity)
     */
    void setParametricCoordinates( int gdim, int gtag, value_type u, value_type v )
    {
        parametric_node_type uv(2);
        uv[0] = u;
        uv[1] = v;
        M_uv = std::make_tuple( gdim,gtag,uv );
    }

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
        ar& boost::serialization::base_object<super>( *this );
        //ar & M_is_vertex;
        ar & *(M_node);
        /*
            ar & M_gdim;
            ar & M_gtag;
            ar & M_uv;
            */
    }

  private:
    bool M_isView;
    node_type * M_node;

    index_type M_master_id;
    self_type const* M_master_vertex;

    bool M_is_vertex;

    // mesh to which the geond element belongs to
    meshbase_type const* M_mesh;

    std::optional<std::tuple<int,int,parametric_node_type>> M_uv;
};

// Alias for Geo0D<3>
typedef Geo0D<3> Point;

/*--------------------------------------------------------------
  Geo0D
  ---------------------------------------------------------------*/

template <uint16_type Dim, typename T, typename IndexT>
Geo0D<Dim, T, IndexT>::Geo0D( index_type id, value_type x, value_type y, value_type z, bool boundary, bool is_vertex )
    : super( id, MESH_ENTITY_INTERNAL ),
      M_isView( false ),
      M_node( new node_type( Dim ) ),
      M_master_id( id ),
      M_is_vertex( is_vertex ),
      M_mesh( nullptr )
{
    this->operator[]( 0 ) = x;

    if ( Dim >= 2 )
        this->operator[]( 1 ) = y;

    if ( Dim == 3 )
        this->operator[]( 2 ) = z;

    this->setOnBoundary( boundary );
}

template <uint16_type Dim, typename T, typename IndexT>
Geo0D<Dim, T, IndexT>::Geo0D( index_type id, const node_type * __p, bool boundary, bool is_vertex, bool isView )
    : super( id, MESH_ENTITY_INTERNAL ),
      M_isView( isView ),
      M_node( (isView)? const_cast<node_type*>(__p) : new node_type(*__p) ),
      M_master_id( id ),
      M_is_vertex( is_vertex ),
      M_mesh( nullptr )
{

    DCHECK( __p && __p->size() == Dim ) << "invalid node dimension : " << __p->size() << " and should be " << Dim;
    this->setOnBoundary( boundary );
}

template <uint16_type Dim, typename T, typename IndexT>
std::ostream&
Geo0D<Dim, T, IndexT>::showMe( bool /*verbose*/, std::ostream& out ) const
{
    out.setf( std::ios::scientific, std::ios::floatfield );
    out << "----- BEGIN of Geo0D ---\n";
    out << "id = " << this->id() << " node:" << this->node() << "\n";
    out << "is a vertex = " << isVertex() << "\n";
    out << "----- END OF Geo0D ---\n";
    return out;
}

template <uint16_type Dim, typename T, typename IndexT>
inline std::ostream&
operator<<( std::ostream& __os, Geo0D<Dim, T, IndexT> const& __n )
{
    __os.setf( std::ios::scientific, std::ios::floatfield );
    __os << "----- BEGIN of Geo0D ---\n";
    __os << "id = " << __n.id() << " node:" << (ublas::vector<double>&)__n.node() << "\n";
    __os << "is a vertex = " << __n.isVertex() << "\n";
    __os << "----- END OF Geo0D ---\n";
    return __os;
}

template <uint16_type Dim, typename T, typename IndexT>
inline FEELPP_DEPRECATED DebugStream&
operator<<( DebugStream& __os, Geo0D<Dim, T, IndexT> const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n.showMe( true, __str );

        __os << __str.str() << "\n";
    }

    return __os;
}
template <uint16_type Dim, typename T, typename IndexT>
inline FEELPP_DEPRECATED NdebugStream&
operator<<( NdebugStream& __os, Geo0D<Dim, T, IndexT> const& __n )
{
    return __os;
}

template <typename T, typename IndexT>
inline T
    distance( Geo0D<1, T, IndexT> const& p1, Geo0D<1, T, IndexT> const& p2 )
{
    return ublas::norm_2( p1.node() - p2.node() );
}

template <typename T, typename IndexT>
inline T
distance( Geo0D<2, T, IndexT> const& p1, Geo0D<2, T, IndexT> const& p2 )
{
    return ublas::norm_2( p1.node() - p2.node() );
}

template <typename T, typename IndexT>
inline T
distance( Geo0D<3, T, IndexT> const& p1, Geo0D<3, T, IndexT> const& p2 )
{
    return ublas::norm_2( p1.node() - p2.node() );
}
template <typename T, typename IndexT>
inline Geo0D<1, T, IndexT>
middle( Geo0D<1, T, IndexT> const& p1, Geo0D<1, T, IndexT> const& p2 )
{
    return ( p1.node() + p2.node() ) / 2;
}

template <typename T, typename IndexT>
inline Geo0D<2, T, IndexT>
middle( Geo0D<2, T, IndexT> const& p1, Geo0D<2, T, IndexT> const& p2 )
{
    return ( p1.node() + p2.node() ) / 2;
}

template <typename T, typename IndexT>
inline Geo0D<3, T, IndexT>
middle( Geo0D<3, T, IndexT> const& p1, Geo0D<3, T, IndexT> const& p2 )
{
    return ( p1.node() + p2.node() ) / 2;
}

template <typename E1, typename E2>
inline ublas::vector<double>
cross( ublas::vector_expression<E1> _p1,
       ublas::vector_expression<E2> _p2 )
{
    ublas::vector<double> v( 3 );
    ublas::vector<double> p1( _p1 );
    ublas::vector<double> p2( _p2 );
    v( 0 ) = p1( 1 ) * p2( 2 ) - p1( 2 ) * p2( 1 );
    v( 1 ) = p1( 2 ) * p2( 0 ) - p1( 0 ) * p2( 2 );
    v( 2 ) = p1( 0 ) * p2( 1 ) - p1( 1 ) * p2( 0 );
    return v;
}

template <typename T, typename IndexT>
inline ublas::vector<double>
cross( Geo0D<3, T, IndexT> p1,
       Geo0D<3, T, IndexT> p2 )
{
    ublas::vector<double> v( 3 );
    v( 0 ) = p1( 1 ) * p2( 2 ) - p1( 2 ) * p2( 1 );
    v( 1 ) = p1( 2 ) * p2( 0 ) - p1( 0 ) * p2( 2 );
    v( 2 ) = p1( 0 ) * p2( 1 ) - p1( 1 ) * p2( 0 );
    return v;
}

} // namespace Feel

#endif
