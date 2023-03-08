/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
   Copyright (C) 2006-2010 Université Joseph Fourier (Grenoble I)

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
/*!
  \file geoElement.hpp
  \brief Geometric elements
  Introduces all the geometric elements
*/

#ifndef _GEOELEMENT_HH_
#define _GEOELEMENT_HH_

#include <boost/tuple/tuple.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/matrix.hpp>
#include <feel/feelmesh/marker.hpp>
#include <feel/feelmesh/geond.hpp>

//#include <feel/feelalg/lu.hpp>

namespace Feel
{
template<int Dim = 1, typename IndexT = uint32_type>
class SubFaceOfNone
{
public:
    static const uint16_type nDim = Dim;
    using index_type = IndexT;
    using size_type = index_type;
    template<typename ET>
    struct Element
    {
        typedef ET type;
    };
    typedef boost::tuple<size_type, size_type, uint16_type> element_connectivity_type;

    boost::none_t element( uint16_type /* e */ ) const
    {
        return boost::none;
    }

    SubFaceOfNone() {}
    SubFaceOfNone( SubFaceOfNone const& ) {}
    SubFaceOfNone( SubFaceOfNone && ) = default;
    SubFaceOfNone& operator=( SubFaceOfNone const& ) = default;
    SubFaceOfNone& operator=( SubFaceOfNone && ) = default;

    virtual ~SubFaceOfNone() {}

    template<typename SFO>
    SubFaceOfNone( SFO const& /*sf*/ )
    {
    }
    bool
    isGhostFace( rank_type /*p*/ ) const
    {
        return false;
    }
    bool
    isInterProcessDomain( rank_type /*p*/ ) const
    {
        return false;
    }
    rank_type partition1( rank_type /*p*/ ) const { return invalid_rank_type_value; }
    rank_type partition2( rank_type /*p*/ ) const { return invalid_rank_type_value; }
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
        }
};

/**
 * @brief Tag base class
 * @ingroup Mesh
 */
struct SubFaceOfBase {};

/**
 * @brief description of a subface or facet (topological d-1) of an element of topologicql dimension d
 * @ingroup Mesh
 * 
 * @tparam ElementType type of the element 
 */
template<typename ElementType>
class SubFaceOf : public SubFaceOfBase
{
public:
    static const uint16_type nDim = ElementType::nDim;
    using index_type = typename ElementType::index_type;
    using size_type = typename ElementType::size_type;
    template<typename ET>
    struct Element
    {
        typedef ElementType type;
    };
    typedef ElementType entity_type;
    typedef boost::tuple<ElementType const*, uint16_type> element_connectivity_type;

    SubFaceOf()
        :
        M_element0( 0, invalid_uint16_type_value ),
        M_element1( 0, invalid_uint16_type_value )
    {}

    SubFaceOf( element_connectivity_type const& connect0 )
        :
        M_element0( connect0 ),
        M_element1( 0, invalid_uint16_type_value )
    {}
    SubFaceOf( element_connectivity_type const& connect0,
               element_connectivity_type const& connect1 )
        :
        M_element0( connect0 ),
        M_element1( connect1 )
    {}

    SubFaceOf( SubFaceOf const& sf )
        :
        M_element0( sf.M_element0 ),
        M_element1( sf.M_element1 )
    {
    }
    SubFaceOf( SubFaceOf && sf )
        :
        M_element0( std::move( sf.M_element0 ) ),
        M_element1( std::move( sf.M_element1 ) )
    {
     }

        template <typename TheEltType >
        /*explicit*/ SubFaceOf( SubFaceOf<TheEltType> const& /*sf*/,
                                std::enable_if_t< !std::is_same<TheEltType, entity_type>::value >* = nullptr )
        :
        M_element0( 0, invalid_uint16_type_value ),
        M_element1( 0, invalid_uint16_type_value )
    {
    }

    template <int SFoD >
    explicit SubFaceOf( SubFaceOfNone<SFoD> const& /*sf*/,
                        std::enable_if_t<SFoD == nDim || SFoD == 0>* = nullptr )
        :
        M_element0( 0, invalid_uint16_type_value ),
        M_element1( 0, invalid_uint16_type_value )
    {
    }
    virtual ~SubFaceOf() {}

    SubFaceOf& operator=( SubFaceOf const& sf ) = default;

    SubFaceOf& operator=( SubFaceOf && sf )
        {
            M_element0 = std::move(sf.M_element0);
            M_element1 = std::move(sf.M_element1);
            ///std::cout << "move assigned SubFaceOf\n";
            return *this;
        }

    entity_type const& element( uint16_type e ) const
    {
        if ( e == 0 )
            return *boost::get<0>( M_element0 );

        else
            return *boost::get<0>( M_element1 );
    }

    entity_type const& element0() const
    {
        return *boost::get<0>( M_element0 );
    }
    entity_type const& element1() const
    {
        return *boost::get<0>( M_element1 );
    }
    size_type idElement0() const { return this->element0().id(); }
    uint16_type idInElement0() const { return boost::get<1>( M_element0 ); }
    rank_type pidElement0() const { return this->element0().processId(); }

    size_type idElement1() const { return this->element1().id(); }
    uint16_type idInElement1() const { return boost::get<1>( M_element1 ); }
    rank_type pidElement1() const { return this->element1().processId(); }

    size_type ad_first() const
    {
        return this->element0().id();
    }
    uint16_type pos_first() const
    {
        return boost::get<1>( M_element0 );
    }
    rank_type proc_first() const
    {
        return this->element0().processId();
    }
    rank_type partition1( rank_type p ) const
    {
        return p;
    }

    size_type ad_second() const
    {
        return this->element1().id();
    }
    uint16_type pos_second() const
    {
        return boost::get<1>( M_element1 );
    }
    rank_type proc_second() const
    {
        return this->element1().processId();
    }
    rank_type partition2( size_type p ) const
    {
        return ( p == this->proc_first() )? this->proc_second() : this->proc_first();
    }

    element_connectivity_type const& connection0() const
    {
        return M_element0;
    }
    element_connectivity_type const& connection1() const
    {
        return M_element1;
    }

    void setConnection( uint16_type f, element_connectivity_type const& connect )
    {
        if ( f == 0 )
            M_element0 = connect;

        else
            M_element1 = connect;

    }

    void setConnection0( element_connectivity_type const& connect )
    {
        M_element0 = connect;
    }
    void setConnection1( element_connectivity_type const& connect )
    {
        M_element1 = connect;
    }

    bool isConnected() const { return isConnectedTo0() && isConnectedTo1(); }

    bool isConnectedTo0() const
    {
        return ( boost::get<0>( M_element0 ) != 0 );
    }
    bool isConnectedTo1() const
    {
        return ( boost::get<0>( M_element1 ) != 0 );
    }

    /**
     * @brief say if the face is a ghost or not for process id @p p
     * 
     * a face is a ghost face on process rank @p p if the process id is greater than @p p
     * 
     * @param p the rank of the communicator
     * @return true if the face is a ghost face 
     * @return false otherwise
     */
    bool
    isGhostFace( rank_type p ) const
    {
        return ( this->isConnectedTo0() && this->isConnectedTo1() &&
                 ( ( ( this->pidElement0() == p ) && ( this->pidElement1() < p ) ) ||
                   ( ( this->pidElement0() < p ) && ( this->pidElement1() == p ) ) ) );
    }

    bool
    isInterProcessDomain( rank_type p ) const
    {
        return ( this->isConnectedTo0() && this->isConnectedTo1() &&
                 ( ( this->pidElement0() == p ) || ( this->pidElement1() == p ) ) &&
                 ( this->pidElement0() != this->pidElement1() ) );
    }
    bool
    isIntraProcessDomain( rank_type p ) const
    {
        bool hasConnect0 = this->isConnectedTo0();
        bool hasConnect1 = this->isConnectedTo1();
        if ( hasConnect0 && hasConnect1 )
            return ( ( this->pidElement0() == p ) && ( this->pidElement1() == p ) );
        else if ( hasConnect0 && !hasConnect1 )
            return ( this->pidElement0() == p );
        else if ( !hasConnect0 && hasConnect1 )
            return ( this->pidElement1() == p );
        return false;
    }

    void disconnect0()
    {
        M_element0 = boost::make_tuple( ( ElementType const* )0, invalid_uint16_type_value );
    }

    void disconnect1()
    {
        M_element1 = boost::make_tuple( ( ElementType const* )0, invalid_uint16_type_value );
    }

    void disconnect()
    {
        disconnect0();
        disconnect1();
    }

    void disconnect( ElementType const& elem )
    {
        if(boost::get<0>( M_element0 ) == boost::addressof(elem))
        {
            DVLOG(2) << "connecting 1 to 0 and disconnecting 1..\n";
            M_element0 = M_element1;
            disconnect1();
        }
        else
        {
            DVLOG(2) << "disconnecting 1..\n";
            disconnect1();
        }
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
#if 1
            ar & M_element0.template get<1>();
            M_element0.template get<0>() = 0;

            ar & M_element1.template get<1>();
            M_element1.template get<0>() = 0;
#endif
        }
private:

    element_connectivity_type M_element0;
    element_connectivity_type M_element1;

};
#if 0
template<typename ElementType>
class SubFaceOfMany
{
public:
    static const uint16_type nDim = ElementType::nDim;
    static const uint16_type nRealDim = ElementType::nRealDim;
    template<typename ET>
    struct Element
    {
        typedef ElementType type;
    };
    typedef ElementType entity_type;
    typedef boost::tuple<ElementType const*, size_type, uint16_type, size_type> element_connectivity_type;

    SubFaceOfMany()
        :
        M_elements()
        {}

    SubFaceOfMany( SubFaceOfMany const& sf )
        :
        M_elements( sf.M_elements )
        {
        }
    SubFaceOfMany( SubFaceOfNone<nDim> const& /*sf*/ )
        :
        M_elements()
        {
        }
    virtual ~SubFaceOfMany() {}

    SubFaceOfMany& operator=( SubFaceOfMany const& sf )
        {
            if ( this != &sf )
            {
                M_elements = sf.M_elements;
            }

            return *this;
        }
    SubFaceOfMany& operator=( SubFaceOfNone<nDim> const& /*sf*/ )
    {
        return *this;
    }

    entity_type const& element0() const
    {
        return *boost::get<0>( *M_elements.begin() );
    }
    entity_type const& element1() const
    {
        return *boost::get<0>( *boost::next(M_elements.begin()) );
    }

    void setConnection( element_connectivity_type const& connect )
    {
        this->insert( connect );
    }

    size_type ad_first() const
    {
        return boost::get<1>( *M_elements.begin() );
    }
    uint16_type pos_first() const
    {
        return boost::get<2>( *M_elements.begin() );
    }
    size_type proc_first() const
    {
        return boost::get<3>( *M_elements.begin() );
    }

    size_type ad_second() const
    {
        return boost::get<1>( *boost::next(M_elements.begin()) );
    }
    uint16_type pos_second() const
    {
        return boost::get<2>( *boost::next(M_elements.begin()) );
    }
    size_type proc_second() const
    {
        return boost::get<3>( *boost::next(M_elements.begin()) );
    }


    void setConnection0( element_connectivity_type const& connect )
    {
        M_elements.insert( connect );
    }
    void setConnection1( element_connectivity_type const& connect )
    {
        M_elements.insert( connect );
    }

    element_connectivity_type const& connection0() const
    {
        return *M_elements.begin();
    }
    element_connectivity_type const& connection1() const
    {
        return *boost::next(M_elements.begin());
    }

    bool isConnectedTo0() const
    {
        return ( boost::get<1>( *M_elements.begin() ) != invalid_v<size_type> &&
                 boost::get<2>( *M_elements.begin() ) != invalid_uint16_type_value &&
                 boost::get<3>( *M_elements.begin() ) != invalid_v<size_type> );
    }
    bool isConnectedTo1() const
    {
        return ( boost::get<1>( *boost::next(M_elements.begin()) ) != invalid_v<size_type> &&
                 boost::get<2>( *boost::next(M_elements.begin()) ) != invalid_uint16_type_value &&
                 boost::get<3>( *boost::next(M_elements.begin()) ) != invalid_v<size_type> );
    }
    bool
    isInterProcessDomain( size_type p ) const
    {
        return false;
    }
    bool
    isIntraProcessDomain( size_type p ) const
    {
        return true;
    }

    entity_type const& element( uint16_type e ) const
    {
        return *boost::get<0>( *M_elements.begin() );
    }
    void disconnect()
    {
        M_elements.clear();
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
        }
private:

    std::set<element_connectivity_type> M_elements;
};
#else
//template<typename ElementType> class SubFaceOfMany: public SubFaceOf<ElementType> {};
#define SubFaceOfMany SubFaceOf

#endif


/**
 * Class for Points and Vertices
 * @ingroup Mesh
 */
template <uint16_type Dim,
          typename SubFace = SubFaceOfNone<Dim>,
          typename T = double,
          typename IndexT = uint32_type>
class GeoElement0D
    :
    public Geo0D<Dim,T,IndexT>,
    public SubFace
{
public:

    static inline const uint16_type nDim = 0;
    static inline const uint16_type nRealDim = Dim;
    static inline const bool is_simplex = true;

    typedef Geo0D<Dim,T,IndexT> geo0d_type;
    using index_type = typename geo0d_type::index_type;
    using size_type = typename geo0d_type::size_type;
    typedef typename geo0d_type::node_type node_type;

    typedef geo0d_type super;
    typedef SubFace super2;

    typedef GeoElement0D<Dim,SubFace,T,IndexT> self_type;
#if 0
    using element_type = typename mpl::if_<mpl::bool_<SubFace::nDim==0>,
                                           mpl::identity<self_type>,
                                           mpl::identity<typename SubFace::template Element<self_type>::type> >::type::type ;
#else
    using element_type = self_type;
    using gm_type = GT_Lagrange<0, 1, nRealDim, Simplex, T>;
    using gm1_type = gm_type;
#endif
    typedef self_type point_type;

    typedef typename super::matrix_node_type matrix_node_type;

    static inline const uint16_type numLocalVertices = super::numVertices;

    GeoElement0D() = default;

    //! Declares item id and if it is on boundary
    GeoElement0D( size_type id, bool boundary = false )
        :
        super( id, boundary ),
        super2()
    {}

    template <typename GeoNodeType>
    GeoElement0D( size_type id, GeoNodeType const& n,  bool boundary = false )
        :
        super( id, n, boundary ),
        super2()
        {}
    GeoElement0D( size_type id, geo0d_type const& n, bool boundary, bool isView )
        :
        super( id, n, boundary, false, isView ),
        super2()
    {}

    //! Declares item id and if it is on boundary, and provides coordinate
    //! data.
    GeoElement0D( size_type id, Real x, Real y, Real z, bool boundary = false )
        :
        super( id, x, y, z, boundary ),
        super2()
    {}

    explicit GeoElement0D( geo0d_type const& p )
        :
        super( p ),
        super2()
        {}
    GeoElement0D( GeoElement0D const & g ) = default;
    GeoElement0D( GeoElement0D && g ) = default;

    template<typename SF>
    GeoElement0D( GeoElement0D<Dim,SF,T> const & g )
        :
        super( g ),
        super2( g )
    {
    }

    ~GeoElement0D() override
    {}

    GeoElement0D & operator = ( GeoElement0D const& g ) = default;
    GeoElement0D & operator = ( GeoElement0D && g ) = default;

    template<typename SF>
    GeoElement0D & operator = ( GeoElement0D<Dim,SF,T> const & g )
    {
        super::operator=( g );
        super2::operator=( g );
        return *this;
    }

    //void setMesh( MeshBase<> const* m ) { super::setMesh( m ); }
    MeshBase<> const* mesh() const
    {
        return super::mesh();
    }

    /**
     * \return id
     */
    size_type id() const
    {
        return super::id();
    }

    /**
     * \return process id
     */
    rank_type processId() const
    {
        return super::processId();
    }

    /**
     * \return process id
     */
    rank_type partition1() const
    {
        return super2::partition1( super::processId() );
    }

    /**
     * \return process id
     */
    rank_type partition2() const
    {
        return super2::partition2( super::processId() );
    }

    bool isGhostFace() const
    {
        return super2::isGhostFace( super::processId()  );
    }

    /**
     * \return \p true if interprocess domain face, \p false otherwise
     */
    bool isInterProcessDomain() const
    {
        return super2::isInterProcessDomain( super::processId()  );
    }

    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    bool isOnBoundary() const
    {
        return super::isOnBoundary();
    }

    /**
     * \return maximum \c dimension of the sub-entity touching the boundary of the element
     */
    uint16_type boundaryEntityDimension() const
    {
        return 0;
    }

    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const
    {
        return super::isGhostCell();
    }

    /**
     * \return the point associated to the face
     */
    geo0d_type const& point( uint16_type /*i*/ ) const
    {
        return *this;
    }

    /**
     * \return the point associated to the face
     */
    geo0d_type & point( uint16_type /*i*/ )
    {
        return *this;
    }

    /**
     * set the geometrical point associated to the face
     */
    void setPoint( uint16_type /*i*/, geo0d_type & e )
    {
        //M_facept = e;
        //M_facept= const_cast<geo0d_type *>( &e );
        //M_facept= std::addressof( e );
    }

    matrix_node_type /*const&*/ G() const
    {
        return super::G();
    }

    matrix_node_type /*const&*/ vertices() const
    {
        return super::vertices();
    }

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & boost::serialization::base_object<super>( *this );
            ar & boost::serialization::base_object<super2>( *this );
        }


};


/**
 * \class GeoElement1D
 * @ingroup Mesh
 * \brief class for 1D elements
 *
 * In the 2D case, we store the size_types of the adjacent 2D elements
 * and their relative position.
 */
template<uint16_type Dim,
         typename GEOSHAPE,
         typename SubFace = SubFaceOfNone<0>,
         typename T = double,
         typename IndexT = uint32_type,
         bool PointTypeIsSubFaceOf = false,
         bool UseMeasuresStorage = false >
class GeoElement1D
    :
    public GeoND<Dim, GEOSHAPE, T, IndexT,
                 typename mpl::if_<mpl::bool_<PointTypeIsSubFaceOf>,
                                   mpl::identity< GeoElement0D<Dim, SubFaceOf<GeoElement1D<Dim, GEOSHAPE, SubFace, T, IndexT, PointTypeIsSubFaceOf, UseMeasuresStorage> >, T, IndexT> >,
                                   mpl::identity< GeoElement0D<Dim, SubFaceOfNone<0>, T, IndexT> > >::type::type,
                 UseMeasuresStorage
                 >,
    public SubFace
{
public:

    typedef GeoND<Dim, GEOSHAPE, T, IndexT,
                  typename mpl::if_<mpl::bool_<PointTypeIsSubFaceOf>,
                                    mpl::identity< GeoElement0D<Dim, SubFaceOf<GeoElement1D<Dim, GEOSHAPE, SubFace, T, IndexT, PointTypeIsSubFaceOf, UseMeasuresStorage> >, T, IndexT> >,
                                    mpl::identity< GeoElement0D<Dim, SubFaceOfNone<0>, T, IndexT> > >::type::type,
                  UseMeasuresStorage
                  > super;

    typedef SubFace super2;

    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;
    static inline const uint16_type nRealDim = super::nRealDim;

    static inline const bool condition = ( Dim==nRealDim );
    BOOST_MPL_ASSERT_MSG( ( condition ), INVALID_ELEMENT_REAL_DIMENSION, ( mpl::int_<Dim>, mpl::int_<nRealDim>, GEOSHAPE ) );

    using index_type = typename super::index_type;
    using size_type = typename super::size_type;
    typedef GEOSHAPE GeoShape;
    typedef GeoElement1D<Dim, GEOSHAPE, SubFace, T, IndexT, PointTypeIsSubFaceOf, UseMeasuresStorage> self_type;
    //typedef typename SubFace::template Element<self_type>::type element_type;
    typedef self_type element_type;
#if 0
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nRealDim>,mpl::int_<1> >,
                              mpl::identity<GeoElement0D<Dim, SubFaceOf<self_type>, T, IndexT> >,
                              mpl::identity<GeoElement0D<Dim, SubFaceOfMany<self_type>, T, IndexT> > >::type::type point_type;
    typedef point_type GeoBElement;
#endif
    typedef typename super::point_type point_type;

    static inline const uint16_type numLocalVertices = super::numLocalVertices;
    static inline const uint16_type numLocalEdges = super::numEdges;
    static inline const uint16_type numLocalFaces = super::numLocalVertices;


    typedef typename super::node_type node_type;
    typedef typename super::vertex_permutation_type vertex_permutation_type;
    typedef typename super::edge_permutation_type edge_permutation_type;
    typedef typename super::face_permutation_type face_permutation_type;
    typedef typename super::permutation_type permutation_type;

    /**
     * default constructor, make it explicit to avoid implicit
     * inversion to \c size_type
     */
    explicit GeoElement1D( size_type id = 0 )
        :
        super( id ),
        super2(),
        M_vertices()
        //M_vertex_permutation( numLocalVertices )
    {
        M_vertices.fill( nullptr );
    }
    /**
     * copy constructor
     */
    GeoElement1D( GeoElement1D const& g ) = default;
    GeoElement1D( GeoElement1D && g )
        :
        super( std::move(g) ),
        super2( std::move(g) ),
        M_map( std::move( g.M_map ) ),
        M_vertices( std::move( g.M_vertices ) )
        {
            //std::cout << "GeoElement1D move ctor\n";
        }
           
        

    /**
     * destructor
     */
    ~GeoElement1D() override
    {}

    /**
     * copy operator
     */
    GeoElement1D& operator=( GeoElement1D const& g ) = default;
    GeoElement1D& operator=( GeoElement1D && g )
        {
            super::operator=( std::move(g) );
            super2::operator=( std::move(g) );
            M_map = std::move( g.M_map );
            M_vertices = std::move( g.M_vertices );
            //std::cout << "GeoElement1D move op\n";
            return *this;
        }


    //void setMesh( MeshBase<> const* m ) { super::setMesh( m ); }
    MeshBase<> const* mesh() const
    {
        return super::mesh();
    }


    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    size_type id() const
    {
        return super::id();
    }

    bool isGhostFace() const
    {
        return super2::isGhostFace( super::processId()  );
    }


    /**
     * \return \p true if interprocess domain face, \p false otherwise
     */
    bool isInterProcessDomain() const
    {
        return super2::isInterProcessDomain( super::processId() );
    }

    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    bool isOnBoundary() const
    {
        return super::isOnBoundary();
    }

    /**
     * \return maximum \c dimension of the sub-entity touching the boundary of the element
     */
    uint16_type boundaryEntityDimension() const
    {
        return 0;
    }

    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const
    {
        return super::isGhostCell();
    }

    /**
     * \return process id
     */
    rank_type processId() const
    {
        return super::processId();
    }

    /**
     * \return process id
     */
    rank_type partition1() const
    {
        return super2::partition1( super::processId() );
    }

    /**
     * \return process id
     */
    rank_type partition2() const
    {
        return super2::partition2( super::processId() );
    }

    void setMap( uint8_type k_1, uint8_type k_2 )
    {
        M_map[k_1] = k_2;
    }

    uint8_type map( uint8_type k_1 ) const
    {
        return M_map[ k_1 ];
    }

    /**
     * Inserts a point as face of the edge geometric element
     */
    void setFace( uint16_type const i, point_type const & p )
    {
        FEELPP_ASSERT( i < numLocalVertices )( i ).error( "invalid local point index" );
        M_vertices[i] = const_cast<point_type*>( boost::addressof( p ) );
    }

    edge_permutation_type permutation( uint16_type /*i*/ ) const
    {
        return edge_permutation_type();
    }
    //!
    //! @return true if GeoElement1D is connected to a face 
    //!
    bool hasFace( uint16_type i ) const
        {
            if ( i >= numLocalFaces )
                return false;
            return M_vertices[i] != nullptr;
        }
    point_type const& edge( uint16_type i ) const
    {
        return *M_vertices[i];
    }
    point_type const& face( uint16_type i ) const
    {
        return *M_vertices[i];
    }
    point_type & face( uint16_type i )
    {
        return *M_vertices[i];
    }
    point_type const* facePtr( uint16_type i ) const
    {
        FEELPP_ASSERT( i < numLocalVertices )( this->id() )( i ).error( "invalid local vertex index" );
        return M_vertices[i];
    }
    point_type* facePtr( uint16_type i )
    {
        FEELPP_ASSERT( i < numLocalVertices )( this->id() )( i ).error( "invalid local vertex index" );
        return M_vertices[i];
    }

    typedef typename std::array<point_type*, numLocalVertices>::iterator face_iterator;
    typedef typename std::array<point_type*, numLocalVertices>::const_iterator face_const_iterator;

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_iterator,face_iterator>
    faces()
    {
        return std::make_pair( M_vertices.begin(), M_vertices.end() );
    }

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_const_iterator,face_const_iterator>
    faces() const
    {
        return std::make_pair( M_vertices.begin(), M_vertices.end() );
    }
    std::vector<size_type> facesId() const
        {
            std::vector<size_type> fid;
            std::for_each( M_vertices.begin(), M_vertices.end(),
                           [&fid]( auto const& f ) { fid.push_back(f->id()); } );
            return fid;
        }
    /**
     * \sa edgePermutation(), permutation()
     */
    vertex_permutation_type facePermutation( uint16_type /*i*/ ) const
    {
        return edge_permutation_type();

        //FEELPP_ASSERT( i < numLocalVertices )( i )( numLocalVertices ).error( "invalid local vertex index" );
        //return M_vertex_permutation[i];
    }
    vertex_permutation_type permutation( uint16_type, mpl::int_<1> ) const override
    {
        return vertex_permutation_type();
    }

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            DVLOG(2) << "Serializing Geoelement1D id: " << this->id() << "...\n";
            ar & boost::serialization::base_object<super>( *this );
            ar & boost::serialization::base_object<super2>( *this );
        }



private:

    std::vector<uint8_type> M_map;
    std::array<point_type*, numLocalVertices> M_vertices;
    //ublas::bounded_array<vertex_permutation_type, numLocalVertices> M_vertex_permutation;

};

/**
 * \class GeoElement2D
 * @ingroup Mesh
 * \brief  Class for 2D elements.
 *
 * In the 3D case, we store the size_types of the adjacent 3D elements
 * and their relative position.
 */
template<uint16_type Dim,
         typename GEOSHAPE,
         typename SubFace = SubFaceOfNone<0>,
         typename T = double,
         typename IndexT = uint32_type,
         bool UseMeasuresStorage = false >
class GeoElement2D
    :
        public GeoND<Dim, GEOSHAPE, T, IndexT, GeoElement0D<Dim, SubFaceOfNone<0>, T, IndexT>, UseMeasuresStorage >,
public SubFace
{
public:


    typedef GeoND<Dim, GEOSHAPE, T, IndexT, GeoElement0D<Dim, SubFaceOfNone<0>, T, IndexT>, UseMeasuresStorage > super;
    typedef SubFace super2;

    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;
    static inline const uint16_type nRealDim = super::nRealDim;

    static inline const bool condition = ( Dim==nRealDim );
    BOOST_MPL_ASSERT_MSG( ( condition ), INVALID_ELEMENT_REAL_DIMENSION, ( mpl::int_<Dim>, mpl::int_<nRealDim>, GEOSHAPE ) );

    //! Number of element edges
    static inline const uint16_type numLocalEdges = super::numEdges;
    static inline const uint16_type numLocalFaces = super::numFaces;

    using index_type = typename super::index_type;
    using size_type = typename super::size_type;
    typedef GEOSHAPE GeoShape;
    typedef typename super::face_type entity_face_type;
    typedef GeoElement2D<Dim, GEOSHAPE,SubFace, T, IndexT, UseMeasuresStorage> self_type;
    //typedef typename SubFace::template Element<self_type>::type element_type;
    typedef self_type element_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nRealDim>,mpl::int_<2> >,
                              mpl::identity<GeoElement1D<Dim, entity_face_type, SubFaceOf<self_type>, T, IndexT> >,
                              mpl::identity<GeoElement1D<Dim, entity_face_type, SubFaceOfMany<self_type>, T, IndexT> > >::type::type edge_type;
    //typedef GeoElement1D<Dim, entity_face_type, SubFaceOf<self_type>, T > edge_type;
    typedef GeoElement0D<Dim, SubFaceOfNone<0>, T, IndexT> point_type;
#if 0
    BOOST_MPL_ASSERT_MSG( ( boost::is_same<point_type,typename edge_type::point_type>::value ),
                          INCOMPATIBLE_POINT_TYPE,
                          ( point_type, typename edge_type::point_type, edge_type, element_type, self_type ) );
    BOOST_STATIC_ASSERT( ( boost::is_same<point_type,typename edge_type::point_type>::value ) );
#endif // 0
    typedef typename super::node_type node_type;

    typedef typename super::vertex_permutation_type vertex_permutation_type;
    typedef typename super::edge_permutation_type edge_permutation_type;
    typedef typename super::edge_permutation_type permutation_type;
    //typedef typename super::face_permutation_type face_permutation_type;
    using face_permutation_type = edge_permutation_type;
    typedef typename super2::element_connectivity_type element_connectivity_type;

    /**
     * default constructor, make it explicit to avoid implicit
     * inversion to \c size_type
     */
    explicit GeoElement2D( size_type id = 0 )
        :
        super( id ),
        super2(),
        M_edges( numLocalEdges, nullptr ),
        M_edge_permutation( numLocalEdges, edge_permutation_type( edge_permutation_type::IDENTITY ) )
    {
    }

    /**
     * copy constructor
     */
    GeoElement2D( GeoElement2D const& g ) = default;
    GeoElement2D( GeoElement2D && g ) = default;

    /**
     * destructor
     */
    ~GeoElement2D() override
    {}

    /**
     * copy operator
     */
    GeoElement2D& operator=( GeoElement2D const& g ) = default;
    GeoElement2D& operator=( GeoElement2D && g ) = default;

    //void setMesh( MeshBase<> const* m ) { super::setMesh( m ); }
    MeshBase<> const* mesh() const
    {
        return super::mesh();
    }
    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    size_type id() const
    {
        return super::id();
    }

    bool isGhostFace() const
    {
        return super2::isGhostFace( super::processId()  );
    }


    /**
     * \return \p true if interprocess domain face, \p false otherwise
     */
    bool isInterProcessDomain() const
    {
        return super2::isInterProcessDomain( super::processId() );
    }

    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    bool isOnBoundary() const
    {
        return super::isOnBoundary();
    }

    /**
     * \return maximum \c dimension of the sub-entity touching the boundary of the element
     */
    uint16_type boundaryEntityDimension() const
    {
        return super::boundaryEntityDimension();
    }

    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const
    {
        return super::isGhostCell();
    }

    /**
     * \return process id
     */
    rank_type processId() const
    {
        return super::processId();
    }

    /**
     * \return process id
     */
    rank_type partition1() const
    {
        return super2::partition1( super::processId() );
    }

    /**
     * \return process id
     */
    rank_type partition2() const
    {
        return super2::partition2( super::processId() );
    }


    /**
     * \return true if element have an edge connected
     */
    bool hasEdge( uint16_type i ) const
    {
        if ( i >= numLocalEdges )
            return false;
        return M_edges[i] != nullptr;
    }

    /**
     * \return true if element have a face connected
     */
    bool hasFace( uint16_type i ) const
    {
        return this->hasEdge( i );
    }

    /**
     * \sa face()
     */
    edge_type const& edge( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();
        return boost::cref( *M_edges[i] );
    }

    /**
     * \sa face()
     */
    edge_type& edge( uint16_type i )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();
        return boost::ref( *M_edges[i] );
    }

    edge_type & face( uint16_type i )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return boost::ref( *M_edges[i] );
    }

    /**
     * \sa edge()
     */
    edge_type const& face( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return boost::cref( *M_edges[i] );
    }

    edge_type const* facePtr( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        return M_edges[i];
    }

    /**
     * Inserts an edge.
     * \sa setEdge()
     */
    void setFace( uint16_type const i, edge_type const & p )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;

        M_edges[i] = const_cast<edge_type*>( boost::addressof( p ) );
    }

    /**
     * \sa facePermutation(), permutation()
     */
    edge_permutation_type edgePermutation( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;

        return M_edge_permutation[i];
    }
    /**
     * \sa edgePermutation(), permutation()
     */
    edge_permutation_type facePermutation( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;

        return M_edge_permutation[i];
    }

    /**
     * \sa edgePermutation(), facePermutation()
     */
    edge_permutation_type permutation( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges << " in element id " << this->id();
        return M_edge_permutation[i];
    }

    edge_permutation_type permutation( uint16_type i, mpl::int_<1> ) const override
    {
        return this->edgePermutation( i );
    }
    vertex_permutation_type permutation( uint16_type, mpl::int_<2> ) const
    {
        return vertex_permutation_type();
    }

    /**
     * Inserts a point.  Uses point references put point
     * \sa setFace()
     */
    void setEdge( uint16_type i, edge_type const & p )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        M_edges[i] = const_cast<edge_type*>( boost::addressof( p ) );
    }

    void setEdgePermutation( uint16_type i, edge_permutation_type o )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;

        M_edge_permutation[i] = o;
    }
    void setFacePermutation( uint16_type i, edge_permutation_type o ) { this->setEdgePermutation( i,o ); }

    typedef typename std::vector<edge_type*>::iterator face_iterator;
    typedef typename std::vector<edge_type*>::const_iterator face_const_iterator;

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_iterator,face_iterator>
    faces()
    {
        return std::make_pair( M_edges.begin(), M_edges.end() );
    }

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_const_iterator,face_const_iterator>
    faces() const
    {
        return std::make_pair( M_edges.begin(), M_edges.end() );
    }

    std::vector<size_type> facesId() const
        {
            std::vector<size_type> fid;
            std::for_each( M_edges.begin(), M_edges.end(),
                           [&fid]( auto const& f ) { fid.push_back(f->id()); } );
            return fid;
        }
    void disconnectSubEntities()
    {
        for(unsigned int i = 0; i<numLocalEdges;++i)
        {
            M_edges[i]->disconnect(*this);
        }
    }

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            DVLOG(2) << "Serializing Geoelement2D id: " << this->id() << "...\n";
            ar & boost::serialization::base_object<super>( *this );
#if 0
            ar & boost::serialization::base_object<super2>( *this );
            ar & M_edges;
#endif
        }


private:

    std::vector<edge_type*> M_edges;
    std::vector<edge_permutation_type> M_edge_permutation;
};

/*-------------------------------------------------------------------------
  GeoElement2D
  --------------------------------------------------------------------------*/



/**
 * \class GeoElement3D
 * @ingroup Mesh
 * \brief Class for 3D elements
 *
 */
template<uint16_type Dim,
         typename GEOSHAPE,
         typename T = double,
         typename IndexT = uint32_type,
         bool UseMeasuresStorage = false >
class GeoElement3D
    :
        public GeoND<Dim, GEOSHAPE, T, IndexT, GeoElement0D<Dim, SubFaceOfNone<0>, T, IndexT>,UseMeasuresStorage >,
public SubFaceOfNone<0>
{
public:

    static inline const uint16_type nDim = Dim;

    typedef GeoND<Dim, GEOSHAPE, T, IndexT, GeoElement0D<Dim, SubFaceOfNone<0>, T, IndexT>, UseMeasuresStorage > super;
    typedef SubFaceOfNone<0> super2;

    using index_type = typename super::index_type;
    using size_type = typename super::size_type;
    typedef GEOSHAPE GeoShape;

    typedef typename super::face_type entity_face_type;

    typedef GeoElement3D<Dim, GEOSHAPE,T,IndexT,UseMeasuresStorage> self_type;
    typedef self_type element_type;
    typedef GeoElement2D<Dim, entity_face_type, SubFaceOf<self_type>, T, IndexT > face_type;
    typedef GeoElement1D<Dim, typename entity_face_type::topological_face_type, SubFaceOfMany<face_type>, T, IndexT> edge_type;
    typedef GeoElement0D<Dim, SubFaceOfNone<0>, T, IndexT> point_type;

    typedef typename super::node_type node_type;

    typedef typename super::vertex_permutation_type vertex_permutation_type;
    typedef typename super::edge_permutation_type edge_permutation_type;
    typedef typename super::face_permutation_type face_permutation_type;
    typedef typename super::face_permutation_type permutation_type;

    //! Number of local Vertices
    static inline const uint16_type numLocalVertices = super::numVertices;
    //! Number of local Faces
    static inline const uint16_type numLocalFaces = super::numFaces;
    //! Number of local Edges (using Euler Formula)
    static inline const uint16_type numLocalEdges = super::numEdges;

    /**
     *
     *
     * @param id identifier of the element
     */
    explicit GeoElement3D( size_type id = 0 )
        :
        super( id ),
        super2(),
        M_edges( numLocalEdges, nullptr ),
        M_faces( numLocalFaces, nullptr ),
        M_edge_permutation( numLocalEdges, edge_permutation_type( edge_permutation_type::IDENTITY ) ),
        M_face_permutation( numLocalFaces, face_permutation_type( face_permutation_type::IDENTITY ) )
    {
        //std::fill( M_faces.begin(), M_faces.end(), ( face_type* )0 );
        //std::fill( M_face_permutation.begin(), M_face_permutation.end(), face_permutation_type( face_permutation_type::IDENTITY ) );
    }

    /**
     * copy/move constructors
     */
    GeoElement3D( GeoElement3D const& g ) = default;
    GeoElement3D( GeoElement3D && g ) = default;

    /**
     * destructor
     */
    ~GeoElement3D() override
    {}

    /**
     * copy operator
     */
    GeoElement3D& operator=( GeoElement3D const& g ) = default;
    GeoElement3D& operator=( GeoElement3D && g ) = default;

    //void setMesh( MeshBase<> const* m ) { super::setMesh( m ); }
    MeshBase<> const* mesh() const
    {
        return super::mesh();
    }

    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    size_type id() const noexcept
    {
        return super::id();
    }

    bool isGhostFace() const
    {
        return super2::isGhostFace( super::processId()  );
    }


    /**
     * \return \p true if interprocess domain face, \p false otherwise
     */
    bool isInterProcessDomain() const
    {
        return super2::isInterProcessDomain( super::processId() );
    }

    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    bool isOnBoundary() const noexcept
    {
        return super::isOnBoundary();
    }

    /**
     * \return maximum \c dimension of the sub-entity touching the boundary of the element
     */
    uint16_type boundaryEntityDimension() const noexcept
    {
        return super::boundaryEntityDimension();
    }

    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const noexcept
    {
        return super::isGhostCell();
    }

    /**
     * \return process id
     */
    rank_type processId() const noexcept
    {
        return super::processId();
    }

    /**
     * \return process id
     */
    rank_type partition1() const
    {
        return super2::partition1( super::processId() );
    }

    /**
     * \return process id
     */
    rank_type partition2() const
    {
        return super2::partition2( super::processId() );
    }


    size_type ad_first() const
    {
        return invalid_v<size_type>;
    }
    uint16_type pos_first() const
    {
        return invalid_uint16_type_value;
    }
    size_type ad_second() const
    {
        return invalid_v<size_type>;
    }
    uint16_type pos_second() const
    {
        return invalid_uint16_type_value;
    }


    /**
     * \return true if GeoElement3D is connected to an edge
     */
    bool hasEdge( uint16_type i ) const
    {
        if ( i >= numLocalEdges )
            return false;
        return M_edges[i] != nullptr;
    }

    /**
     * \return true if GeoElement3D is connected to a face
     */
    bool hasFace( uint16_type i ) const
    {
        if ( i >= numLocalFaces )
            return false;
        return M_faces[i] != nullptr;
    }


    edge_type const& edge( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return *M_edges[i];
    }

    edge_type& edge( uint16_type i )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return *M_edges[i];
    }

    edge_type const* edgePtr( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        //DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return M_edges[i];
    }

    edge_permutation_type edgePermutation( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();


        return M_edge_permutation[i];
    }

    /**
     * Inserts an edge
     */
    void setEdge( uint16_type const i, edge_type const & p )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( boost::addressof( p ) ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();
        M_edges[i] = const_cast<edge_type*>( boost::addressof( p ) );
    }

    void setEdgePermutation( uint16_type i, edge_permutation_type o )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;

        M_edge_permutation[i] = o;
    }

    face_type const& face( uint16_type i ) const
    {
        DCHECK( i < numLocalFaces ) << "invalid local edge index elt " << this->id() << " face " << i;
        DCHECK( M_faces[i] ) << "invalid edge (null pointer) elt " << this->id() << " face " << i;
        return *M_faces[i];
    }

    face_type& face( uint16_type i )
    {
        DCHECK( i < numLocalFaces ) << "invalid local edge index elt " << this->id() << " face " << i;
        DCHECK( M_faces[i] ) << "invalid edge (null pointer) elt " << this->id() << " face " << i;
        return *M_faces[i];
    }

    face_type const* facePtr( uint16_type i ) const
    {
        DCHECK( i < numLocalFaces ) << "invalid local edge index elt " << this->id() << " face " << i;
        //FEELPP_ASSERT( M_faces[i] )( i ).error( "invalid edge (null pointer)" );
        return M_faces[i];
    }

    face_permutation_type facePermutation( uint16_type i ) const
    {
        DCHECK( i < numLocalFaces ) <<  "invalid local face index in elt " << this->id() << " face " << i;
        DCHECK( M_faces[i] ) <<  "invalid face (null pointer in elt " << this->id() << " face " << i;
        return M_face_permutation[i];
    }
    face_permutation_type permutation( uint16_type i ) const
    {
        DCHECK( i < numLocalFaces ) <<  "invalid local face index in elt " << this->id() << " face " << i;
        DCHECK( M_faces[i] ) <<  "invalid face (null pointer in elt " << this->id() << " face " << i;
        return M_face_permutation[i];
    }
    face_permutation_type permutation( uint16_type i, mpl::int_<1> ) const override
    {
        return this->facePermutation( i );
    }
    edge_permutation_type permutation( uint16_type i, mpl::int_<2> ) const
    {
        return this->edgePermutation( i );
    }
    vertex_permutation_type permutation( uint16_type, mpl::int_<3> ) const
    {
        return vertex_permutation_type();
    }

    /**
     * Inserts a face.
     */
    void setFace( uint16_type const i, face_type const & p )
    {
        M_faces[i] = const_cast<face_type*>( boost::addressof( p ) );
    }

    void setFacePermutation( uint16_type i, face_permutation_type o )
    {
        FEELPP_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local face index" );
        M_face_permutation[i] = o;
    }

    typedef typename ublas::bounded_array<face_type*, numLocalFaces>::iterator face_iterator;
    typedef typename ublas::bounded_array<face_type*, numLocalFaces>::const_iterator face_const_iterator;

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_iterator,face_iterator>
    faces()
    {
        return std::make_pair( M_faces.begin(), M_faces.end() );
    }

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_const_iterator,face_const_iterator>
    faces() const
    {
        return std::make_pair( M_faces.begin(), M_faces.end() );
    }
    std::vector<size_type> facesId() const
        {
            std::vector<size_type> fid;
            std::for_each( M_faces.begin(), M_faces.end(),
                           [&fid]( auto const& f ) { fid.push_back(f->id()); } );
            return fid;
        }
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & boost::serialization::base_object<super>( *this );
        }


private:

    ublas::bounded_array<edge_type*, numLocalEdges> M_edges;
    ublas::bounded_array<face_type*, numLocalFaces> M_faces;

    ublas::bounded_array<edge_permutation_type, numLocalEdges> M_edge_permutation;
    ublas::bounded_array<face_permutation_type, numLocalFaces> M_face_permutation;
};
/*@}*/




/*-------------------------------------------------------------------------
  GeoElement3D
  --------------------------------------------------------------------------*/


template <uint16_type Dim, typename SubFace, typename T, typename IndexT>
struct is_geoelement<GeoElement0D<Dim,SubFace,T,IndexT>>: std::true_type {};
template<uint16_type Dim,
         typename GEOSHAPE,
         typename SubFace,
         typename T,
         typename IndexT,
         bool PointTypeIsSubFaceOf,
         bool UseMeasuresStorage >
struct is_geoelement<GeoElement1D<Dim,GEOSHAPE,SubFace,T,IndexT,PointTypeIsSubFaceOf,UseMeasuresStorage>>: std::true_type {};
template<uint16_type Dim,
         typename GEOSHAPE,
         typename SubFace,
         typename T,
         typename IndexT,
         bool UseMeasuresStorage >
struct is_geoelement<GeoElement2D<Dim,GEOSHAPE,SubFace,T,IndexT,UseMeasuresStorage>>: std::true_type {};

template<uint16_type Dim,
         typename GEOSHAPE,
         typename T,
         typename IndexT,
         bool UseMeasuresStorage >
struct is_geoelement<GeoElement3D<Dim,GEOSHAPE,T,IndexT,UseMeasuresStorage>>: std::true_type {};


/**
 * @brief get if a face of an element has the marker @p flag
 * 
 * @tparam EltType type of mesh element
 * @return true if the element \p e has a face with \p flag, false otherwise
 */
template<typename EltType>
bool
hasFaceWithMarker( EltType const& e, boost::any const& flag,
                   std::enable_if_t<(dimension_v<EltType> > 0) && is_geoelement_v<EltType>>* = nullptr )
{
    flag_type theflag = e.mesh()->markerId( flag );
    // for( auto const& f : e.faces() )
    auto [ fbegin, fend ] = e.faces();
    for( auto f = fbegin; f != fend; ++f )
    {
        if ( *f && (*f)->hasMarker() )
        {
            if ( (*f)->marker().value() == theflag )
                return true;
        }
    }
    return false;
}

/**
 * @brief check if a element as faces with any of the string markers
 * 
 * @tparam EltType element type to be checked
 * @param e element to be checked
 * @param flags vector of strings containing the markers
 * @return bool true if has face with at least one of the markers, false otherwise
 */
template<typename EltType>
bool
hasFaceWithAnyOfTheMarkers( EltType const& e, std::vector<std::string> const& flags,
                            std::enable_if_t<(dimension_v<EltType> > 0) && is_geoelement_v<EltType>>* = nullptr )
{
    //flag_type theflag = e.mesh()->markerId( flag );
    // for( auto const& f : e.faces() )
    auto [ fbegin, fend ] = e.faces();
    for( auto f = fbegin; f != fend; ++f )
    {
        if ( *f && (*f)->hasMarker() )
        {
            if ( auto it = std::find( flags.begin(), flags.end(), e.mesh()->markerName( (*f)->marker().value() ) ); it !=  flags.end() )
            {
                return true;
            }
        }
    }
    return false;
}

} // Feel
#endif
