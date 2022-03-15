/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-09-03

  Copyright (C) 2005,2006 EPFL

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
   \file points.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-09-03
 */
#ifndef FEELPP_MESH_POINTS_HPP
#define FEELPP_MESH_POINTS_HPP

#include <unordered_map>

#include <feel/feelcore/commobject.hpp>
#include <feel/feelmesh/geoelement.hpp>

namespace Feel
{
/// \cond detail
/*!
  \class Points
  \brief Points container class

  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  @see
*/
template <uint16_type nDim, typename T = double, typename IndexT = uint32_type, typename SubFace = SubFaceOfNone<0> >
class Points
{
  public:
    /** @name Typedefs
     */
    //@{
    using index_type = IndexT;
    using size_type = index_type;
    typedef GeoElement0D<nDim, SubFace, T, IndexT> point_type;

    typedef std::unordered_map<size_type, point_type> points_type;

    typedef typename points_type::iterator point_iterator;
    typedef typename points_type::const_iterator point_const_iterator;

    typedef std::vector<boost::reference_wrapper<point_type const>> points_reference_wrapper_type;
    typedef std::shared_ptr<points_reference_wrapper_type> points_reference_wrapper_ptrtype;
    typedef typename points_reference_wrapper_type::iterator point_reference_wrapper_iterator;
    typedef typename points_reference_wrapper_type::const_iterator point_reference_wrapper_const_iterator;

    typedef std::vector<boost::reference_wrapper<point_type>> ordered_points_reference_wrapper_type;
    typedef typename ordered_points_reference_wrapper_type::iterator ordered_point_reference_wrapper_iterator;
    typedef typename ordered_points_reference_wrapper_type::const_iterator ordered_point_reference_wrapper_const_iterator;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Points( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        : M_worldComm( worldComm ),
          M_points(),
          M_needToOrderPoints( false )
    {
    }

    Points( Points const& f )
        : M_worldComm( f.M_worldComm ),
          M_points( f.M_points ),
          M_needToOrderPoints( false )
    {
        this->buildOrderedPoints();
    }

    virtual ~Points() {}

    void clear()
    {
        DVLOG( 1 ) << "deleting points...\n";
        M_orderedPoints.clear();
        M_points.clear();
        M_needToOrderPoints = false;
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    Points& operator=( Points const& e )
    {
        if ( this != &e )
        {
            M_worldComm = e.M_worldComm;
            M_points = e.M_points;
            this->buildOrderedPoints();
        }

        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the points container
     */
    points_type& points()
    {
        return M_points;
    }

    /**
     * \return the points container
     */
    points_type const& points() const
    {
        return M_points;
    }

    virtual bool isEmpty() const
    {
        return M_points.empty();
    }
    bool isBoundaryPoint( point_type const& e ) const
    {
        return e.isOnBoundary();
    }
    bool isBoundaryPoint( size_type const& id ) const
    {
        auto itFindPt = M_points.find( id );
        if ( itFindPt == M_points.end() )
            return false;
        return itFindPt->isOnBoundary();
    }

    point_type& point( size_type i )
    {
        auto itFindPt = M_points.find( i );
        CHECK( itFindPt != M_points.end() ) << " point " << i << "does not found";
        return itFindPt->second;
    }
    point_type const& point( size_type i ) const
    {
        auto itFindPt = M_points.find( i );
        CHECK( itFindPt != M_points.end() ) << " point " << i << "does not found";
        return itFindPt->second;
    }

    point_const_iterator pointIterator( size_type i ) const
    {
        return M_points.find( i );
    }
    point_iterator pointIterator( size_type i )
    {
        return M_points.find( i );
    }

    bool hasPoint( size_type i ) const
    {
        return M_points.find( i ) != M_points.end();
    }

    point_iterator beginPoint()
    {
        return M_points.begin();
    }
    point_const_iterator beginPoint() const
    {
        return M_points.begin();
    }
    point_iterator endPoint()
    {
        return M_points.end();
    }
    point_const_iterator endPoint() const
    {
        return M_points.end();
    }

    ordered_point_reference_wrapper_iterator beginOrderedPoint()
    {
        return M_orderedPoints.begin();
    }
    ordered_point_reference_wrapper_const_iterator beginOrderedPoint() const
    {
        return M_orderedPoints.begin();
    }
    ordered_point_reference_wrapper_iterator endOrderedPoint()
    {
        return M_orderedPoints.end();
    }
    ordered_point_reference_wrapper_const_iterator endOrderedPoint() const
    {
        return M_orderedPoints.end();
    }

    /**
     * \return the range of iterator \c (begin,end) over the points
     * with marker \p m on processor \p p
     */
    std::pair<point_iterator, point_iterator>
    pointsRange()
    {
        return std::make_pair( M_points.begin(), M_points.end() );
    }

    /**
     * \return the range of iterator \c (begin,end) over the points
     * with marker \p m on processor \p p
     */
    std::pair<point_const_iterator, point_const_iterator>
    pointsRange() const
    {
        return std::make_pair( M_points.begin(), M_points.end() );
    }

    /**
     * \return the range of iterator \c (begin,end) over the points
     * with \c Id \p m
     */
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithId( size_type m ) const
    {
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        if ( this->hasPoint( m ) )
            mypoints->push_back( boost::cref( this->point( m ) ) );
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }

    /**
     * \return iterator over marked points
     */
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithMarkerByType( uint16_type markerType, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = ( p == invalid_rank_type_value ) ? this->worldCommPoints().localRank() : p;
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        auto it = this->beginOrderedPoint();
        auto en = this->endOrderedPoint();
        for ( ; it != en; ++it )
        {
            auto const& point = unwrap_ref( *it );
            if ( point.processId() != part )
                continue;
            if ( !point.hasMarkerType( markerType ) )
                continue;
            if ( point.marker().isOff() )
                continue;
            mypoints->push_back( boost::cref( point ) );
        }
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }

    /**
     * \return iterator over marked points
     */
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = ( p == invalid_rank_type_value ) ? this->worldCommPoints().localRank() : p;
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        auto it = this->beginOrderedPoint();
        auto en = this->endOrderedPoint();
        for ( ; it != en; ++it )
        {
            auto const& point = unwrap_ref( *it );
            if ( point.processId() != part )
                continue;
            if ( !point.hasMarkerType( markerType ) )
                continue;
            if ( !point.marker( markerType ).hasOneOf( markerFlags ) )
                continue;
            mypoints->push_back( boost::cref( point ) );
        }
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }

    /**
     * \return iterator over marked points
     */
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithMarkerByType( uint16_type markerType, flag_type m, rank_type p = invalid_rank_type_value ) const
    {
        if ( m == invalid_flag_type_value )
            return this->pointsWithMarkerByType( markerType, p );
        else
            return this->pointsWithMarkerByType( markerType, std::set<flag_type>( {m} ), p );
    }

    /**
     * \return iterator over marked points
     */
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithMarker( flag_type m = invalid_flag_type_value, rank_type p = invalid_rank_type_value ) const
    {
        return this->pointsWithMarkerByType( 1, m, p );
    }

    /**
     * get iterator over internal points
     *
     *
     * @return iterator over internal points
     */
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    internalPoints( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = ( p == invalid_rank_type_value ) ? this->worldCommPoints().localRank() : p;
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        auto it = this->beginOrderedPoint();
        auto en = this->endOrderedPoint();
        for ( ; it != en; ++it )
        {
            auto const& point = unwrap_ref( *it );
            if ( point.processId() != part )
                continue;
            if ( !point.isInternal() )
                continue;
            mypoints->push_back( boost::cref( point ) );
        }
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }
    /**
     * get iterator over boundary points
     *
     *
     * @return iterator over boundary points
     */
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    boundaryPoints( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = ( p == invalid_rank_type_value ) ? this->worldCommPoints().localRank() : p;
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        auto it = this->beginOrderedPoint();
        auto en = this->endOrderedPoint();
        for ( ; it != en; ++it )
        {
            auto const& point = unwrap_ref( *it );
            if ( point.processId() != part )
                continue;
            if ( !point.isOnBoundary() )
                continue;
            mypoints->push_back( boost::cref( point ) );
        }
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }

    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = ( p == invalid_rank_type_value ) ? this->worldCommPoints().localRank() : p;
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        auto it = this->beginOrderedPoint();
        auto en = this->endOrderedPoint();
        for ( ; it != en; ++it )
        {
            auto const& point = unwrap_ref( *it );
            if ( point.processId() == part )
                mypoints->push_back( boost::cref( point ) );
        }

        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }

    template <typename SF = SubFace>
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    interProcessPoints( rank_type neighbor_pid = invalid_rank_type_value,
                        std::enable_if_t< std::is_base_of<SubFaceOfBase, SF>::value >* = nullptr ) const
        {
            bool allNeighbor = ( neighbor_pid == invalid_rank_type_value );
            const rank_type part = this->worldCommPoints().localRank();
            points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
            auto it = this->beginOrderedPoint();
            auto en = this->endOrderedPoint();
            for ( ; it!=en;++it )
            {
                auto const& point = unwrap_ref( *it );
                if ( !point.isInterProcessDomain() )
                    continue;
                if ( point.partition1() != part )
                    continue;
                if ( !allNeighbor && point.partition2() != neighbor_pid )
                    continue;
                mypoints->push_back( boost::cref( point ) );
            }
            return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
        }


    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //!
    //! reserve size of container
    //! @param nPoint : the size reserved
    //!
    void reserveNumberOfPoint( size_type nPoint )
    {
        M_points.reserve( nPoint );
        M_orderedPoints.reserve( nPoint );
    }

    /**
     * add a new point in the mesh
     * @param f a new point
     * @return the new point from the list
     */
    std::pair<point_iterator,bool> addPoint( point_type const& f )
    {
        //return M_points.insert( std::make_pair( f.id(), f ) ).first->second;
        auto ret = M_points.emplace( std::make_pair( f.id(), f ) );

        auto& newPoint = ret.first->second;
        if ( ret.second )
        {
            if ( !M_needToOrderPoints && !M_orderedPoints.empty() && unwrap_ref( M_orderedPoints.back() ).id() > newPoint.id() )
                M_needToOrderPoints = true;
            M_orderedPoints.push_back( boost::ref( newPoint ) );
        }
        return ret;
    }

    /**
     * add a new point in the mesh
     * @param f a new point
     * @return the new point from the list
     */
    std::pair<point_iterator,bool> addPoint( point_type&& f )
    {
        //return M_points.insert( std::make_pair( f.id(), f ) ).first->second;
        auto ret = M_points.emplace( std::make_pair( f.id(), f ) );

        auto& newPoint = ret.first->second;
        if ( ret.second )
        {
            if ( !M_needToOrderPoints && !M_orderedPoints.empty() && unwrap_ref( M_orderedPoints.back() ).id() > newPoint.id() )
                M_needToOrderPoints = true;
            M_orderedPoints.push_back( boost::ref( newPoint ) );
        }
        return ret;
    }

    /**
     * erase point at position \p position
     *
     * @param position \p position is a valid dereferenceable iterator of the index.
     *
     * @return An iterator pointing to the point immediately
     * following the one that was deleted, or \c end() if no such point
     * exists.
     */
    point_iterator erasePoint( point_iterator it )
    {
        size_type erasedId = it->first;
        auto itret = M_points.erase( it );
        auto itOrdered = std::find_if( M_orderedPoints.begin(), M_orderedPoints.end(),
                                       [&erasedId]( auto& pointWrap ) { return unwrap_ref( pointWrap ).id() == erasedId; } );
        M_orderedPoints.erase( itOrdered );
        return itret;
    }

    WorldComm const& worldCommPoints() const
    {
        return *M_worldComm;
    }

    void setWorldCommPoints( worldcomm_ptr_t const& _worldComm )
    {
        M_worldComm = _worldComm;
    }

    //@}

    void updateOrderedPoints()
    {
        if ( !M_needToOrderPoints )
            return;
        std::sort( M_orderedPoints.begin(), M_orderedPoints.end(),
                   []( auto const& a, auto const& b ) -> bool {
                       return unwrap_ref( a ).id() < unwrap_ref( b ).id();
                   } );
        M_needToOrderPoints = false;
    }

  private:
    void buildOrderedPoints()
    {
        M_orderedPoints.clear();
        auto it = beginPoint(), en = endPoint();
        size_type nPoint = std::distance( it, en );
        M_orderedPoints.reserve( nPoint );
        for ( ; it != en; ++it )
            M_orderedPoints.push_back( boost::ref( it->second ) );
        M_needToOrderPoints = true;
        this->updateOrderedPoints();
    }

    friend class boost::serialization::access;
    template <class Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
        if ( Archive::is_loading::value )
        {
            M_points.clear();
            M_orderedPoints.clear();
            M_needToOrderPoints = false;
            size_type nPoints = 0;
            ar& BOOST_SERIALIZATION_NVP( nPoints );
            point_type newPoint;
            for ( size_type k = 0; k < nPoints; ++k )
            {
                ar& boost::serialization::make_nvp( "point", newPoint );
                this->addPoint( std::move( newPoint ) );
            }
        }
        else
        {
            auto it = beginOrderedPoint(), en = endOrderedPoint();
            size_type nPoints = std::distance( it, en );
            ar& BOOST_SERIALIZATION_NVP( nPoints );
            for ( ; it != en; ++it )
                ar& boost::serialization::make_nvp( "point", unwrap_ref( *it ) );
        }
    }

  private:
    worldcomm_ptr_t M_worldComm;

    points_type M_points;
    ordered_points_reference_wrapper_type M_orderedPoints;
    bool M_needToOrderPoints;
};
/// \endcond
} // namespace Feel
#endif /* __points_H */
