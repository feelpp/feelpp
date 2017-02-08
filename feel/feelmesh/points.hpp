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
#ifndef __points_H
#define __points_H 1

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include <feel/feelmesh/geoelement.hpp>

namespace Feel
{
namespace multi_index = boost::multi_index;

/// \cond detail
/*!
  \class Points
  \brief Points container class

  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  @see
*/
template<uint16_type nDim,typename T = double>
class Points
{
public:


    /** @name Typedefs
     */
    //@{

    typedef GeoElement0D<nDim,SubFaceOfNone,T> point_type;
    typedef multi_index::multi_index_container<
    point_type,
    multi_index::indexed_by<
    // sort by employee::operator<
        multi_index::ordered_unique<multi_index::identity<point_type> >
#if 0
        // sort by less<int> on marker
        multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_marker>,
                                        multi_index::composite_key<point_type,
                                                                   multi_index::const_mem_fun<point_type,
                                                                                              Marker1 const&,
                                                                                              &point_type::marker>,
                                                                   multi_index::const_mem_fun<point_type,
                                                                                              rank_type,
                                                                                              &point_type::processId> > >,
        // sort by less<int> on processId
        multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_pid>,
                                        multi_index::const_mem_fun<point_type,
                                                                   rank_type,
                                                                   &point_type::processId> >,

        // sort by less<int> on boundary
        multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_location>,
                                        multi_index::const_mem_fun<point_type,
                                                                   bool,
                                                                   &point_type::isOnBoundary> >
#endif
    >
    > points_type;


    typedef typename points_type::iterator point_iterator;
    typedef typename points_type::const_iterator point_const_iterator;

    typedef std::vector<boost::reference_wrapper<point_type const> > points_reference_wrapper_type;
    typedef std::shared_ptr<points_reference_wrapper_type> points_reference_wrapper_ptrtype;
    typedef typename points_reference_wrapper_type::iterator point_reference_wrapper_iterator;
    typedef typename points_reference_wrapper_type::const_iterator point_reference_wrapper_const_iterator;


#if 0
    typedef typename points_type::template index<Feel::detail::by_marker>::type marker_points;
    typedef typename marker_points::iterator marker_point_iterator;
    typedef typename marker_points::const_iterator marker_point_const_iterator;

    typedef typename points_type::template index<Feel::detail::by_pid>::type pid_points;
    typedef typename pid_points::iterator pid_point_iterator;
    typedef typename pid_points::const_iterator pid_point_const_iterator;

    typedef typename points_type::template index<Feel::detail::by_location>::type location_points;
    typedef typename location_points::iterator location_point_iterator;
    typedef typename location_points::const_iterator location_point_const_iterator;
#endif

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Points( WorldComm const& worldComm = Environment::worldComm() )
        :
        M_worldCommPoints( worldComm ),
        M_points()
    {}

    Points( Points const & f )
        :
        M_worldCommPoints( f.M_worldCommPoints ),
        M_points( f.M_points )
    {}

    virtual ~Points()
        {
            this->clear();
        }
    void clear()
    {
        VLOG(1) << "deleting points...\n";
        M_points.clear();
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    Points& operator=( Points const& e )
    {
        if ( this != &e )
        {
            M_worldCommPoints = e.M_worldCommPoints;
            M_points = e.M_points;
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
    points_type & points()
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
    bool isBoundaryPoint( point_type const & e ) const
    {
        return M_points.find( e )->isOnBoundary();
    }
    bool isBoundaryPoint( size_type const & id ) const
    {
        return M_points.find( point_type( id ) )->isOnBoundary();
    }


    point_type const& point( size_type i ) const
    {
        return *M_points.find( point_type( i ) );
    }

    point_iterator pointIterator( size_type i ) const
    {
        return  M_points.find( point_type( i ) );
    }

    bool hasPoint( size_type i ) const
    {
        return M_points.find( point_type( i ) ) != M_points.end();
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
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
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
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    pointsWithMarkerByType( uint16_type markerType, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommPoints().localRank() : p;
            points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
            auto it = this->beginPoint();
            auto en = this->endPoint();
            for ( ; it!=en;++it )
            {
                auto const& point = *it;
                if ( point.processId() != part )
                    continue;
                if ( !point.hasMarker( markerType ) )
                    continue;
                if ( point.marker().isOff() )
                    continue;
                mypoints->push_back(boost::cref(point));
            }
            return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
        }

    /**
     * \return iterator over marked points
     */
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    pointsWithMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommPoints().localRank() : p;
            points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
            auto it = this->beginPoint();
            auto en = this->endPoint();
            for ( ; it!=en;++it )
            {
                auto const& point = *it;
                if ( point.processId() != part )
                    continue;
                if ( !point.hasMarker( markerType ) )
                    continue;
                if ( point.marker( markerType ).isOff() )
                    continue;
                if ( markerFlags.find( point.marker( markerType ).value() ) == markerFlags.end() )
                    continue;
                mypoints->push_back(boost::cref(point));
            }
            return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
        }

    /**
     * \return iterator over marked points
     */
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    pointsWithMarkerByType( uint16_type markerType, flag_type m, rank_type p = invalid_rank_type_value ) const
        {
            if ( m == invalid_flag_type_value )
                return this->pointsWithMarkerByType( markerType, p );
            else
                return this->pointsWithMarkerByType( markerType, std::set<flag_type>( { m } ), p );
        }

    /**
     * \return iterator over marked points
     */
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
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
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    internalPoints( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommPoints().localRank() : p;
            points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
            auto it = this->beginPoint();
            auto en = this->endPoint();
            for ( ; it!=en;++it )
            {
                auto const& point = *it;
                if ( point.processId() != part )
                    continue;
                if ( !point.isInternal() )
                    continue;
                mypoints->push_back(boost::cref(point));
            }
            return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
        }
    /**
     * get iterator over boundary points
     *
     *
     * @return iterator over boundary points
     */
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    boundaryPoints( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommPoints().localRank() : p;
            points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
            auto it = this->beginPoint();
            auto en = this->endPoint();
            for ( ; it!=en;++it )
            {
                auto const& point = *it;
                if ( point.processId() != part )
                    continue;
                if ( !point.isOnBoundary() )
                    continue;
                mypoints->push_back(boost::cref(point));
            }
            return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
        }



    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    pointsWithProcessId( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommPoints().localRank() : p;
            points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
            auto it = this->beginPoint();
            auto en = this->endPoint();
            for ( ; it!=en;++it )
            {
                auto const& point = *it;
                if ( point.processId() == part )
                    mypoints->push_back(boost::cref(point));
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

    /**
     * add a new point in the mesh
     * @param f a new point
     * @return the new point from the list
     */
    point_type const& addPoint( point_type const& f )
    {
        return *M_points.insert( f ).first;
    }

    WorldComm const& worldCommPoints() const
    {
        return M_worldCommPoints;
    }

    void setWorldCommPoints( WorldComm const& _worldComm )
    {
        M_worldCommPoints = _worldComm;
    }

    //@}

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & M_points;
        }

private:
    WorldComm M_worldCommPoints;

    points_type M_points;
};
/// \endcond
} // Feel
#endif /* __points_H */
