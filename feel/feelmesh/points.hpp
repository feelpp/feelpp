/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
template<uint16_type nDim>
class Points
{
public:


    /** @name Typedefs
     */
    //@{

    typedef GeoElement0D<nDim> point_type;
    typedef multi_index::multi_index_container<
    point_type,
    multi_index::indexed_by<
    // sort by employee::operator<
        multi_index::ordered_unique<multi_index::identity<point_type> >,
        // sort by less<int> on marker
        multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_marker>,
                                        multi_index::const_mem_fun<point_type,
                                                                   Marker1 const&,
                                                                   &point_type::marker> >,

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
    >
    > points_type;


    typedef typename points_type::iterator point_iterator;
    typedef typename points_type::const_iterator point_const_iterator;

    typedef typename points_type::template index<Feel::detail::by_marker>::type marker_points;
    typedef typename marker_points::iterator marker_point_iterator;
    typedef typename marker_points::const_iterator marker_point_const_iterator;

    typedef typename points_type::template index<Feel::detail::by_pid>::type pid_points;
    typedef typename pid_points::iterator pid_point_iterator;
    typedef typename pid_points::const_iterator pid_point_const_iterator;

    typedef typename points_type::template index<Feel::detail::by_location>::type location_points;
    typedef typename location_points::iterator location_point_iterator;
    typedef typename location_points::const_iterator location_point_const_iterator;



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
    };

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


    marker_point_iterator beginPointWithMarker( size_type m )
    {
        return M_points.template get<Feel::detail::by_marker>().lower_bound( Marker1(m) );
    }
    marker_point_const_iterator beginPointWithMarker( size_type m ) const
    {
        return M_points.template get<Feel::detail::by_marker>().lower_bound( Marker1(m) );
    }
    marker_point_iterator endPointWithMarker( size_type m )
    {
        return M_points.template get<Feel::detail::by_marker>().upper_bound( Marker1(m) );
    }
    marker_point_const_iterator endPointWithMarker( size_type m ) const
    {
        return M_points.template get<Feel::detail::by_marker>().upper_bound( Marker1(m) );
    }

    point_iterator pointIterator( size_type i ) const
    {
        return  M_points.find( point_type( i ) );
    }

    bool hasPoint( size_type i ) const
        {
            return M_points.find( point_type( i ) ) != M_points.end();
        }
    /**
     * get the points container by id
     *
     *
     * @return the point container by id
     */
    typename points_type::template nth_index<0>::type &
    pointsById()
    {
        return M_points.template get<0>();
    }

    /**
     * get the points container by id
     *
     *
     * @return the point container by id
     */
    typename points_type::template nth_index<0>::type const&
    pointsById() const
    {
        return M_points.template get<0>();
    }

    /**
     * get the points container using the marker view
     *
     *
     * @return the point container using marker view
     */
    marker_points &
    pointsByMarker()
    {
        return M_points.template get<Feel::detail::by_marker>();
    }

    /**
     * get the points container using the marker view
     *
     *
     * @return the point container using marker view
     */
    marker_points const&
    pointsByMarker() const
    {
        return M_points.template get<Feel::detail::by_marker>();
    }
    /**
     * get the points container using the location view
     *
     *
     * @return the point container using location view
     */
    location_points &
    pointsByLocation()
    {
        return M_points.template get<Feel::detail::by_location>();
    }

    /**
     * get the points container using the location view
     *
     *
     * @return the point container using location view
     */
    location_points const&
    pointsByLocation() const
    {
        return M_points.template get<Feel::detail::by_location>();
    }

    /**
     * get the begin() iterator on all the internal points
     *
     * @return the begin() iterator on all the internal points
     */
    location_point_iterator beginInternalPoint()
    {
        return M_points.template get<Feel::detail::by_location>().lower_bound( INTERNAL );
    }
    /**
     * get the end() iterator on all the internal points
     *
     * @return the end() iterator on all the internal points
     */
    location_point_iterator endInternalPoint()
    {
        return M_points.template get<Feel::detail::by_location>().upper_bound( INTERNAL );
    }

    /**
     * get the begin() iterator on all the internal points
     *
     * @return the begin() iterator on all the internal points
     */
    location_point_const_iterator beginInternalPoint() const
    {
        return M_points.template get<Feel::detail::by_location>().lower_bound( INTERNAL );
    }

    /**
     * get the end() iterator on all the internal points
     *
     * @return the end() iterator on all the internal points
     */
    location_point_const_iterator endInternalPoint() const
    {
        return M_points.template get<Feel::detail::by_location>().upper_bound( INTERNAL );
    }

    /**
     * get the begin() iterator on all the boundary points
     *
     * @return the begin() iterator on all the boundary points
     */
    location_point_iterator beginPointOnBoundary()
    {
        return M_points.template get<Feel::detail::by_location>().lower_bound( ON_BOUNDARY );
    }
    /**
     * get the end() iterator on all the boundary points
     *
     * @return the end() iterator on all the boundary points
     */
    location_point_iterator endPointOnBoundary()
    {
        return M_points.template get<Feel::detail::by_location>().upper_bound( ON_BOUNDARY );
    }

    /**
     * get the begin() iterator on all the boundary points
     *
     * @return the begin() iterator on all the boundary points
     */
    location_point_const_iterator beginPointOnBoundary() const
    {
        return M_points.template get<Feel::detail::by_location>().lower_bound( ON_BOUNDARY );
    }

    /**
     * get the end() iterator on all the boundary points
     *
     * @return the end() iterator on all the boundary points
     */
    location_point_const_iterator endPointOnBoundary() const
    {
        return M_points.template get<Feel::detail::by_location>().upper_bound( ON_BOUNDARY );
    }


    pid_point_iterator beginPointWithProcessId( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommPoints().localRank() : p;
        return M_points.template get<Feel::detail::by_pid>().lower_bound( /*boost::make_tuple( part )*/ part );
    }
    pid_point_const_iterator beginPointWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommPoints().localRank() : p;
        return M_points.template get<Feel::detail::by_pid>().lower_bound( /*boost::make_tuple( part )*/ part );
    }
    pid_point_iterator endPointWithProcessId( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommPoints().localRank() : p;
        return M_points.template get<Feel::detail::by_pid>().upper_bound( /*boost::make_tuple( part )*/ part );
    }
    pid_point_const_iterator endPointWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommPoints().localRank() : p;
        return M_points.template get<Feel::detail::by_pid>().upper_bound( /*boost::make_tuple( part )*/ part );
    }


    std::pair<pid_point_iterator, pid_point_iterator>
    pointsWithProcessId( size_type p ) const
    {
        return M_points.template get<Feel::detail::by_pid>().equal_range( p );
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
