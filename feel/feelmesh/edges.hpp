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
   \file edges.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-09-03
 */
#ifndef __edges_H
#define __edges_H 1


#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include <feel/feelmesh/geoelement.hpp>

namespace Feel
{
namespace multi_index = boost::multi_index;

/// \cond \detail
/*!
  \class Edges
  \brief Edges container class

  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  @see
*/
template<typename EdgeType,typename FaceType>
class Edges
{
public:


    /** @name Typedefs
     */
    //@{

    typedef typename mpl::if_<mpl::equal_to<mpl::int_<EdgeType::nRealDim>, mpl::int_<3> >,
                              mpl::identity<GeoElement1D<3, EdgeType,SubFaceOfMany<FaceType> > >,
                              mpl::identity<boost::none_t> >::type::type edge_type;


    typedef multi_index::multi_index_container<
        edge_type,
        multi_index::indexed_by<
            // sort by employee::operator<
            multi_index::ordered_unique<multi_index::identity<edge_type> >,
            // sort by less<int> on marker
            multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_marker>,
                                            multi_index::composite_key<
                                            edge_type,
                                            multi_index::const_mem_fun<edge_type,
                                                                       Marker1 const&,
                                                                       &edge_type::marker>,
                                            multi_index::const_mem_fun<edge_type,
                                                                       rank_type,
                                                                       &edge_type::processId>
                                            > >,

            // sort by less<int> on processId
            multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_pid>,
                                            multi_index::const_mem_fun<edge_type,
                                                                       rank_type,
                                                                       &edge_type::processId> >,

            // sort by less<int> on boundary
            multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_location>,
                                            multi_index::const_mem_fun<edge_type,
                                                                       bool,
                                                                       &edge_type::isOnBoundary> >
            >
        > edges_type;


    typedef typename edges_type::iterator edge_iterator;
    typedef typename edges_type::const_iterator edge_const_iterator;
    typedef typename edges_type::template index<Feel::detail::by_marker>::type marker_edges;

    typedef typename marker_edges::iterator marker_edge_iterator;
    typedef typename marker_edges::const_iterator marker_edge_const_iterator;

    typedef typename edges_type::template index<Feel::detail::by_pid>::type pid_edges;
    typedef typename pid_edges::iterator pid_edge_iterator;
    typedef typename pid_edges::const_iterator pid_edge_const_iterator;

    typedef typename edges_type::template index<Feel::detail::by_location>::type location_edges;
    typedef typename location_edges::iterator location_edge_iterator;
    typedef typename location_edges::const_iterator location_edge_const_iterator;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Edges( WorldComm const& worldComm = Environment::worldComm() )
        :
        M_worldCommEdges( worldComm ),
        M_edges()
    {}

    Edges( Edges const & f )
        :
        M_worldCommEdges( f.M_worldCommEdges ),
        M_edges( f.M_edges )
    {}

    ~Edges()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the points container
     */
    edges_type & edges()
    {
        return M_edges;
    }

    /**
     * \return the edges container
     */
    edges_type const& edges() const
    {
        return M_edges;
    }

    /**
     * \return the world communicatior
     */
    WorldComm const& worldCommEdges() const
    {
        return M_worldCommEdges;
    }

    /**
     * \return true if container is empty, false otherwise
     */
    bool isEmpty() const
    {
        return M_edges.empty();
    }

    bool isBoundaryEdge( edge_type const & e ) const
    {
        return M_edges.find( e )->isOnBoundary();
    }
    bool isBoundaryEdge( size_type const & id ) const
    {
        return M_edges.find( edge_type( id ) )->isOnBoundary();
    }

    edge_type const& edge( size_type i ) const
    {
        return *M_edges.find( edge_type( i ) );
    }

    edge_iterator edgeIterator( size_type i ) const
    {
        return  M_edges.find( edge_type( i ) );
    }

    edge_iterator beginEdge()
    {
        return M_edges.begin();
    }
    edge_const_iterator beginEdge() const
    {
        return M_edges.begin();
    }
    edge_iterator endEdge()
    {
        return M_edges.end();
    }
    edge_const_iterator endEdge() const
    {
        return M_edges.end();
    }

    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with marker \p m on processor \p p
     */
    std::pair<marker_edge_iterator, marker_edge_iterator>
    edgesWithMarker( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommEdges().localRank() : p;
        return M_edges.template get<Feel::detail::by_marker>().equal_range( boost::make_tuple( Marker1( m ), part ) );
    }

    marker_edge_iterator beginEdgeWithMarker( size_type m )
    {
        return M_edges.template get<Feel::detail::by_marker>().lower_bound( Marker1( m ) );
    }
    marker_edge_const_iterator beginEdgeWithMarker( size_type m ) const
    {
        return M_edges.template get<Feel::detail::by_marker>().lower_bound( Marker1( m ) );
    }
    marker_edge_iterator endEdgeWithMarker( size_type m )
    {
        return M_edges.template get<Feel::detail::by_marker>().upper_bound( Marker1( m ) );
    }
    marker_edge_const_iterator endEdgeWithMarker( size_type m ) const
    {
        return M_edges.template get<Feel::detail::by_marker>().upper_bound( Marker1( m ) );
    }

    /**
     * get the edges container by id
     *
     *
     * @return the edge container by id
     */
    typename edges_type::template nth_index<0>::type &
    edgesById()
    {
        return M_edges.template get<0>();
    }

    /**
     * get the edges container by id
     *
     *
     * @return the edge container by id
     */
    typename edges_type::template nth_index<0>::type const&
    edgesById() const
    {
        return M_edges.template get<0>();
    }

    /**
     * get the edges container using the marker view
     *
     *
     * @return the edge container using marker view
     */
    marker_edges &
    edgesByMarker()
    {
        return M_edges.template get<Feel::detail::by_marker>();
    }

    /**
     * get the edges container using the marker view
     *
     *
     * @return the edge container using marker view
     */
    marker_edges const&
    edgesByMarker() const
    {
        return M_edges.template get<Feel::detail::by_marker>();
    }
    /**
     * get the edges container using the location view
     *
     *
     * @return the edge container using location view
     */
    location_edges &
    edgesByLocation()
    {
        return M_edges.template get<Feel::detail::by_location>();
    }

    /**
     * get the edges container using the location view
     *
     *
     * @return the edge container using location view
     */
    location_edges const&
    edgesByLocation() const
    {
        return M_edges.template get<Feel::detail::by_location>();
    }

    /**
     * get the begin() iterator on all the internal edges
     *
     * @return the begin() iterator on all the internal edges
     */
    location_edge_iterator beginInternalEdge()
    {
        return M_edges.template get<Feel::detail::by_location>().lower_bound( INTERNAL );
    }
    /**
     * get the end() iterator on all the internal edges
     *
     * @return the end() iterator on all the internal edges
     */
    location_edge_iterator endInternalEdge()
    {
        return M_edges.template get<Feel::detail::by_location>().upper_bound( INTERNAL );
    }

    /**
     * get the begin() iterator on all the internal edges
     *
     * @return the begin() iterator on all the internal edges
     */
    location_edge_const_iterator beginInternalEdge() const
    {
        return M_edges.template get<Feel::detail::by_location>().lower_bound( INTERNAL );
    }

    /**
     * get the end() iterator on all the internal edges
     *
     * @return the end() iterator on all the internal edges
     */
    location_edge_const_iterator endInternalEdge() const
    {
        return M_edges.template get<Feel::detail::by_location>().upper_bound( INTERNAL );
    }

    /**
     * get the begin() iterator on all the boundary edges
     *
     * @return the begin() iterator on all the boundary edges
     */
    location_edge_iterator beginEdgeOnBoundary()
    {
        return M_edges.template get<Feel::detail::by_location>().lower_bound( ON_BOUNDARY );
    }
    /**
     * get the end() iterator on all the boundary edges
     *
     * @return the end() iterator on all the boundary edges
     */
    location_edge_iterator endEdgeOnBoundary()
    {
        return M_edges.template get<Feel::detail::by_location>().upper_bound( ON_BOUNDARY );
    }

    /**
     * get the begin() iterator on all the boundary edges
     *
     * @return the begin() iterator on all the boundary edges
     */
    location_edge_const_iterator beginEdgeOnBoundary() const
    {
        return M_edges.template get<Feel::detail::by_location>().lower_bound( ON_BOUNDARY );
    }

    /**
     * get the end() iterator on all the boundary edges
     *
     * @return the end() iterator on all the boundary edges
     */
    location_edge_const_iterator endEdgeOnBoundary() const
    {
        return M_edges.template get<Feel::detail::by_location>().upper_bound( ON_BOUNDARY );
    }

    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with processor \p p
     */
    std::pair<pid_edge_iterator, pid_edge_iterator>
    edgesWithProcessId( rank_type p )
    {
        return M_edges.template get<Feel::detail::by_pid>().equal_range( p );
    }
    std::pair<pid_edge_const_iterator, pid_edge_const_iterator>
    edgesWithProcessId( rank_type p ) const
    {
        return M_edges.template get<Feel::detail::by_pid>().equal_range( p );
    }

    pid_edge_iterator beginEdgeWithProcessId( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommEdges().localRank() : p;
        return M_edges.template get<Feel::detail::by_pid>().lower_bound( /*boost::make_tuple( part )*/ part );
    }
    pid_edge_const_iterator beginEdgeWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommEdges().localRank() : p;
        return M_edges.template get<Feel::detail::by_pid>().lower_bound( /*boost::make_tuple( part )*/ part );
    }
    pid_edge_iterator endEdgeWithProcessId( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommEdges().localRank() : p;
        return M_edges.template get<Feel::detail::by_pid>().upper_bound( /*boost::make_tuple( part )*/ part );
    }
    pid_edge_const_iterator endEdgeWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommEdges().localRank() : p;
        return M_edges.template get<Feel::detail::by_pid>().upper_bound( /*boost::make_tuple( part )*/ part );
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
     * add a new edge in the mesh
     * @param f a new edge
     * @return the new edge from the list
     */
    edge_type const& addEdge( edge_type& f )
    {
        f.setId( M_edges.size() );
        return *M_edges.insert( f ).first;
    }

    void setWorldCommEdges( WorldComm const& _worldComm )
    {
        M_worldCommEdges = _worldComm;
    }

    //@}

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & M_edges;
        }

private:
    WorldComm M_worldCommEdges;
    edges_type M_edges;
};
/// \endcond
} // Feel
#endif /* __edges_H */
