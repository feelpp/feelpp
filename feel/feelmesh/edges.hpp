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
            multi_index::ordered_unique<multi_index::identity<edge_type> >
#if 0
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
#endif
            >
        > edges_type;


    typedef typename edges_type::iterator edge_iterator;
    typedef typename edges_type::const_iterator edge_const_iterator;

    typedef std::vector<boost::reference_wrapper<edge_type const> > edges_reference_wrapper_type;
    typedef std::shared_ptr<edges_reference_wrapper_type> edges_reference_wrapper_ptrtype;
    typedef typename edges_reference_wrapper_type::iterator edge_reference_wrapper_iterator;
    typedef typename edges_reference_wrapper_type::const_iterator edge_reference_wrapper_const_iterator;

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

    /**
     * \return \c true if element with id \p i is found, \c false otherwise
     */
    bool hasEdge( size_type i ) const
    {
        return M_edges.template get<0>().find( edge_type( i ) ) !=
               M_edges.template get<0>().end();
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
     * \return the range of iterator \c (begin,end) over the edges
     * with \c Marker1 \p m on processor \p p
     */
    std::tuple<edge_reference_wrapper_const_iterator,edge_reference_wrapper_const_iterator,edges_reference_wrapper_ptrtype>
    edgesWithMarkerByType( uint16_type markerType, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommEdges().localRank() : p;
            edges_reference_wrapper_ptrtype myedges( new edges_reference_wrapper_type );
            auto it = this->beginEdge();
            auto en = this->endEdge();
            for ( ; it!=en;++it )
            {
                auto const& edge = *it;
                if ( edge.processId() != part )
                    continue;
                if ( !edge.hasMarker( markerType ) )
                    continue;
                if ( edge.marker().isOff() )
                    continue;
                myedges->push_back( boost::cref( edge ) );
            }
            return std::make_tuple( myedges->begin(), myedges->end(), myedges );
        }
    /**
     * \return the range of iterator \c (begin,end) over the edges
     * with \c Marker1 \p m on processor \p p
     */
    std::tuple<edge_reference_wrapper_const_iterator,edge_reference_wrapper_const_iterator,edges_reference_wrapper_ptrtype>
    edgesWithMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommEdges().localRank() : p;
            edges_reference_wrapper_ptrtype myedges( new edges_reference_wrapper_type );
            auto it = this->beginEdge();
            auto en = this->endEdge();
            for ( ; it!=en;++it )
            {
                auto const& edge = *it;
                if ( edge.processId() != part )
                    continue;
                if ( !edge.hasMarker( markerType ) )
                    continue;
                if ( edge.marker( markerType ).isOff() )
                    continue;
                if ( markerFlags.find( edge.marker( markerType ).value() ) == markerFlags.end() )
                    continue;
                myedges->push_back( boost::cref( edge ) );
            }
            return std::make_tuple( myedges->begin(), myedges->end(), myedges );
        }

    /**
     * \return the range of iterator \c (begin,end) over the edges
     * with \c Marker1 \p m on processor \p p
     */
    std::tuple<edge_reference_wrapper_const_iterator,edge_reference_wrapper_const_iterator,edges_reference_wrapper_ptrtype>
    edgesWithMarkerByType( uint16_type markerType, flag_type m, rank_type p = invalid_rank_type_value ) const
        {
            if ( m == invalid_flag_type_value )
                return this->edgesWithMarkerByType( markerType, p );
            else
                return this->edgesWithMarkerByType( markerType, std::set<flag_type>( { m } ), p );
        }

    /**
     * \return the range of iterator \c (begin,end) over the edges
     * with \c Marker1 \p m on processor \p p
     */
    std::tuple<edge_reference_wrapper_const_iterator,edge_reference_wrapper_const_iterator,edges_reference_wrapper_ptrtype>
    edgesWithMarker( flag_type m = invalid_flag_type_value, rank_type p = invalid_rank_type_value ) const
        {
            return this->edgesWithMarkerByType( 1, m, p );
        }

    /**
     * \return the range of iterator \c (begin,end) over the boundary
     *  edges on processor \p p
     */
    std::tuple<edge_reference_wrapper_const_iterator,edge_reference_wrapper_const_iterator,edges_reference_wrapper_ptrtype>
    edgesOnBoundary( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommEdges().localRank() : p;
            edges_reference_wrapper_ptrtype myedges( new edges_reference_wrapper_type );
            auto it = this->beginEdge();
            auto en = this->endEdge();
            for ( ; it!=en;++it )
            {
                auto const& edge = *it;
                if ( edge.processId() != part )
                    continue;
                if ( !edge.isOnBoundary() )
                    continue;
                myedges->push_back( boost::cref( edge ) );
            }
            return std::make_tuple( myedges->begin(), myedges->end(), myedges );
        }

    /**
     * \return the range of iterator \c (begin,end) over the internal edges
     * on processor \p p
     */
    std::tuple<edge_reference_wrapper_const_iterator,edge_reference_wrapper_const_iterator,edges_reference_wrapper_ptrtype>
    internalEdges( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommEdges().localRank() : p;
            edges_reference_wrapper_ptrtype myedges( new edges_reference_wrapper_type );
            auto it = this->beginEdge();
            auto en = this->endEdge();
            for ( ; it!=en;++it )
            {
                auto const& edge = *it;
                if ( edge.processId() != part )
                    continue;
                if ( !edge.isInternal() )
                    continue;
                myedges->push_back( boost::cref( edge ) );
            }
            return std::make_tuple( myedges->begin(), myedges->end(), myedges );
        }

    /**
     * \return the range of iterator \c (begin,end) over the edges
     * on processor \p p
     */
    std::tuple<edge_reference_wrapper_const_iterator,edge_reference_wrapper_const_iterator,edges_reference_wrapper_ptrtype>
    edgesWithProcessId( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommEdges().localRank() : p;
            edges_reference_wrapper_ptrtype myedges( new edges_reference_wrapper_type );
            auto it = this->beginEdge();
            auto en = this->endEdge();
            for ( ; it!=en;++it )
            {
                auto const& edge = *it;
                if ( edge.processId() != part )
                    continue;
                myedges->push_back( boost::cref( edge ) );
            }
            return std::make_tuple( myedges->begin(), myedges->end(), myedges );
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
