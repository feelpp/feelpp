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

#include <unordered_map>

#include <feel/feelmesh/geoelement.hpp>

namespace Feel
{

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
    typedef typename FaceType::value_type value_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<EdgeType::nRealDim>, mpl::int_<3> >,
                              mpl::identity<GeoElement1D<3, EdgeType,SubFaceOfMany<FaceType>,value_type > >,
                              mpl::identity<boost::none_t> >::type::type edge_type;

    typedef std::unordered_map<size_type,edge_type> edges_type;

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
        return e.isOnBoundary();
    }
    bool isBoundaryEdge( size_type const & id ) const
    {
        return M_edges.find( id )->second.isOnBoundary();
    }

    /**
     * \return \c true if element with id \p i is found, \c false otherwise
     */
    bool hasEdge( size_type i ) const
    {
        return M_edges.find( i ) != M_edges.end();
    }

    edge_type const& edge( size_type i ) const
    {
        auto itFindEdge = M_edges.find( i );
        CHECK( itFindEdge != M_edges.end() ) << " edge " << i << "does not found";
        return itFindEdge->second;
    }

    edge_const_iterator edgeIterator( size_type i ) const
    {
        return  M_edges.find( i );
    }
    edge_iterator edgeIterator( size_type i )
    {
        return  M_edges.find( i );
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
                auto const& edge = it->second;
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
                auto const& edge = it->second;
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
                auto const& edge = it->second;
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
                auto const& edge = it->second;
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
                auto const& edge = it->second;
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
        //f.setId( M_edges.size() );
        return M_edges.emplace( std::make_pair( f.id(),f ) ).first->second;
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
