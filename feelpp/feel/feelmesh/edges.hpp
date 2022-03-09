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
#ifndef FEELPP_MESH_EDGES_HPP
#define FEELPP_MESH_EDGES_HPP

#include <unordered_map>
#include <feel/feelcore/commobject.hpp>
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
    using index_type = typename FaceType::index_type;
    using size_type = typename FaceType::size_type;
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

    typedef std::vector<boost::reference_wrapper<edge_type> > ordered_edges_reference_wrapper_type;
    typedef typename ordered_edges_reference_wrapper_type::iterator ordered_edge_reference_wrapper_iterator;
    typedef typename ordered_edges_reference_wrapper_type::const_iterator ordered_edge_reference_wrapper_const_iterator;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Edges( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        :
        M_worldComm( worldComm ),
        M_edges(),
        M_needToOrderEdges( false )
    {}

    Edges( Edges const & f )
        :
        M_worldComm( f.M_worldComm ),
        M_edges( f.M_edges ),
        M_needToOrderEdges( false )
    {
        this->buildOrderedEdges();
    }

    virtual ~Edges() {}

    void clear()
        {
            DVLOG(1) << "deleting edges...\n";
            M_orderedEdges.clear();
            M_edges.clear();
            M_needToOrderEdges = false;
        }

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
        return *M_worldComm;
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

    ordered_edge_reference_wrapper_iterator beginOrderedEdge()
        {
            return M_orderedEdges.begin();
        }
    ordered_edge_reference_wrapper_const_iterator beginOrderedEdge() const
        {
            return M_orderedEdges.begin();
        }
    ordered_edge_reference_wrapper_iterator endOrderedEdge()
        {
            return M_orderedEdges.end();
        }
    ordered_edge_reference_wrapper_const_iterator endOrderedEdge() const
        {
            return M_orderedEdges.end();
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
            auto it = this->beginOrderedEdge();
            auto en = this->endOrderedEdge();
            for ( ; it!=en;++it )
            {
                auto const& edge = unwrap_ref( *it );
                if ( edge.processId() != part )
                    continue;
                if ( !edge.hasMarkerType( markerType ) )
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
            auto it = this->beginOrderedEdge();
            auto en = this->endOrderedEdge();
            for ( ; it!=en;++it )
            {
                auto const& edge = unwrap_ref( *it );
                if ( edge.processId() != part )
                    continue;
                if ( !edge.hasMarkerType( markerType ) )
                    continue;
                if ( !edge.marker( markerType ).hasOneOf( markerFlags ) )
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
            auto it = this->beginOrderedEdge();
            auto en = this->endOrderedEdge();
            for ( ; it!=en;++it )
            {
                auto const& edge = unwrap_ref( *it );
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
            auto it = this->beginOrderedEdge();
            auto en = this->endOrderedEdge();
            for ( ; it!=en;++it )
            {
                auto const& edge = unwrap_ref( *it );
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
            auto it = this->beginOrderedEdge();
            auto en = this->endOrderedEdge();
            for ( ; it!=en;++it )
            {
                auto const& edge = unwrap_ref( *it );
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
    std::pair<edge_iterator,bool>  addEdge( edge_type& f )
    {
        //f.setId( M_edges.size() );
        std::pair<edge_iterator,bool> ret = M_edges.emplace( std::make_pair( f.id(),f ) );

        if ( ret.second )
        {
            auto & newEdge = ret.first->second;
            if ( !M_needToOrderEdges && !M_orderedEdges.empty() && unwrap_ref( M_orderedEdges.back() ).id() > newEdge.id() )
                M_needToOrderEdges = true;
            M_orderedEdges.push_back( boost::ref( newEdge ) );
        }

        return ret;
    }

    std::pair<edge_iterator,bool>  addEdge( edge_type&& f )
        {
            //f.setId( M_edges.size() );
            std::pair<edge_iterator,bool> ret = M_edges.emplace( std::make_pair( f.id(),f ) );

            if ( ret.second )
            {
                auto & newEdge = ret.first->second;
                if ( !M_needToOrderEdges && !M_orderedEdges.empty() && unwrap_ref( M_orderedEdges.back() ).id() > newEdge.id() )
                    M_needToOrderEdges = true;
                M_orderedEdges.push_back( boost::ref( newEdge ) );
            }

            return ret;
        }

    /**
     * erase edge at position \p position
     *
     * @param position \p position is a valid dereferenceable iterator of the index.
     *
     * @return An iterator pointing to the edge immediately
     * following the one that was deleted, or \c end() if no such edge
     * exists.
     */
    edge_iterator eraseEdge( edge_iterator it )
        {
            size_type erasedId = it->first;
            auto itret = M_edges.erase( it );
            auto itOrdered = std::find_if( M_orderedEdges.begin(), M_orderedEdges.end(),
                                           [&erasedId]( auto & edgeWrap ) { return unwrap_ref( edgeWrap ).id() == erasedId; } );
            M_orderedEdges.erase( itOrdered );
            return itret;
        }


    void setWorldCommEdges( worldcomm_ptr_t const& _worldComm )
    {
        M_worldComm  = _worldComm;
    }

    void updateOrderedEdges()
        {
            if ( !M_needToOrderEdges )
                return;
            std::sort( M_orderedEdges.begin(), M_orderedEdges.end(),
                       []( auto const& a, auto const& b) -> bool
                       {
                           return unwrap_ref( a ).id() < unwrap_ref( b ).id();
                       });
            M_needToOrderEdges = false;
        }
    //@}

private:

    void buildOrderedEdges()
        {
            M_orderedEdges.clear();
            auto it = beginEdge(), en = endEdge();
            size_type nEdge = std::distance( it, en );
            M_orderedEdges.reserve( nEdge );
            for ( ; it != en ; ++it )
                M_orderedEdges.push_back( boost::ref( it->second ) );
            M_needToOrderEdges = true;
            this->updateOrderedEdges();
        }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            if ( Archive::is_loading::value )
            {
                M_edges.clear();
                M_orderedEdges.clear();
                M_needToOrderEdges = false;
                size_type nEdges = 0;
                ar & BOOST_SERIALIZATION_NVP( nEdges );
                edge_type newEdge;
                for ( size_type k=0 ; k<nEdges ; ++k )
                {
                    ar & boost::serialization::make_nvp( "edge", newEdge );
                    this->addEdge( std::move( newEdge ) );
                }
            }
            else
            {
                auto it = beginOrderedEdge(), en = endOrderedEdge();
                size_type nEdges = std::distance( it, en );
                ar & BOOST_SERIALIZATION_NVP( nEdges );
                for ( ; it != en ; ++it )
                    ar & boost::serialization::make_nvp( "edge", unwrap_ref( *it ) );
            }
        }

private:
    worldcomm_ptr_t M_worldComm;
    edges_type M_edges;
    ordered_edges_reference_wrapper_type M_orderedEdges;
    bool M_needToOrderEdges;
};
/// \endcond
} // Feel
#endif /* __edges_H */
