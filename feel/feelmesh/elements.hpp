/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-09-03

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007,2008,2009,2010 Université Joseph Fourier (Grenoble I)

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
   \file elements.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-09-03
 */
#ifndef FEELPP_ELEMENTS_HPP
#define FEELPP_ELEMENTS_HPP 1


#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

#include <feel/feelmesh/geoelement.hpp>
#include <feel/feelmesh/filters.hpp>

namespace Feel
{
namespace multi_index = boost::multi_index;
/// \cond detail

namespace detail
{
    template <typename EltType >
    void
    updateElementGhostConnectEdgeToElement( EltType& e, int i, mpl::int_<1> /**/)
    {}
    template <typename EltType >
    void
    updateElementGhostConnectEdgeToElement( EltType& e, int i, mpl::int_<2> /**/)
    {}
    template <typename EltType >
    void
    updateElementGhostConnectEdgeToElement( EltType& e, int i, mpl::int_<3> /**/)
    {
        e.edge( i ).addElementGhost( e.processId(),e.id() );
    }
}


/*!
  \class Elements
  \brief Elements container class

  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  @see
*/
template<typename ElementType>
class Elements
{
public:


    /** @name Typedefs
     */
    //@{

    /**
     * Element type depending on the dimension, @see geoelement.hpp
     * \note Elements have their topological dimension equal to the
     * dimension of the geometric space.
     */
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<ElementType::nDim>, mpl::int_<3> >,
            mpl::identity<GeoElement3D<ElementType::nRealDim, ElementType> >,
            typename mpl::if_<mpl::equal_to<mpl::int_<ElementType::nDim>, mpl::int_<2> >,
            mpl::identity<GeoElement2D<ElementType::nRealDim, ElementType> >,
            typename mpl::if_<mpl::equal_to<mpl::int_<ElementType::nDim>, mpl::int_<1> >,
            mpl::identity<GeoElement1D<ElementType::nRealDim, ElementType> >,
            mpl::identity<GeoElement0D<ElementType::nRealDim, ElementType> > >::type>::type>::type::type element_type;

    /**
     * multi-indexed element container
     */
    typedef multi_index::multi_index_container<
        element_type,
        multi_index::indexed_by<
            //multi_index::random_access<>,
            // sort by less<int> on id() + pid()
            multi_index::ordered_unique<
                multi_index::composite_key<element_type,
                                           multi_index::const_mem_fun<element_type,
                                                                      rank_type,
                                                                      &element_type::processId>,
                                           multi_index::const_mem_fun<element_type,
                                                                      size_type,
                                                                      &element_type::id> > >,
            // sort by less<int> on marker
            multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_marker>,
                                            multi_index::composite_key<element_type,
                                                                       multi_index::const_mem_fun<element_type,
                                                                                                  Marker1 const&,
                                                                                                  &element_type::marker>,
                                                                       multi_index::const_mem_fun<element_type,
                                                                                                  rank_type,
                                                                                                  &element_type::processId> > >,
            // sort by less<int> on marker
            multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_marker2>,
                                            multi_index::composite_key<element_type,
                                                                       multi_index::const_mem_fun<element_type,
                                                                                                  Marker2 const&,
                                                                                                  &element_type::marker2>,
                                                                       multi_index::const_mem_fun<element_type,
                                                                                                  rank_type,
                                                                                                  &element_type::processId> > >,

            // sort by less<int> on marker
            multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_marker3>,
                                            multi_index::composite_key<element_type,
                                                                       multi_index::const_mem_fun<element_type,
                                                                                                  Marker3 const&,
                                                                                                  &element_type::marker3>,
                                                                       multi_index::const_mem_fun<element_type,
                                                                                                  rank_type,
                                                                                                  &element_type::processId> > >,

            // sort by less<int> on boundary
            multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_location>,
                                            multi_index::composite_key<element_type,
                                                                       multi_index::const_mem_fun<element_type,
                                                                                                  rank_type,
                                                                                                  &element_type::processId>,
                                                                       multi_index::const_mem_fun<element_type,
                                                                                                  bool,
                                                                                                  &element_type::isOnBoundary>,
                                                                       multi_index::const_mem_fun<element_type,
                                                                                                  uint16_type,
                                                                                                  &element_type::boundaryEntityDimension> > >,


            // sort by less<int> on processId
            multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_pid>,
                                            multi_index::const_mem_fun<element_type,
                                                                       rank_type,
                                                                       &element_type::processId> >,

            // sort by less<int> on processId
            multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_ghostcell>,
                                            multi_index::const_mem_fun<element_type,
                                                                       bool,
                                                                       &element_type::isGhostCell> >


            > > elements_type;


    typedef typename elements_type::iterator element_iterator;
    typedef typename elements_type::const_iterator element_const_iterator;

    // marker
    typedef typename elements_type::template index<Feel::detail::by_marker>::type marker_elements;
    typedef typename marker_elements::iterator marker_element_iterator;
    typedef typename marker_elements::const_iterator marker_element_const_iterator;

    // marker2
    typedef typename elements_type::template index<Feel::detail::by_marker2>::type marker2_elements;
    typedef typename marker2_elements::iterator marker2_element_iterator;
    typedef typename marker2_elements::const_iterator marker2_element_const_iterator;

    // marker3
    typedef typename elements_type::template index<Feel::detail::by_marker3>::type marker3_elements;
    typedef typename marker3_elements::iterator marker3_element_iterator;
    typedef typename marker3_elements::const_iterator marker3_element_const_iterator;

    typedef typename elements_type::template index<Feel::detail::by_pid>::type pid_elements;
    typedef typename pid_elements::iterator pid_element_iterator;
    typedef typename pid_elements::const_iterator pid_element_const_iterator;


    typedef typename elements_type::template index<Feel::detail::by_location>::type location_elements;
    typedef typename location_elements::iterator location_element_iterator;
    typedef typename location_elements::const_iterator location_element_const_iterator;

    typedef typename elements_type::template index<Feel::detail::by_ghostcell>::type ghostcell_elements;
    typedef typename ghostcell_elements::iterator ghostcell_element_iterator;
    typedef typename ghostcell_elements::const_iterator ghostcell_element_const_iterator;

    typedef std::map<int, size_type> parts_map_type;
    typedef typename parts_map_type::const_iterator parts_const_iterator_type;

    /// \cond disabled
    struct update_element_neighbor_type
    {
        update_element_neighbor_type( rank_type n, size_type id )
            :
            M_pos_neigh( n ),
            M_neigh_id( id )
        {}

        void operator()( element_type& e )
        {
            e.setNeighbor( M_pos_neigh, M_neigh_id );
        }

    private:
        rank_type M_pos_neigh;
        size_type M_neigh_id;
    };

    /**
     * @class ElementUpdatePoint
     * @brief update point data structure in element
     *
     * this structure is to be used by \p modify from
     * \c boost::multi_index
     *
     */
    struct ElementUpdatePoint
    {
        /**
         * @param index local index of the point
         * @param pt point to update
         */
        ElementUpdatePoint( uint16_type index, typename element_type::PointType const& pt )
            :
            M_index( index ),
            M_pt( pt )
        {}

        /**
         * update the \p M_index point info in element \p e using
         * \p M_pt
         */
        void operator()( element_type& e )
        {
            e.setPoint( M_index, M_pt );
        }
    private:
        uint16_type M_index;
        typename element_type::PointType const& M_pt;

    };
    /**
     * @class ElementConnectPointToElement
     * @brief connect point to element
     *
     */
    struct ElementConnectPointToElement
    {
        void operator()( element_type& e )
        {
            for ( int i = 0; i < e.numPoints; ++i )
                e.point( i ).addElement( e.id() );
        }
    };

    /**
     * @class ElementConnectPointToElement
     * @brief connect point to element
     *
     */
    struct ElementGhostConnectPointToElement
    {
        void operator()( element_type& e )
        {
            for ( int i = 0; i < e.numPoints; ++i )
            {
                e.point( i ).addElementGhost( e.processId(),e.id() );
                // only if point is on interprocess
                if ( e.point( i ).processId()!=invalid_rank_type_value )
                    e.point( i ).addNeighborPartitionId( e.processId() );
            }
        }
    };

    /**
     * @class ElementGhostConnectEdgeToElement
     * @brief connect edge to element
     *
     */
    struct ElementGhostConnectEdgeToElement
    {
        void operator()( element_type& e )
        {
            for ( int i = 0; i < e.numEdges; ++i )
                Feel::detail::updateElementGhostConnectEdgeToElement(e,i,mpl::int_<element_type::nDim>());
        }
    };


    /// \endcond

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Elements( WorldComm const& worldComm = Environment::worldComm() )
        :
        M_worldCommElements(worldComm),
        M_elements()
    {}

    Elements( Elements const & f )
        :
        M_worldCommElements( f.worldCommElements() ),
        M_elements( f.M_elements )
    {}

    virtual ~Elements()
    {
        this->clear();
    }

    void clear()
        {
            VLOG(1) << "deleting elements...\n";
            M_elements.clear();
        }
    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * copy operator
     */
    Elements& operator=( Elements const& e )
    {
        if ( this != &e )
        {
            M_worldCommElements = e.M_worldCommElements;
            M_elements = e.M_elements;
        }

        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the elements container
     */
    elements_type const& elements() const
    {
        return M_elements;
    }

    /**
     * \return the elements container
     */
    elements_type &      elements()
    {
        return M_elements;
    }

    /**
     * \return \p true if the container is empty, \p false otherwise
     */
    virtual bool isEmpty() const
    {
        return M_elements.empty();
    }

    WorldComm const& worldCommElements() const
    {
        return M_worldCommElements;
    }


    bool isBoundaryElement( element_type const & e ) const
    {
        return elementIterator(e.id())->isOnBoundary();
    }
    bool isBoundaryElement( size_type const & id ) const
    {
        return elementIterator(id)->isOnBoundary();
    }

    element_iterator elementIterator( size_type i ) const
    {
        return  M_elements.template get<0>().find( boost::make_tuple( this->worldCommElements().localRank(), i ) );
    };

    element_iterator elementIterator( size_type i, rank_type p ) const
    {
        return  M_elements.template get<0>().find( boost::make_tuple( p, i ) );
    };

    element_iterator elementIterator( element_type const& elt ) const
    {
        return elementIterator( elt.id(), elt.processId() );
    };

    element_type const& element( size_type i ) const
    {
        return *M_elements.template get<0>().find( boost::make_tuple( this->worldCommElements().localRank(), i ) );
    };

    element_type const& element( size_type i, rank_type p ) const
    {
        return *M_elements.template get<0>().find( boost::make_tuple( p, i ) );
    };

    /**
     * \return \c true if element with id \p i is found, \c false otherwise
     */
    bool hasElement( size_type i, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<0>().find( boost::make_tuple( part, i ) ) !=
               M_elements.template get<0>().end();
    }

    element_iterator beginElement()
    {
        return M_elements.begin();
    }
    element_const_iterator beginElement() const
    {
        return M_elements.begin();
    }
    element_iterator endElement()
    {
        return M_elements.end();
    }
    element_const_iterator endElement() const
    {
        return M_elements.end();
    }


    parts_const_iterator_type beginParts() const
    {
        return M_parts.begin();
    }
    parts_const_iterator_type endParts() const
    {
        return M_parts.end();
    }

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with marker \p m on processor \p p
     */
    std::pair<element_iterator, element_iterator>
    elementsRange()
    {
        return std::make_pair( M_elements.begin(), M_elements.end() );
    }

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with marker \p m on processor \p p
     */
    std::pair<element_const_iterator, element_const_iterator>
    elementsRange() const
    {
        return std::make_pair( M_elements.begin(), M_elements.end() );
    }

    element_iterator beginElementWithId( size_type m )
    {
        return M_elements.template get<0>().lower_bound( boost::make_tuple( this->worldCommElements().localRank(), m ) );
    }
    element_const_iterator beginElementWithId( size_type m ) const
    {
        return M_elements.template get<0>().lower_bound( boost::make_tuple( this->worldCommElements().localRank(), m ) );
    }
    element_iterator endElementWithId( size_type m )
    {
        return M_elements.template get<0>().upper_bound( boost::make_tuple( this->worldCommElements().localRank(), m ) );
    }
    element_const_iterator endElementWithId( size_type m ) const
    {
        return M_elements.template get<0>().upper_bound( boost::make_tuple( this->worldCommElements().localRank(), m ) );
    }

    /**
     * \return the iterator \c begin over the elements with \c Marker1 \p m
     */
    marker_element_const_iterator beginElementWithMarker( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_marker>().equal_range( boost::make_tuple( Marker1( m ), part ) ).first;
    }

    /**
     * \return the iterator \c begin over the elements with \c Marker2 \p m
     */
    marker2_element_const_iterator beginElementWithMarker2( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_marker2>().equal_range( boost::make_tuple( Marker2( m ), part ) ).first;
    }

    /**
     * \return the iterator \c begin over the elements with \c Marker3 \p m
     */
    marker3_element_const_iterator beginElementWithMarker3( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_marker3>().equal_range( boost::make_tuple( Marker3( m ), part ) ).first;
    }

    /**
     * \return the iterator \c end over the elements with \c Marker1 \p m
     */
    marker_element_const_iterator endElementWithMarker( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_marker>().equal_range( boost::make_tuple( Marker1( m ), part ) ).second;
    }

    /**
     * \return the iterator \c end over the elements with \c Marker2 \p m
     */
    marker2_element_const_iterator endElementWithMarker2( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_marker2>().equal_range( boost::make_tuple( Marker2( m ), part ) ).second;
    }

    /**
     * \return the iterator \c end over the elements with \c Marker3 \p m
     */
    marker3_element_const_iterator endElementWithMarker3( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_marker3>().equal_range( boost::make_tuple( Marker3( m ), part ) ).second;
    }

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker1 \p m on processor \p p
     */
    std::pair<marker_element_const_iterator, marker_element_const_iterator>
    elementsWithMarker( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_marker>().equal_range( boost::make_tuple( Marker1( m ), part ) );
    }


    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker2 \p m on processor \p p
     */
    std::pair<marker2_element_const_iterator, marker2_element_const_iterator>
    elementsWithMarker2( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_marker2>().equal_range( boost::make_tuple( Marker2( m ), part ) );
    }

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker3 \p m on processor \p p
     */
    std::pair<marker3_element_const_iterator, marker3_element_const_iterator>
    elementsWithMarker3( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_marker3>().equal_range( boost::make_tuple( Marker3( m ), part ) );
    }

    element_iterator beginElementWithProcessId( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<0>().lower_bound( boost::make_tuple( part ) );
    }
    element_const_iterator beginElementWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<0>().lower_bound( boost::make_tuple( part ) );
    }
    element_iterator endElementWithProcessId( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<0>().upper_bound( boost::make_tuple( part ) );
    }
    element_const_iterator endElementWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<0>().upper_bound( boost::make_tuple( part ) );
    }

    std::pair<element_const_iterator, element_const_iterator>
    elementsWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<0>().equal_range( boost::make_tuple( part ) );
    }

    std::pair<element_iterator, element_iterator>
    elementsWithProcessId( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<0>().equal_range( boost::make_tuple( part ) );
    }

    /**
     * get the elements container by id
     *
     *
     * @return the element container by id
     */
    typename elements_type::template nth_index<0>::type &
    elementsById()
    {
        return M_elements.template get<0>();
    }

    /**
     * get the elements container by id
     *
     *
     * @return the element container by id
     */
    typename elements_type::template nth_index<0>::type const&
    elementsById() const
    {
        return M_elements.template get<0>();
    }

    /**
     * get the elements container using the \c Marker1 view
     *
     *
     * @return the element container using \c Marker1 view
     */
    marker_elements &
    elementsByMarker()
    {
        return M_elements.template get<Feel::detail::by_marker>();
    }

    /**
     * get the elements container using the \c Marker2 view
     *
     *
     * @return the element container using \c Marker2 view
     */
    marker2_elements &
    elementsByMarker2()
    {
        return M_elements.template get<Feel::detail::by_marker2>();
    }

    /**
     * get the elements container using the \c Marker3 view
     *
     *
     * @return the element container using \c Marker3 view
     */
    marker3_elements &
    elementsByMarker3()
    {
        return M_elements.template get<Feel::detail::by_marker3>();
    }

    /**
     * get the elements container using the \c Marker1 view
     *
     *
     * @return the element container using \c Marker1 view
     */
    marker_elements const&
    elementsByMarker() const
    {
        return M_elements.template get<Feel::detail::by_marker>();
    }

    /**
     * get the elements container using the \c Marker2 view
     *
     *
     * @return the element container using \c Marker2 view
     */
    marker2_elements const&
    elementsByMarker2() const
    {
        return M_elements.template get<Feel::detail::by_marker2>();
    }

    /**
     * get the elements container using the \c Marker3 view
     *
     *
     * @return the element container using \c Marker3 view
     */
    marker3_elements const&
    elementsByMarker3() const
    {
        return M_elements.template get<Feel::detail::by_marker3>();
    }

    /**
     * get the elements container using the process id view
     *
     *
     * @return the element container using process id view
     */
    pid_elements &
    elementsByProcessId()
    {
        return M_elements.template get<Feel::detail::by_pid>();
    }

    /**
     * get the elements container using the process id view
     *
     *
     * @return the element container using marker view
     */
    pid_elements const&
    elementsByProcessId() const
    {
        return M_elements.template get<Feel::detail::by_pid>();
    }

    /**
     * \return the range of iterator \c (begin,end) over the boundary
     *  element on processor \p p
     */
    std::pair<location_element_const_iterator, location_element_const_iterator>
    boundaryElements( uint16_type entity_min_dim, uint16_type entity_max_dim, rank_type p = invalid_rank_type_value  ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        auto lower = M_elements.template get<Feel::detail::by_location>().lower_bound( boost::make_tuple( part, bool(ON_BOUNDARY), entity_min_dim ) );
        auto upper = M_elements.template get<Feel::detail::by_location>().upper_bound( boost::make_tuple( part, bool(ON_BOUNDARY), entity_max_dim ) );
        return std::make_pair( lower, upper );

    }

    /**
     * \return the range of iterator \c (begin,end) over the boundary
     *  element on processor \p p
     */
    std::pair<location_element_const_iterator, location_element_const_iterator>
    boundaryElements( rank_type p = invalid_rank_type_value  ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return boundaryElements( 0, 2, part );
        //auto lower = boost::make_tuple( this->worldCommElements().localRank(), bool(ON_BOUNDARY), 0);
        //auto upper = boost::make_tuple( this->worldCommElements().localRank(), bool(ON_BOUNDARY), 2);
        //return M_elements.template get<Feel::detail::by_location>().range( lower, upper );

    }


    /**
     * \return the range of iterator \c (begin,end) over the internal
     *  element on processor \p p
     */
    std::pair<location_element_const_iterator, location_element_const_iterator>
    internalElements( rank_type p = invalid_rank_type_value  ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        auto lower = M_elements.template get<Feel::detail::by_location>().lower_bound( boost::make_tuple( part, bool(INTERNAL), 0 ) );
        auto upper = M_elements.template get<Feel::detail::by_location>().upper_bound( boost::make_tuple( part, bool(INTERNAL), invalid_uint16_type_value ) );
        return std::make_pair( lower, upper );
    }


    /**
     * get the elements container using the location view
     *
     *
     * @return the element container using location view
     */
    location_elements &
    elementsByLocation()
    {
        return M_elements.template get<Feel::detail::by_location>();
    }

    /**
     * get the elements container using the location view
     *
     *
     * @return the element container using location view
     */
    location_elements const&
    elementsByLocation() const
    {
        return M_elements.template get<Feel::detail::by_location>();
    }

    /**
     * get the begin() iterator on all the internal elements
     *
     * @return the begin() iterator on all the internal elements
     */
    location_element_iterator beginInternalElement( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_location>().equal_range( boost::make_tuple( part, INTERNAL, invalid_uint16_type_value ) ).first;
    }
    /**
     * get the end() iterator on all the internal elements
     *
     * @return the end() iterator on all the internal elements
     */
    location_element_iterator endInternalElement( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_location>().equal_range( boost::make_tuple( part, INTERNAL, invalid_uint16_type_value ) ).second;
    }

    /**
     * get the begin() iterator on all the internal elements
     *
     * @return the begin() iterator on all the internal elements
     */
    location_element_const_iterator beginInternalElement( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_location>().equal_range( boost::make_tuple( part, INTERNAL, invalid_uint16_type_value ) ).first;
    }

    /**
     * get the end() iterator on all the internal elements
     *
     * @return the end() iterator on all the internal elements
     */
    location_element_const_iterator endInternalElement( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_location>().equal_range( boost::make_tuple( part, INTERNAL, invalid_uint16_type_value ) ).second;
    }

    /**
     * get the begin() iterator on all the boundary elements
     *
     * @return the begin() iterator on all the boundary elements
     */
    location_element_iterator beginElementOnBoundary( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_location>().lower_bound( boost::make_tuple( part, ON_BOUNDARY, 0 ) );
    }
    /**
     * get the end() iterator on all the boundary elements
     *
     * @return the end() iterator on all the boundary elements
     */
    location_element_iterator endElementOnBoundary( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_location>().upper_bound( boost::make_tuple( part, ON_BOUNDARY, 2 ) );
    }

    /**
     * get the begin() iterator on all the boundary elements
     *
     * @return the begin() iterator on all the boundary elements
     */
    location_element_const_iterator beginElementOnBoundary( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_location>().lower_bound( boost::make_tuple( part, ON_BOUNDARY, 0 ) );
    }

    /**
     * get the end() iterator on all the boundary elements
     *
     * @return the end() iterator on all the boundary elements
     */
    location_element_const_iterator endElementOnBoundary( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<Feel::detail::by_location>().upper_bound( boost::make_tuple( part, ON_BOUNDARY, 2 ) );
    }

    /**
     * get the begin() iterator on all ghost elements
     *
     * @return the begin() iterator on all ghost elements
     */
    ghostcell_element_iterator beginGhostElement()
    {
        //return M_elements.template get<Feel::detail::by_ghostcell>().equal_range(boost::make_tuple(true)).first;
        return M_elements.template get<Feel::detail::by_ghostcell>().equal_range( true ).first;
        //return M_elements.template get<Feel::detail::by_ghostcell>().begin();
    }

    /**
     * get the end() iterator on all ghost elements
     *
     * @return the end() iterator on all ghost elements
     */
    ghostcell_element_iterator endGhostElement()
    {
        //return M_elements.template get<Feel::detail::by_ghostcell>().equal_range(boost::make_tuple(true)).second;
        return M_elements.template get<Feel::detail::by_ghostcell>().equal_range( true ).second;
        //return M_elements.template get<Feel::detail::by_ghostcell>().end();
    }

    /**
     * get the begin() iterator on all ghost elements
     *
     * @return the begin() iterator on all ghost elements
     */
    ghostcell_element_const_iterator beginGhostElement() const
    {
        //return M_elements.template get<Feel::detail::by_ghostcell>().equal_range(boost::make_tuple(true)).first;
        return M_elements.template get<Feel::detail::by_ghostcell>().equal_range( true ).first;
        //return M_elements.template get<Feel::detail::by_ghostcell>().begin();
    }

    /**
     * get the end() iterator on all ghost elements
     *
     * @return the end() iterator on all ghost elements
     */
    ghostcell_element_const_iterator endGhostElement() const
    {
        //return M_elements.template get<Feel::detail::by_ghostcell>().equal_range(boost::make_tuple(true)).second;
        return M_elements.template get<Feel::detail::by_ghostcell>().equal_range( true ).second;
        //return M_elements.template get<Feel::detail::by_ghostcell>().end();
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
     * add a new element in the mesh
     * @param f a new point
     * @return the new point from the list
     */
    element_type const& addElement( element_type& f )
    {
        M_parts[f.marker().value()]++;
        f.setId( M_elements.size() );
        return *M_elements.insert( f ).first;
        //M_elements.push_back( f );
        //return M_elements.back();

    }

    template<typename ElementVecType>
    void updateMarker2( ElementVecType const& evec )
    {
        EntityProcessType entityProcess = (evec.functionSpace()->dof()->buildDofTableMPIExtended())? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
        auto rangeElt = Feel::elements( evec.mesh(), entityProcess );
        auto it = rangeElt.template get<1>();
        auto en = rangeElt.template get<2>();

        for ( ; it != en; ++it )
            M_elements.modify( this->elementIterator( boost::unwrap_ref(*it).id(), boost::unwrap_ref(*it).processId() ), [&evec]( element_type& e )
        {
            e.setMarker2(  evec.localToGlobal( e.id(), 0, 0 ) );
        } );
    }

    template<typename ElementVecType>
    void updateMarker3( ElementVecType const& evec )
    {
        EntityProcessType entityProcess = (evec.functionSpace()->dof()->buildDofTableMPIExtended())? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
        auto rangeElt = Feel::elements( evec.mesh(), entityProcess );
        auto it = rangeElt.template get<1>();
        auto en = rangeElt.template get<2>();

        for ( ; it != en; ++it )
            M_elements.modify( this->elementIterator( boost::unwrap_ref(*it).id(), boost::unwrap_ref(*it).processId() ), [&evec]( element_type& e )
        {
            e.setMarker3(  evec.localToGlobal( e.id(), 0, 0 ) );
        } );
    }

    template<typename IteratorRange>
    void updateMarker2WithRangeElements( IteratorRange const& range, flag_type flag )
    {
        typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_range_type;
        iterator_range_type it, en;
        boost::tie( boost::tuples::ignore, it, en ) = range;

        for (  ; it != en; ++it )
            M_elements.modify( this->elementIterator( it->id() ), [&flag]( element_type& e )
        {
            e.setMarker2( flag );
        } );
    }

    template<typename IteratorRange>
    void updateMarker3WithRangeElements( IteratorRange const& range, flag_type flag )
    {
        typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_type;
        iterator_type it, en;
        boost::tie( boost::tuples::ignore, it, en ) = range;

        for ( ; it != en; ++it )
            M_elements.modify( this->elementIterator( it->id() ), [&flag]( element_type& e )
        {
            e.setMarker3( flag );
        } );
    }


    /**
     * update the elements markers by setting them from the face markers associated to the elements
     * warning the marker2 and marker3 must be > 0. if 2 several markers are find in elt, the element take
     * the last find as marker
     */
    void updateMarkersFromFaces()
    {

        auto it = beginElement(), en = endElement();

        for (  ; it != en; ++it )
            M_elements.modify( it,
                             []( element_type& e )
        {
            int newtag2=0, newtag3=0;
            for (uint16_type f=0;f<e.numTopologicalFaces; ++f)
                {
                    int tag2 = e.face(f).marker2().value();
                    int tag3 = e.face(f).marker3().value();
                    if (tag2>0) newtag2=tag2;
                    if (tag3>0) newtag3=tag3;
                }

            if (newtag2>0) e.setMarker2( newtag2 );
            if (newtag3>0) e.setMarker3( newtag3 );

        } );

    }



    void setWorldCommElements( WorldComm const& _worldComm )
    {
        M_worldCommElements = _worldComm;
    }

    //@}

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & M_elements;
            ar & M_parts;
        }


private:
    WorldComm M_worldCommElements;
    elements_type M_elements;
    parts_map_type M_parts;
};
/// \endcond
} // Feel
#endif /* FEELPP_ELEMENTS_HPP */
