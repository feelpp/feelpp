/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-09-03
 */
#ifndef __elements_H
#define __elements_H 1


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
/*!
  \class Elements
  \brief Elements container class

  @author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
    uint16_type,
    &element_type::processId>,
    multi_index::const_mem_fun<element_type,
    size_type,
    &element_type::id> > >,
    // sort by less<int> on marker
    multi_index::ordered_non_unique<multi_index::tag<detail::by_marker>,
    multi_index::composite_key<element_type,
    multi_index::const_mem_fun<element_type,
    Marker1 const&,
    &element_type::marker>,
    multi_index::const_mem_fun<element_type,
    uint16_type,
    &element_type::processId> > >,
    // sort by less<int> on marker
    multi_index::ordered_non_unique<multi_index::tag<detail::by_marker2>,
    multi_index::composite_key<element_type,
    multi_index::const_mem_fun<element_type,
    Marker2 const&,
    &element_type::marker2>,
    multi_index::const_mem_fun<element_type,
    uint16_type,
    &element_type::processId> > >,

    // sort by less<int> on marker
    multi_index::ordered_non_unique<multi_index::tag<detail::by_marker3>,
    multi_index::composite_key<element_type,
    multi_index::const_mem_fun<element_type,
    Marker3 const&,
    &element_type::marker3>,
    multi_index::const_mem_fun<element_type,
    uint16_type,
    &element_type::processId> > >,

    // sort by less<int> on boundary
    multi_index::ordered_non_unique<multi_index::tag<detail::by_location>,
    multi_index::composite_key<element_type,
    multi_index::const_mem_fun<element_type,
    bool,
    &element_type::isOnBoundary>,
    multi_index::const_mem_fun<element_type,
    uint16_type,
    &element_type::processId> > >,

    // sort by less<int> on processId
    multi_index::ordered_non_unique<multi_index::tag<detail::by_pid>,
    multi_index::const_mem_fun<element_type,
    uint16_type,
    &element_type::processId> >,

    // sort by less<int> on processId
    multi_index::ordered_non_unique<multi_index::tag<detail::by_ghostcell>,
    multi_index::const_mem_fun<element_type,
    bool,
    &element_type::isGhostCell> >


    > > elements_type;


    typedef typename elements_type::iterator element_iterator;
    typedef typename elements_type::const_iterator element_const_iterator;

    // marker
    typedef typename elements_type::template index<detail::by_marker>::type marker_elements;
    typedef typename marker_elements::iterator marker_element_iterator;
    typedef typename marker_elements::const_iterator marker_element_const_iterator;

    // marker2
    typedef typename elements_type::template index<detail::by_marker2>::type marker2_elements;
    typedef typename marker2_elements::iterator marker2_element_iterator;
    typedef typename marker2_elements::const_iterator marker2_element_const_iterator;

    // marker3
    typedef typename elements_type::template index<detail::by_marker3>::type marker3_elements;
    typedef typename marker3_elements::iterator marker3_element_iterator;
    typedef typename marker3_elements::const_iterator marker3_element_const_iterator;

    typedef typename elements_type::template index<detail::by_pid>::type pid_elements;
    typedef typename pid_elements::iterator pid_element_iterator;
    typedef typename pid_elements::const_iterator pid_element_const_iterator;


    typedef typename elements_type::template index<detail::by_location>::type location_elements;
    typedef typename location_elements::iterator location_element_iterator;
    typedef typename location_elements::const_iterator location_element_const_iterator;

    typedef typename elements_type::template index<detail::by_ghostcell>::type ghostcell_elements;
    typedef typename ghostcell_elements::iterator ghostcell_element_iterator;
    typedef typename ghostcell_elements::const_iterator ghostcell_element_const_iterator;

    typedef std::map<int, size_type> parts_map_type;
    typedef typename parts_map_type::const_iterator parts_const_iterator_type;

    /// \cond disabled
    struct update_element_neighbor_type
    {
        update_element_neighbor_type( uint16_type n, size_type id )
            :
            _M_pos_neigh( n ),
            _M_neigh_id( id )
        {}

        void operator()( element_type& e )
        {
            e.setNeighbor( _M_pos_neigh, _M_neigh_id );
        }

    private:
        uint16_type _M_pos_neigh;
        size_type _M_neigh_id;
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
            _M_index( index ),
            _M_pt( pt )
        {}

        /**
         * update the \p _M_index point info in element \p e using
         * \p _M_pt
         */
        void operator()( element_type& e )
        {
            e.setPoint( _M_index, _M_pt );
        }
    private:
        uint16_type _M_index;
        typename element_type::PointType const& _M_pt;

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
    /// \endcond

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Elements( WorldComm const& worldComm = Environment::worldComm() )
        :
        _M_worldCommElements(worldComm),
        _M_elements()
    {}

    Elements( Elements const & f )
        :
        _M_worldCommElements( f.worldCommElements() ),
        _M_elements( f._M_elements )
    {}

    virtual ~Elements()
    {}

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
            _M_worldCommElements = e._M_worldCommElements;
            _M_elements = e._M_elements;
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
        return _M_elements;
    }

    /**
     * \return the elements container
     */
    elements_type &      elements()
    {
        return _M_elements;
    }

    /**
     * \return \p true if the container is empty, \p false otherwise
     */
    virtual bool isEmpty() const
    {
        return _M_elements.empty();
    }

    WorldComm const& worldCommElements() const
    {
        return _M_worldCommElements;
    }


    bool isBoundaryElement( element_type const & e ) const
    {
        return _M_elements.find( e )->isOnBoundary();
    }
    bool isBoundaryElement( size_type const & id ) const
    {
        return _M_elements.find( element_type( id ) )->isOnBoundary();
    }

    element_iterator elementIterator( size_type i ) const
    {
        return  _M_elements.template get<0>().find( boost::make_tuple( this->worldCommElements().localRank(), i ) );
    };

    element_iterator elementIterator( size_type i, size_type p ) const
    {
        return  _M_elements.template get<0>().find( boost::make_tuple( p, i ) );
    };

    element_type const& element( size_type i ) const
    {
        return *_M_elements.template get<0>().find( boost::make_tuple( this->worldCommElements().localRank(), i ) );
    };

    element_type const& element( size_type i, size_type p ) const
    {
        return *_M_elements.template get<0>().find( boost::make_tuple( p, i ) );
    };

    /**
     * \return \c true if element with id \p i is found, \c false otherwise
     */
    bool hasElement( size_type i ) const
    {
        return _M_elements.template get<0>().find( boost::make_tuple( this->worldCommElements().localRank(), i ) ) !=
               _M_elements.template get<0>().end();
    }

    element_iterator beginElement()
    {
        return _M_elements.begin();
    }
    element_const_iterator beginElement() const
    {
        return _M_elements.begin();
    }
    element_iterator endElement()
    {
        return _M_elements.end();
    }
    element_const_iterator endElement() const
    {
        return _M_elements.end();
    }


    parts_const_iterator_type beginParts() const
    {
        return _M_parts.begin();
    }
    parts_const_iterator_type endParts() const
    {
        return _M_parts.end();
    }

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with marker \p m on processor \p p
     */
    std::pair<element_iterator, element_iterator>
    elementsRange()
    {
        return std::make_pair( _M_elements.begin(), _M_elements.end() );
    }

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with marker \p m on processor \p p
     */
    std::pair<element_const_iterator, element_const_iterator>
    elementsRange() const
    {
        return std::make_pair( _M_elements.begin(), _M_elements.end() );
    }

    element_iterator beginElementWithId( size_type m )
    {
        return _M_elements.template get<0>().lower_bound( boost::make_tuple( this->worldCommElements().localRank(), m ) );
    }
    element_const_iterator beginElementWithId( size_type m ) const
    {
        return _M_elements.template get<0>().lower_bound( boost::make_tuple( this->worldCommElements().localRank(), m ) );
    }
    element_iterator endElementWithId( size_type m )
    {
        return _M_elements.template get<0>().upper_bound( boost::make_tuple( this->worldCommElements().localRank(), m ) );
    }
    element_const_iterator endElementWithId( size_type m ) const
    {
        return _M_elements.template get<0>().upper_bound( boost::make_tuple( this->worldCommElements().localRank(), m ) );
    }

    /**
     * \return the iterator \c begin over the elements with \c Marker1 \p m
     */
    marker_element_const_iterator beginElementWithMarker( size_type m ) const
    {
        return _M_elements.template get<detail::by_marker>().equal_range( boost::make_tuple( Marker1( m ), this->worldCommElements().localRank() ) ).first;
    }

    /**
     * \return the iterator \c begin over the elements with \c Marker2 \p m
     */
    marker2_element_const_iterator beginElementWithMarker2( size_type m ) const
    {
        return _M_elements.template get<detail::by_marker2>().equal_range( boost::make_tuple( Marker2( m ), this->worldCommElements().localRank() ) ).first;
    }

    /**
     * \return the iterator \c begin over the elements with \c Marker3 \p m
     */
    marker3_element_const_iterator beginElementWithMarker3( size_type m ) const
    {
        return _M_elements.template get<detail::by_marker3>().equal_range( boost::make_tuple( Marker3( m ), this->worldCommElements().localRank() ) ).first;
    }

    /**
     * \return the iterator \c end over the elements with \c Marker1 \p m
     */
    marker_element_const_iterator endElementWithMarker( size_type m ) const
    {
        return _M_elements.template get<detail::by_marker>().equal_range( boost::make_tuple( Marker1( m ), this->worldCommElements().localRank() ) ).second;
    }

    /**
     * \return the iterator \c end over the elements with \c Marker2 \p m
     */
    marker2_element_const_iterator endElementWithMarker2( size_type m ) const
    {
        return _M_elements.template get<detail::by_marker2>().equal_range( boost::make_tuple( Marker2( m ), this->worldCommElements().localRank() ) ).second;
    }

    /**
     * \return the iterator \c end over the elements with \c Marker3 \p m
     */
    marker3_element_const_iterator endElementWithMarker3( size_type m ) const
    {
        return _M_elements.template get<detail::by_marker3>().equal_range( boost::make_tuple( Marker3( m ), this->worldCommElements().localRank() ) ).second;
    }

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker1 \p m on processor \p p
     */
    std::pair<marker_element_const_iterator, marker_element_const_iterator>
    elementsWithMarker( size_type m, size_type p ) const
    {
        return _M_elements.template get<detail::by_marker>().equal_range( boost::make_tuple( Marker1( m ), this->worldCommElements().localRank() ) );
    }


    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker2 \p m on processor \p p
     */
    std::pair<marker2_element_const_iterator, marker2_element_const_iterator>
    elementsWithMarker2( size_type m, size_type p ) const
    {
        return _M_elements.template get<detail::by_marker2>().equal_range( boost::make_tuple( Marker2( m ), this->worldCommElements().localRank() ) );
    }

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker3 \p m on processor \p p
     */
    std::pair<marker3_element_const_iterator, marker3_element_const_iterator>
    elementsWithMarker3( size_type m, size_type p ) const
    {
        return _M_elements.template get<detail::by_marker3>().equal_range( boost::make_tuple( Marker3( m ), this->worldCommElements().localRank() ) );
    }

    element_iterator beginElementWithProcessId( size_type m )
    {
        return _M_elements.template get<0>().lower_bound( boost::make_tuple( m ) );
    }
    element_const_iterator beginElementWithProcessId( size_type m ) const
    {
        return _M_elements.template get<0>().lower_bound( boost::make_tuple( m ) );
    }
    element_iterator endElementWithProcessId( size_type m )
    {
        return _M_elements.template get<0>().upper_bound( boost::make_tuple( m ) );
    }
    element_const_iterator endElementWithProcessId( size_type m ) const
    {
        return _M_elements.template get<0>().upper_bound( boost::make_tuple( m ) );
    }

    std::pair<element_const_iterator, element_const_iterator>
    elementsWithProcessId( size_type m ) const
    {
        return _M_elements.template get<0>().equal_range( boost::make_tuple( m ) );
    }

    std::pair<element_iterator, element_iterator>
    elementsWithProcessId( size_type m )
    {
        return _M_elements.template get<0>().equal_range( boost::make_tuple( m ) );
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
        return _M_elements.template get<0>();
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
        return _M_elements.template get<0>();
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
        return _M_elements.template get<detail::by_marker>();
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
        return _M_elements.template get<detail::by_marker2>();
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
        return _M_elements.template get<detail::by_marker3>();
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
        return _M_elements.template get<detail::by_marker>();
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
        return _M_elements.template get<detail::by_marker2>();
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
        return _M_elements.template get<detail::by_marker3>();
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
        return _M_elements.template get<detail::by_pid>();
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
        return _M_elements.template get<detail::by_pid>();
    }

    /**
     * \return the range of iterator \c (begin,end) over the boundary
     *  element on processor \p p
     */
    std::pair<location_element_const_iterator, location_element_const_iterator>
    boundaryElements( size_type p  ) const
    {
        return _M_elements.template get<detail::by_location>().equal_range( boost::make_tuple( ON_BOUNDARY, this->worldCommElements().localRank() ) );
    }

    /**
     * \return the range of iterator \c (begin,end) over the internal
     *  element on processor \p p
     */
    std::pair<location_element_const_iterator, location_element_const_iterator>
    internalElements( size_type p  ) const
    {
        return _M_elements.template get<detail::by_location>().equal_range( boost::make_tuple( INTERNAL, this->worldCommElements().localRank() ) );
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
        return _M_elements.template get<detail::by_location>();
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
        return _M_elements.template get<detail::by_location>();
    }

    /**
     * get the begin() iterator on all the internal elements
     *
     * @return the begin() iterator on all the internal elements
     */
    location_element_iterator beginInternalElement()
    {
        return _M_elements.template get<detail::by_location>().equal_range( boost::make_tuple( INTERNAL,this->worldCommElements().localRank() ) ).first;
    }
    /**
     * get the end() iterator on all the internal elements
     *
     * @return the end() iterator on all the internal elements
     */
    location_element_iterator endInternalElement()
    {
        return _M_elements.template get<detail::by_location>().equal_range( boost::make_tuple( INTERNAL,this->worldCommElements().localRank() ) ).second;
    }

    /**
     * get the begin() iterator on all the internal elements
     *
     * @return the begin() iterator on all the internal elements
     */
    location_element_const_iterator beginInternalElement() const
    {
        return _M_elements.template get<detail::by_location>().equal_range( boost::make_tuple( INTERNAL,this->worldCommElements().localRank() ) ).first;
    }

    /**
     * get the end() iterator on all the internal elements
     *
     * @return the end() iterator on all the internal elements
     */
    location_element_const_iterator endInternalElement() const
    {
        return _M_elements.template get<detail::by_location>().equal_range( boost::make_tuple( INTERNAL,this->worldCommElements().localRank() ) ).second;
    }

    /**
     * get the begin() iterator on all the boundary elements
     *
     * @return the begin() iterator on all the boundary elements
     */
    location_element_iterator beginElementOnBoundary()
    {
        return _M_elements.template get<detail::by_location>().equal_range( boost::make_tuple( ON_BOUNDARY,this->worldCommElements().localRank() ) ).first;
    }
    /**
     * get the end() iterator on all the boundary elements
     *
     * @return the end() iterator on all the boundary elements
     */
    location_element_iterator endElementOnBoundary()
    {
        return _M_elements.template get<detail::by_location>().equal_range( boost::make_tuple( ON_BOUNDARY,this->worldCommElements().localRank() ) ).second;
    }

    /**
     * get the begin() iterator on all the boundary elements
     *
     * @return the begin() iterator on all the boundary elements
     */
    location_element_const_iterator beginElementOnBoundary() const
    {
        return _M_elements.template get<detail::by_location>().equal_range( boost::make_tuple( ON_BOUNDARY,this->worldCommElements().localRank() ) ).first;
    }

    /**
     * get the end() iterator on all the boundary elements
     *
     * @return the end() iterator on all the boundary elements
     */
    location_element_const_iterator endElementOnBoundary() const
    {
        return _M_elements.template get<detail::by_location>().equal_range( boost::make_tuple( ON_BOUNDARY,this->worldCommElements().localRank() ) ).second;
    }

    /**
     * get the begin() iterator on all ghost elements
     *
     * @return the begin() iterator on all ghost elements
     */
    ghostcell_element_iterator beginGhostElement()
    {
        //return _M_elements.template get<detail::by_ghostcell>().equal_range(boost::make_tuple(true)).first;
        return _M_elements.template get<detail::by_ghostcell>().equal_range( true ).first;
        //return _M_elements.template get<detail::by_ghostcell>().begin();
    }

    /**
     * get the end() iterator on all ghost elements
     *
     * @return the end() iterator on all ghost elements
     */
    ghostcell_element_iterator endGhostElement()
    {
        //return _M_elements.template get<detail::by_ghostcell>().equal_range(boost::make_tuple(true)).second;
        return _M_elements.template get<detail::by_ghostcell>().equal_range( true ).second;
        //return _M_elements.template get<detail::by_ghostcell>().end();
    }

    /**
     * get the begin() iterator on all ghost elements
     *
     * @return the begin() iterator on all ghost elements
     */
    ghostcell_element_const_iterator beginGhostElement() const
    {
        //return _M_elements.template get<detail::by_ghostcell>().equal_range(boost::make_tuple(true)).first;
        return _M_elements.template get<detail::by_ghostcell>().equal_range( true ).first;
        //return _M_elements.template get<detail::by_ghostcell>().begin();
    }

    /**
     * get the end() iterator on all ghost elements
     *
     * @return the end() iterator on all ghost elements
     */
    ghostcell_element_const_iterator endGhostElement() const
    {
        //return _M_elements.template get<detail::by_ghostcell>().equal_range(boost::make_tuple(true)).second;
        return _M_elements.template get<detail::by_ghostcell>().equal_range( true ).second;
        //return _M_elements.template get<detail::by_ghostcell>().end();
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
        _M_parts[f.marker().value()]++;
        f.setId( _M_elements.size() );
        return *_M_elements.insert( f ).first;
        //_M_elements.push_back( f );
        //return _M_elements.back();

    }

    /**
     * erase element at position \p position
     *
     * @param position \p position is a valid dereferenceable iterator of the index.
     *
     * @return An iterator pointing to the element immediately
     * following the one that was deleted, or \c end() if no such element
     * exists.
     */
    element_iterator eraseElement( element_iterator position )
    {
        return _M_elements.erase( position );
    }

    template<typename ElementVecType>
    void updateMarker2( ElementVecType const& evec )
    {
        element_iterator it;
        element_iterator en;
        boost::tie( it, en ) = elementsWithProcessId( this->worldCommElements().localRank() );

        for ( ; it != en; ++it )
            _M_elements.modify( it, [&evec]( element_type& e )
        {
            e.setMarker2(  evec.localToGlobal( e.id(), 0, 0 ) );
        } );
    }

    template<typename ElementVecType>
    void updateMarker3( ElementVecType const& evec )
    {
        element_iterator it;
        element_iterator en;
        boost::tie( it, en ) = elementsWithProcessId( this->worldCommElements().localRank() );

        for ( ; it != en; ++it )
            _M_elements.modify( it, [&evec]( element_type& e )
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
            _M_elements.modify( this->elementIterator( it->id() ), [&flag]( element_type& e )
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
            _M_elements.modify( this->elementIterator( it->id() ), [&flag]( element_type& e )
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
            _M_elements.modify( it,
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
        _M_worldCommElements = _worldComm;
    }

    //@}

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & _M_elements;
            ar & _M_parts;
        }


private:
    WorldComm _M_worldCommElements;
    elements_type _M_elements;
    parts_map_type _M_parts;
};
/// \endcond
} // Feel
#endif /* __elements_H */
