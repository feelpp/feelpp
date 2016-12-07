/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-09-03

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007,2008,2009,2010 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2016 Feel++ Consortium

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
        if ( e.edgePtr(i) )
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
#if 0
            multi_index::ordered_unique<multi_index::identity<element_type> >
#else
            multi_index::ordered_unique<
                multi_index::composite_key<element_type,
                                           multi_index::const_mem_fun<element_type,
                                                                      rank_type,
                                                                      &element_type::processId>,
                                           multi_index::const_mem_fun<element_type,
                                                                      size_type,
                                                                      &element_type::id> > >
#endif

#if 0
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
            multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_ghostcell>,
                                            multi_index::const_mem_fun<element_type,
                                                                       bool,
                                                                       &element_type::isGhostCell> >
#endif


            > > elements_type;


    typedef typename elements_type::iterator element_iterator;
    typedef typename elements_type::const_iterator element_const_iterator;

    typedef std::vector<boost::reference_wrapper<element_type const> > elements_reference_wrapper_type;
    typedef std::shared_ptr<elements_reference_wrapper_type> elements_reference_wrapper_ptrtype;
    typedef typename elements_reference_wrapper_type::iterator element_reference_wrapper_iterator;
    typedef typename elements_reference_wrapper_type::const_iterator element_reference_wrapper_const_iterator;

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
        // return M_elements.find( element_type( i ) );
    };

    element_iterator elementIterator( size_type i, rank_type p ) const
    {
        return  M_elements.template get<0>().find( boost::make_tuple( p, i ) );
        // return M_elements.find( element_type( i ) );
    };

    element_iterator elementIterator( element_type const& elt ) const
    {
        // return elementIterator( elt.id() );//, elt.processId() );
        return elementIterator( elt.id(), elt.processId() );
    };

    element_type const& element( size_type i ) const
    {
        return *M_elements.template get<0>().find( boost::make_tuple( this->worldCommElements().localRank(), i ) );
        // return *M_elements.find( element_type( i ) );
    };

    element_type const& element( size_type i, rank_type p ) const
    {
        return *M_elements.template get<0>().find( boost::make_tuple( p, i ) );
        // return *M_elements.find( element_type( i ) );
    };

    /**
     * \return \c true if element with id \p i is found, \c false otherwise
     */
    bool hasElement( size_type i, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return M_elements.template get<0>().find( boost::make_tuple( part, i ) ) !=
               M_elements.template get<0>().end();
#if 0
        return M_elements.find( element_type( i ) ) != M_elements.end();
        if ( M_elements.find( element_type( i ) ) != M_elements.end() )
            return false;
        if ( this->element( i ).processId() != part )
            return false;
        return true;
#endif
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

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Id \p m
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithId( size_type m ) const
        {
            elements_reference_wrapper_ptrtype myelements( new elements_reference_wrapper_type );
            if ( this->hasElement( m ) )
                myelements->push_back( boost::cref( this->element( m ) ) );
            return std::make_tuple( myelements->begin(), myelements->end(), myelements );
        }

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker1 \p m on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithMarker( size_type m, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
            elements_reference_wrapper_ptrtype myelements( new elements_reference_wrapper_type );
            auto it = this->beginElement();
            auto en = this->endElement();
            for ( ; it!=en;++it )
            {
                auto const& elt = *it;
                if ( elt.processId() != part )
                    continue;
                if ( elt.marker().isOff() || elt.marker().value() != m )
                    continue;
                myelements->push_back(boost::cref(elt));
            }
            return std::make_tuple( myelements->begin(), myelements->end(), myelements );
        }
    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker2 \p m on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithMarker2( size_type m, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
            elements_reference_wrapper_ptrtype myelements( new elements_reference_wrapper_type );
            auto it = this->beginElement();
            auto en = this->endElement();
            for ( ; it!=en;++it )
            {
                auto const& elt = *it;
                if ( elt.processId() != part )
                    continue;
                if ( elt.marker2().isOff() || elt.marker2().value() != m )
                    continue;
                myelements->push_back(boost::cref(elt));
            }
            return std::make_tuple( myelements->begin(), myelements->end(), myelements );
        }
    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker3 \p m on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithMarker3( size_type m, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
            elements_reference_wrapper_ptrtype myelements( new elements_reference_wrapper_type );
            auto it = this->beginElement();
            auto en = this->endElement();
            for ( ; it!=en;++it )
            {
                auto const& elt = *it;
                if ( elt.processId() != part )
                    continue;
                if ( elt.marker3().isOff() || elt.marker3().value() != m )
                    continue;
                myelements->push_back(boost::cref(elt));
            }
            return std::make_tuple( myelements->begin(), myelements->end(), myelements );
        }
    /**
     * \return the range of iterator \c (begin,end) over the elements
     * on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        elements_reference_wrapper_ptrtype myelements( new elements_reference_wrapper_type );
        auto it = this->beginElement();
        auto en = this->endElement();
        for ( ; it!=en;++it )
        {
            auto const& elt = *it;
            if ( elt.processId() != part )
                continue;
            myelements->push_back(boost::cref(elt));
        }
        return std::make_tuple( myelements->begin(), myelements->end(), myelements );
    }

    /**
     * \return the first iterator over the elements
     * on processor \p p if exist else return endElement iterator
     */
    element_iterator
    firstElementIteratorWithProcessId( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
            auto it = this->beginElement();
            auto en = this->endElement();
            for ( ; it!=en;++it )
            {
                auto const& elt = *it;
                if ( elt.processId() == part )
                    return it;
            }
            return en;
        }

    /**
     * \return the range of iterator \c (begin,end) over the boundary
     *  element on processor \p p which share a subentity of minDim<= dim <= maxDim on boundary
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    boundaryElements( uint16_type entity_min_dim, uint16_type entity_max_dim, rank_type p = invalid_rank_type_value  ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        elements_reference_wrapper_ptrtype myelements( new elements_reference_wrapper_type );
        auto it = this->beginElement();
        auto en = this->endElement();
        for ( ; it!=en;++it )
        {
            auto const& elt = *it;
            if ( elt.processId() != part )
                continue;
            if ( !elt.isOnBoundary() )
                continue;
            if ( elt.boundaryEntityDimension() < entity_min_dim )
                continue;
            if ( elt.boundaryEntityDimension() > entity_max_dim )
                continue;
            myelements->push_back(boost::cref(elt));
        }
        return std::make_tuple( myelements->begin(), myelements->end(), myelements );
    }

    /**
     * \return the range of iterator \c (begin,end) over the boundary
     *  element on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    boundaryElements( rank_type p = invalid_rank_type_value  ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        return boundaryElements( 0, 2, part );
    }

    /**
     * \return the range of iterator \c (begin,end) over the internal
     *  element on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    internalElements( rank_type p = invalid_rank_type_value  ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
        elements_reference_wrapper_ptrtype myelements( new elements_reference_wrapper_type );
        auto it = this->beginElement();
        auto en = this->endElement();
        for ( ; it!=en;++it )
        {
            auto const& elt = *it;
            if ( elt.processId() != part )
                continue;
            if ( !elt.isInternal() )
                continue;
            myelements->push_back(boost::cref(elt));
        }
        return std::make_tuple( myelements->begin(), myelements->end(), myelements );
    }


    /**
     * \return the range of iterator \c (begin,end) over the internal
     *  element on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    ghostElements() const
    {
        elements_reference_wrapper_ptrtype myelements( new elements_reference_wrapper_type );
        auto it = this->beginElement();
        auto en = this->endElement();
        for ( ; it!=en;++it )
        {
            auto const& elt = *it;
            if ( !elt.isGhostCell() )
                continue;
            myelements->push_back(boost::cref(elt));
        }
        return std::make_tuple( myelements->begin(), myelements->end(), myelements );
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
    element_type const& addElement( element_type& f, bool setid = true )
    {
        M_parts[f.marker().value()]++;
        if ( setid )
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
