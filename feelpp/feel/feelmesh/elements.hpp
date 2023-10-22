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
#ifndef FEELPP_MESH_ELEMENTS_HPP
#define FEELPP_MESH_ELEMENTS_HPP

#include <unordered_map>

#include <feel/feelcore/commobject.hpp>
#include <feel/feelmesh/geoelement.hpp>
#include <feel/feelmesh/filters.hpp>

namespace Feel
{

/// \cond detail

namespace detail
{
    template <typename EltType >
    void
    updateElementGhostConnectEdgeToElement( EltType& e, uint16_type i, mpl::int_<1> /**/)
    {}
    template <typename EltType >
    void
    updateElementGhostConnectEdgeToElement( EltType& e, uint16_type i, mpl::int_<2> /**/)
    {}
    template <typename EltType >
    void
    updateElementGhostConnectEdgeToElement( EltType& e, uint16_type i, mpl::int_<3> /**/)
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
template<typename ElementType, typename T = double, typename IndexT = uint32_type>
class Elements 
{
public:


    /** @name Typedefs
     */
    //@{

    using index_type = IndexT;
    using size_type = index_type;
    
    /**
     * Element type depending on the dimension, @see geoelement.hpp
     * \note Elements have their topological dimension equal to the
     * dimension of the geometric space.
     */
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<ElementType::nDim>, mpl::int_<3> >,
                              mpl::identity<GeoElement3D<ElementType::nRealDim, ElementType, T, IndexT, true> >,
            typename mpl::if_<mpl::equal_to<mpl::int_<ElementType::nDim>, mpl::int_<2> >,
                              mpl::identity<GeoElement2D<ElementType::nRealDim, ElementType, SubFaceOfNone<ElementType::nDim,IndexT>, T, IndexT, true> >,
            typename mpl::if_<mpl::equal_to<mpl::int_<ElementType::nDim>, mpl::int_<1> >,
                              mpl::identity<GeoElement1D<ElementType::nRealDim, ElementType, SubFaceOfNone<ElementType::nDim, IndexT>, T, IndexT, true, true> >,
                              mpl::identity<GeoElement0D<ElementType::nRealDim, SubFaceOfNone<ElementType::nDim, IndexT>/*ElementType*/, T, IndexT> > >::type>::type>::type::type element_type;


    typedef std::unordered_map<size_type,element_type> elements_type;

    typedef typename elements_type::iterator element_iterator;
    typedef typename elements_type::const_iterator element_const_iterator;

    typedef std::vector<boost::reference_wrapper<element_type const> > elements_reference_wrapper_type;
    typedef std::shared_ptr<elements_reference_wrapper_type> elements_reference_wrapper_ptrtype;
    typedef typename elements_reference_wrapper_type::iterator element_reference_wrapper_iterator;
    typedef typename elements_reference_wrapper_type::const_iterator element_reference_wrapper_const_iterator;

    typedef std::vector<boost::reference_wrapper<element_type> > ordered_elements_reference_wrapper_type;
    typedef typename ordered_elements_reference_wrapper_type::iterator ordered_element_reference_wrapper_iterator;
    typedef typename ordered_elements_reference_wrapper_type::const_iterator ordered_element_reference_wrapper_const_iterator;

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
            for ( uint16_type i = 0; i < e.numPoints; ++i )
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
            for ( uint16_type i = 0; i < e.numPoints; ++i )
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
            for ( uint16_type i = 0; i < e.numEdges; ++i )
                Feel::detail::updateElementGhostConnectEdgeToElement(e,i,mpl::int_<element_type::nDim>());
        }
    };


    /// \endcond

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Elements( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        :
        M_worldComm( worldComm ),
        M_elements(),
        M_needToOrderElements( false )
    {}

    Elements( Elements const & f )
        :
        M_worldComm( f.M_worldComm ),
        M_elements( f.M_elements ),
        M_needToOrderElements( false )
    {
        this->buildOrderedElements();
    }

    virtual ~Elements() {}

    void clear()
        {
            DVLOG(1) << "deleting elements...\n";
            M_orderedElements.clear();
            M_elements.clear();
            M_needToOrderElements = false;
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
            M_worldComm = e.M_worldComm;
            M_elements = e.M_elements;
            this->buildOrderedElements();
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

    WorldComm & worldCommElements() 
        {
            return *M_worldComm;
        }
    WorldComm const& worldCommElements() const
    {
        return *M_worldComm;
    }


    bool isBoundaryElement( element_type const & e ) const
    {
        return e.isOnBoundary();
    }
    bool isBoundaryElement( size_type const & id ) const
    {
        auto itFindElt = M_elements.find( id );
        if ( itFindElt == M_elements.end() )
            return false;
        return itFindElt->second.isOnBoundary();
    }

    element_const_iterator elementIterator( size_type i ) const
    {
        return M_elements.find( i );
    };
    element_iterator elementIterator( size_type i )
    {
        return M_elements.find( i );
    };

    FEELPP_DEPRECATED
    element_const_iterator elementIterator( size_type i, rank_type p ) const
    {
        return M_elements.find( i );
    };

    element_const_iterator elementIterator( element_type const& elt ) const
    {
        return elementIterator( elt.id() );
    };
    element_iterator elementIterator( element_type const& elt )
    {
        return elementIterator( elt.id() );
    };

    element_type const& element( size_type i ) const
    {
        auto itFindElt = M_elements.find( i );
        CHECK( itFindElt != M_elements.end() ) << " element " << i << "does not found";
        return itFindElt->second;
    };

    FEELPP_DEPRECATED
    element_type const& element( size_type i, rank_type p ) const
    {
        return this->element( i );
    };

    /**
     * \return \c true if element with id \p i is found, \c false otherwise
     */
    bool hasElement( size_type i ) const
    {
        return M_elements.find( i ) != M_elements.end();
    }

    bool hasElement( size_type i, rank_type p ) const
    {
        auto itFindElt = M_elements.find( i );
        if ( itFindElt == M_elements.end() )
            return false;
        if ( itFindElt->second.processId() != p )
            return false;
        return true;
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

    ordered_element_reference_wrapper_iterator beginOrderedElement()
        {
            return M_orderedElements.begin();
        }
    ordered_element_reference_wrapper_const_iterator beginOrderedElement() const
        {
            return M_orderedElements.begin();
        }
    ordered_element_reference_wrapper_iterator endOrderedElement()
        {
            return M_orderedElements.end();
        }
    ordered_element_reference_wrapper_const_iterator endOrderedElement() const
        {
            return M_orderedElements.end();
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
     * with any \c Marker1 \p on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithMarkerByType( uint16_type markerType, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
            elements_reference_wrapper_ptrtype myelements( new elements_reference_wrapper_type );
            auto it = this->beginOrderedElement();
            auto en = this->endOrderedElement();
            for ( ; it!=en;++it )
            {
                auto const& elt = unwrap_ref( *it );
                if ( elt.processId() != part )
                    continue;
                if ( !elt.hasMarkerType( markerType ) )
                    continue;
                if ( elt.marker( markerType ).isOff() )
                    continue;
                myelements->push_back(boost::cref(elt));
            }
            myelements->shrink_to_fit();
            return std::make_tuple( myelements->begin(), myelements->end(), myelements );
        }
    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker1 \p markerFlags on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
            elements_reference_wrapper_ptrtype myelements( new elements_reference_wrapper_type );
            auto it = this->beginOrderedElement();
            auto en = this->endOrderedElement();
            for ( ; it!=en;++it )
            {
                auto const& elt = unwrap_ref( *it );
                if ( elt.processId() != part )
                    continue;
                if ( !elt.hasMarkerType( markerType ) )
                    continue;
                if ( !elt.marker( markerType ).hasOneOf( markerFlags ) )
                    continue;
                myelements->push_back(boost::cref(elt));
            }
            myelements->shrink_to_fit();
            return std::make_tuple( myelements->begin(), myelements->end(), myelements );
        }
    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker1 \p m on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithMarkerByType( uint16_type markerType, flag_type m, rank_type p = invalid_rank_type_value ) const
        {
            if ( m == invalid_flag_type_value )
                return this->elementsWithMarkerByType( markerType, p );
            else
                return this->elementsWithMarkerByType( markerType, std::set<flag_type>( { m } ), p );

        }

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker1 \p m on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithMarker( flag_type m = invalid_flag_type_value, rank_type p = invalid_rank_type_value ) const
        {
            return this->elementsWithMarkerByType( 1, m, p );
        }
    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker2 \p m on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithMarker2( flag_type m, rank_type p = invalid_rank_type_value ) const
        {
            return this->elementsWithMarkerByType( 2, m, p );
        }
    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker3 \p m on processor \p p
     */
    std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>
    elementsWithMarker3( flag_type m, rank_type p = invalid_rank_type_value ) const
        {
            return this->elementsWithMarkerByType( 3, m, p );
        }


    /**
     * \return the range of iterator \c (begin,end) over the elements
     * with \c Marker1 \p markerFlags on processor \p p
     *  - TheType=0 : only actives elements
     * - TheType=1 : only ghosts elements
     * - TheType=2 : actives and ghosts elements (container splited)
     */
    template <int TheType=0>
    auto
    collectionOfElementsWithMarkerByType( uint16_type markerType,  std::map<int,std::set<flag_type>> const& collectionOfMarkerFlagSet, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;

            std::map<int,elements_reference_wrapper_ptrtype> collectionOfElements, collectionOfGhostElements;
            for ( auto const& [part,markersFlag] : collectionOfMarkerFlagSet )
            {
                collectionOfElements[part].reset( new elements_reference_wrapper_type );
                if constexpr ( TheType == 2 )
                    collectionOfGhostElements[part].reset( new elements_reference_wrapper_type );
            }
            auto it = this->beginOrderedElement();
            auto en = this->endOrderedElement();
            for ( ; it!=en;++it )
            {
                auto const& elt = unwrap_ref( *it );
                bool isActiveElt = elt.processId() == part;
                if constexpr ( TheType == 0 )
                    if ( !isActiveElt )
                        continue;
                if constexpr ( TheType == 1 )
                    if ( isActiveElt )
                        continue;

                if ( !elt.hasMarkerType( markerType ) )
                    continue;
                if ( elt.marker( markerType ).isOff() )
                    continue;
                for ( auto const& [part,markersFlag] : collectionOfMarkerFlagSet )
                {
                    if ( !elt.marker( markerType ).hasOneOf( markersFlag ) )
                        continue;

                    if constexpr ( TheType == 0 || TheType == 1 )
                    {
                        collectionOfElements[part]->push_back(boost::cref(elt));
                    }
                    else
                    {
                        if ( isActiveElt )
                            collectionOfElements[part]->push_back(boost::cref(elt));
                        else
                            collectionOfGhostElements[part]->push_back(boost::cref(elt));
                    }
                    break;
                }
            }

            if constexpr ( TheType == 0 || TheType == 1 )
            {
                std::map<int, std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype> > collectionOfRangeElement;
                for ( auto & [part,myelements] : collectionOfElements )
                {
                    myelements->shrink_to_fit();
                    collectionOfRangeElement[part] = std::make_tuple( myelements->begin(), myelements->end(), myelements );
                }
                return collectionOfRangeElement;
            }
            else
            {
                using range_base_type = std::tuple<element_reference_wrapper_const_iterator,element_reference_wrapper_const_iterator,elements_reference_wrapper_ptrtype>;
                std::map<int, std::tuple<range_base_type,range_base_type> > collectionOfRangeElement;
                for ( auto & [part,myelements] : collectionOfElements )
                {
                    myelements->shrink_to_fit();
                    auto & myghostelements = collectionOfGhostElements.at(part);
                    myghostelements->shrink_to_fit();
                    collectionOfRangeElement[part] = std::make_tuple(  std::make_tuple( myelements->begin(), myelements->end(), myelements ),
                                                                       std::make_tuple( myghostelements->begin(), myghostelements->end(), myghostelements ) );
                }
                return collectionOfRangeElement;
            }
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
        auto it = this->beginOrderedElement();
        auto en = this->endOrderedElement();
        for ( ; it!=en;++it )
        {
            auto const& elt = unwrap_ref( *it );
            if ( elt.processId() != part )
                continue;
            myelements->push_back(boost::cref(elt));
        }
        myelements->shrink_to_fit();
        return std::make_tuple( myelements->begin(), myelements->end(), myelements );
    }

    /**
     * \return the first iterator over the elements
     * on processor \p p if exist else return endElement iterator
     */
    element_const_iterator
    firstElementIteratorWithProcessId( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommElements().localRank() : p;
            auto it = this->beginElement();
            auto en = this->endElement();
            for ( ; it!=en;++it )
            {
                auto const& elt = it->second;
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
        auto it = this->beginOrderedElement();
        auto en = this->endOrderedElement();
        for ( ; it!=en;++it )
        {
            auto const& elt = unwrap_ref( *it );
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
        myelements->shrink_to_fit();
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
        auto it = this->beginOrderedElement();
        auto en = this->endOrderedElement();
        for ( ; it!=en;++it )
        {
            auto const& elt = unwrap_ref( *it );
            if ( elt.processId() != part )
                continue;
            if ( !elt.isInternal() )
                continue;
            myelements->push_back(boost::cref(elt));
        }
        myelements->shrink_to_fit();
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
        auto it = this->beginOrderedElement();
        auto en = this->endOrderedElement();
        for ( ; it!=en;++it )
        {
            auto const& elt = unwrap_ref( *it );
            if ( !elt.isGhostCell() )
                continue;
            myelements->push_back(boost::cref(elt));
        }
        myelements->shrink_to_fit();
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

    //!
    //! reserve size of container
    //! @param nElt : the size reserved
    //!
    void reserveNumberOfElement( size_type nElt )
        {
            M_elements.reserve( nElt );
            M_orderedElements.reserve( nElt );
        }

    //!
    //! add a new element in the mesh
    //! @param f a new point
    //! @return the new point from the list
    //!
    std::pair<element_iterator,bool> addElement( element_type& f, bool setid = true )
    {
        if ( setid )
            f.setId( M_elements.size() );
        auto ret = M_elements.emplace( std::make_pair( f.id(), f ) );

        if ( ret.second )
        {
            auto & newElement = ret.first->second;
            if ( !M_needToOrderElements && !M_orderedElements.empty() && unwrap_ref( M_orderedElements.back() ).id() > newElement.id() )
                M_needToOrderElements = true;
            M_orderedElements.push_back( boost::ref( newElement ) );
        }
        return ret;
    }
    //!
    //! move an element into the mesh
    //! @param f a new point
    //! @return the new point from the list
    //!
    std::pair<element_iterator,bool> addElement( element_type&& f )
        {
            auto ret = M_elements.emplace( std::make_pair( f.id(), f ) );

            if ( ret.second )
            {
                auto & newElement = ret.first->second;
                if ( !M_needToOrderElements && !M_orderedElements.empty() && unwrap_ref( M_orderedElements.back() ).id() > newElement.id() )
                    M_needToOrderElements = true;
                M_orderedElements.push_back( boost::ref( newElement ) );
            }
            return ret;
        }


    element_iterator eraseElementOnly( element_iterator it )
        {
            size_type erasedId = it->first;
            auto itret = M_elements.erase( it );
            auto itOrdered = std::find_if( M_orderedElements.begin(), M_orderedElements.end(),
                                           [&erasedId]( auto & eltWrap ) { return unwrap_ref( eltWrap ).id() == erasedId; } );
            M_orderedElements.erase( itOrdered );
            return itret;
        }


    template<typename ElementVecType>
    void updateMarker( uint16_type markerType, ElementVecType const& evec )
    {
        EntityProcessType entityProcess = (evec.functionSpace()->dof()->buildDofTableMPIExtended())? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
        auto rangeElt = Feel::elements( evec.mesh(), entityProcess );
        auto it = rangeElt.template get<1>();
        auto en = rangeElt.template get<2>();
        for ( ; it != en; ++it )
        {
            size_type eltId = boost::unwrap_ref(*it).id();
            auto & eltModified = this->elementIterator( eltId )->second;
            eltModified.setMarker( markerType, evec.localToGlobal( eltId, 0, 0 ) );
        }
    }

    template<typename ElementVecType>
    void updateMarker2( ElementVecType const& evec )
    {
        this->updateMarker( 2, evec );
    }

    template<typename ElementVecType>
    void updateMarker3( ElementVecType const& evec )
    {
        this->updateMarker( 3, evec );
    }

    template<typename IteratorRange>
    void updateMarkerWithRangeElements( uint16_type markerType, IteratorRange const& range, flag_type flag )
    {
        for ( auto const& elt : range )
            this->elementIterator( boost::unwrap_ref( elt ) )->second.setMarker( markerType,flag );
    }

    template<typename IteratorRange>
    void updateMarkerWithRangeElements( IteratorRange const& range, flag_type flag )
    {
        this->updateMarkerWithRangeElements( 1,range,flag );
    }

    template<typename IteratorRange>
    void updateMarker2WithRangeElements( IteratorRange const& range, flag_type flag )
    {
        this->updateMarkerWithRangeElements( 2,range,flag );
    }

    template<typename IteratorRange>
    void updateMarker3WithRangeElements( IteratorRange const& range, flag_type flag )
    {
        this->updateMarkerWithRangeElements( 3,range,flag );
    }


    /**
     * update the elements markers by setting them from the face markers associated to the elements
     * if 2 several markers are find in elt, the element take the last find as marker
     */
    void updateMarkersFromFaces( std::initializer_list<uint16_type> const& markersType )
    {
        auto it = beginElement(), en = endElement();
        for ( ; it != en; ++it )
        {
            auto const& elt = it->second;

            auto & eltModified = this->elementIterator( elt )->second;

            for (uint16_type f=0;f<elt.numTopologicalFaces; ++f)
            {
                auto const& theface = elt.face(f);
                for ( uint16_type const& markerType : markersType )
                {
                    if ( theface.hasMarkerType( markerType ) )
                        eltModified.addMarker( markerType, theface.marker( markerType ) );
                }
            }
        }
    }
    void updateMarkersFromFaces()
    {
        this->updateMarkersFromFaces( { 2,3 } );
    }


    void setWorldCommElements( worldcomm_ptr_t const& _worldComm )
    {
        M_worldComm = _worldComm;
    }

    void updateOrderedElement()
        {
            if ( !M_needToOrderElements )
                return;
            std::sort( M_orderedElements.begin(), M_orderedElements.end(),
                       []( auto const& a, auto const& b) -> bool
                       {
                           return unwrap_ref( a ).id() < unwrap_ref( b ).id();
                       });
            M_needToOrderElements = false;
        }

    //@}

private:

    void buildOrderedElements()
        {
            M_orderedElements.clear();
            auto it = beginElement(), en = endElement();
            size_type nElement = std::distance( it, en );
            M_orderedElements.reserve( nElement );
            for ( ; it != en ; ++it )
                M_orderedElements.push_back( boost::ref( it->second ) );
            M_needToOrderElements = true;
            this->updateOrderedElements();
        }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            if ( Archive::is_loading::value )
            {
                M_elements.clear();
                M_orderedElements.clear();
                M_needToOrderElements =  false;
                size_type nElements = 0;
                ar & BOOST_SERIALIZATION_NVP( nElements );
                element_type newElt;
                for ( size_type k=0 ; k<nElements ; ++k )
                {
                    ar & boost::serialization::make_nvp( "element", newElt );
                    this->addElement( std::move( newElt ) );
                }
            }
            else
            {
                auto it = beginOrderedElement(), en = endOrderedElement();
                size_type nElements = std::distance( it, en );
                ar & BOOST_SERIALIZATION_NVP( nElements );
                for ( ; it != en ; ++it )
                    ar & boost::serialization::make_nvp( "element", unwrap_ref( *it ) );
            }
        }


private:
    worldcomm_ptr_t M_worldComm;
    elements_type M_elements;
    ordered_elements_reference_wrapper_type M_orderedElements;
    bool M_needToOrderElements;
};
/// \endcond
} // Feel
#endif /* FEELPP_ELEMENTS_HPP */
