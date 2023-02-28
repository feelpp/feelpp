/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-09-03

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007,2010 Université Joseph Fourier

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
   \file faces.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-09-03
 */
#ifndef FEELPP_MESH_FACES_HPP
#define FEELPP_MESH_FACES_HPP

#include <unordered_map>
#include <feel/feelcore/commobject.hpp>
#include <feel/feelmesh/geoelement.hpp>
#include <feel/feelmesh/filters.hpp>

namespace Feel
{

/// \cond detail
/**
 * \class Faces
 * \brief Faces container class
 *
 * @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * @see Elements, Edges, Points
 */
template<typename EntityType, typename ElementType>
class Faces
{
public:


    /** @name Typedefs
     */
    //@{
    
    typedef typename ElementType::value_type value_type;
    using index_type = typename ElementType::index_type;
    using size_type = typename ElementType::size_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<EntityType::nDim>, mpl::int_<EntityType::nRealDim-1> >,
                              mpl::identity<typename mpl::if_<mpl::equal_to<mpl::int_<EntityType::nDim>, mpl::int_<0> >,
                                                              mpl::identity<GeoElement0D<EntityType::nRealDim, SubFaceOf<ElementType>, value_type, index_type > >,
                                                              typename mpl::if_<mpl::equal_to<mpl::int_<EntityType::nDim>, mpl::int_<1> >,
                                                                                mpl::identity<GeoElement1D<EntityType::nRealDim, EntityType,  SubFaceOf<ElementType>, value_type, index_type > >,
                                                                                mpl::identity<GeoElement2D<EntityType::nRealDim, EntityType,  SubFaceOf<ElementType>, value_type > >
                                                                                >::type>::type>,
                              mpl::identity<typename mpl::if_<mpl::equal_to<mpl::int_<EntityType::nDim>, mpl::int_<0> >,
                                                              mpl::identity<GeoElement0D<EntityType::nRealDim, SubFaceOfMany<ElementType>, value_type, index_type > >,
                                                              typename mpl::if_<mpl::equal_to<mpl::int_<EntityType::nDim>, mpl::int_<1> >,
                                                                                mpl::identity<GeoElement1D<EntityType::nRealDim, EntityType,  SubFaceOfMany<ElementType>, value_type, index_type > >,
                                                                                mpl::identity<GeoElement2D<EntityType::nRealDim, EntityType,  SubFaceOfMany<ElementType>, value_type, index_type > >
                                                                                >::type>::type> >::type::type::type face_type;


    typedef std::unordered_map<index_type,face_type,std::hash<index_type>> faces_type;

    typedef typename faces_type::iterator face_iterator;
    typedef typename faces_type::const_iterator face_const_iterator;

    typedef std::vector<boost::reference_wrapper<face_type const> > faces_reference_wrapper_type;
    typedef std::shared_ptr<faces_reference_wrapper_type> faces_reference_wrapper_ptrtype;
    typedef typename faces_reference_wrapper_type::iterator face_reference_wrapper_iterator;
    typedef typename faces_reference_wrapper_type::const_iterator face_reference_wrapper_const_iterator;

    typedef std::vector<boost::reference_wrapper<face_type> > ordered_faces_reference_wrapper_type;
    typedef typename ordered_faces_reference_wrapper_type::iterator ordered_face_reference_wrapper_iterator;
    typedef typename ordered_faces_reference_wrapper_type::const_iterator ordered_face_reference_wrapper_const_iterator;

    //@}

    /**
     * @class FaceUpdatePoint
     * @brief update point data structure in face
     *
     * this structure is to be used by \p modify from
     * \c boost::multi_index
     *
     */
    struct FaceUpdatePoint
    {
        /**
         * @param index local index of the point
         * @param pt point to update
         */
        FaceUpdatePoint( uint16_type index, typename face_type::point_type const& pt )
            :
            M_index( index ),
            M_pt( pt )
        {}

        /**
         * update the \p M_index point info in element \p e using
         * \p M_pt
         */
        void operator()( face_type& e )
        {
            //VLOG(1) << "FaceUpdatePoint] update point index " << M_index << " with "<< M_pt.id() << "\n";
            e.setPoint( M_index, M_pt );
            //VLOG(1) << "FaceUpdatePoint] update point "<< e.point(M_index).id() << "\n";
        }
    private:
        uint16_type M_index;
        typename face_type::point_type const& M_pt;

    };

    /** @name Constructors, destructor
     */
    //@{

    explicit Faces( worldcomm_ptr_t  const& worldComm = Environment::worldCommPtr() )
        :
        M_worldComm( worldComm ),
        M_faces(),
        M_needToOrderFaces( false )
    {}

    Faces( Faces const & f )
        :
        M_worldComm( f.M_worldComm ),
        M_faces( f.M_faces ),
        M_needToOrderFaces( false )
    {
        this->buildOrderedFaces();
    }

    virtual ~Faces() {}

    void clear()
        {
            DVLOG(1) << "deleting faces...\n";
            M_orderedFaces.clear();
            M_faces.clear();
            M_needToOrderFaces = false;
        }
    //@}

    /** @name Operator overloads
     */
    //@{

    Faces& operator=( Faces const& e )
    {
        if ( this != &e )
        {
            M_worldComm = e.M_worldComm;
            M_faces = e.M_faces;
            this->buildOrderedFaces();
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
    faces_type & faces()
    {
        return M_faces;
    }

    /**
     * \return the faces container
     */
    faces_type const& faces() const
    {
        return M_faces;
    }

    WorldComm const& worldCommFaces() const
    {
        return *M_worldComm;
    }

    virtual bool isEmpty() const
    {
        return M_faces.empty();
    }
    bool isBoundaryFace( face_type const & e ) const
    {
        return e.isOnBoundary();
    }
    bool isBoundaryFace( size_type const & id ) const
    {
        auto itFindFace = M_faces.find( id );
        if ( itFindFace == M_faces.end() )
            return false;
        return itFindFace->second.isOnBoundary();
    }

    /**
     * \return \c true if element with id \p i is found, \c false otherwise
     */
    bool hasFace( size_type i ) const
    {
        return M_faces.find( i ) != M_faces.end();
    }

    face_type const& face( size_type i ) const
    {
        auto itFindFace = M_faces.find( i );
        CHECK( itFindFace != M_faces.end() ) << " face " << i << " does not exist";
        return itFindFace->second;
    }

    face_const_iterator faceIterator( size_type i ) const
    {
        return  M_faces.find( i );
    }
    face_iterator faceIterator( size_type i )
    {
        return  M_faces.find( i );
    }
    face_const_iterator faceIterator( face_type const& face ) const
    {
        return faceIterator( face.id() );
    }
    face_iterator faceIterator( face_type const& face )
    {
        return faceIterator( face.id() );
    }

    face_iterator beginFace()
    {
        return M_faces.begin();
    }
    face_const_iterator beginFace() const
    {
        return M_faces.begin();
    }
    face_iterator endFace()
    {
        return M_faces.end();
    }
    face_const_iterator endFace() const
    {
        return M_faces.end();
    }

    ordered_face_reference_wrapper_iterator beginOrderedFace()
        {
            return M_orderedFaces.begin();
        }
    ordered_face_reference_wrapper_const_iterator beginOrderedFace() const
        {
            return M_orderedFaces.begin();
        }
    ordered_face_reference_wrapper_iterator endOrderedFace()
        {
            return M_orderedFaces.end();
        }
    ordered_face_reference_wrapper_const_iterator endOrderedFace() const
        {
            return M_orderedFaces.end();
        }


    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with \c Id \p m
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithId( size_type m ) const
        {
            faces_reference_wrapper_ptrtype myfaces( new faces_reference_wrapper_type );
            if ( this->hasFace( m ) )
                myfaces->push_back( boost::cref( this->face( m ) ) );
            return std::make_tuple( myfaces->begin(), myfaces->end(), myfaces );
        }

    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with any \c Marker1 \p on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarkerByType( uint16_type markerType, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
            faces_reference_wrapper_ptrtype myfaces( new faces_reference_wrapper_type );
            auto it = this->beginOrderedFace();
            auto en = this->endOrderedFace();
            for ( ; it!=en;++it )
            {
                auto const& face = unwrap_ref( *it );
                if ( face.processId() != part )
                    continue;
                if ( !face.hasMarkerType( markerType ) )
                    continue;
                if ( face.marker( markerType ).isOff() )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
            myfaces->shrink_to_fit();
            return std::make_tuple( myfaces->begin(), myfaces->end(), myfaces );
        }
    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with \c Marker1 \p markerFlags on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
            faces_reference_wrapper_ptrtype myfaces( new faces_reference_wrapper_type );
            auto it = this->beginOrderedFace();
            auto en = this->endOrderedFace();
            for ( ; it!=en;++it )
            {
                auto const& face = unwrap_ref( *it );
                if ( face.processId() != part )
                    continue;
                if ( !face.hasMarkerType( markerType ) )
                    continue;
                if ( !face.marker( markerType ).hasOneOf( markerFlags ) )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
            myfaces->shrink_to_fit();
            return std::make_tuple( myfaces->begin(), myfaces->end(), myfaces );
        }
    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with \c Marker1 \p m on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarkerByType( uint16_type markerType, flag_type m, rank_type p = invalid_rank_type_value ) const
        {
            if ( m == invalid_flag_type_value )
                return this->facesWithMarkerByType( markerType, p );
            else
                return this->facesWithMarkerByType( markerType, std::set<flag_type>( { m } ), p );
        }

    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with \c Marker1 \p m on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarker( flag_type m = invalid_flag_type_value, rank_type p = invalid_rank_type_value ) const
        {
            return this->facesWithMarkerByType( 1, m, p );
        }

    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with \c Marker2 \p m on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarker2( flag_type m = invalid_flag_type_value, rank_type p = invalid_rank_type_value ) const
        {
            return this->facesWithMarkerByType( 2, m, p );
        }

    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with \c Marker3 \p m on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarker3( flag_type m = invalid_flag_type_value, rank_type p = invalid_rank_type_value ) const
        {
            return this->facesWithMarkerByType( 3, m, p );
        }


    /**
     * \return the range of iterator \c (begin,end) over the boundary
     *  faces on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesOnBoundary( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
            faces_reference_wrapper_ptrtype myfaces( new faces_reference_wrapper_type );
            auto it = this->beginOrderedFace();
            auto en = this->endOrderedFace();
            for ( ; it!=en;++it )
            {
                auto const& face = unwrap_ref( *it );
                if ( face.processId() != part )
                    continue;
                if ( !face.isOnBoundary() )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
            myfaces->shrink_to_fit();
            return std::make_tuple( myfaces->begin(), myfaces->end(), myfaces );
        }

    /**
     * \return the range of iterator \c (begin,end) over the internal faces
     * on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    internalFaces( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
            faces_reference_wrapper_ptrtype myfaces( new faces_reference_wrapper_type );
            auto it = this->beginOrderedFace();
            auto en = this->endOrderedFace();
            for ( ; it!=en;++it )
            {
                auto const& face = unwrap_ref( *it );
                if ( face.processId() != part )
                    continue;
                if ( !face.isInternal() )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
            myfaces->shrink_to_fit();
            return std::make_tuple( myfaces->begin(), myfaces->end(), myfaces );
        }

    /**
     * \return the range of iterator \c (begin,end) over the inter-process domain faces
     * on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    interProcessFaces( rank_type neighbor_pid = invalid_rank_type_value ) const
        {
            // const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
            bool allNeighbor = ( neighbor_pid == invalid_rank_type_value );
            const rank_type part = this->worldCommFaces().localRank();
            faces_reference_wrapper_ptrtype myfaces( new faces_reference_wrapper_type );
            auto it = this->beginOrderedFace();
            auto en = this->endOrderedFace();
            for ( ; it!=en;++it )
            {
                auto const& face = unwrap_ref( *it );
                if ( !face.isInterProcessDomain() )
                    continue;
                if ( face.partition1() != part )
                    continue;
                if ( !allNeighbor && face.partition2() != neighbor_pid )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
            myfaces->shrink_to_fit();
            return std::make_tuple( myfaces->begin(), myfaces->end(), myfaces );
        }

    /**
     * \return the range of iterator \c (begin,end) over the intra-process domain faces
     * on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    intraProcessFaces( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
            faces_reference_wrapper_ptrtype myfaces( new faces_reference_wrapper_type );
            auto it = this->beginOrderedFace();
            auto en = this->endOrderedFace();
            for ( ; it!=en;++it )
            {
                auto const& face = unwrap_ref( *it );
                if ( face.processId() != part )
                    continue;
                if ( !face.isIntraProcessDomain( part ) )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
            myfaces->shrink_to_fit();
            return std::make_tuple( myfaces->begin(), myfaces->end(), myfaces );
        }

    /**
     * \return the range of iterator \c (begin,end) over the faces
     * on processor \p p
     */
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithProcessId( rank_type p = invalid_rank_type_value ) const
        {
            const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
            faces_reference_wrapper_ptrtype myfaces( new faces_reference_wrapper_type );
            auto it = this->beginOrderedFace();
            auto en = this->endOrderedFace();
            for ( ; it!=en;++it )
            {
                auto const& face = unwrap_ref( *it );
                if ( face.processId() != part )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
            myfaces->shrink_to_fit();
            return std::make_tuple( myfaces->begin(), myfaces->end(), myfaces );
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
    //! @brief add a new face in the mesh
    //!  @param f a new point
    //! @return the new point from the list
    //! 
    std::pair<face_iterator,bool> addFace( face_type& f )
    {
        std::pair<face_iterator,bool> ret =  M_faces.emplace/*insert*/( std::make_pair( f.id(),f ) );
        DLOG_IF(WARNING, ret.second == false )
            << "addFace failed, face not added to container : "
            << ret.first->second.id() << " face id:"
            << f.id();

        if ( ret.second )
        {
            auto & newFace = ret.first->second;
            if ( !M_needToOrderFaces && !M_orderedFaces.empty() && unwrap_ref( M_orderedFaces.back() ).id() > newFace.id() )
                M_needToOrderFaces = true;
            M_orderedFaces.push_back( boost::ref( newFace ) );
        }

        return ret;
    }
    //!
    //! @brief move a new face into the mesh
    //! @param f a new point
    //! @return the new point from the list
    //!
    std::pair<face_iterator,bool> addFace( face_type&& f )
        {
            std::pair<face_iterator,bool> ret =  M_faces.emplace/*insert*/( std::make_pair( f.id(),f ) );
            DLOG_IF(WARNING, ret.second == false )
                << "addFace failed, face not added to container : "
                << ret.first->second.id() << " face id:"
                << f.id();

            if ( ret.second )
            {
                auto & newFace = ret.first->second;
                if ( !M_needToOrderFaces && !M_orderedFaces.empty() && unwrap_ref( M_orderedFaces.back() ).id() > newFace.id() )
                    M_needToOrderFaces = true;
                M_orderedFaces.push_back( boost::ref( newFace ) );
            }

            return ret;
        }
#if 0
    //!
    //! @brief move a new face into the mesh
    //! @param f a new point
    //! @param pos position hint where to move
    //! @return the new point from the list
    //!
    template<typename... Args>
    face_iterator emplaceFace( Args&&... f )
        {
            return M_faces.emplace( f... );
        }
#endif
    /**
     * erase face at position \p position
     *
     * @param position \p position is a valid dereferenceable iterator of the index.
     *
     * @return An iterator pointing to the face immediately
     * following the one that was deleted, or \c end() if no such face
     * exists.
     */
    face_iterator eraseFace( face_iterator it )
    {
        size_type erasedId = it->first;
        auto itOrdered = std::find_if( M_orderedFaces.begin(), M_orderedFaces.end(),
                                       [&erasedId]( auto & faceWrap ) { return unwrap_ref( faceWrap ).id() == erasedId; } );
        auto itret = M_faces.erase( it );
        M_orderedFaces.erase( itOrdered );
        return itret;
    }

    /**
     * update the faces markers by setting them from the elements markers associated to the face
     */
    void updateMarkersFromElements( std::initializer_list<uint16_type> const& markersType )
    {
        auto it = beginFace(), en = endFace();
        for (  ; it != en; ++it )
        {
            auto & faceModified = it->second;
            for ( uint16_type const& markerType : markersType )
            {
                if ( !faceModified.isConnectedTo0() )
                    continue;
                if( !faceModified.isConnectedTo1() )
                {
                    if ( !faceModified.element0().hasMarkerType( markerType ) )
                        continue;
                    flag_type tag_0 = faceModified.element0().marker( markerType ).value();
                    faceModified.setMarker( markerType, tag_0 );
                }
                else
                {
                    bool hasMarkerElt0 = faceModified.element0().hasMarkerType( markerType );
                    bool hasMarkerElt1 = faceModified.element1().hasMarkerType( markerType );
                    flag_type tag_0 = (hasMarkerElt0)? faceModified.element0().marker( markerType ).value() : 0;
                    flag_type tag_1 = (hasMarkerElt1)? faceModified.element1().marker( markerType ).value() : 0;
                    if ( hasMarkerElt0 && hasMarkerElt1 )
                        faceModified.setMarker( markerType, std::max(tag_0,tag_1) );
                    else if ( hasMarkerElt0 && !hasMarkerElt1 )
                        faceModified.setMarker( markerType, tag_0 );
                    else if ( !hasMarkerElt0 && hasMarkerElt1 )
                        faceModified.setMarker( markerType, tag_1 );
                }
            }
        }
    }
    /**
     * update the faces markers by setting them from the elements markers associated to the face
     */
    void updateMarkersFromElements( uint16_type markerType )
        {
            this->updateMarkersFromElements( { markerType } );
        }
    /**
     * update the faces markers by setting them from the elements markers associated to the face
     */
    void updateMarkersFromElements()
    {
        this->updateMarkersFromElements( { 2,3 } );
    }

    /**
     * update faces marker 2 from a vector whose size is exactly the number of
     * faces. This vector can be generated using a P0 discontinuous space
     * associated to a mesh whose elements are the faces
     */
    template<typename ElementVecType>
    void updateFacesMarker( uint16_type markerType, ElementVecType const& evec )
    {
        auto rangeElt = Feel::elements( evec.mesh() );
        auto it = rangeElt.template get<1>();
        auto en = rangeElt.template get<2>();
        size_type id = 0;
        for ( ; it != en; ++it )
        {
            auto const& elt = unwrap_ref( *it );
            id = elt.id();
            size_type fid = evec.mesh()->subMeshToMesh( id );
            auto & faceModified = this->faceIterator( fid )->second;
            auto dof_value = evec.localToGlobal( id, 0, 0 );
            faceModified.setMarker( markerType, dof_value );
        }
    }


    /**
     * update faces marker 2 from a vector whose size is exactly the number of
     * faces. This vector can be generated using a P0 discontinuous space
     * associated to a mesh whose elements are the faces
     */
    template<typename ElementVecType>
    void updateFacesMarker2( ElementVecType const& evec )
    {
        this->updateFacesMarker( 2, evec );
    }

    /**
     * update faces marker 3 from a vector whose size is exactly the number of
     * faces. This vector can be generated using a P0 discontinuous space
     * associated to a mesh whose elements are the faces
     */
    template<typename ElementVecType>
    void updateFacesMarker3( ElementVecType const& evec )
    {
        this->updateFacesMarker( 3, evec );
    }


    template<typename IteratorRange>
    void updateMarkerWithRangeFaces( uint16_type markerType, IteratorRange const& range, flag_type flag )
    {
        auto it = boost::get<1>( range );
        auto en = boost::get<2>( range );
        for (  ; it != en; ++it )
        {
            auto & faceModified = this->faceIterator( boost::unwrap_ref( *it ).id() )->second;
            faceModified.setMarker( markerType, flag );
        }
    }
    template<typename IteratorRange>
    void updateMarker2WithRangeFaces( IteratorRange const& range, flag_type flag )
    {
        this->updateMarkerWithRangeFaces( 2, range, flag );
    }
    template<typename IteratorRange>
    void updateMarker3WithRangeFaces( IteratorRange const& range, flag_type flag )
    {
        this->updateMarkerWithRangeFaces( 3, range, flag );
    }

    void setWorldCommFaces( worldcomm_ptr_t const& _worldComm )
    {
        M_worldComm = _worldComm;
    }

    void updateOrderedFace()
        {
            if ( !M_needToOrderFaces )
                return;
            std::sort( M_orderedFaces.begin(), M_orderedFaces.end(),
                       []( auto const& a, auto const& b) -> bool
                       {
                           return unwrap_ref( a ).id() < unwrap_ref( b ).id();
                       });
            M_needToOrderFaces = false;
        }

    //@}

private:

    void buildOrderedFaces()
        {
            M_orderedFaces.clear();
            auto it = beginFace(), en = endFace();
            size_type nFace = std::distance( it, en );
            M_orderedFaces.reserve( nFace );
            for ( ; it != en ; ++it )
                M_orderedFaces.push_back( boost::ref( it->second ) );
            M_needToOrderFaces = true;
            this->updateOrderedFaces();
        }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            if ( Archive::is_loading::value )
            {
                M_faces.clear();
                M_orderedFaces.clear();
                M_needToOrderFaces = false;
                size_type nFaces = 0;
                ar & BOOST_SERIALIZATION_NVP( nFaces );
                face_type newFace;
                for ( size_type k=0 ; k<nFaces ; ++k )
                {
                    ar & boost::serialization::make_nvp( "face", newFace );
                    this->addFace( std::move( newFace ) );
                }
            }
            else
            {
                auto it = beginOrderedFace(), en = endOrderedFace();
                size_type nFaces = std::distance( it, en );
                ar & BOOST_SERIALIZATION_NVP( nFaces );
                for ( ; it != en ; ++it )
                    ar & boost::serialization::make_nvp( "face", unwrap_ref( *it ) );
            }
        }

private:
    worldcomm_ptr_t M_worldComm;
    faces_type M_faces;
    ordered_faces_reference_wrapper_type M_orderedFaces;
    bool M_needToOrderFaces;
};
/// \endcond
} // Feel
#endif /* __faces_H */
