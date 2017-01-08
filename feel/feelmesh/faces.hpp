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
#ifndef __faces_H
#define __faces_H 1


#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include <feel/feelmesh/geoelement.hpp>
#include <feel/feelmesh/filters.hpp>

namespace Feel
{
namespace multi_index = boost::multi_index;

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
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<EntityType::nDim>, mpl::int_<EntityType::nRealDim-1> >,
                              mpl::identity<typename mpl::if_<mpl::equal_to<mpl::int_<EntityType::nDim>, mpl::int_<0> >,
                                                              mpl::identity<GeoElement0D<EntityType::nRealDim, SubFaceOf<ElementType> > >,
                                                              typename mpl::if_<mpl::equal_to<mpl::int_<EntityType::nDim>, mpl::int_<1> >,
                                                                                mpl::identity<GeoElement1D<EntityType::nRealDim, EntityType,  SubFaceOf<ElementType> > >,
                                                                                mpl::identity<GeoElement2D<EntityType::nRealDim, EntityType,  SubFaceOf<ElementType> > >
                                                                                >::type>::type>,
                              mpl::identity<typename mpl::if_<mpl::equal_to<mpl::int_<EntityType::nDim>, mpl::int_<0> >,
                                                              mpl::identity<GeoElement0D<EntityType::nRealDim, SubFaceOfMany<ElementType> > >,
                                                              typename mpl::if_<mpl::equal_to<mpl::int_<EntityType::nDim>, mpl::int_<1> >,
                                                                                mpl::identity<GeoElement1D<EntityType::nRealDim, EntityType,  SubFaceOfMany<ElementType> > >,
                                                                                mpl::identity<GeoElement2D<EntityType::nRealDim, EntityType,  SubFaceOfMany<ElementType> > >
                                                                                >::type>::type> >::type::type::type face_type;

    typedef multi_index::multi_index_container<
    face_type,
    multi_index::indexed_by<
    // sort by employee::operator<
#if 1
        multi_index::ordered_unique<multi_index::identity<face_type> >
#else
        multi_index::ordered_unique<
            multi_index::composite_key<face_type,
                                       multi_index::const_mem_fun<face_type,
                                                                  rank_type,
                                                                  &face_type::processId>,
                                       multi_index::const_mem_fun<face_type,
                                                                  size_type,
                                                                  &face_type::id> > >,
#endif

#if 0
        // sort by less<int> on marker
        multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_marker>,
                                        multi_index::composite_key<
                                            face_type,
                                            multi_index::const_mem_fun<face_type,
                                                                       Marker1 const&,
                                                                       &face_type::marker>,
                                            multi_index::const_mem_fun<face_type,
                                                                       rank_type,
                                                                       &face_type::processId>
                                            > >,
        multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_marker2>,
                                        multi_index::composite_key<
                                            face_type,
                                            multi_index::const_mem_fun<face_type,
                                                                       Marker2 const&,
                                                                       &face_type::marker2>,
                                            multi_index::const_mem_fun<face_type,
                                                                       rank_type,
                                                                       &face_type::processId>
                                            > >,
        multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_marker3>,
                                        multi_index::composite_key<
                                            face_type,
                                            multi_index::const_mem_fun<face_type,
                                                                       Marker3 const&,
                                                                       &face_type::marker3>,
                                            multi_index::const_mem_fun<face_type,
                                                                       rank_type,
                                                                       &face_type::processId>
                                            > >,
        // sort by less<int> on processId
        multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_pid>,
                                        multi_index::const_mem_fun<face_type,
                                                                   rank_type,
                                                                   &face_type::processId> >,



        // sort by less<int> on boundary
        multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_interprocessdomain>,
                                        multi_index::composite_key<
                                            face_type,
                                            multi_index::const_mem_fun<face_type,
                                                                       bool,
                                                                       &face_type::isInterProcessDomain>,
                                            multi_index::const_mem_fun<face_type,
                                                                       rank_type,
                                                                       &face_type::partition1>,
                                            multi_index::const_mem_fun<face_type,
                                                                       rank_type,
                                                                       &face_type::partition2>
                                            >
                                        >,
        // sort by less<int> on boundary
        multi_index::ordered_non_unique<multi_index::tag<Feel::detail::by_location>,
                                        multi_index::composite_key<
                                            face_type,
                                            multi_index::const_mem_fun<face_type,
                                                                       bool,
                                                                       &face_type::isOnBoundary>,
                                            multi_index::const_mem_fun<face_type,
                                                                       rank_type,
                                                                       &face_type::processId>
                                            >
                                        >
#endif
        >


    > faces_type;


    typedef typename faces_type::iterator face_iterator;
    typedef typename faces_type::const_iterator face_const_iterator;

    typedef std::vector<boost::reference_wrapper<face_type const> > faces_reference_wrapper_type;
    typedef std::shared_ptr<faces_reference_wrapper_type> faces_reference_wrapper_ptrtype;
    typedef typename faces_reference_wrapper_type::iterator face_reference_wrapper_iterator;
    typedef typename faces_reference_wrapper_type::const_iterator face_reference_wrapper_const_iterator;

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

    Faces( WorldComm const& worldComm = Environment::worldComm() )
        :
        M_worldCommFaces( worldComm ),
        M_faces()
    {}

    Faces( Faces const & f )
        :
        M_worldCommFaces( f.M_worldCommFaces ),
        M_faces( f.M_faces )
    {}

    virtual ~Faces()
    {
        this->clear();
    }
    void clear()
        {
            VLOG(1) << "deleting faces...\n";
            M_faces.clear();
        }
    //@}

    /** @name Operator overloads
     */
    //@{

    Faces& operator=( Faces const& e )
    {
        if ( this != &e )
        {
            M_worldCommFaces = e.M_worldCommFaces;
            M_faces = e.M_faces;
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
        return M_worldCommFaces;
    }

    virtual bool isEmpty() const
    {
        return M_faces.empty();
    }
    bool isBoundaryFace( face_type const & e ) const
    {
        return M_faces.find( e )->isOnBoundary();
    }
    bool isBoundaryFace( size_type const & id ) const
    {
        return M_faces.find( face_type( id ) )->isOnBoundary();
    }

    /**
     * \return \c true if element with id \p i is found, \c false otherwise
     */
    bool hasFace( size_type i ) const
    {
        return M_faces.template get<0>().find( face_type( i ) ) !=
               M_faces.template get<0>().end();
    }

    face_type const& face( size_type i ) const
    {
        return *M_faces.find( face_type( i ) );
    }

    face_iterator faceIterator( size_type i ) const
    {
        return  M_faces.find( face_type( i ) );
    }
    face_iterator faceIterator( face_type const& face ) const
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
            auto it = this->beginFace();
            auto en = this->endFace();
            for ( ; it!=en;++it )
            {
                auto const& face = *it;
                if ( face.processId() != part )
                    continue;
                if ( !face.hasMarker( markerType ) )
                    continue;
                if ( face.marker( markerType ).isOff() )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
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
            auto it = this->beginFace();
            auto en = this->endFace();
            for ( ; it!=en;++it )
            {
                auto const& face = *it;
                if ( face.processId() != part )
                    continue;
                if ( !face.hasMarker( markerType ) )
                    continue;
                if ( face.marker( markerType ).isOff() )
                    continue;
                if ( markerFlags.find( face.marker( markerType ).value() ) == markerFlags.end() )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
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
            auto it = this->beginFace();
            auto en = this->endFace();
            for ( ; it!=en;++it )
            {
                auto const& face = *it;
                if ( face.processId() != part )
                    continue;
                if ( !face.isOnBoundary() )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
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
            auto it = this->beginFace();
            auto en = this->endFace();
            for ( ; it!=en;++it )
            {
                auto const& face = *it;
                if ( face.processId() != part )
                    continue;
                if ( !face.isInternal() )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
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
            auto it = this->beginFace();
            auto en = this->endFace();
            for ( ; it!=en;++it )
            {
                auto const& face = *it;
                if ( !face.isInterProcessDomain() )
                    continue;
                if ( face.partition1() != part )
                    continue;
                if ( !allNeighbor && face.partition2() != neighbor_pid )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
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
            auto it = this->beginFace();
            auto en = this->endFace();
            for ( ; it!=en;++it )
            {
                auto const& face = *it;
                if ( face.processId() != part )
                    continue;
                if ( !face.isIntraProcessDomain( part ) )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
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
            auto it = this->beginFace();
            auto en = this->endFace();
            for ( ; it!=en;++it )
            {
                auto const& face = *it;
                if ( face.processId() != part )
                    continue;
                myfaces->push_back( boost::cref( face ) );
            }
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
        std::pair<face_iterator,bool> ret =  M_faces.insert( f );
        DLOG_IF(WARNING, ret.second == false )
            << "addFace failed, face not added to container : "
            << ret.first->id() << " face id:"
            << f.id();
        
        return ret;
    }
    //!
    //! @brief move a new face into the mesh
    //! @param f a new point
    //! @return the new point from the list
    //! 
    std::pair<face_iterator,bool> addFace( face_type&& f )
        {
            std::pair<face_iterator,bool> ret =  M_faces.insert( f );
            DLOG_IF(WARNING, ret.second == false )
                << "addFace failed, face not added to container : "
                << ret.first->id() << " face id:"
                << f.id();
        
            return ret;
        }
    //!
    //! @brief copy a new face into the mesh
    //! @param f a new point
    //! @param pos position hint where to move
    //! @return the new point from the list
    //!
    face_iterator addFace( face_iterator pos, face_type& f )
        {
            return M_faces.insert( pos, f );
        }
    //!
    //! @brief move a new face into the mesh
    //! @param f a new point
    //! @param pos position hint where to move
    //! @return the new point from the list
    //!
    face_iterator addFace( face_iterator pos, face_type&& f )
        {
            return M_faces.insert( pos, f );
        }

    //!
    //! @brief move a new face into the mesh
    //! @param f a new point
    //! @param pos position hint where to move
    //! @return the new point from the list
    //!
    template<typename... Args>
    face_iterator emplaceFace( face_iterator pos, Args&&... f )
        {
            return M_faces.emplace_hint( pos, f... );
        }

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

    /**
     * erase face at position \p position
     *
     * @param position \p position is a valid dereferenceable iterator of the index.
     *
     * @return An iterator pointing to the face immediately
     * following the one that was deleted, or \c end() if no such face
     * exists.
     */
    face_iterator eraseFace( face_iterator position )
    {
        return M_faces.erase( position );
    }

    /**
     * update the faces markers by setting them from the elements markers associated to the face
     */
    void updateMarkersFromElements( std::initializer_list<uint16_type> const& markersType )
    {
        auto it = beginFace(), en = endFace();
        for (  ; it != en; ++it )
            M_faces.modify( it,
                             [&markersType]( face_type& e )
        {
            for ( uint16_type const& markerType : markersType )
            {
                if ( !e.isConnectedTo0() )
                    continue;
                if( !e.isConnectedTo1() )
                {
                    if ( !e.element0().hasMarker( markerType ) )
                        continue;
                    flag_type tag_0 = e.element0().marker( markerType ).value();
                    e.setMarker( markerType, tag_0 );
                }
                else
                {
                    bool hasMarkerElt0 = e.element0().hasMarker( markerType );
                    bool hasMarkerElt1 = e.element1().hasMarker( markerType );
                    flag_type tag_0 = (hasMarkerElt0)? e.element0().marker( markerType ).value() : 0;
                    flag_type tag_1 = (hasMarkerElt1)? e.element1().marker( markerType ).value() : 0;
                    if ( hasMarkerElt0 && hasMarkerElt1 )
                        e.setMarker( markerType, std::max(tag_0,tag_1) );
                    else if ( hasMarkerElt0 && !hasMarkerElt1 )
                        e.setMarker( markerType, tag_0 );
                    else if ( !hasMarkerElt0 && hasMarkerElt1 )
                        e.setMarker( markerType, tag_1 );
                }
            }
        } );
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
     * update faces marker 2 from a vector whose size is exactely the number of
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
        auto update_marker = [&markerType,&evec,&id]( face_type& e )
            {
                auto dof_value = evec.localToGlobal( id, 0, 0 );
                e.setMarker( markerType, dof_value );
            };
        for ( ; it != en; ++it )
        {
            id = boost::unwrap_ref(*it).id();
            auto const& theface = face( evec.mesh()->subMeshToMesh( id ) );
            auto fid = evec.mesh()->subMeshToMesh( id );
            auto it = this->faceIterator( fid );
            bool r = M_faces.modify( it, update_marker );
            DLOG_IF(WARNING, r == false ) << "update marker2 failed for element id " << id << " face id " << fid;
        }
    }


    /**
     * update faces marker 2 from a vector whose size is exactely the number of
     * faces. This vector can be generated using a P0 discontinuous space
     * associated to a mesh whose elements are the faces
     */
    template<typename ElementVecType>
    void updateFacesMarker2( ElementVecType const& evec )
    {
        this->updateFacesMarker( 2, evec );
    }
    
    /**
     * update faces marker 3 from a vector whose size is exactely the number of
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
            M_faces.modify( this->faceIterator( boost::unwrap_ref( *it ).id() ), [&markerType,&flag]( face_type& e )
        {
            e.setMarker( markerType, flag );
        } );
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

    void setWorldCommFaces( WorldComm const& _worldComm )
    {
        M_worldCommFaces = _worldComm;
    }

    //@}

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & M_faces;
        }

private:
    WorldComm M_worldCommFaces;
    faces_type M_faces;
};
/// \endcond
} // Feel
#endif /* __faces_H */
