/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
        multi_index::ordered_unique<multi_index::identity<face_type> >,
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
                                        > >
    > faces_type;


    typedef typename faces_type::iterator face_iterator;
    typedef typename faces_type::const_iterator face_const_iterator;
    typedef typename faces_type::template index<Feel::detail::by_marker>::type marker_faces;
    typedef typename faces_type::template index<Feel::detail::by_marker2>::type marker2_faces;
    typedef typename faces_type::template index<Feel::detail::by_marker3>::type marker3_faces;
    typedef typename marker_faces::iterator marker_face_iterator;
    typedef typename marker_faces::const_iterator marker_face_const_iterator;
    typedef typename marker2_faces::iterator marker2_face_iterator;
    typedef typename marker2_faces::const_iterator marker2_face_const_iterator;
    typedef typename marker3_faces::iterator marker3_face_iterator;
    typedef typename marker3_faces::const_iterator marker3_face_const_iterator;

    typedef typename faces_type::template index<Feel::detail::by_location>::type location_faces;
    typedef typename location_faces::iterator location_face_iterator;
    typedef typename location_faces::const_iterator location_face_const_iterator;

    typedef typename faces_type::template index<Feel::detail::by_pid>::type pid_faces;
    typedef typename pid_faces::iterator pid_face_iterator;
    typedef typename pid_faces::const_iterator pid_face_const_iterator;

    typedef typename faces_type::template index<Feel::detail::by_interprocessdomain>::type interprocess_faces;
    typedef typename interprocess_faces::iterator interprocess_face_iterator;
    typedef typename interprocess_faces::const_iterator interprocess_face_const_iterator;

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


    marker_face_iterator beginFaceWithMarker()
    {
        return M_faces.template get<Feel::detail::by_marker>().begin();
    }
    marker_face_const_iterator beginFaceWithMarker() const
    {
        return M_faces.template get<Feel::detail::by_marker>().begin();
    }
    marker_face_iterator endFaceWithMarker()
    {
        return M_faces.template get<Feel::detail::by_marker>().end();
    }
    marker_face_const_iterator endFaceWithMarker() const
    {
        return M_faces.template get<Feel::detail::by_marker>().end();
    }

    marker_face_iterator beginFaceWithMarker( size_type m )
    {
        return M_faces.template get<Feel::detail::by_marker>().lower_bound( Marker1( m ) );
    }
    marker_face_const_iterator beginFaceWithMarker( size_type m ) const
    {
        return M_faces.template get<Feel::detail::by_marker>().lower_bound( Marker1( m ) );
    }
    marker_face_iterator endFaceWithMarker( size_type m )
    {
        return M_faces.template get<Feel::detail::by_marker>().upper_bound( Marker1( m ) );
    }
    marker_face_const_iterator endFaceWithMarker( size_type m ) const
    {
        return M_faces.template get<Feel::detail::by_marker>().upper_bound( Marker1( m ) );
    }

    marker2_face_iterator beginFaceWithMarker2()
    {
        return M_faces.template get<Feel::detail::by_marker2>().begin();
    }
    marker2_face_const_iterator beginFaceWithMarker2() const
    {
        return M_faces.template get<Feel::detail::by_marker2>().begin();
    }
    marker2_face_iterator endFaceWithMarker2()
    {
        return M_faces.template get<Feel::detail::by_marker2>().end();
    }
    marker2_face_const_iterator endFaceWithMarker2() const
    {
        return M_faces.template get<Feel::detail::by_marker2>().end();
    }

    marker2_face_iterator beginFaceWithMarker2( size_type m )
    {
        return M_faces.template get<Feel::detail::by_marker2>().lower_bound( Marker2( m ) );
    }
    marker2_face_const_iterator beginFaceWithMarker2( size_type m ) const
    {
        return M_faces.template get<Feel::detail::by_marker2>().lower_bound( Marker2( m ) );
    }
    marker2_face_iterator endFaceWithMarker2( size_type m )
    {
        return M_faces.template get<Feel::detail::by_marker2>().upper_bound( Marker2( m ) );
    }
    marker2_face_const_iterator endFaceWithMarker2( size_type m ) const
    {
        return M_faces.template get<Feel::detail::by_marker2>().upper_bound( Marker2( m ) );
    }

    marker3_face_iterator beginFaceWithMarker3()
    {
        return M_faces.template get<Feel::detail::by_marker3>().begin();
    }
    marker3_face_const_iterator beginFaceWithMarker3() const
    {
        return M_faces.template get<Feel::detail::by_marker3>().begin();
    }
    marker3_face_iterator endFaceWithMarker3()
    {
        return M_faces.template get<Feel::detail::by_marker3>().end();
    }
    marker3_face_const_iterator endFaceWithMarker3() const
    {
        return M_faces.template get<Feel::detail::by_marker3>().end();
    }

    marker3_face_iterator beginFaceWithMarker3( size_type m )
    {
        return M_faces.template get<Feel::detail::by_marker3>().lower_bound( Marker3( m ) );
    }
    marker3_face_const_iterator beginFaceWithMarker3( size_type m ) const
    {
        return M_faces.template get<Feel::detail::by_marker3>().lower_bound( Marker3( m ) );
    }
    marker3_face_iterator endFaceWithMarker3( size_type m )
    {
        return M_faces.template get<Feel::detail::by_marker3>().upper_bound( Marker3( m ) );
    }
    marker3_face_const_iterator endFaceWithMarker3( size_type m ) const
    {
        return M_faces.template get<Feel::detail::by_marker3>().upper_bound( Marker3( m ) );
    }

    face_iterator beginFaceWithId( size_type m )
    {
        return M_faces.lower_bound( face_type( m ) );
    }
    face_const_iterator beginFaceWithId( size_type m ) const
    {
        return M_faces.lower_bound( face_type( m ) );
    }
    face_iterator endFaceWithId( size_type m )
    {
        return M_faces.upper_bound( face_type( m ) );
    }
    face_const_iterator endFaceWithId( size_type m ) const
    {
        return M_faces.upper_bound( face_type( m ) );
    }

    pid_face_iterator beginFaceWithProcessId( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_pid>().lower_bound( /*boost::make_tuple( part )*/ part );
    }
    pid_face_const_iterator beginFaceWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_pid>().lower_bound( /*boost::make_tuple( part )*/ part );
    }
    pid_face_iterator endFaceWithProcessId( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_pid>().upper_bound( /*boost::make_tuple( part )*/ part );
    }
    pid_face_const_iterator endFaceWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_pid>().upper_bound( /*boost::make_tuple( part )*/ part );
    }

    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with marker \p m on processor \p p
     */
    std::pair<marker_face_iterator, marker_face_iterator>
    facesWithMarker( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_marker>().equal_range( boost::make_tuple( Marker1( m ), part ) );
    }

    std::pair<marker2_face_iterator, marker2_face_iterator>
    facesWithMarker2( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_marker2>().equal_range( boost::make_tuple( Marker2( m ), part ) );
    }

    std::pair<marker3_face_iterator, marker3_face_iterator>
    facesWithMarker3( size_type m, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_marker3>().equal_range( boost::make_tuple( Marker3( m ), part ) );
    }


    /**
     * \return the range of iterator \c (begin,end) over the faces
     * with marker \p m on processor \p p
     */
    std::pair<location_face_iterator, location_face_iterator>
    facesOnBoundary() const
    {
        return M_faces.template get<Feel::detail::by_location>().equal_range( boost::make_tuple( ON_BOUNDARY ) );
    }

    /**
     * \return the range of iterator \c (begin,end) over the boundary
     *  faces on processor \p p
     */
    std::pair<location_face_iterator, location_face_iterator>
    facesOnBoundary( rank_type p  ) const
    {
        return M_faces.template get<Feel::detail::by_location>().equal_range( boost::make_tuple( ON_BOUNDARY, p ) );
    }

    /**
     * \return the range of iterator \c (begin,end) over the internal faces
     * on processor \p p
     */
    std::pair<location_face_iterator, location_face_iterator>
    internalFaces( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_location>().equal_range( boost::make_tuple( INTERNAL, part ) );
    }

    /**
     * \return the range of iterator \c (begin,end) over the inter-process domain faces
     * on processor \p p
     */
    std::pair<interprocess_face_iterator, interprocess_face_iterator>
    interProcessFaces( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part =  this->worldCommFaces().localRank();
        if ( p != invalid_rank_type_value )
            return M_faces.template get<Feel::detail::by_interprocessdomain>().equal_range( boost::make_tuple( true, part, p ) );
        else
            return M_faces.template get<Feel::detail::by_interprocessdomain>().equal_range( boost::make_tuple( true, part ) );
    }

#if 0
    /**
     * \return the range of iterator \c (begin,end) over the intra-process domain faces
     * on processor \p p
     */
    std::pair<interprocess_face_iterator, interprocess_face_iterator>
    intraProcessFaces( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_interprocessdomain>().equal_range( boost::make_tuple( false, part ) );
    }
#endif

    /**
     * \return the range of iterator \c (begin,end) over the elements
     * on processor \p p
     */
    std::pair<pid_face_iterator, pid_face_iterator>
    facesWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_pid>().equal_range( part );
    }

    /**
     * get the faces container by id
     *
     *
     * @return the face container by id
     */
    typename faces_type::template nth_index<0>::type &
    facesById()
    {
        return M_faces.template get<0>();
    }

    /**
     * get the faces container by id
     *
     *
     * @return the face container by id
     */
    typename faces_type::template nth_index<0>::type const&
    facesById() const
    {
        return M_faces.template get<0>();
    }

    /**
     * get the faces container using the marker view
     *
     *
     * @return the face container using marker view
     */
    marker_faces &
    facesByMarker()
    {
        return M_faces.template get<Feel::detail::by_marker>();
    }

    /**
     * get the faces container using the marker view
     *
     *
     * @return the face container using marker view
     */
    marker_faces const&
    facesByMarker() const
    {
        return M_faces.template get<Feel::detail::by_marker>();
    }

    /**
     * get the faces container using the marker view
     *
     *
     * @return the face container using marker view
     */
    marker2_faces &
    facesByMarker2()
    {
        return M_faces.template get<Feel::detail::by_marker2>();
    }

    /**
     * get the faces container using the marker view
     *
     *
     * @return the face container using marker view
     */
    marker2_faces const&
    facesByMarker2() const
    {
        return M_faces.template get<Feel::detail::by_marker2>();
    }

    /**
     * get the faces container using the marker view
     *
     *
     * @return the face container using marker view
     */
    marker3_faces &
    facesByMarker3()
    {
        return M_faces.template get<Feel::detail::by_marker3>();
    }

    /**
     * get the faces container using the marker view
     *
     *
     * @return the face container using marker view
     */
    marker3_faces const&
    facesByMarker3() const
    {
        return M_faces.template get<Feel::detail::by_marker3>();
    }


    /**
     * get the faces container using the location view
     *
     *
     * @return the face container using location view
     */
    location_faces &
    facesByLocation()
    {
        return M_faces.template get<Feel::detail::by_location>();
    }

    /**
     * get the faces container using the location view
     *
     *
     * @return the face container using location view
     */
    location_faces const&
    facesByLocation() const
    {
        return M_faces.template get<Feel::detail::by_location>();
    }

    /**
     * get the begin() iterator on all the internal faces
     *
     * @return the begin() iterator on all the internal faces
     */
    location_face_iterator beginInternalFace( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_location>().lower_bound( boost::make_tuple( INTERNAL, part ) );
    }
    /**
     * get the end() iterator on all the internal faces
     *
     * @return the end() iterator on all the internal faces
     */
    location_face_iterator endInternalFace( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_location>().upper_bound( boost::make_tuple( INTERNAL, part ) );
    }

    /**
     * get the begin() iterator on all the internal faces
     *
     * @return the begin() iterator on all the internal faces
     */
    location_face_const_iterator beginInternalFace( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_location>().lower_bound( boost::make_tuple( INTERNAL, part ) );
    }

    /**
     * get the end() iterator on all the internal faces
     *
     * @return the end() iterator on all the internal faces
     */
    location_face_const_iterator endInternalFace( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_location>().upper_bound( boost::make_tuple( INTERNAL, part ) );
    }

    /**
     * get the begin() iterator on all the boundary faces
     *
     * @return the begin() iterator on all the boundary faces
     */
    location_face_iterator beginFaceOnBoundary( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_location>().lower_bound( boost::make_tuple( ON_BOUNDARY, part ) );
    }
    /**
     * get the end() iterator on all the boundary faces
     *
     * @return the end() iterator on all the boundary faces
     */
    location_face_iterator endFaceOnBoundary( rank_type p = invalid_rank_type_value )
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_location>().upper_bound( boost::make_tuple( ON_BOUNDARY, part ) );
    }

    /**
     * get the begin() iterator on all the boundary faces
     *
     * @return the begin() iterator on all the boundary faces
     */
    location_face_const_iterator beginFaceOnBoundary( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_location>().lower_bound( boost::make_tuple( ON_BOUNDARY, part ) );
    }

    /**
     * get the end() iterator on all the boundary faces
     *
     * @return the end() iterator on all the boundary faces
     */
    location_face_const_iterator endFaceOnBoundary( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = (p==invalid_rank_type_value)? this->worldCommFaces().localRank() : p;
        return M_faces.template get<Feel::detail::by_location>().upper_bound( boost::make_tuple( ON_BOUNDARY, part ) );
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
     * add a new face in the mesh
     * @param f a new point
     * @return the new point from the list
     */
    std::pair<face_iterator,bool> addFace( face_type& f )
    {
        std::pair<face_iterator,bool> ret =  M_faces.insert( f );
        FEELPP_ASSERT( ret.second )( ret.second )( ret.first->id() )( f.id() ).warn( "face not added to container" );
        return ret;
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
    face_iterator eraseFace( face_iterator position )
    {
        return M_faces.erase( position );
    }

    /**
     * update the faces markers by setting them from the elements markers associated to the face
     */
    void updateMarkersFromElements()
    {
        //pid_face_iterator it;
        //pid_face_iterator en
        //face_iterator it;
        //face_iterator en;
        //boost::tie( it, en ) = facesWithProcessId( this->worldCommFaces().localRank() );

        auto it = beginFace(), en = endFace();

        for (  ; it != en; ++it )
            M_faces.modify( it,
                             []( face_type& e )
        {
            int tag2_0 = e.isConnectedTo0()?e.element0().marker2().value():-1;
            int tag2_1 = e.isConnectedTo1()?e.element1().marker2().value():-1;
            int tag3_0 = e.isConnectedTo0()?e.element0().marker3().value():-1;
            int tag3_1 = e.isConnectedTo1()?e.element1().marker3().value():-1;

            if ( ( tag2_0 != -1 && tag2_0 == tag2_1 ) || e.isOnBoundary() )
                e.setMarker2( tag2_0 );
            else if ( tag2_0 != -1 && tag2_1 != -1 )
                e.setMarker2( std::max(tag2_0,tag2_1) );
            else
                e.setMarker2( 0 );

            if ( ( tag3_0 != -1 && tag3_0 == tag3_1 ) || e.isOnBoundary() )
                e.setMarker3( tag3_0 );
            else if ( tag3_0 != -1 && tag3_1 != -1 )
                e.setMarker3( std::max(tag3_0,tag3_1) );
            else
                e.setMarker3( 0 );
        } );

    }


    template<typename IteratorRange>
    void updateMarker2WithRangeFaces( IteratorRange const& range, flag_type flag )
    {
        typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_range_type;
        iterator_range_type it, en;
        boost::tie( boost::tuples::ignore, it, en ) = range;

        for (  ; it != en; ++it )
            M_faces.modify( this->faceIterator( it->id() ), [&flag]( face_type& e )
        {
            e.setMarker2( flag );
        } );
    }

    template<typename IteratorRange>
    void updateMarker3WithRangeFaces( IteratorRange const& range, flag_type flag )
    {
        typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_range_type;
        iterator_range_type it, en;
        boost::tie( boost::tuples::ignore, it, en ) = range;

        for (  ; it != en; ++it )
            M_faces.modify( this->faceIterator( it->id() ), [&flag]( face_type& e )
        {
            e.setMarker3( flag );
        } );
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
