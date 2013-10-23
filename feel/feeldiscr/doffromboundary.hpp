/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-10-23

  Copyright (C) 2013 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file dofboundary.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-10-23
 */
#ifndef FEELPP_DofFromBoundary_H
#define FEELPP_DofFromBoundary_H 1

namespace Feel
{
/**
 * \class DofFromBoundary
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class DofFromBoundary
{
public:


  /** @name Constants
*/
 //@{


//@}

/** @name Typedefs
*/
 //@{


//@}

/** @name Constructors, destructor
*/
//@{

//! default constructor
DofFromBoundary();
  //! copy constructor
DofFromBoundary( DofFromBoundary const & );
  //! destructor
~DofFromBoundary();

//@}

/** @name Operator overloads
*/
//@{

//! copy operator
DofFromBoundary& operator=( DofFromBoundary const & o)
{
 if (this != &o )
 {
 }
 return *this;
}
//@}

/** @name Accessors
*/
//@{


//@}

/** @name  Mutators
*/
//@{


//@}

/** @name  Methods
*/
//@{

    void add( face_iterator it )
        {
            uint16_type lcVertex = 0;
            uint16_type lcEdge = 0;
            uint16_type lcFace = 0;

            addVertexBoundaryDof( it, lcVertex );
            addEdgeBoundaryDof( it, lcEdge );
            addFaceBoundaryDof( it, lcFace );
        }

//@}



protected:

private:

    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator __face_it, uint16_type& lc )
    {
        addVertexBoundaryDof( __face_it, lc, mpl::bool_<(fe_type::nDofPerVertex>0)>(), mpl::int_<nDim>() );
    }
    template<typename FaceIterator> void addVertexBoundaryDof( FaceIterator /*__face_it*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<1> ) {}
    template<typename FaceIterator> void addVertexBoundaryDof( FaceIterator /*__face_it*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<2> ) {}
    template<typename FaceIterator> void addVertexBoundaryDof( FaceIterator /*__face_it*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<3> ) {}
    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator __face_it, uint16_type& lc, mpl::bool_<true>, mpl::int_<1>  )
    {
        BOOST_STATIC_ASSERT( face_type::numVertices );
#if !defined(FEELPP_ENABLE_MPI_MODE)
        // id of the element adjacent to the face
        // \warning NEED TO INVESTIGATE THIS
        size_type iElAd = __face_it->ad_first();
        FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

        // local id of the face in its adjacent element
        uint16_type iFaEl = __face_it->pos_first();
        FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
#else // MPI
        uint16_type iFaEl;
        size_type iElAd;

        if ( __face_it->processId() == __face_it->proc_first() )
        {
            iElAd = __face_it->ad_first();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = __face_it->pos_first();
            FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
        }

        else
        {
            iElAd = __face_it->ad_second();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = __face_it->pos_second();
            FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
        }

#endif
        // Loop number of Dof per vertex
        const int ncdof = is_product?nComponents:1;

        for ( int c = 0; c < ncdof; ++c )
        {
            for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l )
            {
                auto temp= this->localToGlobal( iElAd,
                                                iFaEl * fe_type::nDofPerVertex + l,
                                                c );
                M_face_l2g[ __face_it->id()][ lc++ ] = boost::make_tuple( boost::get<0>( temp ),boost::get<1>( temp ),boost::get<2>( temp ),
                                                        iFaEl * fe_type::nDofPerVertex + l );
            }
        }
    }
    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator __face_it, uint16_type& lc, mpl::bool_<true>, mpl::int_<2>  )
    {
        BOOST_STATIC_ASSERT( face_type::numVertices );
#if !defined(FEELPP_ENABLE_MPI_MODE)
        // id of the element adjacent to the face
        // \warning NEED TO INVESTIGATE THIS
        size_type iElAd = __face_it->ad_first();
        FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

        // local id of the face in its adjacent element
        uint16_type iFaEl = __face_it->pos_first();
        FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
#else // MPI
        uint16_type iFaEl;
        size_type iElAd;

        if ( __face_it->processId() == __face_it->proc_first() )
        {
            iElAd = __face_it->ad_first();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = __face_it->pos_first();
            FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
        }

        else
        {
            iElAd = __face_it->ad_second();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = __face_it->pos_second();
            FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
        }

#endif

        size_type ndofF = ( face_type::numVertices * fe_type::nDofPerVertex +
                            face_type::numEdges * fe_type::nDofPerEdge +
                            face_type::numFaces * fe_type::nDofPerFace );


        //M_dof2elt[gDof].push_back( boost::make_tuple( iElAd, lc-1, 48, 0 ) );
        // loop on face vertices
        const int ncdof = is_product?nComponents:1;

        for ( int c = 0; c < ncdof; ++c )
        {
            uint16_type lcc=c*ndofF;

            for ( uint16_type iVeFa = 0; iVeFa < face_type::numVertices; ++iVeFa )
            {
                // local vertex number (in element)
                uint16_type iVeEl = element_type::fToP( iFaEl, iVeFa );

                FEELPP_ASSERT( iVeEl != invalid_uint16_type_value ).error( "invalid local dof" );

                // Loop number of Dof per vertex
                for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l )
                {
                    auto temp = this->localToGlobal( iElAd,
                                                     iVeEl * fe_type::nDofPerVertex + l,
                                                     c );
                    M_face_l2g[ __face_it->id()][ lcc++ ] = boost::make_tuple( boost::get<0>( temp ),boost::get<1>( temp ),boost::get<2>( temp ),
                            iVeEl * fe_type::nDofPerVertex + l );
                }
            }
        }
    }
    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator __face_it, uint16_type& lc, mpl::bool_<true>, mpl::int_<3>  )
    {
        addVertexBoundaryDof( __face_it, lc, mpl::bool_<true>(), mpl::int_<2>() );
    }
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator __face_it, uint16_type& lc )
    {
        static const bool cond = fe_type::nDofPerEdge*face_type::numEdges > 0;
        addEdgeBoundaryDof( __face_it, lc, mpl::bool_<cond>(), mpl::int_<nDim>() );
    }
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator /*__face_it*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<1> ) {}
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator /*__face_it*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<2> ) {}
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator /*__face_it*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<3> ) {}
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator __face_it, uint16_type& lc, mpl::bool_<true>, mpl::int_<2> )
    {
#if !defined(FEELPP_ENABLE_MPI_MODE)
        // id of the element adjacent to the face
        // \warning NEED TO INVESTIGATE THIS
        size_type iElAd = __face_it->ad_first();
        FEELPP_ASSERT( iElAd != invalid_size_type_value )
        ( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

        // local id of the face in its adjacent element
        uint16_type iFaEl = __face_it->pos_first();
#else // MPI
        uint16_type iFaEl;
        size_type iElAd;

        if ( __face_it->processId() == __face_it->proc_first() )
        {
            iElAd = __face_it->ad_first();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )
            ( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = __face_it->pos_first();
        }

        else
        {
            iElAd = __face_it->ad_second();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )
            ( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = __face_it->pos_second();
        }

#endif


        FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
#if !defined(NDEBUG)
        DVLOG(4) << " local face id : " << iFaEl << "\n";
#endif
        size_type nVerticesF = face_type::numVertices * fe_type::nDofPerVertex;
        size_type ndofF = ( face_type::numVertices * fe_type::nDofPerVertex +
                            face_type::numEdges * fe_type::nDofPerEdge +
                            face_type::numFaces * fe_type::nDofPerFace );


        const int ncdof = is_product?nComponents:1;

        for ( int c = 0; c < ncdof; ++c )
        {
            uint16_type lcc=nVerticesF+c*ndofF;

            // Loop number of Dof per edge
            for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l )
            {
                auto temp = this->localToGlobal( iElAd,
                                                 element_type::numVertices*fe_type::nDofPerVertex +
                                                 iFaEl * fe_type::nDofPerEdge + l,
                                                 c );
                M_face_l2g[ __face_it->id()][ lcc++ ] = boost::make_tuple( boost::get<0>( temp ),boost::get<1>( temp ),boost::get<2>( temp ),
                        element_type::numVertices*fe_type::nDofPerVertex +
                        iFaEl * fe_type::nDofPerEdge + l );
            }
        }
    }
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator __face_it, uint16_type& lc, mpl::bool_<true>, mpl::int_<3> )
    {
        //BOOST_STATIC_ASSERT( face_type::numEdges );
#if !defined(FEELPP_ENABLE_MPI_MODE)
        // id of the element adjacent to the face
        // \warning NEED TO INVESTIGATE THIS
        size_type iElAd = __face_it->ad_first();
        FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

        // local id of the face in its adjacent element
        uint16_type iFaEl = __face_it->pos_first();
        FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
#else // MPI
        uint16_type iFaEl;
        size_type iElAd;

        if ( __face_it->processId() == __face_it->proc_first() )
        {
            iElAd = __face_it->ad_first();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )
            ( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = __face_it->pos_first();
        }

        else
        {
            iElAd = __face_it->ad_second();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )
            ( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = __face_it->pos_second();
        }

#endif

#if !defined(NDEBUG)
        DVLOG(4) << " local face id : " << iFaEl << "\n";
#endif
        size_type nVerticesF = face_type::numVertices * fe_type::nDofPerVertex;
        size_type ndofF = ( face_type::numVertices * fe_type::nDofPerVertex +
                            face_type::numEdges * fe_type::nDofPerEdge +
                            face_type::numFaces * fe_type::nDofPerFace );

        const int ncdof = is_product?nComponents:1;

        for ( int c = 0; c < ncdof; ++c )
        {
            uint16_type lcc=nVerticesF+c*ndofF;

            // loop on face vertices
            for ( uint16_type iEdFa = 0; iEdFa < face_type::numEdges; ++iEdFa )
            {
                // local edge number (in element)
                uint16_type iEdEl = element_type::fToE( iFaEl, iEdFa );

                FEELPP_ASSERT( iEdEl != invalid_uint16_type_value ).error( "invalid local dof" );

                // Loop number of Dof per edge
                for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l )
                {
                    auto temp = this->localToGlobal( iElAd,
                                                     element_type::numVertices*fe_type::nDofPerVertex +
                                                     iEdEl * fe_type::nDofPerEdge + l,
                                                     c );
                    M_face_l2g[ __face_it->id()][ lcc++ ] = boost::make_tuple( boost::get<0>( temp ),boost::get<1>( temp ),boost::get<2>( temp ),
                            element_type::numVertices*fe_type::nDofPerVertex +
                            iEdEl * fe_type::nDofPerEdge + l );
                }
            }
        }
    }


    template<typename FaceIterator>
    void addFaceBoundaryDof( FaceIterator __face_it, uint16_type& lc )
    {
        addFaceBoundaryDof( __face_it, lc, mpl::bool_<(face_type::numFaces*fe_type::nDofPerFace > 0)>() );
    }
    template<typename FaceIterator>
    void addFaceBoundaryDof( FaceIterator /*__face_it*/, uint16_type& /*lc*/, mpl::bool_<false> )
    {
    }
    template<typename FaceIterator>
    void addFaceBoundaryDof( FaceIterator __face_it, uint16_type& lc, mpl::bool_<true> )
    {
#if !defined(FEELPP_ENABLE_MPI_MODE)
        // id of the element adjacent to the face
        // \warning NEED TO INVESTIGATE THIS
        size_type iElAd = __face_it->ad_first();
        FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

        // local id of the face in its adjacent element
        uint16_type iFaEl = __face_it->pos_first();
        FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
#else // MPI
        uint16_type iFaEl;
        size_type iElAd;

        if ( __face_it->processId() == __face_it->proc_first() )
        {
            iElAd = __face_it->ad_first();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )
            ( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = __face_it->pos_first();
        }

        else
        {
            iElAd = __face_it->ad_second();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )
            ( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = __face_it->pos_second();
        }

#endif

#if !defined(NDEBUG)
        DVLOG(4) << " local face id : " << iFaEl << "\n";
#endif
        size_type nVerticesAndEdgeF = ( face_type::numVertices * fe_type::nDofPerVertex +
                                        face_type::numEdges * fe_type::nDofPerEdge );
        size_type ndofF = ( face_type::numVertices * fe_type::nDofPerVertex +
                            face_type::numEdges * fe_type::nDofPerEdge +
                            face_type::numFaces * fe_type::nDofPerFace );

        const int ncdof = is_product?nComponents:1;

        for ( int c = 0; c < ncdof; ++c )
        {
            uint16_type lcc=nVerticesAndEdgeF+c*ndofF;

            // Loop on number of Dof per face
            for ( uint16_type l = 0; l < fe_type::nDofPerFace; ++l )
            {
                auto temp = this->localToGlobal( iElAd,
                                                 element_type::numVertices*fe_type::nDofPerVertex +
                                                 element_type::numEdges*fe_type::nDofPerEdge +
                                                 iFaEl * fe_type::nDofPerFace + l,
                                                 c );
                M_face_l2g[ __face_it->id()][ lcc++ ] = boost::make_tuple( boost::get<0>( temp ),boost::get<1>( temp ),boost::get<2>( temp ),
                        element_type::numVertices*fe_type::nDofPerVertex +
                        element_type::numEdges*fe_type::nDofPerEdge +
                        iFaEl * fe_type::nDofPerFace + l );
            }
        }
    }

};
#endif /* FEELPP_DofFromBoundary_H */
