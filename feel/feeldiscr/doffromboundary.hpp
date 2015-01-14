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
 * \brief Local Dof contribution from boundary dof
 *
 * @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * @see DofTable, Dof, DofFromElement
 */
template <typename DofTableType, typename FEType>
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

    typedef DofTableType doftable_type;
    typedef typename doftable_type::mesh_type mesh_type;
    typedef typename doftable_type::element_type element_type;
    typedef typename doftable_type::face_type face_type;
    typedef typename doftable_type::ref_shift_type ref_shift_type;
    typedef typename doftable_type::localdof_type localdof_type;
    typedef FEType fe_type;

    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;

    static const uint16_type nOrder = fe_type::nOrder;
    static const uint16_type nDim = mesh_type::nDim;
    static const uint16_type nRealDim = mesh_type::nRealDim;
    static const uint16_type Shape = mesh_type::Shape;
    static const uint16_type nComponents = fe_type::nComponents;
    static const uint16_type nComponents1 = fe_type::nComponents1;
    static const uint16_type nComponents2 = fe_type::nComponents2;


    static const bool is_continuous = fe_type::isContinuous;
    static const bool is_discontinuous_locally = fe_type::continuity_type::is_discontinuous_locally;
    static const bool is_discontinuous_totally = fe_type::continuity_type::is_discontinuous_totally;

    static const bool is_scalar = fe_type::is_scalar;
    static const bool is_vectorial = fe_type::is_vectorial;
    static const bool is_tensor2 = fe_type::is_tensor2;
    static const bool is_modal = fe_type::is_modal;
    static const bool is_product = fe_type::is_product;

    static const bool is_p0_continuous = ( ( nOrder == 0 ) && is_continuous );

    static const uint16_type nDofPerElement = mpl::if_<mpl::bool_<is_product>, mpl::int_<fe_type::nLocalDof*nComponents1>, mpl::int_<fe_type::nLocalDof> >::type::value;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    DofFromBoundary( doftable_type* doftable, fe_type const& fe )
        :
        M_doftable( doftable ),
        M_fe( fe )
        {}

    //! destructor
    ~DofFromBoundary() {}

    //@}

    /** @name Operator overloads
     */
    //@{

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

    template<typename FaceIterator>
    void add( FaceIterator it )
        {
            bool useConnection0 = it->processId() == it->proc_first();
            if ( it->isGhostCell() )
            {
                if ( M_doftable->isElementDone( it->ad_first() ) )
                {
                    useConnection0 = true;
                }
                else
                {
                    CHECK( it->isConnectedTo1() ) << "no connection1";
                    CHECK( M_doftable->isElementDone( it->ad_second() ) ) << " no dof table define on this elt " << it->ad_second() << "\n";
                    useConnection0 = false;
                }
            }

            uint16_type lcVertex = 0;
            uint16_type lcEdge = 0;
            uint16_type lcFace = 0;

            addVertexBoundaryDof( it, useConnection0, lcVertex );
            addEdgeBoundaryDof( it, useConnection0, lcEdge );
            addFaceBoundaryDof( it, useConnection0, lcFace );
        }

    //@}

protected:

private:
    doftable_type * M_doftable;
    fe_type const& M_fe;

private:

    //! default constructor
    DofFromBoundary();
    //! copy constructor
    DofFromBoundary( DofFromBoundary const & );
    //! copy operator
    DofFromBoundary& operator=( DofFromBoundary const & o)
        {
            if (this != &o )
            {
            }
            return *this;
        }

    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator face_it, bool useConnection0, uint16_type& lc )
    {
        addVertexBoundaryDof( face_it, useConnection0, lc, mpl::bool_<(fe_type::nDofPerVertex>0)>(), mpl::int_<nDim>() );
    }
    template<typename FaceIterator> void addVertexBoundaryDof( FaceIterator /*face_it*/, bool /*useConnection0*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<1> ) {}
    template<typename FaceIterator> void addVertexBoundaryDof( FaceIterator /*face_it*/, bool /*useConnection0*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<2> ) {}
    template<typename FaceIterator> void addVertexBoundaryDof( FaceIterator /*face_it*/, bool /*useConnection0*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<3> ) {}
    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator face_it, bool useConnection0, uint16_type& lc, mpl::bool_<true>, mpl::int_<1>  )
    {
        BOOST_STATIC_ASSERT( face_type::numVertices );

        uint16_type iFaEl;
        size_type iElAd;

        if ( useConnection0 )
        {
            iElAd = face_it->ad_first();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )( face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = face_it->pos_first();
            FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
        }

        else
        {
            iElAd = face_it->ad_second();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )( face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = face_it->pos_second();
            FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
        }

        // Loop number of Dof per vertex
        const int ncdof = is_product?nComponents:1;

        for ( int c = 0; c < ncdof; ++c )
        {
            for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l, ++lc )
            {
                uint16_type ldinelt = iFaEl * fe_type::nDofPerVertex + l;
                auto const& temp= M_doftable->localToGlobal( iElAd, ldinelt, c );
                M_doftable->M_face_l2g[ face_it->id()][ lc ] = FaceDof( temp, lc, ldinelt );
            }
        }
    }
    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator face_it, bool useConnection0, uint16_type& lc, mpl::bool_<true>, mpl::int_<2>  )
    {
        BOOST_STATIC_ASSERT( face_type::numVertices );

        uint16_type iFaEl;
        size_type iElAd;

        if ( useConnection0 )
        {
            iElAd = face_it->ad_first();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )( face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = face_it->pos_first();
            FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
        }
        else
        {
            iElAd = face_it->ad_second();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )( face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = face_it->pos_second();
            FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
        }
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
                for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l, ++lcc )
                {
                    uint16_type ldinelt = iVeEl * fe_type::nDofPerVertex + l;
                    auto const& temp = M_doftable->localToGlobal( iElAd, ldinelt, c );
                    M_doftable->M_face_l2g[ face_it->id()][ lcc ] = FaceDof( temp, lcc, ldinelt );
                }
            }
        }
    }
    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator face_it, bool useConnection0, uint16_type& lc, mpl::bool_<true>, mpl::int_<3>  )
    {
        addVertexBoundaryDof( face_it, useConnection0, lc, mpl::bool_<true>(), mpl::int_<2>() );
    }
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator face_it, bool useConnection0, uint16_type& lc )
    {
        static const bool cond = fe_type::nDofPerEdge*face_type::numEdges > 0;
        addEdgeBoundaryDof( face_it, useConnection0, lc, mpl::bool_<cond>(), mpl::int_<nDim>() );
    }
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator /*face_it*/, bool /*useConnection0*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<1> ) {}
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator /*face_it*/, bool /*useConnection0*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<2> ) {}
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator /*face_it*/, bool /*useConnection0*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<3> ) {}
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator face_it, bool useConnection0, uint16_type& lc, mpl::bool_<true>, mpl::int_<2> )
    {
        uint16_type iFaEl;
        size_type iElAd;



        if ( useConnection0 )
        {
            iElAd = face_it->ad_first();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )
            ( face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = face_it->pos_first();
        }

        else
        {
            iElAd = face_it->ad_second();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )
            ( face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = face_it->pos_second();
        }

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
            for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l, ++lcc )
            {
                uint16_type ldinelt = element_type::numVertices*fe_type::nDofPerVertex +
                    iFaEl * fe_type::nDofPerEdge + l ;
                auto const& temp = M_doftable->localToGlobal( iElAd,ldinelt, c );
                M_doftable->M_face_l2g[ face_it->id()][ lcc ] = FaceDof( temp, lcc, ldinelt );
            }
        }
    }
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator face_it, bool useConnection0, uint16_type& lc, mpl::bool_<true>, mpl::int_<3> )
    {
        //BOOST_STATIC_ASSERT( face_type::numEdges );
        uint16_type iFaEl;
        size_type iElAd;

        if ( useConnection0 )
        {
            iElAd = face_it->ad_first();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )
            ( face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = face_it->pos_first();
        }

        else
        {
            iElAd = face_it->ad_second();
            FEELPP_ASSERT( iElAd != invalid_size_type_value )
            ( face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

            // local id of the face in its adjacent element
            iFaEl = face_it->pos_second();
        }

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
                for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l, ++lcc )
                {
                    uint16_type ldinelt = element_type::numVertices*fe_type::nDofPerVertex +
                        iEdEl * fe_type::nDofPerEdge + l;
                    auto const& temp = M_doftable->localToGlobal( iElAd, ldinelt, c );
                    M_doftable->M_face_l2g[ face_it->id()][ lcc ] = FaceDof( temp, lcc, ldinelt );
                }
            }
        }
    }


    template<typename FaceIterator>
    void addFaceBoundaryDof( FaceIterator face_it, bool useConnection0, uint16_type& lc )
    {
        addFaceBoundaryDof( face_it, useConnection0, lc, mpl::bool_<(face_type::numFaces*fe_type::nDofPerFace > 0)>() );
    }
    template<typename FaceIterator>
    void addFaceBoundaryDof( FaceIterator /*face_it*/, bool /*useConnection0*/, uint16_type& /*lc*/, mpl::bool_<false> )
    {
    }
    template<typename FaceIterator>
    void addFaceBoundaryDof( FaceIterator face_it, bool useConnection0, uint16_type& lc, mpl::bool_<true> )
        {
            uint16_type iFaEl;
            size_type iElAd;

            if ( useConnection0 )
            {
                iElAd = face_it->ad_first();
                FEELPP_ASSERT( iElAd != invalid_size_type_value )
                    ( face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

                // local id of the face in its adjacent element
                iFaEl = face_it->pos_first();
            }

            else
            {
                iElAd = face_it->ad_second();
                FEELPP_ASSERT( iElAd != invalid_size_type_value )
                    ( face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

                // local id of the face in its adjacent element
                iFaEl = face_it->pos_second();
            }

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
                for ( uint16_type l = 0; l < fe_type::nDofPerFace; ++l, ++lcc )
                {
                    uint16_type ldinelt = element_type::numVertices*fe_type::nDofPerVertex +
                        element_type::numEdges*fe_type::nDofPerEdge +
                        iFaEl * fe_type::nDofPerFace + l;
                    auto const& temp = M_doftable->localToGlobal( iElAd, ldinelt, c );
                    M_doftable->M_face_l2g[ face_it->id()][ lcc ] = FaceDof( temp, lcc, ldinelt );
                }
            }
        }

};
}
#endif /* FEELPP_DofFromBoundary_H */
