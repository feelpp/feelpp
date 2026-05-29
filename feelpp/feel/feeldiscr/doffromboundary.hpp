/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-10-23
 */
#ifndef FEELPP_DofFromBoundary_H
#define FEELPP_DofFromBoundary_H 1

#include <algorithm>
#include <cmath>
#include <vector>
#include <feel/feeldiscr/doflayout.hpp>
#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/hdivpolynomialset.hpp>
#include <feel/feelpoly/hcurlpolynomialset.hpp>

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
    using global_dof_from_entity_type = typename doftable_type::global_dof_from_entity_type;
    typedef FEType fe_type;

    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;
    using size_type = typename mesh_type::size_type;
    static inline const uint16_type nOrder = fe_type::nOrder;
    static inline const uint16_type nDim = mesh_type::nDim;
    static inline const uint16_type nRealDim = mesh_type::nRealDim;
    static inline const uint16_type Shape = mesh_type::Shape;
    static inline const uint16_type nComponents = fe_type::nComponents;
    static inline const uint16_type nComponents1 = fe_type::nComponents1;
    static inline const uint16_type nComponents2 = fe_type::nComponents2;


    static inline const bool is_continuous = fe_type::isContinuous;
    static inline const bool is_discontinuous_locally = fe_type::continuity_type::is_discontinuous_locally;
    static inline const bool is_discontinuous_totally = fe_type::continuity_type::is_discontinuous_totally;

    static inline const bool is_scalar = fe_type::is_scalar;
    static inline const bool is_vectorial = fe_type::is_vectorial;
    static inline const bool is_tensor2 = fe_type::is_tensor2;
    static inline const bool is_modal = fe_type::is_modal;
    static inline const bool is_product = fe_type::is_product;

    static inline const bool is_p0_continuous = ( ( nOrder == 0 ) && is_continuous );

    static inline const uint16_type nDofPerElement = mpl::if_<mpl::bool_<is_product>, mpl::int_<fe_type::nLocalDof*nComponents>, mpl::int_<fe_type::nLocalDof> >::type::value;

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

    bool add( face_type const& face )
        {
            uint8_type connectionId = invalid_v<uint8_type>;
            if ( face.isConnectedTo0() && M_doftable->isElementDone( face.ad_first() ) )
                connectionId = 0;
            else if ( face.isConnectedTo1() && M_doftable->isElementDone( face.ad_second() ) )
                connectionId = 1;
            if ( connectionId == invalid_v<uint8_type> )
                return false;

            this->add( face, connectionId );
            return true;
        }

    void add( face_type const& face, uint8_type connectionId )
        {
            bool useConnection0 = (connectionId == 0);

            if ( this->addBoundaryDofUsingFiniteElementLayout( face, useConnection0 ) )
                return;

            if ( !this->hasCompleteFiniteElementLayout() &&
                 this->addBoundaryDofUsingFiniteElementOrdering( face, useConnection0 ) )
                return;

            uint16_type lcVertex = 0;
            uint16_type lcEdge = 0;
            uint16_type lcFace = 0;

            addVertexBoundaryDof( face, useConnection0, lcVertex );
            addEdgeBoundaryDof( face, useConnection0, lcEdge );
            addFaceBoundaryDof( face, useConnection0, lcFace );
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

    static constexpr bool is_order_dynamic = orderIsDynamic<fe_type>;

    [[nodiscard]] uint16_type runtimeDofPerVertex() const noexcept
    {
        if constexpr ( is_order_dynamic )
        {
            if constexpr ( requires( fe_type const& fe ) { fe.dofPerVertex(); } )
                return M_fe.dofPerVertex();
            else
                return M_fe.runtimeDofPerVertex();
        }
        else
            return fe_type::nDofPerVertex;
    }
    [[nodiscard]] uint16_type runtimeDofPerEdge() const noexcept
    {
        if constexpr ( is_order_dynamic )
        {
            if constexpr ( requires( fe_type const& fe ) { fe.dofPerEdge(); } )
                return M_fe.dofPerEdge();
            else
                return M_fe.runtimeDofPerEdge();
        }
        else
            return fe_type::nDofPerEdge;
    }
    [[nodiscard]] uint16_type runtimeDofPerFace() const noexcept
    {
        if constexpr ( is_order_dynamic )
        {
            if constexpr ( requires( fe_type const& fe ) { fe.dofPerFace(); } )
                return M_fe.dofPerFace();
            else
                return M_fe.runtimeDofPerFace();
        }
        else
            return fe_type::nDofPerFace;
    }
    [[nodiscard]] uint16_type runtimeLocalDof() const noexcept
    {
        if constexpr ( is_order_dynamic )
        {
            if constexpr ( requires( fe_type const& fe ) { fe.localDof(); } )
                return M_fe.localDof();
            else
                return M_fe.runtimeLocalDof();
        }
        else
            return fe_type::nLocalDof;
    }
    [[nodiscard]] size_type runtimeNDofOnFace() const noexcept
    {
        return face_type::numVertices * runtimeDofPerVertex() +
               face_type::numEdges * runtimeDofPerEdge() +
               face_type::numFaces * runtimeDofPerFace();
    }

    [[nodiscard]] bool hasCompleteFiniteElementLayout() const
    {
        if constexpr ( !FiniteElementDofLayoutProvider<fe_type> )
            return false;
        else
        {
            const uint16_type nLocalDof = finiteElementLocalDofCount( M_fe );
            if ( nLocalDof == 0 )
                return false;

            for ( uint16_type localDof = 0; localDof < nLocalDof; ++localDof )
                if ( !M_fe.localDofLayout( localDof ).attachment.isValid() )
                    return false;
            return true;
        }
    }

    void getAdjacentElementAndLocalFace( face_type const& face, bool useConnection0, size_type& iElAd, uint16_type& iFaEl ) const
    {
        if ( useConnection0 )
        {
            iElAd = face.ad_first();
            FEELPP_ASSERT( iElAd != invalid_v<size_type> )( face.id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );
            iFaEl = face.pos_first();
        }
        else
        {
            iElAd = face.ad_second();
            FEELPP_ASSERT( iElAd != invalid_v<size_type> )( face.id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );
            iFaEl = face.pos_second();
        }
        FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error( "invalid element index in face" );
    }

    template<typename MatA, typename MatB>
    static bool pointsMatch( MatA const& a, uint16_type ia, MatB const& b, uint16_type ib, double tol = 1e-12 )
    {
        CHECK( a.size1() == b.size1() ) << "incompatible point dimensions " << a.size1() << " vs " << b.size1();
        for ( uint16_type d = 0; d < static_cast<uint16_type>( a.size1() ); ++d )
        {
            if ( std::abs( a( d, ia ) - b( d, ib ) ) > tol )
                return false;
        }
        return true;
    }

    [[nodiscard]] uint16_type localFaceDofIndexFromAttachment( uint16_type iFaEl,
                                                               typename fe_type::DofAttachment const& attachment,
                                                               uint16_type nDofPerVertex,
                                                               uint16_type nDofPerEdge,
                                                               uint16_type nDofPerFace ) const
    {
        const uint16_type invalid = invalid_uint16_type_value;
        if ( !attachment.isValid() )
            return invalid;

        if ( attachment.entityDim == 0 )
        {
            if ( nDofPerVertex == 0 || attachment.ordinal >= nDofPerVertex )
                return invalid;
            for ( uint16_type iVeFa = 0; iVeFa < face_type::numVertices; ++iVeFa )
            {
                const uint16_type iVeEl = element_type::fToP( iFaEl, iVeFa );
                if ( iVeEl == attachment.entityId )
                    return static_cast<uint16_type>( iVeFa * nDofPerVertex + attachment.ordinal );
            }
            return invalid;
        }

        if ( attachment.entityDim == 1 )
        {
            if ( nDofPerEdge == 0 || attachment.ordinal >= nDofPerEdge )
                return invalid;

            if constexpr ( nDim == 2 )
            {
                if ( attachment.entityId != iFaEl )
                    return invalid;
                return static_cast<uint16_type>( face_type::numVertices * nDofPerVertex + attachment.ordinal );
            }
            else if constexpr ( nDim == 3 )
            {
                for ( uint16_type iEdFa = 0; iEdFa < face_type::numEdges; ++iEdFa )
                {
                    const uint16_type iEdEl = element_type::fToE( iFaEl, iEdFa );
                    if ( iEdEl == attachment.entityId )
                    {
                        return static_cast<uint16_type>( face_type::numVertices * nDofPerVertex +
                                                         iEdFa * nDofPerEdge + attachment.ordinal );
                    }
                }
            }
            return invalid;
        }

        if ( attachment.entityDim == 2 )
        {
            if ( nDofPerFace == 0 || attachment.ordinal >= nDofPerFace )
                return invalid;

            if constexpr ( nDim == 3 )
            {
                if ( attachment.entityId != iFaEl )
                    return invalid;
                return static_cast<uint16_type>( face_type::numVertices * nDofPerVertex +
                                                 face_type::numEdges * nDofPerEdge +
                                                 attachment.ordinal );
            }
            return invalid;
        }

        return invalid;
    }

    bool addBoundaryDofUsingFiniteElementLayout( face_type const& face, bool useConnection0 )
    {
        if constexpr ( !FiniteElementDofLayoutProvider<fe_type> )
            return false;

        if constexpr ( fe_type::is_modal )
            return false;

        if ( !this->hasCompleteFiniteElementLayout() )
            return false;

        size_type iElAd;
        uint16_type iFaEl;
        getAdjacentElementAndLocalFace( face, useConnection0, iElAd, iFaEl );

        const uint16_type nLocalDof = finiteElementLocalDofCount( M_fe );
        const uint16_type nDofPerVertex = runtimeDofPerVertex();
        const uint16_type nDofPerEdge = runtimeDofPerEdge();
        const uint16_type nDofPerFace = runtimeDofPerFace();

        const size_type ndofFPerComponent = M_fe.localDofCountOnFacet( iFaEl, true );
        if ( nLocalDof == 0 || ndofFPerComponent == 0 )
            return false;

        std::vector<global_dof_from_entity_type> faceDofs;
        std::vector<bool> hasFaceDof;

        for ( uint16_type localDof = 0; localDof < nLocalDof; ++localDof )
        {
            auto const layout = M_fe.localDofLayout( localDof );
            auto const localIndexOnFace = localFaceDofIndexFromAttachment( iFaEl, layout.attachment,
                                                                           nDofPerVertex, nDofPerEdge, nDofPerFace );
            if ( localIndexOnFace == invalid_uint16_type_value )
                continue;

            const uint16_type lcc = static_cast<uint16_type>( layout.component * ndofFPerComponent + localIndexOnFace );
            if ( faceDofs.size() <= lcc )
            {
                faceDofs.resize( lcc + 1 );
                hasFaceDof.resize( lcc + 1, false );
            }

            const uint16_type ldinelt = M_doftable->localDofId( layout.parentLocalDofId, layout.component );
            auto const& temp = M_doftable->localToGlobal( iElAd, layout.parentLocalDofId, layout.component );
            faceDofs[lcc] = FaceDof( temp, lcc, ldinelt );
            hasFaceDof[lcc] = true;
        }

        if ( faceDofs.empty() ||
             std::any_of( hasFaceDof.begin(), hasFaceDof.end(), []( bool x ) { return !x; } ) )
            return false;

        for ( uint16_type lcc = 0; lcc < static_cast<uint16_type>( faceDofs.size() ); ++lcc )
            M_doftable->M_face_l2g[face.id()][lcc] = faceDofs[lcc];

        return true;
    }

    bool addBoundaryDofUsingFiniteElementOrdering( face_type const& face, bool useConnection0 )
    {
        if constexpr ( is_hdiv_conforming_v<fe_type> || is_hcurl_conforming_v<fe_type> )
            return false;

        if constexpr ( fe_type::is_modal )
            return false;

        size_type iElAd;
        uint16_type iFaEl;
        getAdjacentElementAndLocalFace( face, useConnection0, iElAd, iFaEl );

        auto const& facePts = M_fe.points( iFaEl );
        auto const& eltPts = M_fe.points();

        const uint16_type nFacePts = static_cast<uint16_type>( facePts.size2() );
        const uint16_type nEltPts = static_cast<uint16_type>( eltPts.size2() );
        if ( nFacePts == 0 || nEltPts == 0 )
            return false;

        std::vector<uint16_type> faceToElt( nFacePts, invalid_uint16_type_value );
        for ( uint16_type q = 0; q < nFacePts; ++q )
        {
            for ( uint16_type l = 0; l < nEltPts; ++l )
            {
                if ( pointsMatch( facePts, q, eltPts, l ) )
                {
                    faceToElt[q] = l;
                    break;
                }
            }
            CHECK( faceToElt[q] != invalid_uint16_type_value )
                << "failed to map face dof point " << q << " on face " << iFaEl
                << " of element " << iElAd;
        }

        size_type ndofF = runtimeNDofOnFace();
        CHECK( nFacePts == ndofF ) << "invalid face dof count mismatch " << nFacePts << " vs " << ndofF;

        const int ncdof = is_product ? nComponents : 1;
        for ( int c = 0; c < ncdof; ++c )
        {
            const uint16_type cOffset = static_cast<uint16_type>( c * ndofF );
            for ( uint16_type q = 0; q < nFacePts; ++q )
            {
                const uint16_type lcc = cOffset + q;
                const uint16_type ldinelt = faceToElt[q];
                auto const& temp = M_doftable->localToGlobal( iElAd, ldinelt, c );
                M_doftable->M_face_l2g[face.id()][lcc] = FaceDof( temp, lcc, M_doftable->localDofId( ldinelt, c ) );
            }
        }
        return true;
    }

    void addVertexBoundaryDof( face_type const& face, bool useConnection0, uint16_type& lc )
    {
        const uint16_type nDofPerVertex = runtimeDofPerVertex();
        if ( nDofPerVertex == 0 )
            return;

        size_type iElAd;
        uint16_type iFaEl;
        getAdjacentElementAndLocalFace( face, useConnection0, iElAd, iFaEl );

        const int ncdof = is_product ? nComponents : 1;
        if constexpr ( nDim == 1 )
        {
            for ( int c = 0; c < ncdof; ++c )
            {
                for ( uint16_type l = 0; l < nDofPerVertex; ++l, ++lc )
                {
                    uint16_type ldinelt = iFaEl * nDofPerVertex + l;
                    auto const& temp = M_doftable->localToGlobal( iElAd, ldinelt, c );
                    M_doftable->M_face_l2g[face.id()][lc] = FaceDof( temp, lc, M_doftable->localDofId( ldinelt, c ) );
                }
            }
            return;
        }

        size_type ndofF = runtimeNDofOnFace();
        for ( int c = 0; c < ncdof; ++c )
        {
            uint16_type lcc = c * ndofF;
            for ( uint16_type iVeFa = 0; iVeFa < face_type::numVertices; ++iVeFa )
            {
                uint16_type iVeEl = element_type::fToP( iFaEl, iVeFa );
                FEELPP_ASSERT( iVeEl != invalid_uint16_type_value ).error( "invalid local dof" );
                for ( uint16_type l = 0; l < nDofPerVertex; ++l, ++lcc )
                {
                    uint16_type ldinelt = iVeEl * nDofPerVertex + l;
                    auto const& temp = M_doftable->localToGlobal( iElAd, ldinelt, c );
                    M_doftable->M_face_l2g[face.id()][lcc] = FaceDof( temp, lcc, M_doftable->localDofId( ldinelt, c ) );
                }
            }
        }
        lc = static_cast<uint16_type>( ncdof * ndofF );
    }

    void addEdgeBoundaryDof( face_type const& face, bool useConnection0, uint16_type& lc )
    {
        const uint16_type nDofPerEdge = runtimeDofPerEdge();
        if ( nDofPerEdge == 0 || face_type::numEdges == 0 )
            return;

        size_type iElAd;
        uint16_type iFaEl;
        getAdjacentElementAndLocalFace( face, useConnection0, iElAd, iFaEl );

        const uint16_type nDofPerVertex = runtimeDofPerVertex();
        size_type nVerticesF = face_type::numVertices * nDofPerVertex;
        size_type ndofF = runtimeNDofOnFace();
        const int ncdof = is_product ? nComponents : 1;

        if constexpr ( nDim == 2 )
        {
            for ( int c = 0; c < ncdof; ++c )
            {
                uint16_type lcc = nVerticesF + c * ndofF;
                for ( uint16_type l = 0; l < nDofPerEdge; ++l, ++lcc )
                {
                    uint16_type ldinelt = element_type::numVertices * nDofPerVertex + iFaEl * nDofPerEdge + l;
                    auto const& temp = M_doftable->localToGlobal( iElAd, ldinelt, c );
                    M_doftable->M_face_l2g[face.id()][lcc] = FaceDof( temp, lcc, M_doftable->localDofId( ldinelt, c ) );
                }
            }
        }
        else if constexpr ( nDim == 3 )
        {
            for ( int c = 0; c < ncdof; ++c )
            {
                uint16_type lcc = nVerticesF + c * ndofF;
                for ( uint16_type iEdFa = 0; iEdFa < face_type::numEdges; ++iEdFa )
                {
                    uint16_type iEdEl = element_type::fToE( iFaEl, iEdFa );
                    FEELPP_ASSERT( iEdEl != invalid_uint16_type_value ).error( "invalid local dof" );
                    for ( uint16_type l = 0; l < nDofPerEdge; ++l, ++lcc )
                    {
                        uint16_type ldinelt = element_type::numVertices * nDofPerVertex + iEdEl * nDofPerEdge + l;
                        auto const& temp = M_doftable->localToGlobal( iElAd, ldinelt, c );
                        M_doftable->M_face_l2g[face.id()][lcc] = FaceDof( temp, lcc, M_doftable->localDofId( ldinelt, c ) );
                    }
                }
            }
        }
        lc = static_cast<uint16_type>( nVerticesF + ncdof * face_type::numEdges * nDofPerEdge );
    }

    void addFaceBoundaryDof( face_type const& face, bool useConnection0, uint16_type& lc )
    {
        const uint16_type nDofPerFace = runtimeDofPerFace();
        if ( nDofPerFace == 0 || face_type::numFaces == 0 )
            return;

        size_type iElAd;
        uint16_type iFaEl;
        getAdjacentElementAndLocalFace( face, useConnection0, iElAd, iFaEl );

        const uint16_type nDofPerVertex = runtimeDofPerVertex();
        const uint16_type nDofPerEdge = runtimeDofPerEdge();
        size_type nVerticesAndEdgeF = face_type::numVertices * nDofPerVertex + face_type::numEdges * nDofPerEdge;
        size_type ndofF = runtimeNDofOnFace();
        const int ncdof = is_product ? nComponents : 1;

        for ( int c = 0; c < ncdof; ++c )
        {
            uint16_type lcc = nVerticesAndEdgeF + c * ndofF;
            for ( uint16_type l = 0; l < nDofPerFace; ++l, ++lcc )
            {
                uint16_type ldinelt = element_type::numVertices * nDofPerVertex +
                                      element_type::numEdges * nDofPerEdge +
                                      iFaEl * nDofPerFace + l;
                auto const& temp = M_doftable->localToGlobal( iElAd, ldinelt, c );
                M_doftable->M_face_l2g[face.id()][lcc] = FaceDof( temp, lcc, M_doftable->localDofId( ldinelt, c ) );
            }
        }
        lc = static_cast<uint16_type>( nVerticesAndEdgeF + ncdof * face_type::numFaces * nDofPerFace );
    }

};
}
#endif /* FEELPP_DofFromBoundary_H */
