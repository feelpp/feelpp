//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 23 Mar 2013
//! @copyright 2013-2017 Feel++ Consortium
//!
#ifndef FEELPP_DOFFROMELEMENT_HPP
#define FEELPP_DOFFROMELEMENT_HPP 1

#include <feel/feeldiscr/dof.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feelpoly/hcurlpolynomialset.hpp>
#include <feel/feelpoly/hdivpolynomialset.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelpoly/order.hpp>

#include <feel/feelmesh/marker.hpp>
#include <vector>

namespace Feel
{
/**
 * \brief local dof contribution from an element
 *
 * @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * @see DofTable, Dof
 */
template <typename DofTableType, typename FEType>
class DofFromElement
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
    typedef typename doftable_type::dof_relation dof_relation;
    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;

    using size_type = typename mesh_type::size_type;
    using mesh_marker_type = typename doftable_type::mesh_marker_type;

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
    static inline const bool is_tensor2symm = fe_type::is_tensor2 && is_symm_v<fe_type>;
    static inline const bool is_modal = fe_type::is_modal;
    static inline const bool is_product = fe_type::is_product;

    //! @brief True if polynomial order is determined at runtime
    static constexpr bool is_order_dynamic = orderIsDynamic<fe_type>;

    static inline const bool is_p0_continuous = ( ( nOrder == 0 ) && is_continuous );

    static constexpr int numberOfDofPerElement()
    {
        if constexpr ( is_product )
            return fe_type::nLocalDof * nComponents;
        else
            return fe_type::nLocalDof;
    }
    static inline const uint16_type nDofPerElement = numberOfDofPerElement();

    //@}

    /** @name Constructors, destructor
     */
    //@{

    DofFromElement( doftable_type* doftable, fe_type const& fe )
        : M_doftable( doftable ), M_fe( fe ) {}

    //! destructor
    ~DofFromElement() {}

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

    void add( element_type const& __elt,
              size_type& next_free_dof,
              rank_type processor = 0,
              size_type shift = 0 );

    //@}

  protected:
  private:
    doftable_type* M_doftable;
    fe_type const& M_fe;
    mesh_marker_type M_emptyMarker;

    struct DofInsertEntry
    {
        uint16_type localDof = 0;
        uint16_type localEntity = 0;
        uint16_type entityDim = 0;
        size_type globalDof = 0;
        int32_type sign = 1;
        bool hasLocGlobSign = false;
        int32_type locGlobSign = 1;
        mesh_marker_type marker;
    };

    //! @name Runtime DOF count accessors
    //! @{
    //! @brief Get DOFs per vertex (runtime-aware)
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

    //! @brief Get DOFs per edge (runtime-aware)
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

    //! @brief Get DOFs per face (runtime-aware)
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

    //! @brief Get DOFs per volume (runtime-aware)
    [[nodiscard]] uint16_type runtimeDofPerVolume() const noexcept
    {
        if constexpr ( is_order_dynamic )
        {
            if constexpr ( requires( fe_type const& fe ) { fe.dofPerVolume(); } )
                return M_fe.dofPerVolume();
            else
                return M_fe.runtimeDofPerVolume();
        }
        else
            return fe_type::nDofPerVolume;
    }

    //! @brief Get total local DOFs (runtime-aware)
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
    //! @}

    [[nodiscard]] size_type canonicalPointId( size_type id ) const
    {
        auto* mesh = M_doftable->mesh();
        return mesh ? mesh->canonicalPointId( id ) : id;
    }
    [[nodiscard]] size_type canonicalEdgeId( size_type id ) const
    {
        auto* mesh = M_doftable->mesh();
        return mesh ? mesh->canonicalEdgeId( id ) : id;
    }
    [[nodiscard]] size_type canonicalFaceId( size_type id ) const
    {
        auto* mesh = M_doftable->mesh();
        return mesh ? mesh->canonicalFaceId( id ) : id;
    }

    bool addDofUsingFiniteElementLayout( element_type const& __elt,
                                         rank_type processor,
                                         size_type& next_free_dof,
                                         size_type shift )
    {
        if constexpr ( !( is_hdiv_conforming_v<fe_type> || is_hcurl_conforming_v<fe_type> ) )
            return false;

        if constexpr ( fe_type::is_modal )
            return false;

        const uint16_type nLocalDof = runtimeLocalDof();
        if ( nLocalDof == 0 )
            return false;

        const uint16_type nDofPerVertex = runtimeDofPerVertex();
        const uint16_type nDofPerEdge = runtimeDofPerEdge();
        const uint16_type nDofPerFace = runtimeDofPerFace();
        const uint16_type nDofPerVolume = runtimeDofPerVolume();

        std::vector<DofInsertEntry> entries;
        entries.reserve( nLocalDof );

        for ( uint16_type parentLid = 0; parentLid < nLocalDof; ++parentLid )
        {
            auto const layout = M_fe.localDofLayout( parentLid );
            auto const& attachment = layout.attachment;

            if ( !attachment.isValid() )
                return false;

            DofInsertEntry entry;
            entry.localDof = M_doftable->localDofId( layout.parentLocalDofId, layout.component );
            entry.entityDim = static_cast<uint16_type>( attachment.entityDim );

            switch ( attachment.entityDim )
            {
            case 0:
            {
                if ( nDofPerVertex == 0 || attachment.entityId >= element_type::numVertices ||
                     attachment.ordinal >= nDofPerVertex )
                    return false;

                auto const& point = __elt.point( attachment.entityId );
                entry.localEntity = attachment.entityId;
                entry.globalDof = this->canonicalPointId( point.id() ) * nDofPerVertex + attachment.ordinal;
                entry.marker = point.hasMarker() ? point.marker() : M_emptyMarker;
                break;
            }
            case 1:
            {
                if ( nDofPerEdge == 0 || attachment.ordinal >= nDofPerEdge )
                    return false;

                if constexpr ( nDim == 1 )
                {
                    entry.localEntity = attachment.ordinal;
                    entry.globalDof = is_p0_continuous ? attachment.ordinal : __elt.id() * nDofPerEdge + attachment.ordinal;
                    entry.marker = __elt.hasMarker() ? __elt.marker() : M_emptyMarker;
                }
                else
                {
                    if ( attachment.entityId >= element_type::numEdges )
                        return false;

                    entry.localEntity = attachment.entityId;
                    entry.globalDof = this->canonicalEdgeId( __elt.edge( attachment.entityId ).id() ) * nDofPerEdge;
                    entry.marker = __elt.edge( attachment.entityId ).hasMarker() ? __elt.edge( attachment.entityId ).marker() : M_emptyMarker;

                    if ( __elt.edgePermutation( attachment.entityId ).value() == edge_permutation_type::IDENTITY )
                    {
                        entry.globalDof += attachment.ordinal;
                        if constexpr ( is_hdiv_conforming_v<fe_type> || is_hcurl_conforming_v<fe_type> )
                        {
                            entry.hasLocGlobSign = true;
                            entry.locGlobSign = 1;
                        }
                    }
                    else if ( __elt.edgePermutation( attachment.entityId ).value() == edge_permutation_type::REVERSE_PERMUTATION )
                    {
                        entry.globalDof += nDofPerEdge - 1 - attachment.ordinal;
                        if constexpr ( is_hdiv_conforming_v<fe_type> || is_hcurl_conforming_v<fe_type> )
                        {
                            entry.sign = -1;
                            entry.hasLocGlobSign = true;
                            entry.locGlobSign = -1;
                        }
                    }
                    else
                        return false;
                }
                break;
            }
            case 2:
            {
                if ( nDofPerFace == 0 || attachment.ordinal >= nDofPerFace )
                    return false;

                if constexpr ( nDim == 2 )
                {
                    entry.localEntity = attachment.ordinal;
                    entry.globalDof = is_p0_continuous ? attachment.ordinal : __elt.id() * nDofPerFace + attachment.ordinal;
                    entry.marker = __elt.hasMarker() ? __elt.marker() : M_emptyMarker;
                }
                else if constexpr ( nDim == 3 )
                {
                    if ( attachment.entityId >= element_type::numFaces )
                        return false;

                    face_permutation_type permutation = __elt.facePermutation( attachment.entityId );
                    if ( permutation == face_permutation_type( 0 ) )
                        return false;

                    entry.localEntity = attachment.entityId;
                    entry.globalDof = this->canonicalFaceId( __elt.face( attachment.entityId ).id() ) * nDofPerFace;
                    entry.marker = __elt.face( attachment.entityId ).hasMarker() ? __elt.face( attachment.entityId ).marker() : M_emptyMarker;

                    if ( nDofPerFace == 1 || permutation == face_permutation_type( face_permutation_type::IDENTITY ) )
                        entry.globalDof += attachment.ordinal;
                    else
                    {
                        if ( !M_doftable->hasValidFacePermutation( permutation, nDofPerFace ) )
                            return false;
                        auto const& perm = M_doftable->facePermutationVector( permutation, nDofPerFace );
                        entry.globalDof += perm( attachment.ordinal );
                    }

                    if constexpr ( is_hdiv_conforming_v<fe_type> || is_hcurl_conforming_v<fe_type> )
                    {
                        entry.hasLocGlobSign = true;
                        entry.locGlobSign = ( __elt.face( attachment.entityId ).ad_first() == __elt.id() ) ? 1 : -1;
                    }
                }
                else
                    return false;
                break;
            }
            case 3:
            {
                if ( nDofPerVolume == 0 || attachment.ordinal >= nDofPerVolume )
                    return false;

                entry.localEntity = attachment.ordinal;
                entry.globalDof = is_p0_continuous ? attachment.ordinal : __elt.id() * nDofPerVolume + attachment.ordinal;
                entry.marker = __elt.hasMarker() ? __elt.marker() : M_emptyMarker;
                break;
            }
            default:
                return false;
            }

            entries.push_back( std::move( entry ) );
        }

        const size_type ie = __elt.id();
        for ( auto const& entry : entries )
        {
            M_doftable->insertDof( ie, entry.localDof, entry.localEntity,
                                   std::make_tuple( entry.entityDim, entry.globalDof ),
                                   processor, next_free_dof, entry.sign, false, shift, entry.marker );
            if ( entry.hasLocGlobSign )
                M_doftable->M_locglob_signs[ie][entry.localDof] = entry.locGlobSign;
        }

        return true;
    }

  private:
    //! default constructor
    DofFromElement();
    //! copy constructor
    DofFromElement( DofFromElement const& );
    //! copy operator
    DofFromElement& operator=( DofFromElement const& o )
    {
        if ( this != &o )
        {
        }
        return *this;
    }

    void addVertexDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
                       ref_shift_type& shifts )
    {
        // Early return if no vertex DOFs
        if constexpr ( !is_order_dynamic && fe_type::nDofPerVertex == 0 )
            return;

        // Runtime check for dynamic order
        const uint16_type nDofPerVertex = runtimeDofPerVertex();
        if ( nDofPerVertex == 0 )
            return;

        auto& [local_shift, global_shift] = shifts;

        size_type ie = __elt.id();

        uint16_type lc = local_shift;

        for ( uint16_type i = 0; i < element_type::numVertices; ++i )
        {
            auto const& thepoint = __elt.point( i );
            mesh_marker_type const& pointMarker = thepoint.hasMarker() ? thepoint.marker() : M_emptyMarker;
            for ( uint16_type l = 0; l < nDofPerVertex; ++l, ++lc )
            {
                const size_type gDof = this->canonicalPointId( thepoint.id() ) * nDofPerVertex + l;
                M_doftable->insertDof( ie, lc, i, std::make_tuple( 0, gDof ),
                                       processor, next_free_dof, 1, false, global_shift, pointMarker );
            }
        }

        // update shifts
        local_shift = lc;

#if !defined( NDEBUG )
        DVLOG( 4 ) << "[Dof::updateVolumeDof(addVertexDof] vertex proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }
    void addEdgeDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts )
    {
        // Early return if no edge DOFs
        if constexpr ( !is_order_dynamic && fe_type::nDofPerEdge == 0 )
            return;

        // Runtime check for dynamic order
        const uint16_type nDofPerEdge = runtimeDofPerEdge();
        if ( nDofPerEdge == 0 )
            return;

        if constexpr ( nDim == 1 )
        {
            auto& [local_shift, global_shift] = shifts;

            size_type ie = __elt.id();
            uint16_type lc = local_shift;
            mesh_marker_type const& eltMarker = __elt.hasMarker() ? __elt.marker() : M_emptyMarker;

            for ( uint16_type l = 0; l < nDofPerEdge; ++l, ++lc )
            {
                const size_type gDof = is_p0_continuous ? l : ie * nDofPerEdge + l;
                M_doftable->insertDof( ie, lc, l, std::make_tuple( 1, gDof ), processor, next_free_dof, 1, false, global_shift, eltMarker );
            }

            // update shifts
            local_shift = lc;
#if !defined( NDEBUG )
            DVLOG( 4 ) << "[Dof::addEdgeDof(1)] element proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        }
        else if constexpr ( nDim == 2 )
        {
            auto& [local_shift, global_shift] = shifts;

            size_type ie = __elt.id();
            uint16_type lc = local_shift;

            /** The boundary dofs are constructed in the same way if the basis is modal **/

            for ( uint16_type i = 0; i < element_type::numEdges; ++i )
            {
                mesh_marker_type const& edgeMarker = __elt.edge( i ).hasMarker() ? __elt.edge( i ).marker() : M_emptyMarker;

                for ( uint16_type l = 0; l < nDofPerEdge; ++l, ++lc )
                {
                    size_type gDof = this->canonicalEdgeId( __elt.edge( i ).id() ) * nDofPerEdge;
                    int32_type sign = 1;

                    if ( __elt.edgePermutation( i ).value() == edge_permutation_type::IDENTITY )
                    {
                        gDof += l; // both nodal and modal case
                        if ( is_hdiv_conforming<fe_type>::value || is_hcurl_conforming<fe_type>::value )
                        {

                            M_doftable->M_locglob_signs[ie][lc] = 1;
                        }
                    }

                    else if ( __elt.edgePermutation( i ).value() == edge_permutation_type::REVERSE_PERMUTATION )
                    {

                        if ( fe_type::is_modal )
                        {
                            //only half of the modes (odd polynomial order) are negative.
                            sign = ( l % 2 ) ? ( -1 ) : ( 1 );
                            gDof += l;
                        }
                        else
                            gDof += nDofPerEdge - 1 - l;
                        if ( is_hdiv_conforming<fe_type>::value || is_hcurl_conforming<fe_type>::value )
                        {
                            sign = -1;
                            M_doftable->M_locglob_signs[ie][lc] = -1;
                        }
                    }

                    else
                        FEELPP_ASSERT( 0 ).error( "invalid edge permutation" );

                    M_doftable->insertDof( ie, lc, i, std::make_tuple( 1, gDof ), processor, next_free_dof, sign, false, global_shift, edgeMarker );
                }
            }

            // update shifts
            local_shift = lc;
#if !defined( NDEBUG )
            DVLOG( 4 ) << "[Dof::addEdgeDof] edge proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        }
        else if constexpr ( nDim == 3 )
        {
            auto& [local_shift, global_shift] = shifts;

            size_type ie = __elt.id();
            uint16_type lc = local_shift;

            for ( uint16_type i = 0; i < element_type::numEdges; ++i )
            {
                mesh_marker_type const& edgeMarker = __elt.edge( i ).hasMarker() ? __elt.edge( i ).marker() : M_emptyMarker;
                for ( uint16_type l = 0; l < nDofPerEdge; ++l, ++lc )
                {
                    size_type gDof = this->canonicalEdgeId( __elt.edge( i ).id() ) * nDofPerEdge;

                    int32_type sign = 1;

                    if ( __elt.edgePermutation( i ).value() == edge_permutation_type::IDENTITY )
                    {
                        gDof += l; // both nodal and modal case
                        if ( is_hcurl_conforming<fe_type>::value )
                        {
                            M_doftable->M_locglob_signs[ie][lc] = 1;
                        }
                    }

                    else if ( __elt.edgePermutation( i ).value() == edge_permutation_type::REVERSE_PERMUTATION )
                    {

                        if constexpr ( fe_type::is_modal )
                            {
                                //only half of the modes (odd polynomial order) are negative.
                                sign = ( l % 2 ) ? ( -1 ) : ( 1 );
                                gDof += l;
                            }
                        else
                            gDof += nDofPerEdge - 1 - l;
                        if constexpr ( is_hcurl_conforming_v<fe_type> )
                        {
                            sign = -1;
                            M_doftable->M_locglob_signs[ie][lc] = -1;
                        }
                    }

                    else
                        FEELPP_ASSERT( 0 ).error( "invalid edge permutation" );

                    M_doftable->insertDof( ie, lc, i, std::make_tuple( 1, gDof ), processor, next_free_dof, sign, false, global_shift, edgeMarker );
                }
            }

            // update shifts
            local_shift = lc;
#if !defined( NDEBUG )
            DVLOG( 4 ) << "[Dof::addEdgeDof] edge proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        } // 3D
    }         // addEdgeDof
    void addFaceDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts )
    {
        // Early return if no face DOFs
        if constexpr ( !is_order_dynamic && fe_type::nDofPerFace == 0 )
            return;

        // Runtime check for dynamic order
        const uint16_type nDofPerFace = runtimeDofPerFace();
        if ( nDofPerFace == 0 )
            return;

        if constexpr ( nDim == 2 )
        {
            auto& [local_shift, global_shift] = shifts;

            size_type ie = __elt.id();
            uint16_type lc = local_shift;
            mesh_marker_type const& eltMarker = __elt.hasMarker() ? __elt.marker() : M_emptyMarker;

            for ( uint16_type l = 0; l < nDofPerFace; ++l, ++lc )
            {
                const size_type gDof = is_p0_continuous ? l : ie * nDofPerFace + l;
                M_doftable->insertDof( ie, lc, l, std::make_tuple( 2, gDof ), processor, next_free_dof, 1, false, global_shift, eltMarker );
            }

            // update shifts
            local_shift = lc;
#if !defined( NDEBUG )
            DVLOG( 4 ) << "[Dof::addFaceDof(2,true)] face proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        }
        else if constexpr ( nDim == 3 )
        {
            auto& [local_shift, global_shift] = shifts;

            size_type ie = __elt.id();

            uint16_type lc = local_shift;

            for ( uint16_type i = 0; i < element_type::numFaces; ++i )
            {
                face_permutation_type permutation = __elt.facePermutation( i );
                mesh_marker_type const& faceMarker = __elt.face( i ).hasMarker() ? __elt.face( i ).marker() : M_emptyMarker;

                DCHECK( permutation != face_permutation_type( 0 ) ) << "invalid face permutation";

                // Polynomial order in each direction
                uint16_type p = 1;
                uint16_type q = 0;

                // MaxOrder = Order - 2
                int MaxOrder = int( ( 3 + std::sqrt( 1 + 8 * nDofPerFace ) ) / 2 ) - 2;

                for ( uint16_type l = 0; l < nDofPerFace; ++l, ++lc )
                {

                    // TODO: orient the dof indices such
                    // that they match properly the faces
                    // dof of the connected faces. There
                    // are a priori many permutations of
                    // the dof face indices
                    size_type gDof = this->canonicalFaceId( __elt.face( i ).id() ) * nDofPerFace;
                    int32_type sign = 1;

                    q = q + 1;

                    if ( q > MaxOrder )
                    {
                        q = 1;
                        p = p + 1;
                        MaxOrder = MaxOrder - 1;
                    }

                    if constexpr ( !fe_type::is_modal )
                    {
                        if constexpr ( is_hdiv_conforming_v<fe_type> || is_hcurl_conforming_v<fe_type> )
                        {
                            // no need of permutation if identity or only one dof on face
                            if ( nDofPerFace == 1 || permutation == face_permutation_type( face_permutation_type::IDENTITY ) )
                                gDof += l;
                            else
                            {
                                auto const& perm = M_doftable->facePermutationVector( permutation, nDofPerFace );
                                gDof += perm( l );
                            }

                            if ( __elt.face( i ).ad_first() == __elt.id() )
                                M_doftable->M_locglob_signs[ie][lc] = 1;
                            else
                                M_doftable->M_locglob_signs[ie][lc] = -1;
                        }
                        else
                        {
                            // no need of permutation if identity or only one dof on face
                            if ( nDofPerFace == 1 || permutation == face_permutation_type( face_permutation_type::IDENTITY ) )
                                gDof += l;
                            else
                            {
                                auto const& perm = M_doftable->facePermutationVector( permutation, nDofPerFace );
                                gDof += perm( l );
                            }
                        }
                    }

                    else
                    {
                        gDof += l;

                        if ( permutation == face_permutation_type( 2 ) )
                        {
                            // Reverse sign if polynomial order in
                            // eta_1 direction is odd

                            if ( p % 2 == 0 )
                                sign = -1;
                        }
                    }

                    M_doftable->insertDof( ie, lc, i, std::make_tuple( 2, gDof ), processor, next_free_dof, sign, false, global_shift, faceMarker );
                }
            }

            // update shifts
            local_shift = lc;
#if !defined( NDEBUG )
            DVLOG( 4 ) << "[Dof::addFaceDof<3>] face proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        }
    }
    void addVolumeDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
                       ref_shift_type& shifts )
    {
        // Early return if no volume DOFs (compile-time check for static order)
        if constexpr ( !is_order_dynamic && fe_type::nDofPerVolume == 0 )
            return;

        // For elements that don't have volumes (like 2D triangles), skip
        if constexpr ( element_type::numVolumes == 0 )
            return;

        // Runtime check for dynamic order
        const uint16_type nDofPerVolume = runtimeDofPerVolume();
        if ( nDofPerVolume == 0 )
            return;

        auto& [local_shift, global_shift] = shifts;

        size_type ie = __elt.id();
        uint16_type lc = local_shift;
        mesh_marker_type const& eltMarker = __elt.hasMarker() ? __elt.marker() : M_emptyMarker;

        for ( uint16_type l = 0; l < nDofPerVolume; ++l, ++lc )
        {
            const size_type gDof = is_p0_continuous ? l : ie * nDofPerVolume + l;
            M_doftable->insertDof( ie, lc, l, std::make_tuple( 3, gDof ), processor, next_free_dof, 1, false, global_shift, eltMarker );
        }

        // update shifts
        local_shift = lc;
#if !defined( NDEBUG )
        DVLOG( 4 ) << "[Dof::updateVolumeDof(<2>)] element proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }
};

template <typename DofTableType, typename FEType>
void DofFromElement<DofTableType, FEType>::add( element_type const& __elt,
                                                size_type& next_free_dof,
                                                rank_type processor,
                                                size_type shift )
{

    //size_type tndof = M_doftable->nLocalDof(  );
    size_type nldof = M_doftable->nLocalDof( true );

    DVLOG( 3 ) << "adding dof from element " << __elt.id() << "\n";
    size_type gdofcount = shift;
    DVLOG( 3 ) << "next_free_dof " << next_free_dof << "\n";
    DVLOG( 3 ) << "current dof " << M_doftable->dofIndex( next_free_dof ) << "\n";

    /*
     * Only in the continuous , we need to have the ordering [vertex,edge,face,volume]
     */
    if constexpr ( is_continuous || is_discontinuous_locally )
    {
        if ( this->addDofUsingFiniteElementLayout( __elt, processor, next_free_dof, shift ) )
            return;

        /* idem as above but for local element
           numbering except that it is
           reset to 0 after each element */
        uint16_type ldofcount = 0;

        /* pack the shifts into a tuple */
        std::tuple<uint16_type&, size_type&> shifts{ std::ref( ldofcount ), std::ref( gdofcount ) };

        /* \warning: the order of function calls is
           crucial here we order the degrees of freedom
           wrt the topological entities of the mesh
           elements from lowest dimension (vertex) to
           highest dimension (element)
        */
        addVertexDof( __elt, processor, next_free_dof, shifts );
        addEdgeDof( __elt, processor, next_free_dof, shifts );
        addFaceDof( __elt, processor, next_free_dof, shifts );
        addVolumeDof( __elt, processor, next_free_dof, shifts );
    }

    else
    {

        size_type ie = __elt.id();

        const int ncdof = is_product ? nComponents : 1;
        const uint16_type nLocalDof = runtimeLocalDof();

        for ( uint16_type l = 0; l < nldof; ++l )
        {
            if constexpr ( is_tensor2symm )
            {
                for ( int c1 = 0; c1 < nComponents1; ++c1 )
                {
                    for ( int c2 = 0; c2 < c1; ++c2, ++next_free_dof )
                    {
                        M_doftable->M_el_l2g.insert( dof_relation( localdof_type( ie, nLocalDof * ( nComponents1 * c1 + c2 ) + l ),
                                                                   Dof( ( M_doftable->dofIndex( next_free_dof ) ), 1, false ) ) );
                        M_doftable->M_el_l2g.insert( dof_relation( localdof_type( ie, nLocalDof * ( nComponents1 * c2 + c1 ) + l ),
                                                                   Dof( ( M_doftable->dofIndex( next_free_dof ) ), 1, false ) ) );
                    }
                    M_doftable->M_el_l2g.insert( dof_relation( localdof_type( ie, nLocalDof * ( nComponents1 * c1 + c1 ) + l ),
                                                               Dof( ( M_doftable->dofIndex( next_free_dof ) ), 1, false ) ) );
                    ++next_free_dof;
                }
            }
            else
            {
                for ( int c = 0; c < ncdof; ++c, ++next_free_dof )
                {
                    M_doftable->M_el_l2g.insert( dof_relation( localdof_type( ie, nLocalDof * c + l ),
                                                               Dof( ( M_doftable->dofIndex( next_free_dof ) ), 1, false ) ) );
                }
            }
        }
    }
}
} // namespace Feel
#endif /* FEELPP_DofFromElement_H */
