/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2015-03-06

  Copyright (C) 2015 Feel++ Consortium

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
#ifndef FEELPP_DOFFROMEDGE_HPP
#define FEELPP_DOFFROMEDGE_HPP 1

#include <algorithm>
#include <vector>
#include <feel/feelpoly/order.hpp>

namespace Feel
{
/**
 * \brief Local Dof contribution from edge dof
 *
 * @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * @see DofTable, Dof, DofFromElement
 */
template <typename DofTableType, typename FEType>
class DofFromEdge
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
    using mesh_type = typename doftable_type::mesh_type;
    typedef typename doftable_type::element_type element_type;
    typedef typename doftable_type::face_type face_type;
    typedef typename doftable_type::edge_type edge_type;
    typedef typename doftable_type::ref_shift_type ref_shift_type;
    typedef typename doftable_type::localdof_type localdof_type;
    typedef FEType fe_type;
    
    
    typedef typename element_type::edge_permutation_type edge_permutation_type;
    
    using global_dof_from_entity_type = typename doftable_type::global_dof_from_entity_type;
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

    static constexpr bool is_order_dynamic = orderIsDynamic<fe_type>;
    static inline const bool is_p0_continuous = ( ( nOrder == 0 ) && is_continuous );

    static inline const uint16_type nDofPerElement = mpl::if_<mpl::bool_<is_product>, mpl::int_<fe_type::nLocalDof*nComponents1>, mpl::int_<fe_type::nLocalDof> >::type::value;

    //@}

    /** @name Constructors, destructor
     */
    //@{
    DofFromEdge() = delete;
    //! copy constructor
    DofFromEdge( DofFromEdge const & ) = default;
    //! copy operator
    DofFromEdge& operator=( DofFromEdge const & o) = default;

    
    DofFromEdge( doftable_type* doftable )
        :
        M_doftable( doftable ),
        M_fe( doftable->fe() )
        {}

    //! destructor
    ~DofFromEdge() {}

    //@}

    /** @name  Methods
     */
    //@{

    std::vector<global_dof_from_entity_type> operator()( size_type elid, uint16_type edge_id ) const;

    template<typename Iterator>
    std::vector<global_dof_from_entity_type> operator()( Iterator it ) const
        {
            CHECK( it->elements().size() ) << "Invalid call, there are no elements connected to the edge, edges must be marked to operate on them";
            size_type eid = it->elements().begin()->first;
            size_type edgeid_in_element = it->elements().begin()->second;
            return this->operator()( eid, edgeid_in_element );
        }

    
    //@}

protected:

private:
    doftable_type * M_doftable;
    fe_type const& M_fe;

private:
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
    [[nodiscard]] size_type runtimeNDofOnEdge() const noexcept
    {
        return edge_type::numVertices * runtimeDofPerVertex() +
               edge_type::numEdges * runtimeDofPerEdge();
    }

    [[nodiscard]] uint16_type localEdgeDofIndexFromAttachment( uint16_type edge_id,
                                                               typename fe_type::DofAttachment const& attachment,
                                                               uint16_type nDofPerVertex,
                                                               uint16_type nDofPerEdge ) const
    {
        const uint16_type invalid = invalid_uint16_type_value;
        if ( !attachment.isValid() )
            return invalid;

        if ( attachment.entityDim == 0 )
        {
            if ( nDofPerVertex == 0 || attachment.ordinal >= nDofPerVertex )
                return invalid;

            for ( uint16_type iVeEd = 0; iVeEd < edge_type::numVertices; ++iVeEd )
            {
                const uint16_type iVeEl = element_type::eToP( edge_id, iVeEd );
                if ( iVeEl == attachment.entityId )
                    return static_cast<uint16_type>( iVeEd * nDofPerVertex + attachment.ordinal );
            }
            return invalid;
        }

        if ( attachment.entityDim == 1 )
        {
            if ( nDofPerEdge == 0 || attachment.ordinal >= nDofPerEdge )
                return invalid;

            if ( attachment.entityId == edge_id )
                return static_cast<uint16_type>( edge_type::numVertices * nDofPerVertex + attachment.ordinal );
        }

        return invalid;
    }

    bool addEdgeDofUsingFiniteElementLayout( size_type el_id, uint16_type edge_id,
                                             std::vector<global_dof_from_entity_type>& edge_dof ) const
    {
        if constexpr ( !( is_hdiv_conforming_v<fe_type> || is_hcurl_conforming_v<fe_type> ) )
            return false;

        if constexpr ( fe_type::is_modal )
            return false;

        const uint16_type nLocalDof = runtimeLocalDof();
        const uint16_type nDofPerVertex = runtimeDofPerVertex();
        const uint16_type nDofPerEdge = runtimeDofPerEdge();
        const size_type ndofF = runtimeNDofOnEdge();
        if ( nLocalDof == 0 || ndofF == 0 )
            return false;

        const int ncdof = is_product ? nComponents : 1;
        std::vector<global_dof_from_entity_type> edgeDofs( ncdof * ndofF );
        std::vector<bool> hasEdgeDof( ncdof * ndofF, false );

        for ( uint16_type parentLid = 0; parentLid < nLocalDof; ++parentLid )
        {
            auto const layout = M_fe.localDofLayout( parentLid );
            const uint16_type localIndexOnEdge =
                localEdgeDofIndexFromAttachment( edge_id, layout.attachment, nDofPerVertex, nDofPerEdge );
            if ( localIndexOnEdge == invalid_uint16_type_value )
                continue;

            for ( int c = 0; c < ncdof; ++c )
            {
                const uint16_type lcc = static_cast<uint16_type>( c * ndofF + localIndexOnEdge );
                const uint16_type ldinelt = M_doftable->localDofId( parentLid, c );
                auto const& temp = M_doftable->localToGlobal( el_id, parentLid, c );
                edgeDofs[lcc] = global_dof_from_entity_type( temp, lcc, ldinelt );
                hasEdgeDof[lcc] = true;
            }
        }

        if ( std::any_of( hasEdgeDof.begin(), hasEdgeDof.end(), []( bool x ) { return !x; } ) )
            return false;

        edge_dof = std::move( edgeDofs );
        return true;
    }

    template<typename Container>
    void addVertexEdgeDof( size_type el_id, uint16_type edge_id, uint16_type& lc, std::back_insert_iterator<Container> d ) const
    {
        if constexpr (fe_type::nDofPerVertex>0)
        {
            if constexpr ( nDim == 1 )
            {
                BOOST_STATIC_ASSERT( edge_type::numVertices );

                // Loop number of Dof per vertex
                const int ncdof = is_product?nComponents:1;

                for ( int c = 0; c < ncdof; ++c )
                {
                    for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l, ++lc, ++d )
                    {
                        uint16_type ldinelt = edge_id * fe_type::nDofPerVertex + l;
                        auto const& temp= M_doftable->localToGlobal( el_id, ldinelt, c );
                        *d = std::move( global_dof_from_entity_type( temp, lc, M_doftable->localDofId( ldinelt, c ) ) );
                    }
                }
            }
            if constexpr ( nDim > 1 )
            {
                BOOST_STATIC_ASSERT( edge_type::numVertices );

                size_type ndofF = ( edge_type::numVertices * fe_type::nDofPerVertex +
                                    edge_type::numEdges * fe_type::nDofPerEdge );
                
                // loop on edge vertices
                const int ncdof = is_product?nComponents:1;

                for ( int c = 0; c < ncdof; ++c )
                {
                    for ( uint16_type iVeFa = 0; iVeFa < edge_type::numVertices; ++iVeFa )
                    {
                        // local vertex number (in element)
                        uint16_type iVeEl = element_type::eToP( edge_id, iVeFa );

                        DCHECK( iVeEl != invalid_uint16_type_value ) <<  "invalid local dof";

                        // Loop number of Dof per vertex
                        for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l, ++lc, ++d )
                        {
                            uint16_type ldinelt = iVeEl * fe_type::nDofPerVertex + l;
                            auto const& temp = M_doftable->localToGlobal( el_id, ldinelt, c );
                            *d = std::move( global_dof_from_entity_type( temp, lc, M_doftable->localDofId( ldinelt, c ) ) );
                        }
                    }
                }
            }
        }
    }
    template<typename Container>
    void addEdgeEdgeDof( size_type el_id, uint16_type edge_id, uint16_type& lc, std::back_insert_iterator<Container> d ) const
    {
        constexpr bool cond = fe_type::nDofPerEdge*face_type::numEdges > 0;
        if constexpr ( cond )
        {
            size_type nVerticesF = edge_type::numVertices * fe_type::nDofPerVertex;
            size_type ndofF = ( edge_type::numVertices * fe_type::nDofPerVertex +
                                edge_type::numEdges * fe_type::nDofPerEdge );

            const int ncdof = is_product?nComponents:1;

            for ( int c = 0; c < ncdof; ++c )
            {
                uint16_type lcc=nVerticesF+c*ndofF;

                // Loop number of Dof per edge
                for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l, ++lcc, ++lc, ++d )
                {
                    uint16_type ldinelt = element_type::numVertices*fe_type::nDofPerVertex +
                        edge_id * fe_type::nDofPerEdge + l ;
                    auto const& temp = M_doftable->localToGlobal( el_id,ldinelt, c );
                    *d = std::move( global_dof_from_entity_type( temp, lcc, M_doftable->localDofId( ldinelt, c ) ) );
                }
            }
        }
    }
}; // class

template <typename DofTableType, typename FEType>
std::vector<typename DofFromEdge<DofTableType,FEType>::global_dof_from_entity_type> 
DofFromEdge<DofTableType,FEType>::operator()( size_type elid, uint16_type edge_id ) const
{
    std::vector<global_dof_from_entity_type> edge_dof;
    if ( this->addEdgeDofUsingFiniteElementLayout( elid, edge_id, edge_dof ) )
        return edge_dof;

    size_type ndofF = ( edge_type::numVertices * runtimeDofPerVertex() +
                        edge_type::numEdges * runtimeDofPerEdge() );
    edge_dof.reserve( ndofF );
    
    uint16_type lcVertex = 0;
    uint16_type lcEdge = 0;

    addVertexEdgeDof( elid, edge_id, lcVertex, std::back_inserter( edge_dof ) );
    DVLOG(3) << "n local dof vertex for edge (" << elid << "," << edge_id << ")=" << lcVertex;
    addEdgeEdgeDof( elid, edge_id, lcEdge, std::back_inserter(edge_dof) );
    DVLOG(3) << "n local dof edge for edge (" << elid << "," << edge_id << ")=" << lcVertex;
    return edge_dof;
}
}
#endif /* FEELPP_DofFromEdge_H */
