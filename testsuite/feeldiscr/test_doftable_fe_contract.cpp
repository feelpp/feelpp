/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2026-05-09

  Copyright (C) 2026 Feel++ Consortium

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
 * @file test_doftable_fe_contract.cpp
 * @brief Compatibility tests for the public FE/DofTable contract.
 */
#define BOOST_TEST_MODULE test_doftable_fe_contract
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/doflayout.hpp>
#include <feel/feeldiscr/doftablebase.hpp>
#include <feel/feeldiscr/neh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/brezzidouglasmarini.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

#include <array>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
using contract_mesh_type = Mesh<Simplex<2, 1>>;
using contract_mesh_ptrtype = std::shared_ptr<contract_mesh_type>;

contract_mesh_ptrtype
makeContractMesh()
{
    return unitSquare( 0.35 );
}

auto
CRh( contract_mesh_ptrtype const& mesh )
{
    using space_type = FunctionSpace<contract_mesh_type, bases<CrouzeixRaviart<1>>>;
    return space_type::New( mesh );
}

template<typename SpacePtrType>
void
checkElementPublicContract( SpacePtrType const& Xh, std::string const& label )
{
    using space_type = typename SpacePtrType::element_type;
    using size_type = typename space_type::size_type;
    using doftable_base_type = DofTableBase<size_type>;

    auto const mesh = Xh->mesh();
    auto const dof = Xh->dof();
    auto const& fe = *Xh->basis();
    auto const& dofBase = static_cast<doftable_base_type const&>( *dof );

    size_type checkedElements = 0;
    auto rangeElements = elements( mesh, entity_process_t::LOCAL_ONLY );
    for ( auto const& eltWrap : rangeElements )
    {
        auto const& elt = unwrap_ref( eltWrap );
        auto const elid = elt.id();
        auto const nLocalPerComponent = static_cast<uint16_type>( dof->nLocalDof( true ) );
        auto const nLocalWithComponents = dof->getIndicesSize( elid );

        BOOST_TEST_CONTEXT( label << " element " << elid )
        {
            BOOST_REQUIRE_GT( nLocalPerComponent, 0 );
            BOOST_REQUIRE_GE( nLocalWithComponents, nLocalPerComponent );
            BOOST_CHECK_EQUAL( nLocalPerComponent, fe.localDofCount( true ) );
            BOOST_CHECK_EQUAL( nLocalWithComponents, fe.localDofCount() );
            BOOST_CHECK_EQUAL( dof->nLocalDof(), fe.localDofCount() );

            auto const indices = dof->getIndices( elid );
            BOOST_REQUIRE_EQUAL( indices.size(), nLocalWithComponents );

            auto const clusterIndices = dof->getIndicesOnGlobalCluster( elid );
            BOOST_REQUIRE( clusterIndices.empty() || clusterIndices.size() == nLocalWithComponents );

            auto const& localToGlobalIndices = dof->localToGlobalIndices( elid );
            BOOST_REQUIRE_EQUAL( static_cast<size_type>( localToGlobalIndices.size() ), nLocalWithComponents );

            auto const& localToGlobalSigns = dof->localToGlobalSigns( elid );
            BOOST_REQUIRE_GE( static_cast<size_type>( localToGlobalSigns.size() ), nLocalWithComponents );

            for ( uint16_type localDof = 0; localDof < nLocalPerComponent; ++localDof )
            {
                auto const& fromId = dof->localToGlobal( elid, localDof, 0 );
                auto const fromElement = dof->localToGlobal( elt, localDof, 0 );
                auto const& fromBase = dofBase.localToGlobal( elid, localDof, 0 );
                auto const onCluster = dof->localToGlobalOnCluster( elid, localDof, 0 );

                BOOST_CHECK_EQUAL( fromId.index(), fromElement.index() );
                BOOST_CHECK_EQUAL( fromId.index(), fromBase.index() );
                BOOST_CHECK_LT( fromId.index(), dof->nLocalDofWithGhost() );
                BOOST_CHECK_EQUAL( localToGlobalIndices( localDof ), static_cast<int>( fromId.index() ) );
                BOOST_CHECK_EQUAL( onCluster.index(), dof->mapGlobalProcessToGlobalCluster( fromId.index() ) );
                BOOST_CHECK_LT( onCluster.index(), dof->nDof() );
            }

            auto const localDofs = dof->localDof( elid );
            for ( auto it = localDofs.first; it != localDofs.second; ++it )
            {
                auto const localDof = it->first.localDof();
                auto const globalProcessDof = it->second.index();
                BOOST_CHECK_LT( localDof, nLocalWithComponents );
                BOOST_CHECK_LT( globalProcessDof, dof->nLocalDofWithGhost() );
                BOOST_CHECK_LT( dof->mapGlobalProcessToGlobalCluster( globalProcessDof ), dof->nDof() );
            }
        }
        ++checkedElements;
    }

    mpi::all_reduce( Environment::worldComm(), mpi::inplace( checkedElements ), std::plus<size_type>() );
    BOOST_CHECK_GT( checkedElements, 0 );

    if ( !dof->mapGDof().empty() )
    {
        size_type flatKeyEntries = 0;
        std::set<std::tuple<uint16_type,size_type,uint16_type>> keys;
        BOOST_CHECK( dof->mapGDofUsesFlatLocalDofEntries() );
        dof->forEachGlobalDofKeyEntry( [&]( auto const& entry )
        {
            BOOST_CHECK_LT( entry.localDofIndex, dof->nLocalDofWithGhost() );
            keys.emplace( std::get<0>( entry.key ), std::get<1>( entry.key ), entry.component );
            ++flatKeyEntries;
        } );
        BOOST_CHECK_EQUAL( flatKeyEntries, static_cast<size_type>( dof->mapGDof().size() ) );
        BOOST_CHECK_EQUAL( static_cast<size_type>( keys.size() ), flatKeyEntries );
    }
}

template<typename SpacePtrType>
void
checkFacePublicContract( SpacePtrType const& Xh, std::string const& label )
{
    using space_type = typename SpacePtrType::element_type;
    using size_type = typename space_type::size_type;
    using doftable_base_type = DofTableBase<size_type>;

    auto const mesh = Xh->mesh();
    auto const dof = Xh->dof();
    auto const& dofBase = static_cast<doftable_base_type const&>( *dof );
    auto const nFaceDofPerComponent = static_cast<uint16_type>( dof->nLocalDofOnFace( true ) );

    if ( nFaceDofPerComponent == 0 )
        return;

    BOOST_CHECK_EQUAL( dof->nLocalDofOnFacet( true ), dof->nLocalDofOnFace( true ) );

    size_type checkedFaces = 0;
    for ( auto fit = mesh->beginFace(), fen = mesh->endFace(); fit != fen; ++fit )
    {
        auto const& face = fit->second;
        auto const faceDofs = dof->faceLocalDof( face.id() );
        auto const facetDofs = dof->facetLocalDof( face.id() );
        if ( faceDofs.first == faceDofs.second )
            continue;

        BOOST_TEST_CONTEXT( label << " face " << face.id() )
        {
            BOOST_REQUIRE_GE( static_cast<size_type>( std::distance( faceDofs.first, faceDofs.second ) ),
                              static_cast<size_type>( nFaceDofPerComponent ) );
            BOOST_CHECK_EQUAL( std::distance( faceDofs.first, faceDofs.second ),
                               std::distance( facetDofs.first, facetDofs.second ) );

            for ( uint16_type localFaceDof = 0; localFaceDof < nFaceDofPerComponent; ++localFaceDof )
            {
                auto const& byId = dof->faceLocalToGlobal( face.id(), localFaceDof, 0 );
                auto const& byFacetId = dof->facetLocalToGlobal( face.id(), localFaceDof, 0 );
                auto const byEntity = dof->localToGlobal( face, localFaceDof, 0 );
                auto const& byBase = dofBase.faceLocalToGlobal( face.id(), localFaceDof, 0 );
                auto const& byBaseFacet = dofBase.facetLocalToGlobal( face.id(), localFaceDof, 0 );

                BOOST_CHECK_EQUAL( byId.index(), byEntity.index() );
                BOOST_CHECK_EQUAL( byId.index(), byFacetId.index() );
                BOOST_CHECK_EQUAL( byId.index(), byBase.index() );
                BOOST_CHECK_EQUAL( byId.index(), byBaseFacet.index() );
                BOOST_CHECK_EQUAL( byId.localDofInEntity(), localFaceDof );
                BOOST_CHECK_LT( byId.index(), dof->nLocalDofWithGhost() );
                BOOST_CHECK_LT( dof->mapGlobalProcessToGlobalCluster( byId.index() ), dof->nDof() );
            }
        }
        ++checkedFaces;
    }

    mpi::all_reduce( Environment::worldComm(), mpi::inplace( checkedFaces ), std::plus<size_type>() );
    BOOST_CHECK_GT( checkedFaces, 0 );
}

template<typename SpacePtrType>
void
checkPublicDofContract( SpacePtrType const& Xh, std::string const& label )
{
    checkElementPublicContract( Xh, label );
    checkFacePublicContract( Xh, label );
}

template<typename SpacePtrType>
void
checkFeOwnedLocalCardinality( SpacePtrType const& Xh, std::string const& label )
{
    auto const dof = Xh->dof();
    auto const& fe = *Xh->basis();

    BOOST_TEST_CONTEXT( label )
    {
        BOOST_CHECK_EQUAL( dof->nLocalDof( true ), fe.localDofCount( true ) );
        BOOST_CHECK_EQUAL( dof->nLocalDof(), fe.localDofCount() );
        BOOST_CHECK_EQUAL( dof->nLocalDofOnFacet( true ), fe.localDofCountOnFacet( 0, true ) );
        BOOST_CHECK_EQUAL( dof->nLocalDofOnFacet(), fe.localDofCountOnFacet( 0 ) );
    }
}

template<typename SpacePtrType>
void
checkFlatDescriptorInsertion( SpacePtrType const& Xh, std::string const& label )
{
    using space_type = typename SpacePtrType::element_type;
    using size_type = typename space_type::size_type;

    auto const mesh = Xh->mesh();
    auto const dof = Xh->dof();
    auto const& fe = *Xh->basis();
    auto const nFlatLocalDof = fe.localDofCount();

    size_type checkedElements = 0;
    auto rangeElements = elements( mesh, entity_process_t::LOCAL_ONLY );
    for ( auto const& eltWrap : rangeElements )
    {
        auto const& elt = unwrap_ref( eltWrap );
        auto const elid = elt.id();
        auto const indices = dof->getIndices( elid );
        auto const& localToGlobalIndices = dof->localToGlobalIndices( elid );
        auto const localDofs = dof->localDof( elid );
        std::vector<bool> seen( nFlatLocalDof, false );

        BOOST_TEST_CONTEXT( label << " element " << elid )
        {
            BOOST_REQUIRE_EQUAL( dof->getIndicesSize( elid ), nFlatLocalDof );
            BOOST_REQUIRE_EQUAL( indices.size(), nFlatLocalDof );
            BOOST_REQUIRE_EQUAL( static_cast<size_type>( localToGlobalIndices.size() ), nFlatLocalDof );

            for ( auto it = localDofs.first; it != localDofs.second; ++it )
            {
                auto const localDof = it->first.localDof();
                BOOST_REQUIRE_LT( localDof, nFlatLocalDof );
                seen[localDof] = true;
            }

            for ( uint16_type localDof = 0; localDof < nFlatLocalDof; ++localDof )
            {
                auto const layout = fe.localDofLayout( localDof );
                auto const& fromElement = dof->localToGlobal( elid, layout.parentLocalDofId, layout.component );
                BOOST_CHECK( seen[localDof] );
                BOOST_CHECK_EQUAL( layout.localDofId, localDof );
                BOOST_CHECK_EQUAL( localToGlobalIndices( localDof ), static_cast<int>( fromElement.index() ) );
            }
        }
        ++checkedElements;
    }

    mpi::all_reduce( Environment::worldComm(), mpi::inplace( checkedElements ), std::plus<size_type>() );
    BOOST_CHECK_GT( checkedElements, 0 );
}

template<typename FEType>
void
checkLocalDescriptorContract( FEType const& fe, std::string const& label )
{
    static_assert( FiniteElementDofLayoutProvider<FEType> );

    auto const descriptors = makeLocalDofDescriptors( fe );
    auto const nLocalDof = finiteElementLocalDofCount( fe );
    auto const validation = validateFiniteElementDofLayout( fe );
    std::array<uint16_type, 4> countByDim = { 0, 0, 0, 0 };

    BOOST_REQUIRE_EQUAL( descriptors.size(), nLocalDof );
    BOOST_CHECK_EQUAL( fe.localDofCount(), nLocalDof );
    BOOST_CHECK_EQUAL( fe.localDofCount( true ), fe.localDofPerComponent() );
    BOOST_REQUIRE_MESSAGE( static_cast<bool>( validation ),
                           label << " invalid FE local layout failure="
                                 << static_cast<int>( validation.failure )
                                 << " localDof=" << validation.localDof
                                 << " entityDim=" << static_cast<int>( validation.entityDim )
                                 << " localEntity=" << validation.localEntity
                                 << " ordinal=" << validation.ordinal );
    BOOST_CHECK( finiteElementDofLayoutIsComplete( fe ) );

    for ( uint16_type localDof = 0; localDof < nLocalDof; ++localDof )
    {
        auto const layout = fe.localDofLayout( localDof );
        auto const& descriptor = descriptors[localDof];
        auto const expectedKind = dofEntityKindFromTopologicalDim( static_cast<int8_type>( layout.attachment.entityDim ),
                                                                   FEType::nDim );

        BOOST_TEST_CONTEXT( label << " local dof " << localDof )
        {
            BOOST_CHECK_EQUAL( descriptor.localDof, layout.localDofId );
            BOOST_CHECK_EQUAL( descriptor.localDof, localDof );
            BOOST_CHECK_EQUAL( descriptor.parentLocalDof, layout.parentLocalDofId );
            BOOST_CHECK_EQUAL( descriptor.parentLocalDof, fe.dofParent( localDof ) );
            BOOST_CHECK_EQUAL( descriptor.component, layout.component );
            BOOST_CHECK_EQUAL( descriptor.component, fe.component( localDof ) );
            BOOST_CHECK( descriptor.attachment.isValid() );
            BOOST_CHECK( descriptor.attachment.kind == expectedKind );
            BOOST_CHECK_EQUAL( static_cast<int>( descriptor.attachment.topologicalDim ),
                               static_cast<int>( layout.attachment.entityDim ) );
            BOOST_CHECK_EQUAL( descriptor.attachment.localEntity, layout.attachment.entityId );
            BOOST_CHECK_EQUAL( descriptor.attachment.ordinal, layout.attachment.ordinal );
            BOOST_CHECK_EQUAL( descriptor.attachment.shared,
                               static_cast<uint16_type>( layout.attachment.entityDim ) < FEType::nDim );
            BOOST_CHECK_EQUAL( static_cast<uint16_type>( descriptor.functional ),
                               static_cast<uint16_type>( finiteElementDofFunctionalKind( fe, localDof ) ) );

            if ( descriptor.attachment.topologicalDim >= 0 &&
                 descriptor.attachment.topologicalDim + 1 == FEType::nDim )
                BOOST_CHECK( descriptor.attachment.kind == DofEntityKind::Facet );

            BOOST_REQUIRE_GE( layout.attachment.entityDim, 0 );
            BOOST_REQUIRE_LT( static_cast<std::size_t>( layout.attachment.entityDim ), countByDim.size() );
            countByDim[layout.attachment.entityDim] += 1;
        }
    }

    auto const nDofPerVertex = [&fe]() -> uint16_type
    {
        if constexpr ( requires { fe.runtimeDofPerVertex(); } )
            return fe.runtimeDofPerVertex();
        else
            return FEType::nDofPerVertex;
    }();
    auto const nDofPerEdge = [&fe]() -> uint16_type
    {
        if constexpr ( requires { fe.runtimeDofPerEdge(); } )
            return fe.runtimeDofPerEdge();
        else
            return FEType::nDofPerEdge;
    }();
    auto const nDofPerFace = [&fe]() -> uint16_type
    {
        if constexpr ( requires { fe.runtimeDofPerFace(); } )
            return fe.runtimeDofPerFace();
        else
            return FEType::nDofPerFace;
    }();
    auto const nDofPerVolume = [&fe]() -> uint16_type
    {
        if constexpr ( requires { fe.runtimeDofPerVolume(); } )
            return fe.runtimeDofPerVolume();
        else
            return FEType::nDofPerVolume;
    }();

    BOOST_CHECK_EQUAL( countByDim[0],
                       static_cast<uint16_type>( FEType::reference_convex_type::numVertices * nDofPerVertex ) );
    BOOST_CHECK_EQUAL( countByDim[1],
                       static_cast<uint16_type>( FEType::reference_convex_type::numEdges * nDofPerEdge ) );
    BOOST_CHECK_EQUAL( countByDim[2],
                       static_cast<uint16_type>( FEType::reference_convex_type::numFaces * nDofPerFace ) );
    BOOST_CHECK_EQUAL( countByDim[3], nDofPerVolume );

    if constexpr ( FEType::nDim > 1 )
    {
        using face_type = typename FEType::face_type;
        const uint16_type expectedFacetDofPerComponent =
            static_cast<uint16_type>( face_type::numVertices * nDofPerVertex +
                                      face_type::numEdges * nDofPerEdge +
                                      face_type::numFaces * nDofPerFace );
        for ( uint16_type localFacet = 0; localFacet < FEType::reference_convex_type::numFaces; ++localFacet )
        {
            BOOST_CHECK_EQUAL( fe.localDofCountOnFacet( localFacet, true ),
                               expectedFacetDofPerComponent );
            if constexpr ( FEType::is_product )
                BOOST_CHECK_EQUAL( fe.localDofCountOnFacet( localFacet ),
                                   static_cast<uint16_type>( FEType::nComponents * expectedFacetDofPerComponent ) );
            else
                BOOST_CHECK_EQUAL( fe.localDofCountOnFacet( localFacet ),
                                   expectedFacetDofPerComponent );
        }
    }
}

template<typename FEType>
uint16_type
testDofPerVertex( FEType const& fe )
{
    if constexpr ( requires { fe.runtimeDofPerVertex(); } )
        return fe.runtimeDofPerVertex();
    else
        return FEType::nDofPerVertex;
}

template<typename FEType>
uint16_type
testDofPerEdge( FEType const& fe )
{
    if constexpr ( requires { fe.runtimeDofPerEdge(); } )
        return fe.runtimeDofPerEdge();
    else
        return FEType::nDofPerEdge;
}

template<typename FEType>
uint16_type
testDofPerFace( FEType const& fe )
{
    if constexpr ( requires { fe.runtimeDofPerFace(); } )
        return fe.runtimeDofPerFace();
    else
        return FEType::nDofPerFace;
}

template<typename FEType>
uint16_type
testDofPerVolume( FEType const& fe )
{
    if constexpr ( requires { fe.runtimeDofPerVolume(); } )
        return fe.runtimeDofPerVolume();
    else
        return FEType::nDofPerVolume;
}

template<typename FEType>
bool
testRuntimeP0Continuous( FEType const& fe )
{
    if constexpr ( !FEType::isContinuous )
        return false;
    else if constexpr ( requires { fe.order(); } )
        return fe.order() == 0;
    else
        return FEType::nOrder == 0;
}

template<typename FEType>
constexpr uint16_type
testFamilyTag()
{
    if constexpr ( requires { FEType::TAG; } )
        return FEType::TAG;
    else
        return 0;
}

template<typename ElementType, typename FaceType, typename AttachmentType>
uint16_type
testLocalFacetDofIndexFromAttachment( uint16_type localFace,
                                      AttachmentType const& attachment,
                                      uint16_type nDofPerVertex,
                                      uint16_type nDofPerEdge,
                                      uint16_type nDofPerFace )
{
    if ( !attachment.isValid() )
        return invalid_uint16_type_value;

    if ( attachment.entityDim == 0 )
    {
        if ( nDofPerVertex == 0 || attachment.ordinal >= nDofPerVertex )
            return invalid_uint16_type_value;

        for ( uint16_type localVertexInFace = 0; localVertexInFace < FaceType::numVertices; ++localVertexInFace )
        {
            auto const localVertexInElement = ElementType::fToP( localFace, localVertexInFace );
            if ( localVertexInElement == attachment.entityId )
                return static_cast<uint16_type>( localVertexInFace * nDofPerVertex + attachment.ordinal );
        }
        return invalid_uint16_type_value;
    }

    if ( attachment.entityDim == 1 )
    {
        if ( nDofPerEdge == 0 || attachment.ordinal >= nDofPerEdge )
            return invalid_uint16_type_value;

        if constexpr ( ElementType::nDim == 2 )
        {
            if ( attachment.entityId == localFace )
                return static_cast<uint16_type>( FaceType::numVertices * nDofPerVertex + attachment.ordinal );
        }
        else if constexpr ( ElementType::nDim == 3 )
        {
            for ( uint16_type localEdgeInFace = 0; localEdgeInFace < FaceType::numEdges; ++localEdgeInFace )
            {
                auto const localEdgeInElement = ElementType::fToE( localFace, localEdgeInFace );
                if ( localEdgeInElement == attachment.entityId )
                    return static_cast<uint16_type>( FaceType::numVertices * nDofPerVertex +
                                                     localEdgeInFace * nDofPerEdge + attachment.ordinal );
            }
        }
        return invalid_uint16_type_value;
    }

    if ( attachment.entityDim == 2 )
    {
        if ( nDofPerFace == 0 || attachment.ordinal >= nDofPerFace )
            return invalid_uint16_type_value;

        if constexpr ( ElementType::nDim == 3 )
            if ( attachment.entityId == localFace )
                return static_cast<uint16_type>( FaceType::numVertices * nDofPerVertex +
                                                 FaceType::numEdges * nDofPerEdge + attachment.ordinal );
    }

    return invalid_uint16_type_value;
}

template<typename ElementType, typename EdgeType, typename AttachmentType>
uint16_type
testLocalEdgeDofIndexFromAttachment( uint16_type localEdge,
                                     AttachmentType const& attachment,
                                     uint16_type nDofPerVertex,
                                     uint16_type nDofPerEdge )
{
    if ( !attachment.isValid() )
        return invalid_uint16_type_value;

    if ( attachment.entityDim == 0 )
    {
        if ( nDofPerVertex == 0 || attachment.ordinal >= nDofPerVertex )
            return invalid_uint16_type_value;

        for ( uint16_type localVertexInEdge = 0; localVertexInEdge < EdgeType::numVertices; ++localVertexInEdge )
        {
            auto const localVertexInElement = ElementType::eToP( localEdge, localVertexInEdge );
            if ( localVertexInElement == attachment.entityId )
                return static_cast<uint16_type>( localVertexInEdge * nDofPerVertex + attachment.ordinal );
        }
        return invalid_uint16_type_value;
    }

    if ( attachment.entityDim == 1 )
    {
        if ( nDofPerEdge == 0 || attachment.ordinal >= nDofPerEdge )
            return invalid_uint16_type_value;

        if ( attachment.entityId == localEdge )
            return static_cast<uint16_type>( EdgeType::numVertices * nDofPerVertex + attachment.ordinal );
    }

    return invalid_uint16_type_value;
}

template<typename SpacePtrType>
void
checkDescriptorBackedEntityViews( SpacePtrType const& Xh, std::string const& label )
{
    using space_type = typename SpacePtrType::element_type;
    using dof_type = typename space_type::dof_type;
    using size_type = typename space_type::size_type;
    using element_type = typename dof_type::element_type;
    using face_type = typename dof_type::face_type;
    using edge_type = typename dof_type::edge_type;
    using fe_type = std::remove_cvref_t<decltype( *Xh->basis() )>;

    auto const mesh = Xh->mesh();
    auto const dof = Xh->dof();
    auto const& fe = *Xh->basis();

    auto const nDofPerVertex = testDofPerVertex( fe );
    auto const nDofPerEdge = testDofPerEdge( fe );
    auto const nDofPerFace = testDofPerFace( fe );
    auto const nFacetDof = static_cast<uint16_type>( dof->nLocalDofOnFacet( true ) );
    auto const nEdgeDof = static_cast<uint16_type>( edge_type::numVertices * nDofPerVertex +
                                                   edge_type::numEdges * nDofPerEdge );

    size_type checkedFacets = 0;
    for ( auto fit = mesh->beginFace(), fen = mesh->endFace(); fit != fen; ++fit )
    {
        auto const& face = fit->second;
        size_type elementId = invalid_v<size_type>;
        uint16_type localFace = invalid_uint16_type_value;

        if ( face.isConnectedTo0() && dof->isElementDone( face.ad_first() ) )
        {
            elementId = face.ad_first();
            localFace = face.pos_first();
        }
        else if ( face.isConnectedTo1() && dof->isElementDone( face.ad_second() ) )
        {
            elementId = face.ad_second();
            localFace = face.pos_second();
        }

        if ( elementId == invalid_v<size_type> || localFace == invalid_uint16_type_value )
            continue;

        auto const facetDofs = dof->facetLocalDof( face.id() );
        if ( facetDofs.first == facetDofs.second )
            continue;

        std::vector<uint16_type> expected( nFacetDof, invalid_uint16_type_value );
        for ( uint16_type parentLid = 0; parentLid < fe.localDofPerComponent(); ++parentLid )
        {
            auto const layout = fe.localDofLayout( parentLid );
            auto const localIndex = testLocalFacetDofIndexFromAttachment<element_type, face_type>(
                localFace, layout.attachment, nDofPerVertex, nDofPerEdge, nDofPerFace );
            if ( localIndex != invalid_uint16_type_value )
                expected[localIndex] = dof->localDofId( parentLid, 0 );
        }

        for ( uint16_type localFacetDof = 0; localFacetDof < nFacetDof; ++localFacetDof )
        {
            BOOST_REQUIRE_MESSAGE( expected[localFacetDof] != invalid_uint16_type_value,
                                   label << " missing descriptor for facet dof " << localFacetDof );
            auto const& facetDof = dof->facetLocalToGlobal( face.id(), localFacetDof, 0 );
            auto const& elementDof = dof->localToGlobal( elementId, expected[localFacetDof], 0 );
            BOOST_CHECK_EQUAL( facetDof.localDofInElement(), expected[localFacetDof] );
            BOOST_CHECK_EQUAL( facetDof.index(), elementDof.index() );
        }
        ++checkedFacets;
    }

    size_type checkedEdges = 0;
    auto rangeElements = elements( mesh, entity_process_t::LOCAL_ONLY );
    for ( auto const& eltWrap : rangeElements )
    {
        auto const& elt = unwrap_ref( eltWrap );
        for ( uint16_type localEdge = 0; localEdge < element_type::numEdges; ++localEdge )
        {
            auto const edgeDofs = dof->edgeLocalDof( elt.id(), localEdge );
            if ( edgeDofs.empty() )
                continue;

            BOOST_REQUIRE_GE( static_cast<size_type>( edgeDofs.size() ), static_cast<size_type>( nEdgeDof ) );

            std::vector<uint16_type> expected( nEdgeDof, invalid_uint16_type_value );
            for ( uint16_type parentLid = 0; parentLid < fe.localDofPerComponent(); ++parentLid )
            {
                auto const layout = fe.localDofLayout( parentLid );
                auto const localIndex = testLocalEdgeDofIndexFromAttachment<element_type, edge_type>(
                    localEdge, layout.attachment, nDofPerVertex, nDofPerEdge );
                if ( localIndex != invalid_uint16_type_value )
                    expected[localIndex] = dof->localDofId( parentLid, 0 );
            }

            for ( uint16_type localEdgeDof = 0; localEdgeDof < nEdgeDof; ++localEdgeDof )
            {
                BOOST_REQUIRE_MESSAGE( expected[localEdgeDof] != invalid_uint16_type_value,
                                       label << " missing descriptor for edge dof " << localEdgeDof );
                auto const& edgeDof = edgeDofs[localEdgeDof];
                auto const& elementDof = dof->localToGlobal( elt.id(), expected[localEdgeDof], 0 );
                BOOST_CHECK_EQUAL( edgeDof.localDofInElement(), expected[localEdgeDof] );
                BOOST_CHECK_EQUAL( edgeDof.index(), elementDof.index() );
            }
            ++checkedEdges;
        }
    }

    mpi::all_reduce( Environment::worldComm(), mpi::inplace( checkedFacets ), std::plus<size_type>() );
    mpi::all_reduce( Environment::worldComm(), mpi::inplace( checkedEdges ), std::plus<size_type>() );
    BOOST_CHECK_GT( checkedFacets, 0 );
    BOOST_CHECK_GT( checkedEdges, 0 );
}

template<typename SpacePtrType>
void
checkDescriptorKeyLocalToGlobalConsistency( SpacePtrType const& Xh, std::string const& label )
{
    using space_type = typename SpacePtrType::element_type;
    using dof_type = typename space_type::dof_type;
    using size_type = typename space_type::size_type;
    using element_type = typename dof_type::element_type;
    using edge_permutation_type = typename element_type::edge_permutation_type;
    using face_permutation_type = typename element_type::face_permutation_type;
    using fe_type = std::remove_cvref_t<decltype( *Xh->basis() )>;

    auto const mesh = Xh->mesh();
    auto const dof = Xh->dof();
    auto const& fe = *Xh->basis();

    auto const nDofPerVertex = testDofPerVertex( fe );
    auto const nDofPerEdge = testDofPerEdge( fe );
    auto const nDofPerFace = testDofPerFace( fe );
    auto const nDofPerVolume = testDofPerVolume( fe );
    auto const p0Continuous = testRuntimeP0Continuous( fe );

    std::map<DofKey<size_type>, size_type> keyToIndex;
    std::map<size_type, DofKey<size_type>> indexToKey;
    size_type checkedDofs = 0;

    auto makeKey = []( uint16_type topologicalDim,
                       size_type canonicalEntityId,
                       uint16_type ordinal,
                       uint16_type component,
                       uint16_type functional )
    {
        return DofKey<size_type>{
            .topologicalDim = static_cast<uint8_type>( topologicalDim ),
            .canonicalEntityId = canonicalEntityId,
            .ordinal = ordinal,
            .component = component,
            .familyTag = testFamilyTag<fe_type>(),
            .variant = functional
        };
    };

    auto rangeElements = elements( mesh, entity_process_t::LOCAL_ONLY );
    for ( auto const& eltWrap : rangeElements )
    {
        auto const& elt = unwrap_ref( eltWrap );
        for ( uint16_type parentLid = 0; parentLid < fe.localDofPerComponent(); ++parentLid )
        {
            auto const layout = fe.localDofLayout( parentLid );
            auto const& attachment = layout.attachment;
            BOOST_REQUIRE( attachment.isValid() );

            DofKey<size_type> key;
            switch ( attachment.entityDim )
            {
            case 0:
                BOOST_REQUIRE_GT( nDofPerVertex, 0 );
                BOOST_REQUIRE_LT( attachment.entityId, element_type::numVertices );
                BOOST_REQUIRE_LT( attachment.ordinal, nDofPerVertex );
                key = makeKey( 0,
                               mesh->canonicalPointId( elt.point( attachment.entityId ).id() ),
                               attachment.ordinal,
                               layout.component,
                               attachment.kind );
                break;
            case 1:
                BOOST_REQUIRE_GT( nDofPerEdge, 0 );
                BOOST_REQUIRE_LT( attachment.ordinal, nDofPerEdge );
                if constexpr ( element_type::nDim == 1 )
                {
                    key = makeKey( 1,
                                   p0Continuous ? 0 : elt.id(),
                                   attachment.ordinal,
                                   layout.component,
                                   attachment.kind );
                }
                else
                {
                    BOOST_REQUIRE_LT( attachment.entityId, element_type::numEdges );
                    auto ordinal = attachment.ordinal;
                    if ( elt.edgePermutation( attachment.entityId ).value() == edge_permutation_type::REVERSE_PERMUTATION )
                        ordinal = static_cast<uint16_type>( nDofPerEdge - 1 - attachment.ordinal );
                    key = makeKey( 1,
                                   mesh->canonicalEdgeId( elt.edge( attachment.entityId ).id() ),
                                   ordinal,
                                   layout.component,
                                   attachment.kind );
                }
                break;
            case 2:
                BOOST_REQUIRE_GT( nDofPerFace, 0 );
                BOOST_REQUIRE_LT( attachment.ordinal, nDofPerFace );
                if constexpr ( element_type::nDim == 2 )
                {
                    key = makeKey( 2,
                                   p0Continuous ? 0 : elt.id(),
                                   attachment.ordinal,
                                   layout.component,
                                   attachment.kind );
                }
                else if constexpr ( element_type::nDim == 3 )
                {
                    BOOST_REQUIRE_LT( attachment.entityId, element_type::numFaces );
                    auto ordinal = attachment.ordinal;
                    auto const permutation = elt.facePermutation( attachment.entityId );
                    if ( nDofPerFace != 1 && permutation != face_permutation_type( face_permutation_type::IDENTITY ) )
                        ordinal = static_cast<uint16_type>( dof->facePermutationVector( permutation, nDofPerFace )( attachment.ordinal ) );
                    key = makeKey( 2,
                                   mesh->canonicalFaceId( elt.face( attachment.entityId ).id() ),
                                   ordinal,
                                   layout.component,
                                   attachment.kind );
                }
                break;
            case 3:
                BOOST_REQUIRE_GT( nDofPerVolume, 0 );
                BOOST_REQUIRE_LT( attachment.ordinal, nDofPerVolume );
                key = makeKey( 3,
                               p0Continuous ? 0 : elt.id(),
                               attachment.ordinal,
                               layout.component,
                               attachment.kind );
                break;
            default:
                BOOST_FAIL( "invalid descriptor entity dimension" );
            }

            BOOST_REQUIRE( key.isValid() );
            auto const index = dof->localToGlobal( elt.id(), parentLid, 0 ).index();
            auto const [keyIt, insertedKey] = keyToIndex.emplace( key, index );
            if ( !insertedKey )
                BOOST_CHECK_EQUAL( keyIt->second, index );

            auto const [indexIt, insertedIndex] = indexToKey.emplace( index, key );
            if ( !insertedIndex )
                BOOST_CHECK( indexIt->second == key );

            BOOST_TEST_CONTEXT( label << " element " << elt.id() << " local dof " << parentLid )
            {
                BOOST_CHECK_LT( index, dof->nLocalDofWithGhost() );
            }
            ++checkedDofs;
        }
    }

    mpi::all_reduce( Environment::worldComm(), mpi::inplace( checkedDofs ), std::plus<size_type>() );
    BOOST_CHECK_GT( checkedDofs, 0 );
}

template<typename SpacePtrType>
void
checkDescriptorKeyMpiOwnership( SpacePtrType const& Xh, std::string const& label )
{
    using space_type = typename SpacePtrType::element_type;
    using dof_type = typename space_type::dof_type;
    using size_type = typename space_type::size_type;
    using element_type = typename dof_type::element_type;
    using edge_permutation_type = typename element_type::edge_permutation_type;
    using face_permutation_type = typename element_type::face_permutation_type;
    using fe_type = std::remove_cvref_t<decltype( *Xh->basis() )>;

    if ( Environment::worldComm().localSize() <= 1 )
        return;

    auto const mesh = Xh->mesh();
    auto const dof = Xh->dof();
    auto const& fe = *Xh->basis();
    auto const myRank = Environment::worldComm().localRank();

    auto const nDofPerVertex = testDofPerVertex( fe );
    auto const nDofPerEdge = testDofPerEdge( fe );
    auto const nDofPerFace = testDofPerFace( fe );

    auto familyTag = []() constexpr -> uint16_type
    {
        if constexpr ( requires { fe_type::TAG; } )
            return fe_type::TAG;
        else
            return 0;
    };

    auto partitionCanonicalEntity = []( auto const& entity )
    {
        rank_type canonicalPartition = entity.pidInPartition();
        if ( canonicalPartition == invalid_rank_type_value )
            canonicalPartition = entity.processId();

        size_type canonicalId = entity.id();
        if ( canonicalPartition != entity.pidInPartition() )
        {
            auto const itId = entity.idInOthersPartitions().find( canonicalPartition );
            if ( itId != entity.idInOthersPartitions().end() && itId->second != invalid_v<size_type> )
                canonicalId = itId->second;
        }

        for ( auto const& [pid,entityId] : entity.idInOthersPartitions() )
        {
            if ( entityId != invalid_v<size_type> && pid < canonicalPartition )
            {
                canonicalPartition = pid;
                canonicalId = entityId;
            }
        }
        return std::make_pair( canonicalPartition, canonicalId );
    };

    auto makeKey = [&familyTag]( uint16_type topologicalDim,
                                 rank_type canonicalPartition,
                                 size_type canonicalEntityId,
                                 uint16_type ordinal,
                                 uint16_type component,
                                 uint16_type functional )
    {
        return DofKey<size_type>{
            .topologicalDim = static_cast<uint8_type>( topologicalDim ),
            .canonicalPartition = canonicalPartition,
            .canonicalEntityId = canonicalEntityId,
            .ordinal = ordinal,
            .component = component,
            .familyTag = familyTag(),
            .variant = functional
        };
    };

    auto packKey = []( DofKey<size_type> const& key, std::vector<size_type>& payload )
    {
        payload.push_back( static_cast<size_type>( key.topologicalDim ) );
        payload.push_back( static_cast<size_type>( key.canonicalPartition ) );
        payload.push_back( key.canonicalEntityId );
        payload.push_back( static_cast<size_type>( key.ordinal ) );
        payload.push_back( static_cast<size_type>( key.component ) );
        payload.push_back( static_cast<size_type>( key.familyTag ) );
        payload.push_back( static_cast<size_type>( key.variant ) );
    };

    auto unpackKey = []( std::vector<size_type> const& payload, std::size_t offset )
    {
        return DofKey<size_type>{
            .topologicalDim = static_cast<uint8_type>( payload[offset] ),
            .canonicalPartition = static_cast<rank_type>( payload[offset+1] ),
            .canonicalEntityId = payload[offset+2],
            .ordinal = static_cast<uint16_type>( payload[offset+3] ),
            .component = static_cast<uint16_type>( payload[offset+4] ),
            .familyTag = static_cast<uint16_type>( payload[offset+5] ),
            .variant = static_cast<uint16_type>( payload[offset+6] )
        };
    };

    auto makeSharedKey = [&]( element_type const& elt,
                              uint16_type parentLid,
                              DofKey<size_type>& key ) -> bool
    {
        auto const layout = fe.localDofLayout( parentLid );
        auto const& attachment = layout.attachment;
        if ( !attachment.isValid() )
            return false;

        switch ( attachment.entityDim )
        {
        case 0:
        {
            if ( nDofPerVertex == 0 || attachment.entityId >= element_type::numVertices ||
                 attachment.ordinal >= nDofPerVertex )
                return false;

            auto const& point = elt.point( attachment.entityId );
            if ( !mesh->isInterprocessPoints( point.id() ) && point.idInOthersPartitions().empty() )
                return false;

            auto const [entityPartition,entityId] = partitionCanonicalEntity( point );
            key = makeKey( 0, entityPartition, entityId,
                           attachment.ordinal, layout.component, attachment.kind );
            return true;
        }
        case 1:
        {
            if ( nDofPerEdge == 0 || attachment.ordinal >= nDofPerEdge )
                return false;

            if constexpr ( element_type::nDim == 1 )
                return false;
            else
            {
                if ( attachment.entityId >= element_type::numEdges )
                    return false;

                auto ordinal = attachment.ordinal;
                auto const edgePermutation = elt.edgePermutation( attachment.entityId );
                if ( edgePermutation.value() == edge_permutation_type::REVERSE_PERMUTATION )
                    ordinal = static_cast<uint16_type>( nDofPerEdge - 1 - attachment.ordinal );
                else if ( edgePermutation.value() != edge_permutation_type::IDENTITY )
                    return false;

                if constexpr ( element_type::nDim == 2 )
                {
                    auto const facePtr = elt.facePtr( attachment.entityId );
                    if ( !facePtr || ( !facePtr->isInterProcessDomain() && facePtr->idInOthersPartitions().empty() ) )
                        return false;

                    auto const [entityPartition,entityId] = partitionCanonicalEntity( *facePtr );
                    key = makeKey( 1, entityPartition, entityId,
                                   ordinal, layout.component, attachment.kind );
                    return true;
                }
                else
                {
                    auto const edgePtr = elt.edgePtr( attachment.entityId );
                    if ( !edgePtr || ( !mesh->isInterprocessEdges( edgePtr->id() ) && edgePtr->idInOthersPartitions().empty() ) )
                        return false;

                    auto const [entityPartition,entityId] = partitionCanonicalEntity( *edgePtr );
                    key = makeKey( 1, entityPartition, entityId,
                                   ordinal, layout.component, attachment.kind );
                    return true;
                }
            }
        }
        case 2:
        {
            if ( nDofPerFace == 0 || attachment.ordinal >= nDofPerFace )
                return false;

            if constexpr ( element_type::nDim == 3 )
            {
                if ( attachment.entityId >= element_type::numFaces )
                    return false;

                auto const facePtr = elt.facePtr( attachment.entityId );
                if ( !facePtr || ( !facePtr->isInterProcessDomain() && facePtr->idInOthersPartitions().empty() ) )
                    return false;

                auto ordinal = attachment.ordinal;
                auto const permutation = elt.facePermutation( attachment.entityId );
                if ( permutation == face_permutation_type( 0 ) )
                    return false;
                if ( nDofPerFace != 1 && permutation != face_permutation_type( face_permutation_type::IDENTITY ) )
                    ordinal = static_cast<uint16_type>( dof->facePermutationVector( permutation, nDofPerFace )( attachment.ordinal ) );

                auto const [entityPartition,entityId] = partitionCanonicalEntity( *facePtr );
                key = makeKey( 2, entityPartition, entityId,
                               ordinal, layout.component, attachment.kind );
                return true;
            }
            else
                return false;
        }
        default:
            return false;
        }
    };

    std::map<DofKey<size_type>, size_type> localSharedKeyToDof;
    auto rangeElements = elements( mesh, entity_process_t::LOCAL_ONLY );
    for ( auto const& eltWrap : rangeElements )
    {
        auto const& elt = unwrap_ref( eltWrap );
        for ( uint16_type parentLid = 0; parentLid < fe.localDofPerComponent(); ++parentLid )
        {
            DofKey<size_type> key;
            if ( !makeSharedKey( elt, parentLid, key ) )
                continue;

            auto const layout = fe.localDofLayout( parentLid );
            auto const localDofIndex = dof->localToGlobal( elt.id(), layout.parentLocalDofId, layout.component ).index();
            auto const [it, inserted] = localSharedKeyToDof.emplace( key, localDofIndex );
            if ( !inserted )
                BOOST_CHECK_EQUAL( it->second, localDofIndex );
        }
    }

    std::vector<size_type> localPayload;
    localPayload.reserve( localSharedKeyToDof.size()*10 );
    for ( auto const& [key,localDofIndex] : localSharedKeyToDof )
    {
        packKey( key, localPayload );
        localPayload.push_back( localDofIndex );
        localPayload.push_back( dof->mapGlobalProcessToGlobalCluster( localDofIndex ) );
        localPayload.push_back( static_cast<size_type>( dof->dofGlobalProcessIsGhost( localDofIndex ) ? 1 : 0 ) );
    }

    std::vector<std::vector<size_type>> gatheredPayloads;
    mpi::all_gather( Environment::worldComm(), localPayload, gatheredPayloads );

    struct key_observation
    {
        rank_type rank;
        size_type localDof;
        size_type globalClusterDof;
        bool isGhost;
    };

    std::map<DofKey<size_type>, std::vector<key_observation>> observationsByKey;
    for ( rank_type rank = 0; rank < static_cast<rank_type>( gatheredPayloads.size() ); ++rank )
    {
        auto const& payload = gatheredPayloads[rank];
        BOOST_REQUIRE_MESSAGE( payload.size() % 10 == 0,
                               label << " invalid phase4 payload size " << payload.size() << " from rank " << rank );
        for ( std::size_t k = 0; k < payload.size(); k += 10 )
        {
            auto const key = unpackKey( payload, k );
            observationsByKey[key].push_back( key_observation{
                .rank = rank,
                .localDof = payload[k+7],
                .globalClusterDof = payload[k+8],
                .isGhost = payload[k+9] != 0
            } );
        }
    }

    size_type checkedSharedKeys = 0;
    for ( auto const& [key,observations] : observationsByKey )
    {
        std::set<rank_type> ranks;
        for ( auto const& observation : observations )
            ranks.insert( observation.rank );

        if ( ranks.size() <= 1 )
            continue;

        auto const expectedOwner = *ranks.begin();
        auto const expectedGlobalClusterDof = observations.front().globalClusterDof;

        BOOST_TEST_CONTEXT( label << " key dim=" << static_cast<int>( key.topologicalDim )
                            << " part=" << key.canonicalPartition
                            << " id=" << key.canonicalEntityId
                            << " ordinal=" << key.ordinal )
        {
            BOOST_REQUIRE_NE( expectedGlobalClusterDof, invalid_v<size_type> );
            BOOST_CHECK_EQUAL( dof->procOnGlobalCluster( expectedGlobalClusterDof ), expectedOwner );

            for ( auto const& observation : observations )
            {
                BOOST_CHECK_EQUAL( observation.globalClusterDof, expectedGlobalClusterDof );
                BOOST_CHECK_EQUAL( observation.isGhost, observation.rank != expectedOwner );
                if ( observation.rank == myRank )
                    BOOST_CHECK_LT( observation.localDof, dof->nLocalDofWithGhost() );
            }
        }
        ++checkedSharedKeys;
    }

    BOOST_CHECK_GT( checkedSharedKeys, 0 );
}

template<typename SpacePtrType>
void
checkTransformProjectionContract( SpacePtrType const& Xh, std::string const& label )
{
    using space_type = typename SpacePtrType::element_type;
    using dof_type = typename space_type::dof_type;
    using size_type = typename space_type::size_type;
    using element_type = typename dof_type::element_type;
    using fe_type = std::remove_cvref_t<decltype( *Xh->basis() )>;

    static_assert( FiniteElementDofTransformProvider<fe_type, element_type> );

    auto const mesh = Xh->mesh();
    auto const dof = Xh->dof();
    auto const& fe = *Xh->basis();

    size_type checkedDofs = 0;
    size_type negativeTransforms = 0;
    auto rangeElements = elements( mesh, entity_process_t::LOCAL_ONLY );
    for ( auto const& eltWrap : rangeElements )
    {
        auto const& elt = unwrap_ref( eltWrap );
        auto const& transforms = dof->localToGlobalTransforms( elt.id() );
        auto const& signs = dof->localToGlobalSigns( elt.id() );
        BOOST_REQUIRE_GE( transforms.size(), static_cast<std::size_t>( fe.localDofPerComponent() ) );
        BOOST_REQUIRE_GE( static_cast<size_type>( signs.size() ), static_cast<size_type>( fe.localDofPerComponent() ) );

        for ( uint16_type parentLid = 0; parentLid < fe.localDofPerComponent(); ++parentLid )
        {
            auto const layout = fe.localDofLayout( parentLid );
            auto const localDof = dof->localDofId( layout.parentLocalDofId, layout.component );
            auto const expected = fe.dofTransform( elt, parentLid );
            auto const& actual = dof->localToGlobalTransform( elt.id(), localDof );
            auto const actualSign = dof->dofTransformSignProjection( actual );
            auto const expectedSign = dof->dofTransformSignProjection( expected );

            BOOST_TEST_CONTEXT( label << " element " << elt.id() << " local dof " << localDof )
            {
                BOOST_CHECK_EQUAL( static_cast<int>( actual.kind ), static_cast<int>( expected.kind ) );
                BOOST_CHECK_EQUAL( actual.sign, expected.sign );
                BOOST_CHECK_EQUAL( actualSign, expectedSign );
                BOOST_CHECK_EQUAL( signs( localDof ), actualSign );
            }

            if ( actualSign < 0 )
                ++negativeTransforms;
            ++checkedDofs;
        }
    }

    mpi::all_reduce( Environment::worldComm(), mpi::inplace( checkedDofs ), std::plus<size_type>() );
    mpi::all_reduce( Environment::worldComm(), mpi::inplace( negativeTransforms ), std::plus<size_type>() );
    BOOST_CHECK_GT( checkedDofs, 0 );
    BOOST_CHECK_GT( negativeTransforms, 0 );
}

template<typename SpacePtrType>
void
checkRepresentativePointContract( SpacePtrType const& Xh,
                                  std::string const& label,
                                  bool expectRepresentativePoints,
                                  DofFunctionalKind expectedSharedMomentKind = DofFunctionalKind::Other )
{
    using space_type = typename SpacePtrType::element_type;
    using size_type = typename space_type::size_type;
    using fe_type = std::remove_cvref_t<decltype( *Xh->basis() )>;

    static_assert( FiniteElementDofRepresentativePointProvider<fe_type> );
    static_assert( FiniteElementDofFunctionalProvider<fe_type> );

    auto const mesh = Xh->mesh();
    auto const dof = Xh->dof();
    auto const& fe = *Xh->basis();
    auto const& dofPoints = dof->dofPoints();

    std::set<size_type> checkedGlobalDofs;
    size_type checkedLocalDofs = 0;
    size_type checkedPointDofs = 0;

    auto rangeElements = elements( mesh, entity_process_t::LOCAL_ONLY );
    for ( auto const& eltWrap : rangeElements )
    {
        auto const& elt = unwrap_ref( eltWrap );
        auto const localDofs = dof->localDof( elt.id() );
        for ( auto it = localDofs.first; it != localDofs.second; ++it )
        {
            auto const localDof = it->first.localDof();
            auto const globalProcessDof = it->second.index();
            auto const layout = fe.localDofLayout( localDof );
            auto const functional = finiteElementDofFunctionalKind( fe, localDof );

            BOOST_TEST_CONTEXT( label << " element " << elt.id() << " local dof " << localDof )
            {
                BOOST_CHECK_EQUAL( fe.dofHasRepresentativePoint( localDof ), expectRepresentativePoints );
                BOOST_CHECK_EQUAL( static_cast<uint16_type>( functional ), fe.dofFunctionalKind( localDof ) );

                if ( expectRepresentativePoints )
                {
                    auto const pointId = fe.dofRepresentativePointIndex( localDof );
                    BOOST_CHECK_LT( static_cast<int>( pointId ), static_cast<int>( fe.points().size2() ) );
                    BOOST_REQUIRE( dof->hasDofPoint( globalProcessDof ) );

                    if ( checkedGlobalDofs.insert( globalProcessDof ).second )
                    {
                        auto const& dofPoint = dof->dofPoint( globalProcessDof );
                        BOOST_CHECK_EQUAL( boost::get<1>( dofPoint ), globalProcessDof );
                        BOOST_CHECK_EQUAL( boost::get<2>( dofPoint ), layout.component );
                        ++checkedPointDofs;
                    }
                }
                else
                {
                    BOOST_CHECK( !dof->hasDofPoint( globalProcessDof ) );
                    if ( layout.attachment.isValid() )
                    {
                        if ( static_cast<uint16_type>( layout.attachment.entityDim ) == fe_type::nDim )
                            BOOST_CHECK( functional == DofFunctionalKind::InteriorMoment );
                        else if ( expectedSharedMomentKind != DofFunctionalKind::Other )
                            BOOST_CHECK( functional == expectedSharedMomentKind );
                    }
                }
            }
            ++checkedLocalDofs;
        }
    }

    auto nDofPoints = static_cast<size_type>( dofPoints.size() );
    mpi::all_reduce( Environment::worldComm(), mpi::inplace( checkedLocalDofs ), std::plus<size_type>() );
    mpi::all_reduce( Environment::worldComm(), mpi::inplace( checkedPointDofs ), std::plus<size_type>() );
    mpi::all_reduce( Environment::worldComm(), mpi::inplace( nDofPoints ), std::plus<size_type>() );

    BOOST_CHECK_GT( checkedLocalDofs, 0 );
    if ( expectRepresentativePoints )
    {
        BOOST_CHECK_GT( checkedPointDofs, 0 );
        BOOST_CHECK_GT( nDofPoints, 0 );
    }
    else
    {
        BOOST_CHECK_EQUAL( checkedPointDofs, 0 );
        BOOST_CHECK_EQUAL( nDofPoints, 0 );
    }
}

} // namespace

BOOST_AUTO_TEST_SUITE( doftable_fe_contract_suite )

BOOST_AUTO_TEST_CASE( phase0_lagrange_public_contract )
{
    auto mesh = makeContractMesh();

    checkPublicDofContract( Pch<1>( mesh ), "Pch<1>" );
    checkPublicDofContract( Pch<Dynamic>( mesh, RuntimeOrder{ 1 } ), "Pch<Dynamic,P1>" );
    checkPublicDofContract( Pdh<1>( mesh ), "Pdh<1>" );
    checkPublicDofContract( CRh( mesh ), "CrouzeixRaviart<1>" );
}

BOOST_AUTO_TEST_CASE( phase0_hdiv_hcurl_public_contract )
{
    auto mesh = makeContractMesh();

    checkPublicDofContract( RTh<0>( mesh ), "RTh<0>" );
    checkPublicDofContract( RTh<1>( mesh ), "RTh<1>" );
    checkPublicDofContract( RTh<Dynamic>( mesh, RuntimeOrder{ 1 } ), "RTh<Dynamic,P1>" );
    checkPublicDofContract( BDMh<0>( mesh ), "BDMh<0>" );
    checkPublicDofContract( BDMh<Dynamic>( mesh, RuntimeOrder{ 1 } ), "BDMh<Dynamic,P1>" );
    checkPublicDofContract( Neh<0>( mesh ), "Neh<0>" );
}

BOOST_AUTO_TEST_CASE( phase1_fe_descriptor_contract )
{
    using lagrange_p2_type = typename Lagrange<2, Scalar>::template apply<2, 2, double, Simplex<2>>::type;
    using lagrange_dynamic_type = typename Lagrange<Dynamic, Scalar>::template apply<2, 2, double, Simplex<2>>::type;
    using rt0_type = typename RaviartThomas<0>::template apply<2, 2, double, Simplex<2>>::type;
    using rt1_type = typename RaviartThomas<1>::template apply<2, 2, double, Simplex<2>>::type;
    using rt_dynamic_type = typename RaviartThomas<Dynamic>::template apply<2, 2, double, Simplex<2>>::type;
    using nedelec0_type = typename Nedelec<0, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;
    using bdm0_type = typename BrezziDouglasMarini<0>::template apply<2, 2, double, Simplex<2>>::type;
    using bdm_dynamic_type = typename BrezziDouglasMarini<Dynamic>::template apply<2, 2, double, Simplex<2>>::type;
    using cr1_type = typename CrouzeixRaviart<1>::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( FiniteElementDofLayoutProvider<lagrange_p2_type> );
    static_assert( FiniteElementDofLayoutProvider<lagrange_dynamic_type> );
    static_assert( FiniteElementDofLayoutProvider<rt0_type> );
    static_assert( FiniteElementDofLayoutProvider<rt1_type> );
    static_assert( FiniteElementDofLayoutProvider<rt_dynamic_type> );
    static_assert( FiniteElementDofLayoutProvider<nedelec0_type> );
    static_assert( FiniteElementDofLayoutProvider<bdm0_type> );
    static_assert( FiniteElementDofLayoutProvider<bdm_dynamic_type> );
    static_assert( FiniteElementDofLayoutProvider<cr1_type> );
    static_assert( FiniteElementRuntimeSizedDofProvider<lagrange_p2_type> );
    static_assert( FiniteElementRuntimeSizedDofProvider<lagrange_dynamic_type> );
    static_assert( FiniteElementRuntimeSizedDofProvider<rt0_type> );
    static_assert( FiniteElementRuntimeSizedDofProvider<rt1_type> );
    static_assert( FiniteElementRuntimeSizedDofProvider<rt_dynamic_type> );
    static_assert( FiniteElementRuntimeSizedDofProvider<nedelec0_type> );
    static_assert( FiniteElementRuntimeSizedDofProvider<bdm0_type> );
    static_assert( FiniteElementRuntimeSizedDofProvider<bdm_dynamic_type> );
    static_assert( FiniteElementRuntimeSizedDofProvider<cr1_type> );
    static_assert( FiniteElementDofCardinalityProvider<lagrange_p2_type> );
    static_assert( FiniteElementDofCardinalityProvider<lagrange_dynamic_type> );
    static_assert( FiniteElementDofCardinalityProvider<rt0_type> );
    static_assert( FiniteElementDofCardinalityProvider<rt_dynamic_type> );
    static_assert( FiniteElementDofCardinalityProvider<nedelec0_type> );
    static_assert( FiniteElementDofCardinalityProvider<bdm0_type> );
    static_assert( FiniteElementDofCardinalityProvider<bdm_dynamic_type> );
    static_assert( FiniteElementDofCardinalityProvider<cr1_type> );
    static_assert( FiniteElementDofRepresentativePointProvider<lagrange_p2_type> );
    static_assert( FiniteElementDofRepresentativePointProvider<rt0_type> );
    static_assert( FiniteElementDofFunctionalProvider<lagrange_p2_type> );
    static_assert( FiniteElementDofFunctionalProvider<rt0_type> );
    static_assert( !FiniteElementDofTransformProvider<lagrange_p2_type, contract_mesh_type::element_type> );
    static_assert( FiniteElementDofTransformProvider<rt0_type, contract_mesh_type::element_type> );
    static_assert( FiniteElementDofTransformProvider<rt1_type, contract_mesh_type::element_type> );
    static_assert( FiniteElementDofTransformProvider<rt_dynamic_type, contract_mesh_type::element_type> );
    static_assert( FiniteElementDofTransformProvider<nedelec0_type, contract_mesh_type::element_type> );
    static_assert( FiniteElementDofTransformProvider<bdm0_type, contract_mesh_type::element_type> );
    static_assert( FiniteElementDofTransformProvider<bdm_dynamic_type, contract_mesh_type::element_type> );
    static_assert( !FiniteElementDofKeyProvider<lagrange_p2_type, contract_mesh_type::element_type> );

    lagrange_p2_type lagrangeP2;
    lagrange_dynamic_type lagrangeDynamicP2{ RuntimeOrder{ 2 } };
    rt0_type rt0;
    rt1_type rt1;
    rt_dynamic_type rtDynamic0{ RuntimeOrder{ 0 } };
    rt_dynamic_type rtDynamic1{ RuntimeOrder{ 1 } };
    nedelec0_type nedelec0;
    bdm0_type bdm0;
    bdm_dynamic_type bdmDynamic1{ RuntimeOrder{ 1 } };
    cr1_type cr1;

    checkLocalDescriptorContract( lagrangeP2, "Lagrange<2>" );
    checkLocalDescriptorContract( lagrangeDynamicP2, "Lagrange<Dynamic,P2>" );
    checkLocalDescriptorContract( rt0, "RaviartThomas<0>" );
    checkLocalDescriptorContract( rt1, "RaviartThomas<1>" );
    checkLocalDescriptorContract( rtDynamic0, "RaviartThomas<Dynamic,P0>" );
    checkLocalDescriptorContract( rtDynamic1, "RaviartThomas<Dynamic,P1>" );
    checkLocalDescriptorContract( nedelec0, "Nedelec<0>" );
    checkLocalDescriptorContract( bdm0, "BrezziDouglasMarini<0>" );
    checkLocalDescriptorContract( bdmDynamic1, "BrezziDouglasMarini<Dynamic,P1>" );
    checkLocalDescriptorContract( cr1, "CrouzeixRaviart<1>" );

    DofKey<> key{
        .topologicalDim = 1,
        .canonicalEntityId = 42,
        .ordinal = 0,
        .component = 0,
        .familyTag = 1,
        .variant = 0
    };
    auto nextKey = key;
    nextKey.ordinal = 1;
    BOOST_CHECK( key.isValid() );
    BOOST_CHECK( key < nextKey );

    DofTransform identity;
    DofTransform positiveSign{ .kind = DofTransformKind::Sign, .sign = 1 };
    DofTransform negativeSign{ .kind = DofTransformKind::Sign, .sign = -1 };
    BOOST_CHECK( identity.isIdentity() );
    BOOST_CHECK( positiveSign.isIdentity() );
    BOOST_CHECK( !negativeSign.isIdentity() );
}

BOOST_AUTO_TEST_CASE( phase2_descriptor_backed_entity_views )
{
    auto mesh = makeContractMesh();

    checkDescriptorBackedEntityViews( Pch<2>( mesh ), "Pch<2>" );
    checkDescriptorBackedEntityViews( Pch<Dynamic>( mesh, RuntimeOrder{ 2 } ), "Pch<Dynamic,P2>" );
    checkDescriptorBackedEntityViews( CRh( mesh ), "CrouzeixRaviart<1>" );
    checkDescriptorBackedEntityViews( RTh<1>( mesh ), "RTh<1>" );
    checkDescriptorBackedEntityViews( BDMh<0>( mesh ), "BDMh<0>" );
    checkDescriptorBackedEntityViews( Neh<0>( mesh ), "Neh<0>" );
}

BOOST_AUTO_TEST_CASE( phase3_descriptor_backed_local_to_global_keys )
{
    auto mesh = makeContractMesh();

    checkDescriptorKeyLocalToGlobalConsistency( Pch<2>( mesh ), "Pch<2>" );
    checkDescriptorKeyLocalToGlobalConsistency( Pch<Dynamic>( mesh, RuntimeOrder{ 2 } ), "Pch<Dynamic,P2>" );
    checkDescriptorKeyLocalToGlobalConsistency( CRh( mesh ), "CrouzeixRaviart<1>" );
    checkDescriptorKeyLocalToGlobalConsistency( RTh<1>( mesh ), "RTh<1>" );
    checkDescriptorKeyLocalToGlobalConsistency( BDMh<0>( mesh ), "BDMh<0>" );
    checkDescriptorKeyLocalToGlobalConsistency( Neh<0>( mesh ), "Neh<0>" );

    if ( Environment::worldComm().localSize() == 1 )
    {
        checkDescriptorKeyLocalToGlobalConsistency( Pch<3>( mesh ), "Pch<3>" );
        checkDescriptorKeyLocalToGlobalConsistency( Pch<Dynamic>( mesh, RuntimeOrder{ 3 } ), "Pch<Dynamic,P3>" );
    }
}

BOOST_AUTO_TEST_CASE( phase4_key_based_mpi_ownership )
{
    auto mesh = makeContractMesh();

    checkDescriptorKeyMpiOwnership( Pch<2>( mesh ), "Pch<2>" );
    checkDescriptorKeyMpiOwnership( Pch<Dynamic>( mesh, RuntimeOrder{ 2 } ), "Pch<Dynamic,P2>" );
    checkDescriptorKeyMpiOwnership( CRh( mesh ), "CrouzeixRaviart<1>" );
    checkDescriptorKeyMpiOwnership( RTh<1>( mesh ), "RTh<1>" );
    checkDescriptorKeyMpiOwnership( BDMh<0>( mesh ), "BDMh<0>" );
    checkDescriptorKeyMpiOwnership( Neh<0>( mesh ), "Neh<0>" );
}

BOOST_AUTO_TEST_CASE( phase5_hdiv_hcurl_transform_projection )
{
    auto mesh = makeContractMesh();

    checkTransformProjectionContract( RTh<0>( mesh ), "RTh<0>" );
    checkTransformProjectionContract( RTh<1>( mesh ), "RTh<1>" );
    checkTransformProjectionContract( RTh<Dynamic>( mesh, RuntimeOrder{ 0 } ), "RTh<Dynamic,P0>" );
    checkTransformProjectionContract( BDMh<0>( mesh ), "BDMh<0>" );
    checkTransformProjectionContract( Neh<0>( mesh ), "Neh<0>" );
}

BOOST_AUTO_TEST_CASE( phase5b_fe_owned_local_cardinality )
{
    auto mesh = makeContractMesh();

    checkFeOwnedLocalCardinality( Pch<2>( mesh ), "Pch<2>" );
    checkFeOwnedLocalCardinality( Pch<Dynamic>( mesh, RuntimeOrder{ 2 } ), "Pch<Dynamic,P2>" );
    checkFeOwnedLocalCardinality( CRh( mesh ), "CrouzeixRaviart<1>" );
    checkFeOwnedLocalCardinality( RTh<1>( mesh ), "RTh<1>" );
    checkFeOwnedLocalCardinality( BDMh<0>( mesh ), "BDMh<0>" );
    checkFeOwnedLocalCardinality( Neh<0>( mesh ), "Neh<0>" );
}

BOOST_AUTO_TEST_CASE( phase7a_runtime_sized_fe_instance_contract )
{
    auto mesh = makeContractMesh();

    auto lagrangeDynamicP2 = Pch<Dynamic>( mesh, RuntimeOrder{ 2 } );
    auto const& lagrangeFe = *lagrangeDynamicP2->basis();
    BOOST_CHECK_EQUAL( lagrangeFe.order(), 2 );
    BOOST_CHECK_EQUAL( lagrangeFe.runtimeOrder(), 2 );
    BOOST_CHECK_EQUAL( lagrangeDynamicP2->dof()->nLocalDof(), lagrangeFe.localDofCount() );
    BOOST_CHECK_EQUAL( lagrangeDynamicP2->dof()->nLocalDof( true ), lagrangeFe.localDofCount( true ) );
    BOOST_CHECK_GT( lagrangeFe.localDofCountOnEntity( 1, 0, true ), 0 );
    BOOST_REQUIRE( static_cast<bool>( validateFiniteElementDofLayout( lagrangeFe ) ) );

    auto rtDynamic0 = RTh<Dynamic>( mesh, RuntimeOrder{ 0 } );
    auto const& rtFe = *rtDynamic0->basis();
    BOOST_CHECK_EQUAL( rtFe.order(), 0 );
    BOOST_CHECK_EQUAL( rtFe.runtimeOrder(), 0 );
    BOOST_CHECK_EQUAL( rtDynamic0->dof()->nLocalDof(), rtFe.localDofCount() );
    BOOST_CHECK_EQUAL( rtDynamic0->dof()->nLocalDof( true ), rtFe.localDofCount( true ) );
    BOOST_CHECK_GT( rtFe.localDofCountOnFacet( 0, true ), 0 );
    BOOST_CHECK( !rtFe.dofHasRepresentativePoint( 0 ) );
    BOOST_REQUIRE( static_cast<bool>( validateFiniteElementDofLayout( rtFe ) ) );

    auto rtDynamic1 = RTh<Dynamic>( mesh, RuntimeOrder{ 1 } );
    auto const& rtFe1 = *rtDynamic1->basis();
    BOOST_CHECK_EQUAL( rtFe1.order(), 1 );
    BOOST_CHECK_EQUAL( rtFe1.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( rtFe1.localDofCount(), 8 );
    BOOST_CHECK_EQUAL( rtDynamic1->dof()->nLocalDof(), rtFe1.localDofCount() );
    BOOST_CHECK_GT( rtFe1.localDofCountOnFacet( 0, true ), 0 );
    BOOST_REQUIRE( static_cast<bool>( validateFiniteElementDofLayout( rtFe1 ) ) );

    auto bdmDynamic1 = BDMh<Dynamic>( mesh, RuntimeOrder{ 1 } );
    auto const& bdmFe1 = *bdmDynamic1->basis();
    BOOST_CHECK_EQUAL( bdmFe1.order(), 1 );
    BOOST_CHECK_EQUAL( bdmFe1.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( bdmFe1.localDofCount(), 12 );
    BOOST_CHECK_EQUAL( bdmDynamic1->dof()->nLocalDof(), bdmFe1.localDofCount() );
    BOOST_CHECK_GT( bdmFe1.localDofCountOnFacet( 0, true ), 0 );
    BOOST_REQUIRE( static_cast<bool>( validateFiniteElementDofLayout( bdmFe1 ) ) );
}

BOOST_AUTO_TEST_CASE( phase5c_descriptor_driven_flat_insertion )
{
    auto mesh = makeContractMesh();

    checkFlatDescriptorInsertion( Pch<2>( mesh ), "Pch<2>" );
    checkFlatDescriptorInsertion( Pch<Dynamic>( mesh, RuntimeOrder{ 2 } ), "Pch<Dynamic,P2>" );
    checkFlatDescriptorInsertion( Pchv<1>( mesh ), "Pchv<1>" );
    checkFlatDescriptorInsertion( CRh( mesh ), "CrouzeixRaviart<1>" );
    checkFlatDescriptorInsertion( RTh<1>( mesh ), "RTh<1>" );
    checkFlatDescriptorInsertion( BDMh<0>( mesh ), "BDMh<0>" );
    checkFlatDescriptorInsertion( Neh<0>( mesh ), "Neh<0>" );
}

BOOST_AUTO_TEST_CASE( phase6_dof_representative_points )
{
    auto mesh = makeContractMesh();

    checkRepresentativePointContract( Pch<2>( mesh ), "Pch<2>", true );
    checkRepresentativePointContract( Pch<Dynamic>( mesh, RuntimeOrder{ 2 } ), "Pch<Dynamic,P2>", true );
    checkRepresentativePointContract( Pchv<1>( mesh ), "Pchv<1>", true );
    checkRepresentativePointContract( CRh( mesh ), "CrouzeixRaviart<1>", true );
    checkRepresentativePointContract( RTh<0>( mesh ), "RTh<0>", false, DofFunctionalKind::NormalMoment );
    checkRepresentativePointContract( RTh<1>( mesh ), "RTh<1>", false, DofFunctionalKind::NormalMoment );
    checkRepresentativePointContract( RTh<Dynamic>( mesh, RuntimeOrder{ 0 } ), "RTh<Dynamic,P0>", false, DofFunctionalKind::NormalMoment );
    checkRepresentativePointContract( BDMh<0>( mesh ), "BDMh<0>", false, DofFunctionalKind::NormalMoment );
    checkRepresentativePointContract( Neh<0>( mesh ), "Neh<0>", false, DofFunctionalKind::TangentialMoment );
}

BOOST_AUTO_TEST_SUITE_END()
