/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#pragma once

#include "natural_second_operator_checks.hpp"
#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/ned2h.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>

namespace Feel::differential_interpolation_checks
{
using namespace natural_second_checks;

template <typename Operator, typename Field, typename ImageSpace, typename Exact>
void checkMatrix( Operator&& interpolation, Field const& field,
                  ImageSpace const& image, Exact const& exact )
{
    BOOST_TEST_CONTEXT( "operand=" << interpolation.interpolationType().interpolationOperand()
                                   << " domain=" << field.functionSpace()->basisName() << field.functionSpace()->order()
                                   << " image=" << image->basisName() << image->order() )
    {
        auto actual = image->element(), expected = image->element();
        interpolation.apply( field, actual );
        expected.on( _range = elements( image->mesh() ), _expr = exact,
                     _geomap = GeomapStrategyType::GEOMAP_HO );
        BOOST_CHECK_SMALL( normL2( _range = elements( image->mesh() ),
                                   _expr = idv( actual ) - idv( expected ), _quad = 12,
                                   _geomap = GeomapStrategyType::GEOMAP_HO ),
                           2e-8 );
    }
}

template <typename SpacePtr, typename MeshPtr>
void checkSpace( SpacePtr const& space, MeshPtr const& imageMesh, bool polynomial )
{
    using element_type = typename SpacePtr::element_type::element_type;
    constexpr int Dim = element_type::nDim;
    auto field = space->element();
    field.on( _range = elements( space->mesh() ), _expr = quadraticField<Dim>(),
              _geomap = GeomapStrategyType::GEOMAP_HO );
    auto scalarImage = Pdh<Dynamic>( imageMesh, RuntimeOrder{ 1 } );
    auto vectorImage = Pdhv<Dynamic>( imageMesh, RuntimeOrder{ 1 } );
    auto constantImage = Pdhv<Dynamic>( imageMesh, RuntimeOrder{ 0 } );
    BOOST_TEST_CONTEXT( space->basisName() << " FE" << space->order() << " G" << space->mesh()->order()
                                           << " related=" << imageMesh->isRelatedTo( space->mesh() ) )
    {
        auto divergence = Div( _domainSpace = space, _imageSpace = scalarImage );
        if ( polynomial )
            checkMatrix( divergence, field, scalarImage, exactDivergence<Dim>() );
        else
            checkMatrix( divergence, field, scalarImage, divv( field ) );

        auto checkCurl = [&]( auto const& image )
        {
            // An explicit policy must override the default through the shared factory.
            auto curl = Curl( _domainSpace = space, _imageSpace = image,
                              _type = makeCurlInterpolation( conforming_t{} ) );
            if ( polynomial )
                checkMatrix( curl, field, image, exactCurl<Dim>() );
            else
                checkMatrix( curl, field, image, curlv( field ) );
        };
        if constexpr ( Dim == 2 )
            checkCurl( scalarImage );
        else
            checkCurl( vectorImage );

        {
            auto gradDiv = GradDiv( _domainSpace = space, _imageSpace = constantImage );
            if ( polynomial )
                checkMatrix( gradDiv, field, constantImage, exactGradDiv<Dim>() );
            else
                checkMatrix( gradDiv, field, constantImage, graddivv( field ) );
        }
        {
            auto curlCurl = CurlCurl( _domainSpace = space, _imageSpace = constantImage );
            if ( polynomial )
                checkMatrix( curlCurl, field, constantImage, exactCurlCurl<Dim>() );
            else
                checkMatrix( curlCurl, field, constantImage, curlcurlv( field ) );
        }
    }
}

template <int Dim, typename MeshPtr>
void checkAllFamilies( MeshPtr const& mesh, MeshPtr const& imageMesh, bool polynomial )
{
    checkSpace( Pchv<Dynamic>( mesh, RuntimeOrder{ 2 } ), imageMesh, polynomial );
    checkSpace( RTh<Dynamic>( mesh, RuntimeOrder{ 2 } ), imageMesh, polynomial );
    checkSpace( BDMh<Dynamic>( mesh, RuntimeOrder{ 2 } ), imageMesh, polynomial );
    checkSpace( Ned1h<Dynamic>( mesh, RuntimeOrder{ 2 } ), imageMesh, polynomial );
    checkSpace( Ned2h<Dynamic>( mesh, RuntimeOrder{ 2 } ), imageMesh, polynomial );
}

template <int Dim>
void checkDimension()
{
    for ( uint16_type order : { 1, 2, 3 } )
    {
        auto mesh = box<Dim>( order );
        checkAllFamilies<Dim>( mesh, mesh, true );
    }
    for ( uint16_type order : { 2, 3 } )
    {
        auto mesh = curved<Dim>( order );
        checkAllFamilies<Dim>( mesh, mesh, false );
    }
    auto mesh = box<Dim>( 1 );
    checkSpace( RTh<2>( mesh ), mesh, true );
    checkSpace( Ned1h<2>( mesh ), mesh, true );

    auto hdiv = RTh<Dynamic>( mesh, RuntimeOrder{ 2 } );
    auto hcurl = Ned1h<Dynamic>( mesh, RuntimeOrder{ 2 } );
    auto divField = hdiv->element();
    auto curlField = hcurl->element();
    divField.on( _range = elements( mesh ), _expr = quadraticField<Dim>() );
    curlField.on( _range = elements( mesh ), _expr = quadraticField<Dim>() );
    auto divImage = RTh<Dynamic>( mesh, RuntimeOrder{ 0 } );
    auto curlImage = Ned1h<Dynamic>( mesh, RuntimeOrder{ 0 } );
    checkMatrix( I( _domainSpace = hdiv, _imageSpace = divImage ), divField, divImage, idv( divField ) );
    checkMatrix( I( _domainSpace = hcurl, _imageSpace = curlImage ), curlField, curlImage, idv( curlField ) );
    checkMatrix( GradDiv( _domainSpace = hdiv, _imageSpace = divImage ), divField, divImage, exactGradDiv<Dim>() );
    checkMatrix( CurlCurl( _domainSpace = hcurl, _imageSpace = curlImage ), curlField, curlImage, exactCurlCurl<Dim>() );

    auto potentialSpace = Pch<Dynamic>( mesh, RuntimeOrder{ 2 } );
    auto potential = potentialSpace->element();
    if constexpr ( Dim == 2 )
    {
        potential.on( _range = elements( mesh ), _expr = Px() * Px() + Px() * Py() + Py() * Py() );
        checkMatrix( Grad( _domainSpace = potentialSpace, _imageSpace = hcurl ), potential, hcurl,
                     vec( 2 * Px() + Py(), Px() + 2 * Py() ) );
    }
    else
    {
        potential.on( _range = elements( mesh ), _expr = Px() * Px() + Px() * Py() + Py() * Py() + Pz() * Pz() );
        checkMatrix( Grad( _domainSpace = potentialSpace, _imageSpace = hcurl ), potential, hcurl,
                     vec( 2 * Px() + Py(), Px() + 2 * Py(), 2 * Pz() ) );
        // The mapped de Rham composition must also vanish on genuinely curved cells.
        for ( uint16_type order : { 1, 2, 3 } )
        {
            auto chainMesh = order == 1 ? mesh : curved<Dim>( order );
            auto source = Ned1h<Dynamic>( chainMesh, RuntimeOrder{ 2 } );
            auto intermediate = RTh<Dynamic>( chainMesh, RuntimeOrder{ 2 } );
            auto image = Pdh<Dynamic>( chainMesh, RuntimeOrder{ 1 } );
            auto u = source->element();
            auto curlU = intermediate->element();
            auto divCurlU = image->element();
            u.on( _range = elements( chainMesh ), _expr = quadraticField<Dim>(), _geomap = GeomapStrategyType::GEOMAP_HO );
            Curl( _domainSpace = source, _imageSpace = intermediate ).apply( u, curlU );
            Div( _domainSpace = intermediate, _imageSpace = image ).apply( curlU, divCurlU );
            BOOST_CHECK_SMALL( normL2( _range = elements( chainMesh ), _expr = idv( divCurlU ),
                                       _quad = 12, _geomap = GeomapStrategyType::GEOMAP_HO ),
                               2e-7 );
        }
    }
    auto unrelatedMesh = box<Dim>( 2 );
    BOOST_REQUIRE( !unrelatedMesh->isRelatedTo( mesh ) );
    if ( Environment::worldComm().size() == 1 )
    {
        // A second mesh object exercises physical localization rather than related-cell pairing.
        checkAllFamilies<Dim>( mesh, unrelatedMesh, true );
        auto divImage = RTh<Dynamic>( unrelatedMesh, RuntimeOrder{ 0 } );
        auto curlImage = Ned1h<Dynamic>( unrelatedMesh, RuntimeOrder{ 0 } );
        checkMatrix( I( _domainSpace = hdiv, _imageSpace = divImage ), divField, divImage, idv( divField ) );
        checkMatrix( I( _domainSpace = hcurl, _imageSpace = curlImage ), curlField, curlImage, idv( curlField ) );
        checkMatrix( GradDiv( _domainSpace = hdiv, _imageSpace = divImage ), divField, divImage, exactGradDiv<Dim>() );
        checkMatrix( CurlCurl( _domainSpace = hcurl, _imageSpace = curlImage ), curlField, curlImage, exactCurlCurl<Dim>() );
    }
    else
    {
        auto image = Pdh<Dynamic>( unrelatedMesh, RuntimeOrder{ 1 } );
        BOOST_CHECK_THROW( ( Div( _domainSpace = hdiv, _imageSpace = image ) ), std::logic_error );
    }
}
} // namespace Feel::differential_interpolation_checks
