/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#pragma once

#include <feel/feeldiscr/stencil.hpp>

#include <type_traits>
#include <vector>

namespace Feel
{
namespace detail
{
template<typename SpaceT>
inline constexpr bool is_product_space_entry_v =
    std::is_base_of_v<ProductSpaceBase, decay_type<SpaceT>>;

template<typename SpaceT>
int
csrGraphBlockCount( SpaceT const& space )
{
    if constexpr ( is_product_space_entry_v<SpaceT> )
        return space->numberOfSpaces();
    else
        return 1;
}

template<typename PatternVectorT>
uint32_type
csrGraphBlockPattern( PatternVectorT const& patterns, int row, int col, int nCol )
{
    return uint32_type( patterns[row*nCol + col] );
}

template<typename TestSpaceT, typename TrialSpaceT, typename PatternVectorT, typename RangeMapT>
void
setCsrGraphBlock( BlocksBaseGraphCSR& graph,
                  TestSpaceT const& testSpace,
                  TrialSpaceT const& trialSpace,
                  int row, int col,
                  int nTrialBlocks,
                  PatternVectorT const& patterns,
                  RangeMapT const& range )
{
    LOG(INFO) << "filling out stencil (" << row << "," << col << ")\n";
    graph( row, col ) =
        stencil( _test=testSpace,
                 _trial=trialSpace,
                 _pattern=csrGraphBlockPattern( patterns, row, col, nTrialBlocks ),
                 _range=range,
                 _diag_is_nonzero=false,
                 _close=false )->graph();
}

template<typename TestSpaceT, typename TrialSpaceT, typename PatternVectorT, typename RangeMapT>
void
setCsrGraphBlocks( BlocksBaseGraphCSR& graph,
                   TestSpaceT const& testSpace,
                   TrialSpaceT const& trialSpace,
                   int rowStart, int colStart,
                   int nTrialBlocks,
                   PatternVectorT const& patterns,
                   RangeMapT const& range )
{
    if constexpr ( is_product_space_entry_v<TestSpaceT> && is_product_space_entry_v<TrialSpaceT> )
    {
        for ( int i = 0; i < testSpace->numberOfSpaces(); ++i )
            for ( int j = 0; j < trialSpace->numberOfSpaces(); ++j )
                setCsrGraphBlock( graph, (*testSpace)[i], (*trialSpace)[j],
                                  rowStart + i, colStart + j, nTrialBlocks, patterns, range );
    }
    else if constexpr ( is_product_space_entry_v<TestSpaceT> )
    {
        for ( int i = 0; i < testSpace->numberOfSpaces(); ++i )
            setCsrGraphBlock( graph, (*testSpace)[i], trialSpace,
                              rowStart + i, colStart, nTrialBlocks, patterns, range );
    }
    else if constexpr ( is_product_space_entry_v<TrialSpaceT> )
    {
        for ( int j = 0; j < trialSpace->numberOfSpaces(); ++j )
            setCsrGraphBlock( graph, testSpace, (*trialSpace)[j],
                              rowStart, colStart + j, nTrialBlocks, patterns, range );
    }
    else
    {
        setCsrGraphBlock( graph, testSpace, trialSpace,
                          rowStart, colStart, nTrialBlocks, patterns, range );
    }
}

template<typename PatternVectorT>
void
checkCsrGraphBlockPatterns( PatternVectorT const& patterns, int nRow, int nCol )
{
    CHECK_EQ( patterns.size(), static_cast<std::size_t>( nRow*nCol ) )
        << "invalid block pattern size for product space graph: got "
        << patterns.size() << ", expected " << nRow*nCol;
}
} // namespace detail

/**
 * Build rectangular blocks of CSR graphs from distinct row/test and column/trial
 * product spaces.
 */
template<typename TestPS, typename TrialPS, typename PatternSizeT, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpacesBase, std::remove_reference_t<TestPS>>
          && std::is_base_of_v<ProductSpacesBase, std::remove_reference_t<TrialPS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( TestPS&& testPs,
                TrialPS&& trialPs,
                std::vector<PatternSizeT> const& patterns,
                RangeMapT range = stencilRangeMap() )
{
    int nTestBlocks = testPs.numberOfSpaces();
    int nTrialBlocks = trialPs.numberOfSpaces();
    Feel::detail::checkCsrGraphBlockPatterns( patterns, nTestBlocks, nTrialBlocks );
    BlocksBaseGraphCSR graph( nTestBlocks, nTrialBlocks );

    int rowStart = 0;
    auto const& testTuple = testPs.tupleSpaces();
    auto const& trialTuple = trialPs.tupleSpaces();
    hana::for_each( testTuple, [&]( auto const& testSpace )
                    {
                        int colStart = 0;
                        hana::for_each( trialTuple, [&]( auto const& trialSpace )
                                        {
                                            Feel::detail::setCsrGraphBlocks( graph, testSpace, trialSpace,
                                                                             rowStart, colStart, nTrialBlocks, patterns, range );
                                            colStart += Feel::detail::csrGraphBlockCount( trialSpace );
                                        } );
                        rowStart += Feel::detail::csrGraphBlockCount( testSpace );
                    } );
    return graph;
}

template<typename TestPS, typename TrialPS, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpacesBase, std::remove_reference_t<TestPS>>
          && std::is_base_of_v<ProductSpacesBase, std::remove_reference_t<TrialPS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( TestPS&& testPs,
                TrialPS&& trialPs,
                uint32_type pattern = Pattern::COUPLED,
                RangeMapT range = stencilRangeMap() )
{
    return csrGraphBlocks( testPs, trialPs,
                           std::vector<size_type>( testPs.numberOfSpaces()*trialPs.numberOfSpaces(), pattern ),
                           range );
}

template<typename TestPS, typename TrialPS, typename PatternSizeT, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpacesBase, std::remove_reference_t<TestPS>>
          && std::is_base_of_v<ProductSpaceBase, std::remove_reference_t<TrialPS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( TestPS&& testPs,
                TrialPS&& trialPs,
                std::vector<PatternSizeT> const& patterns,
                RangeMapT range = stencilRangeMap() )
{
    int nTestBlocks = testPs.numberOfSpaces();
    int nTrialBlocks = trialPs.numberOfSpaces();
    Feel::detail::checkCsrGraphBlockPatterns( patterns, nTestBlocks, nTrialBlocks );
    BlocksBaseGraphCSR graph( nTestBlocks, nTrialBlocks );

    int rowStart = 0;
    auto const& testTuple = testPs.tupleSpaces();
    hana::for_each( testTuple, [&]( auto const& testSpace )
                    {
                        for ( int col = 0; col < nTrialBlocks; ++col )
                            Feel::detail::setCsrGraphBlocks( graph, testSpace, trialPs[col],
                                                             rowStart, col, nTrialBlocks, patterns, range );
                        rowStart += Feel::detail::csrGraphBlockCount( testSpace );
                    } );
    return graph;
}

template<typename TestPS, typename TrialPS, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpacesBase, std::remove_reference_t<TestPS>>
          && std::is_base_of_v<ProductSpaceBase, std::remove_reference_t<TrialPS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( TestPS&& testPs,
                TrialPS&& trialPs,
                uint32_type pattern = Pattern::COUPLED,
                RangeMapT range = stencilRangeMap() )
{
    return csrGraphBlocks( testPs, trialPs,
                           std::vector<size_type>( testPs.numberOfSpaces()*trialPs.numberOfSpaces(), pattern ),
                           range );
}

template<typename TestPS, typename TrialPS, typename PatternSizeT, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpaceBase, std::remove_reference_t<TestPS>>
          && std::is_base_of_v<ProductSpacesBase, std::remove_reference_t<TrialPS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( TestPS&& testPs,
                TrialPS&& trialPs,
                std::vector<PatternSizeT> const& patterns,
                RangeMapT range = stencilRangeMap() )
{
    int nTestBlocks = testPs.numberOfSpaces();
    int nTrialBlocks = trialPs.numberOfSpaces();
    Feel::detail::checkCsrGraphBlockPatterns( patterns, nTestBlocks, nTrialBlocks );
    BlocksBaseGraphCSR graph( nTestBlocks, nTrialBlocks );

    auto const& trialTuple = trialPs.tupleSpaces();
    for ( int row = 0; row < nTestBlocks; ++row )
    {
        int colStart = 0;
        hana::for_each( trialTuple, [&]( auto const& trialSpace )
                        {
                            Feel::detail::setCsrGraphBlocks( graph, testPs[row], trialSpace,
                                                             row, colStart, nTrialBlocks, patterns, range );
                            colStart += Feel::detail::csrGraphBlockCount( trialSpace );
                        } );
    }
    return graph;
}

template<typename TestPS, typename TrialPS, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpaceBase, std::remove_reference_t<TestPS>>
          && std::is_base_of_v<ProductSpacesBase, std::remove_reference_t<TrialPS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( TestPS&& testPs,
                TrialPS&& trialPs,
                uint32_type pattern = Pattern::COUPLED,
                RangeMapT range = stencilRangeMap() )
{
    return csrGraphBlocks( testPs, trialPs,
                           std::vector<size_type>( testPs.numberOfSpaces()*trialPs.numberOfSpaces(), pattern ),
                           range );
}

template<typename TestPS, typename TrialPS, typename PatternSizeT, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpaceBase, std::remove_reference_t<TestPS>>
          && std::is_base_of_v<ProductSpaceBase, std::remove_reference_t<TrialPS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( TestPS&& testPs,
                TrialPS&& trialPs,
                std::vector<PatternSizeT> const& patterns,
                RangeMapT range = stencilRangeMap() )
{
    int nTestBlocks = testPs.numberOfSpaces();
    int nTrialBlocks = trialPs.numberOfSpaces();
    Feel::detail::checkCsrGraphBlockPatterns( patterns, nTestBlocks, nTrialBlocks );
    BlocksBaseGraphCSR graph( nTestBlocks, nTrialBlocks );

    for ( int row = 0; row < nTestBlocks; ++row )
        for ( int col = 0; col < nTrialBlocks; ++col )
            Feel::detail::setCsrGraphBlock( graph, testPs[row], trialPs[col],
                                            row, col, nTrialBlocks, patterns, range );
    return graph;
}

template<typename TestPS, typename TrialPS, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpaceBase, std::remove_reference_t<TestPS>>
          && std::is_base_of_v<ProductSpaceBase, std::remove_reference_t<TrialPS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( TestPS&& testPs,
                TrialPS&& trialPs,
                uint32_type pattern = Pattern::COUPLED,
                RangeMapT range = stencilRangeMap() )
{
    return csrGraphBlocks( testPs, trialPs,
                           std::vector<size_type>( testPs.numberOfSpaces()*trialPs.numberOfSpaces(), pattern ),
                           range );
}

template<typename PS, typename PatternSizeT, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpacesBase, decay_type<PS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( PS&& ps,
                std::vector<PatternSizeT> const& patterns,
                RangeMapT range = stencilRangeMap() )
{
    auto&& productSpace = remove_shared_ptr_f( std::forward<PS>( ps ) );
    return csrGraphBlocks( productSpace, productSpace, patterns, range );
}

template<typename PS, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpacesBase, decay_type<PS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( PS&& ps,
                uint32_type pattern = Pattern::COUPLED,
                RangeMapT range = stencilRangeMap() )
{
    auto&& productSpace = remove_shared_ptr_f( std::forward<PS>( ps ) );
    return csrGraphBlocks( productSpace, productSpace, pattern, range );
}

template<typename PS, typename PatternSizeT, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpaceBase, decay_type<PS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( PS&& ps,
                std::vector<PatternSizeT> const& patterns,
                RangeMapT range = stencilRangeMap() )
{
    auto&& productSpace = remove_shared_ptr_f( std::forward<PS>( ps ) );
    return csrGraphBlocks( productSpace, productSpace, patterns, range );
}

template<typename PS, typename RangeMapT = StencilRangeMap0Type>
    requires std::is_base_of_v<ProductSpaceBase, decay_type<PS>>
          && std::is_base_of_v<StencilRangeMapTypeBase, RangeMapT>
BlocksBaseGraphCSR
csrGraphBlocks( PS&& ps,
                uint32_type pattern = Pattern::COUPLED,
                RangeMapT range = stencilRangeMap() )
{
    auto&& productSpace = remove_shared_ptr_f( std::forward<PS>( ps ) );
    return csrGraphBlocks( productSpace, productSpace, pattern, range );
}
}
