/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-12-23

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
#ifndef __galerkingraph_H
#define __galerkingraph_H 1

#define FEELPP_EXPORT_GRAPH 0
#if FEELPP_EXPORT_GRAPH
#include <feel/feelfilters/exporter.hpp>
#endif

#include <feel/feelvf/pattern.hpp>
#include <feel/feelvf/block.hpp>
#include <feel/feelalg/graphcsr.hpp>
#include <feel/feeldiscr/functionspace.hpp>


#if 1
namespace Feel
{
namespace detail
{
template<typename BFType, typename Space1Type>
struct compute_graph3
{
    compute_graph3( BFType* bf, boost::shared_ptr<Space1Type> const& space1, size_type trial_index, size_type hints  )
        :
        M_stencil( bf ),
        M_space1( space1 ),
        M_test_index( 0 ),
        M_trial_index( trial_index ),
        M_hints( hints )
    {}

    template <typename Space2>
    void operator()( boost::shared_ptr<Space2> const& space2 ) const
    {
        if ( M_stencil->testSpace()->worldsComm()[M_test_index].isActive() )
        {
            if ( M_stencil->isBlockPatternZero( M_test_index,M_trial_index ) )
            {
#if !defined(FEELPP_ENABLE_MPI_MODE)
                const size_type proc_id           = M_stencil->testSpace()->mesh()->comm().rank();
                const size_type n1_dof_on_proc    = space2->nLocalDof();
                const size_type first1_dof_on_proc = space2->dof()->firstDof( proc_id );
                const size_type last1_dof_on_proc = space2->dof()->lastDof( proc_id );
                const size_type first2_dof_on_proc = M_space1->dof()->firstDof( proc_id );
                const size_type last2_dof_on_proc = M_space1->dof()->lastDof( proc_id );
#else
                const size_type proc_id           = M_stencil->testSpace()->worldsComm()[M_test_index].globalRank();
                const size_type n1_dof_on_proc    = space2->nLocalDof();
                const size_type first1_dof_on_proc = space2->dof()->firstDofGlobalCluster( proc_id );
                const size_type last1_dof_on_proc = space2->dof()->lastDofGlobalCluster( proc_id );
                const size_type first2_dof_on_proc = M_space1->dof()->firstDofGlobalCluster( proc_id );
                const size_type last2_dof_on_proc = M_space1->dof()->lastDofGlobalCluster( proc_id );
#endif
                typename BFType::graph_ptrtype zerograph( new typename BFType::graph_type( n1_dof_on_proc,
                        first1_dof_on_proc, last1_dof_on_proc,
                        first2_dof_on_proc, last2_dof_on_proc,
                        space2->worldComm() ) );
                zerograph->zero();
                M_stencil->mergeGraph( M_stencil->testSpace()->nDofStart( M_test_index ), M_stencil->trialSpace()->nDofStart( M_trial_index ) , zerograph );
            }

            else
            {
                auto thestencil = stencil( _test=space2, _trial=M_space1,
                                           _pattern=M_stencil->blockPattern( M_test_index,M_trial_index ),
                                           _pattern_block=M_stencil->blockPattern(),
                                           _diag_is_nonzero=false,
                                           _collect_garbage=false );

                if ( M_stencil->testSpace()->worldComm().globalSize()>1 && M_stencil->testSpace()->hasEntriesForAllSpaces() )
                    M_stencil->mergeGraphMPI( M_test_index, M_trial_index,
                                              space2->mapOn(), M_space1->mapOn(),
                                              thestencil->graph() );
                else
                    M_stencil->mergeGraph( M_stencil->testSpace()->nDofStart( M_test_index ),
                                           M_stencil->trialSpace()->nDofStart( M_trial_index ),
                                           thestencil->graph() );
            }
        } // if ( M_stencil->testSpace()->worldsComm()[M_test_index].isActive() )

        ++M_test_index;
    }


    mutable BFType* M_stencil;
    boost::shared_ptr<Space1Type> const& M_space1;
    mutable size_type M_test_index;
    size_type M_trial_index;
    size_type M_hints;
};


template<typename BFType, typename Space1Type>
struct compute_graph2
{
    compute_graph2( BFType* bf, boost::shared_ptr<Space1Type> const& space1, size_type test_index, size_type hints  )
        :
        M_stencil( bf ),
        M_space1( space1 ),
        M_test_index( test_index ),
        M_trial_index( 0 ),
        M_hints( hints )
    {}

    template <typename Space2>
    void operator()( boost::shared_ptr<Space2> const& space2 ) const
    {
        if ( M_stencil->testSpace()->worldsComm()[M_test_index].isActive() )
        {
            if ( M_stencil->isBlockPatternZero( M_test_index,M_trial_index ) )
            {
#if !defined(FEELPP_ENABLE_MPI_MODE)
                const size_type proc_id           = M_stencil->testSpace()->template mesh<0>()->comm().rank();
                const size_type n1_dof_on_proc    = M_space1->nLocalDof();
                const size_type first1_dof_on_proc = M_space1->dof()->firstDof( proc_id );
                const size_type last1_dof_on_proc = M_space1->dof()->lastDof( proc_id );
                const size_type first2_dof_on_proc = space2->dof()->firstDof( proc_id );
                const size_type last2_dof_on_proc = space2->dof()->lastDof( proc_id );
#else
                const size_type proc_id           = M_stencil->testSpace()->worldsComm()[M_test_index].globalRank();
                const size_type n1_dof_on_proc    = M_space1->nLocalDof();
                const size_type first1_dof_on_proc = M_space1->dof()->firstDofGlobalCluster( proc_id );
                const size_type last1_dof_on_proc = M_space1->dof()->lastDofGlobalCluster( proc_id );
                const size_type first2_dof_on_proc = space2->dof()->firstDofGlobalCluster( proc_id );
                const size_type last2_dof_on_proc = space2->dof()->lastDofGlobalCluster( proc_id );
#endif

                typename BFType::graph_ptrtype zerograph( new typename BFType::graph_type( n1_dof_on_proc,
                        first1_dof_on_proc, last1_dof_on_proc,
                        first2_dof_on_proc, last2_dof_on_proc,
                        M_space1->worldComm() ) );
                zerograph->zero();
                M_stencil->mergeGraph( M_stencil->testSpace()->nDofStart( M_test_index ),
                                       M_stencil->trialSpace()->nDofStart( M_trial_index ),
                                       zerograph );
            }

            else
            {

                auto thestencil = stencil( _test=M_space1, _trial=space2,
                                           _pattern=M_stencil->blockPattern( M_test_index,M_trial_index ),
                                           _pattern_block=M_stencil->blockPattern(),
                                           _diag_is_nonzero=false,
                                           _collect_garbage=false );

                if ( M_stencil->testSpace()->worldComm().globalSize()>1 && M_stencil->testSpace()->hasEntriesForAllSpaces() )
                    M_stencil->mergeGraphMPI( M_test_index, M_trial_index,
                                              M_space1->mapOn(), space2->mapOn(),
                                              thestencil->graph() );
                else
                    M_stencil->mergeGraph( M_stencil->testSpace()->nDofStart( M_test_index ),
                                           M_stencil->trialSpace()->nDofStart( M_trial_index ),
                                           thestencil->graph() );
            }
        } // if ( M_stencil->testSpace()->worldsComm()[M_test_index].isActive() )

        ++M_trial_index;
    }


    mutable BFType* M_stencil;
    boost::shared_ptr<Space1Type> const& M_space1;
    size_type M_test_index;
    mutable size_type M_trial_index;
    size_type M_hints;
};


template<typename BFType>
struct compute_graph1
{
    compute_graph1( BFType* bf, size_type hints )
        :
        M_stencil( bf ),
        M_test_index( 0 ),
        M_hints( hints )
    {}

    template <typename Space1>
    void operator()( boost::shared_ptr<Space1> const& space1 ) const
    {
        fusion::for_each( M_stencil->trialSpace()->functionSpaces(),
                          compute_graph2<BFType,Space1>( M_stencil, space1, M_test_index, M_hints ) );
        ++M_test_index;
    }
    mutable BFType* M_stencil;
    mutable size_type M_test_index;
    size_type M_hints;
};

}
template<typename X1, typename X2>
class Stencil
{
public:
    typedef X1 test_space_ptrtype;
    typedef X2 trial_space_ptrtype;
    typedef typename X1::value_type test_space_type;
    typedef typename X2::value_type trial_space_type;
    typedef GraphCSR graph_type;
    typedef boost::shared_ptr<graph_type> graph_ptrtype;
    typedef Stencil<X1,X2> self_type;

    Stencil( test_space_ptrtype Xh, trial_space_ptrtype Yh,
             size_type graph_hints,
             BlocksStencilPattern block_pattern=BlocksStencilPattern(1,1,Pattern::HAS_NO_BLOCK_PATTERN),
             bool diag_is_nonzero=false )
        :
        _M_X1( Xh ),
        _M_X2( Yh ),
#if !defined(FEELPP_ENABLE_MPI_MODE)
        M_graph( new graph_type( Xh->nLocalDof(),
                                 Xh->nDofStart(), Xh->nDofStart()+ Xh->nLocalDof()-1,
                                 Yh->nDofStart(), Yh->nDofStart()+ Yh->nLocalDof()-1 ) ),
#else
        M_graph( new graph_type( Xh->nLocalDof(),
                                 Xh->dof()->firstDofGlobalCluster( Xh->worldComm().globalRank() ), Xh->dof()->lastDofGlobalCluster( Xh->worldComm().globalRank() ),
                                 Yh->dof()->firstDofGlobalCluster( Yh->worldComm().globalRank() ), Yh->dof()->lastDofGlobalCluster( Yh->worldComm().globalRank() ),
                                 Xh->worldComm() ) ),
#endif
        M_block_pattern( block_pattern )
    {
        // init block_pattern if empty
        uint16_type nbSubSpace1 = _M_X1->nSubFunctionSpace();
        uint16_type nbSubSpace2 = _M_X2->nSubFunctionSpace();

        if ( this->isBlockPatternNoPattern( 0,0 ) )
        {
            M_block_pattern = BlocksStencilPattern(nbSubSpace1,nbSubSpace2,graph_hints);
        }

        M_graph = this->computeGraph( graph_hints );

        if ( diag_is_nonzero && _M_X1->nLocalDofWithoutGhost()>0 && _M_X2->nLocalDofWithoutGhost()>0 ) M_graph->addMissingZeroEntriesDiagonal();

        M_graph->close();
    }

    Stencil( test_space_ptrtype Xh, trial_space_ptrtype Yh, size_type graph_hints, graph_ptrtype g )
        :
        _M_X1( Xh ),
        _M_X2( Yh ),
        M_graph( g ),
        M_block_pattern(Xh->nSubFunctionSpace(),Yh->nSubFunctionSpace(),size_type( graph_hints/*Pattern::HAS_NO_BLOCK_PATTERN*/ ))
    {}


    BlocksStencilPattern blockPattern() const
    {
        return  M_block_pattern;
    }

    size_type blockPattern( uint16_type i,uint16_type j ) const
    {
        return  M_block_pattern(i,j);
    }

    bool isBlockPatternNoPattern( uint16_type i,uint16_type j ) const
    {
        Feel::Context ctx( M_block_pattern(i,j) );
        return ctx.test( HAS_NO_BLOCK_PATTERN );

    }
    bool isBlockPatternZero( uint16_type i,uint16_type j ) const
    {
        Feel::Context ctx( M_block_pattern(i,j) );
        return ctx.test( ZERO );
    }

    graph_ptrtype computeGraph( size_type hints );

    graph_ptrtype computeGraph( size_type hints, mpl::bool_<true> );
    graph_ptrtype computeGraph( size_type hints, mpl::bool_<false> );
    graph_ptrtype computeGraph( size_type hints, mpl::bool_<true>, mpl::bool_<true> );
    graph_ptrtype computeGraph( size_type hints, mpl::bool_<false>, mpl::bool_<true> );
    graph_ptrtype computeGraph( size_type hints, mpl::bool_<true>, mpl::bool_<false> );

    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true> );
    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false> );
    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<true> );
    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false>, mpl::bool_<true> );
    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<false> );


    void mergeGraph( int row, int col, graph_ptrtype g );
    void mergeGraphMPI( size_type test_index, size_type trial_index,
                        DataMap const& mapOnTest, DataMap const& mapOnTrial,
                        graph_ptrtype g);

    test_space_ptrtype testSpace() const
    {
        return _M_X1;
    }
    trial_space_ptrtype trialSpace() const
    {
        return _M_X2;
    }
    graph_ptrtype graph() const
    {
        return M_graph;
    }
    graph_ptrtype graph()
    {
        return M_graph;
    }

private:

    test_space_ptrtype _M_X1;
    trial_space_ptrtype _M_X2;
    graph_ptrtype M_graph;
    BlocksStencilPattern M_block_pattern;

};
namespace detail
{
template<typename Args>
struct compute_stencil_type
{
    typedef typename remove_pointer_const_reference_type<Args,tag::test>::type _test_type;
    typedef typename remove_pointer_const_reference_type<Args,tag::trial>::type _trial_type;
    typedef Stencil<_test_type, _trial_type> type;
    typedef boost::shared_ptr<type> ptrtype;

};

}

class StencilManagerImpl:
    public std::map<boost::tuple<boost::shared_ptr<FunctionSpaceBase>,
    boost::shared_ptr<FunctionSpaceBase>,
    size_type,
    std::vector<size_type>,
    bool >, boost::shared_ptr<GraphCSR> >,
public boost::noncopyable
{
public:
    typedef boost::shared_ptr<GraphCSR> graph_ptrtype;
    typedef boost::tuple<boost::shared_ptr<FunctionSpaceBase>,
            boost::shared_ptr<FunctionSpaceBase>,
            size_type,
            std::vector<size_type>,
            bool > key_type;
    typedef std::map<key_type, graph_ptrtype> graph_manager_type;

};

typedef Feel::Singleton<StencilManagerImpl> StencilManager;

//! function that cleans up the StencilManager any time \c stencil() is called
void stencilManagerGarbageCollect();
//! function that cleans up one entry the StencilManager
void stencilManagerGarbage(StencilManagerImpl::key_type const& key);
//! function that add an entry in the StencilManager
void stencilManagerAdd(StencilManagerImpl::key_type const& key,StencilManagerImpl::graph_ptrtype graph);
//! print entries stored in stencil manager
void stencilManagerPrint();

extern BlocksStencilPattern default_block_pattern;

BOOST_PARAMETER_FUNCTION(
    ( typename detail::compute_stencil_type<Args>::ptrtype ), // 1. return type
    stencil,                                       // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required                                  // 4. one required parameter, and
      ( test,             *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
      ( trial,            *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
    )
    ( optional                                  //    four optional parameters, with defaults
      ( pattern,          *( boost::is_integral<mpl::_> ), Pattern::COUPLED )
      ( pattern_block,    *, default_block_pattern )
      ( diag_is_nonzero,  *( boost::is_integral<mpl::_> ), false )
      ( collect_garbage, *( boost::is_integral<mpl::_> ), true )
    )
)
{
    if ( collect_garbage )
    {
        // cleanup memory before doing anything
        stencilManagerGarbageCollect();
    }

    Feel::detail::ignore_unused_variable_warning( args );
    typedef typename detail::compute_stencil_type<Args>::ptrtype stencil_ptrtype;
    typedef typename detail::compute_stencil_type<Args>::type stencil_type;

    // we look into the spaces dictionary for existing graph
    auto git = StencilManager::instance().find( boost::make_tuple( test, trial, pattern, pattern_block.getSetOfBlocks(), diag_is_nonzero ) );

    if (  git != StencilManager::instance().end() )
    {
        //std::cout << "Found a  stencil in manager (" << test.get() << "," << trial.get() << "," << pattern << ")\n";
        auto s = stencil_ptrtype( new stencil_type( test, trial, pattern, git->second ) );
        return s;
    }

    else
    {
        // look for transposed stencil if it exist and transpose it to get the stencil
        auto git_trans = StencilManager::instance().find( boost::make_tuple( trial, test, pattern, pattern_block.transpose().getSetOfBlocks(), diag_is_nonzero ) );

        if ( git_trans != StencilManager::instance().end() )
        {
            auto g = git_trans->second->transpose(test->mapOn());
            //auto g = git_trans->second->transpose();
            StencilManager::instance().operator[]( boost::make_tuple( test, trial, pattern, pattern_block.getSetOfBlocks(), diag_is_nonzero ) ) = g;
            auto s = stencil_ptrtype( new stencil_type( test, trial, pattern, g ) );
            //std::cout << "Found a  transposed stencil in manager (" << test.get() << "," << trial.get() << "," << pattern << ")\n";
            return s;
        }

        else
        {
            //std::cout << "Creating a new stencil in manager (" << test.get() << "," << trial.get() << "," << pattern << ")\n";
            auto s = stencil_ptrtype( new stencil_type( test, trial, pattern, pattern_block, diag_is_nonzero ) );
            StencilManager::instance().operator[]( boost::make_tuple( test, trial, pattern, pattern_block.getSetOfBlocks(), diag_is_nonzero ) ) = s->graph();
            return s;
        }
    }
}


namespace detail
{
template<typename BidirectionalIterator>
inline
void
sortSparsityRow ( const BidirectionalIterator begin,
                  BidirectionalIterator       middle,
                  const BidirectionalIterator end )
{
    if ( ( begin == middle ) || ( middle == end ) ) return;

    assert ( std::distance ( begin,  middle ) > 0 );
    assert ( std::distance ( middle, end )    > 0 );
    FEELPP_ASSERT( std::unique ( begin,  middle ) == middle )
    ( *begin )( *middle ).error( "duplicate dof(begin,middle)" );
    FEELPP_ASSERT ( std::unique ( middle, end )    == end )
    ( *begin )( *middle ).error( "duplicate dof (middle,end)" );

    while ( middle != end )
    {
        BidirectionalIterator
        b = middle,
        a = b-1;

        // Bubble-sort the middle value downward
        while ( !( *a < *b ) ) // *a & *b are less-than comparable, so use <
        {
            std::swap ( *a, *b );

            if ( a == begin ) break;

            b=a;
            --a;
        }

        ++middle;
    }

    // Assure the algorithm worked if we are in DEBUG mode
#ifdef DEBUG
    {
        // SGI STL extension!
        // assert (std::is_sorted(begin,end));

        BidirectionalIterator
        prev  = begin,
        first = begin;

        for ( ++first; first != end; prev=first, ++first )
            if ( *first < *prev )
                assert( false );
    }
#endif

    // Make sure the two ranges did not contain any common elements
    assert ( std::unique ( begin, end ) == end );
} //

}
template<typename X1,  typename X2>
void
Stencil<X1,X2>::mergeGraph( int row, int col, graph_ptrtype g )
{
    boost::timer tim;
    Debug( 5050 ) << "[merge graph] for composite bilinear form\n";
    Debug( 5050 ) << "[mergeGraph] row = " << row << "\n";
    Debug( 5050 ) << "[mergeGraph] col = " << col << "\n";

    // nothing yet in store
    //if ( !M_graph || M_graph->empty() )
    if ( 0 )
    {
        Debug( 5050 ) << "[merge graph] nothing yet in store, copy graph\n";
        M_graph = g;

#if 0
        typename graph_type::const_iterator it = g->begin();
        typename graph_type::const_iterator en = g->end();

        for ( ; it != en; ++it )
        {
            std::vector<size_type> const& ivec = boost::get<2>( it->second );

            for ( int i = 0; i < ivec.size(); ++i )
            {
                Debug( 5050 ) << "[mergeGraph] ivec[" << i << "] = " << ivec[i] << "\n";
            }
        }

#endif // 9
    }

    else
    {
        //std::cout << "\n row " << row << " col " << col << " with god rank" <<  this->testSpace()->worldComm().godRank() << std::endl;
        Debug( 5050 ) << "[merge graph] already something in store\n";
        typename graph_type::const_iterator it = g->begin();
        typename graph_type::const_iterator en = g->end();

        for ( ; it != en; ++it )
        {
            int theglobalrow = row+it->first;
            int thelocalrow;// warning : not the same in parallel

            if ( this->testSpace()->worldComm().globalSize()>1 )
                thelocalrow = /*row +*/ boost::get<1>( it->second );

            else
                thelocalrow = row + boost::get<1>( it->second );

            //auto row1_entries = boost::unwrap_ref( boost::ref( M_graph->row(theglobalrow).template get<2>() ) );
            std::set<size_type>& row1_entries = M_graph->row( theglobalrow ).template get<2>();
            std::set<size_type> const& row2_entries = boost::get<2>( it->second );

            Debug( 5050 ) << "[mergeGraph] adding information to global row [" << theglobalrow << "], localrow=" << thelocalrow << "\n";
            M_graph->row( theglobalrow ).template get<1>() = thelocalrow;
            M_graph->row( theglobalrow ).template get<0>() = this->testSpace()->worldComm().mapLocalRankToGlobalRank()[it->second.get<0>()];

            if ( !row2_entries.empty() )
            {

                if ( row1_entries.empty() )
                {
                    // if row is empty then no need to shift the dof in
                    // composite case since the merge in done block-row-wise
                    //row1_entries = row2_entries;
                    if ( col==0 ) row1_entries = row2_entries;

                    else
                    {
                        for ( auto it = row2_entries.begin(), en = row2_entries.end() ; it!=en; ++it ) row1_entries.insert( *it+col );
                    }

                }

                else
                {
                    // ensure unique sorted ids
                    auto itg = boost::prior( row1_entries.end() );
                    // shift dofs in case of composite spaces
                    std::for_each( row2_entries.begin(), row2_entries.end(),[&]( size_type o )
                    {
                        itg = row1_entries.insert( itg, o+col );
                    } );
                }
            }
        } // for( ; it != en; ++it )
    }

    Debug( 5050 ) << " -- merge_graph (" << row << "," << col << ") in " << tim.elapsed() << "\n";
    Debug( 5050 ) << "merge graph for composite bilinear form done\n";
}

template<typename X1,  typename X2>
void
Stencil<X1,X2>::mergeGraphMPI( size_type test_index, size_type trial_index,
                               DataMap const& mapOnTest, DataMap const& mapOnTrial,
                               graph_ptrtype g )
{

    const int row = ( this->testSpace()->dof()->firstDofGlobalCluster()  +  this->testSpace()->nLocalDofWithoutGhostStart( test_index ) ) - mapOnTest.firstDofGlobalCluster();
    const int col = ( this->trialSpace()->dof()->firstDofGlobalCluster() +  this->trialSpace()->nLocalDofWithoutGhostStart( trial_index ) ) - mapOnTrial.firstDofGlobalCluster();
    typename graph_type::const_iterator it = g->begin();
    typename graph_type::const_iterator en = g->end();

    for ( ; it != en; ++it )
        {
            int theglobalrow = row+it->first;
            int thelocalrow = this->testSpace()->nLocalDofWithoutGhostOnProcStart(this->testSpace()->worldComm().globalRank(), test_index ) + it->second.get<1>();

            if (it->second.get<0>()!=g->worldComm().globalRank() )
                {
                    const int proc = it->second.get<0>();
                    const size_type realrow = this->testSpace()->dof()->firstDofGlobalCluster(proc)
                        + this->testSpace()->nLocalDofWithoutGhostOnProcStart(proc, test_index )
                        - mapOnTest.firstDofGlobalCluster(proc);
                    theglobalrow = realrow+it->first;
                    thelocalrow = this->testSpace()->nLocalDofWithoutGhostOnProcStart(proc, test_index ) + it->second.get<1>();
                }

            std::set<size_type>& row1_entries = M_graph->row( theglobalrow ).template get<2>();
            std::set<size_type> const& row2_entries = boost::get<2>( it->second );

            Debug( 5050 ) << "[mergeGraph] adding information to global row [" << theglobalrow << "], localrow=" << thelocalrow << "\n";
            M_graph->row( theglobalrow ).template get<1>() = thelocalrow;
            M_graph->row( theglobalrow ).template get<0>() = this->testSpace()->worldComm().mapLocalRankToGlobalRank()[it->second.get<0>()];

            if ( !row2_entries.empty() )
                {
                    for ( auto it = row2_entries.begin(), en = row2_entries.end() ; it!=en; ++it )
                        {
                            const auto dofcol = *it+col;
                            if (mapOnTrial.dofGlobalClusterIsOnProc(*it))
                                row1_entries.insert( dofcol );
                            else
                                {
                                    const int realproc = mapOnTrial.procOnGlobalCluster(*it);
                                    const size_type realcol = this->trialSpace()->dof()->firstDofGlobalCluster(realproc)
                                        + this->trialSpace()->nLocalDofWithoutGhostOnProcStart( realproc, trial_index )
                                        - mapOnTrial.firstDofGlobalCluster(realproc);

                                    row1_entries.insert(*it+realcol);
                                }
                        }
                }

        } // for( ; it != en; ++it )


}

template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraph( size_type hints )
{
    if ( (is_shared_ptr<typename test_space_type::mesh_ptrtype>::value && is_shared_ptr<typename trial_space_type::mesh_ptrtype>::value ) &&
         dynamic_cast<void*>( _M_X1->template mesh<0>().get() ) == dynamic_cast<void*>( _M_X2->template mesh<0>().get() ) )
        return this->computeGraph( hints, mpl::bool_<mpl::and_< mpl::bool_< ( test_space_type::nSpaces == 1 )>,
                                                                mpl::bool_< ( trial_space_type::nSpaces == 1 )> >::type::value >() );

    else
        return this->computeGraphInCaseOfInterpolate( hints, mpl::bool_<mpl::and_< mpl::bool_< ( test_space_type::nSpaces == 1 )>,
                                                                                   mpl::bool_< ( trial_space_type::nSpaces == 1 )> >::type::value >() );
}


template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraph( size_type hints, mpl::bool_<false> )
{
    boost::timer t;
    Debug( 5050 ) << "compute graph for composite bilinear form with interpolation\n";

    auto graph = computeGraph( hints, mpl::bool_< ( test_space_type::nSpaces > 1 )>(), mpl::bool_< ( trial_space_type::nSpaces > 1 )>() );

    Debug( 5050 ) << "closing graph for composite bilinear form with interpolation done in " << t.elapsed() << "s\n";
    t.restart();
    //graph->close();
    Debug( 5050 ) << "compute graph for composite bilinear form done in " << t.elapsed() << "s\n";

    return graph;
}

template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraph( size_type hints, mpl::bool_<true>, mpl::bool_<true> )
{
    fusion::for_each( _M_X1->functionSpaces(),
                      detail::compute_graph1<self_type>( this, hints ) );
    return M_graph;
}

template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraph( size_type hints, mpl::bool_<true>, mpl::bool_<false> )
{
    fusion::for_each( _M_X1->functionSpaces(),
                      detail::compute_graph3<self_type,trial_space_type>( this, _M_X2, 0, hints ) );
    return M_graph;
}

template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraph( size_type hints, mpl::bool_<false>, mpl::bool_<true> )
{
    fusion::for_each( _M_X2->functionSpaces(),
                      detail::compute_graph2<self_type,test_space_type>( this, _M_X1, 0, hints ) );
    return M_graph;
}

template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false> )
{
    boost::timer t;
    Debug( 5050 ) << "compute graph for composite bilinear form with interpolation\n";
    auto graph = computeGraphInCaseOfInterpolate( hints,
                                                  mpl::bool_< ( test_space_type::nSpaces > 1 )>(),
                                                  mpl::bool_< ( trial_space_type::nSpaces > 1 )>() );

    Debug( 5050 ) << "closing graph for composite bilinear form with interpolation done in " << t.elapsed() << "s\n";
    return graph;
}

template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<true> )
{
    fusion::for_each( _M_X1->functionSpaces(),
                      detail::compute_graph1<self_type>( this, hints ) );
    return M_graph;
}

template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<false> )
{
    fusion::for_each( _M_X1->functionSpaces(),
                      detail::compute_graph3<self_type,trial_space_type>( this, _M_X2, 0, hints ) );
    return M_graph;
}

template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false>, mpl::bool_<true> )
{
    fusion::for_each( _M_X2->functionSpaces(),
                      detail::compute_graph2<self_type,test_space_type>( this, _M_X1, 0, hints ) );
    return M_graph;
}




#if 0
template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraph( size_type hints, mpl::bool_<true> )
{
    boost::timer t;
    // Compute the sparsity structure of the global matrix.  This can be
    // fed into a PetscMatrix to allocate exacly the number of nonzeros
    // necessary to store the matrix.  This algorithm should be linear
    // in the (# of elements)*(# nodes per element)
    const size_type proc_id           = _M_X1->mesh()->comm().rank();
    const size_type n1_dof_on_proc    = _M_X1->nLocalDof();
    //const size_type n2_dof_on_proc    = _M_X2->nLocalDof();
    const size_type first1_dof_on_proc = _M_X1->dof()->firstDof( proc_id );
    const size_type last1_dof_on_proc = _M_X1->dof()->lastDof( proc_id );
    const size_type first2_dof_on_proc = _M_X2->dof()->firstDof( proc_id );
    const size_type last2_dof_on_proc = _M_X2->dof()->lastDof( proc_id );
    graph_ptrtype sparsity_graph( new graph_type( n1_dof_on_proc,
                                  first1_dof_on_proc, last1_dof_on_proc,
                                  first2_dof_on_proc, last2_dof_on_proc ) );

    typedef typename mesh_type::element_const_iterator mesh_element_const_iterator;

    mesh_element_const_iterator       elem_it  = _M_X1->mesh()->beginElementWithProcessId( proc_id );
    const mesh_element_const_iterator elem_en  = _M_X1->mesh()->endElementWithProcessId( proc_id );

    Feel::Context graph( hints );
    // If the user did not explicitly specify the DOF coupling
    // then all the DOFS are coupled to each other.  Furthermore,
    // we can take a shortcut and do this more quickly here.  So
    // we use an if-test.
    Debug( 5050 ) << "[computeGraph] test : " << ( graph.test ( Pattern::COUPLED ) || graph.test ( Pattern::EXTENDED ) ) << "\n";
    Debug( 5050 ) << "[computeGraph]  : graph.test ( Pattern::COUPLED )=" <<  graph.test ( Pattern::COUPLED ) << "\n";
    Debug( 5050 ) << "[computeGraph]  : graph.test ( Pattern::EXTENDED)=" <<  graph.test ( Pattern::EXTENDED ) << "\n";
#if 0

    if ( graph.test ( Pattern::COUPLED ) ||
            graph.test ( Pattern::EXTENDED ) )
#else
    if ( 1 )
#endif
    {
        Debug( 5050 ) << "[computeGraph] test (Pattern::COUPLED || Pattern::EXTENDED) ok\n";
        std::vector<size_type>
        element_dof1,
        element_dof2,
        neighbor_dof,
        dof_to_add;

        for ( ; elem_it != elem_en; ++elem_it )
        {
#if !defined(NDEBUG)
            Debug( 5050 ) << "[Stencil::computePatter] element " << elem_it->id() << " on proc " << elem_it->processId() << "\n";
#endif /* NDEBUG */
            mesh_element_type const& elem = *elem_it;

            // Get the global indices of the DOFs with support on this element
            element_dof1 = _M_X1->dof()->getIndices( elem.id() );
            element_dof2 = _M_X2->dof()->getIndices( elem.id() );

            // We can be more efficient if we sort the element DOFs
            // into increasing order
            std::sort( element_dof1.begin(), element_dof1.end() );
            std::sort( element_dof2.begin(), element_dof2.end() );

            const uint16_type  n1_dof_on_element = element_dof1.size();
            const uint16_type  n2_dof_on_element = element_dof2.size();

            for ( size_type i=0; i<n1_dof_on_element; i++ )
            {
                const size_type ig1 = element_dof1[i];

                // Only bother if this matrix row will be stored
                // on this processor.
                //if ((ig1 >= first1_dof_on_proc) &&
                //(ig1 <= last1_dof_on_proc))
                {
                    // This is what I mean
                    // assert ((ig - first_dof_on_proc) >= 0);
                    // but do the test like this because ig and
                    // first_dof_on_proc are size_types
#if 0
                    FEELPP_ASSERT ( ig1 >= first1_dof_on_proc )( ig1 )( first1_dof_on_proc ).error ( "invalid dof index" );
                    FEELPP_ASSERT ( ( ig1 - first1_dof_on_proc ) < sparsity_graph->size() )
                    ( ig1 )( first1_dof_on_proc )( sparsity_graph->size() ).error( "invalid dof index" );
#endif
                    graph_type::row_type& row = sparsity_graph->row( ig1 );
                    bool is_on_proc = ( ig1 >= first1_dof_on_proc ) && ( ig1 <= last1_dof_on_proc );
                    row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                    row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                    Debug( 5051 ) << "work with row " << ig1 << " local index " << ig1 - first1_dof_on_proc << "\n";

                    // If the row is empty we will add *all* the element DOFs,
                    // so just do that.
                    if ( row.get<2>().empty() )
                    {
                        row.get<2>() = element_dof2;
                    }

                    else
                    {
                        // Build a list of the DOF indices not found in the
                        // sparsity graph
                        dof_to_add.clear();

                        // Cache iterators.  Low will move forward, subsequent
                        // searches will be on smaller ranges
                        std::vector<size_type>::iterator
                        low  = std::lower_bound ( row.get<2>().begin(), row.get<2>().end(), element_dof2.front() ),
                        high = std::upper_bound ( low,         row.get<2>().end(), element_dof2.back() );

                        for ( size_type j=0; j<n2_dof_on_element; j++ )
                        {
                            const size_type jg = element_dof2[j];
                            //VLOG(1) << "computeGraph : ig:" << ig1 << ", lig: " << ig1-first1_dof_on_proc  << ", jg = " << jg << "\n";
#if 0
                            // See if jg is in the sorted range
                            std::pair<std::vector<size_type>::iterator,
                                std::vector<size_type>::iterator>
                                pos = std::equal_range ( low, high, jg );

                            // Must add jg if it wasn't found
                            if ( pos.first == pos.second )
                                dof_to_add.push_back( jg );

                            // pos.first is now a valid lower bound for any
                            // remaining element Dof. (That's why we sorted them.)
                            // Use it for the next search
                            low = pos.first;
#else
                            // See if jg is in the sorted range
                            std::pair<std::vector<size_type>::iterator,
                                std::vector<size_type>::iterator>
                                pos = std::equal_range ( row.get<2>().begin(), row.get<2>().end(), jg );

                            // Insert jg if it wasn't found
                            if ( pos.first == pos.second )
                                dof_to_add.push_back( jg );

#endif
                        }

                        // Add to the sparsity graph
                        if ( !dof_to_add.empty() )
                        {
                            const size_type old_size = row.get<2>().size();

                            row.get<2>().insert ( row.get<2>().end(),
                                                  dof_to_add.begin(),
                                                  dof_to_add.end() );

                            //std::inplace_merge (row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end());
                            sortSparsityRow ( row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end() );
                        }

                    }

                    // Now (possibly) add dof from neighboring elements
                    //if ( graph.test( Pattern::EXTENDED ) )
                    for ( uint16_type ms=0; ms < elem.nNeighbors(); ms++ )
                    {
                        mesh_element_type const* neighbor = NULL;
                        size_type neighbor_id = elem.neighbor( ms ).first;
                        size_type neighbor_process_id = elem.neighbor( ms ).second;

                        if ( neighbor_id != invalid_size_type_value )
                            //&& neighbor_process_id != proc_id )
                        {

#if 0
                            VLOG(1) << "element id " << elem.id()
                                    << ", element neighbor id " << neighbor_id
                                    << " in proc " << neighbor_process_id << "\n";
#endif
                            neighbor = boost::addressof( _M_X1->mesh()->element( neighbor_id,
                                                         neighbor_process_id ) );

#if 0
                            VLOG(1) << "found neighbor of element id " << elem.id()
                                    << ", element neighbor id " << neighbor->id()
                                    << " in proc " << neighbor->processId() << "\n";
#endif

                            if ( neighbor_id == neighbor->id()  )
                            {
                                neighbor_dof = _M_X2->dof()->getIndices( neighbor->id() );

                                const size_type n_dof_on_neighbor = neighbor_dof.size();
#if 0

                                for ( size_type j=0; j<n_dof_on_neighbor; j++ )
                                {
                                    Debug( 5051 ) << "neighbor elem id: " << neighbor->id() << " dof " << neighbor_dof[j] << "\n";
                                }

                                Debug( 5051 ) << "looking for dof " << ig1  << "\n";
#endif
#if 0
                                std::pair<std::vector<size_type>::iterator,
                                    std::vector<size_type>::iterator>
                                    posig = std::equal_range ( neighbor_dof.begin(), neighbor_dof.end(), ig1 );
#else
                                std::vector<size_type>::iterator posig = std::find( neighbor_dof.begin(), neighbor_dof.end(), ig1 );
#endif

                                // Insert jg if it wasn't found
                                //if (posig.first != posig.second)
                                if ( posig != neighbor_dof.end() ||
                                        graph.test ( Pattern::EXTENDED ) )

                                {
                                    //VLOG(1) << "found element in proc " << neighbor_process_id << " that shares dof\n";
                                    for ( size_type j=0; j<n_dof_on_neighbor; j++ )
                                    {
                                        const size_type jg = neighbor_dof[j];

#if 0
                                        // See if jg is in the sorted range
                                        std::pair<std::vector<size_type>::iterator,
                                            std::vector<size_type>::iterator>
                                            pos = std::equal_range ( row.get<2>().begin(), row.get<2>().end(), jg );
#else
                                        std::vector<size_type>::iterator pos = std::find( row.get<2>().begin(), row.get<2>().end(), jg );

#endif

                                        // Insert jg if it wasn't found
                                        if ( pos == row.get<2>().end() )
                                        {
                                            const size_type old_size = row.get<2>().size();
                                            row.get<2>().push_back ( jg );
                                            //std::inplace_merge (row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end());
                                            sortSparsityRow ( row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end() );
                                        }
                                    }
                                }
                            }
                        }

                    } // neighbor graph

#if 0

                    for ( int k = 0; k < row.get<2>().size(); ++k )
                        VLOG(1) << "row[ " << ig1 - first1_dof_on_proc << ","<< k << " ]=" << row.get<2>()[k] << "\n";

#endif
                } // only dof on proc

            }// dof loop
        } // element iterator loop
    }

    else
    {}

    Debug( 5050 ) << "[computeGraph<true>] before calling close in " << t.elapsed() << "s\n";
    //sparsity_graph->close();
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    return sparsity_graph;
}
#else
template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraph( size_type hints, mpl::bool_<true> )
{
    boost::timer t;
    // Compute the sparsity structure of the global matrix.  This can be
    // fed into a PetscMatrix to allocate exacly the number of nonzeros
    // necessary to store the matrix.  This algorithm should be linear
    // in the (# of elements)*(# nodes per element)
#if !defined(FEELPP_ENABLE_MPI_MODE) // NOT MPI
    const size_type proc_id           = _M_X1->mesh()->comm().rank();
    const size_type n1_dof_on_proc    = _M_X1->nLocalDof();
    //const size_type n2_dof_on_proc    = _M_X2->nLocalDof();
    const size_type first1_dof_on_proc = _M_X1->dof()->firstDof( proc_id );
    const size_type last1_dof_on_proc = _M_X1->dof()->lastDof( proc_id );
    const size_type first2_dof_on_proc = _M_X2->dof()->firstDof( proc_id );
    const size_type last2_dof_on_proc = _M_X2->dof()->lastDof( proc_id );
#else // MPI
    const size_type proc_id           = _M_X1->worldsComm()[0].localRank();
    const size_type n1_dof_on_proc    = _M_X1->nLocalDof();
    //const size_type n2_dof_on_proc    = _M_X2->nLocalDof();
    const size_type first1_dof_on_proc = _M_X1->dof()->firstDofGlobalCluster( proc_id );
    const size_type last1_dof_on_proc = _M_X1->dof()->lastDofGlobalCluster( proc_id );
    const size_type first2_dof_on_proc = _M_X2->dof()->firstDofGlobalCluster( proc_id );
    const size_type last2_dof_on_proc = _M_X2->dof()->lastDofGlobalCluster( proc_id );
#endif

    graph_ptrtype sparsity_graph( new graph_type( n1_dof_on_proc,
                                  first1_dof_on_proc, last1_dof_on_proc,
                                  first2_dof_on_proc, last2_dof_on_proc,
                                  _M_X1->worldComm() ) );

    if (_M_X1->nLocalDofWithoutGhost()==0 && _M_X2->nLocalDofWithoutGhost()==0 ) return sparsity_graph;

    auto elem_it  = _M_X1->mesh()->beginElementWithProcessId( _M_X1->mesh()->worldComm().localRank() /*proc_id*/ );
    auto elem_en  = _M_X1->mesh()->endElementWithProcessId( _M_X1->mesh()->worldComm().localRank() /*proc_id*/ );

    Feel::Context graph( hints );
    // If the user did not explicitly specify the DOF coupling
    // then all the DOFS are coupled to each other.  Furthermore,
    // we can take a shortcut and do this more quickly here.  So
    // we use an if-test.
    Debug( 5050 ) << "[computeGraph]  : graph.test ( Pattern::DEFAULT )=" <<  graph.test ( Pattern::DEFAULT ) << "\n";
    Debug( 5050 ) << "[computeGraph]  : graph.test ( Pattern::COUPLED )=" <<  graph.test ( Pattern::COUPLED ) << "\n";
    Debug( 5050 ) << "[computeGraph]  : graph.test ( Pattern::EXTENDED)=" <<  graph.test ( Pattern::EXTENDED ) << "\n";
    bool do_less =  ( ( graph.test( Pattern::DEFAULT ) &&
                        ( _M_X1->dof()->nComponents ==
                          _M_X2->dof()->nComponents ) ) &&
                      !graph.test( Pattern::COUPLED ) );
    std::vector<size_type>
    element_dof2( _M_X2->dof()->getIndicesSize() ),
                  neighbor_dof;

    for ( ; elem_it != elem_en; ++elem_it )
    {
#if !defined(NDEBUG)
        Debug( 5050 ) << "[Stencil::computePatter] element " << elem_it->id() << " on proc " << elem_it->processId() << "\n";
#endif /* NDEBUG */
        const auto & elem = *elem_it;

        // Get the global indices of the DOFs with support on this element
        //element_dof1 = _M_X1->dof()->getIndices( elem.id() );
#if !defined(FEELPP_ENABLE_MPI_MODE) // NOT MPI
        _M_X2->dof()->getIndicesSet( elem.id(), element_dof2 );
#else // MPI
        _M_X2->dof()->getIndicesSetOnGlobalCluster( elem.id(), element_dof2 );
#endif
        // We can be more efficient if we sort the element DOFs
        // into increasing order
        //std::sort(element_dof1.begin(), element_dof1.end());
        std::sort( element_dof2.begin(), element_dof2.end() );

        //const uint16_type  n1_dof_on_element = element_dof1.size();
        const uint16_type  n1_dof_on_element = _M_X1->dof()->getIndicesSize();
        const uint16_type  n2_dof_on_element = element_dof2.size();

        for ( size_type i=0; i<n1_dof_on_element; i++ )
            //BOOST_FOREACH( auto ig1, _M_X1->dof()->getIndices( elem.id() ) )
        {
#if !defined(FEELPP_ENABLE_MPI_MODE) // NOT MPI
            const size_type ig1 = _M_X1->dof()->localToGlobalId( elem.id(), i );
#else // MPI
            const size_type ig1 = _M_X1->dof()->mapGlobalProcessToGlobalCluster()[_M_X1->dof()->localToGlobalId( elem.id(), i )];
            auto theproc = _M_X1->dof()->procOnGlobalCluster( ig1 );
            // numLocal without ghosts ! very important for the graph with petsc
            const size_type il1 = ig1 - _M_X1->dof()->firstDofGlobalCluster( theproc );
#endif
            //const size_type ig1 = element_dof1[i];
            const int ndofpercomponent1 = n1_dof_on_element / _M_X1->dof()->nComponents;
            const int ncomp1 = i / ndofpercomponent1;
            const int ndofpercomponent2 = n2_dof_on_element / _M_X2->dof()->nComponents;

            {
                // This is what I mean
                // assert ((ig - first_dof_on_proc) >= 0);
                // but do the test like this because ig and
                // first_dof_on_proc are size_types
#if 0
                FEELPP_ASSERT ( ig1 >= first1_dof_on_proc )( ig1 )( first1_dof_on_proc ).error ( "invalid dof index" );
                FEELPP_ASSERT ( ( ig1 - first1_dof_on_proc ) < sparsity_graph->size() )
                ( ig1 )( first1_dof_on_proc )( sparsity_graph->size() ).error( "invalid dof index" );
#endif
                graph_type::row_type& row = sparsity_graph->row( ig1 );
#if !defined(FEELPP_ENABLE_MPI_MODE) // NOT MPI
                bool is_on_proc = ( ig1 >= first1_dof_on_proc ) && ( ig1 <= last1_dof_on_proc );
                row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
#else // MPI
                row.get<0>() = theproc ;
                row.get<1>() = il1;
#endif
                Debug( 5051 ) << "work with row " << ig1 << " local index " << ig1 - first1_dof_on_proc << "\n";

                if ( do_less )
                {
                    if ( ncomp1 == ( _M_X2->dof()->nComponents-1 ) )
                        row.get<2>().insert( element_dof2.begin()+ncomp1*ndofpercomponent2,
                                             element_dof2.end() );

                    else
                        row.get<2>().insert( element_dof2.begin()+ncomp1*ndofpercomponent2,
                                             element_dof2.begin()+( ncomp1+1 )*ndofpercomponent2 );
                }

                else
                {
                    row.get<2>().insert( element_dof2.begin(), element_dof2.end() );
                }

                // Now (possibly) add dof from neighboring elements
                if ( graph.test( Pattern::EXTENDED ) )
                {
                    for ( uint16_type ms=0; ms < elem.nNeighbors(); ms++ )
                    {
                        const auto * neighbor = boost::addressof( elem );
                        size_type neighbor_id = elem.neighbor( ms ).first;
                        size_type neighbor_process_id = elem.neighbor( ms ).second;

                        if ( neighbor_id != invalid_size_type_value )
                            //&& neighbor_process_id != proc_id )
                        {

                            neighbor = boost::addressof( _M_X1->mesh()->element( neighbor_id,
                                                         neighbor_process_id ) );

                            if ( neighbor_id == neighbor->id()  )
                            {
                                neighbor_dof = _M_X2->dof()->getIndices( neighbor->id() );

                                if ( do_less )
                                {
                                    if ( ncomp1 == ( _M_X2->dof()->nComponents-1 ) )
                                        row.get<2>().insert( neighbor_dof.begin()+ncomp1*ndofpercomponent2,
                                                             neighbor_dof.end() );

                                    else
                                        row.get<2>().insert( neighbor_dof.begin()+ncomp1*ndofpercomponent2,
                                                             neighbor_dof.begin()+( ncomp1+1 )*ndofpercomponent2 );

                                }

                                else
                                {
                                    row.get<2>().insert( neighbor_dof.begin(), neighbor_dof.end() );
                                }

                            } // neighbor_id
                        }

                    } // neighbor graph
                }
            } // only dof on proc

        }// dof loop
    } // element iterator loop

    Debug( 5050 )<< "[computeGraph<true>] before calling close in " << t.elapsed() << "s\n";
    //sparsity_graph->close();
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    return sparsity_graph;
}
#endif


template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true> )
{
    //std::cout << "\n start graphInterp "<< std::endl;

    typedef mpl::int_<20/*50*/> order_1d_type;
    typedef mpl::int_<12/*20*/> order_2d_type;
    typedef mpl::int_<8/*10*/> order_3d_type;

    typedef typename test_space_type::mesh_type test_mesh_type;
    typedef typename trial_space_type::mesh_type trial_mesh_type;
    typedef typename mpl::if_< mpl::equal_to<mpl::int_<test_mesh_type::nDim>,mpl::int_<1> >,
            order_1d_type,
            typename mpl::if_<mpl::equal_to<mpl::int_<test_mesh_type::nDim>,mpl::int_<2> >,
            order_2d_type,
            order_3d_type >::type>::type order_used_type;
    typedef typename mpl::if_<mpl::bool_<test_mesh_type::element_type::is_simplex>,
            mpl::identity<typename _Q<order_used_type::value>::template apply<test_mesh_type::element_type::nDim, typename test_mesh_type::value_type, Simplex>::type >,
                mpl::identity<typename _Q<order_used_type::value>::template apply<test_mesh_type::element_type::nDim, typename test_mesh_type::value_type, Hypercube>::type >
    >::type::type theim_type;

    typedef typename test_mesh_type::gm_type::precompute_type thepc_type;
    typedef typename test_mesh_type::gm_type::precompute_ptrtype thepc_ptrtype;
    typedef typename test_mesh_type::gm_type::template Context<vm::POINT, typename test_mesh_type::element_type> thegmc_type;
    typedef boost::shared_ptr<thegmc_type> thegmc_ptrtype;

    typedef typename test_mesh_type::Localization::matrix_node_type matrix_node_type;

    //-----------------------------------------------------------------------//

    //std::cout << "\n OrderUse " << order_used_type::value << std::endl;

    const size_type proc_id           = _M_X1->mesh()->comm().rank();
    const size_type n1_dof_on_proc    = _M_X1->nLocalDof();
    //const size_type n2_dof_on_proc    = _M_X2->nLocalDof();
    const size_type first1_dof_on_proc = _M_X1->dof()->firstDof( proc_id );
    const size_type last1_dof_on_proc = _M_X1->dof()->lastDof( proc_id );
    const size_type first2_dof_on_proc = _M_X2->dof()->firstDof( proc_id );
    const size_type last2_dof_on_proc = _M_X2->dof()->lastDof( proc_id );

    graph_ptrtype sparsity_graph( new graph_type( n1_dof_on_proc,
                                  first1_dof_on_proc, last1_dof_on_proc,
                                  first2_dof_on_proc, last2_dof_on_proc ) );

    typedef typename test_mesh_type::element_const_iterator mesh_element_const_iterator;
    mesh_element_const_iterator       elem_it  = _M_X1->mesh()->beginElementWithProcessId( proc_id );
    const mesh_element_const_iterator elem_en  = _M_X1->mesh()->endElementWithProcessId( proc_id );


    //-----------------------------------------------------------------------//
    // init localisation tools
    auto locToolForXh2 = _M_X2->mesh()->tool_localization();
    locToolForXh2->updateForUse();
    bool doExtrapolationAtStartXh2 = locToolForXh2->doExtrapolation();
    if (doExtrapolationAtStartXh2) locToolForXh2->setExtrapolation( false );
    //locTool->kdtree()->nbNearNeighbor(_M_X2->mesh()->numElements() );
    auto locToolForXh1 = _M_X1->mesh()->tool_localization();
    locToolForXh1->updateForUse();
    bool doExtrapolationAtStartXh1 = locToolForXh1->doExtrapolation();
    if (doExtrapolationAtStartXh1) locToolForXh1->setExtrapolation( false );


    matrix_node_type ptsReal( elem_it->vertices().size1(), 1 );
    size_type IdEltInXh2 = invalid_size_type_value;
    //node_type trialNodeRef,testNodeRef;

#if FEELPP_EXPORT_GRAPH
    std::map<size_type,std::list<size_type> > mapBetweenMeshes;
#endif

    std::vector<size_type> element_dof1, element_dof2;
    std::set<size_type> neighLocalizedInXh1;

    std::set<size_type > listTup;

    theim_type theim;
    thepc_ptrtype geopc( new thepc_type( elem_it->gm(), theim.points() ) );
    thegmc_ptrtype gmc( new thegmc_type( elem_it->gm(), *elem_it, geopc ) );


    //-----------------------------------------------------------------------//

    if ( _M_X1->nDof()>1 )
    {
        for ( ; elem_it != elem_en; ++elem_it )
        {
            auto const& elem = *elem_it;

            // Get the global indices of the DOFs with support on this element
            element_dof1 = _M_X1->dof()->getIndices( elem.id() );

            const uint16_type n1_dof_on_element = element_dof1.size();

            std::vector<boost::tuple<bool,size_type> > hasFinds( n1_dof_on_element,boost::make_tuple( false,invalid_size_type_value ) );

            for ( size_type i=0; i<n1_dof_on_element; i++ )
            {
                const size_type ig1 = element_dof1[i];
                auto const ptRealDof = boost::get<0>( _M_X1->dof()->dofPoint( ig1 ) );

                ublas::column(ptsReal,0 ) = ptRealDof;
                auto resLocalisationInXh2 = locToolForXh2->run_analysis(ptsReal,IdEltInXh2,elem_it->vertices(),mpl::int_<0>());
                IdEltInXh2 = resLocalisationInXh2.template get<1>();
                bool hasFind = resLocalisationInXh2.template get<0>()[0];

                listTup.clear();

                if ( hasFind )
                {
                    listTup.insert( IdEltInXh2 );
                    hasFinds[i] = boost::make_tuple( true,IdEltInXh2 );
                    // maybe is on boundary->more elts
                    //size_type idElt1 = elem.id();
                    //size_type idElt2 = resTemp.template get<1>();
                    auto const& geoelt2 = _M_X2->mesh()->element( IdEltInXh2 );
                    std::vector<size_type> neighbor_ids;

                    for ( uint16_type ms=0; ms < geoelt2.nNeighbors(); ms++ )
                    {
                        size_type neighbor_id = geoelt2.neighbor( ms ).first;

                        if ( neighbor_id!=invalid_size_type_value )
                            neighbor_ids.push_back( neighbor_id );
                    }

                    auto resIsIn = locToolForXh2->isIn( neighbor_ids,ptRealDof );
                    uint16_type cpt=0;

                    for ( auto it=resIsIn.template get<1>().begin(),en=resIsIn.template get<1>().end(); it<en; ++it,++cpt )
                    {
                        if ( *it )
                        {
                            listTup.insert( neighbor_ids[cpt] );
                        }
                    }

                    auto res_it = listTup.begin();
                    auto res_en = listTup.end();
                    for ( ; res_it != res_en ; ++res_it )
                    {
#if FEELPP_EXPORT_GRAPH
                        mapBetweenMeshes[*res_it].push_back( elem.id() );
                        //std::cout << " test id " << *res_it << " trial id " << elem.id() << std::endl;
#endif
                        // not efficient but sometimes necessary
                        for ( size_type ii=0; ii<n1_dof_on_element; ii++ )
                            {
                                const size_type ig1ongraph = element_dof1[ii];

                                element_dof2 = _M_X2->dof()->getIndices( *res_it );

                                graph_type::row_type& row = sparsity_graph->row( ig1ongraph );
                                bool is_on_proc = ( ig1ongraph >= first1_dof_on_proc ) && ( ig1ongraph <= last1_dof_on_proc );
                                row.template get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                                row.template get<1>() = is_on_proc?ig1ongraph - first1_dof_on_proc:invalid_size_type_value;
                                //if ( do_less ) {}
                                row.template get<2>().insert( element_dof2.begin(), element_dof2.end() );
                            }
                    }//res
                } // if (hasFind)

                else
                {
#if 0
                    //std::cout << "\n not find"<<std::endl;
                    // row empty
                    graph_type::row_type& row = sparsity_graph->row( ig1 );
                    bool is_on_proc = ( ig1 >= first1_dof_on_proc ) && ( ig1 <= last1_dof_on_proc );
                    row.template get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                    row.template get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                    row.template get<2>().clear();
#endif
                }
            } // for (size_type i=0; i<n1_dof_on_element; i++)

            uint16_type nbFind = hasFinds.size();
            bool doQ=false;
            size_type id1,id2;

            for ( uint16_type nF = 0 ; nF<nbFind && !doQ ; ++nF )
                if ( hasFinds[nF].template get<0>() )
                {
                    id1=hasFinds[nF].template get<1>();

                    for ( uint16_type nF2 = 0 ; nF2<nbFind && !doQ  ; ++nF2 )
                    {
                        if ( hasFinds[nF2].template get<0>() )
                        {
                            id2=hasFinds[nF2].template get<1>();

                            if ( id1!=id2 )
                                doQ=true; //if( !_M_X2->mesh()->element(id1).isNeighbor(_M_X2->mesh()->element(id2) ) )  doQ=true;
                        }
                    }
                }


            if ( doQ )
            {

                //theim_type theim;
                //thepc_ptrtype geopc( new thepc_type( elem.gm(), theim.points() ) );
                //thegmc_ptrtype gmc( new thegmc_type( elem.gm(), elem, geopc ) );
                gmc->update( elem );

                for ( int q = 0; q <  gmc->nPoints(); ++ q )
                {
                    ublas::column(ptsReal,0 ) = gmc->xReal( q );
                    //auto const resQuad = locToolForXh2->run_analysis(ptsReal,IdEltInXh2,elem_it->vertices(),mpl::int_<0>());
                    auto const resQuad = locToolForXh2->searchElement( gmc->xReal( q ) );

                    if ( resQuad.template get<0>() )
                    {
                        IdEltInXh2 = resQuad.template get<1>();
#if FEELPP_EXPORT_GRAPH
                        mapBetweenMeshes[IdEltInXh2].push_back( elem.id() );
#endif
                        element_dof2 = _M_X2->dof()->getIndices( IdEltInXh2 );

                        for ( size_type i=0; i<n1_dof_on_element; i++ )
                        {
                            const size_type ig1 = element_dof1[i];
                            graph_type::row_type& row = sparsity_graph->row( ig1 );
                            bool is_on_proc = ( ig1 >= first1_dof_on_proc ) && ( ig1 <= last1_dof_on_proc );
                            row.template get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                            row.template get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                            //if ( do_less ) {}
                            row.template get<2>().insert( element_dof2.begin(), element_dof2.end() );
                        }
                    }
                } // for ( int q = 0; q <  gmc->nPoints(); ++ q )
            } // if (doQ)
        } // for ( ; elem_it ... )
    } //Xh1->nof >1

    else //case Xh1->nof == 1
    {
        for ( ; elem_it != elem_en; ++elem_it )
        {
            /*mesh_element_1_type*/auto const& elem = *elem_it;

            // Get the global indices of the DOFs with support on this element
            element_dof1 = _M_X1->dof()->getIndices( elem.id() );

            const uint16_type nPtGeo = elem.G().size2();
            std::vector<boost::tuple<bool,size_type> > hasFinds( nPtGeo,boost::make_tuple( false,invalid_size_type_value ) );

            for ( size_type i=0; i<nPtGeo; i++ )
            {
                const size_type ig1 = element_dof1[0];
                //const int ndofpercomponent1 = n1_dof_on_element / _M_X1->dof()->nComponents;
                //const int ncomp1 = i / ndofpercomponent1;
                //const int ndofpercomponent2 = n2_dof_on_element / _M_X2->dof()->nComponents;

                typename matrix_node<typename test_mesh_type::value_type>::type ptReal( test_mesh_type::nRealDim , 1 );
                ublas::column( ptReal ,0 ) = ublas::column( elem.G(),i );

                //auto res = locTool->searchElements(ublas::column( elem.G(),i )/*ptReal*/);
                //auto hasFind = res.template get<0>();

                auto resTemp = locToolForXh2->searchElement( ublas::column( elem.G(),i ) );
                bool hasFind = resTemp.template get<0>();
                std::set<size_type > listTup;

                if ( hasFind )
                {
                    listTup.insert( resTemp.template get<1>() );
                    hasFinds[i] = boost::make_tuple( true,resTemp.template get<1>() );
                    // maybe is on boundary->more elts
                    //size_type idElt1 = elem.id();
                    size_type idElt2 = resTemp.template get<1>();
                    auto const& geoelt2 = _M_X2->mesh()->element( idElt2 );
                    std::vector<size_type> neighbor_ids( geoelt2.nNeighbors() ); //neighbor_ids.clear();//(geoelt2.nNeighbors());

                    for ( uint16_type ms=0; ms < geoelt2.nNeighbors(); ms++ )
                    {
                        size_type neighbor_id = geoelt2.neighbor( ms ).first;

                        if ( neighbor_id!=invalid_size_type_value ) neighbor_ids.push_back( neighbor_id );

                        //neighbor_ids[ms]=neighbor_id;
                    }

                    auto resIsIn = locToolForXh2->isIn( neighbor_ids,ublas::column( elem.G(),i ) );
                    uint16_type cpt=0;

                    for ( auto it=resIsIn.template get<1>().begin(),en=resIsIn.template get<1>().end(); it<en; ++it,++cpt )
                    {
                        if ( *it )
                        {
                            listTup.insert( neighbor_ids[cpt] );
                        }

                        else if ( test_mesh_type::nDim==trial_mesh_type::nDim ) // pt n'est pas chez le voison, peut-etre les pts G()(si dim=dim!!!)
                        {
                            for ( auto it_neigh = neighbor_ids.begin(), en_neigh = neighbor_ids.end() ; it_neigh != en_neigh ; ++it_neigh )
                            {
                                if ( ( sparsity_graph->row( ig1 ) ).template get<2>().find( *it_neigh )==( sparsity_graph->row( ig1 ) ).template get<2>().end() )
                                {
                                    if ( neighLocalizedInXh1.find( *it_neigh )==neighLocalizedInXh1.end() )
                                    {
                                        neighLocalizedInXh1.insert( *it_neigh );
                                        bool findNeih=false;
                                        auto const& geoeltNEW = _M_X2->mesh()->element( *it_neigh );
                                        const uint16_type nPtGeoBis = _M_X2->mesh()->element( *it_neigh ).G().size2();

                                        for ( size_type iii=0; iii<nPtGeoBis && !findNeih ; iii++ )
                                        {
                                            auto  resBis = locToolForXh1->searchElement( ublas::column( geoeltNEW.G(),iii ) );

                                            if ( resBis.template get<0>() )
                                            {
                                                listTup.insert( *it_neigh );
                                                findNeih=true;
                                            }
                                        }
                                    }
                                }
                            }
                        } //else
                    } // for (auto it=resIsIn...

                    auto res_it = listTup.begin();
                    auto res_en = listTup.end();
                    for ( ; res_it != res_en ; ++res_it )
                    {
#if FEELPP_EXPORT_GRAPH
                        //mapBetweenMeshes[elem.id()].push_back(*res_it);
                        mapBetweenMeshes[*res_it/*->get<0>()*/].push_back( elem.id() );
#endif


                        element_dof2 = _M_X2->dof()->getIndices( *res_it/*->get<0>()*/ );

                        graph_type::row_type& row = sparsity_graph->row( ig1 );
                        bool is_on_proc = ( ig1 >= first1_dof_on_proc ) && ( ig1 <= last1_dof_on_proc );
                        row.template get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                        row.template get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                        //if ( do_less ) {}
                        row.template get<2>().insert( element_dof2.begin(), element_dof2.end() );
                    }
                } // hasFind
            } // for (size_type i=0; i<nPtGeo; i++)

            uint16_type nbFind = hasFinds.size();
            bool doQ=false;
            size_type id1,id2;

            for ( uint16_type nF = 0 ; nF<nbFind && !doQ ; ++nF )
                if ( hasFinds[nF].template get<0>() )
                {
                    id1=hasFinds[nF].template get<1>();

                    for ( uint16_type nF2 = 0 ; nF2<nbFind && !doQ  ; ++nF2 )
                    {
                        if ( hasFinds[nF2].template get<0>() )
                        {
                            id2=hasFinds[nF2].template get<1>();

                            if ( id1!=id2 )
                                doQ=true; //if( !_M_X2->mesh()->element(id1).isNeighbor(_M_X2->mesh()->element(id2) ) )  doQ=true;
                        }
                    }
                }

            if ( doQ )
            {
                //theim_type theim;
                //thepc_ptrtype geopc( new thepc_type( elem.gm(), theim.points() ) );
                //thegmc_ptrtype gmc( new thegmc_type( elem.gm(), elem, geopc ) );
                gmc->update( elem );

                for ( int q = 0; q <  gmc->nPoints(); ++ q )
                {
                    auto resQuad = locToolForXh2->searchElement( gmc->xReal( q ) );

                    if ( resQuad.template get<0>() )
                    {
#if FEELPP_EXPORT_GRAPH
                        mapBetweenMeshes[resQuad.template get<1>()].push_back( elem.id() );
#endif
                        element_dof2 = _M_X2->dof()->getIndices( resQuad.template get<1>() );
                        //for (size_type i=0; i<n1_dof_on_element; i++)
                        //   {
                        const size_type ig1 = element_dof1[0];
                        graph_type::row_type& row = sparsity_graph->row( ig1 );
                        bool is_on_proc = ( ig1 >= first1_dof_on_proc ) && ( ig1 <= last1_dof_on_proc );
                        row.template get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                        row.template get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                        //if ( do_less ) {}
                        row.template get<2>().insert( element_dof2.begin(), element_dof2.end() );
                        //    }
                    }
                }
            } // if (doQ)


        } //  for ( ; elem_it != elem_en; ++elem_it)
    } // M_X1->nDof==1

    if (doExtrapolationAtStartXh1) locToolForXh1->setExtrapolation( true );
    if (doExtrapolationAtStartXh2) locToolForXh2->setExtrapolation( true );
    //locTool->kdtree()->nbNearNeighbor( 15 );

    //sparsity_graph->close();

#if FEELPP_EXPORT_GRAPH
#if 0
    typedef mesh_1_type mesh_export_type;
    typedef FunctionSpace<mesh_1_type, bases<Lagrange<0, Scalar,Discontinuous> > > space_disc_type;
    auto spaceGraphProj = space_disc_type::New( _M_X1->mesh() );
    auto elem_itt  = _M_X1->mesh()->beginElementWithProcessId( proc_id );
    auto elem_ent  = _M_X1->mesh()->endElementWithProcessId( proc_id );
#else
    typedef trial_mesh_type mesh_export_type;
    typedef FunctionSpace<mesh_export_type, bases<Lagrange<0, Scalar,Discontinuous> > > space_disc_type;
    auto spaceGraphProj = space_disc_type::New( _M_X2->mesh() );
    auto elem_itt  = _M_X2->mesh()->beginElementWithProcessId( proc_id );
    auto elem_ent  = _M_X2->mesh()->endElementWithProcessId( proc_id );
#endif
    auto graphProj = spaceGraphProj->element();
    graphProj.zero();

    for ( ; elem_itt != elem_ent; ++elem_itt )
    {
        if ( mapBetweenMeshes.find( elem_itt->id() ) != mapBetweenMeshes.end() )
        {
            if ( mapBetweenMeshes.find( elem_itt->id() )->second.size()>0 )
            {
                element_dof1 = spaceGraphProj->dof()->getIndices( elem_itt->id() );
                const uint16_type n1_dof_on_element = element_dof1.size();
                //std::cout << "\nn1_dof_on_element " << n1_dof_on_element << std::endl;
                for ( uint16_type i=0; i<n1_dof_on_element; i++ )
                {
                    const size_type ig1 = element_dof1[i];
                    //std::cout << "ig1 "<<ig1 << " graphProj.nDof() "<< graphProj.nDof() <<std::endl;
                    graphProj( ig1 ) = 1; //graphProj.set(0, 1. );
                }
            }
        }
    }

    auto exporter = boost::shared_ptr<Feel::Exporter<mesh_export_type> >( Feel::Exporter<mesh_export_type>::New( "ensight", "ExportGraph" ) );
    exporter->step( 0 )->setMesh( graphProj.mesh() );
    exporter->step( 0 )->add( "graphProj", graphProj );
    exporter->save();
#endif
    //std::cout << "\n finish graphInterp "<< std::endl;
    return sparsity_graph;
}


}
#endif
#endif /* __galerkingraph_H */
