/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
            auto thestencil = stencil(_test=space2, _trial=M_space1, _pattern=M_hints );

            M_stencil->mergeGraph( M_stencil->testSpace()->nDofStart( M_test_index ), M_stencil->trialSpace()->nDofStart( M_trial_index ) , thestencil->graph() );

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
            auto thestencil = stencil(_test=M_space1, _trial=space2, _pattern=M_hints );

            M_stencil->mergeGraph( M_stencil->testSpace()->nDofStart( M_test_index ), M_stencil->trialSpace()->nDofStart( M_trial_index ) , thestencil->graph() );

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

    Stencil( test_space_ptrtype Xh, trial_space_ptrtype Yh, size_type graph_hints  )
        :
        _M_X1( Xh ),
        _M_X2( Yh ),
        M_graph( new graph_type( Xh->nLocalDof(),
                                 Xh->nDofStart(), Xh->nDofStart()+Xh->nLocalDof(),
                                 Yh->nDofStart(), Yh->nDofStart()+Yh->nLocalDof() ) )
        {
            const size_type n1_dof_on_proc = _M_X1->nLocalDof();

            boost::timer t;
            //std::cout << "compute graph\n";

            if ( dynamic_cast<void*>( _M_X1->mesh().get()) == dynamic_cast<void*>( _M_X2->mesh().get()) )
                M_graph = computeGraph( graph_hints, mpl::bool_<mpl::and_< mpl::bool_< (test_space_type::nSpaces == 1)>,
                                                                         mpl::bool_< (trial_space_type::nSpaces == 1)> >::type::value >() );
            else
                M_graph = computeGraphInCaseOfInterpolate( graph_hints, mpl::bool_<mpl::and_< mpl::bool_< (test_space_type::nSpaces == 1)>,
                                                                                            mpl::bool_< (trial_space_type::nSpaces == 1)> >::type::value >() );
            M_graph->close();
            //std::cout << "computed graph in " << t.elapsed() << "s\n"; t.restart();
        }
    Stencil( test_space_ptrtype Xh, trial_space_ptrtype Yh, size_type graph_hints, graph_ptrtype g )
        :
        _M_X1( Xh ),
        _M_X2( Yh ),
        M_graph( g )
        {}

    graph_ptrtype computeGraph( size_type hints, mpl::bool_<true> );
    graph_ptrtype computeGraph( size_type hints, mpl::bool_<false> );
    graph_ptrtype computeGraph( size_type hints, mpl::bool_<true>, mpl::bool_<true> );
    graph_ptrtype computeGraph( size_type hints, mpl::bool_<false>, mpl::bool_<true> );
    graph_ptrtype computeGraph( size_type hints, mpl::bool_<true>, mpl::bool_<false> );


    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true> );

    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false> )
        {
            boost::timer t;
            Debug( 5050 ) << "compute graph for composite bilinear form with interpolation\n";

            auto graph = computeGraphInCaseOfInterpolate( hints, mpl::bool_< ( test_space_type::nSpaces > 1)>(), mpl::bool_< ( trial_space_type::nSpaces > 1)>() );

            Debug( 5050 ) << "closing graph for composite bilinear form with interpolation done in " << t.elapsed() << "s\n"; t.restart();
            graph->close();
            Debug( 5050 ) << "compute graph for composite bilinear form done in " << t.elapsed() << "s\n";

            return graph;
        }


    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<true> )
        {
            fusion::for_each( _M_X1->functionSpaces(), detail::compute_graph1<self_type>( this, hints ) );

            return M_graph;
        }
    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false>, mpl::bool_<true> )
        {
            fusion::for_each( _M_X2->functionSpaces(),
                              detail::compute_graph2<self_type,test_space_type>( this, _M_X1, 0, hints ) );
            return M_graph;
        }

    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<false> )
        {
            fusion::for_each( _M_X1->functionSpaces(),
                              detail::compute_graph3<self_type,trial_space_type>( this, _M_X2, 0, hints ) );
            return M_graph;
        }


    void mergeGraph( int row, int col, graph_ptrtype g );

    test_space_ptrtype testSpace() const { return _M_X1; }
    trial_space_ptrtype trialSpace() const { return _M_X1; }
    graph_ptrtype graph() const { return M_graph; }
    graph_ptrtype graph() { return M_graph; }
private:
    test_space_ptrtype _M_X1;
    trial_space_ptrtype _M_X2;
    graph_ptrtype M_graph;
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
    public std::map<boost::tuple<boost::shared_ptr<FunctionSpaceBase>,boost::shared_ptr<FunctionSpaceBase>,size_type>, boost::shared_ptr<GraphCSR> >,
    public boost::noncopyable
{
public:
    typedef boost::shared_ptr<GraphCSR> graph_ptrtype;
    typedef boost::tuple<boost::shared_ptr<FunctionSpaceBase>,boost::shared_ptr<FunctionSpaceBase>,size_type> key_type;
    typedef std::map<key_type, graph_ptrtype> graph_manager_type;

};

typedef Feel::Singleton<StencilManagerImpl> StencilManager;

BOOST_PARAMETER_FUNCTION(
    (typename detail::compute_stencil_type<Args>::ptrtype), // 1. return type
    stencil,                                       // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    (required                                   // 4. one required parameter, and
     (test,             *(boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> >))
     (trial,             *(boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> >))
        )
    (optional                                   //    four optional parameters, with defaults
     (pattern,             *(boost::is_integral<mpl::_>), Pattern::COUPLED )
        )
    )
{
    Feel::detail::ignore_unused_variable_warning(args);
    typedef typename detail::compute_stencil_type<Args>::ptrtype stencil_ptrtype;
    typedef typename detail::compute_stencil_type<Args>::type stencil_type;

    // we look into the spaces dictionary for existing graph
    auto git = StencilManager::instance().find( boost::make_tuple( trial, test, pattern ) );
    if (  git != StencilManager::instance().end() )
    {
        //std::cout << "Found a  stencil in manager (" << test.get() << "," << trial.get() << "," << pattern << ")\n";
        auto s = stencil_ptrtype( new stencil_type( test, trial, pattern, git->second ) );
        return s;
    }
    else
    {
        //std::cout << "Creating a new stencil in manager (" << test.get() << "," << trial.get() << "," << pattern << ")\n";
        auto s = stencil_ptrtype( new stencil_type( test, trial, pattern ) );
        StencilManager::instance().operator[](boost::make_tuple( trial, test, pattern )) = s->graph();
        return s;
    }
}


namespace detail
{
template<typename BidirectionalIterator>
inline
void
sortSparsityRow (const BidirectionalIterator begin,
                 BidirectionalIterator       middle,
                 const BidirectionalIterator end)
{
    if ((begin == middle) || (middle == end)) return;

    assert (std::distance (begin,  middle) > 0);
    assert (std::distance (middle, end)    > 0);
    FEEL_ASSERT( std::unique (begin,  middle) == middle )
        ( *begin )( *middle ).error( "duplicate dof(begin,middle)" );
    FEEL_ASSERT (std::unique (middle, end)    == end)
        (*begin)( *middle ).error( "duplicate dof (middle,end)" );

    while (middle != end)
        {
            BidirectionalIterator
                b = middle,
                a = b-1;

            // Bubble-sort the middle value downward
            while (!(*a < *b)) // *a & *b are less-than comparable, so use <
                {
                    std::swap (*a, *b);

                    if (a == begin) break;

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

        for (++first; first != end; prev=first, ++first)
            if (*first < *prev)
                assert(false);
    }
#endif

    // Make sure the two ranges did not contain any common elements
    assert (std::unique (begin, end) == end);
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
            for( ; it != en; ++it )
                {
                    std::vector<size_type> const& ivec = boost::get<2>( it->second );
                    for( int i = 0;i < ivec.size(); ++i )
                        {
                            Debug( 5050 ) << "[mergeGraph] ivec[" << i << "] = " << ivec[i] << "\n";
                        }
                }
#endif // 9
        }
    else
        {
            M_graph->setLastRowEntryOnProc( row + g->lastRowEntryOnProc() );
            M_graph->setLastColEntryOnProc( col + g->lastColEntryOnProc() );

            Debug( 5050 ) << "[merge graph] already something in store\n";
            typename graph_type::const_iterator it = g->begin();
            typename graph_type::const_iterator en = g->end();
            for( ; it != en; ++it )
            {
                int theglobalrow = row+it->first;
                int thelocalrow = row + boost::get<1>( it->second );
                //auto row1_entries = boost::unwrap_ref( boost::ref( M_graph->row(theglobalrow).template get<2>() ) );
                std::set<size_type>& row1_entries = M_graph->row(theglobalrow).template get<2>();
                std::set<size_type> const& row2_entries = boost::get<2>( it->second );

                Debug( 5050 ) << "[mergeGraph] adding information to global row [" << theglobalrow << "], localrow=" << thelocalrow << "\n";
                M_graph->row(theglobalrow).template get<1>() = thelocalrow;

                if ( row1_entries.empty() )
                {
                    // if row is empty then no need to shift the dof in
                    // composite case since the merge in done block-row-wise
                    row1_entries = row2_entries;
                }
                else
                {
                    // ensure unique sorted ids
                    auto itg = boost::prior(row1_entries.end());
                    // shift dofs in case of composite spaces
                    std::for_each( row2_entries.begin(), row2_entries.end(),[&]( size_type o ){ itg = row1_entries.insert( itg, o+col); });
                }
            }
        }
    Debug( 5050 ) << " -- merge_graph (" << row << "," << col << ") in " << tim.elapsed() << "\n";
    Debug( 5050 ) << "merge graph for composite bilinear form done\n";
}




template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraph( size_type hints, mpl::bool_<false> )
{
    boost::timer t;
    Debug( 5050 ) << "compute graph for composite bilinear form with interpolation\n";

    auto graph = computeGraph( hints, mpl::bool_< ( test_space_type::nSpaces > 1)>(), mpl::bool_< ( trial_space_type::nSpaces > 1)>() );

    Debug( 5050 ) << "closing graph for composite bilinear form with interpolation done in " << t.elapsed() << "s\n"; t.restart();
    graph->close();
    Debug( 5050 ) << "compute graph for composite bilinear form done in " << t.elapsed() << "s\n";

    return graph;
}

template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraph( size_type hints, mpl::bool_<true>, mpl::bool_<true> )
{
    fusion::for_each( _M_X1->functionSpaces(), detail::compute_graph1<self_type>( this, hints ) );

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

#if 0 //vincent
template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraph( size_type hints, mpl::bool_<false> )
{
    boost::timer t;
    Debug( 5050 ) << "compute graph for composite bilinear form\n";
    fusion::for_each( _M_X1->functionSpaces(), detail::compute_graph1<self_type>( this, hints ) );
    Debug( 5050 ) << "closing graph for composite bilinear form done in " << t.elapsed() << "s\n"; t.restart();

#if 0
    typename graph_type::const_iterator it = M_graph->begin();
    typename graph_type::const_iterator en = M_graph->end();
    for( ; it != en; ++it )
        {
            std::cout << "row " << it->first << ", " <<  boost::get<1>( it->second ) << ": ";
            std::copy( boost::get<2>( it->second ).begin(), boost::get<2>( it->second ).end(), std::ostream_iterator<int>( std::cout, " " ) );
            std::cout << "\n";
        }
#endif
    M_graph->close();
    Debug( 5050 ) << "compute graph for composite bilinear form done in " << t.elapsed() << "s\n";


    return M_graph;
}
#endif

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
    const size_type nprocs           = _M_X1->mesh()->comm().size();
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
        if (1)
#endif
        {
            Debug( 5050 ) << "[computeGraph] test (Pattern::COUPLED || Pattern::EXTENDED) ok\n";
            std::vector<size_type>
                element_dof1,
                element_dof2,
                neighbor_dof,
                dof_to_add;

            for ( ; elem_it != elem_en; ++elem_it)
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
                    std::sort(element_dof1.begin(), element_dof1.end());
                    std::sort(element_dof2.begin(), element_dof2.end());

                    const uint16_type  n1_dof_on_element = element_dof1.size();
                    const uint16_type  n2_dof_on_element = element_dof2.size();

                    for (size_type i=0; i<n1_dof_on_element; i++)
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
                                    FEEL_ASSERT (ig1 >= first1_dof_on_proc )( ig1 )( first1_dof_on_proc ).error ("invalid dof index");
                                    FEEL_ASSERT ((ig1 - first1_dof_on_proc) < sparsity_graph->size() )
                                        ( ig1 )( first1_dof_on_proc )( sparsity_graph->size() ).error( "invalid dof index" );
#endif
                                    graph_type::row_type& row = sparsity_graph->row(ig1);
                                    bool is_on_proc = ( ig1 >= first1_dof_on_proc) && (ig1 <= last1_dof_on_proc);
                                    row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                                    row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                                    Debug( 5051 ) << "work with row " << ig1 << " local index " << ig1 - first1_dof_on_proc << "\n";

                                    // If the row is empty we will add *all* the element DOFs,
                                    // so just do that.
                                    if (row.get<2>().empty())
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
                                                low  = std::lower_bound (row.get<2>().begin(), row.get<2>().end(), element_dof2.front()),
                                                high = std::upper_bound (low,         row.get<2>().end(), element_dof2.back());

                                            for (size_type j=0; j<n2_dof_on_element; j++)
                                                {
                                                    const size_type jg = element_dof2[j];
                                                    //Debug() << "computeGraph : ig:" << ig1 << ", lig: " << ig1-first1_dof_on_proc  << ", jg = " << jg << "\n";
#if 0
                                                    // See if jg is in the sorted range
                                                    std::pair<std::vector<size_type>::iterator,
                                                        std::vector<size_type>::iterator>
                                                        pos = std::equal_range (low, high, jg);

                                                    // Must add jg if it wasn't found
                                                    if (pos.first == pos.second)
                                                        dof_to_add.push_back(jg);

                                                    // pos.first is now a valid lower bound for any
                                                    // remaining element Dof. (That's why we sorted them.)
                                                    // Use it for the next search
                                                    low = pos.first;
#else
                                                    // See if jg is in the sorted range
                                                    std::pair<std::vector<size_type>::iterator,
                                                        std::vector<size_type>::iterator>
                                                        pos = std::equal_range (row.get<2>().begin(), row.get<2>().end(), jg);

                                                    // Insert jg if it wasn't found
                                                    if (pos.first == pos.second)
                                                        dof_to_add.push_back(jg);
#endif
                                                }

                                            // Add to the sparsity graph
                                            if (!dof_to_add.empty())
                                                {
                                                    const size_type old_size = row.get<2>().size();

                                                    row.get<2>().insert (row.get<2>().end(),
                                                                dof_to_add.begin(),
                                                                dof_to_add.end());

                                                    //std::inplace_merge (row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end());
                                                    sortSparsityRow (row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end());
                                                }

                                        }

                                    // Now (possibly) add dof from neighboring elements
                                    //if ( graph.test( Pattern::EXTENDED ) )
                                        for (uint16_type ms=0; ms < elem.nNeighbors(); ms++)
                                            {
                                                mesh_element_type const* neighbor = NULL;
                                                size_type neighbor_id = elem.neighbor(ms).first;
                                                size_type neighbor_process_id = elem.neighbor(ms).second;
                                                if ( neighbor_id != invalid_size_type_value )
                                                    //&& neighbor_process_id != proc_id )
                                                    {

#if 0
                                                        Debug() << "element id " << elem.id()
                                                                << ", element neighbor id " << neighbor_id
                                                                << " in proc " << neighbor_process_id << "\n";
#endif
                                                        neighbor = boost::addressof( _M_X1->mesh()->element( neighbor_id,
                                                                                                             neighbor_process_id ) );

#if 0
                                                        Debug() << "found neighbor of element id " << elem.id()
                                                                << ", element neighbor id " << neighbor->id()
                                                                << " in proc " << neighbor->processId() << "\n";
#endif

                                                        if ( neighbor_id == neighbor->id()  )
                                                            {
                                                                neighbor_dof = _M_X2->dof()->getIndices( neighbor->id() );

                                                                const size_type n_dof_on_neighbor = neighbor_dof.size();
#if 0
                                                                for (size_type j=0; j<n_dof_on_neighbor; j++)
                                                                    {
                                                                        Debug( 5051 ) << "neighbor elem id: " << neighbor->id() << " dof " << neighbor_dof[j] << "\n";
                                                                    }
                                                                Debug( 5051 ) << "looking for dof " << ig1  << "\n";
#endif
#if 0
                                                                std::pair<std::vector<size_type>::iterator,
                                                                    std::vector<size_type>::iterator>
                                                                    posig = std::equal_range (neighbor_dof.begin(), neighbor_dof.end(), ig1 );
#else
                                                                std::vector<size_type>::iterator posig = std::find( neighbor_dof.begin(), neighbor_dof.end(), ig1 );
#endif
                                                                // Insert jg if it wasn't found
                                                                //if (posig.first != posig.second)
                                                                if ( posig != neighbor_dof.end() ||
                                                                     graph.test ( Pattern::EXTENDED ) )

                                                                    {
                                                                        //Debug() << "found element in proc " << neighbor_process_id << " that shares dof\n";
                                                                        for (size_type j=0; j<n_dof_on_neighbor; j++)
                                                                            {
                                                                                const size_type jg = neighbor_dof[j];

#if 0
                                                                                // See if jg is in the sorted range
                                                                                std::pair<std::vector<size_type>::iterator,
                                                                                    std::vector<size_type>::iterator>
                                                                                    pos = std::equal_range (row.get<2>().begin(), row.get<2>().end(), jg);
#else
                                                                                std::vector<size_type>::iterator pos = std::find( row.get<2>().begin(), row.get<2>().end(), jg );

#endif
                                                                                // Insert jg if it wasn't found
                                                                                if (pos == row.get<2>().end() )
                                                                                    {
                                                                                        const size_type old_size = row.get<2>().size();
                                                                                        row.get<2>().push_back (jg);
                                                                                        //std::inplace_merge (row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end());
                                                                                        sortSparsityRow (row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end());
                                                                                    }
                                                                            }
                                                                    }
                                                            }
                                                    }

                                            } // neighbor graph
#if 0
                                        for( int k = 0; k < row.get<2>().size(); ++k )
                                            Debug() << "row[ " << ig1 - first1_dof_on_proc << ","<< k << " ]=" << row.get<2>()[k] << "\n";
#endif
                                } // only dof on proc

                        }// dof loop
                } // element iterator loop
        }
    else
        {}

    Debug( 5050 ) << "[computeGraph<true>] before calling close in " << t.elapsed() << "s\n";
    sparsity_graph->close();
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
    const size_type nprocs           = _M_X1->mesh()->comm().size();
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

    auto elem_it  = _M_X1->mesh()->beginElementWithProcessId( proc_id );
    auto elem_en  = _M_X1->mesh()->endElementWithProcessId( proc_id );

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
        element_dof2(_M_X2->dof()->getIndicesSize()),
        neighbor_dof;

    for ( ; elem_it != elem_en; ++elem_it)
    {
#if !defined(NDEBUG)
        Debug( 5050 ) << "[Stencil::computePatter] element " << elem_it->id() << " on proc " << elem_it->processId() << "\n";
#endif /* NDEBUG */
        const auto & elem = *elem_it;

        // Get the global indices of the DOFs with support on this element
        //element_dof1 = _M_X1->dof()->getIndices( elem.id() );
        _M_X2->dof()->getIndicesSet( elem.id(), element_dof2 );

        // We can be more efficient if we sort the element DOFs
        // into increasing order
        //std::sort(element_dof1.begin(), element_dof1.end());
        std::sort(element_dof2.begin(), element_dof2.end());

        //const uint16_type  n1_dof_on_element = element_dof1.size();
        const uint16_type  n1_dof_on_element = _M_X1->dof()->getIndicesSize();
        const uint16_type  n2_dof_on_element = element_dof2.size();

        for (size_type i=0; i<n1_dof_on_element; i++)
        //BOOST_FOREACH( auto ig1, _M_X1->dof()->getIndices( elem.id() ) )
        {
            const size_type ig1 = _M_X1->dof()->localToGlobalId( elem.id(), i );
            //const size_type ig1 = element_dof1[i];
            const int ndofpercomponent1 = n1_dof_on_element / _M_X1->dof()->nComponents;
            const int ncomp1 = i / ndofpercomponent1;
            const int ndofpercomponent2 = n2_dof_on_element / _M_X2->dof()->nComponents;

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
                FEEL_ASSERT (ig1 >= first1_dof_on_proc )( ig1 )( first1_dof_on_proc ).error ("invalid dof index");
                FEEL_ASSERT ((ig1 - first1_dof_on_proc) < sparsity_graph->size() )
                    ( ig1 )( first1_dof_on_proc )( sparsity_graph->size() ).error( "invalid dof index" );
#endif
                graph_type::row_type& row = sparsity_graph->row(ig1);
                bool is_on_proc = ( ig1 >= first1_dof_on_proc) && (ig1 <= last1_dof_on_proc);
                row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                Debug( 5051 ) << "work with row " << ig1 << " local index " << ig1 - first1_dof_on_proc << "\n";

                if ( do_less )
                {
                    if ( ncomp1 == (_M_X2->dof()->nComponents-1) )
                        row.get<2>().insert( element_dof2.begin()+ncomp1*ndofpercomponent2,
                                             element_dof2.end() );
                    else
                        row.get<2>().insert( element_dof2.begin()+ncomp1*ndofpercomponent2,
                                             element_dof2.begin()+(ncomp1+1)*ndofpercomponent2 );
                }
                else
                {
                    row.get<2>().insert( element_dof2.begin(), element_dof2.end() );
                }
                // Now (possibly) add dof from neighboring elements
                if ( graph.test( Pattern::EXTENDED ) )
                {
                    for (uint16_type ms=0; ms < elem.nNeighbors(); ms++)
                    {
                        const auto * neighbor = boost::addressof( elem );
                        size_type neighbor_id = elem.neighbor(ms).first;
                        size_type neighbor_process_id = elem.neighbor(ms).second;
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
                                    if ( ncomp1 == (_M_X2->dof()->nComponents-1) )
                                        row.get<2>().insert( neighbor_dof.begin()+ncomp1*ndofpercomponent2,
                                                             neighbor_dof.end() );
                                    else
                                        row.get<2>().insert( neighbor_dof.begin()+ncomp1*ndofpercomponent2,
                                                             neighbor_dof.begin()+(ncomp1+1)*ndofpercomponent2 );

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
    sparsity_graph->close();
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    return sparsity_graph;
}
#endif




template<typename X1,  typename X2>
typename Stencil<X1,X2>::graph_ptrtype
Stencil<X1,X2>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true> )
{

    const size_type nprocs           = _M_X1->mesh()->comm().size();
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

    auto elem_it  = _M_X1->mesh()->beginElementWithProcessId( proc_id );
    auto elem_en  = _M_X1->mesh()->endElementWithProcessId( proc_id );

    auto locTool = _M_X2->mesh()->tool_localization();
    locTool->updateForUse();
    locTool->setExtrapolation(false);
    //locTool->kdtree()->nbNearNeighbor(_M_X2->mesh()->numElements() );

    std::vector<size_type>
        element_dof1,
        element_dof2;

    for ( ; elem_it != elem_en; ++elem_it)
    {
        const auto& elem = *elem_it;

        // Get the global indices of the DOFs with support on this element
        element_dof1 = _M_X1->dof()->getIndices( elem.id() );

        const uint16_type  n1_dof_on_element = element_dof1.size();
        //const uint16_type  n2_dof_on_element = element_dof2.size();

        for (size_type i=0; i<n1_dof_on_element; i++)
            {
                const size_type ig1 = element_dof1[i];
                const int ndofpercomponent1 = n1_dof_on_element / _M_X1->dof()->nComponents;
                const int ncomp1 = i / ndofpercomponent1;
                //const int ndofpercomponent2 = n2_dof_on_element / _M_X2->dof()->nComponents;

                auto ptRealDof = boost::get<0>(_M_X1->dof()->dofPoint(ig1));

                //std::cout << "Pt dof " << ptRealDof<<std::endl;

                auto res = locTool->searchElements(ptRealDof);
                auto hasFind = res.get<0>();

                if ( hasFind )
                    {
                        //std::cout << "\n  find"<<std::endl;

                        auto res_it = res.get<1>().begin();
                        auto res_en = res.get<1>().end();
                        for ( ; res_it != res_en ; ++res_it)
                            {
                                element_dof2 = _M_X2->dof()->getIndices( res_it->get<0>() );

                                graph_type::row_type& row = sparsity_graph->row(ig1);
                                bool is_on_proc = ( ig1 >= first1_dof_on_proc) && (ig1 <= last1_dof_on_proc);
                                row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                                row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;

                                //if ( do_less ) {}
                                //else
                                row.get<2>().insert( element_dof2.begin(), element_dof2.end() );


                                auto elem = _M_X2->mesh()->element(res_it->get<0>());
                                if ( /*graph.test( DOF_PATTERN_NEIGHBOR )*/false )
                                    {
                                        for (uint16_type ms=0; ms < elem.nNeighbors(); ms++)
                                            {
                                                const auto* neighbor = boost::addressof( elem );
                                                size_type neighbor_id = elem.neighbor(ms).first;
                                                size_type neighbor_process_id = elem.neighbor(ms).second;
                                                if ( neighbor_id != invalid_size_type_value )
                                                    //&& neighbor_process_id != proc_id )
                                                    {
                                                        neighbor = boost::addressof( _M_X2->mesh()->element( neighbor_id,
                                                                                                             neighbor_process_id ) );
                                                        if ( neighbor_id == neighbor->id()  )
                                                            {
                                                                auto neighbor_dof = _M_X2->dof()->getIndices( neighbor->id() );
#if 0
                                                                if ( do_less )
                                                                    {
                                                                        if ( ncomp1 == (_M_X2->dof()->nComponents-1) )
                                                                            row.get<2>().insert( neighbor_dof.begin()+ncomp1*ndofpercomponent2,
                                                                                                 neighbor_dof.end() );
                                                                        else
                                                                            row.get<2>().insert( neighbor_dof.begin()+ncomp1*ndofpercomponent2,
                                                                                                 neighbor_dof.begin()+(ncomp1+1)*ndofpercomponent2 );
                                                                    }
                                                                else
                                                                    {
#endif
                                                                        row.get<2>().insert( neighbor_dof.begin(), neighbor_dof.end() );
                                                                        //}

                                                            } // neighbor_id
                                                    } // neighbor_id

                                            } // ms
                                    } // true


                            } // res

                    } // if (hasFind)
                else
                    {
                        //std::cout << "\n not find"<<std::endl;
                        // row empty
                        graph_type::row_type& row = sparsity_graph->row(ig1);
                        bool is_on_proc = ( ig1 >= first1_dof_on_proc) && (ig1 <= last1_dof_on_proc);
                        row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                        row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                        row.get<2>().clear();
                    }

            } // for (size_type i=0; i<n1_dof_on_element; i++)



    } // for ( ; elem_it ... )

    locTool->setExtrapolation(true);
    //locTool->kdtree()->nbNearNeighbor( 15 );

    sparsity_graph->close();

    return sparsity_graph;
}


}
#endif
#endif /* __galerkingraph_H */
