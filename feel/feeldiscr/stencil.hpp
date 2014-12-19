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
        static const uint16_type tag1 = Space1Type::basis_type::TAG;
        static const uint16_type tag2 = Space2::basis_type::TAG;
        typedef mpl::bool_<BFType::template rangeiteratorType<tag2,tag1>::hasnotfindrange_type::value> hasnotfindrange_type;
        typedef mpl::bool_<BFType::template rangeExtendedIteratorType<tag2,tag1>::hasnotfindrange_type::value> hasnotfindrange_extended_type;

        if ( M_stencil->testSpace()->worldsComm()[M_test_index].isActive() )
        {
            if ( M_stencil->isBlockPatternZero( M_test_index,M_trial_index ) )
            {
                typename BFType::graph_ptrtype zerograph( new typename BFType::graph_type( space2->dof(), M_space1->dof() ) );
                zerograph->zero();
                M_stencil->mergeGraph( M_stencil->testSpace()->nDofStart( M_test_index ), M_stencil->trialSpace()->nDofStart( M_trial_index ) , zerograph );
            }

            else
            {
                auto thestencil = stencil( _test=space2, _trial=M_space1,
                                           _pattern=M_stencil->blockPattern( M_test_index,M_trial_index ),
                                           _pattern_block=M_stencil->blockPattern(),
                                           _diag_is_nonzero=false,
                                           _collect_garbage=false,
                                           _close=false,
                                           _range=M_stencil->template subRangeIterator<tag2,tag1>( hasnotfindrange_type() ),
                                           _range_extended=M_stencil->template subRangeExtendedIterator<tag2,tag1>( hasnotfindrange_extended_type() ),
                                           _quad=typename BFType::nonstandard_quadset_type()
                                           );

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
        static const uint16_type tag1 = Space1Type::basis_type::TAG;
        static const uint16_type tag2 = Space2::basis_type::TAG;
        typedef mpl::bool_<BFType::template rangeiteratorType<tag1,tag2>::hasnotfindrange_type::value> hasnotfindrange_type;
        typedef mpl::bool_<BFType::template rangeExtendedIteratorType<tag1,tag2>::hasnotfindrange_type::value> hasnotfindrange_extended_type;

        if ( M_stencil->testSpace()->worldsComm()[M_test_index].isActive() )
        {
            if ( M_stencil->isBlockPatternZero( M_test_index,M_trial_index ) )
            {
                typename BFType::graph_ptrtype zerograph( new typename BFType::graph_type(  M_space1->dof(), space2->dof() ) );
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
                                           _collect_garbage=false,
                                           _close=false,
                                           _range=M_stencil->template subRangeIterator<tag1,tag2>( hasnotfindrange_type() ),
                                           _range_extended=M_stencil->template subRangeExtendedIterator<tag2,tag1>( hasnotfindrange_extended_type() ),
                                           _quad=typename BFType::nonstandard_quadset_type()
                                           );

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



template <typename RangeType>
struct stencilrangetype
{
    typedef typename mpl::if_< boost::is_std_list<RangeType>,
                               mpl::identity<RangeType>,
                               mpl::identity<std::list<RangeType> > >::type::type::value_type type;
    typedef std::list<type> list_type;
};

template <int I,int J,typename IteratorRange>
fusion::pair< fusion::pair<mpl::int_<I>,mpl::int_<J> >,typename stencilrangetype<IteratorRange>::list_type >
stencilRange( IteratorRange const& r, mpl::true_ /**/ )
{
    return fusion::make_pair< fusion::pair<mpl::int_<I>,mpl::int_<J> > >( r);
}
template <int I,int J,typename IteratorRange>
fusion::pair< fusion::pair<mpl::int_<I>,mpl::int_<J> >,typename stencilrangetype<IteratorRange>::list_type >
stencilRange( IteratorRange const& r, mpl::false_ /**/)
{
    typename stencilrangetype<IteratorRange>::list_type res;
    res.push_back( r );
    return fusion::make_pair< fusion::pair<mpl::int_<I>,mpl::int_<J> > >( res );
}
template <int I,int J,typename IteratorRange>
fusion::pair< fusion::pair<mpl::int_<I>,mpl::int_<J> >,typename stencilrangetype<IteratorRange>::list_type >
stencilRange( IteratorRange const& r)
{
    return stencilRange<I,J,IteratorRange>(r, mpl::bool_< boost::is_std_list<IteratorRange>::value >() );
}



struct stencilRangeMapTypeBase {};

struct stencilRangeMap0Type
    :
    public stencilRangeMapTypeBase,
    public fusion::map<>
{
    typedef fusion::map<> super_type;

    stencilRangeMap0Type()
        :
        super_type()
    {}

    static bool isNullRange() { return true; }
};

template <typename ThePair1Type>
struct stencilRangeMap1Type
    :
    public stencilRangeMapTypeBase,
    public fusion::map< fusion::pair< typename ThePair1Type::first_type, typename ThePair1Type::second_type > >
{
    typedef typename ThePair1Type::first_type key1_type;
    typedef fusion::map< fusion::pair< key1_type, typename ThePair1Type::second_type > > super_type;

    stencilRangeMap1Type( ThePair1Type const& p )
        :
        super_type( p )
    {}

    static bool isNullRange() { return false; }
};

template <typename ThePair1Type,typename ThePair2Type>
struct stencilRangeMap2Type
    :
    public stencilRangeMapTypeBase,
    public fusion::map< fusion::pair< typename ThePair1Type::first_type, typename ThePair1Type::second_type >,
                        fusion::pair< typename ThePair2Type::first_type, typename ThePair2Type::second_type > >
{
    typedef typename ThePair1Type::first_type key1_type;
    typedef typename ThePair2Type::first_type key2_type;
    typedef fusion::map< fusion::pair< key1_type, typename ThePair1Type::second_type >,
                         fusion::pair< key2_type, typename ThePair2Type::second_type > > super_type;

    stencilRangeMap2Type( ThePair1Type const& p1, ThePair2Type const& p2 )
        :
        super_type( p1,p2 )
    {}

    static bool isNullRange() { return false; }
};

stencilRangeMap0Type stencilRangeMap();

template <typename ThePair1Type>
stencilRangeMap1Type<ThePair1Type>
stencilRangeMap( ThePair1Type const& p)
{
    return stencilRangeMap1Type<ThePair1Type>( p.second );
}

template <typename ThePair1Type,typename ThePair2Type>
stencilRangeMap2Type<ThePair1Type,ThePair2Type>
stencilRangeMap( ThePair1Type const& p1, ThePair2Type const& p2)
{
    return stencilRangeMap2Type<ThePair1Type,ThePair2Type>( p1.second, p2.second );
}

/**
 * define the quadrature order use with non standard stencil
 */
struct stencilQuadSetBase {};

template <int QuadOrder1d=20, int QuadOrder2d=12, int QuadOrder3d=8 >
struct stencilQuadSet :
        public stencilQuadSetBase,
        public fusion::vector< _Q<QuadOrder1d>,_Q<QuadOrder2d>,_Q<QuadOrder3d> >
{
    typedef fusion::vector< _Q<QuadOrder1d>,_Q<QuadOrder2d>,_Q<QuadOrder3d> > super_type;

    stencilQuadSet()
        :
        super_type( _Q<QuadOrder1d>(), _Q<QuadOrder2d>(), _Q<QuadOrder3d>() )
    {}
};


template<typename X1, typename X2,
         typename RangeIteratorTestType = stencilRangeMap0Type,
         typename RangeExtendedIteratorType = stencilRangeMap0Type,
         typename QuadSetType = stencilQuadSet<> >
class Stencil
{
public:
    typedef X1 test_space_ptrtype;
    typedef X2 trial_space_ptrtype;
    typedef typename X1::element_type test_space_type;
    typedef typename X2::element_type trial_space_type;
    typedef GraphCSR graph_type;
    typedef boost::shared_ptr<graph_type> graph_ptrtype;
    typedef Stencil<X1,X2,RangeIteratorTestType,RangeExtendedIteratorType,QuadSetType> self_type;

    typedef RangeIteratorTestType rangeiterator_test_type;
    typedef RangeExtendedIteratorType rangeiterator_extended_type;
    typedef QuadSetType nonstandard_quadset_type;

    Stencil( test_space_ptrtype Xh, trial_space_ptrtype Yh,
             size_type graph_hints,
             BlocksStencilPattern block_pattern=BlocksStencilPattern(1,1,Pattern::HAS_NO_BLOCK_PATTERN),
             bool diag_is_nonzero=false,
             bool close=true,
             rangeiterator_test_type r=rangeiterator_test_type(),
             rangeiterator_extended_type rangeExtended=rangeiterator_extended_type() )
        :
        _M_X1( Xh ),
        _M_X2( Yh ),
        M_graph( new graph_type( Xh->dof(),Yh->dof() ) ),
        M_block_pattern( block_pattern ),
        M_rangeIteratorTest( r ),
        M_rangeIteratorExtended( rangeExtended )
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

        if ( close ) M_graph->close();
    }

    Stencil( test_space_ptrtype Xh, trial_space_ptrtype Yh, size_type graph_hints, graph_ptrtype g,
             rangeiterator_test_type r=rangeiterator_test_type(),
             rangeiterator_extended_type rangeExtended=rangeiterator_extended_type() )
        :
        _M_X1( Xh ),
        _M_X2( Yh ),
        M_graph( g ),
        M_block_pattern(Xh->nSubFunctionSpace(),Yh->nSubFunctionSpace(),size_type( graph_hints/*Pattern::HAS_NO_BLOCK_PATTERN*/ )),
        M_rangeIteratorTest( r ),
        M_rangeIteratorExtended( rangeExtended )
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

    template<typename EltType>
    std::set<std::pair<size_type,rank_type> >
    testElementIdFromRange( mpl::size_t<MESH_ELEMENTS> /**/, EltType const& elem )
    {
        const bool test_related_to_range = _M_X1->mesh()->isSubMeshFrom( elem.mesh() );
        const bool range_related_to_test = elem.mesh()->isSubMeshFrom( _M_X1->mesh() );
        std::set<std::pair<size_type,rank_type> > res;
        if ( test_related_to_range )
        {
            const size_type test_eid = _M_X1->mesh()->meshToSubMesh( elem.id() );
            if ( test_eid != invalid_size_type_value )
                res.insert( std::make_pair( test_eid,elem.processId() ) );
        }
        else if ( range_related_to_test )
        {
            const size_type test_eid = elem.mesh()->subMeshToMesh( elem.id() );
            if ( test_eid != invalid_size_type_value )
                res.insert( std::make_pair( test_eid,elem.processId() ) );
        }
        else // same mesh
        {
            res.insert( std::make_pair( elem.id(),elem.processId() ) );
        }
        return res;
    }

    template<typename FaceType>
    std::set<std::pair<size_type,rank_type> >
    testElementIdFromRange( mpl::size_t<MESH_FACES> /**/, FaceType const& theface )
    {
        std::set<std::pair<size_type,rank_type> > res;
        if ( theface.isConnectedTo0() && !theface.element0().isGhostCell() )
        {
            auto resElt0 = testElementIdFromRange( mpl::size_t<MESH_ELEMENTS>(), theface.element( 0 ) );
            for ( std::pair<size_type,rank_type> const& idElt0 : resElt0 )
                res.insert( idElt0 );
        }
        if ( theface.isConnectedTo1() && !theface.element1().isGhostCell() )
        {
            auto resElt1 = testElementIdFromRange( mpl::size_t<MESH_ELEMENTS>(), theface.element( 1 ) );
            for ( std::pair<size_type,rank_type> const& idElt1 : resElt1 )
                res.insert( idElt1 );
        }
        return res;
    }

    std::set<size_type> trialElementId( size_type test_eid, mpl::int_<0> /**/ )
        {
            std::set<size_type> idsFind;
            const bool test_related_to_trial = _M_X1->mesh()->isSubMeshFrom( _M_X2->mesh() );
            const bool trial_related_to_test = _M_X2->mesh()->isSubMeshFrom( _M_X1->mesh() );
            if ( test_related_to_trial )
            {
                const size_type domain_eid = _M_X1->mesh()->subMeshToMesh( test_eid );
                DVLOG(2) << "[test_related_to_trial] test element id: "  << test_eid << " trial element id : " << domain_eid << "\n";
                if ( domain_eid != invalid_size_type_value ) idsFind.insert( domain_eid );
            }
            else if( trial_related_to_test )
            {
                const size_type domain_eid = _M_X2->mesh()->meshToSubMesh( test_eid );
                DVLOG(2) << "[trial_related_to_test] test element id: "  << test_eid << " trial element id : " << domain_eid << "\n";
                if ( domain_eid != invalid_size_type_value ) idsFind.insert( domain_eid );
            }
            else // same mesh
            {
                idsFind.insert(test_eid);
            }

            return idsFind;
        }
    std::set<size_type> trialElementId( size_type test_eid, mpl::int_<1> /**/ )
        {
            std::set<size_type> idsFind;
            const bool test_related_to_trial = _M_X1->mesh()->isSubMeshFrom( _M_X2->mesh() );
            const bool trial_related_to_test = _M_X2->mesh()->isSubMeshFrom( _M_X1->mesh() );
            const bool trial_sibling_of_test = _M_X2->mesh()->isSiblingOf( _M_X1->mesh() );
            if ( test_related_to_trial )
            {
                auto const& theface = _M_X2->mesh()->face(_M_X1->mesh()->subMeshToMesh( test_eid ));
                size_type domain_eid = invalid_size_type_value;
                if ( !theface.element0().isGhostCell() )
                    domain_eid = theface.element0().id();
                else if ( theface.isConnectedTo1() && !theface.element1().isGhostCell() )
                    domain_eid = theface.element1().id();
                else
                    CHECK(false) << " error : maybe the faces is not on partition or invalid connection\n";

                DVLOG(2) << "[test_related_to_trial<1>] test element id: "  << test_eid << " trial element id : " << domain_eid << "\n";
                if ( domain_eid != invalid_size_type_value ) idsFind.insert( domain_eid );
            }
            else if( trial_related_to_test )
            {
                auto const& eltTest = _M_X1->mesh()->element(test_eid);
                for (uint16_type f=0;f< _M_X1->mesh()->numLocalFaces();++f)
                    {
                        const size_type idFind = _M_X2->mesh()->meshToSubMesh( eltTest.face(f).id() );
                        if ( idFind != invalid_size_type_value ) idsFind.insert( idFind );
                    }
                DVLOG(2) << "[trial_related_to_test<1>] test element id: "  << test_eid << " idsFind.size() "<< idsFind.size() << "\n";
            }
            else if ( trial_sibling_of_test )
            {
                DVLOG(1) << "test_eid = " << test_eid;
                size_type id_in_sibling = _M_X2->mesh()->meshToSubMesh( _M_X1->mesh(), test_eid );
                DVLOG(1) << "id_in_sibling = " << id_in_sibling;
                size_type domain_eid = invalid_size_type_value;
                if ( id_in_sibling!=invalid_size_type_value)
                {
                    domain_eid = _M_X2->mesh()->face(id_in_sibling).element0().id();
                    DVLOG(1) << "[test_sibling_of_trial<1>] test element id: "  << test_eid << " trial element id : " << domain_eid << "\n";
                    idsFind.insert( domain_eid );
                }
            }
            else
            {
                CHECK ( false ) << "[trial_related_to_test<1>] : test and trial mesh can not be the same here\n";
            }
            return idsFind;
        }
    std::set<size_type> trialElementId( size_type test_eid, mpl::int_<2> /**/ )
    {
        CHECK ( false ) << "[trial_related_to_test<2>] : submesh relation with codim=2 is not implement\n";
        return std::set<size_type>();
    }


public :
    template <int I,int J>
    struct rangeiteratorType
    {
        typedef typename fusion::result_of::find<rangeiterator_test_type,fusion::pair<mpl::int_<I>,mpl::int_<J> > >::type resultfindrange_it_type;
        typedef typename boost::is_same<resultfindrange_it_type, typename fusion::result_of::end<rangeiterator_test_type>::type> hasnotfindrange_type;

        typedef typename boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                                      typename MeshTraits<typename test_space_type::mesh_type>::element_const_iterator,
                                      typename MeshTraits<typename test_space_type::mesh_type>::element_const_iterator> defaultrange_type;

        // fix compilation from boost 1.55
        // if not find in fusion map else there is a problem now with result_of::value_of
        // so we create an bid iterator to fix that
        typedef fusion::map< fusion::pair< fusion::pair<mpl::int_<-1>,mpl::int_<-1> >, defaultrange_type > > mapbid;
        typedef typename fusion::result_of::find< mapbid,fusion::pair<mpl::int_<-1>,mpl::int_<-1> > >::type resultfindrange_it_bid_type;
        typedef typename mpl::if_< hasnotfindrange_type,
                                   resultfindrange_it_bid_type,
                                   resultfindrange_it_type >::type resultfindrange_it2_type;
        typedef typename fusion::result_of::value_of<resultfindrange_it2_type>::type resultfindrange_type;

        typedef typename mpl::if_< hasnotfindrange_type,
                                   mpl::identity< defaultrange_type >,
                                   mpl::identity< resultfindrange_type >
                                  >::type::type type;
    };

    template <int I,int J>
    struct rangeExtendedIteratorType
    {
        typedef typename fusion::result_of::find<rangeiterator_extended_type,fusion::pair<mpl::int_<I>,mpl::int_<J> > >::type resultfindrange_it_type;
        typedef typename boost::is_same<resultfindrange_it_type, typename fusion::result_of::end<rangeiterator_extended_type>::type> hasnotfindrange_type;

        typedef typename boost::tuple<mpl::size_t<MESH_FACES>,
                                      typename MeshTraits<typename test_space_type::mesh_type>::location_face_const_iterator,
                                      typename MeshTraits<typename test_space_type::mesh_type>::location_face_const_iterator> defaultrange_type;

        // fix compilation from boost 1.55
        // if not find in fusion map else there is a problem now with result_of::value_of
        // so we create an bid iterator to fix that
        typedef fusion::map< fusion::pair< fusion::pair<mpl::int_<-1>,mpl::int_<-1> >, defaultrange_type > > mapbid;
        typedef typename fusion::result_of::find< mapbid,fusion::pair<mpl::int_<-1>,mpl::int_<-1> > >::type resultfindrange_it_bid_type;
        typedef typename mpl::if_< hasnotfindrange_type,
                                   resultfindrange_it_bid_type,
                                   resultfindrange_it_type >::type resultfindrange_it2_type;
        typedef typename fusion::result_of::value_of<resultfindrange_it2_type>::type resultfindrange_type;

        typedef typename mpl::if_< hasnotfindrange_type,
                                   mpl::identity< defaultrange_type >,
                                   mpl::identity< resultfindrange_type >
                                  >::type::type type;
    };

    /**
     * range/sub-range for standard stencil
     */
#if 0
    template <int I,int J>
    typename rangeiteratorType<I,J>::type
    rangeiterator() const
    {
        return rangeiterator<I,J>( mpl::bool_<rangeiteratorType<I,J>::hasnotfindrange_type::value>() );
    }
#endif
    template <int I,int J>
    std::list< typename  rangeiteratorType<I,J>::defaultrange_type>
    rangeiterator(mpl::bool_<true> /**/) const
    {
        std::list<typename rangeiteratorType<I,J>::defaultrange_type> res;
        res.push_back( elements( _M_X1->mesh() ) );
        return res;
    }
    template <int I,int J>
    typename rangeiteratorType<I,J>::resultfindrange_type::second_type
    rangeiterator(mpl::bool_<false> /**/) const
    {
        typedef fusion::pair<mpl::int_<I>,mpl::int_<J> > key_type;
        return fusion::at_key< key_type >( M_rangeIteratorTest );
    }
    template <int I,int J>
    stencilRangeMap0Type
    subRangeIterator( mpl::bool_<true> /**/ )
    {
        return stencilRangeMap0Type();
    }
    template <int I,int J>
    stencilRangeMap1Type< fusion::pair<  fusion::pair<mpl::int_<0>,mpl::int_<0> >, typename rangeiteratorType<I,J>::resultfindrange_type::second_type > >
    subRangeIterator( mpl::bool_<false> /**/ )
    {
        return stencilRangeMap( stencilRange<0,0>( rangeiterator<I,J>(mpl::bool_<false>()) ) );
    }

    /**
     * range/sub-range for extended stencil
     */
    template <int I,int J>
    std::list<typename rangeExtendedIteratorType<I,J>::defaultrange_type>
    rangeExtendedIterator(mpl::bool_<true> /**/) const
    {
        std::list<typename rangeExtendedIteratorType<I,J>::defaultrange_type> res;
        res.push_back( internalfaces( _M_X1->mesh() ) );
        return res;
    }
    template <int I,int J>
    typename rangeExtendedIteratorType<I,J>::resultfindrange_type::second_type
    rangeExtendedIterator(mpl::bool_<false> /**/) const
    {
        typedef fusion::pair<mpl::int_<I>,mpl::int_<J> > key_type;
        return fusion::at_key< key_type >( M_rangeIteratorExtended );
    }
    template <int I,int J>
    stencilRangeMap0Type
    subRangeExtendedIterator( mpl::bool_<true> /**/ )
    {
        return stencilRangeMap0Type();
    }
    template <int I,int J>
    stencilRangeMap1Type< fusion::pair< fusion::pair<mpl::int_<0>,mpl::int_<0> >, typename rangeExtendedIteratorType<I,J>::resultfindrange_type::second_type > >
    subRangeExtendedIterator( mpl::bool_<false> /**/ )
    {
        return stencilRangeMap( stencilRange<0,0>( rangeExtendedIterator<I,J>( mpl::bool_<false>() ) ) );
    }

private:

    test_space_ptrtype _M_X1;
    trial_space_ptrtype _M_X2;
    graph_ptrtype M_graph;
    BlocksStencilPattern M_block_pattern;
    rangeiterator_test_type M_rangeIteratorTest;
    rangeiterator_extended_type M_rangeIteratorExtended;
};
namespace detail
{
template<typename Args>
struct compute_stencil_type
{
    typedef typename remove_pointer_const_reference_type<Args,tag::test>::type _test_type;
    typedef typename remove_pointer_const_reference_type<Args,tag::trial>::type _trial_type;
    typedef typename remove_pointer_const_reference_default_type<Args,tag::range, stencilRangeMap0Type >::type _range_type;
    typedef typename remove_pointer_const_reference_default_type<Args,tag::range_extended, stencilRangeMap0Type >::type _range_extended_type;
    typedef typename remove_pointer_const_reference_default_type<Args,tag::quad, stencilQuadSet<> >::type _quad_type;
    typedef Stencil<_test_type, _trial_type, _range_type, _range_extended_type, _quad_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};

}

class StencilManagerImpl:
    public std::map<boost::tuple<boost::weak_ptr<FunctionSpaceBase>,
    boost::weak_ptr<FunctionSpaceBase>,
    size_type,
    std::vector<size_type>,
    bool >, boost::shared_ptr<GraphCSR> >,
public boost::noncopyable
{
public:
    typedef boost::shared_ptr<GraphCSR> graph_ptrtype;
    typedef boost::tuple<boost::weak_ptr<FunctionSpaceBase>,
            boost::weak_ptr<FunctionSpaceBase>,
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
    ( typename Feel::detail::compute_stencil_type<Args>::ptrtype ), // 1. return type
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
      ( collect_garbage,  *( boost::is_integral<mpl::_> ), true )
      ( close,            *( boost::is_integral<mpl::_> ), true )
      ( range,            *( boost::is_convertible<mpl::_, stencilRangeMapTypeBase>) , stencilRangeMap0Type() )
      ( range_extended,   *( boost::is_convertible<mpl::_, stencilRangeMapTypeBase>) , stencilRangeMap0Type() )
      ( quad,             *( boost::is_convertible<mpl::_, stencilQuadSetBase>), stencilQuadSet<>() )
    )
)
{
    if ( collect_garbage )
    {
        // cleanup memory before doing anything
        stencilManagerGarbageCollect();
    }

    Feel::detail::ignore_unused_variable_warning( args );
    typedef typename Feel::detail::compute_stencil_type<Args>::ptrtype stencil_ptrtype;
    typedef typename Feel::detail::compute_stencil_type<Args>::type stencil_type;

    // we look into the spaces dictionary for existing graph
    auto git = StencilManager::instance().find( boost::make_tuple( test, trial, pattern, pattern_block.getSetOfBlocks(), diag_is_nonzero ) );

    if ( git != StencilManager::instance().end() && range.isNullRange() && range_extended.isNullRange() )
    {
        //std::cout << "Found a  stencil in manager (" << test.get() << "," << trial.get() << "," << pattern << ")\n";
        auto s = stencil_ptrtype( new stencil_type( test, trial, pattern, git->second, range, range_extended ) );
        return s;
    }

    else
    {
        // look for transposed stencil if it exist and transpose it to get the stencil
        auto git_trans = StencilManager::instance().find( boost::make_tuple( trial, test, pattern, pattern_block.transpose().getSetOfBlocks(), diag_is_nonzero ) );

        if ( git_trans != StencilManager::instance().end() && range.isNullRange() && range_extended.isNullRange() )
        {
            auto g = git_trans->second->transpose(close);
            //auto g = git_trans->second->transpose();
            stencilManagerAdd( boost::make_tuple( test, trial, pattern, pattern_block.getSetOfBlocks(), diag_is_nonzero ), g );

            auto s = stencil_ptrtype( new stencil_type( test, trial, pattern, g, range, range_extended ) );
            //std::cout << "Found a  transposed stencil in manager (" << test.get() << "," << trial.get() << "," << pattern << ")\n";
            return s;
        }

        else
        {
            //std::cout << "Creating a new stencil in manager (" << test.get() << "," << trial.get() << "," << pattern << ")\n";
            auto s = stencil_ptrtype( new stencil_type( test, trial, pattern, pattern_block, diag_is_nonzero, close, range, range_extended ) );
            if ( range.isNullRange() && range_extended.isNullRange() )
                stencilManagerAdd( boost::make_tuple( test, trial, pattern, pattern_block.getSetOfBlocks(), diag_is_nonzero ), s->graph() );
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
template<typename X1,  typename X2, typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
void
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::mergeGraph( int row, int col, graph_ptrtype g )
{
    boost::timer tim;
    DVLOG(2) << "[merge graph] for composite bilinear form\n";
    DVLOG(2) << "[mergeGraph] row = " << row << "\n";
    DVLOG(2) << "[mergeGraph] col = " << col << "\n";

    // nothing yet in store
    //if ( !M_graph || M_graph->empty() )
    if ( 0 )
    {
        DVLOG(2) << "[merge graph] nothing yet in store, copy graph\n";
        M_graph = g;

#if 0
        typename graph_type::const_iterator it = g->begin();
        typename graph_type::const_iterator en = g->end();

        for ( ; it != en; ++it )
        {
            std::vector<size_type> const& ivec = boost::get<2>( it->second );

            for ( int i = 0; i < ivec.size(); ++i )
            {
                DVLOG(2) << "[mergeGraph] ivec[" << i << "] = " << ivec[i] << "\n";
            }
        }

#endif // 9
    }

    else
    {
        //std::cout << "\n row " << row << " col " << col << " with god rank" <<  this->testSpace()->worldComm().godRank() << std::endl;
        DVLOG(2) << "[merge graph] already something in store\n";
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

            DVLOG(2) << "[mergeGraph] adding information to global row [" << theglobalrow << "], localrow=" << thelocalrow << "\n";
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
                        for ( auto itcol = row2_entries.begin(), encol = row2_entries.end() ; itcol!=encol; ++itcol ) row1_entries.insert( *itcol+col );
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

    DVLOG(2) << " -- merge_graph (" << row << "," << col << ") in " << tim.elapsed() << "\n";
    DVLOG(2) << "merge graph for composite bilinear form done\n";
}

template<typename X1,  typename X2, typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
void
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::mergeGraphMPI( size_type test_index, size_type trial_index,
                               DataMap const& mapOnTest, DataMap const& mapOnTrial,
                               graph_ptrtype g )
{
    const int myrank = this->testSpace()->worldComm().globalRank();

    const size_type globalDofRowStart = this->testSpace()->dof()->firstDofGlobalCluster()  +  this->testSpace()->nLocalDofWithoutGhostOnProcStart( myrank, test_index );
    const size_type globalDofColStart = this->trialSpace()->dof()->firstDofGlobalCluster()  +  this->trialSpace()->nLocalDofWithoutGhostOnProcStart( myrank, trial_index );
    const size_type locdofStart = this->testSpace()->nLocalDofWithGhostOnProcStart( myrank, test_index );

    typename graph_type::const_iterator it = g->begin();
    typename graph_type::const_iterator en = g->end();
    for ( ; it != en; ++it )
    {
        size_type theglobalrow = globalDofRowStart + ( it->first - mapOnTest.firstDofGlobalCluster() );
        const size_type thelocalrow = locdofStart + it->second.get<1>();

        if (it->second.get<0>()!=g->worldComm().globalRank() )
        {
            const int proc = it->second.get<0>();
            const size_type realrowStart = this->testSpace()->dof()->firstDofGlobalCluster(proc)
                + this->testSpace()->nLocalDofWithoutGhostOnProcStart(proc, test_index );
            theglobalrow = realrowStart+(it->first-mapOnTest.firstDofGlobalCluster(proc));
        }

        std::set<size_type>& row1_entries = M_graph->row( theglobalrow ).template get<2>();
        std::set<size_type> const& row2_entries = boost::get<2>( it->second );

        DVLOG(2) << "[mergeGraph] adding information to global row [" << theglobalrow << "], localrow=" << thelocalrow << "\n";
        M_graph->row( theglobalrow ).template get<1>() = thelocalrow;
        M_graph->row( theglobalrow ).template get<0>() = this->testSpace()->worldComm().mapLocalRankToGlobalRank()[it->second.get<0>()];

        if ( !row2_entries.empty() )
        {
            for ( auto itcol = row2_entries.begin(), encol = row2_entries.end() ; itcol!=encol; ++itcol )
            {
                if (mapOnTrial.dofGlobalClusterIsOnProc(*itcol))
                {
                    const size_type dofcol = globalDofColStart + (*itcol-mapOnTrial.firstDofGlobalCluster());
                    row1_entries.insert( dofcol );
                }
                else
                {
                    const int realproc = mapOnTrial.procOnGlobalCluster(*itcol);
                    const size_type realcolStart = this->trialSpace()->dof()->firstDofGlobalCluster(realproc)
                        + this->trialSpace()->nLocalDofWithoutGhostOnProcStart( realproc, trial_index );
                    const size_type dofcol = realcolStart + (*itcol - mapOnTrial.firstDofGlobalCluster(realproc));
                    row1_entries.insert(dofcol);
                }
            }
        }

    } // for( ; it != en; ++it )

}

template<typename X1,  typename X2, typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraph( size_type hints )
{
    VLOG(2) << "computeGraph: deciding whether the mesh are related to optimize the stencil\n";
    //if ( (is_shared_ptr<typename test_space_type::mesh_ptrtype>::value && is_shared_ptr<typename trial_space_type::mesh_ptrtype>::value ) &&
    //dynamic_cast<void*>( _M_X1->template mesh<0>().get() ) == dynamic_cast<void*>( _M_X2->template mesh<0>().get() ) )
    if ( _M_X1->template mesh<0>()->isRelatedTo( _M_X2->template mesh<0>() ) )
    {
        VLOG(2) << "computeGraph: meshes are related\n";
        return this->computeGraph( hints, mpl::bool_<mpl::and_< mpl::bool_< ( test_space_type::nSpaces == 1 )>,
                                                                mpl::bool_< ( trial_space_type::nSpaces == 1 )> >::type::value >() );
    }
    else
    {
        VLOG(2) << "computeGraph: meshes are not related\n";
        return this->computeGraphInCaseOfInterpolate( hints, mpl::bool_<mpl::and_< mpl::bool_< ( test_space_type::nSpaces == 1 )>,
                                                                                   mpl::bool_< ( trial_space_type::nSpaces == 1 )> >::type::value >() );
    }
}


template<typename X1,  typename X2, typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraph( size_type hints, mpl::bool_<false> )
{
    boost::timer t;
    DVLOG(2) << "compute graph for composite bilinear form with interpolation\n";

    auto graph = computeGraph( hints, mpl::bool_< ( test_space_type::nSpaces > 1 )>(), mpl::bool_< ( trial_space_type::nSpaces > 1 )>() );

    DVLOG(2) << "closing graph for composite bilinear form with interpolation done in " << t.elapsed() << "s\n";
    t.restart();
    //graph->close();
    DVLOG(2) << "compute graph for composite bilinear form done in " << t.elapsed() << "s\n";

    return graph;
}

 template<typename X1,  typename X2,typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
 typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
 Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraph( size_type hints, mpl::bool_<true>, mpl::bool_<true> )
{
    fusion::for_each( _M_X1->functionSpaces(),
                      Feel::detail::compute_graph1<self_type>( this, hints ) );
    return M_graph;
}

template<typename X1,  typename X2,typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraph( size_type hints, mpl::bool_<true>, mpl::bool_<false> )
{
    fusion::for_each( _M_X1->functionSpaces(),
                      Feel::detail::compute_graph3<self_type,trial_space_type>( this, _M_X2, 0, hints ) );
    return M_graph;
}

template<typename X1,  typename X2,typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraph( size_type hints, mpl::bool_<false>, mpl::bool_<true> )
{
    fusion::for_each( _M_X2->functionSpaces(),
                      Feel::detail::compute_graph2<self_type,test_space_type>( this, _M_X1, 0, hints ) );
    return M_graph;
}

template<typename X1,  typename X2,typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false> )
{
    boost::timer t;
    DVLOG(2) << "compute graph for composite bilinear form with interpolation\n";
    auto graph = computeGraphInCaseOfInterpolate( hints,
                                                  mpl::bool_< ( test_space_type::nSpaces > 1 )>(),
                                                  mpl::bool_< ( trial_space_type::nSpaces > 1 )>() );

    DVLOG(2) << "closing graph for composite bilinear form with interpolation done in " << t.elapsed() << "s\n";
    return graph;
}

template<typename X1,  typename X2,typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<true> )
{
    fusion::for_each( _M_X1->functionSpaces(),
                      Feel::detail::compute_graph1<self_type>( this, hints ) );
    return M_graph;
}

template<typename X1,  typename X2,typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<false> )
{
    fusion::for_each( _M_X1->functionSpaces(),
                      Feel::detail::compute_graph3<self_type,trial_space_type>( this, _M_X2, 0, hints ) );
    return M_graph;
}

template<typename X1,  typename X2,typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false>, mpl::bool_<true> )
{
    fusion::for_each( _M_X2->functionSpaces(),
                      Feel::detail::compute_graph2<self_type,test_space_type>( this, _M_X1, 0, hints ) );
    return M_graph;
}




#if 0
template<typename X1,  typename X2,typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraph( size_type hints, mpl::bool_<true> )
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
    DVLOG(2) << "[computeGraph] test : " << ( graph.test ( Pattern::COUPLED ) || graph.test ( Pattern::EXTENDED ) ) << "\n";
    DVLOG(2) << "[computeGraph]  : graph.test ( Pattern::COUPLED )=" <<  graph.test ( Pattern::COUPLED ) << "\n";
    DVLOG(2) << "[computeGraph]  : graph.test ( Pattern::EXTENDED)=" <<  graph.test ( Pattern::EXTENDED ) << "\n";
#if 0

    if ( graph.test ( Pattern::COUPLED ) ||
            graph.test ( Pattern::EXTENDED ) )
#else
    if ( 1 )
#endif
    {
        DVLOG(2) << "[computeGraph] test (Pattern::COUPLED || Pattern::EXTENDED) ok\n";
        std::vector<size_type>
        element_dof1,
        element_dof2,
        neighbor_dof,
        dof_to_add;

        for ( ; elem_it != elem_en; ++elem_it )
        {
#if !defined(NDEBUG)
            DVLOG(2) << "[Stencil::computePatter] element " << elem_it->id() << " on proc " << elem_it->processId() << "\n";
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
                    DVLOG(2) << "work with row " << ig1 << " local index " << ig1 - first1_dof_on_proc << "\n";

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
                                    DVLOG(2) << "neighbor elem id: " << neighbor->id() << " dof " << neighbor_dof[j] << "\n";
                                }

                                DVLOG(2) << "looking for dof " << ig1  << "\n";
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

    DVLOG(2) << "[computeGraph<true>] before calling close in " << t.elapsed() << "s\n";
    //sparsity_graph->close();
    DVLOG(2) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    DVLOG(2) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    return sparsity_graph;
}
#else
template<typename X1,  typename X2,typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraph( size_type hints, mpl::bool_<true> )
{
    static const bool hasNotFindRangeStandard = rangeiteratorType<0,0>::hasnotfindrange_type::value;
    static const bool hasNotFindRangeExtended = rangeExtendedIteratorType<0,0>::hasnotfindrange_type::value;

    boost::timer t;
    // Compute the sparsity structure of the global matrix.  This can be
    // fed into a PetscMatrix to allocate exacly the number of nonzeros
    // necessary to store the matrix.  This algorithm should be linear
    // in the (# of elements)*(# nodes per element)
    const size_type proc_id           = _M_X1->worldsComm()[0].localRank();
    const size_type n1_dof_on_proc    = _M_X1->nLocalDof();
    //const size_type n2_dof_on_proc    = _M_X2->nLocalDof();
    const size_type first1_dof_on_proc = _M_X1->dof()->firstDofGlobalCluster( proc_id );
    const size_type last1_dof_on_proc = _M_X1->dof()->lastDofGlobalCluster( proc_id );
    const size_type first2_dof_on_proc = _M_X2->dof()->firstDofGlobalCluster( proc_id );
    const size_type last2_dof_on_proc = _M_X2->dof()->lastDofGlobalCluster( proc_id );

    graph_ptrtype sparsity_graph( new graph_type( _M_X1->dof(),_M_X2->dof() ) );

    Feel::Context graph( hints );
    if ( graph.test( Pattern::ZERO ) ) { sparsity_graph->zero(); return sparsity_graph; }

    //if (_M_X1->nLocalDofWithoutGhost()==0 && _M_X2->nLocalDofWithoutGhost()==0 ) return sparsity_graph;

    //auto elem_it  = _M_X1->mesh()->beginElementWithProcessId( _M_X1->mesh()->worldComm().localRank() /*proc_id*/ );
    //auto elem_en  = _M_X1->mesh()->endElementWithProcessId( _M_X1->mesh()->worldComm().localRank() /*proc_id*/ );

    // If the user did not explicitly specify the DOF coupling
    // then all the DOFS are coupled to each other.  Furthermore,
    // we can take a shortcut and do this more quickly here.  So
    // we use an if-test.
    DVLOG(2) << "[computeGraph]  : graph.test ( Pattern::DEFAULT )=" <<  graph.test ( Pattern::DEFAULT ) << "\n";
    DVLOG(2) << "[computeGraph]  : graph.test ( Pattern::COUPLED )=" <<  graph.test ( Pattern::COUPLED ) << "\n";
    DVLOG(2) << "[computeGraph]  : graph.test ( Pattern::EXTENDED)=" <<  graph.test ( Pattern::EXTENDED ) << "\n";
    bool do_less =  ( ( graph.test( Pattern::DEFAULT ) &&
                        ( _M_X1->dof()->nComponents ==
                          _M_X2->dof()->nComponents ) ) &&
                      !graph.test( Pattern::COUPLED ) );
    std::vector<size_type>
    element_dof2( _M_X2->dof()->getIndicesSize() ),
                  neighbor_dof;

    static const uint16_type nDimTest = test_space_type::mesh_type::nDim;
    static const uint16_type nDimTrial = trial_space_type::mesh_type::nDim;
    static const uint16_type nDimDiffBetweenTestTrial = ( nDimTest > nDimTrial )? nDimTest-nDimTrial : nDimTrial-nDimTest;

    auto rangeListTest = this->rangeiterator<0,0>( mpl::bool_<hasNotFindRangeStandard>() );

    for ( auto const& rangeTest : rangeListTest )
    {
    auto iDimRange = rangeTest.template get<0>();
    auto elem_it = rangeTest.template get<1>();
    auto elem_en = rangeTest.template get<2>();


    for ( ; elem_it != elem_en; ++elem_it )
    {
        DVLOG(4) << "[Stencil::computePattern] element " << elem_it->id() << " on proc " << elem_it->processId() << "\n";

        const std::set<std::pair<size_type,rank_type> > infoTestElts = testElementIdFromRange( iDimRange,*elem_it );
        for ( std::pair<size_type,rank_type> const& infoTestElt : infoTestElts )
        {
        size_type idTestElt = infoTestElt.first;
        rank_type pidTestElt = infoTestElt.second;
        auto const& elem = _M_X1->mesh()->element( idTestElt,pidTestElt );

        auto const domains_eid_set = trialElementId( elem.id(), mpl::int_<nDimDiffBetweenTestTrial>() );
        //const uint16_type  n1_dof_on_element = element_dof1.size();
        const uint16_type  n1_dof_on_element = _M_X1->dof()->getIndicesSize(elem.id());

        for ( const size_type domain_eid : domains_eid_set )
        {
            // Get the global indices of the DOFs with support on this element
            _M_X2->dof()->getIndicesSetOnGlobalCluster( domain_eid, element_dof2 );

            // We can be more efficient if we sort the element DOFs
            // into increasing order
            //std::sort(element_dof1.begin(), element_dof1.end());
            std::sort( element_dof2.begin(), element_dof2.end() );

            const uint16_type  n2_dof_on_element = element_dof2.size();

            for ( size_type i=0; i<n1_dof_on_element; i++ )
            {
                // numLocal without ghosts ! very important for the graph with petsc
                const size_type il1 = _M_X1->dof()->localToGlobalId( elem.id(), i );// ig1 - _M_X1->dof()->firstDofGlobalCluster( theproc );
                const size_type ig1 = _M_X1->dof()->mapGlobalProcessToGlobalCluster()[il1];
                auto theproc = _M_X1->dof()->procOnGlobalCluster( ig1 );

                //const size_type ig1 = element_dof1[i];
                const int ndofpercomponent1 = n1_dof_on_element / _M_X1->dof()->nComponents;
                const int ncomp1 = i / ndofpercomponent1;
                const int ndofpercomponent2 = n2_dof_on_element / _M_X2->dof()->nComponents;

                {
                    graph_type::row_type& row = sparsity_graph->row( ig1 );
                    row.get<0>() = theproc ;
                    row.get<1>() = il1;

                    DVLOG(4) << "work with row " << ig1 << " local index " << il1 << " proc " << theproc << "\n";

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
                    if ( graph.test( Pattern::EXTENDED ) && hasNotFindRangeExtended )
                    {
                        for ( uint16_type ms=0; ms < elem.nNeighbors(); ms++ )
                        {
                            const auto * neighbor = boost::addressof( *_M_X1->mesh()->beginElementWithProcessId() /*elem*/ );
                            size_type neighbor_id = elem.neighbor( ms ).first;
                            size_type neighbor_process_id = elem.neighbor( ms ).second;

                            // warning ! the last condition is a temporary solution
                            if ( neighbor_id != invalid_size_type_value )
                            {
                                if ( neighbor_process_id != proc_id )
                                    CHECK( ( _M_X1->dof()->buildDofTableMPIExtended() &&
                                             _M_X2->dof()->buildDofTableMPIExtended() ) )
                                        << "Both spaces must have the extended dof table and none of them should be P0 Continuous to build the matrix stencil. Use block pattern construction instead!";

                                neighbor = boost::addressof( _M_X1->mesh()->element( neighbor_id,
                                                                                     neighbor_process_id ) );

                                if ( neighbor_id == neighbor->id()  )
                                {
                                    auto const domainsExtended_eid_set = trialElementId( neighbor_id/*elem.id()*/, mpl::int_<nDimDiffBetweenTestTrial>() );
                                    for ( const size_type neighborEltIdTrial : domainsExtended_eid_set )
                                    {

                                    neighbor_dof = _M_X2->dof()->getIndicesOnGlobalCluster( neighborEltIdTrial/*neighbor->id()*/ );

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
                                    }
                                } // neighbor_id
                            }
                        } // neighbor graph
                    } // if ( graph.test( Pattern::EXTENDED ) )
                } // only dof on proc
            } // dof loop
        } // trial id loop
        } // for ( size_type idTestElt : idsTestElt )
    } // element iterator loop
    } // rangeListTest

    if ( graph.test( Pattern::EXTENDED ) && !hasNotFindRangeExtended )
    {
        //auto rangeExtended = this->rangeExtendedIterator<0,0>( mpl::bool_<hasNotFindRangeExtended>() );
        auto rangeListExtended = this->rangeExtendedIterator<0,0>( mpl::bool_<hasNotFindRangeExtended>() );
        auto rangeExtended = rangeListExtended.front();
        auto iDimRangeExtended = rangeExtended.template get<0>();
        if ( iDimRangeExtended == ElementsType::MESH_ELEMENTS )
        {
            CHECK( false ) << "a range with MESH_ELEMENTS is not implemented";
        }
        else if ( iDimRangeExtended == ElementsType::MESH_FACES )
        {
            auto faceExtended_it = rangeExtended.template get<1>();
            auto faceExtended_en = rangeExtended.template get<2>();
            for ( ; faceExtended_it != faceExtended_en ;++faceExtended_it )
            {
                if ( !faceExtended_it->isConnectedTo0() || !faceExtended_it->isConnectedTo1() ) continue;

                if ( faceExtended_it->isInterProcessDomain() )
                    CHECK( ( _M_X1->dof()->buildDofTableMPIExtended() &&
                             _M_X2->dof()->buildDofTableMPIExtended() )  )
                        << "Both spaces must have the extended dof table and none of them should be P0 Continuous to build the matrix stencil. Use block pattern construction instead!";
#if 0
                    CHECK( _M_X1->dof()->buildDofTableMPIExtended() &&
                           _M_X2->dof()->buildDofTableMPIExtended() )
                        << "DofTableMPIExtended is not built!";
#endif

                auto const& elt0 = faceExtended_it->element0();
                auto const& elt1 = faceExtended_it->element1();

                const uint16_type  n1_dof_on_element = _M_X1->dof()->getIndicesSize();
                const uint16_type  n2_dof_on_element = _M_X2->dof()->getIndicesSize();

                neighbor_dof = _M_X2->dof()->getIndicesOnGlobalCluster( elt1.id() );

                for ( size_type i=0; i<n1_dof_on_element; i++ )
                {
                    const size_type ig1 = _M_X1->dof()->mapGlobalProcessToGlobalCluster()[_M_X1->dof()->localToGlobalId( elt0.id(), i )];
                    auto theproc = _M_X1->dof()->procOnGlobalCluster( ig1 );
                    const size_type il1 = _M_X1->dof()->localToGlobalId( elt0.id(), i );

                    const int ndofpercomponent1 = n1_dof_on_element / _M_X1->dof()->nComponents;
                    const int ncomp1 = i / ndofpercomponent1;
                    const int ndofpercomponent2 = n2_dof_on_element / _M_X2->dof()->nComponents;

                    graph_type::row_type& row = sparsity_graph->row( ig1 );
                    row.get<0>() = theproc ;
                    row.get<1>() = il1;

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
                }

                neighbor_dof = _M_X2->dof()->getIndicesOnGlobalCluster( elt0.id() );
                for ( size_type i=0; i<n1_dof_on_element; i++ )
                {
                    const size_type ig1 = _M_X1->dof()->mapGlobalProcessToGlobalCluster()[_M_X1->dof()->localToGlobalId( elt1.id(), i )];
                    auto theproc = _M_X1->dof()->procOnGlobalCluster( ig1 );
                    const size_type il1 = _M_X1->dof()->localToGlobalId( elt1.id(), i );

                    const int ndofpercomponent1 = n1_dof_on_element / _M_X1->dof()->nComponents;
                    const int ncomp1 = i / ndofpercomponent1;
                    const int ndofpercomponent2 = n2_dof_on_element / _M_X2->dof()->nComponents;

                    graph_type::row_type& row = sparsity_graph->row( ig1 );
                    row.get<0>() = theproc ;
                    row.get<1>() = il1;

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

                }

            } // for ( ; faceExtended_it != faceExtended_en ;++faceExtended_it )
        } // if ( iDimRangeExtended == ElementsType::MESH_FACES )
    } // if ( graph.test( Pattern::EXTENDED ) && !hasNotFindRangeExtended )

    DVLOG(2) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    return sparsity_graph;
}
#endif
namespace detail
{
template<typename EltType>
struct gmcDefStencil
{
    typedef typename EltType::gm_type::precompute_type pc_type;
    typedef typename EltType::gm_type::precompute_ptrtype pc_ptrtype;
    typedef typename EltType::gm_type::template Context<vm::POINT, EltType> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
};

template<typename FaceType>
struct gmcDefFaceStencil
{
    typedef typename boost::remove_reference<FaceType>::type FaceType1;
    typedef typename boost::remove_const<FaceType1>::type FaceType2;

    typedef typename FaceType2::super2::template Element<FaceType2>::type EltType;
    typedef typename gmcDefStencil<EltType>::pc_type pc_type;
    typedef typename gmcDefStencil<EltType>::pc_ptrtype pc_ptrtype;
    typedef typename gmcDefStencil<EltType>::gmc_type gmc_type;
    typedef typename gmcDefStencil<EltType>::gmc_ptrtype gmc_ptrtype;
};

template<typename ImType,typename EltType>
typename gmcDefStencil<EltType>::gmc_ptrtype
gmcStencil( mpl::size_t<MESH_ELEMENTS> /**/, EltType const& elem )
{
    typedef typename gmcDefStencil<EltType>::pc_type pc_type;
    typedef typename gmcDefStencil<EltType>::pc_ptrtype pc_ptrtype;
    typedef typename gmcDefStencil<EltType>::gmc_type gmc_type;
    typedef typename gmcDefStencil<EltType>::gmc_ptrtype gmc_ptrtype;

    ImType theim;
    pc_ptrtype geopc( new pc_type( elem.gm(), theim.points() ) );
    gmc_ptrtype gmc( new gmc_type( elem.gm(), elem, geopc ) );
    return gmc;
}

template<typename ImType,typename FaceType>
typename gmcDefFaceStencil<FaceType>::gmc_ptrtype
gmcStencil( mpl::size_t<MESH_FACES> /**/, FaceType const& theface )
{
    typedef typename gmcDefFaceStencil<FaceType>::pc_type pc_type;
    typedef typename gmcDefFaceStencil<FaceType>::pc_ptrtype pc_ptrtype;
    typedef typename gmcDefFaceStencil<FaceType>::gmc_type gmc_type;
    typedef typename gmcDefFaceStencil<FaceType>::gmc_ptrtype gmc_ptrtype;

    typedef typename QuadMapped<ImType>::permutation_type permutation_type;
    typedef typename QuadMapped<ImType>::permutation_points_type permutation_points_type;

    ImType theim;
    QuadMapped<ImType> qm;
    permutation_points_type ppts( qm(theim) );

    std::vector<std::map<permutation_type, pc_ptrtype> > __geopc( theim.nFaces() );
    for ( uint16_type __f = 0; __f < theim.nFaces(); ++__f )
    {
        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            __geopc[__f][__p] = pc_ptrtype(  new pc_type( theface.element( 0 ).gm(), ppts[__f].find( __p )->second ) );
        }
    }

    uint16_type __face_id_in_elt_0 = theface.pos_first();
    gmc_ptrtype gmc( new gmc_type( theface.element( 0 ).gm(),
                                   theface.element( 0 ),
                                   __geopc,
                                   __face_id_in_elt_0 ) );
    return gmc;
}
template<typename EltType>
void
gmcUpdateStencil( mpl::size_t<MESH_ELEMENTS> /**/, EltType const& elem, typename gmcDefStencil<EltType>::gmc_ptrtype &gmc )
{
    gmc->update( elem );
}
template<typename FaceType>
void
gmcUpdateStencil( mpl::size_t<MESH_FACES> /**/, FaceType const& theface, typename gmcDefFaceStencil<FaceType>::gmc_ptrtype &gmc )
{
    gmc->update( theface.element( 0 ), theface.pos_first() );
}


template<typename EltType>
std::set<size_type>
idEltStencil( mpl::size_t<MESH_ELEMENTS> /**/, EltType const& elem )
{
    std::set<size_type> res;
    res.insert( elem.id() );
    return res;
}
template<typename FaceType>
std::set<size_type>
idEltStencil( mpl::size_t<MESH_FACES> /**/, FaceType const& theface )
{
    std::set<size_type> res;
    if ( theface.isConnectedTo0() )
        res.insert( theface.element( 0 ).id() );
    if ( theface.isConnectedTo1() )
        res.insert( theface.element( 1 ).id() );
    return res;
}

} // namespace detail


template<typename X1,  typename X2,typename RangeItTestType, typename RangeExtendedItType, typename QuadSetType>
typename Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::graph_ptrtype
Stencil<X1,X2,RangeItTestType,RangeExtendedItType,QuadSetType>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true> )
{
    //std::cout << "\n start graphInterp "<< std::endl;
    typedef mpl::int_< boost::remove_reference<typename fusion::result_of::at_c<QuadSetType,0>::type>::type::order > order_1d_type;
    typedef mpl::int_< boost::remove_reference<typename fusion::result_of::at_c<QuadSetType,1>::type>::type::order > order_2d_type;
    typedef mpl::int_< boost::remove_reference<typename fusion::result_of::at_c<QuadSetType,2>::type>::type::order > order_3d_type;

    typedef typename test_space_type::mesh_type test_mesh_type;
    typedef typename trial_space_type::mesh_type trial_mesh_type;
    typedef typename mpl::if_< mpl::equal_to<mpl::int_<test_mesh_type::nDim>,mpl::int_<1> >,
            order_1d_type,
            typename mpl::if_<mpl::equal_to<mpl::int_<test_mesh_type::nDim>,mpl::int_<2> >,
            order_2d_type,
            order_3d_type >::type>::type order_used_type;
    typedef typename mpl::if_<mpl::bool_<test_mesh_type::element_type::is_simplex>,
            mpl::identity<typename _Q<order_used_type::value>::template applyIMGeneral<test_mesh_type::element_type::nDim, typename test_mesh_type::value_type, Simplex>::type >,
                mpl::identity<typename _Q<order_used_type::value>::template applyIMGeneral<test_mesh_type::element_type::nDim, typename test_mesh_type::value_type, Hypercube>::type >
    >::type::type theim_type;

    typedef typename test_mesh_type::Localization::matrix_node_type matrix_node_type;

    //-----------------------------------------------------------------------//

    const size_type proc_id           = _M_X1->mesh()->comm().rank();
    const size_type n1_dof_on_proc    = _M_X1->nLocalDof();
    //const size_type n2_dof_on_proc    = _M_X2->nLocalDof();
    const size_type first1_dof_on_proc = _M_X1->dof()->firstDof( proc_id );
    const size_type last1_dof_on_proc = _M_X1->dof()->lastDof( proc_id );
    const size_type first2_dof_on_proc = _M_X2->dof()->firstDof( proc_id );
    const size_type last2_dof_on_proc = _M_X2->dof()->lastDof( proc_id );

    graph_ptrtype sparsity_graph( new graph_type( _M_X1->dof(),_M_X2->dof() ) );

    //-----------------------------------------------------------------------//
    // init localisation tools
    auto locToolForXh2 = _M_X2->mesh()->tool_localization();
    locToolForXh2->updateForUse();
    bool doExtrapolationAtStartXh2 = locToolForXh2->doExtrapolation();
    if (doExtrapolationAtStartXh2) locToolForXh2->setExtrapolation( false );
    const auto nbNearNeighborAtStartTrial = locToolForXh2->kdtree()->nPtMaxNearNeighbor();
    bool notUseOptLocTrial = trial_mesh_type::nDim!=trial_mesh_type::nRealDim;
    if (notUseOptLocTrial) locToolForXh2->kdtree()->nbNearNeighbor(trial_mesh_type::element_type::numPoints);

    //locTool->kdtree()->nbNearNeighbor(_M_X2->mesh()->numElements() );
    auto locToolForXh1 = _M_X1->mesh()->tool_localization();
    locToolForXh1->updateForUse();
    bool doExtrapolationAtStartXh1 = locToolForXh1->doExtrapolation();
    if (doExtrapolationAtStartXh1) locToolForXh1->setExtrapolation( false );


    size_type IdEltInXh2 = invalid_size_type_value;
    //node_type trialNodeRef,testNodeRef;

#if FEELPP_EXPORT_GRAPH
    std::map<size_type,std::list<size_type> > mapBetweenMeshes;
#endif

    std::vector<size_type> element_dof1_range, element_dof1, element_dof2;
    std::set<size_type> neighLocalizedInXh1;

    std::set<size_type > listTup;

    //-----------------------------------------------------------------------//

    static const bool hasNotFindRangeStandard = rangeiteratorType<0,0>::hasnotfindrange_type::value;
    auto rangeListTest = this->rangeiterator<0,0>( mpl::bool_<hasNotFindRangeStandard>() );
    for ( auto const& rangeTest : rangeListTest )
    {
    auto iDimRange = rangeTest.template get<0>();
    auto elem_it = rangeTest.template get<1>();
    auto elem_en = rangeTest.template get<2>();

    if ( elem_it==elem_en ) return sparsity_graph;

    CHECK( elem_it->mesh()->isSameMesh( _M_X1->mesh() ) ) << "case not take into account : mesh range is not the same that test mesh";

    matrix_node_type ptsReal( elem_it->vertices().size1(), 1 );
    auto gmc = Feel::detail::gmcStencil<theim_type>( iDimRange,*elem_it );

    if ( _M_X1->nDof()>1 )
    {
        for ( ; elem_it != elem_en; ++elem_it )
        {
            auto const& elem = *elem_it;
            // Get the global indices of the DOFs with support on this element
            element_dof1_range = _M_X1->dof()->getIndices( elem.id(), iDimRange );

            const std::set<size_type> idsElt = Feel::detail::idEltStencil( iDimRange,elem );
            if ( idsElt.empty() ) continue;
            element_dof1 = _M_X1->dof()->getIndices( *(idsElt.begin()) );
            const uint16_type n1_dof_on_element_range = element_dof1_range.size();
            const uint16_type n1_dof_on_element = element_dof1.size();

            std::vector<boost::tuple<bool,size_type> > hasFinds( n1_dof_on_element_range,boost::make_tuple( false,invalid_size_type_value ) );

            for ( size_type i=0; i<n1_dof_on_element_range; i++ )
            {
                const size_type ig1 = element_dof1_range[i];
                auto const ptRealDof = boost::get<0>( _M_X1->dof()->dofPoint( ig1 ) );

                ublas::column(ptsReal,0 ) = ptRealDof;
                if (notUseOptLocTrial) IdEltInXh2=invalid_size_type_value;
                auto resLocalisationInXh2 = locToolForXh2->run_analysis(ptsReal,IdEltInXh2,elem_it->vertices(),mpl::int_<0>());
                IdEltInXh2 = resLocalisationInXh2.template get<1>();
                bool hasFind = resLocalisationInXh2.template get<0>()[0];

                listTup.clear();

                if ( hasFind )
                {
                    listTup.insert( IdEltInXh2 );
                    hasFinds[i] = boost::make_tuple( true,IdEltInXh2 );
                    // maybe is on boundary->more elts
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
            } // for (size_type i=0; i<n1_dof_on_element_range; i++)

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
                Feel::detail::gmcUpdateStencil( iDimRange,elem,gmc );

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
            auto const& elem = *elem_it;

            // Get the global indices of the DOFs with support on this element
            const std::set<size_type> idsElt = Feel::detail::idEltStencil( iDimRange,elem );
            if ( idsElt.empty() ) continue;
            element_dof1 = _M_X1->dof()->getIndices( *(idsElt.begin()) );

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
                    auto const& geoelt2 = _M_X2->mesh()->element( idElt2 );// warning miss processId
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
                                        auto const& geoeltNEW = _M_X2->mesh()->element( *it_neigh );// warning miss processId
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
                Feel::detail::gmcUpdateStencil(iDimRange,elem,gmc);

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

    } // rangeListTest

    if (doExtrapolationAtStartXh1) locToolForXh1->setExtrapolation( true );
    if (doExtrapolationAtStartXh2) locToolForXh2->setExtrapolation( true );
    if (notUseOptLocTrial) locToolForXh2->kdtree()->nbNearNeighbor(nbNearNeighborAtStartTrial);

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
