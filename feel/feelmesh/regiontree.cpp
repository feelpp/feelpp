/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-01

  Copyright (C) 2005,2006 EPFL

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
   \file RegionTree.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-02-01
 */
#include <boost/timer.hpp>

#include <feel/feelmesh/regiontree.hpp>

namespace Feel
{
/**
 *
 */
typedef node<double>::type node_type;
/// \cond detail
struct RegionTree::element_base
{
    enum { RECTS_PER_LEAF=8 };
    bool isleaf() const
    {
        return isleaf_;
    }
    element_base( bool leaf,
                  const node_type& rmin_,
                  const node_type& rmax_ )
        :
        isleaf_( leaf ),
        rmin( rmin_ ),
        rmax( rmax_ )
    {}

    bool isleaf_;
    node_type rmin, rmax;
};

struct tree_node
        :
    public RegionTree::element_base
{
    tree_node( const node_type& bmin,
               const node_type& bmax,
               RegionTree::element_base *left_,
               RegionTree::element_base *right_ )
        :
        RegionTree::element_base( false, bmin, bmax ),
        left( left_ ),
        right( right_ )
    {}

    element_base *left;
    element_base *right;
};

struct leaf
        :
    public RegionTree::element_base
{
    leaf( const node_type& bmin,
          const node_type& bmax,
          RegionTree::pbox_container_type& lst_ )
        :
        RegionTree::element_base( true, bmin, bmax )
    {
        lst.swap( lst_ );
    }

    RegionTree::pbox_container_type lst;
};

/* enlarge box to hold [a..b] */
FEELPP_NO_EXPORT
void
updateBox( node_type& bmin, node_type& bmax,
           const node_type& a, const node_type& b )
{
    node_type::iterator itmin = bmin.begin();
    node_type::iterator itmax = bmax.begin();

    for ( size_type i=0; i < a.size(); ++i )
    {
        bmin[i] = std::min( bmin[i], a[i] );
        bmax[i] = std::max( bmax[i], b[i] );
    }
}

FEELPP_NO_EXPORT
bool
r1_ge_r2( const node_type& min1, const node_type& max1,
          const node_type& min2, const node_type& max2 )
{
    for ( size_type i=0; i < min1.size(); ++i )
    {
        if ( !( min1[i] <= min2[i] && max1[i] >= max2[i] ) )
            return false;
    }

    return true;
}

FEELPP_NO_EXPORT
bool
r1_inter_r2( const node_type& min1, const node_type& max1,
             const node_type& min2, const node_type& max2 )
{
    for ( size_type i=0; i < min1.size(); ++i )
    {
        if ( max1[i] < min2[i] || min1[i] > max2[i] )
            return false;
    }

    return true;
}

/* some predicates for searches */
struct FEELPP_NO_EXPORT intersection_p
{
    intersection_p( const node_type& min_, const node_type& max_ )
        :
        min( min_ ),
        max( max_ )
    {}
    bool operator()( const node_type& min2, const node_type& max2 )
    {
        return r1_inter_r2( min,max,min2,max2 );
    }

    const node_type min,max;

};

/* match boxes containing [min..max] */
struct FEELPP_NO_EXPORT contains_p
{
    contains_p( const node_type& min_, const node_type& max_ )
        :
        min( min_ ),
        max( max_ )
    {}

    bool operator()( const node_type& min2, const node_type& max2 )
    {
        return r1_ge_r2( min2,max2,min,max );
    }

    const node_type min,max;
};

/* match boxes contained in [min..max] */
struct FEELPP_NO_EXPORT contained_p
{
    contained_p( const node_type& min_, const node_type& max_ )
        :
        min( min_ ),
        max( max_ )
    {}
    bool operator()( const node_type& min2, const node_type& max2 )
    {
        FEELPP_ASSERT( min.size() == min2.size() &&
                       min.size() == max2.size() &&
                       max.size() == min2.size() &&
                       max.size() == max2.size() )

        ( min.size() )( max.size() )( min2.size() )( max2.size() ).error( "invalid box size" );

        return r1_ge_r2( min,max,min2,max2 );
    }

    const node_type min,max;
};

struct FEELPP_NO_EXPORT has_point_p
{
    has_point_p( const node_type& P_ ) : P( P_ ) {}

    bool operator()( const node_type& min2, const node_type& max2 )
    {
        FEELPP_ASSERT( P.size() == min2.size() &&
                       P.size() == max2.size() )
        ( P.size() )( min2.size() )( max2.size() ).error( "invalid point size" );

        for ( size_type i=0; i < P.size(); ++i )
        {
            if ( P[i] < min2[i] || P[i] > max2[i] )
                return false;
        }

        return true;
    }

    const node_type P;
};


template <typename Predicate>
FEELPP_NO_EXPORT void findMatchingBoxes( RegionTree::element_base *n,
        RegionTree::pbox_set_type& boxlst,
        Predicate p )
{
    DVLOG(2) << "find_matching_boxes_: "
                  << n->rmin << ".."
                  << n->rmax << "\n";

    if ( n->isleaf() )
    {
        DVLOG(2) << "findMatchingBoxes in leaf\n";

        const leaf *rl = static_cast<leaf*>( n );

        for ( RegionTree::pbox_container_type::const_iterator it = rl->lst.begin();
                it != rl->lst.end();
                ++it )
        {
            DVLOG(2) << "  ->match(" << ( *it )->id << "="
                          << ( *it )->min << "," << ( *it )->max << " -> "
                          << p( ( *it )->min, ( *it )->max ) << "\n";

            if ( p( ( *it )->min, ( *it )->max ) )
            {
                boxlst.insert( *it );
            }
        }
    }

    else
    {
        DVLOG(2) << "findMatchingBoxes in branch\n";

        const tree_node *rn = static_cast<tree_node*>( n );

        if ( p( rn->left->rmin,rn->left->rmax ) )
        {
            findMatchingBoxes( rn->left, boxlst, p );
        }

        if ( p( rn->right->rmin,rn->right->rmax ) )
        {
            findMatchingBoxes( rn->right, boxlst, p );
        }
    }
}

/// \endcond detail
void
RegionTree::addBox( node_type min, node_type max, size_type id )
{
    box_index_type bi;
    bi.min = min;
    bi.max = max;
    bi.id = ( id + 1 ) ? id : M_boxes.size();
    M_boxes.push_back( bi );
}
void
RegionTree::findIntersectingBoxes( const node_type& bmin,
                                   const node_type& bmax,
                                   pbox_set_type& boxes )
{
    boxes.clear();

    if ( !M_root )
        build();

    findMatchingBoxes( M_root, boxes, intersection_p( bmin,bmax ) );
    DVLOG(2)<< "findIntersectingBoxes : found " << boxes.size() << " matches\n";
}
void
RegionTree::findIntersectingBoxes( const node_type& bmin,
                                   const node_type& bmax,
                                   std::vector<size_type>& idvec )
{
    pbox_set_type bs;
    findIntersectingBoxes( bmin, bmax, bs );
    toIdVector( bs, idvec );
}
void
RegionTree::findContainingBoxes( const node_type& bmin,
                                 const node_type& bmax,
                                 pbox_set_type& boxes )
{
    boxes.clear();

    if ( !M_root )
        build();

    findMatchingBoxes( M_root, boxes, contains_p( bmin,bmax ) );
    DVLOG(2)<< "findContainingBoxes : found " << boxes.size() << " matches\n";
}
void
RegionTree::findContainingBoxes( const node_type& bmin,
                                 const node_type& bmax,
                                 std::vector<size_type>& idvec )
{
    pbox_set_type bs;
    findContainingBoxes( bmin, bmax, bs );
    toIdVector( bs, idvec );
}
void
RegionTree::findContainedBoxes( const node_type& bmin,
                                const node_type& bmax,
                                pbox_set_type& boxes )
{
    boxes.clear();

    if ( !M_root )
        build();

    findMatchingBoxes( M_root, boxes, contained_p( bmin,bmax ) );
    DVLOG(2)<< "findContainedBoxes : found " << boxes.size() << " matches\n";
}
void
RegionTree::findContainedBoxes( const node_type& bmin,
                                const node_type& bmax,
                                std::vector<size_type>& idvec )
{
    pbox_set_type bs;
    findContainedBoxes( bmin, bmax, bs );
    toIdVector( bs, idvec );
}
void
RegionTree::findBoxesAtPoint( const node_type& P,
                              pbox_set_type& boxes )
{
    boxes.clear();

    if ( !M_root )
        build();

    findMatchingBoxes( M_root, boxes, has_point_p( P ) );
    DVLOG(2)<< "findBoxesAtPointb : found " << boxes.size() << " matches\n";
}

void
RegionTree::findBoxesAtPoint( const node_type& P,
                              std::vector<size_type>& idvec )
{
    pbox_set_type bs;
    findBoxesAtPoint( P, bs );
    toIdVector( bs, idvec );
}

void
RegionTree::toIdVector( pbox_set_type const& bs, std::vector<size_type>& idvec )
{
    idvec.reserve( bs.size() );
    idvec.resize( 0 );

    for ( pbox_set_type::const_iterator it=bs.begin(); it != bs.end(); ++it )
    {
        idvec.push_back( ( *it )->id );
    }
}


/*
  try to split at the approximate center of the box. Could be much more
  sophisticated
*/
FEELPP_NO_EXPORT
bool
splitTest( const RegionTree::pbox_container_type& b,
           const node_type& bmin,
           const node_type& bmax,
           size_type dir,
           scalar_type& split_v )
{
    scalar_type v = bmin[dir] + ( bmax[dir] - bmin[dir] )/2;
    split_v = v;
    size_type cnt = 0;

    DVLOG(2) << "[enter]Split_test: dir=" << dir
                  << ", split_v=" << v
                  << ", bmin=" << bmin
                  << ", bmax=" << bmax << "\n";

    RegionTree::pbox_container_type::const_iterator it = b.begin();
    RegionTree::pbox_container_type::const_iterator en = b.end();

    while ( it != en )
    {
        if ( ( *it )->max[dir] < v )
        {
            if ( cnt == 0 )
            {
                split_v = ( *it )->max[dir];
            }

            else
            {
                split_v = std::max( ( *it )->max[dir],split_v );
            }

            cnt++;
        }

        ++it;
    }

    DVLOG(2) << "[exit] Split_test cnt = " << cnt
                  << ", b.size()=" << b.size()
                  << ", split_v=" << split_v << "\n";

    return ( cnt > 0 && cnt < b.size() );
}
/*
  there are many flavors of rtree ... this one is more or less a quadtree
  where splitting does not occurs at predefined positions (hence the split_test
  function above). Regions of the tree do not overlap (box are splitted).
*/
FEELPP_NO_EXPORT
RegionTree::element_base*
build( RegionTree::pbox_container_type& b,
       const node_type& bmin,
       const node_type& bmax,
       size_type last_dir )
{
    size_type N=bmin.size();
    scalar_type split_v;
    size_type split_dir = ( last_dir+1 )%N;

    DVLOG(2) << " build_tree_ [b.size=" << b.size() << "],"
                  << "bmin=" << bmin << ", "
                  << "bmax=" << bmax << "\n";

    bool split_ok = false;

    if ( b.size() > RegionTree::element_base::RECTS_PER_LEAF )
    {
        for ( size_type itry=0; itry < N; ++itry )
        {
            DVLOG(2) << "split_test: dir=" << split_dir << "\n";

            if ( splitTest( b, bmin, bmax, split_dir, split_v ) )
            {
                split_ok = true;
                break;
            }

            split_dir = ( split_dir+1 )%N;
        }
    }

    if ( split_ok )
    {
        size_type cnt1=0,cnt2=0;

        DVLOG(2) << "splitting with v=" << split_v << "\n";

        RegionTree::pbox_container_type::const_iterator it = b.begin();
        RegionTree::pbox_container_type::const_iterator en = b.end();

        while ( it != en )
        {
            DVLOG(2) << " . test box" << ( *it )->min[split_dir] << ".." << ( *it )->max[split_dir] << "\n";

            if ( ( *it )->min[split_dir] < split_v )
            {
                cnt1++;
            }

            if ( ( *it )->max[split_dir] > split_v )
            {
                cnt2++;
            }

            ++it;
        }

        DVLOG(2) << "  -> left : " << cnt1 << " boxes, right : " << cnt2 << " boxes\n";

        FEELPP_ASSERT( cnt1 )( cnt1 ).error( "counter 1 is 0" );
        FEELPP_ASSERT( cnt2 )( cnt2 ).error( "counter 2 is 0" );
        FEELPP_ASSERT( cnt1+cnt2 >= b.size() )( cnt1 )( cnt2 )( cnt1+cnt2 )( b.size() ).error( "counter sum should be greater or equal to the number of boxes" );

        RegionTree::pbox_container_type v1( cnt1 );
        RegionTree::pbox_container_type v2( cnt2 );
        node_type bmin1( bmax ), bmax1( bmin );
        node_type bmin2( bmax ), bmax2( bmin );
        cnt1 = cnt2 = 0;

        it = b.begin();
        en = b.end();

        while ( it != en )
        {
            if ( ( *it )->min[split_dir] < split_v )
            {
                v1[cnt1++] = *it;
                updateBox( bmin1, bmax1,( *it )->min,( *it )->max );
                DVLOG(2) << "update_box bmin1=" << bmin1 << ", bmax1=" << bmax1 << "\n";
            }

            if ( ( *it )->max[split_dir] > split_v )
            {
                v2[cnt2++] = *it;
                updateBox( bmin2,bmax2,( *it )->min,( *it )->max );
            }

            ++it;
        }

        for ( size_type k=0; k < N; ++k )
        {
            bmin1[k] = std::max( bmin1[k],bmin[k] );
            bmax1[k] = std::min( bmax1[k],bmax[k] );
            bmin2[k] = std::max( bmin2[k],bmin[k] );
            bmax2[k] = std::min( bmax2[k],bmax[k] );
        }

        bmax1[split_dir] = std::min( bmax1[split_dir], split_v );
        bmin2[split_dir] = std::max( bmin2[split_dir], split_v );
        FEELPP_ASSERT( cnt1 == v1.size() )( cnt1 )( v1.size() ).error( "sizes should be equal" );;
        FEELPP_ASSERT( cnt2 == v2.size() )( cnt2 )( v2.size() ).error( "sizes should be equal" );;
        return new tree_node( bmin,bmax,
                              build( v1, bmin1, bmax1, split_dir ),
                              build( v2, bmin2, bmax2, split_dir ) );
    }

    else
    {
        return new leaf( bmin,bmax,b );
    }
}

void
RegionTree::build()
{
    DVLOG(2) << "build tree\n";
    boost::timer __timer;

    if ( M_boxes.size() == 0 )
        return;

    FEELPP_ASSERT( M_root == 0 ).error( "the tree has already been built" );

    pbox_container_type b( M_boxes.size() );
    pbox_container_type::iterator b_it = b.begin();

    node_type bmin( M_boxes.front().min );
    node_type bmax( M_boxes.front().max );

    // build the bounding box of all the boxes stored
    // and store the boxes address in the pbox_container_type
    box_container_type::const_iterator it=M_boxes.begin();
    box_container_type::const_iterator en=M_boxes.end();

    while ( it != en )
    {
        updateBox( bmin, bmax,
                   ( *it ).min,( *it ).max );
        *b_it++ = &( *it );

        ++it;
    }

    // build a tree view of the boxes
    M_root = Feel::build( b, bmin, bmax, 0 );

    DVLOG(2) << "build tree done in " << __timer.elapsed() << "s\n";
}


FEELPP_NO_EXPORT
void
dump( RegionTree::element_base *p, int level, size_type& count )
{
    if ( !p )
        return;

    std::cout << level << "|";

    for ( int i=0; i < level; ++i )
        std::cout << "  ";

    std::cout << "span=" << p->rmin << ".." << p->rmax << " ";

    if ( p->isleaf() )
    {
        leaf *rl = static_cast<leaf*>( p );

        std::cout << "Leaf [" << rl->lst.size() << " elts] = ";

        for ( size_type i=0; i < rl->lst.size(); ++i )
            std::cout << " " << rl->lst[i]->id;

        std::cout << "\n";

        count += rl->lst.size();
    }

    else
    {
        std::cout << "Node\n";
        const tree_node *rn = static_cast<tree_node*>( p );

        if ( rn->left )
        {
            dump( rn->left, level+1, count );
        }

        if ( rn->right )
        {
            dump( rn->right, level+1, count );
        }
    }
}

void
RegionTree::dump()
{
    std::cout << "tree dump follows\n";

    if ( !M_root )
        build();

    size_type count = 0;
    Feel::dump( M_root, 0, count );

    std::cout << " --- end of tree dump, nb of boxes: " << M_boxes.size()
              << ", rectangle ref in tree: " << count << "\n";
}


FEELPP_NO_EXPORT
void
destroy( RegionTree::element_base *n )
{
    if ( n->isleaf() )
    {
        delete static_cast<leaf*>( n );
    }

    else
    {
        const tree_node *rn = static_cast<tree_node*>( n );

        if ( rn->left )
        {
            destroy( rn->left );
        }

        if ( rn->right )
        {
            destroy( rn->right );
        }

        delete rn;
    }
}

void
RegionTree::destroy()
{
    if (  M_root )
        Feel::destroy( M_root );
}

}
