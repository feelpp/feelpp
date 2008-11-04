/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-07

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file kdtree.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-06-07
 */
#include <life/lifemesh/kdtree.hpp>

namespace Life
{

struct KDTree::Element
{
    enum { POINTS_PER_LEAF = 8 };

    /* 0 => is a tree node, != 0 => tree leaf storing n points */
    size_type n;

    virtual ~Element() {}
    bool isleaf() const { return (n != 0); }
    Element( size_type n_ ) : n(n_) {}
};

struct KDTree::Leaf : public KDTree::Element
{
    KDTree::points_const_iterator it;
    Leaf( KDTree::points_const_iterator begin,
          KDTree::points_const_iterator end)
        :
        Element( std::distance(begin,end) )
    {
        LIFE_ASSERT( begin != end )( std::distance(begin,end) ).error( "invalid iterator" );
        it = begin;
    }
};

struct KDTree::Node : public KDTree::Element
{
    scalar_type split_v;

    /* left: <=v, right: >v */
    Element *left, *right;

    Node(scalar_type v, Element *left_, Element *right_)
        :
        Element(0),
        split_v(v),
        left(left_),
        right(right_)
    {}
};

namespace detail
{
/**
 *sorting of points_type with respect to a component of the points
 */
struct component_sort
{
    component_sort(unsigned d)
        :
        dir(d)
    {}
    bool operator()( KDTree::index_node_type const& a,
                     KDTree::index_node_type const& b )
    {
        return boost::get<0>(a)(dir) < boost::get<0>(b)(dir);
    }
    unsigned dir;
};

static void
swap( KDTree::index_node_type& a,
      KDTree::index_node_type& b )
{
    KDTree::index_node_type tmp( a );
    a = b;
    b = tmp;
}
static
KDTree::points_iterator
partition( KDTree::points_iterator begin,
           KDTree::points_iterator end,
           unsigned dir,
           scalar_type median )
{
    --end;
    do
        {
            while ( (begin <= end) &&
                    (boost::get<0>(*begin)(dir) <= median) )
                ++begin;
            while ( (begin <= end) &&
                    ( boost::get<0>(*end)(dir) > median) )
                --end;
            if (begin < end)
                {
                    //begin->swap(*end);
                    swap( *begin, *end );
                    ++begin;
                    --end;
                }
            else break;
        }
    while (1);
    return begin;
}
/*
  build (recursively) a complete tree for points in the interval [begin, end[
  dir is the splitting direction
*/
static
KDTree::Element*
build_tree( KDTree::points_iterator begin,
            KDTree::points_iterator end,
            unsigned dir )
{
    if (begin == end)
        return 0;

    size_type npts = std::distance(begin,end);
    Debug( 4011 ) << "[KDTree::build_tree] nbpts = " << npts << "\n";

    if (npts > KDTree::Element::POINTS_PER_LEAF)
        {
            Debug( 4011 ) << "[KDTree::build_tree] split\n";
            KDTree::points_iterator itmedian;
            scalar_type median;
            size_type N = boost::get<0>(*begin).size();
            unsigned ndir_tests = dir/N;
            dir %= N;
            if (npts > 50)
                {
                    /* too much points for an exact median: estimation of the median .. */
                    std::vector<KDTree::index_node_type> v(30);
                    //random_sample(begin,end,v.begin(),v.end());
                    for (size_type i=0; i < v.size(); ++i)
                        v[i] = begin[rand() % npts];
                    std::sort(v.begin(), v.end(), component_sort(dir));
                    median = ( boost::get<0>( v[v.size()/2-1] )(dir)+
                               boost::get<0>( v[v.size()/2] )(dir))/2;
                    itmedian = partition(begin,end,dir,median);
                }
            else
                {
                    /* exact median computation */
                    std::sort(begin, end, component_sort(dir));
                    itmedian = begin + npts/2 - 1;
                    median = boost::get<0>(*itmedian)[dir];
                    while ( itmedian < end &&
                            boost::get<0>(*itmedian)[dir] == median )
                        itmedian++;
                }
            /* could not split the set (all points have same value for component 'dir' !) */
            if (itmedian == end)
                {
                    /* tested all N direction ? so all points are strictly the same */
                    if (ndir_tests == N-1)
                        return new KDTree::Leaf(begin,end);
                    else
                        return new KDTree::Node(median,
                                                build_tree(begin,
                                                           itmedian,
                                                           (dir+1)%N + (ndir_tests+1)*N),
                                                0);
                }
            else
                {
                    /* the general case */
                    LIFE_ASSERT( (boost::get<0>(*itmedian)[dir] > median &&
                                 boost::get<0>(*(itmedian-1))[dir] <= median ) )
                        ( dir )
                        ( boost::get<0>(*itmedian)[dir])
                        (boost::get<0>(*(itmedian-1))[dir])
                        (median).error( "invalid median iterator" );

                    return new KDTree::Node(median,
                                            build_tree(begin, itmedian, (dir+1)%N),
                                            build_tree(itmedian,end, (dir+1)%N));
                }
        }
    else
        {
            Debug( 4011 ) << "[KDTree::build_tree] return leaf\n";
            return new KDTree::Leaf(begin,end);
        }
}

static void destroy_tree(KDTree::Element *t)
{
    if (t == 0)
        return;
    if (!t->isleaf())
        {
            KDTree::Node *tn = static_cast<KDTree::Node*>(t);
            destroy_tree(tn->right);
            destroy_tree(tn->left);
            delete tn;
        }
    else
        {
            KDTree::Leaf *tl = static_cast<KDTree::Leaf*>(t);
            delete tl;
        }
}

/* avoid pushing too much arguments on the stack for points_in_box_ */
struct points_in_box_data
{
    KDTree::node_type::const_iterator bmin;
    KDTree::node_type::const_iterator bmax;
    KDTree::points_type *ipts;
    size_type N;
};

/* recursive lookup for points inside a given box */
static
void
points_in_box(const points_in_box_data& p,
              const KDTree::Element *t,
              unsigned dir)
{
    if (!t->isleaf())
        {
            const KDTree::Node *tn = static_cast<const KDTree::Node*>(t);
            if (p.bmin[dir] <= tn->split_v && tn->left)
                {
                    points_in_box( p,
                                   tn->left,
                                   (dir+1)%p.N );
                }
            if (p.bmax[dir] > tn->split_v && tn->right)
                {
                    points_in_box( p,
                                   tn->right,
                                   (dir+1)%p.N );
                }
        }
    else
        {
            const KDTree::Leaf *tl = static_cast<const KDTree::Leaf*>(t);
            KDTree::points_const_iterator itpt = tl->it;
            for (size_type i=tl->n;i; --i, ++itpt)
                {
                    Debug( 4011 ) <<  "i = " << i << " pt = " << boost::get<0>(*itpt) << " ptindex = " << boost::get<1>(*itpt)<< "\n";
                    bool is_in = true;
                    KDTree::node_type::const_iterator it = boost::get<0>(*itpt).begin();
                    KDTree::node_type::const_iterator en = boost::get<0>(*itpt).end();
                    //LIFE_ASSERT( it != en  )( i ).error( "invalid iterator");
                    if ( it == en )
                        continue;


                    for (size_type k=0; k < p.N; ++k)
                        {
                            // debugging output
#if 1
                            Debug( 4011 ) << "test: k=" << k << ", "
                                    << *it
                                    << ", p.bmin[k]=" << p.bmin[k]
                                    << ", p.bmax[k]=" << p.bmax[k]
                                    << ", isin = " << bool(( it[k] > p.bmin[k]) || ( it[k] < p.bmax[k]))
                                    << "\n";
#endif
                            if ( ( it[k] < p.bmin[k]) || ( it[k] > p.bmax[k]) )
                                {
                                    is_in = false;
                                    break;
                                }
                        }
                    if (is_in)
                        {
                            Debug( 4011 ) << "new point in box pt = " << boost::get<0>(*itpt) << " ptindex = " << boost::get<1>(*itpt)<< "\n";
                            p.ipts->push_back(*itpt);
                            Debug( 4011 ) << "size p.ipts = " << p.ipts->size() << "\n";
                        }
                }
        }
}



} // detail namespace

void
KDTree::clearTree()
{
    M_pts.resize( 0 );
    detail::destroy_tree( M_tree );
    M_tree = 0;
}

void
KDTree::pointsInBox( points_type &inpts,
                     const node_type &min,
                     const node_type &max)
{
    inpts.resize(0);

    // construct the tree from the points set
    if (M_tree == 0)
        {
            M_tree = detail::build_tree( M_pts.begin(),
                                         M_pts.end(),
                                         0 );
            if (!M_tree)
                return;
        }

    // check that we have a box
    node_type bmin(min);
    node_type bmax(max);

    for (size_type i=0; i < bmin.size(); ++i)
        if (bmin[i] > bmax[i])
            return;

    detail::points_in_box_data p;
    p.bmin = bmin.begin();
    p.bmax = bmax.begin();

    p.ipts = &inpts;
    p.N = boost::get<0>( *M_pts.begin() ).size();
    detail::points_in_box(p, M_tree, 0);

    Debug( 4011 ) << "size inpts = " << inpts.size() << "\n";

}


} //Life
