/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-06-07

  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)

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
   \file kdtree.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-06-07
 */
#include <feel/feelmesh/kdtree.hpp>
#include <fstream>

namespace Feel
{
/// \cond detail
struct KDTree::Element
{
    enum { POINTS_PER_LEAF = 16 };

    /* 0 => is a tree node, != 0 => tree leaf storing n points */
    size_type n;

    virtual ~Element() {}
    bool isleaf() const
    {
        return ( n != 0 );
    }
    Element( size_type n_ ) : n( n_ ) {}
};

struct KDTree::Leaf : public KDTree::Element
{
    KDTree::points_const_iterator it;
    Leaf( KDTree::points_const_iterator begin,
          KDTree::points_const_iterator end )
        :
        Element( std::distance( begin,end ) )
    {
        FEELPP_ASSERT( begin != end )( std::distance( begin,end ) ).error( "invalid iterator" );
        it = begin;
    }
};

struct KDTree::Node : public KDTree::Element
{
    scalar_type split_v;

    //pts decrivant le decoupage du noeud
    node_type ptmin,ptmax;

    /* left: <=v, right: >v */
    Element *left, *right;

    Node( scalar_type v, node_type min, node_type max, Element *left_, Element *right_ )
        :
        Element( 0 ),
        split_v( v ),
        ptmin( min ),
        ptmax( max ),
        left( left_ ),
        right( right_ )
    {}
};

namespace detail
{
/**
 *sorting of points_type with respect to a component of the points
 */
struct component_sort
{
    component_sort( unsigned d )
        :
        dir( d )
    {}
    bool operator()( KDTree::index_node_type const& a,
                     KDTree::index_node_type const& b )
    {
        return boost::get<0>( a )( dir ) < boost::get<0>( b )( dir );
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
        while ( ( begin <= end ) &&
                ( boost::get<0>( *begin )( dir ) <= median ) )
            ++begin;

        while ( ( begin <= end ) &&
                ( boost::get<0>( *end )( dir ) > median ) )
            --end;

        if ( begin < end )
        {
            //begin->swap(*end);
            swap( *begin, *end );
            ++begin;
            --end;
        }

        else break;
    }
    while ( 1 );

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
    if ( begin == end )
        return 0;

    size_type npts = std::distance( begin,end );
    DVLOG(2) << "[KDTree::build_tree] nbpts = " << npts << "\n";

    if ( npts > KDTree::Element::POINTS_PER_LEAF )
    {
        DVLOG(2) << "[KDTree::build_tree] split\n";
        KDTree::points_iterator itmedian;
        scalar_type median;
        size_type N = boost::get<0>( *begin ).size();
        KDTree::node_type ptmin( N ), ptmax( N );
        unsigned ndir_tests = dir/N;
        dir %= N;

        if ( npts > 50 )
        {
            /* too much points for an exact median: estimation of the median .. */
            std::vector<KDTree::index_node_type> v( 30 );

            //random_sample(begin,end,v.begin(),v.end());
            for ( size_type i=0; i < v.size(); ++i )
                v[i] = begin[rand() % npts];

            std::sort( v.begin(), v.end(), component_sort( dir ) );
            median = ( boost::get<0>( v[v.size()/2-1] )( dir )+
                       boost::get<0>( v[v.size()/2] )( dir ) )/2;
            itmedian = partition( begin,end,dir,median );
            ptmin( dir )=median;
            ptmax( dir )=median;

            if ( N==2 )
            {
                //allows to confine the median plane
                std::vector<KDTree::index_node_type> v2( 30 );

                for ( size_type i=0; i < v2.size(); ++i )
                    v2[i] = begin[rand() % npts];

                std::sort( v2.begin(), v2.end(), component_sort( ( dir+1 )%N ) );
                ptmin( ( dir+1 )%N )= boost::get<0>( *( v2.begin() ) ) ( ( dir+1 )%N ) ;
                ptmax( ( dir+1 )%N )= boost::get<0>( *( --v2.end() ) ) ( ( dir+1 )%N ) ;
            }

            else if ( N==3 )
            {
                //allows to confine the median plane
                std::vector<KDTree::index_node_type> v2( 30 );

                for ( size_type i=0; i < v2.size(); ++i )
                    v2[i] = begin[rand() % npts];

                std::sort( v2.begin(), v2.end(), component_sort( ( dir+1 )%N ) );
                ptmin( ( dir+1 )%N )= boost::get<0>( *( v2.begin() ) ) ( ( dir+1 )%N ) ;
                ptmax( ( dir+1 )%N )= boost::get<0>( *( --v2.end() ) ) ( ( dir+1 )%N ) ;

                std::sort( v2.begin(), v2.end(), component_sort( ( dir+2 )%N ) );
                ptmin( ( dir+2 )%N )= boost::get<0>( *( v2.begin() ) ) ( ( dir+2 )%N ) ;
                ptmax( ( dir+2 )%N )= boost::get<0>( *( --v2.end() ) ) ( ( dir+2 )%N ) ;
            }
        }

        else
        {
            /* exact median computation */
            std::sort( begin, end, component_sort( dir ) );
            itmedian = begin + npts/2 - 1;
            median = boost::get<0>( *itmedian )[dir];

            while ( itmedian < end &&
                    boost::get<0>( *itmedian )[dir] == median )
                itmedian++;

            ptmin( dir )=median;
            ptmax( dir )=median;

            if ( N==2 )
            {
                //allows to confine the median plane
                std::vector<KDTree::index_node_type> v2( npts );

                for ( size_type i=0; i < v2.size(); ++i )
                    v2[i] = begin[rand() % npts];

                std::sort( v2.begin(), v2.end(), component_sort( ( dir+1 )%N ) );
                ptmin( ( dir+1 )%N )= boost::get<0>( *( v2.begin() ) ) ( ( dir+1 )%N ) ;
                ptmax( ( dir+1 )%N )= boost::get<0>( *( --v2.end() ) ) ( ( dir+1 )%N ) ;
            }

            else if ( N==3 )
            {
                //allows to confine the median plane
                std::vector<KDTree::index_node_type> v2( 30 );

                for ( size_type i=0; i < v2.size(); ++i )
                    v2[i] = begin[rand() % npts];

                std::sort( v2.begin(), v2.end(), component_sort( ( dir+1 )%N ) );
                ptmin( ( dir+1 )%N )= boost::get<0>( *( v2.begin() ) ) ( ( dir+1 )%N ) ;
                ptmax( ( dir+1 )%N )= boost::get<0>( *( --v2.end() ) ) ( ( dir+1 )%N ) ;

                std::sort( v2.begin(), v2.end(), component_sort( ( dir+2 )%N ) );
                ptmin( ( dir+2 )%N )= boost::get<0>( *( v2.begin() ) ) ( ( dir+2 )%N ) ;
                ptmax( ( dir+2 )%N )= boost::get<0>( *( --v2.end() ) ) ( ( dir+2 )%N ) ;
            }
        }

        /* could not split the set (all points have same value for component 'dir' !) */
        if ( itmedian == end )
        {
            /* tested all N direction ? so all points are strictly the same */
            if ( ndir_tests == N-1 )
                return new KDTree::Leaf( begin,end );

            else
                return new KDTree::Node( median, ptmin, ptmax,
                                         build_tree( begin,
                                                     itmedian,
                                                     ( dir+1 )%N + ( ndir_tests+1 )*N ),
                                         0 );
        }

        else
        {
            /* the general case */
            FEELPP_ASSERT( ( boost::get<0>( *itmedian )[dir] > median &&
                             boost::get<0>( *( itmedian-1 ) )[dir] <= median ) )
            ( dir )
            ( boost::get<0>( *itmedian )[dir] )
            ( boost::get<0>( *( itmedian-1 ) )[dir] )
            ( median ).error( "invalid median iterator" );

            return new KDTree::Node( median, ptmin, ptmax,
                                     build_tree( begin, itmedian, ( dir+1 )%N ),
                                     build_tree( itmedian,end, ( dir+1 )%N ) );
        }
    }

    else
    {
        DVLOG(2) << "[KDTree::build_tree] return leaf\n";
        return new KDTree::Leaf( begin,end );
    }
}

static void destroy_tree( KDTree::Element *t )
{
    if ( t == 0 )
        return;

    if ( !t->isleaf() )
    {
        KDTree::Node *tn = static_cast<KDTree::Node*>( t );
        destroy_tree( tn->right );
        destroy_tree( tn->left );
        delete tn;
    }

    else
    {
        KDTree::Leaf *tl = static_cast<KDTree::Leaf*>( t );
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
points_in_box( const points_in_box_data& p,
               const KDTree::Element *t,
               unsigned dir )
{
    if ( !t->isleaf() )
    {
        const KDTree::Node *tn = static_cast<const KDTree::Node*>( t );

        if ( p.bmin[dir] <= tn->split_v && tn->left )
        {
            points_in_box( p,
                           tn->left,
                           ( dir+1 )%p.N );
        }

        if ( p.bmax[dir] > tn->split_v && tn->right )
        {
            points_in_box( p,
                           tn->right,
                           ( dir+1 )%p.N );
        }
    }

    else
    {
        const KDTree::Leaf *tl = static_cast<const KDTree::Leaf*>( t );
        KDTree::points_const_iterator itpt = tl->it;

        for ( size_type i=tl->n; i; --i, ++itpt )
        {
            DVLOG(2) <<  "i = " << i << " pt = " << boost::get<0>( *itpt ) << " ptindex = " << boost::get<1>( *itpt )<< "\n";
            bool is_in = true;
            KDTree::node_type::const_iterator it = boost::get<0>( *itpt ).begin();
            KDTree::node_type::const_iterator en = boost::get<0>( *itpt ).end();

            //FEELPP_ASSERT( it != en  )( i ).error( "invalid iterator");
            if ( it == en )
                continue;


            for ( size_type k=0; k < p.N; ++k )
            {
                // debugging output
#if 1
                DVLOG(2) << "test: k=" << k << ", "
                              << *it
                              << ", p.bmin[k]=" << p.bmin[k]
                              << ", p.bmax[k]=" << p.bmax[k]
                              << ", isin = " << bool( ( it[k] > p.bmin[k] ) || ( it[k] < p.bmax[k] ) )
                              << "\n";
#endif

                if ( ( it[k] < p.bmin[k] ) || ( it[k] > p.bmax[k] ) )
                {
                    is_in = false;
                    break;
                }
            }

            if ( is_in )
            {
                DVLOG(2) << "new point in box pt = " << boost::get<0>( *itpt ) << " ptindex = " << boost::get<1>( *itpt )<< "\n";
                p.ipts->push_back( *itpt );
                DVLOG(2) << "size p.ipts = " << p.ipts->size() << "\n";
            }
        }
    }
}

/*
 * Compute the square of the distance between two nodes
 */
double
distanceNodes( const KDTree::node_type & p1, const KDTree::node_type & p2 )
{

    double res=0.0;

    for ( uint16_type i=0; i<p1.size(); ++i )
    {
        res+=( p2( i )-p1( i ) )*( p2( i )-p1( i ) );
    }

    //return math::sqrt(res);
    return res;
}


double
restriction( const double & min, const double & max, const double & x )
{
    if ( x < min )
    {
        return min;
    }

    else if ( x > max )
    {
        return max;
    }

    else return x;
}

/**
* Compute the orthogonal project of a node \c node_search on the segment [p1,p2]
* if the projection is not in the segment, the closest point in the segment is
* returned
*/

void
projection2d( KDTree::node_type & res,
              const KDTree::node_type & pMin,
              const KDTree::node_type & pMax,
              const KDTree::node_type & node_search )
{
    // computation of the projection
    if ( pMin( 0 )==pMax( 0 ) )
    {
        res( 0 )=pMin( 0 );
        res( 1 )=restriction( pMin( 1 ),pMax( 1 ),node_search( 1 ) );
    }

    else
    {
        res( 1 )=pMin( 1 );
        res( 0 )=restriction( pMin( 0 ),pMax( 0 ),node_search( 0 ) );
    }

    //return res;

}

void
projection3d( KDTree::node_type & res,
              const KDTree::node_type & pMin,
              const KDTree::node_type & pMax,
              const KDTree::node_type & node_search )
{
    // computation of the projection
    if ( pMin( 0 )==pMax( 0 ) )
    {
        res( 0 )=pMin( 0 );
        res( 1 )=restriction( pMin( 1 ),pMax( 1 ),node_search( 1 ) );
        res( 2 )=restriction( pMin( 2 ),pMax( 2 ),node_search( 2 ) );
    }

    else if ( pMin( 1 )==pMax( 1 ) )
    {
        res( 0 )=restriction( pMin( 0 ),pMax( 0 ),node_search( 0 ) );
        res( 1 )=pMin( 1 );
        res( 2 )=restriction( pMin( 2 ),pMax( 2 ),node_search( 2 ) );
    }

    else if ( pMin( 2 )==pMax( 2 ) )
    {
        res( 0 )=restriction( pMin( 0 ),pMax( 0 ),node_search( 0 ) );
        res( 1 )=restriction( pMin( 1 ),pMax( 1 ),node_search( 1 ) );
        res( 2 )=pMin( 2 );
    }

}

/**
 * Recursive function wich wrote in latex format the hyperplan associated with a node
 */
void
writeDecompositionData( KDTree::Element * tree, boost::shared_ptr<std::ostringstream> __ostr )
{
    if ( ! tree->isleaf() )
    {
        const KDTree::Node *tn = static_cast<const KDTree::Node*>( tree );

        if ( tn->ptmin.size()==2 ) //2d
        {
            *__ostr << "\\draw[-,gray!100,dashed] ("
                    << tn->ptmin( 0 ) << "," << tn->ptmin( 1 )
                    << ") -- ("
                    << tn->ptmax( 0 ) << "," << tn->ptmax( 1 )
                    << ");" << "\n";
        }

        else if ( tn->ptmin.size()==3 ) //3d
        {
            if ( tn->ptmin( 2 )==tn->ptmax( 2 ) )
                *__ostr << "\\draw[-,gray!100,fill=gray!20,opacity=0.5] ("
                        << tn->ptmin( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmin( 2 )
                        << ") -- ("
                        << tn->ptmin( 0 ) << "," << tn->ptmax( 1 ) << "," << tn->ptmin( 2 )
                        << ") -- ("
                        << tn->ptmax( 0 ) << "," << tn->ptmax( 1 ) << "," << tn->ptmax( 2 )
                        << ") -- ("
                        << tn->ptmax( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmin( 2 )
                        << ") -- ("
                        << tn->ptmin( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmin( 2 )
                        << ");" << "\n";

            else if ( tn->ptmin( 1 )==tn->ptmax( 1 ) )
                *__ostr << "\\draw[-,gray!100,fill=gray!20,opacity=0.5] ("
                        << tn->ptmin( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmin( 2 )
                        << ") -- ("
                        << tn->ptmin( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmax( 2 )
                        << ") -- ("
                        << tn->ptmax( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmax( 2 )
                        << ") -- ("
                        << tn->ptmax( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmin( 2 )
                        << ") -- ("
                        << tn->ptmin( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmin( 2 )
                        << ");" << "\n";

            else if ( tn->ptmin( 0 )==tn->ptmax( 0 ) )
                *__ostr << "\\draw[-,gray!100,fill=gray!20,opacity=0.5] ("
                        << tn->ptmin( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmin( 2 )
                        << ") -- ("
                        << tn->ptmin( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmax( 2 )
                        << ") -- ("
                        << tn->ptmin( 0 ) << "," << tn->ptmax( 1 ) << "," << tn->ptmax( 2 )
                        << ") -- ("
                        << tn->ptmin( 0 ) << "," << tn->ptmax( 1 ) << "," << tn->ptmin( 2 )
                        << ") -- ("
                        << tn->ptmin( 0 ) << "," << tn->ptmin( 1 ) << "," << tn->ptmin( 2 )
                        << ");" << "\n";
        }

        if ( tn->right ) writeDecompositionData( tn->right,__ostr );

        if ( tn->left ) writeDecompositionData( tn->left,__ostr );
    }
}


} // detail namespace
/// \endcond detail

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
                     const node_type &max )
{
    inpts.resize( 0 );

    // construct the tree from the points set
    if ( M_tree == 0 )
    {
        M_tree = detail::build_tree( M_pts.begin(),
                                     M_pts.end(),
                                     0 );

        if ( !M_tree )
            return;
    }

    // check that we have a box
    node_type bmin( min );
    node_type bmax( max );

    for ( size_type i=0; i < bmin.size(); ++i )
        if ( bmin[i] > bmax[i] )
            return;

    detail::points_in_box_data p;
    p.bmin = bmin.begin();
    p.bmax = bmax.begin();

    p.ipts = &inpts;
    p.N = boost::get<0>( *M_pts.begin() ).size();
    detail::points_in_box( p, M_tree, 0 );

    DVLOG(2) << "size inpts = " << inpts.size() << "\n";

}


void
KDTree::search( const node_type & node_ )
{

    M_node_search = node_;

    // construct the tree from the points set
    if ( M_tree == 0 )
    {
        M_tree = detail::build_tree( M_pts.begin(),
                                     M_pts.end(),
                                     0 );

        if ( !M_tree )
            return;
    }

    // clean the old research
    M_PtsNearest.clear();
    M_distanceMax = INT_MAX;

    // run search
    run_search( M_tree,0 );

}

void
KDTree::showResultSearch()
{

    std::cout<<"\n["<<M_PtsNearest.size()<<"]\n";

    points_search_const_iterator itpts=M_PtsNearest.begin();
    points_search_const_iterator itpts_end=M_PtsNearest.end();

    for ( ; itpts<itpts_end; ++itpts )
    {
        std::cout<<"("
                 <<boost::get<0>( *itpts )<<" "
                 <<boost::get<1>( *itpts )<<" "
                 <<boost::get<2>( *itpts )<<" "
                 <<boost::get<3>( *itpts )<<" "
                 <<boost::get<4>( *itpts )
                 <<")\n";
    }

}

void
KDTree::run_search( KDTree::Element * tree, uint16_type iter )
{

    if ( ! tree->isleaf() )
    {
        bool aGauche=false;
        const KDTree::Node *tn = static_cast<const KDTree::Node*>( tree );
        size_type N = tn->ptmin.size();

        if ( ( iter%N )==0 )
        {
            if ( M_node_search( 0 )<tn->split_v && tn->left )
            {
                run_search( tn->left,iter+1 );
                aGauche=true;
            }

            else if ( tn->right )
            {
                run_search( tn->right,iter+1 );
                aGauche=false;
            }
        }

        else if ( ( iter%N )==1 )
        {
            if ( M_node_search( 1 )<tn->split_v && tn->left )
            {
                run_search( tn->left,iter+1 );
                aGauche=true;
            }

            else if ( tn->right )
            {
                run_search( tn->right,iter+1 );
                aGauche=false;
            }
        }

        else if ( ( iter%N )==2 )
        {
            if ( M_node_search( 2 )<tn->split_v && tn->left )
            {
                run_search( tn->left,iter+1 );
                aGauche=true;
            }

            else if ( tn->right )
            {
                run_search( tn->right,iter+1 );
                aGauche=false;
            }
        }

        //compute the distance from the point of research with the border cutting kd-tree
        KDTree::node_type proj( N );

        switch ( N )
        {
        case 1:
        {
            proj=tn->ptmin;
            break;
        }

        case 2:
        {
            detail::projection2d( proj,tn->ptmin, tn->ptmax,M_node_search );
            break;
        }

        case 3:
        {
            detail::projection3d( proj,tn->ptmin, tn->ptmax,M_node_search );
            break;
        }
        }

        //if the border is close enough runs on the other side of the graph
        if ( detail::distanceNodes( proj,M_node_search )<( M_distanceMax ) )
        {
            if ( aGauche )
            {
                if ( tn->right ) run_search( tn->right,iter+1 );
            }

            else if ( tn->left ) run_search( tn->left,iter+1 );
        }

    }

    else
    {
        const KDTree::Leaf *tl = static_cast<const KDTree::Leaf*>( tree );

        KDTree::points_const_iterator itpt = tl->it;

        for ( size_type i=0; i<tl->n; ++i,++itpt )
            update_Pts_search( *itpt );

    }


}

void
KDTree::update_Pts_search( const index_node_type & p )
{

    points_search_iterator itpts;
    points_search_iterator itpts_end;

    double d=detail::distanceNodes( boost::get<0>( p ) , this->M_node_search );

    if ( M_PtsNearest.size()<this->M_nbPtMax )
    {
        index_node_search_type newEl( boost::make_tuple( boost::get<0>( p ),
                                      boost::get<1>( p ),
                                      boost::get<2>( p ),
                                      boost::get<3>( p ),
                                      d  )
                                    );
        itpts=M_PtsNearest.begin();
        itpts_end=M_PtsNearest.end();

        while ( itpts!=itpts_end && d > boost::get<4>( *itpts ) )
        {
            ++itpts;
        }

        M_PtsNearest.insert( itpts,newEl );
    }

    else if ( d<this->M_distanceMax )
    {
        itpts=M_PtsNearest.begin();
        itpts_end=M_PtsNearest.end();

        for ( ; itpts<itpts_end; ++itpts )
        {
            if ( d < boost::get<4>( *itpts ) )
            {
                index_node_search_type newEl( boost::make_tuple( boost::get<0>( p ),
                                              boost::get<1>( p ),
                                              boost::get<2>( p ),
                                              boost::get<3>( p ),
                                              d  )
                                            );
                M_PtsNearest.insert( itpts,newEl );
                M_PtsNearest.pop_back();
                M_distanceMax= boost::get<4>( M_PtsNearest.back() );
                itpts=itpts_end;
            }

        }

    }

}



void
KDTree::writeLatexData( std::string __nameFile )
{
    std::ofstream __file( __nameFile.c_str() );
    __file << "\\begin{figure}[h]" << "\n"
           << "\\begin{center}" << "\n"
           << "\\subfigure[Organisation des points de l'espace] {" << "\n"
           << "\\begin{tikzpicture}[scale=1]" << "\n";

    points_const_iterator itpts = M_pts.begin();
    points_const_iterator itpts_end = M_pts.end();
    uint16_type __dim = boost::get<0>( *itpts ).size();

    if ( __dim==2 )
        for ( ; itpts!=itpts_end ; ++itpts )
        {
            __file << "\\draw ("
                   << boost::get<0>( *itpts )( 0 ) << ","
                   << boost::get<0>( *itpts )( 1 ) << ") node {$\\bullet$};"
                   << "\n";
        }

    else if ( __dim==3 )
        for ( ; itpts!=itpts_end ; ++itpts )
        {
            __file << "\\draw ("
                   << boost::get<0>( *itpts )( 0 ) << ","
                   << boost::get<0>( *itpts )( 1 ) << ","
                   << boost::get<0>( *itpts )( 2 ) << ") node {$\\bullet$};"
                   << "\n";
        }


    // construct the tree from the points set
    if ( M_tree == 0 )
    {
        M_tree = detail::build_tree( M_pts.begin(),
                                     M_pts.end(),
                                     0 );

        if ( !M_tree )
            return;
    }

    boost::shared_ptr<std::ostringstream> __ostr( new std::ostringstream() );
    detail::writeDecompositionData( M_tree,__ostr );
    __file << __ostr->str();


    __file <<"\\end{tikzpicture}" << "\n"
           << "}" <<"\n"
           << "\\end{center}" << "\n"
           << "\\end{figure}" << "\n";


}


} //Feel
