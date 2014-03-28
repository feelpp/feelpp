/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-06-07

  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

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
   \file kdtree.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-06-07
 */
#ifndef __KDTree_H
#define __KDTree_H 1

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>


namespace Feel
{


/**
 * \class KDTree
 * \brief KDTree class
 *
 * In computer science, a kd-tree (short for k-dimensional tree) is a
 * space-partitioning data structure for organizing points in a k-dimensional
 * space. kd-trees are a useful data structure for several applications, such as
 * searches involving a multidimensional search key (e.g. range searches and
 * nearest neighbor searches). kd-trees are a special case of BSP trees.
 *
 * @author Christophe Prud'homme
 * @see RegionTree
 */
class KDTree
{
public:


    /** @name Typedefs
     */
    //@{
    struct Element;
    struct Leaf;
    struct Node;

    typedef ublas::vector<double> node_type;
    //in the following template, the last size_type corresponds to the global index of node in the mesh
    typedef boost::tuple<node_type, size_type, uint16_type, size_type> index_node_type;
    typedef std::vector<index_node_type> points_type;
    typedef points_type::iterator points_iterator;
    typedef points_type::const_iterator points_const_iterator;

    //here, the double corresponds at the distance with the node that is search
    typedef boost::tuple<node_type, size_type, uint16_type, size_type, double > index_node_search_type;
    typedef std::vector<index_node_search_type> points_search_type;
    typedef points_search_type::iterator points_search_iterator;
    typedef points_search_type::const_iterator points_search_const_iterator;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    KDTree()
        :
        M_tree( 0 ),
        M_pts(),
        M_node_search(),
        M_PtsNearest(),
        M_distanceMax( INT_MAX ),
        M_nbPtMax( 4 )
    {}
    KDTree( KDTree const & tree )
        :
        M_tree( tree.M_tree ),
        M_pts( tree.M_pts ),
        M_node_search( tree.M_node_search ),
        M_PtsNearest( tree.M_PtsNearest ),
        M_distanceMax( tree.M_distanceMax ),
        M_nbPtMax( tree.M_nbPtMax )
    {}

    ~KDTree()
    {
        // destroy the tree
        clearTree();
    }

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * get the number of points in the tree
     */
    size_type nPoints() const
    {
        return M_pts.size();
    }

    /**
     * get the points set
     */
    const points_type &points() const
    {
        return M_pts;
    }

    /**
     * get the points Near Neighbor set
     */
    const points_search_type &pointsNearNeighbor() const
    {
        return M_PtsNearest;
    }


    size_type nPtMaxNearNeighbor()
    {
        return M_nbPtMax;
    }


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * reset the tree, remove all points
     */
    void clear()
    {
        clearTree();
        M_pts.clear();
    }

    /**
     * reserve memory for the index/node pair
     */
    void reserve( size_type n )
    {
        M_pts.reserve( n );
    }

    /**
     * define the max number of point for the research( Default is 4 )
     */
    void nbNearNeighbor( size_type n )
    {
        M_nbPtMax=n;
    }

    /**
     * insert a new point in the tree
     * @return  the index of the point
     */
    size_type addPoint( node_type const& n, size_type indice_global=0 )
    {
        size_type i = M_pts.size();
        addPointWithId( n,i,0,indice_global );
        return i;
    }
    /**
     * insert a new point, with an associated number.
     */
    void addPointWithId( const node_type& n, size_type i, uint16_type comp, size_type indice_global=0  )
    {
        if ( M_tree )
            clearTree();

        M_pts.push_back( boost::make_tuple( n, i, comp,indice_global ) );
    }
    /**
     * fills ipts with the indexes of points in the box
     * [min,max]
     */
    void pointsInBox( points_type &inpts,
                      const node_type &min,
                      const node_type &max );

    /**
     * search the neighbors points of M_node_search in the kd-tree
     */
    void search( const node_type & node_ );

    /**
     * print the result of the research
     */
    void showResultSearch();

    /**
     * print in a file the scatterplot with the cuts of hyperplan (2d et 3d)
     */
    void writeLatexData( std::string __nameFile = "kdtreeData.tex" );
    //@}


private:

    /**
     * destroy the tree data structure
     */
    void clearTree();

    /**
     * Run the the research of M_node_search neighbors(recursive)
     */
    void run_search( Element * tree, uint16_type iter );

    /**
     * Updating the list of nearest points
     */
    void update_Pts_search( const index_node_type & p );

private:
    Element* M_tree;
    points_type M_pts;

    //the point which we search the nearest neighbors
    node_type M_node_search;
    //vector of the closest neighbors
    points_search_type M_PtsNearest;
    //greater distance from the vector of nearest neighbors
    double M_distanceMax;
    //the maximum number of neighbors points that we want to search
    size_type M_nbPtMax;


};
} // Feel
#endif /* __KDTree_H */
