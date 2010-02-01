/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-07

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-06-07
 */
#ifndef __KDTree_H
#define __KDTree_H 1

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifealg/glas.hpp>


namespace Life
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
    struct Element; struct Leaf; struct Node;

    typedef ublas::vector<double> node_type;
    //le dernier size_type du template correspond a l'indice global du noeud dans le maillage
    typedef boost::tuple<node_type, size_type, uint16_type, size_type> index_node_type;
    typedef std::vector<index_node_type> points_type;
    typedef points_type::iterator points_iterator;
    typedef points_type::const_iterator points_const_iterator;
    
    //ici, le double correspond à la distance avec le noeud que l'on recherche
    typedef boost::tuple<node_type, size_type, uint16_type,double > index_node_search_type;
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
        M_distanceMax(INT_MAX),
        M_nbPtMax(4)
    {}
    KDTree( KDTree const & tree )
        :
        M_tree( tree.M_tree ),
        M_pts( tree.M_pts ),
        M_node_search(tree.M_node_search),
        M_PtsNearest(tree.M_PtsNearest),
        M_distanceMax(tree.M_distanceMax),
        M_nbPtMax(tree.M_nbPtMax)
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
    void reserve(size_type n)
    {
        M_pts.reserve(n);
    }
    /**
     * insert a new point in the tree
     * @return  the index of the point
     */
    size_type addPoint( node_type const& n, size_type indice_global=0 )
    {
        size_type i = M_pts.size();
        addPointWithId(n,i,0,indice_global);
        return i;
    }
    /**
     * insert a new point, with an associated number.
     */
    void addPointWithId(const node_type& n, size_type i, uint16_type comp, size_type indice_global=0  )
    {
        if (M_tree)
            clearTree();
        M_pts.push_back( boost::make_tuple( n, i, comp,indice_global ) );
    }
    /**
     * fills ipts with the indexes of points in the box
     * [min,max]
     */
    void pointsInBox( points_type &inpts,
                      const node_type &min,
                      const node_type &max);
                      
    /**
     * recherche les points voisins de M_node_search dans le kd-tree
     */
    void search(const node_type & node_);
    
    /**
     * affiche le resultat de la recherche
     */
    void showResultSearch();

    //@}


private:

    /**
     * destroy the tree data structure
     */
    void clearTree();
    
    /**
     * Lance la recherche des voisins de M_node_search (recursif)
     */
    void run_search( Element * tree, uint iter);
    
    /**
     * Mise a jour de la liste des points les plus proches
     */
    void update_Pts_search(const index_node_type & p);

private:
    Element* M_tree;
    points_type M_pts;

    //le point dont on cherche les plus proches voisins
    node_type M_node_search;
    //vecteur des plus proches voisins
    points_search_type M_PtsNearest;
    //plus grande distance du vecteur des plus proches voisins
    double M_distanceMax;
    //le nbre max de pt voisin que l'on veut chercher
    uint M_nbPtMax;


};
} // Life
#endif /* __KDTree_H */
