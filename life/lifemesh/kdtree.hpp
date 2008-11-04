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
    typedef boost::tuple<node_type, size_type, uint16_type > index_node_type;
    typedef std::vector<index_node_type> points_type;
    typedef points_type::iterator points_iterator;
    typedef points_type::const_iterator points_const_iterator;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    KDTree()
        :
        M_tree( 0 ),
        M_pts()
    {}
    KDTree( KDTree const & tree )
        :
        M_tree( tree.M_tree ),
        M_pts( tree.M_pts )
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
    size_type addPoint( node_type const& n )
    {
        size_type i = M_pts.size();
        addPointWithId(n,i,0);
        return i;
    }
    /**
     * insert a new point, with an associated number.
     */
    void addPointWithId(const node_type& n, size_type i, uint16_type comp  )
    {
        if (M_tree)
            clearTree();
        M_pts.push_back( boost::make_tuple( n, i, comp ) );
    }
    /**
     * fills ipts with the indexes of points in the box
     * [min,max]
     */
    void pointsInBox( points_type &inpts,
                      const node_type &min,
                      const node_type &max);

    //@}


private:

    /**
     * destroy the tree data structure
     */
    void clearTree();

private:
    Element* M_tree;
    points_type M_pts;


};
} // Life
#endif /* __KDTree_H */
