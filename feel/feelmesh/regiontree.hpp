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
   \file RegionTree.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-02-01
 */
#ifndef __RegionTree_H
#define __RegionTree_H 1

#include <deque>
#include <set>
#include <vector>

#include <boost/shared_ptr.hpp>


#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/glas.hpp>

namespace Feel
{
/*!
  \class RegionTree
  \brief implements a region-tree for point search in a set of boxes

  @author Christophe Prud'homme
*/
class RegionTree
{
public:


    /** @name Typedefs
     */
    //@{

    typedef node<double>::type node_type;

    struct box_index_type
    {
        size_type id;
        node_type min, max;
    };

    typedef std::deque<box_index_type> box_container_type;
    typedef std::vector<const box_index_type*> pbox_container_type;
    typedef std::set<const box_index_type*> pbox_set_type;

    struct element_base;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    RegionTree()
        :
        M_root( 0 )
    {}

    ~RegionTree()
    {
        destroy();
    }

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    size_type nbBoxes() const
    {
        return M_boxes.size();
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
     * add a box in the region tree
     * \param min min coordinates for bounding box
     * \param max max coordinates for bounding box
     * \param id id of the element stored in the bounding box
     */
    void addBox( node_type min, node_type max, size_type id=size_type( -1 ) );

    /**
       clear the tree
    */
    void clear()
    {
        destroy();
        M_boxes.clear();
    }


    void findIntersectingBoxes( const node_type& bmin, const node_type& bmax, pbox_set_type& boxlst );
    void findContainingBoxes( const node_type& bmin, const node_type& bmax, pbox_set_type& boxlst );
    void findContainedBoxes( const node_type& bmin, const node_type& bmax, pbox_set_type& boxlst );
    void findBoxesAtPoint( const node_type& P, pbox_set_type& boxlst );

    void findIntersectingBoxes( const node_type& bmin, const node_type& bmax, std::vector<size_type>& idvec );
    void findContainingBoxes( const node_type& bmin, const node_type& bmax, std::vector<size_type>& idvec );
    void findContainedBoxes( const node_type& bmin, const node_type& bmax, std::vector<size_type>& idvec );
    void findBoxesAtPoint( const node_type& P, std::vector<size_type>& idvec );
    void dump();

    //@}



protected:

private:




    void operator=( const RegionTree& ) {} /* non-copiable */
    RegionTree( const RegionTree& ) {} /*non-copiable */

    void build();
    void destroy();

    static void toIdVector( pbox_set_type const& bs, std::vector<size_type>& idvec );

    element_base* M_root;

    box_container_type M_boxes;
};

/**
   \typedef RegionTree region_tree_type;

   an alias for RegionTree
 */
typedef RegionTree region_tree_type;

/**
   \typedef boost::shared_ptr<RegionTree> region_tree_ptrtype;

   pointer type for Region Tree.
 */
typedef boost::shared_ptr<region_tree_type> region_tree_ptrtype;

}
#endif /* __RegionTree_H */
