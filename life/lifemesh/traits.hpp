/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <prudhomm@zion>
       Date: 2005-08-27

  Copyright (C) 2005,2006 EPFL

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
   \file traits.hpp
   \author Christophe Prud'homme <prudhomm@zion>
   \date 2005-08-27
 */
#ifndef __LIFE_MESH_TRAITS_HPP
#define __LIFE_MESH_TRAITS_HPP 1

#include <life/lifecore/traits.hpp>

namespace Life
{
/**
 * \namespace Life::mesh
 *
 * The \c mesh namespace stores some mesh specific classes that may
 * conflict with classes in the namespace \c Life.
 */
namespace mesh
{
/**
 * \class Traits
 * \brief Traits for meshes
 *
 * @author Christophe Prud'homme
 */
template<typename MeshType>
struct Traits
{
    /** @name Typedefs
     */
    //@{

    typedef Traits<MeshType> self_type;
    typedef typename boost::remove_pointer<typename remove_shared_ptr<MeshType>::type >::type mesh_type;

    typedef typename mesh_type::shape_type element_shape_type;

    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::face_type face_type;

    // element iterators
    typedef typename mesh_type::element_iterator element_iterator;
    typedef typename mesh_type::element_const_iterator element_const_iterator;

    typedef typename mesh_type::marker_element_iterator marker_element_iterator;
    typedef typename mesh_type::marker_element_const_iterator marker_element_const_iterator;

    typedef typename mesh_type::marker2_element_iterator marker2_element_iterator;
    typedef typename mesh_type::marker2_element_const_iterator marker2_element_const_iterator;

    typedef typename mesh_type::marker3_element_iterator marker3_element_iterator;
    typedef typename mesh_type::marker3_element_const_iterator marker3_element_const_iterator;

    typedef typename mesh_type::location_element_iterator location_element_iterator;
    typedef typename mesh_type::location_element_const_iterator location_element_const_iterator;

    typedef typename mesh_type::pid_element_iterator pid_element_iterator;
    typedef typename mesh_type::pid_element_const_iterator pid_element_const_iterator;

    // face iterators
    typedef typename mesh_type::face_iterator face_iterator;
    typedef typename mesh_type::face_const_iterator face_const_iterator;

    typedef typename mesh_type::marker_face_iterator marker_face_iterator;
    typedef typename mesh_type::marker_face_const_iterator marker_face_const_iterator;

    typedef typename mesh_type::location_face_iterator location_face_iterator;
    typedef typename mesh_type::location_face_const_iterator location_face_const_iterator;

    // point iterators
    typedef typename mesh_type::point_iterator point_iterator;
    typedef typename mesh_type::point_const_iterator point_const_iterator;

    //@}
};




} // mesh
} // Life
#endif /* __Traits_H */
