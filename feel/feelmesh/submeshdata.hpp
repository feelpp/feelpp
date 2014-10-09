/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-01-27

  Copyright (C) 2013 Université de Strasbourg

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
   \file submeshdata.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-01-27
 */
#ifndef __SubMeshData_H
#define __SubMeshData_H 1

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wredeclared-class-member"
#endif
#include <boost/bimap.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace Feel
{
class MeshBase;

/**
 * \class SubMeshData
 * \brief data structure storing sub mesh data
 *
 * \author Christophe Prud'homme
 * \see Mesh
 */
class SubMeshData
{
public:
    /** @name Typedefs
     */
    //@{
    typedef const MeshBase mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef boost::bimap< size_type, size_type > bm_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{
    template<typename MeshType>
    SubMeshData( boost::shared_ptr<MeshType> m ) : mesh( m )
        {}

    ~SubMeshData()
        {
            VLOG(2) << "delete sub mesh data\n";
        }

    //MeshBase* meshBase() { return dynamic_cast<MeshBase *>( mesh.get() ); }
    //MeshBase const* meshBase() const { return dynamic_cast<MeshBase const*>( mesh.get() ); }

    //@}


    // parent mesh
    mesh_ptrtype mesh;

    // bi-directional map
    bm_type bm;

};
} // Feel
#endif /* __SubMeshData_H */
