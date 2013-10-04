/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-09

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
   \file meshbase.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-09
 */
#include <feel/feelmesh/meshbase.hpp>

namespace Feel
{

MeshBase::MeshBase( WorldComm const& worldComm )
    :
    M_components( MESH_ALL_COMPONENTS ),
    M_is_updated( false ),
    M_is_parametric( false ),
    M_n_vertices( 0 ),
    M_n_parts( 1 ),
    M_worldComm( worldComm )
{
    DVLOG(2) << "[MeshBase] constructor...\n";
}

MeshBase::MeshBase( MeshBase const& m )
    :
    M_components( m.M_components ),
    M_is_updated( m.M_is_updated ),
    M_is_parametric( m.M_is_parametric ),
    M_n_vertices( m.M_n_vertices ),
    M_n_parts( m.M_n_parts ),
    M_worldComm( m.M_worldComm )
{}

MeshBase::~MeshBase()
{}

MeshBase&
MeshBase::operator=( MeshBase const& m )
{
    if ( this != &m )
    {
        M_components = m.M_components;
        M_is_updated = m.M_is_updated;
        M_is_parametric = m.M_is_parametric;
        M_n_vertices = m.M_n_vertices;
        M_n_parts = m.M_n_parts;
        M_worldComm = m.M_worldComm;
    }

    return *this;
}

void
MeshBase::clear()
{
    M_is_updated = false;

    M_n_vertices = 0;

    // Reset the number of partitions
    M_n_parts = 1;

    M_components = MESH_ALL_COMPONENTS;
}
void
MeshBase::updateForUse( size_type components )
{
    this->setComponents( components );
    this->updateForUse();
}

bool
MeshBase::isPartitioned() const
{
    if ( mpi::environment::initialized() )
        return M_n_parts == M_worldComm.localSize();

    else
        return M_n_parts == 1;
}
} // Feel
