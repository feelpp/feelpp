/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 13 Apr 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_MESHIO_HPP
#define FEELPP_MESHIO_HPP 1

#include <feel/feelfilters/partitionio.hpp>

#if defined( FEELPP_HAS_HDF5 )

namespace Feel
{

template <typename Shape, typename T, int Tag, typename IndexT, bool _EnableSharedFromThis>
void Mesh<Shape, T, Tag, IndexT, _EnableSharedFromThis>::ioHDF5( IOStatus status, std::string const& filename, size_type ctxMeshUpdate, double scale )
{
    PartitionIO<mesh_type> io( filename );
    if ( status == IOStatus::isLoading )
        io.read( this->shared_from_this(), ctxMeshUpdate, scale );
    else if ( status == IOStatus::isSaving )
        io.write( this->shared_from_this(), scale );
}

} // namespace Feel

#endif

#endif
