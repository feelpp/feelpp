//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 16 Apr 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_MAKEMESH_HPP
#define FEELPP_MAKEMESH_HPP 1

#include <feel/feeldiscr/mesh.hpp>


namespace Feel {

//!
//! Create a shared pointer \p Mesh<T> from \p Args
//! @code
//! // create a tetrahedron mesh std::shared_ptr
//! auto mesh = makeSharedMesh<Simplex<3>>();
//! @endcode
//!
template< class T, class... Args >
std::shared_ptr<Mesh<T>>
makeSharedMesh( Args&&... args )
{
    return std::make_shared<Mesh<T>>( args... );
}

//!
//! @see makeSharedMesh
//!
template< class T, class... Args >
decltype(auto)
makeMesh( Args&&... args ) 
{
    return makeSharedMesh<T>( args... );
}

//!
//! Create a unique pointer \p Mesh<T> from \p Args
//! @code
//! // create a tetrahedron mesh std::unique_ptr
//! auto mesh = makeUniqueMesh<Simplex<3>>();
//! @endcode
//!
template< class T, class... Args >
std::unique_ptr<Mesh<T>>
makeUniqueMesh( Args&&... args )
{
    return std::make_unique<Mesh<T>>( args... );
}

} // Feel

#endif
