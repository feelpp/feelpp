/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date:27 Nov 2017

 Copyright (C) 2017 Feel++ Consortium

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

#if !defined(FEELPP_MESH_SUPPORT_BASE_HPP)
#define FEELPP_MESH_SUPPORT_BASE_HPP 1

#include <feel/feelcore/feel.hpp>
#include <unordered_set>

namespace Feel
{

/**
 * \brief Description of a mesh support.
 * A function space can use this object for built a space in a part of the mesh
 */
class MeshSupportBase
{
public :

    virtual ~MeshSupportBase() {}
    
    virtual bool isFullSupport() const = 0;
    virtual bool isPartialSupport() const = 0;

    virtual size_type numElements() const = 0;

    virtual bool hasElement( size_type eltId ) const = 0;

    virtual bool hasGhostElement( size_type eltId ) const = 0;

    virtual std::unordered_set<size_type> const& rangeMeshElementsIdsPartialSupport() const = 0;
    virtual std::unordered_set<size_type> const& rangeMeshElementsGhostIdsPartialSupport() const = 0;

    virtual void resetLocalizationTool() = 0;
};

} // namespace Feel

#endif
