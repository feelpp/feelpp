//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
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
//! @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
//! @date 2023-10-30
//! @copyright 2023 Feel++ Consortium
//! @copyright 2023 Université de Strasbourg
//!
#pragma once

namespace Feel
{
    template <typename GeoShape, typename T = double, int Tag = 0, typename IndexT = uint32_type, bool _EnableSharedFromThis=EnableSharedFromThis>
    class Mesh;

    template <typename GeoShape, typename T = double, typename IndexT = uint32_type>
    class MeshStructured;

    template<typename MeshType>
    class MeshSupport;
} // namespace Feel
