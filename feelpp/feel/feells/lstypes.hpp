/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 07 juil. 2015
 
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
#ifndef FEELPP_LSTYPES_HPP
#define FEELPP_LSTYPES_HPP 1

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>



namespace Feel {


// Instantiate 2d and 3d only when compiling reinit_fms.cpp (where FEELPP_INSTANTIATE_FMS is defined), else "extern" avoid useless instantiation

template<int Dim, int GeoOrder=1,template<uint16_type,int,uint16_type> class Convex=Simplex>
using ls_mesh_type = Feel::Mesh< Convex<Dim,GeoOrder,Dim> >;
template<int Dim, int GeoOrder=1,template<uint16_type,int,uint16_type> class Convex=Simplex>
using ls_mesh_ptrtype = std::shared_ptr<Feel::Mesh< Convex<Dim,GeoOrder,Dim> >>;

template<int Dim,int GeoOrder=1,template<uint16_type,int,uint16_type> class Convex=Simplex>
using ls_space_type = Pch_type<ls_mesh_type<Dim,GeoOrder,Convex>,1>;
template<int Dim,int GeoOrder=1,template<uint16_type,int,uint16_type> class Convex=Simplex>
using ls_space_ptrtype = Pch_ptrtype<ls_mesh_type<Dim,GeoOrder,Convex>,1>;

template<int Dim,int GeoOrder=1,template<uint16_type,int,uint16_type> class Convex=Simplex>
using ls_element_type = typename ls_space_type<Dim,GeoOrder,Convex>::element_type;

template<int Dim, int GeoOrder=1>
using lsh_mesh_type = ls_mesh_type<Dim,GeoOrder,Hypercube>;
template<int Dim,int GeoOrder=1>
using lsh_space_type = ls_space_type<Dim,GeoOrder,Hypercube>;
template<int Dim,int GeoOrder=1>
using lsh_element_type = ls_element_type<Dim,GeoOrder,Hypercube>;


}
#endif
