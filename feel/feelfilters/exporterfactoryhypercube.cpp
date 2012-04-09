/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Goncalo Pena <gpena@mat.uc.pt>
   Date: 2010-04-16

   Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file exporterfactory.cpp
   \author Goncalo Pena <gpena@mat.uc.pt>
   \date 2010-04-16
*/

#include <feel/feelfilters/exporterimpl.hpp>
#include <feel/feelfilters/exporterensight.hpp>
#include <feel/feelfilters/exportergmsh.hpp>

namespace Feel
{
/// \cond detail
namespace detail
{
template<typename MT>
Exporter<MT>* createEnsight()
{
    return new ExporterEnsight<MT>;
}
template<typename MT>
Exporter<MT>* createGmsh()
{
    return new ExporterGmsh<MT>;
}
} // detail

#if 0
//
// Hypercube 1,1
//
typedef Mesh<Hypercube<1,1> > meshsp11_t;
const bool meshsp11e = Exporter<meshsp11_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshsp11_t> );
const bool meshsp11g = Exporter<meshsp11_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshsp11_t> );

//
// Hypercube 2,1
//
typedef Mesh<Hypercube<2,1> > meshsp21_t;
const bool meshsp21e = Exporter<meshsp21_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshsp21_t> );
const bool meshsp21g = Exporter<meshsp21_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshsp21_t> );

//
// Hypercube 2,2
//
typedef Mesh<Hypercube<2,2> > meshsp22_t;
const bool meshsp22e = Exporter<meshsp22_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshsp22_t> );
const bool meshsp22g = Exporter<meshsp22_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshsp22_t> );

//
// Hypercube 2,3
//
typedef Mesh<Hypercube<2,3> > meshsp23_t;
const bool meshsp23e = Exporter<meshsp23_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshsp23_t> );
const bool meshsp23g = Exporter<meshsp23_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshsp23_t> );

//
// Hypercube 3,1
//
typedef Mesh<Hypercube<3,1> > meshsp31_t;
const bool meshsp31e = Exporter<meshsp31_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshsp31_t> );
const bool meshsp31g = Exporter<meshsp31_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshsp31_t> );

typedef Mesh<Hypercube<3,2> > meshsp32_t;
const bool meshsp32e = Exporter<meshsp32_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshsp32_t> );
const bool meshsp32g = Exporter<meshsp32_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshsp32_t> );

#endif
//#if defined( FEELPP_INSTANTIATION_MODE )
//
// explicit instances
//
template class Exporter<Mesh<Hypercube<1,1,1> > >;
template class Exporter<Mesh<Hypercube<2,1,2> > >;
template class Exporter<Mesh<Hypercube<3,1,3> > >;
template class Exporter<Mesh<Hypercube<3,2,3> > >;

template class Exporter<Mesh<Hypercube<2,2,2> > >;
template class Exporter<Mesh<Hypercube<2,3,2> > >;

/// \endcond detail
}
//#endif // FEELPP_INSTANTIATION_MODE
