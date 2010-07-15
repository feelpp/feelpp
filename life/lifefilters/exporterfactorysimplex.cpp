/* -*- mode: c++ -*-

   This file is part of the Life library

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

#include <life/lifefilters/exporterimpl.hpp>
#include <life/lifefilters/exporterensight.hpp>
#include <life/lifefilters/exportergmsh.hpp>

namespace Life
{
/// \cond detail
    namespace detail {
        template<typename MT>
        Exporter<MT>* createEnsight() { return new ExporterEnsight<MT>; }
        template<typename MT>
        Exporter<MT>* createGmsh() { return new ExporterGmsh<MT>; }
    } // detail

    //
    // Simplex 1,1
    //
    typedef Mesh<Simplex<1,1> > meshs11_t;
    const bool meshs11e = Exporter<meshs11_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshs11_t> );
    const bool meshs11g = Exporter<meshs11_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshs11_t> );

    typedef Mesh<Simplex<1,1,2> > meshs112_t;
    const bool meshs112e = Exporter<meshs112_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshs112_t> );
    const bool meshs112g = Exporter<meshs112_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshs112_t> );

    //
    // Simplex 2,1
    //
    typedef Mesh<Simplex<2,1> > meshs21_t;
    const bool meshs21e = Exporter<meshs21_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshs21_t> );
    const bool meshs21g = Exporter<meshs21_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshs21_t> );

    typedef Mesh<Simplex<2,1,3> > meshs213_t;
    const bool meshs213e = Exporter<meshs213_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshs213_t> );
    const bool meshs213g = Exporter<meshs213_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshs213_t> );

    //
    // Simplex 2,2
    //
    typedef Mesh<Simplex<2,2> > meshs22_t;
    const bool meshs22e = Exporter<meshs22_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshs22_t> );
    const bool meshs22g = Exporter<meshs22_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshs22_t> );

    //
    // Simplex 2,3
    //
    typedef Mesh<Simplex<2,3> > meshs23_t;
    const bool meshs23e = Exporter<meshs23_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshs23_t> );
    const bool meshs23g = Exporter<meshs23_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshs23_t> );


    //
    // Simplex 3,1
    //
    typedef Mesh<Simplex<3,1> > meshs31_t;
    const bool meshs31e = Exporter<meshs31_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshs31_t> );
    const bool meshs31g = Exporter<meshs31_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshs31_t> );

    //#if defined( LIFE_INSTANTIATION_MODE )
    //
    // explicit instances
    //
    template class Exporter<Mesh<Simplex<1,1,1> > >;
    template class Exporter<Mesh<Simplex<1,1,2> > >;
    template class Exporter<Mesh<Simplex<2,1,2> > >;
    template class Exporter<Mesh<Simplex<2,2,2> > >;
    template class Exporter<Mesh<Simplex<2,1,3> > >;
    template class Exporter<Mesh<Simplex<3,1,3> > >;

    template class Exporter<Mesh<Simplex<2,3,2> > >;
/// \endcond detail
}
//#endif // LIFE_INSTANTIATION_MODE
