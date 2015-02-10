/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Goncalo Pena <gpena@mat.uc.pt>
   Date: 2010-04-16

   Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

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
#if 0
namespace detail
{
template<typename MT>
Exporter<MT,1>* createEnsight()
{
    return new ExporterEnsight<MT>;
}
template<typename MT, int N>
Exporter<MT,N>* createGmsh()
{
    return new ExporterGmsh<MT,N>;
}
} // detail


# define DIMS BOOST_PP_TUPLE_TO_LIST(3,(1,2,3))
# define ORDERS BOOST_PP_TUPLE_TO_LIST(5,(1,2,3,4,5))
# define ORDERS_FUN_ENSIGHT BOOST_PP_TUPLE_TO_LIST(1,(1))
# define ORDERS_FUN_GMSH BOOST_PP_TUPLE_TO_LIST(5,(1,2,3,4,5))

#define MeshName(name,LDIM,LORDER,ORDERFUN) BOOST_PP_CAT(BOOST_PP_CAT(BOOST_PP_CAT(name,LDIM),LORDER),ORDERFUN)
#define MeshType(name,LDIM,LORDER,ORDERFUN) BOOST_PP_CAT(BOOST_PP_CAT(BOOST_PP_CAT(BOOST_PP_CAT(name,LDIM),LORDER),ORDERFUN),_t)

// Ensight
# define FACTORY_ENSIGHT(LDIM,LORDER,ORDERFUN)                          \
    typedef Mesh<Simplex<LDIM,LORDER> > MeshType(mesh,LDIM,LORDER,ORDERFUN); \
    const bool MeshName(meshespensight,LDIM,LORDER,ORDERFUN) =          \
        Exporter<MeshType(mesh,LDIM,LORDER,ORDERFUN),ORDERFUN>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<MeshType(mesh,LDIM,LORDER,ORDERFUN),ORDERFUN> );


# define FACTORY_ENSIGHT_OP(_, GDO) FACTORY_ENSIGHT GDO

// Gmsh
# define FACTORY_GMSH(LDIM,LORDER,ORDERFUN)                             \
    typedef Mesh<Simplex<LDIM,LORDER> > MeshType(mesh,LDIM,LORDER,ORDERFUN); \
    const bool MeshName(meshexpgmsh,LDIM,LORDER,ORDERFUN) =             \
        Exporter<MeshType(mesh,LDIM,LORDER,ORDERFUN),ORDERFUN>::Factory::type::instance().registerProduct( "ensight", &detail::createGmsh<MeshType(mesh,LDIM,LORDER,ORDERFUN),ORDERFUN> );

# define FACTORY_GMSH_OP(_, GDO) FACTORY_GMSH GDO

// ensight
# define FACTORY(LDIM,LORDER,ORDERFUN) template class Exporter<Mesh<Simplex<LDIM,LORDER,LDIM> >, ORDERFUN >;
# define FACTORY_OP(_, GDO) FACTORY GDO

//BOOST_PP_LIST_FOR_EACH_PRODUCT(FACTORY_ENSIGHT_OP, 3, (DIMS, ORDERS,ORDERS_FUN_ENSIGHT))
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_GMSH_OP, 3, ( DIMS, ORDERS,ORDERS_FUN_GMSH ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_OP, 3, ( DIMS, ORDERS, ORDERS_FUN_GMSH ) )

#if 0
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

//
// Simplex 3,1
//
typedef Mesh<Simplex<3,2> > meshs32_t;
const bool meshs32e = Exporter<meshs32_t>::Factory::type::instance().registerProduct( "ensight", &detail::createEnsight<meshs32_t> );
const bool meshs32g = Exporter<meshs32_t>::Factory::type::instance().registerProduct( "gmsh", &detail::createGmsh<meshs32_t> );

//#if defined( FEELPP_INSTANTIATION_MODE )
//
// explicit instances
//
template class Exporter<Mesh<Simplex<1,1,1> > >;
template class Exporter<Mesh<Simplex<1,1,2> > >;
template class Exporter<Mesh<Simplex<2,1,2> > >;
template class Exporter<Mesh<Simplex<2,2,2> > >;
template class Exporter<Mesh<Simplex<2,1,3> > >;
template class Exporter<Mesh<Simplex<3,1,3> > >;
template class Exporter<Mesh<Simplex<3,2,3> > >;

template class Exporter<Mesh<Simplex<2,3,2> > >;
#endif
#endif
/// \endcond detail
}
//#endif // FEELPP_INSTANTIATION_MODE
