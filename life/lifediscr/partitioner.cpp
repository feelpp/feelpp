 /* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-31

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
   \file partitioner.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-31
 */
#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <life/lifediscr/partitioner.hpp>
#include <life/lifediscr/partitionermetis.hpp>
#include <life/lifediscr/partitionerparmetis.hpp>
#include <life/lifediscr/mesh.hpp>

namespace Life
{

namespace detail
{
template<typename MT>
Partitioner<MT>* createMetis()
{
    Debug( 4020 ) << "creating Metis partitioner for " << typeid( MT ).name() << "\n";
    return new PartitionerMetis<MT>;
}
template<typename MT>
Partitioner<MT>* createParmetis()
{
    Debug( 4020 ) << "creating Parmetis partitioner for " << typeid( MT ).name() << "\n";
    return new PartitionerParmetis<MT>;
}


//
// Simplex 1,1
//
typedef Mesh<Simplex<1,1,1> > meshs11_t;
const bool meshs11e = Partitioner<meshs11_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshs11_t> );
const bool meshs11g = Partitioner<meshs11_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshs11_t> );

typedef Mesh<Simplex<1,1,2> > meshs112_t;
const bool meshs112e = Partitioner<meshs112_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshs112_t> );
const bool meshs112g = Partitioner<meshs112_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshs112_t> );

//
// Simplex 2,1
//
typedef Mesh<Simplex<2,1,2> > meshs21_t;
const bool meshs21e = Partitioner<meshs21_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshs21_t> );
const bool meshs21g = Partitioner<meshs21_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshs21_t> );

typedef Mesh<Simplex<2,2,2> > meshs22_t;
const bool meshs22e = Partitioner<meshs22_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshs22_t> );
const bool meshs22g = Partitioner<meshs22_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshs22_t> );

typedef Mesh<Simplex<2,3,2> > meshs23_t;
const bool meshs23e = Partitioner<meshs23_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshs23_t> );
const bool meshs23g = Partitioner<meshs23_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshs23_t> );

typedef Mesh<Simplex<2,4,2> > meshs24_t;
const bool meshs24e = Partitioner<meshs24_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshs24_t> );
const bool meshs24g = Partitioner<meshs24_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshs24_t> );

typedef Mesh<Simplex<2,5,2> > meshs25_t;
const bool meshs25e = Partitioner<meshs25_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshs25_t> );
const bool meshs25g = Partitioner<meshs25_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshs25_t> );

//
// Simplex 2,1,3
//
typedef Mesh<Simplex<2,1,3> > meshs213_t;
const bool meshs213e = Partitioner<meshs213_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshs213_t> );
const bool meshs213g = Partitioner<meshs213_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshs213_t> );

//
// Simplex 3,1
//
typedef Mesh<Simplex<3,1,3> > meshs31_t;
const bool meshs31e = Partitioner<meshs31_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshs31_t> );
const bool meshs31g = Partitioner<meshs31_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshs31_t> );

typedef Mesh<Simplex<3,2,3> > meshs32_t;
const bool meshs32e = Partitioner<meshs32_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshs32_t> );
const bool meshs32g = Partitioner<meshs32_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshs32_t> );

//
// SimplexProduct 1,1
//
typedef Mesh<SimplexProduct<1,1,1> > meshsp11_t;
const bool meshsp11e = Partitioner<meshsp11_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshsp11_t> );
const bool meshsp11g = Partitioner<meshsp11_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshsp11_t> );

//
// SimplexProduct 2,1
//
typedef Mesh<SimplexProduct<2,1,2> > meshsp21_t;
const bool meshsp21e = Partitioner<meshsp21_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshsp21_t> );
const bool meshsp21g = Partitioner<meshsp21_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshsp21_t> );

//
// SimplexProduct 3,1
//
typedef Mesh<SimplexProduct<3,1,3> > meshsp31_t;
const bool meshsp31e = Partitioner<meshsp31_t>::Factory::type::instance().registerProduct( "metis", &detail::createMetis<meshsp31_t> );
const bool meshsp31g = Partitioner<meshsp31_t>::Factory::type::instance().registerProduct( "parmetis", &detail::createParmetis<meshsp31_t> );


} // detail




template<typename MeshType>
Partitioner<MeshType>::Partitioner()
{
#if defined(LIFE_INSTANTIATION_MODE)
    test_partitioner();
#endif
}



#if defined( LIFE_INSTANTIATION_MODE )
void test_partitioner()
{
    Debug( 4020 ) << "meshsp21e= " << detail::meshsp21e << "\n"
                  << "meshs31e = " << detail::meshs31e << "\n";
}

//
// explicit instances
//
// 1D
template class Partitioner<Mesh<Simplex<1,1,1> > >;
// 1D2D
template class Partitioner<Mesh<Simplex<1,1,2> > >;
// 2D
template class Partitioner<Mesh<Simplex<2,1,2> > >;
#if BOOST_PP_GREATER_EQUAL( LIFE_MESH_MAX_ORDER, 2 )
template class Partitioner<Mesh<Simplex<2,2,2> > >;
#endif
#if BOOST_PP_GREATER_EQUAL( LIFE_MESH_MAX_ORDER, 3 )
template class Partitioner<Mesh<Simplex<2,3,2> > >;
#endif
#if BOOST_PP_GREATER_EQUAL( LIFE_MESH_MAX_ORDER, 4 )
template class Partitioner<Mesh<Simplex<2,4,2> > >;
#endif
#if BOOST_PP_GREATER_EQUAL( LIFE_MESH_MAX_ORDER, 5 )
template class Partitioner<Mesh<Simplex<2,5,2> > >;
#endif
// 2D3D
template class Partitioner<Mesh<Simplex<2,1,3> > >;
// 3D
template class Partitioner<Mesh<Simplex<3,1,3> > >;
#if BOOST_PP_GREATER_EQUAL( LIFE_MESH_MAX_ORDER, 2 )
template class Partitioner<Mesh<Simplex<3,2,3> > >;
#endif
template class Partitioner<Mesh<SimplexProduct<1,1,1> > >;
template class Partitioner<Mesh<SimplexProduct<2,1,2> > >;
template class Partitioner<Mesh<SimplexProduct<3,1,3> > >;
#endif // LIFE_INSTANTIATION_MODE
}
