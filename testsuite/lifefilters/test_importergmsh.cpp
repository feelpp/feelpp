/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-16

  Copyright (C) 2007-2009 Université Joseph Fourier (Grenoble I)

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
/**
   \file test_importergmsh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-06-16
 */
#include <sstream>

#include <life/lifediscr/mesh.hpp>
#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>

using namespace Life;

template<uint16_type Dim, template<uint16_type,uint16_type,uint16_type> class Entity>
struct TestImporterGmsh
{
    typedef Entity<Dim, 1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;
public:
    TestImporterGmsh()
        :
        M_mesh()
    {}
    void test( double hsize, int version  = 1 )
    {
        M_mesh = mesh_ptr_type( new mesh_type );
        Debug() << "testing ImporterGmsh with file format version " << version << "\n";
        std::string fname;
        GmshTensorizedDomain<entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,Entity> td;
        td.setVersion( version );
        td.setCharacteristicLength( hsize );
        fname = td.generate( entity_type::name().c_str() );
        ImporterGmsh<mesh_type> import( fname );
        std::ostringstream ostr;
        ostr << version << ".0";
        import.setVersion( ostr.str() );
        M_mesh->accept( import );
        Debug() << "testing ImporterGmsh with file format version " << version << " done\n";
    }
private:
    mesh_ptr_type M_mesh;
};

int main( int argc, char** argv )
{
    /* assertions handling */
    Life::Assert::setLog( "test_importergmsh.assert");

    double h = 1.0;
    if ( argc == 2 )
        h = ::atof( argv[1] );

    //TestImporterGmsh<1,Simplex> test_simplex_1;
    TestImporterGmsh<2,Simplex> test_simplex_2;
    test_simplex_2.test( h, 1 );
    test_simplex_2.test( h, 2 );
    //TestImporterGmsh<2,Simplex> test_simplex_2;
}
