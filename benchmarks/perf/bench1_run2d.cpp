/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-06-14

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file bench1_run2d.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-06-14
 */
#include <bench1_impl.hpp>

namespace Life
{
void
Bench1::run2d()
{
    using namespace Life;

    typedef Mesh<Simplex<2> > mesh_type;
    boost::shared_ptr<mesh_type> aMesh( new mesh_type );

    GmshTensorizedDomain<2,1,2,Simplex> td;
    td.setCharacteristicLength( meshSize );
    std::string fname = td.generate( "bench12d" );

    ImporterGmsh<mesh_type> import( fname );
    aMesh->accept( import );


    Log() << "run2d starts" << "\n";
    bench1<mesh_type, 1>( aMesh );
    bench1<mesh_type, 2>( aMesh );
    bench1<mesh_type, 5>( aMesh );
    //bench1<mesh_type, 8>( aMesh );
    Log() << "run2d ends" << "\n";

}
}
