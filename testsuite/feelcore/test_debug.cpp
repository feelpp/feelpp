/*
  This file is part of the Feel library.

  Author: Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>

  Copyright (C) 2004 EPFL

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
#include <iostream>
#include <feel/feelcore/feel.hpp>

template<typename streamT>
void
showMe( streamT out )
{
    out << "showMe:: Hello World\n";
    out.flush();
}
int
main()
{

    Feel::Debug() << "Hello\n";


    //
    // To see this message setup the DEBUG variable to 2000
    // export DEBUG=2000
    // then execute test_debug
    //
    Feel::Debug( 2000 ) << "AREA 2000 is now enabled\n";

    showMe( Feel::Debug( 2000 ) );
}

