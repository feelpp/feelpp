/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-02-14

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file adfstream.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */

#include <fstream>
#include <iomanip>
#include <iostream>

namespace Feel
{
class ADOfstream
    :
public std::ofstream
{

public:

    ADOfstream( const char* filename,
                const char* description,
                const char* sourceFile,
                const char* mnemonic )
        :
        std::ofstream( filename )
    {
        ( *this ) <<
                  "/* \n"
                  "   " << filename << "\t" << description << std::endl <<
                  "   This file is part of gstlibs.\n"
                  "\n"
                  "   Copyright (C) 2011 Christophe Prud'homme\n"
                  "\n"
                  "   gstlibs is free software; you can redistribute it and/or modify\n"
                  "   it under the terms of the GNU Lesser General Public License as published by\n"
                  "   the Free Software Foundation; either version 2 of the License, or\n"
                  "   (at your option) any later version.\n"
                  "   \n"
                  "   gstlibs is distributed in the hope that it will be useful,\n"
                  "   but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
                  "   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
                  "   GNU Lesser General Public License for more details.\n"
                  "   \n"
                  "   You should have received a copy of the GNU Lesser General Public License\n"
                  "   along with gstlibs; if not, write to the Free Software\n"
                  "   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n"
                  "*/"
                  << std::endl << std::endl
                  << "// Generated source file.  Do not edit. " << std::endl
                  << "// " << sourceFile << " " << __DATE__ << " " << __TIME__
                  << std::endl << std::endl
                  << "#ifndef " << mnemonic << std::endl
                  << "#define " << mnemonic << std::endl << std::endl;
    }

    void include( const char* filename )
    {
        ( *this ) << "#include <" << filename << ">" << std::endl;
    }

    void beginNamespace( std::string const& __ns )
    {
        ( *this ) << "namespace " << __ns << "\n"
                  << "{\n";
        //(*this) << "BZ_NAMESPACE(blitz)" << std::endl << std::endl;
    }
    void endNamespace()
    {
        ( *this ) << "}\n";
    }
    ~ADOfstream()
    {
        ( *this ) <<  std::endl << "#endif" << std::endl;

    }

};

}

