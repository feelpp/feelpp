/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2018-06-19

  Copyright (C) 2018 Feel++ Consortium

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

#include <feel/feelcore/remotedata.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;

	po::options_description rdoptions( "Remote data options" );
	rdoptions.add_options()
        ( "upload", po::value<std::string>(), "mesh dimension" )
        ( "data", po::value<std::string>(), "mesh dimension" )
		;

    Environment env( _argc=argc, _argv=argv,
                     _desc=rdoptions,
                     _about=about( _name="remotedata" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ),
                     _directory=".",_subdir=false );

    if ( Environment::vm().count("upload") )
    {
        RemoteData rd( soption(_name="upload") );
        if ( !rd.canUpload() )
            return 0;
        if ( !Environment::vm().count("data") )
        {
            Feel::cout << "no data to upload\n";
            return 0;
        }
        std::string data = soption(_name="data");
        rd.upload( data );
    }

    return 0;
}
