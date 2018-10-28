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
        ( "upload", po::value<std::string>(), "upload desc" )
        ( "download", po::value<std::string>(), "download desc" )
        ( "data", po::value<std::string>(), "specify the datas to upload or the download directory" )
        ( "contents", po::value<std::string>(), "contents desc" )
		;

    fs::path initialCurrentPath = fs::current_path();

    Environment env( _argc=argc, _argv=argv,
                     _desc=rdoptions,
                     _about=about( _name="remotedata" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ),
                     _subdir=false );

    if ( Environment::vm().count("upload") )
    {
        RemoteData rd( soption(_name="upload") );
        if ( !rd.canUpload() )
        {
            Feel::cout << "invalid upload\n";
            return 0;
        }
        if ( !Environment::vm().count("data") )
        {
            Feel::cout << "no data to upload\n";
            return 0;
        }
        std::string data = soption(_name="data");
        if ( fs::path(data).is_relative() )
            data = (initialCurrentPath/fs::path(data)).string();
        rd.upload( data );
    }
    else if ( Environment::vm().count("download") )
    {
        RemoteData rd( soption(_name="download") );
        if ( !rd.canDownload() )
        {
            Feel::cout << "invalid download\n";
            return 0;
        }
        std::string dir = Environment::downloadsRepository();
        if ( Environment::vm().count("data") )
            dir = soption(_name="data");
        Feel::cout << "download data in : " << dir << "\n";
        rd.download( dir );
    }
    else if ( Environment::vm().count("contents") )
    {
        RemoteData rd( soption(_name="contents") );
        auto res = rd.contents();
        for ( auto const& folderInfo : std::get<0>( res ) )
            std::cout << "-------------------------------------------------------\n"
                      << folderInfo->print().str() << "\n";
        for ( auto const& itemInfo : std::get<1>( res ) )
            std::cout << "-------------------------------------------------------\n"
                      << itemInfo->print().str() << "\n";
        for ( auto const& fileInfo : std::get<2>( res ) )
            std::cout << "-------------------------------------------------------\n"
                      << fileInfo->print().str() << "\n";
    }
    return 0;
}
