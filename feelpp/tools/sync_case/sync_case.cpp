/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Thomas Saigre <saigre@math.unistra.fr>
 Date: 14 July 2022

 Copyright (C) 2022 Feel++ Consortium

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

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/repository.hpp>
#include <feel/feelcore/table.hpp>
#include <sys/wait.h>

using namespace Feel;

/**
 * Execute the rsync command to synchronize the case
 *
 * @param dir the directory to synchronize
 * @param login the remote user name to synchronize with
 * @param server the host to synchronize with
 * @param local_dir the local directory where file will be synchronized
 * 
 * @return the exit status of the rsync command
 */
int sync(const std::string dir, const std::string login, const std::string server, const std::string local_dir,
         int verbose = 0)
{
    LOG(INFO) << "Sync from " << login << "@" << server << ":" << dir << " to " << local_dir << std::endl; 

    std::string remote;
    const std::string opts = (verbose != 0) ? "-avzP" : "-azP";

    switch (fork())
    {
        case -1:
            Feel::cerr << "fork() failed" << std::endl;
            break;
        case 0:
            remote = login + "@" + server + ":" + dir;
            execlp("rsync", "rsync", opts.c_str(), remote.c_str(), local_dir.c_str(), NULL);
            exit(0);
            break;
        default:
            break;
    }

    int reason = 0, e = 0;
    wait( &reason );
    if ( WIFEXITED(reason) )
    {
        e = WEXITSTATUS( reason );
        if ( e != 0 )
        {
            Feel::cerr << "Error with rsync" << std::endl;
        }
    }

    const std::string full_source = login + "@" + server + ":" + dir;

    Table summary;
    summary.add_row( { "feelpp_sync_case" } );
    summary( 0, 0 ).format().setFontAlign( Font::Align::center );
    Table data;
    data.add_row( { "Source dir", full_source } );
    data.add_row( { "Target dir", local_dir } );
    data.add_row( { "Status", (e == 0) ? "\x1B[32mDone ✓\033[0m" : "\x1B[31mFailed ✘\033[0m" } );
    summary.add_row({data});
    cout << summary << std::endl;

    return e;
}

int main( int argc, char** argv )
{

	po::options_description scoptions( "sync_case options" );
	scoptions.add_options()
        ( "dir", po::value<std::string>()->default_value( "" ), "path in remote server to the cases export directory" )
        ( "login", po::value<std::string>()->default_value( "" ), "name of user on remote server" )
        ( "server", po::value<std::string>()->default_value( "" ), "name of the remote server" )
        ( "local-dir", po::value<std::string>()->default_value( "" ),
                "local path where files are copied (default to ${feeldir})" )
        ( "vb", po::value<int>()->default_value( 0 ), "verbose mode" )
		;

    fs::path initialCurrentPath = fs::current_path();

    Environment env( _argc=argc, _argv=argv,
                     _desc=scoptions,
                     _about=about( _name   = "sync_case" ,
                                   _author = "Feel++ Consortium",
                                   _email  = "feelpp-devel@feelpp.org" ),
                     _subdir=false );

    if (soption("dir") == "")
    {
        LOG(ERROR) << "Directory to sync must me specified" << std::endl;
        Feel::cout << tc::bold << tc::red << "Directory to sync must me specified" << tc::reset << std::endl;
        return -1;
    }


    const auto cfg = Repository::Config();
    const std::string server = (soption("server") == "") ? cfg.dist.server : soption("server");
    const std::string login = (soption("login") == "") ? cfg.dist.login : soption("login");
    const std::string local_dir = (soption("local-dir") == "") ? cfg.global_root.string() + "/" + cfg.dist.name : soption("local-dir");

    if ( server == "" || login == "" )
    {
        LOG(ERROR) << "server or login not set" << std::endl;
        Feel::cout << tc::bold << tc::red << "server or login not set" << tc::reset << std::endl;
        return -1;
    }
    
    int res = sync( soption( "dir" ), login, server, local_dir, ioption( "vb" ) );

    return res;

}

