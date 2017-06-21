/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-06-15

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file crbdb.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-06-15
 */
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
//#include <boost/assign/std/vector.hpp>
#include <boost/algorithm/string.hpp>
#include <feel/feelcrb/crbdb.hpp>

namespace Feel
{
CRBDB::CRBDB( std::string const& name, std::string const& ext, WorldComm const& worldComm )
    :
    CRBDB( name, ext, Environment::randomUUID(true), worldComm )
{}

CRBDB::CRBDB( std::string const& name, std::string const& ext, uuids::uuid const& uuid, WorldComm const& worldComm )
    :
    M_worldComm( worldComm ),
    M_name( algorithm::to_lower_copy(name) ),
    M_ext( ext ),
    M_uuid( uuid ),
    M_isloaded( false )
{
    if ( M_ext.empty() )
        this->setDBFilename( ( boost::format( "%1%.crbdb" ) %M_name%M_ext ).str() );
    else
        this->setDBFilename( ( boost::format( "%1%.%2%.crbdb" ) %M_name%M_ext ).str() );

    this->setDBDirectory( M_uuid );
}

void
CRBDB::setDBDirectory( uuids::uuid const& i )
{
    M_uuid = i;
    std::string database_subdir = ( boost::format( "%1%/%2%" )% M_name % uuids::to_string(M_uuid)).str();
    M_dbDirectory = ( boost::format( "%1%/crbdb/%2%" )
                      % Feel::Environment::rootRepository()
                      % database_subdir ).str();
}
//! destructor
CRBDB::~CRBDB()
{}

fs::path
CRBDB::dbSystemPath() const
{
#if 0
    std::vector<std::string> sysdir{Feel::Info::prefix(), "/usr", "/usr/local", "/opt/local"};
    BOOST_FOREACH( auto dir, sysdir )
    {
        // generate the local repository db path
        std::string syspath = ( boost::format( "%1%/share/feel/db/crb/%2%/" )
                                % dir
                                % M_prefixdir ).str();

        //std::cout << "Look for " << syspath  << "\n";
        if ( fs::exists( syspath ) )
            return syspath;
    }
#endif
    return fs::path();
}
fs::path
CRBDB::dbLocalPath() const
{
#if 0
    int proc_number =  Environment::worldComm().globalRank();
    int nb_proc = Environment::worldComm().globalSize();
    std::string suf;
    if( Environment::vm().count( "crb.results-repo-name" ) )
        {
            std::string database_name = soption(_name="crb.results-repo-name");
            suf = database_name + ( boost::format("_proc%1%on%2%") %proc_number %nb_proc ).str() ;
        }
    else
    {
        std::string database_name = "default_repo";
        suf = M_name + "_" + database_name + ( boost::format("_proc%1%on%2%") %proc_number %nb_proc ).str() ;
    }

    // generate the local repository db path
    std::string localpath = ( boost::format( "%1%/db/crb/%2%/%3%" )
                              % Feel::Environment::rootRepository()
                              % M_prefixdir
                              % suf ).str();
    fs::path rep_path = localpath;
    fs::create_directories( rep_path );
#endif
    fs::path rep_path = M_dbDirectory;
    return rep_path;
}

fs::path
CRBDB::lookForDB() const
{
    //std::cout << "db fdilename=" << this->dbFilename() << "\n";
    // look in local repository $HOME/feel/db/crb/...
    if ( fs::exists( this->dbLocalPath() / this->dbFilename() ) )
    {
        //std::cout << "[CRBDB::lookForDB] found database in " << this->dbLocalPath() << "\n";
        return this->dbLocalPath() / this->dbFilename();
    }

    // then look into the system for install databases
    if ( fs::exists( this->dbSystemPath() / this->dbFilename() ) )
    {
        //std::cout << "[CRBDB::lookForDB] found database in " << this->dbSystemPath() << "\n";
        return this->dbSystemPath() / this->dbFilename();
    }

    return fs::path();
}
void
CRBDB::saveDB()
{
}
bool
CRBDB::loadDB()
{
    return false;
}

}
