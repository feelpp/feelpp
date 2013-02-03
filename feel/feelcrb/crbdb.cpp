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
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcrb/crbdb.hpp>

namespace Feel
{
CRBDB::CRBDB()
    :
    M_name( "noname" ),
    M_vm(),
    M_isloaded( false )
{}

//! constructor from command line options
CRBDB::CRBDB( std::string prefixdir,
              std::string name,
              std::string dbprefix,
              po::variables_map const& vm )
    :
    M_prefixdir( prefixdir ),
    M_name( name ),
    M_vm( vm ),
    M_isloaded( false )
{
    LOG(INFO) << prefixdir << "," << name << "\n";

    this->setDBFilename( ( boost::format( "%1%.crbdb" ) % dbprefix ).str() );
    LOG(INFO) << "database name " << dbFilename() << "\n";


}

//! copy constructor
CRBDB::CRBDB( CRBDB const & o )
    :
    M_name( o.M_name ),
    M_vm( o.M_vm ),
    M_isloaded( o.M_isloaded )
{}

//! destructor
CRBDB::~CRBDB()
{}

fs::path
CRBDB::dbSystemPath() const
{
    std::vector<std::string> sysdir = boost::assign::list_of( Feel::Info::prefix() )( "/usr" )( "/usr/local" )( "/opt/local" );
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
    return fs::path();
}
fs::path
CRBDB::dbLocalPath() const
{

    int proc_number =  Environment::worldComm().globalRank();
    int nb_proc = Environment::worldComm().globalSize();
    std::string suf;
    if( M_vm.count( "crb.results-repo-name" ) )
        {
            std::string database_name = M_vm["crb.results-repo-name"].as<std::string>();
            suf = database_name + ( boost::format("_proc%1%on%2%") %proc_number %nb_proc ).str() ;
        }
    else
        {
            std::string database_name = "default_repo";
            suf = database_name + ( boost::format("_proc%1%on%2%") %proc_number %nb_proc ).str() ;
        }


    // generate the local repository db path
    std::string localpath = ( boost::format( "%1%/db/crb/%2%/%3%" )
                              % Feel::Environment::rootRepository()
                              % M_prefixdir
                              % suf ).str();
    fs::path rep_path = localpath;
    fs::create_directories( rep_path );

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
