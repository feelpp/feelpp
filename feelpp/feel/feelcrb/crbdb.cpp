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

#include <feel/feelcore/disablewarnings.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <feel/feelcore/reenablewarnings.hpp>
//#include <boost/assign/std/vector.hpp>
#include <boost/algorithm/string.hpp>
#include <feel/feelcrb/crbdb.hpp>

namespace Feel
{
CRBDB::CRBDB( std::string const& name, std::string const& ext, worldcomm_ptr_t const& worldComm )
    :
    CRBDB( name, ext, Environment::randomUUID(true), worldComm )
{}

CRBDB::CRBDB( std::string const& name, std::string const& ext, uuids::uuid const& uuid, worldcomm_ptr_t const& worldComm )
    :
    super( worldComm ),
    M_name( algorithm::to_lower_copy(name) ),
    M_ext( ext ),
    M_uuid( uuid ),
    M_isloaded( false )
{
    if ( M_ext.empty() )
        this->setDBFilename( ( boost::format( "%1%.crbdb" ) %M_name ).str() );
    else
        this->setDBFilename( ( boost::format( "%1%.%2%.crbdb" ) %M_name%M_ext ).str() );

    this->setDBDirectory( M_name, M_uuid );
}

void
CRBDB::setDBDirectory( uuids::uuid const& i )
{
    this->setDBDirectory( M_name, i );
}
void
CRBDB::setDBDirectory( std::string const& name, uuids::uuid const& i )
{
    M_uuid = i;
    std::string database_subdir = ( boost::format( "%1%/%2%" )% name % uuids::to_string(M_uuid)).str();
    M_dbDirectory = ( boost::format( "%1%/crbdb/%2%" )
                      % Feel::Environment::rootRepository()
                      % database_subdir ).str();
}
//! destructor
CRBDB::~CRBDB()
{}

uuids::uuid
CRBDB::id( fs::path p ) const
{
    if ( !fs::exists( p ) )
        throw std::invalid_argument( "Invalid filename " + p.string() );
    auto up = p.parent_path().filename();
    return boost::lexical_cast<uuids::uuid>( up.string() );
}

fs::path
CRBDB::dbLocalPath() const
{
    fs::path rep_path = M_dbDirectory;
    if ( !M_dbSubDirectory.empty() )
        rep_path /= M_dbSubDirectory;
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

    return fs::path();
}

bool
CRBDB::findDBUuid(int l, std::string const& uid_file)
{
    uuids::uuid id = uuids::nil_uuid();
    switch( l )
    {
    case 0:
        id = this->id( uid_file );
        break;
    case 1:
        id = this->idFromDBLast( crb::last::created );
        break;
    case 2:
        id = this->idFromDBLast( crb::last::modified );
        break;
    case 3:
        id = this->idFromId( uid_file );
        break;
    default:
        break;
    }
    if( !id.is_nil() )
    {
        this->setDBDirectory(id);
        return true;
    }
    else
        return false;
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

fs::path
CRBDB::db( std::string const& f )  const
{
    auto p = fs::path( f );
    if ( !fs::exists( p ) && ( p.extension().string() == ".json" ) )
        throw std::invalid_argument("Database file " + f + " not found" );
    return p;
}

uuids::uuid
CRBDB::idFromId( std::string const& id ) const
{
    auto dir = ( boost::format( "%1%/crbdb/%2%" )
                 % Feel::Environment::rootRepository()
                 % M_name ).str();
    if( !fs::exists(dir)  && !fs::is_directory(dir) )
        throw std::invalid_argument(std::string("db directory ") + dir + " does not exist");

    std::vector<fs::path> d;

    std::copy(fs::directory_iterator(dir), fs::directory_iterator(), std::back_inserter(d));
    std::sort(d.begin(), d.end());

    for( auto const& dbdir : d )
    {
        if( boost::ends_with( dbdir.string(), id) )
        {
            fs::path dbfilename = dbdir / fs::path(this->dbFilename());
            if( !exists(dbfilename) )
                continue;
            return this->id(dbfilename);
        }
    }
    throw std::invalid_argument(std::string("Database for ") + name() + " with id " + id + " not found");
}

fs::path
CRBDB::dbFromId( std::string const& id, std::string const& root )  const
{
    if ( !fs::exists( root ) )
        throw std::invalid_argument(std::string("root repository ") + root + " does not exist");
    fs::path crbdb = fs::path(root) / "crbdb";
    if ( !fs::exists( crbdb ) )
        throw std::invalid_argument(std::string("crbdb repository ") + crbdb.string() + " does not exist");
    
    fs::path dbbasedir = crbdb / fs::path(name()) ;
    if ( !fs::exists( dbbasedir ) && !fs::is_directory(dbbasedir) )
        throw std::invalid_argument(std::string("db directory ") + dbbasedir.string() + " does not exist");
    // either id provides the full directory or part of it
    // try first full path
    typedef std::vector<fs::path> vec;
    vec d;
    
    std::copy(fs::directory_iterator(dbbasedir), fs::directory_iterator(), std::back_inserter(d));
    std::sort(d.begin(), d.end());
    
    //std::cout << "dbbasedir=" << dbbasedir.string()<< " id=" << id << std::endl;
    for( auto& dbdir : d )
    {
        //std::cout << "dbdir = " << dbdir.string() << std::endl;
        
        if ( boost::ends_with( dbdir.string(), id ) )
        {
            fs::path dbfilename = dbdir / fs::path(name() + ".crb.json");
            if (!fs::exists(dbfilename))
                continue;
            return dbfilename;
        }
    }
    throw std::invalid_argument(std::string("Database for ") + name() + " with id " + id + " not found");
}

void
CRBDB::loadDBFromId( std::string const& id, crb::load l, std::string const& root )
{
    loadDB( dbFromId( id, root ).string(), l );
}

uuids::uuid
CRBDB::idFromDBLast( crb::last last ) const
{
    auto dir = ( boost::format( "%1%/crbdb/%2%" )
                 % Feel::Environment::rootRepository()
                 % M_name ).str();
    if( !fs::exists(dir)  && !fs::is_directory(dir) )
        return uuids::nil_uuid();

    std::vector<fs::path> d;
    typedef std::multimap<std::time_t, fs::path> result_set_t;
    result_set_t result_set;

    for( auto const& dd: boost::make_iterator_range( fs::directory_iterator(dir),{} ) )
    {
        fs::path dbfilename = dd.path() / fs::path(this->dbFilename());
        if (fs::exists( dbfilename ) )
        {
            if ( last == crb::last::created )
            {
                result_set.insert(result_set_t::value_type(fs::last_write_time(dd.path()), dbfilename));
            }
            else if ( last == crb::last::modified )
            {
                result_set.insert(result_set_t::value_type(fs::last_write_time(dbfilename), dbfilename));
            }
        }
    }
    if ( result_set.size() )
    {
        fs::path dbfname =  result_set.rbegin()->second;
        std::cout << "Last " << ((last==crb::last::modified)?"modified":"created") << " db: " << dbfname.string() << std::endl;
        return this->id(dbfname);
    }
    return uuids::nil_uuid();
}

fs::path
CRBDB::dbLast( crb::last last, std::string const& root ) const
{
    if ( !fs::exists( root ) )
        throw std::invalid_argument(std::string("root repository ") + root + " does not exist");
    fs::path crbdb = fs::path(root) / "crbdb";
    if ( !fs::exists( crbdb ) )
        throw std::invalid_argument(std::string("crbdb repository ") + crbdb.string() + " does not exist");

    fs::path dbbasedir = crbdb / fs::path(name()) ;
    if ( !fs::exists( dbbasedir ) && !fs::is_directory(dbbasedir) )
        throw std::invalid_argument(std::string("db directory ") + dbbasedir.string() + " does not exist");
    // either id provides the full directory or part of it
    // try first full path
    typedef std::vector<fs::path> vec;
    vec d;
    typedef std::multimap<std::time_t, fs::path> result_set_t;
    result_set_t result_set;

    for( auto const& dir: boost::make_iterator_range( fs::directory_iterator(dbbasedir),{} ) )
    {
        fs::path dbfilename = dir.path() / fs::path(name() + ".crb.json");
        if (fs::exists( dbfilename ) )
        {
            if ( last == crb::last::created )
            {
                result_set.insert(result_set_t::value_type(fs::last_write_time(dir.path()), dbfilename));
            }
            else if ( last == crb::last::modified )
            {
                result_set.insert(result_set_t::value_type(fs::last_write_time(dbfilename), dbfilename));
            }
        }
    }
    if ( result_set.size() )
    {
        fs::path dbfname =  result_set.rbegin()->second;
        std::cout << "Last " << ((last==crb::last::modified)?"modified":"created") << " db: " << dbfname.string() << std::endl;
        return dbfname;
    }
    throw std::invalid_argument(std::string("Last database for ") + name() + " not found");
}

void
CRBDB::loadDBLast( crb::last last, crb::load l, std::string const& root ) 
{
    this->loadDB( dbLast( last, root ).string(), l );
}
}
