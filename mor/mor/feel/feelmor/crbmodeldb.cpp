/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 09 Oct 2017

 Copyright (C) 2017 Feel++ Consortium

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

#include <feel/feelmor/crbmodeldb.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>


namespace Feel
{

CRBModelDB::CRBModelDB( std::string const& name, uuids::uuid const& uid, std::string const& root )
    :
    M_name( algorithm::to_lower_copy(name ) ),
    M_root( root ),
    M_uuid( uid )
{}

std::string
CRBModelDB::jsonFilename() const
{
    return CRBModelDB::jsonFilename( this->name() );
}

std::string
CRBModelDB::jsonFilename( std::string const& name )
{
    return fmt::format("{}.crb.json",name);
}

std::string
CRBModelDB::dbRepository() const
{
    fs::path dbdir = fs::path(this->rootRepository()) / "crbdb" / this->name() / uuids::to_string(M_uuid);
    return dbdir.string();
}

void
CRBModelDB::updateIdFromDBFilename( std::string const& filename )
{
    uuids::uuid myuid = CRBModelDB::idFromDBFilename( this->name(), filename );
    if ( !myuid.is_nil() )
        M_uuid = myuid;
}
void
CRBModelDB::updateIdFromDBLast( crb::last last )
{
    uuids::uuid myuid = CRBModelDB::idFromDBLast( this->name(), last, this->rootRepository() );
    if ( !myuid.is_nil() )
        M_uuid = myuid;
}
void
CRBModelDB::updateIdFromId( std::string const& uid )
{
    uuids::uuid myuid = CRBModelDB::idFromId( this->name(), uid, this->rootRepository() );
    if ( !myuid.is_nil() )
        M_uuid = myuid;
}

uuids::uuid
CRBModelDB::idFromDBFilename( std::string const& name, std::string const& filename )
{
    if ( !fs::exists( filename ) )
        return uuids::nil_uuid();

    auto json_str_wo_comments = removeComments(readFromFile(filename));
    boost::property_tree::ptree ptree;
    std::istringstream istr( json_str_wo_comments );
    boost::property_tree::read_json( istr, ptree );

    std::string nameReaded;
    if( boost::optional<std::string> d = ptree.get_optional<std::string>("name") )
        nameReaded = *d;
    else if ( auto ptreeCrbModel = ptree.get_child_optional( "crbmodel" ) )
        nameReaded = ptreeCrbModel->get<std::string>("name");

    if ( nameReaded != name )
        return uuids::nil_uuid();

    std::string uidstring = ptree.get<std::string>("uuid");
    return boost::lexical_cast<uuids::uuid>( uidstring );
}

uuids::uuid
CRBModelDB::idFromDBLast( std::string const& name, crb::last last, std::string const& root )
{
    if ( !fs::exists( root ) )
        return uuids::nil_uuid();

    fs::path crbdb = fs::path(root) / "crbdb";
    if ( !fs::exists( crbdb ) )
        return uuids::nil_uuid();

    fs::path dbbasedir = crbdb / fs::path(name) ;
    if ( !fs::exists( dbbasedir ) && !fs::is_directory(dbbasedir) )
        return uuids::nil_uuid();
    //throw std::invalid_argument(std::string("db directory ") + dbbasedir.string() + " does not exist");
    // either id provides the full directory or part of it
    // try first full path
    typedef std::vector<fs::path> vec;
    vec d;
    typedef std::multimap<std::time_t, fs::path> result_set_t;
    result_set_t result_set;

    for( auto const& dir: boost::make_iterator_range( fs::directory_iterator(dbbasedir),{} ) )
    {
        fs::path dbfilename = dir.path() / CRBModelDB::jsonFilename( name );
        if (fs::exists( dbfilename ) )
        {
            fs::path uidpath = dir.path().filename();
            if ( last == crb::last::created )
            {
                result_set.insert(result_set_t::value_type(fs::last_write_time(dir.path()), uidpath ));
            }
            else if ( last == crb::last::modified )
            {
                result_set.insert(result_set_t::value_type(fs::last_write_time(dbfilename), uidpath ));
            }
        }
    }
    if ( result_set.size() )
    {
        std::string uidstring = result_set.rbegin()->second.string();
        //std::cout << "Last " << ((last==crb::last::modified)?"modified":"created") << " db uid: " << uidstring << std::endl;
        return boost::lexical_cast<uuids::uuid>( uidstring );
    }
    //throw std::invalid_argument(std::string("Last database for ") + name() + " not found");

    return uuids::nil_uuid();
}

uuids::uuid
CRBModelDB::idFromId( std::string const& name, std::string const& uid, std::string const& root )
{
    if ( !fs::exists( root ) )
        return uuids::nil_uuid();

    fs::path crbdb = fs::path(root) / "crbdb";
    if ( !fs::exists( crbdb ) )
        return uuids::nil_uuid();

    fs::path dbbasedir = crbdb / fs::path(name) ;
    if ( !fs::exists( dbbasedir ) && !fs::is_directory(dbbasedir) )
        return uuids::nil_uuid();

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

        if ( boost::ends_with( dbdir.string(), uid ) )
        {
            fs::path dbfilename = dbdir / CRBModelDB::jsonFilename( name );
            if (!fs::exists(dbfilename))
                continue;
            return boost::lexical_cast<uuids::uuid>( dbdir.filename().string() );
            //return dbfilename;
        }
    }
    return uuids::nil_uuid();
}

} // namespace Feel
