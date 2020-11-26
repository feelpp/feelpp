//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//! 
//! This file is part of the Feel library
//! 
//! Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! Date: 2020-11-10
//! 
//! Copyright (C) 2020 Feel++ Consortium
//! 
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//! 
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//! 
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//! 
//! 
//! \file environment.hpp
//! \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! \date 2020-11-10
//! 
#include <feel/feelcore/repository.hpp>
#include <feel/feelcore/feelio.hpp>
#include <boost/algorithm/string/replace.hpp>
namespace Feel
{
fs::path findHome()
{
    const char * home = getenv ("HOME");
    if (home == NULL) 
    {
      std::cerr << "error: HOME variable not set." << std::endl;
      throw std::invalid_argument ("error: HOME environment variable not set.");
    }
    return fs::path( home );
}
std::optional<fs::path> findGitDirectory( fs::path p )
{
    //Feel::cout << "path:"  << p << " parent: " << p.parent_path() << std::endl;
    if ( fs::exists( p/".git" ) && ( (fs::status(p).permissions() & fs::perms::owner_all) != fs::perms::no_perms ) )
        return p;
    else if ( !p.parent_path().empty() && (fs::status(p.parent_path()).permissions() & fs::perms::owner_all) != fs::perms::no_perms )
        return findGitDirectory( p.parent_path() );
    else 
        return {};
}
Repository::Repository( Config c )
    :
    config_( c )
{}
void
Repository::configure()
{
    fs::path home = findHome();
    if ( fs::exists( home / "feel" ) && fs::is_directory( home / "feel" ) )
        global_root_ = home / "feel";
    else
        global_root_ = home / "feelppdb";

    if ( config_.directory.is_absolute() )
    {
        root_= config_.directory;
        // force Location::given 
        config_.location = Location::given;
    }
    else
    {
        if ( config_.location == Location::local )
            root_ = fs::current_path() / "feelppdb";
        if ( config_.location == Location::git )
        {
            if ( auto p = findGitDirectory( fs::current_path() ); p )
                root_ = p.value() / "feelppdb";
            else
                throw std::invalid_argument( ("the current directory " + fs::current_path().string()  + "  is not inside a git repository").c_str() );
        }
        if ( config_.location == Location::global )
        {
            root_ = global_root_;
        }
    }
    geo_ = root_ / "geo";

    // read config file in root_/.feelppconfig then try $HOME/.feelppconfig
    if ( fs::exists( root_/".feelppconfig" ) )
    {
        std::ifstream i((root_/".feelppconfig").string());
        i >> config_.data;
    }
    else if ( config_.location != Location::global && fs::exists( global_root_.parent_path()/".feelppconfig" ) )
    {
        Feel::cout << "[feelpp] reading config " << this->globalRoot().parent_path()/".feelppconfig" << std::endl;
        std::ifstream i((this->globalRoot().parent_path()/".feelppconfig").string());
        i >> config_.data;
    }
    std::string e_d = "exprs";
    if ( Environment::vm().count( "subdir.expr" ) )
        e_d = Environment::vm()["subdir.expr"].as<std::string>();
    exprs_ = ( !config_.data.empty() ) ? directoryWithoutAppenders() / config_.data.value("/directory/exprs"_json_pointer,e_d): directoryWithoutAppenders() / e_d;

    if ( Environment::isMasterRank() )
    {
        if ( !fs::exists( root_ ) )
        {
            bool created = fs::create_directories( root_ );
            Feel::cout << "[feelpp] create Feel++ repository " << root_ << " ..." << std::endl;
        }
        if ( !fs::exists( geo_ ) )
        {
            bool created = fs::create_directories( geo_ );
            Feel::cout << "[feelpp] create Feel++ geo repository " << geo_ << " ..." << std::endl;
        }
        if ( !fs::exists( directory() ) )
        {
            bool created = fs::create_directories( directory() );
            Feel::cout << "[feelpp] create Feel++ results directory " << this->directory() << " ..." << std::endl;
        }
        if ( !fs::exists( exprs() ) )
        {
            bool created = fs::create_directories( exprs() );
            Feel::cout << "[feelpp] create Feel++ exprs directory " << ( this->exprs() ) << " ..." << std::endl;
        }
    }
    Environment::worldComm().barrier();
    
}


fs::path
Repository::directory() const 
{ 
    return  root() / relativeDirectory(); 
}
fs::path
Repository::directoryWithoutAppenders() const 
{ 
    if ( isGiven() )
        return config_.directory;
    return root() / config_.directory;
}
fs::path 
Repository::relativeDirectory() const 
{ 
    if ( isGiven() ) return {};
    else
    {
        fs::path p = config_.directory;
        if ( !config_.data.empty() )
        {
            std::cout << config_.data.dump( 1 ) << std::endl;
            bool append_date = config_.data.value( "/directory/append/date"_json_pointer, false );
            if ( append_date )
            {
                using boost::gregorian::day_clock;
                using boost::posix_time::ptime;
                using boost::posix_time::second_clock;
                using boost::posix_time::to_simple_string;

                ptime todayUtc( day_clock::universal_day(), second_clock::universal_time().time_of_day() );
                std::string today = boost::replace_all_copy( boost::replace_all_copy( to_simple_string( todayUtc ), " ", "-" ), ":", "-" );
                p = p / today;
            }
            bool append_np = config_.data.value( "/directory/append/np"_json_pointer, false );
            if ( append_np )
                p /= "np_" + std::to_string( Environment::numberOfProcessors() );
        }
        return p; 
    }        
}

} // namespace Feel
