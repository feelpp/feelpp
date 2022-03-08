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
#include <boost/algorithm/string/replace.hpp>
#include <cpr/cpr.h>
#include <feel/feelcore/feelio.hpp>
#include <feel/feelcore/repository.hpp>
#include <fmt/core.h>
#include <pwd.h>

namespace Feel
{
std::string findUser()
{
    auto pw = getpwuid(geteuid());
    if (pw)
    {
        return std::string { pw->pw_name };
    }
    return std::string{};
}
fs::path findHome()
{
    const char * home = getenv ("HOME");
    if (home == NULL) 
    {
        home = getpwuid( getuid() )->pw_dir;
    }
    return fs::path( home );
}
std::optional<fs::path> findGitDirectory( fs::path p )
{
    //Feel::cout << "path:"  << p << " parent: " << p.parent_path() << std::endl;
    if ( fs::exists( p/".git" ) && ( (fs::status(p).permissions() & fs::perms::owner_all) != fs::perms::none ) )
        return p;
    else if ( !p.parent_path().empty() && (fs::status(p.parent_path()).permissions() & fs::perms::owner_all) != fs::perms::none )
        return findGitDirectory( p.parent_path() );
    else 
        return {};
}
GithubUser::GithubUser()
{
}

bool
GithubUser::update()
{
    fmt::print("githubuser update\n");
    if (login && !login->empty())
    {
        cpr::Response r = cpr::Get( cpr::Url{ fmt::format("https://api.github.com/users/",*login) } );
        if ( r.status_code == 200 )
        {
            nl::json j = nl::json::parse( r.text );
            from_json(j, *this);
            //std::cout << fmt::format( "github : {}\n\n", j.dump() );
            return true;
        }
        else
        {
            LOG( WARNING ) << fmt::format( "invalid request to get the github user information: {}, status: {}", *login, r.status_code );
        }
    }
    return false;
}
Repository::Config::Config()
{
    nl::json jc;
    to_json( jc, *this );
    auto home = findHome();
    VLOG(2) << fmt::format( "home: {}\n", home.string() ) << std::endl;
    if ( fs::exists( home / ".feelppconfig" ) )
    {
        VLOG(2) << fmt::format("loading {}\n", ( home / ".feelppconfig" ).string() );
        std::ifstream i( ( home / ".feelppconfig" ).string() );
        nl::json j;
        i >> j;
        jc.merge_patch( j );
        from_json( jc, *this );
    }
    else 
    {
        if ( fs::exists( home / "feel" ) && fs::is_directory( home / "feel" ) )
            global_root = home / "feel";
        else
            global_root = home / feelppdb;
        owner.name = findUser();
        to_json( jc, *this );
        std::ofstream i( ( home / ".feelppconfig" ).string() );
        i << jc.dump( 1 );
    }
    VLOG(2) << fmt::format( "config: {}", jc.dump(1) ) << std::endl;
#if 0    
    if ( owner.github )
    {
        owner.github->update();
        if ( fs::exists( home / ".feelppconfig" ) )
        {
            std::ofstream i( ( home / ".feelppconfig" ).string() );
            nl::json j;
            to_json( j, *this );
            i << j.dump(1);
        }
    }
#endif    
}
Repository::Config::Config( nl::json const& j ): Config()
{
    nl::json jc;
    to_json( jc, *this );
    jc.merge_patch( j );
    to_json( jc, *this );
}

Repository::Config::Config( nl::json&& j )
    : Config()
{
    nl::json jc;
    to_json( jc, *this );
    jc.merge_patch( std::move(j) );
    to_json( jc, *this );
}
Repository::Config::Config( fs::path d, Location l ): Config()
{
    directory = d;
    location = l;
}
Repository::Config::Config( fs::path d, Location l, nl::json const& dat )
    : Config( dat )
{
    directory = d;
    location = l;
};

Repository::Repository( Config c )
    :
    config_( c )
{}
Repository&
Repository::configure()
{
    if ( config_.directory.is_absolute() )
    {
        root_= config_.directory;
        // force Location::given 
        config_.location = Location::absolute;
    }
    else
    {
        if ( config_.location == Location::relative )
            root_ = fs::current_path() / config_.feelppdb;
        if ( config_.location == Location::git )
        {
            if ( auto p = findGitDirectory( fs::current_path() ); p )
                root_ = p.value() / config_.feelppdb;
            else
                throw std::invalid_argument( ("the current directory " + fs::current_path().string()  + "  is not inside a git repository").c_str() );
        }
        if ( config_.location == Location::global )
        {
            root_ = config_.global_root;
        }
    }
    geo_ = root_ / config_.geos;
    exprs_ = directoryWithoutAppenders() / config_.exprs;
    logs_ = directory() / config_.logs;

    if ( Environment::isMasterRank() )
    {
        auto create_dir = []( fs::path const& d, std::string const& str )
        {
            try
            {
                bool created = fs::create_directories( d );
                Feel::cout << fmt::format("[feelpp] create Feel++ {}: {}\n", str, d.string()) << std::endl;
            }
            catch ( const boost::filesystem::filesystem_error& e )
            {
                std::cerr << fmt::format("[feelpp.boost.filesystem.filesystem_error] cannot create directory {}: {}", str, d.string()) << std::endl;
                throw;
            }
        };
        if ( !fs::exists( root_ ) )
        {
            create_dir(root_, "root repository");
        }
        if ( !fs::exists( geo_ ) )
        {
            create_dir( geo_, "geo repository");
        }
        if ( !fs::exists( directory() ) )
        {
            create_dir( directory(), "results directory" );
        }
        if ( !fs::exists( exprs() ) )
        {
            create_dir( exprs(), "expressions directory" );
        }
        if ( !fs::exists( logs() ) )
        {
            create_dir( logs(), "logs directory" );
        }
    }
    Environment::worldComm().barrier();
    return *this;
}

Repository&
Repository::cd()
{
    ::chdir( directory().string().c_str() );
    return *this;
}

fs::path
Repository::directory() const 
{ 
    return  root() / relativeDirectory(); 
}
fs::path
Repository::directoryWithoutAppenders() const 
{ 
    if ( isAbsolute() )
        return config_.directory;
    return root() / config_.directory;
}
fs::path 
Repository::relativeDirectory() const 
{ 
    fs::path p;
    if ( !isAbsolute() )
        p = config_.directory;
    if ( config_.append_date )
    {
        using boost::gregorian::day_clock;
        using boost::posix_time::ptime;
        using boost::posix_time::second_clock;
        using boost::posix_time::to_simple_string;
        ptime todayUtc( day_clock::universal_day(), second_clock::universal_time().time_of_day() );
        std::string today = boost::replace_all_copy( boost::replace_all_copy( to_simple_string( todayUtc ), " ", "-" ), ":", "-" );
        p = p / today;
    }
    if ( config_.append_np )
    {
        p /= fmt::format("np_{}", std::to_string( Environment::numberOfProcessors() ) );
    }
    return p;        
}

} // namespace Feel
