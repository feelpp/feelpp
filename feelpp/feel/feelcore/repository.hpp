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
#pragma once

#include <string>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/json.hpp>

namespace Feel
{

std::string findUser();
fs::path findHome();

enum class Location {
    unknown=0, ///! unknown directory
    global=10, ///! global repository
    relative, ///! relative to current directory
    absolute, ///! absolute directory given 
    git ///! relative to git repository
};
// map TaskState values to JSON as strings
NLOHMANN_JSON_SERIALIZE_ENUM( Location, { { Location::unknown, nullptr },
                                          { Location::global, "global" },
                                          { Location::relative, "relative" },
                                          { Location::absolute, "absolute" },
                                          { Location::git, "git" } } )
/**
 * @brief get the location strings assocation to enum Location
 * 
 * @return const std::map<location, std::string>& 
 */
inline const std::map<Location, std::string> &location_strings() 
{
    static const std::map<Location, std::string> location_strings = {
        {Location::unknown, ""},
        {Location::global, "global"},
        {Location::relative, "relative"},
        {Location::absolute, "absolute"},
        {Location::git, "git"},
    };
    return location_strings;
}
/**
 * @brief get location string from Location enum
 * the string is empty if the Location is unknown
 * @param location_enum 
 * @return const std::string& 
 */
inline const std::string &location(Location location_enum) noexcept 
{
    auto pos = location_strings().find(location_enum);
    if(pos == location_strings().end()) 
    {
      static std::string empty_string;
      return empty_string;
    }
    return pos->second;
}
/**
 * @brief get the Location enum from string
 * 
 * @param location_string string encoding the Location
 * @return Location 
 */
inline Location location(const std::string &location_string) noexcept 
{
    if(location_string.size() < 3)
      return Location::unknown;
    if ( location_string  == "global" ) return Location::global;
    if ( location_string  == "relative" ) return Location::relative;
    if ( location_string  == "absolute" ) return Location::absolute;
    if ( location_string  == "git" ) return Location::git;
    return Location::unknown;
};
struct GithubUser
{
    GithubUser();
    bool update();
    std::optional<std::string> login;
    std::optional<std::string> name;
    std::optional<std::string> email;
    std::optional<std::string> company;
    std::optional<std::string> location;
    std::optional<std::string> html_url;
    std::optional<std::string> blog;
    NLOHMANN_DEFINE_TYPE_INTRUSIVE( GithubUser, login, name, email, company, location, html_url, blog )
};
/**
 * @ structure to handle repository of results
 * 
 */
class Repository
{

public:
    struct Owner {
        std::string name = {};
        std::string email = {};
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(Owner,name,email)
    };

    struct Dist {
        std::string login = {};
        std::string server = {};
        std::string name = {};
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(Dist,login,server,name)
    };


    struct Config {
        Config();
        Config( nl::json const& j );
        Config( nl::json && j );
        Config( fs::path d, Location l);
        Config( fs::path d, Location l, nl::json const& dat );
        Owner owner;
        Dist dist;
        fs::path feelppdb = "feelppdb";
        Location location = Location::global;
        fs::path global_root;
        fs::path directory;
        fs::path exprs = "exprs";
        fs::path logs = "logs";
        fs::path geos = "geo";
        bool append_date = false;
        bool append_np = true;
        NLOHMANN_DEFINE_TYPE_INTRUSIVE( Config, owner, dist, feelppdb, location, global_root, directory, exprs, logs, geos, append_date, append_np )
    };

    Repository() = default;
    Repository( Config c );

    /**
     * @brief configure the repository directories based on @p Config
     * 
     */
    Repository& configure();

    /**
     * @brief configure for a json configuration
     *
     * @param j json file to configure the repo
     */
    Repository& configure( nl::json const& j )
    {
        config_ = j.get<Config>();
        return configure();
    }
    /**
     * @brief configure for a given directory
     *
     * @param d directory path
     * @param l type of location 
     */
    Repository& configure( fs::path const& d, Location const& l )
    {
        config_.directory = d;
        config_.location = l;
        return configure();
    }

    /**
     * @brief configure for a given directory
     * 
     * @param d directory path
     */
    Repository& configure( fs::path const& d ) { config_.directory = d; return configure(); }

    /**
     * @return get the root of the repository where the results will be stored
     */
    fs::path const& root()  const { return root_; }

    /**
     * @return get the global root repository associated to the global location
     */
    fs::path const& globalRoot()  const { return config_.global_root; }

    /**
     * @brief get the geo directory
     * The geo directory contains  geo files of Gmsh
     *
     * @return get the geo directory of the repository
     */
    fs::path const& geo()  const { return geo_; }

    /**
     * @brief get the expressions directory
     * The expressions directory contains c++ and plugins associated with the expressions
     *
     * @return get the exprs directory of the repository
     */
    fs::path const& exprs()  const { return exprs_; }

    /**
     * @brief get the logs directory
     * The logs directory contains feelpp apps log files
     *
     * @return get the logs directory of the repository
     */
    fs::path const& logs() const { return logs_; }

    /**
     * @brief the result repository absolute directory 
     * 
     * @return fs::path 
     */
    fs::path directory() const; 

    /**
     * @brief the result repository absolute directory without appenders
     * 
     * @return fs::path 
     */
    fs::path directoryWithoutAppenders() const;

    /**
     * @brief the result repository directory relative to root
     * 
     * @return fs::path 
     */
    fs::path relativeDirectory() const;

    /**
     *
     * @return true if repository is local, false otherwise
     */
    bool isLocal() const { return config_.location == Location::relative; }

    /**
     *
     * @return true if repository is local, false otherwise
     */
    bool isRelative() const { return config_.location == Location::relative; }

    /**
     * 
     * @return true if repository is global, false otherwise
     */
    bool isGlobal() const { return config_.location == Location::global; }

    /**
     *
     * @return true if repository is absolute, false otherwise
     */
    bool isAbsolute() const { return config_.location == Location::absolute; }

    /**
     * 
     * @return true if repository is relative to a git repository, false otherwise
     */
    bool isGit() const { return config_.location == Location::git; }

    /**
     * @brief the user name
     * 
     * @return std::string 
     */
    std::string userName() const { return config_.owner.name; }

    /**
     * @brief the user email
     * 
     * @return std::string 
     */
    std::string userEmail() const { return config_.owner.email; }

    /**
     * @brief get the configuration of the repository
     * 
     * @return Config const& 
     */
    Config const& config() const { return config_; }

    /**
     * @brief get the configuration of the repository
     * 
     * @return Config & 
     */
    Config& config() { return config_; }

    /**
     * @brief change directory to current configuration
     * 
     * @return Repository& return the current directory
     */
    Repository& cd();

private:
    Config config_;
    fs::path root_;
    fs::path geo_;
    fs::path exprs_;
    fs::path logs_;
};

inline Repository::Config globalRepository( std::string reldir, nl::json d = {} )
{
    return Repository::Config(fs::path(reldir), Location::global, d);
}
inline Repository::Config localRepository( std::string reldir, nl::json d = {} )
{
    return Repository::Config(fs::path(reldir), Location::relative, d );
}   
inline Repository::Config unknownRepository()
{
    return Repository::Config({}, Location::unknown, {} );
}
} // namespace Feel
