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
#include <stdexcept>
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
    git, ///! relative to git repository
    custom ///! custom location determined by a callback function
};
// map TaskState values to JSON as strings
NLOHMANN_JSON_SERIALIZE_ENUM( Location, { { Location::unknown, nullptr },
                                          { Location::global, "global" },
                                          { Location::relative, "relative" },
                                          { Location::absolute, "absolute" },
                                          { Location::git, "git" },
                                          { Location::custom, "custom" } } )
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
        {Location::custom, "custom"},
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
    if ( location_string  == "custom" ) return Location::custom;
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
        
        /**
         * @brief Construct a Config with a custom location callback
         * 
         * @param d directory path (used as fallback or within callback)
         * @param callback function to compute the repository root after options are processed
         */
        Config( fs::path d, std::function<fs::path()> callback );
        
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
        
        /// @brief Custom callback to determine repository location (used when location == Location::custom)
        std::function<fs::path()> custom_location_callback;
        
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
    Repository& configure( fs::path const& d ) { config_.directory = d; configured_ = true; return configure(); }

    /**
     * @return get the root of the repository where the results will be stored
     */
    fs::path const& root()  const
    {
        if ( !configured_ )
            throw std::logic_error("Repository::root() called before configure()");
        return root_;
    }

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
    fs::path const& geo()  const
    {
        if ( !configured_ )
            throw std::logic_error("Repository::geo() called before configure()");
        return geo_;
    }

    /**
     * @brief get the expressions directory
     * The expressions directory contains c++ and plugins associated with the expressions
     *
     * @return get the exprs directory of the repository
     */
    fs::path const& exprs()  const
    {
        if ( !configured_ )
            throw std::logic_error("Repository::exprs() called before configure()");
        return exprs_;
    }

    /**
     * @brief get the logs directory
     * The logs directory contains feelpp apps log files
     *
     * @return get the logs directory of the repository
     */
    fs::path const& logs() const
    {
        if ( !configured_ )
            throw std::logic_error("Repository::logs() called before configure()");
        return logs_;
    }

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
     * @brief check if the repository is configured
     *
     * @return true if repository is configured, false otherwise
     */
    bool isConfigured() const noexcept { return configured_; }

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
     * 
     * @return true if repository uses a custom location callback, false otherwise
     */
    bool isCustom() const { return config_.location == Location::custom; }

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

    /**
     * @brief Verify that all repository directories exist
     * 
     * This function checks that root, geo, logs, results, and exprs directories
     * have been created. It's useful for debugging to find where directory access
     * happens before creation.
     * 
     * @param caller_info Optional string describing where this check is called from
     * @return true if all directories exist, false otherwise
     */
    bool verifyDirectoriesExist(std::string const& caller_info = "") const;

private:
    bool configured_ = false;
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
/**
 * @brief Create a Repository Config with custom location callback
 * 
 * The callback will be invoked during Repository::configure() after options are processed.
 * This allows for dynamic repository location based on runtime configuration.
 * 
 * @param reldir fallback directory (can be used within callback)
 * @param callback function that returns the computed repository path
 * @return Repository::Config configured with custom location
 */
inline Repository::Config customRepository( std::string reldir, std::function<fs::path()> callback )
{
    return Repository::Config(fs::path(reldir), callback);
}

/**
 * @brief Create repository config for git location
 * 
 * Tells Feel++ to detect the nearest .git directory and use <git-root>/feelppdb
 * 
 * @param start_path starting path for .git search (default: current directory)
 * @return Repository::Config configured with git location
 */
inline Repository::Config gitRepository( fs::path start_path = "." )
{
    return Repository::Config(start_path, Location::git);
}

/**
 * @brief Create repository config for absolute path location
 * 
 * Use a fixed absolute path for the repository root.
 * Useful for temporary directories or per-run storage.
 * 
 * @param abs_path absolute path to use as repository root
 * @return Repository::Config configured with absolute location
 */
inline Repository::Config absoluteRepository( fs::path abs_path )
{
    return Repository::Config(abs_path, Location::absolute);
}
} // namespace Feel
