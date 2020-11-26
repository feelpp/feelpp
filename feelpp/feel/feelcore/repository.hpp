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
enum class Location {
    unknown=0,
    global=10,
    local,
    git,
    given
};
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
        {Location::local, "local"},
        {Location::git, "git"},
        {Location::given, "given"}
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
    if ( location_string  == "local" ) return Location::local;
    if ( location_string  == "git" ) return Location::git;
    if ( location_string  == "given" ) return Location::given;
    return Location::unknown;
}
/**
 * @brief data structure to handle repository of results
 * 
 */
class Repository
{

public:
    struct Config {
        Config() = default;
        Config( fs::path d, Location l): directory(d), location(l) {};
        Config( fs::path d, Location l, nl::json dat ): directory(d), location(l), data(dat) {};
        fs::path directory;
        Location location =  Location::global;
        nl::json data;
    };

    Repository() = default;
    Repository( Config c );

    /**
     * @brief configure the repository directories based on @p Config
     * 
     */
    void configure();
    
    /**
     * @brief configure for a given directory
     * 
     * @param d directory path
     */
    void configure( fs::path const& d ) { config_.directory = d; configure(); }

    /**
     * @return get the root of the repository
     */
    fs::path const& root()  const { return root_; }

    /**
     * @return get the global root repository
     */
    fs::path const& globalRoot()  const { return global_root_; }

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
    bool isLocal() const { return config_.location == Location::local; }

    /**
     * 
     * @return true if repository is global, false otherwise
     */
    bool isGlobal() const { return config_.location == Location::global; }

    /**
     * 
     * @return true if repository is relative to a git repository, false otherwise
     */
    bool isGit() const { return config_.location == Location::git; }

    /**
     * if Config.directory is absolute, then the repository is considered given
     * @return true if repository is given, false otherwise
     */
    bool isGiven() const { return config_.location == Location::given; }

    /**
     * @brief the user name
     * 
     * @return std::string 
     */
    std::string userName() const { return ( config_.data.contains("user") && config_.data["user"].contains("name") )?config_.data["user"]["name"].get<std::string>():""; }

    /**
     * @brief the user email
     * 
     * @return std::string 
     */
    std::string userEmail() const { return ( config_.data.contains("user") && config_.data["user"].contains("email") )?config_.data["user"]["email"].get<std::string>():""; }

    /**
     * @brief data configuration
     * 
     */
    json const& data() const { return config_.data; }
private:
    Config config_;
    fs::path root_;
    fs::path global_root_;
    fs::path geo_;
    fs::path exprs_;
};

inline Repository::Config globalRepository( std::string reldir, nl::json d = {} )
{
    return Repository::Config(fs::path(reldir), Location::global, d);
}
inline Repository::Config localRepository( std::string reldir, nl::json d = {} )
{
    return Repository::Config(fs::path(reldir), Location::local, d );
}   
inline Repository::Config unknownRepository()
{
    return Repository::Config({}, Location::unknown, {} );
}
} // namespace Feel
