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
    if ( fs::exists( p/".git" ) && ( fs::status(p).permissions() == fs::perms::owner_write ) )
        return p;
    else if ( fs::status(p.parent_path()).permissions() == fs::perms::owner_write )
        return findGitDirectory( p.parent_path() );
    else 
        return {};
}
Repository::Repository( Config c )
    :
    config_( c )
{
    if ( config_.directory.is_absolute() )
    {
        root_= config_.directory;
        // force Location::given 
        config_.location = Location::given;
    }
    else
    {
        if ( config_.location == Location::local )
            root_ = fs::current_path() / "feelpp";
        if ( config_.location == Location::git )
        {
            if ( auto p = findGitDirectory( fs::current_path() ); p )
                root_ = p.value();
            else
                throw std::invalid_argument( "the current directory is not inside a git repository" );
        }
        if ( config_.location == Location::global )
        {
            fs::path home = findHome();
            if ( fs::exists( home / "feel" ) && fs::is_directory( home / "feel" ) )
                root_ = home / "feel";
            else
                root_ = home / "feelpp";
        }
    }
    geo_ = root_ / "geo";
}


} // namespace Feel
