/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2014-10-09

  Copyright (C) 2014-2016 Feel++ Consortium

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

#include <feel/feelvf/detail/ginacbuildlibrary.hpp>

#include <feel/feelcore/environment.hpp>

#include <cstdint>
#include <sstream>

namespace Feel
{
namespace vf
{
namespace detail
{
namespace
{
std::string
ginacExprDescHash( std::string const& exprDesc )
{
    std::uint64_t hash = 14695981039346656037ULL;
    for ( unsigned char c : exprDesc )
    {
        hash ^= c;
        hash *= 1099511628211ULL;
    }

    std::ostringstream ostr;
    ostr << std::hex << hash;
    return ostr.str();
}

bool
endsWith( std::string const& value, std::string const& suffix )
{
    return value.size() >= suffix.size() && value.compare( value.size() - suffix.size(), suffix.size(), suffix ) == 0;
}

std::string
ginacExprFileName( std::string const& filename, std::string const& exprDesc )
{
    if ( filename.empty() || exprDesc.empty() )
        return filename;

    std::string exprHash = ginacExprDescHash( exprDesc );
    std::string defaultStem = "ginacExpr_" + exprHash;
    fs::path filenamePath( filename );
    std::string stem = filenamePath.filename().string();
    if ( stem == defaultStem || endsWith( stem, "." + exprHash ) )
        return filename;

    return ( filenamePath.parent_path() / fs::path( stem + "." + exprHash ) ).string();
}
}

FEELPP_EXPORT
void ginacBuildLibrary( GiNaC::lst const& exprs, GiNaC::lst const& syml, std::string const& exprDesc, std::string const& filename,
                        WorldComm const& world,
                        std::shared_ptr<GiNaC::FUNCP_CUBA>& cfun )
{
    // register the GinacExprManager into Feel::Environment so that it gets the
    // GinacExprManager is cleared up when the Environment is deleted
    static bool observed = false;
    if ( !observed )
    {
        Environment::addDeleteObserver( GinacExprManagerDefaultFileNameDeleter::instance() );
        Environment::addDeleteObserver( GinacExprManagerDeleter::instance() );
        observed = true;
    }

    std::string keyExprManager = exprDesc;
    if ( exprDesc.empty() && !filename.empty() )
        keyExprManager = filename;

    bool hasLinked = ( GinacExprManager::instance().find( keyExprManager /*exprDesc*/ /*filename*/ ) != GinacExprManager::instance().end() ) ? true : false;
    if ( hasLinked )
    {
        cfun = GinacExprManager::instance().find( keyExprManager /*exprDesc*/ /*filename*/ )->second;
    }
    else
    {
        fs::path filename_p = fs::path( filename );
        fs::path filename_parent_p = filename_p.parent_path();
        std::string filenameForCompile = filename;
        if ( !filename.empty() && !fs::path(filename).is_absolute() && ( filename_parent_p != Environment::exprRepository() ) )
            filenameForCompile = (fs::path(Environment::exprRepository()) / filename_p).string();

        filenameForCompile = ginacExprFileName( filenameForCompile, exprDesc );
        std::string filenameDescExpr = filenameForCompile.empty() ? std::string() : filenameForCompile + ".desc";
        std::string filenameWithSuffix = filenameForCompile.empty() ? std::string() : filenameForCompile + ".so";
        DVLOG(2) << "filename: " << filename << std::endl;
        DVLOG(2) << "filenameForCompile: " << filenameForCompile << std::endl;
        DVLOG(2) << "filenameWithSuffix: " << filenameWithSuffix << std::endl;

        if ( !filenameForCompile.empty() )
        {
            fs::path filenameForCompileParent = fs::path( filenameForCompile ).parent_path();
            if ( !filenameForCompileParent.empty() && !fs::exists( filenameForCompileParent ) )
                fs::create_directories( filenameForCompileParent );
            if ( !filenameForCompileParent.empty() && !fs::exists( filenameForCompileParent ) )
            {
                using namespace std::string_literals;
                throw std::logic_error( "directories "s + filenameForCompileParent.string() + " not created");
            }
        }

        DVLOG( 2 ) << "GiNaC::compile_ex with filename " << filenameForCompile << "\n";
        GiNaC::compile_ex( exprs, syml, *cfun, filenameForCompile, !exprDesc.empty() );

        hasLinked = true;
        if ( !filename.empty() )
        {
            GinacExprManager::instance().operator[]( keyExprManager /*exprDesc*/ /*filename*/ ) = cfun;

            if ( world.isMasterRank() && !exprDesc.empty() && !filenameDescExpr.empty() )
            {
                std::ofstream file( filenameDescExpr, std::ios::out | std::ios::trunc );
                file << exprDesc;
                file.close();
            }
        }
    }
}

FEELPP_EXPORT std::string
ginacGetDefaultFileName( std::string const& exprDesc, std::string const& dirLibExpr )
{
    std::string res;
    std::string managerKey = exprDesc + "\n" + dirLibExpr;
    if ( GinacExprManagerDefaultFileName::instance().find( managerKey ) != GinacExprManagerDefaultFileName::instance().end() )
        res = GinacExprManagerDefaultFileName::instance().find( managerKey )->second;
    else
    {
        std::string defaultFileNameUsed = "ginacExpr_" + ginacExprDescHash( exprDesc );
        if ( dirLibExpr.empty() )
            res = Environment::exprRepository() + "/" + defaultFileNameUsed;
        else
        {
            fs::path fsdir = fs::path( dirLibExpr );
            if ( fsdir.is_absolute() )
                res = ( fsdir / defaultFileNameUsed ).string();
            else
                res = (fs::path(Environment::exprRepository()) / fsdir / defaultFileNameUsed ).string();
        }
        GinacExprManagerDefaultFileName::instance().operator[]( managerKey ) = res;
    }
    return res;
}

} // namespace detail
} // namespace vf
} // namespace Feel
