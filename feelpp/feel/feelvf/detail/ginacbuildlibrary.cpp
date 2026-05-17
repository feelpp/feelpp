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

namespace Feel
{
namespace vf
{
namespace detail
{
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
        std::string filenameDescExpr = filename + ".desc";
        std::string filenameWithSuffix =  filename + ".so";
        if ( !filename.empty() && !fs::path(filename).is_absolute() && ( filename_parent_p != Environment::exprRepository() ) )
        {
            filenameForCompile = (fs::path(Environment::exprRepository()) / filename_p).string();
            filenameDescExpr = (fs::path(Environment::exprRepository()) / fs::path( filename + ".desc" )).string();
            filenameWithSuffix =  (fs::path(Environment::exprRepository()) / fs::path(filename + ".so")).string();
        }
        DVLOG(2) << "filename: " << filename << std::endl;
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
        GiNaC::compile_ex( exprs, syml, *cfun, filenameForCompile );

        hasLinked = true;
        if ( !filename.empty() )
        {
            GinacExprManager::instance().operator[]( keyExprManager /*exprDesc*/ /*filename*/ ) = cfun;

            if ( world.isMasterRank() && !exprDesc.empty() )
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
    if ( GinacExprManagerDefaultFileName::instance().find( exprDesc ) != GinacExprManagerDefaultFileName::instance().end() )
        res = GinacExprManagerDefaultFileName::instance().find( exprDesc )->second;
    else
    {
        std::string defaultFileNameUsed = ( boost::format( "ginacExprDefaultFileName%1%" ) % GinacExprManagerDefaultFileName::instance().size() ).str();
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
        GinacExprManagerDefaultFileName::instance().operator[]( exprDesc ) = res;
    }
    return res;
}

} // namespace detail
} // namespace vf
} // namespace Feel
