/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2014-10-09

  Copyright (C) 2014 Feel++ Consortium

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
void ginacBuildLibrary( GiNaC::lst const& exprs, GiNaC::lst const& syml, std::string const& exprDesc, std::string const& filename,WorldComm const& world,
                        boost::shared_ptr<GiNaC::FUNCP_CUBA> & cfun )
{
    // register the GinacExprManager into Feel::Environment so that it gets the
    // GinacExprManager is cleared up when the Environment is deleted
    static bool observed=false;
    if ( !observed )
    {
        Environment::addDeleteObserver( GinacExprManagerDefaultFileNameDeleter::instance() );
        Environment::addDeleteObserver( GinacExprManagerDeleter::instance() );
        observed = true;
    }

    std::string keyExprManager = exprDesc;
    if ( exprDesc.empty() && !filename.empty() )
        keyExprManager = filename;

    bool hasLinked = ( GinacExprManager::instance().find( keyExprManager/*exprDesc*//*filename*/ ) != GinacExprManager::instance().end() )? true : false;
    if ( hasLinked )
    {
        cfun = GinacExprManager::instance().find( keyExprManager/*exprDesc*//*filename*/ )->second;
    }
    else
    {
        std::string filenameDescExpr = filename + ".desc";
        std::string filenameWithSuffix = filename + ".so";

        bool doRebuildGinacLib = true;
        // load .desc file and compare if both expr are identical
        if ( world.isMasterRank() && !filename.empty() && fs::exists( filenameDescExpr ) && fs::exists( filenameWithSuffix ) )
        {
            std::string exprInFile;
            std::ifstream file( filenameDescExpr, std::ios::in);
            std::getline(file, exprInFile );
            file.close();
            // if equal else not rebuild ginac
            if ( exprInFile == exprDesc )
                doRebuildGinacLib = false;
        }

        // master rank check if the lib exist and compile this one if not done
        if ( ( world.isMasterRank() && doRebuildGinacLib ) || filename.empty() )
        {
            if ( !filename.empty() && fs::path(filename).is_absolute() && !fs::exists(fs::path(filename).parent_path()) )
                fs::create_directories( fs::path(filename).parent_path() );
            DVLOG(2) << "GiNaC::compile_ex with filenameWithSuffix " << filenameWithSuffix << "\n";
            GiNaC::compile_ex(exprs, syml, *cfun, filename);

            hasLinked=true;
            if ( !filename.empty() )
            {
                GinacExprManager::instance().operator[]( keyExprManager/*exprDesc*//*filename*/ ) = cfun;

                if ( !exprDesc.empty() )
                {
                    std::ofstream file( filenameDescExpr, std::ios::out | std::ios::trunc);
                    file << exprDesc;
                    file.close();
                }
            }
        }
        // wait the lib compilation
        if ( !filename.empty() )
            world.barrier();
        // link with other process
        if ( !hasLinked )
        {
            DVLOG(2) << "GiNaC::link_ex with filenameWithSuffix " << filenameWithSuffix << "\n";
            GiNaC::link_ex(filenameWithSuffix, *cfun);
            if ( !filename.empty() )
                GinacExprManager::instance().operator[]( keyExprManager/*exprDesc*//*filename*/ ) = cfun;
        }
    }

}

std::string
ginacGetDefaultFileName( std::string const& exprDesc )
{
    std::string res;
    if ( GinacExprManagerDefaultFileName::instance().find( exprDesc ) != GinacExprManagerDefaultFileName::instance().end() )
        res = GinacExprManagerDefaultFileName::instance().find( exprDesc )->second;
    else
    {
        std::string defaultFileNameUsed = (boost::format("ginacExprDefaultFileName%1%")%GinacExprManagerDefaultFileName::instance().size()).str();
        res = (fs::current_path()/defaultFileNameUsed).string();
        GinacExprManagerDefaultFileName::instance().operator[]( exprDesc ) = res;
    }
    return res;
}




} // namespace detail
} // namespace vf
} // namespace Feel

