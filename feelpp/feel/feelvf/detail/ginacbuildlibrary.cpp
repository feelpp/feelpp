// SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2014-2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-2.1-or-later

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
/** Hash an expression description deterministically for persistent cache keys. */
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

/** Test whether a cache stem already carries the requested expression hash. */
bool endsWith( std::string const& value, std::string const& suffix )
{
    return value.size() >= suffix.size() && value.compare( value.size() - suffix.size(), suffix.size(), suffix ) == 0;
}

/** Add the expression hash to explicit filenames that are not yet addressed by content. */
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
} // namespace

FEELPP_EXPORT
void ginacBuildLibrary( GiNaC::lst const& exprs, GiNaC::lst const& syml, std::string const& exprDesc, std::string const& filename,
                        WorldComm const& world,
                        std::shared_ptr<GiNaC::FUNCP_CUBA>& cfun )
{
    // Clear callbacks and filename associations before Environment unloads the modules.
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
        if ( !filename.empty() && !fs::path( filename ).is_absolute() && ( filename_parent_p != Environment::exprRepository() ) )
            filenameForCompile = ( fs::path( Environment::exprRepository() ) / filename_p ).string();

        filenameForCompile = ginacExprFileName( filenameForCompile, exprDesc );
        std::string filenameDescExpr = filenameForCompile.empty() ? std::string() : filenameForCompile + ".desc";
        std::string filenameWithSuffix = filenameForCompile.empty() ? std::string() : filenameForCompile + ".so";
        DVLOG( 2 ) << "filename: " << filename << std::endl;
        DVLOG( 2 ) << "filenameForCompile: " << filenameForCompile << std::endl;
        DVLOG( 2 ) << "filenameWithSuffix: " << filenameWithSuffix << std::endl;

        DVLOG( 2 ) << "Feel++ expression compilation with filename " << filenameForCompile << "\n";
        Feel::detail::ginacjit::compile_ex( exprs, syml, *cfun, filenameForCompile, !exprDesc.empty(), world.localComm() );

        hasLinked = true;
        if ( !filename.empty() )
        {
            GinacExprManager::instance().operator[]( keyExprManager /*exprDesc*/ /*filename*/ ) = cfun;

            if ( world.localRank() == 0 && !exprDesc.empty() && !filenameDescExpr.empty() )
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
                res = ( fs::path( Environment::exprRepository() ) / fsdir / defaultFileNameUsed ).string();
        }
        GinacExprManagerDefaultFileName::instance().operator[]( managerKey ) = res;
    }
    return res;
}

} // namespace detail
} // namespace vf
} // namespace Feel
