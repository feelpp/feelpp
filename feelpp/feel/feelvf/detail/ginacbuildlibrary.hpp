// SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2014-2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-2.1-or-later

#ifndef FEELPP_DETAIL_GINACBUILDLIBRARY_HPP
#define FEELPP_DETAIL_GINACBUILDLIBRARY_HPP 1

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feelmacros.hpp>

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <ginac/ginac.h>
extern template GiNaC::registered_class_info GiNaC::container<std::list>::reg_info;
extern template GiNaC::registered_class_info GiNaC::container<std::vector>::reg_info;
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <feel/feelcore/singleton.hpp>
#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelvf/detail/ginacjit.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{

/**
 * Reuse or compile an expression library collectively on the supplied WorldComm.
 * The expression description identifies the cached module; the expression
 * manager retains the resulting function pointer until Environment teardown.
 */
FEELPP_EXPORT void
ginacBuildLibrary( GiNaC::lst const& exprs, GiNaC::lst const& syml, std::string const& exprDesc, std::string const& filename, WorldComm const& world,
                   std::shared_ptr<GiNaC::FUNCP_CUBA>& cfun );

/**
 * Return a deterministic filename for a GiNaC expression library.
 */
FEELPP_EXPORT std::string
ginacGetDefaultFileName( std::string const& exprDesc, std::string const& dirLibExpr = Environment::exprRepository() );

} // namespace detail
} // namespace vf

/** Registry of compiled expression callbacks owned by the Environment. */
class FEELPP_NO_EXPORT GinacExprManagerImpl : public std::map<std::string, std::shared_ptr<GiNaC::FUNCP_CUBA>>,
                                              public boost::noncopyable
{
  public:
    typedef std::shared_ptr<GiNaC::FUNCP_CUBA> value_type;
    typedef std::string key_type;
    typedef std::map<key_type, value_type> ginac_expr_manager_type;
};

typedef Feel::Singleton<GinacExprManagerImpl> GinacExprManager;

/** Clear expression callbacks before unloading their compiled modules. */
struct FEELPP_NO_EXPORT GinacExprManagerDeleterImpl
{
    /** Release all registered callbacks. */
    void operator()() const
    {
        VLOG( 2 ) << "[GinacManagerDeleter] clear GinacExprManager Singleton: " << GinacExprManager::instance().size() << "\n";
        GinacExprManager::instance().clear();
        VLOG( 2 ) << "[GinacManagerDeleter] clear GinacExprManager done\n";
    }
};

typedef Feel::Singleton<GinacExprManagerDeleterImpl> GinacExprManagerDeleter;

/** Registry mapping expression descriptions to their persistent filenames. */
class FEELPP_NO_EXPORT GinacExprManagerDefaultFileNameImpl : public std::map<std::string, std::string>,
                                                             public boost::noncopyable
{
  public:
    typedef std::string value_type;
    typedef std::string key_type;
    typedef std::map<key_type, value_type> ginac_expr_manager_default_filename_type;
};

typedef Feel::Singleton<GinacExprManagerDefaultFileNameImpl> GinacExprManagerDefaultFileName;

/** Clear cached filename associations during Environment teardown. */
struct FEELPP_NO_EXPORT GinacExprManagerDefaultFileNameDeleterImpl
{
    /** Release all registered filename associations. */
    void operator()() const
    {
        VLOG( 2 ) << "[GinacManagerDeleter] clear GinacExprManager Singleton: " << GinacExprManager::instance().size() << "\n";
        GinacExprManagerDefaultFileName::instance().clear();
        VLOG( 2 ) << "[GinacManagerDeleter] clear GinacExprManager done\n";
    }
};

typedef Feel::Singleton<GinacExprManagerDefaultFileNameDeleterImpl> GinacExprManagerDefaultFileNameDeleter;

} // namespace Feel

#endif /* FEELPP_DETAIL_GINACBUILDLIBRARY_HPP */
