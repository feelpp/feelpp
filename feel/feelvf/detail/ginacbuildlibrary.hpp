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

#ifndef FEELPP_DETAIL_GINACBUILDLIBRARY_HPP
#define FEELPP_DETAIL_GINACBUILDLIBRARY_HPP 1

#include <feel/feelcore/feelmacros.hpp>
#include <feel/feelcore/environment.hpp>

// #if defined(__clang__)
// #if  !defined(__APPLE__) && FEELPP_CLANG_AT_LEAST(3,9)
// #pragma clang diagnostic push
// #pragma clang diagnostic ignored "-Wundefined-var-template"
// #endif
// #endif

#include <ginac/ginac.h>
extern template GiNaC::registered_class_info GiNaC::container<std::list>::reg_info;
extern template GiNaC::registered_class_info GiNaC::container<std::vector>::reg_info;

// #if defined(__clang__)
// #if  !defined(__APPLE__) && FEELPP_CLANG_AT_LEAST(3,9)
// #pragma clang diagnostic pop
// #endif
// #endif

#include <feel/feelcore/singleton.hpp>
#include <feel/feelcore/worldcomm.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{

/**
 * @brief
 * - Check if ginac library is alreday done or not ( compare .desc file )
 * - Genereate library (compilation+link) if necessary
 * - Store in Singleton GiNaC::FUNCP_CUBA
 */
FEELPP_EXPORT void
ginacBuildLibrary( GiNaC::lst const& exprs, GiNaC::lst const& syml, std::string const& exprDesc, std::string const& filename,WorldComm const& world,
                   std::shared_ptr<GiNaC::FUNCP_CUBA> & cfun );

/**
 * @brief get a filename for ginac lib define by use a singleton counter
 */
FEELPP_EXPORT std::string
ginacGetDefaultFileName( std::string const& exprDesc, std::string const& dirLibExpr = Environment::exprRepository() );


} // namespace detail
} // namespace vf

class FEELPP_NO_EXPORT GinacExprManagerImpl :
        public std::map<std::string, std::shared_ptr<GiNaC::FUNCP_CUBA> > ,
        public boost::noncopyable
{
public:
    typedef std::shared_ptr<GiNaC::FUNCP_CUBA> value_type;
    typedef std::string key_type;
    typedef std::map<key_type,value_type> ginac_expr_manager_type;
};

typedef Feel::Singleton<GinacExprManagerImpl> GinacExprManager;

struct FEELPP_NO_EXPORT  GinacExprManagerDeleterImpl
{
    void operator()() const
        {
            VLOG(2) << "[GinacManagerDeleter] clear GinacExprManager Singleton: " << GinacExprManager::instance().size() << "\n";
            GinacExprManager::instance().clear();
            VLOG(2) << "[GinacManagerDeleter] clear GinacExprManager done\n";
        }
};
typedef Feel::Singleton<GinacExprManagerDeleterImpl> GinacExprManagerDeleter;


class FEELPP_NO_EXPORT  GinacExprManagerDefaultFileNameImpl :
        public std::map<std::string, std::string >,
        public boost::noncopyable
{
public:
    typedef std::string value_type;
    typedef std::string key_type;
    typedef std::map<key_type,value_type> ginac_expr_manager_default_filename_type;
};

typedef Feel::Singleton<GinacExprManagerDefaultFileNameImpl> GinacExprManagerDefaultFileName;

struct FEELPP_NO_EXPORT GinacExprManagerDefaultFileNameDeleterImpl
{
    void operator()() const
        {
            VLOG(2) << "[GinacManagerDeleter] clear GinacExprManager Singleton: " << GinacExprManager::instance().size() << "\n";
            GinacExprManagerDefaultFileName::instance().clear();
            VLOG(2) << "[GinacManagerDeleter] clear GinacExprManager done\n";
        }
};
typedef Feel::Singleton<GinacExprManagerDefaultFileNameDeleterImpl> GinacExprManagerDefaultFileNameDeleter;

} // namespace Feel

#endif /* FEELPP_DETAIL_GINACBUILDLIBRARY_HPP */
