/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013-2016 Feel++ Consortium

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
/**
   \file crbenums.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_MOR_CRBENUMS_HPP)
#define FEELPP_MOR_CRBENUMS_HPP 1

#include <fmt/core.h>
namespace Feel {

namespace crb {

enum class stage
{
    //! offline mode is triggered
    offline = 0,

    //! online mode is triggered
    online = 1
};

/**
 * @brief get the string representation of the stage
 *
 * @param s the stage
 * @return the string representation of the stage
 */
inline std::string stageToString( stage s )
{
    switch ( s )
    {
    case stage::offline:
        return "offline";
    case stage::online:
        return "online";
    }
    LOG( WARNING ) << fmt::format( "unknown stage {} return offline", s );
    return "offline";
}

/**
 * @brief get the stage from the string representation
 *
 * @param s the string representation of the stage
 * @return the stage
 */
inline stage stringToStage( std::string const& s )
{
    if ( s == "offline" )
        return stage::offline;
    else if ( s == "online" )
        return stage::online;
    else
    {
        LOG(WARNING) << fmt::format("unknown stage {} return offline", s );
        return stage::offline;
    }
}

enum class load
{
    //! load none of the data
    none=0,
    //! load reduced basis data
    rb = 1,
    //! load finite element data
    fe = 2,
    //! load all
    all = 3
};

/**
 * get the string representation of the load enum
 *
 */
inline std::string loadToString( crb::load l )
{
    switch ( l )
    {
    case crb::load::none:
        return "none";
    case crb::load::rb:
        return "rb";
    case crb::load::fe:
        return "fe";
    case crb::load::all:
        return "all";
    default:
        LOG(WARNING) << fmt::format("unknown load enum {} returning rb", static_cast<int>( l ) );
        return "rb";
    }
}

/**
 * get the load enum from the string representation
 *
 */
inline crb::load loadFromString( std::string const& s )
{
    if ( s == "none" )
        return crb::load::none;
    else if ( s == "rb" )
        return crb::load::rb;
    else if ( s == "fe" )
        return crb::load::fe;
    else if ( s == "all" )
        return crb::load::all;
    else
    {
        LOG(WARNING) << fmt::format("unknown load string {} returning rb", s );
        return crb::load::rb;
    }
}

enum class last
{
    //! last created
    created = 1,
    //! last modified
    modified = 2
};

/**
 * @brief get string representation of the last enum
 * 
 * @param l last enum
 * @return std::string representation of the last enum
 */
inline std::string lastToString( crb::last l )
{
    switch ( l )
    {
    case crb::last::created:
        return "created";
    case crb::last::modified:
        return "modified";
    default:
        LOG(WARNING) << fmt::format("unknown last enum {} returning modified", static_cast<int>( l ) );
        return "modified";
    }
}

/**
 * @brief get last enum from the string representation
 * 
 * @param s string representation of the last enum
 * @return crb::last enum
 */
inline crb::last lastFromString( std::string const& s )
{
    if ( s == "created" )
        return crb::last::created;
    else if ( s == "modified" )
        return crb::last::modified;
    else
    {
        LOG(WARNING) << fmt::format("unknown last string {} returning modified", s );
        return crb::last::modified;
    }
}

/**
 * @brief attribute of crb db to load
 *
 * an attribute is used to identify a DB to be loaded
 */
enum class attribute
{
    //! id attribute
    id = 0,
    //! name attribute
    name = 1,
    //! last created attribute
    last_created = 2,
    //! last modified attribute
    last_modified = 3,
};

/**
 * @brief get the string representation of the attribute
 *
 * @param a the attribute
 * @return the string representation of the attribute
 */
inline std::string attributeToString( attribute a )
{
    switch ( a )
    {
    case attribute::id:
        return "id";
    case attribute::name:
        return "name";
    case attribute::last_created:
        return "last_created";
    case attribute::last_modified:
        return "last_modified";
    default:
        LOG(WARNING) << fmt::format("unknown attribute enum {} returning last_modified", static_cast<int>( a ) );
        return "last_modified";
    }
    return "unknown";
}

/**
 * @brief get the attribute from the string representation
 *
 * @param s the string representation of the attribute
 * @return the attribute
 */
inline attribute attributeFromString( std::string const& s )
{
    if ( s == "id" )
        return attribute::id;
    else if ( s == "name" )
        return attribute::name;
    else if ( s == "last_created" )
        return attribute::last_created;
    else if ( s == "last_modified" )
        return attribute::last_modified;
    else
    {
        LOG(WARNING) << fmt::format("unknown attribute string {} returning last_modified", s );
        return attribute::last_modified;
    }

}

} // crb

/**
 * CRBErrorType
 * Determine the type of error estimation used
 * - CRB_RESIDUAL : use the residual error estimation without algorithm SCM
 * - CRB_RESIDUAL_SCM : use the residual error estimation and also algorithm SCM
 * - CRB_NO_RESIDUAL : in this case we don't compute error estimation
 * - CRB_EMPIRICAL : compute |S_n - S_{n-1}| where S_n is the output obtained by using a reduced basis with n elements
 */
enum CRBErrorType { CRB_RESIDUAL = 0, CRB_RESIDUAL_SCM=1, CRB_NO_RESIDUAL=2 , CRB_EMPIRICAL=3};


}

#endif /* FEELPP_CRBENUMS_HPP */
