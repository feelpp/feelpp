// This file is part of the Feel library
//
// Author(s): Feel++ Contortium
//      Date: 2017-07-10
//
// @copyright (C) 2017 University of Strasbourg
// @copyright (C) 2012-2017 Feel++ Consortium
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef FEELPP_MONGODB_HPP
#define FEELPP_MONGODB_HPP 1

#include <sstream>
#include <feel/feelconfig.h>

#if defined(FEELPP_HAS_MONGOCXX)
#include <mongocxx/instance.hpp>

namespace Feel
{

class MongoConfig
{
public:
    std::string name = "feelpp";
    std::string host = "localhost";
    std::string port = "27017";
    std::string user = "";
    std::string password = "";
    std::string authsrc = "admin";
    std::string collection = "";

    // Build the full URI.
    const std::string operator()()
    {
        std::stringstream ss;
        ss << "mongodb://";
        // In case MongoDB auth is enabled.
        if( not user.empty() 
            and not password.empty() )
        {
            ss << user << ":" 
                << password << "@";
        }
        ss << host << ":"
            << port << "/"
            << name << "?authSource="
            << authsrc;
        return ss.str();
    }
};

//! Singleton to guaranty a MongoDB unique instance.
class MongoCxx
{
    public:
        static std::unique_ptr<mongocxx::instance>&& instance()
        {
            if( not S_mongocxx_instance )
                S_mongocxx_instance = std::make_unique<mongocxx::instance>();
            return std::move(S_mongocxx_instance);
        }

        static void reset()
        {
            S_mongocxx_instance.reset();
        }

    private:
        MongoCxx() {}
        ~MongoCxx() {}
        MongoCxx( const MongoCxx& ) = delete ;
        MongoCxx( MongoCxx&& ) = delete ;

        static std::unique_ptr<mongocxx::instance> S_mongocxx_instance;
};

} // Feel namespace

#endif // FEELPP_HAS_MONGOCXX

#endif // FEELPP_MONGODB_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
