/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 09 Oct 2017

 Copyright (C) 2017 Feel++ Consortium

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

#ifndef FEELPP_CRB_CRBMODELDB_H
#define FEELPP_CRB_CRBMODELDB_H 1

//#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcrb/crbenums.hpp>

namespace Feel
{

class CRBModelDB
{
public :
    CRBModelDB( std::string const& name, std::string const& root = Environment::rootRepository() )
        :
        CRBModelDB( name, Environment::randomUUID( true ), root )
        {}
    CRBModelDB( std::string const& name, uuids::uuid const& uid, std::string const& root = Environment::rootRepository() );

    std::string const& name() const { return M_name; }
    std::string const& rootRepository() const { return M_root; }
    uuids::uuid const& uuid() const { return M_uuid; }

    std::string jsonFilename() const;
    std::string dbRepository() const;

    void setName( std::string const& name ) { M_name = algorithm::to_lower_copy(name); }
    void setId( uuids::uuid const& i ) { M_uuid = i; }

    void updateIdFromDBFilename( std::string const& filename );
    void updateIdFromDBLast( crb::last last );
    void updateIdFromId( std::string const& uid );

    static uuids::uuid idFromDBFilename( std::string const& name, std::string const& filename );
    static uuids::uuid idFromDBLast( std::string const& name, crb::last last, std::string const& root = Environment::rootRepository() );
    static uuids::uuid idFromId( std::string const& name, std::string const& uid, std::string const& root = Environment::rootRepository() );

private :
    std::string M_name;
    std::string M_root;
    uuids::uuid M_uuid;
};

} // namespace Feel

#endif
