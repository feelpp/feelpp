/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 9 Feb 2018

 Copyright (C) 2018 Feel++ Consortium

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
#ifndef FEELPP_MODELMODELS_HPP
#define FEELPP_MODELMODELS_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/json.hpp>

namespace Feel {

class FEELPP_EXPORT ModelModel
{
  public :
    ModelModel() = default;
    ModelModel( std::string const& type, std::string const& name, nl::json const& p = {} );
    ModelModel( std::string const& type, nl::json const& p = {} ) : ModelModel( type, "", p ) {}
    ModelModel( ModelModel const& ) = default;
    ModelModel( ModelModel && ) = default;
    std::string type() const { return M_type; }
    std::string const& name() const { return M_name; }
    nl::json const& setup() const { return M_setup; }
    std::map<std::string,std::set<std::string>> const& submodels() const { return M_submodels; }
    std::set<std::string> const& materials() const { return M_materials; }

  private:
    std::string M_type;
    std::string M_name;
    std::map<std::string,std::set<std::string>> M_submodels; // type -> ( name1,name2,..)
    std::set<std::string> M_materials;
    nl::json M_setup;
};

class FEELPP_EXPORT ModelModelsSameType : public std::map<std::string,ModelModel>
{
  public:
    ModelModelsSameType() = default;
    ModelModelsSameType( ModelModelsSameType const& ) = default;
    ModelModelsSameType( ModelModelsSameType && ) = default;

    void setup( std::string type, nl::json const& jarg );
};

class FEELPP_EXPORT ModelModels : public std::map<std::string,ModelModelsSameType>
{
  public :
    ModelModels() = default;
    ModelModels( ModelModels const& ) = default;
    ModelModels( ModelModels && ) = default;
    ModelModelsSameType const& models( std::string const& type ) const { CHECK( this->hasType( type ) ) << "no type registered"; return this->find( type )->second; }
    bool hasType( std::string const& type ) const { return this->find( type ) != this->end(); }

    void setPTree( nl::json const& jarg ) { this->setup(jarg); }
  private :
    void setup( nl::json const& jarg );
};

} // namespace Feel

#endif
