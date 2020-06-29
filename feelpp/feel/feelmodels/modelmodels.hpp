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

#include <boost/property_tree/ptree.hpp>
#include <feel/feelcore/feel.hpp>

namespace Feel {

namespace pt = boost::property_tree;

class FEELPP_EXPORT ModelModel
{
  public :
    ModelModel() = default;
    ModelModel( std::string const& name, pt::ptree const& p, bool addVariants = true );
    ModelModel( ModelModel const& ) = default;
    ModelModel( ModelModel && ) = default;
    std::string const& name() const { return M_name; }
    pt::ptree const& ptree() const { return M_ptree; }
    std::string const& equations() const { return M_equations; }
    std::set<std::string> const& submodels() const { return M_submodels; }
    std::map<std::string,ModelModel> const& variants() const { return M_variants; }

  private:
    std::string M_name;
    pt::ptree M_ptree;
    std::string M_equations;
    std::set<std::string> M_submodels;
    std::map<std::string,ModelModel> M_variants;
};

class FEELPP_EXPORT ModelModels : public std::map<std::string,ModelModel>
{
  public :
    ModelModels();
    ModelModels( ModelModels const& ) = default;
    ModelModels( ModelModels && ) = default;
    pt::ptree const& pTree() const { return M_p; }
    bool useModelName() const { return M_useModelName; }
    ModelModel const& model( std::string const& name = "" ) const;
    bool hasModel( std::string const& name = "" ) const;

    void setPTree( pt::ptree const& _p ) { M_p = _p; setup(); }
  private :
    void setup();
  private :
    pt::ptree M_p;
    bool M_useModelName;
    ModelModel M_emptyModel;
};

} // namespace Feel

#endif
