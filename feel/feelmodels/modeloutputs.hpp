/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Mar 2015

 Copyright (C) 2015 Feel++ Consortium

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
#ifndef FEELPP_MODELOUTPUTS_HPP
#define FEELPP_MODELOUTPUTS_HPP 1

#include <boost/property_tree/ptree.hpp>
#include <feel/feelvf/ginac.hpp>

namespace Feel {

namespace pt =  boost::property_tree;

/** \class ModelOutput
 * \brief class containing an output
 *
 * an output is a name, a type, a range and a dimension
 * "myoutput":
 * {
 *   "type": "mytype",
 *   "range": ["marker1","marker2"],
 *   "topodim": 2
 * }
 */
class FEELPP_EXPORT ModelOutput
{
  public:
    ModelOutput( WorldComm const& worldComm = Environment::worldComm() );
    ModelOutput( ModelOutput const& ) = default;
    ModelOutput( ModelOutput&& ) = default;
    ModelOutput& operator=( ModelOutput const& ) = default;
    ModelOutput& operator=( ModelOutput && ) = default;
    ModelOutput( std::string name, pt::ptree const& p,
                 WorldComm const& worldComm = Environment::worldComm(),
                 std::string const& directoryLibExpr = "" );

    std::string name() const { return M_name; }
    std::string type() const { return M_type; }
    std::set<std::string> range() const { return M_range; }
    int dim() const { return M_dim; }
    std::string getString( std::string const& key ) const;

  private:
    WorldComm const * M_worldComm;
    pt::ptree M_p;
    std::string M_directoryLibExpr;

    std::string M_name;
    std::string M_type;
    std::set<std::string> M_range;
    int M_dim;
};

std::ostream& operator<<( std::ostream& os, ModelOutput const& o );

/** \class ModelOutputs
 * \brief class containing all outputs in the model
 *
 * key: name of the output
 * value: ModelOutput
 * "Outputs":
 * {
 *   "myoutput1":{...},
 *   "myoutput2":{...}
 * }
 */
class FEELPP_EXPORT ModelOutputs: public std::map<std::string,ModelOutput>
{
  public:
    using value_type = std::map<std::string,ModelOutput>::value_type;
    ModelOutputs( WorldComm const& worldComm = Environment::worldComm() );
    ModelOutputs( pt::ptree const& p, WorldComm const& worldComm = Environment::worldComm() );
    virtual ~ModelOutputs() = default;
    void setPTree( pt::ptree const& _p ) { M_p = _p; setup(); }
    ModelOutput loadOutput( std::string const&, std::string const& );
    ModelOutput getOutput( pt::ptree const&, std::string const& );

    ModelOutput const& output( std::string const& o ) const
    {
        auto it = this->find( o );
        if( it == this->end() )
            throw std::invalid_argument( std::string("ModelOutput: Invalid output name ") + o );
        return it->second;
    }

    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

  private:
    void setup();

    WorldComm const* M_worldComm;
    pt::ptree M_p;
    std::string M_directoryLibExpr;

};

}
#endif
