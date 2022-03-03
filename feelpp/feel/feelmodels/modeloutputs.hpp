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

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel {

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
class FEELPP_EXPORT ModelOutput : public CommObject
{
  public:
    using super = CommObject;
    ModelOutput( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ModelOutput( ModelOutput const& ) = default;
    ModelOutput( ModelOutput&& ) = default;
    ModelOutput& operator=( ModelOutput const& ) = default;
    ModelOutput& operator=( ModelOutput && ) = default;
    ModelOutput( std::string name, nl::json const& jarg,
                 worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                 std::string const& directoryLibExpr = "" );

    std::string name() const { return M_name; }
    std::string type() const { return M_type; }
    std::set<std::string> markers() const { return M_markers; }
    std::vector<double> coord() const { return M_coord; }
    double radius() const { return M_radius; }
    int dim() const { return M_dim; }
    std::string getString( std::string const& key ) const;

  private:
    nl::json M_p;
    std::string M_directoryLibExpr;

    std::string M_name;
    std::string M_type;
    ModelMarkers M_markers;
    std::vector<double> M_coord;
    double M_radius;
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
class FEELPP_EXPORT ModelOutputs: public std::map<std::string,ModelOutput>, public CommObject
{
  public:
    using super = CommObject;
    using value_type = std::map<std::string,ModelOutput>::value_type;
    ModelOutputs( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ModelOutputs( nl::json const& jarg, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    virtual ~ModelOutputs() = default;
    void setPTree( nl::json const& jarg ) { M_p = jarg; setup(); }
    ModelOutput loadOutput( std::string const&, std::string const& );
    ModelOutput getOutput( nl::json const& jarg, std::string const& );

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

    nl::json M_p;
    std::string M_directoryLibExpr;

};

}
#endif
