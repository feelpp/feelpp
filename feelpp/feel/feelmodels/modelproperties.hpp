/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 15 Mar 2015

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
#ifndef FEELPP_MODELPROPERTIES_HPP
#define FEELPP_MODELPROPERTIES_HPP 1

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feelmodels/modelmodels.hpp>
#include <feel/feelmodels/modelparameters.hpp>
#include <feel/feelmodels/modelmaterials.hpp>
#include <feel/feelmodels/modelpostprocess.hpp>
#include <feel/feelmodels/modeloutputs.hpp>
#include <feel/feelmodels/modelboundaryconditions.hpp>
#include <feel/feelmodels/modelinitialconditions.hpp>


namespace Feel {

namespace pt =  boost::property_tree;

class FEELPP_EXPORT ModelProperties : public CommObject
{
public:
    using super = CommObject;
    ModelProperties( std::string const& directoryLibExpr = "",
                     worldcomm_ptr_t const& world = Environment::worldCommPtr(),
                     std::string const& prefix="",
                     po::variables_map const& vm = Environment::vm() );
    ModelProperties( std::string const& filename,// = Environment::expand(soption("mod-file")),
                     std::string const& directoryLibExpr = "",
                     worldcomm_ptr_t const& world = Environment::worldCommPtr(),
                     std::string const& prefix="",
                     po::variables_map const& vm = Environment::vm() );

    template <typename T,std::enable_if_t< std::is_same_v<T,nl::json>, bool> = true >
    ModelProperties( T && jarg,
                     std::string const& directoryLibExpr = "",
                     worldcomm_ptr_t const& world = Environment::worldCommPtr(),
                     std::string const& prefix="",
                     po::variables_map const& vm = Environment::vm() )
        :
        ModelProperties( directoryLibExpr,world,prefix,vm )
    {
        M_jsonData = std::forward<T>( jarg );
        this->setupImpl();
    }

    virtual ~ModelProperties() {}

    //! setup from a collection json filename
    void setup( std::vector<std::string> const& filename );
    //! setup from json filename
    void setup( std::string const& filename ) { this->setup( std::vector<std::string>{ filename } ); }
    //! setup from json object
    template <typename T,std::enable_if_t< std::is_same_v<T,nl::json>, bool> = true >
    void setup( T && jarg )
    {
        M_jsonData = std::forward<T>( jarg );
        this->setupImpl();
    }
    //! setup from filename option
    void setupFromFilenameOption( po::variables_map const& vm = Environment::vm() );


    nl::json const& jsonData() const { return M_jsonData; }

    std::string const& name() const {  return M_name; }
    void setName( std::string const& t) { M_name = t; }
    std::string const& shortName() const {  return M_shortname; }
    void setShortName( std::string const& t) { M_shortname = t; }

    std::string const& description() const {  return M_description; }
    void setDescription( std::string const& t) { M_description = t; }

    ModelModels & models() { return M_models; }
    ModelModels const& models() const { return M_models; }

    std::string const& unit() const {  return M_unit; }
    void setUnit( std::string const& t) { M_unit = t; }

    ModelParameters & parameters()  {  return M_params; }
    ModelParameters const& parameters() const {  return M_params; }

    ModelMaterials & materials() {  return M_mat; }
    ModelMaterials const& materials() const {  return M_mat; }

    /**
     * enable BoundaryConditions2 class as a simplified BC class
     * @return ModelProperties reference for chaining
     */
    ModelProperties& enableBoundaryConditions2() { M_bc2_enabled = true; return *this; }

    ModelBoundaryConditions & boundaryConditions2()
    {
      if ( !M_bc2_enabled )
        throw std::logic_error("BoundaryConditions2 are not enabled, call enableBoundaryConditions2() first");
      return M_bc2;
    }
    ModelBoundaryConditions const& boundaryConditions2() const { return M_bc2; }

    ModelBoundaryConditionsNEW const& boundaryConditions() const { return M_bc; }

    ModelInitialConditions & initialConditions() { return M_ic; }
    ModelInitialConditions const& initialConditions() const { return M_ic; }

    ModelPostprocess& postProcess() { return M_postproc; }
    ModelPostprocess const& postProcess() const { return M_postproc; }

    ModelOutputs & outputs() { return M_outputs; }
    ModelOutputs const& outputs() const { return M_outputs; }
  private :
    static nl::json read_json( std::string const& filename, worldcomm_ptr_t const& world );

    void setupImpl();

private:
    std::string M_prefix, M_directoryLibExpr;
    nl::json M_jsonData;
    nl::json M_json_merge_patch;
    nl::json::array_t M_json_patch;

    std::string M_name, M_shortname, M_description, M_unit;
    ModelModels M_models;
    ModelParameters M_params;
    ModelMaterials M_mat;
    ModelInitialConditions M_ic;
    bool M_bc2_enabled = false;
    ModelBoundaryConditions M_bc2;
    ModelBoundaryConditionsNEW M_bc;
    ModelPostprocess M_postproc;
    ModelOutputs M_outputs;
};


}
#endif
