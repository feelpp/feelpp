/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 01 Dec 2019

 Copyright (C) 2019 Feel++ Consortium

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
#ifndef FEELPP_FMUMODEL2_HPP
#define FEELPP_FMUMODEL2_HPP 1

#include <fmilib.h>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfmi/fmumodelbase.hpp>

namespace Feel
{
struct Fmi2Variable;

class FmuModel2 :
        public FmuModelBase
{
public :
    typedef FmuModel2 self_type;
    typedef FmuModelBase super_type;

    typedef typename super_type::callbacks_ptrtype callbacks_ptrtype;
    typedef Fmi2Variable variable_type;
    typedef std::shared_ptr<variable_type> variable_ptrtype;
    typedef std::map<std::string,variable_ptrtype> var_map_type;

    FmuModel2( fmi_import_context_t* context, std::string tmp_dir, callbacks_ptrtype callbacks );
    ~FmuModel2();

    void reset() override;
    void setupExperiment( double const& t_init, double const& t_final, double const& tol ) override;
    void initialize( double t_init ) override;
    void terminate() override;
    void doStep( double t_cur, double step, bool newStep ) override;

    void printVariablesInfo() const override;
    void printInfo() const override;
    void exportValues() const override;

    double defaultStartTime() const override;
    double defaultFinalTime() const override;
    double defaultTolerance() const override;

    void setValue( std::string name, double value ) override;
    void setValue( std::string name, int value ) override;
    void setValue( std::string name, std::string value ) override;
    void setValue( std::string name, bool value ) override;

    void getValue( std::string name, double& value ) const override;
    void getValue( std::string name, int& value ) const override;
    void getValue( std::string name, std::string& value ) const override;
    void getValue( std::string name, bool& value ) const override;

    template <typename VariableType>
    VariableType getValue( std::string const& name ) const
    {
        VariableType value;
        getValue( name, value );
        return value;
    }

private :
    static void fmi2logger(fmi2_component_environment_t env, fmi2_string_t instanceName,
                           fmi2_status_t status, fmi2_string_t category, fmi2_string_t message, ...);
    static void stepFinished(fmi2_component_environment_t env, fmi2_status_t status);

    std::string strValue( std::string const& var );

private :
    fmi2_import_t* M_fmu;
    fmi2_callback_functions_t M_callbackfunctions;
    var_map_type M_v_map;
    std::map< std::string,std::vector<std::string> > M_values;
    fmi2_type_t M_type;

};


} // Feel

#endif
