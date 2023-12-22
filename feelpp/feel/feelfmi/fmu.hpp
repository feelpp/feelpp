//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @date 21 Sep 2018
//! @copyright 2018 Feel++ Consortium
//!


#ifndef FMIWRAPPER_HPP
#define FMIWRAPPER_HPP

#include <feel/feelcore/environment.hpp>


#if FEELPP_HAS_FMILIB
#include <fmilib.h>

#include <feel/feelfmi/fmumodel2.hpp>
#include <feel/feelfmi/cssolver.hpp>

namespace Feel
{

//!
//! @brief Functional Mockup Unit class
//!
class FMU :
        public std::enable_shared_from_this<FMU>
{
public :
    typedef jm_callbacks callbacks_type;
    typedef std::shared_ptr<jm_callbacks> callbacks_ptrtype;

    typedef FmuModelBase fmumodel_type;
    typedef std::shared_ptr<fmumodel_type> fmumodel_ptrtype;

    typedef SolverBase solver_type;
    typedef std::shared_ptr<solver_type> solver_ptrtype;

    typedef std::vector<std::string> var_list_type;
    typedef std::shared_ptr<var_list_type> var_list_ptrtype;


    FMU( std::string prefix="" );
    ~FMU();

    int load( std::string _path="" );
    void simulate( double t_init=-1, double t_final=-1, double tolerance=-1 );
    void reset();
    void initialize( double t_init=-1, double t_final=-1, double tolerance=-1 );
    void doSteps( double t_stop );

    void initSolver();
    void setSolverTimeStep( double const& step );

    void printModelInfo();
    void addExportedVariables( std::vector<std::string> const& var_list );
    void setExportedVariables( std::vector<std::string> const& var_list );
    void setExportDirectory( std::string const& path );

    template <typename VariableType>
    void setValue( std::string var_name, VariableType value )
    {
        CHECK( M_model ) <<"FMU trying to access variable without model\n";
        M_model->setValue( var_name, value );
    }

    double currentTime();

    template <typename VariableType>
    VariableType getValue( std::string var_name )
    {
        VariableType value;
        CHECK( M_model ) <<"FMU trying to access variable without model\n";
        M_model->getValue( var_name, value );
        return value;
    }

private :
    static void importlogger(jm_callbacks* c,
                             jm_string module, jm_log_level_enu_t log_level, jm_string message)
    {
        LOG(INFO) << "FeelFMI["<< module <<"] : " << message << std::endl;
    }


private :
    std::string M_prefix, M_tmp_dir, M_export_directory;
    bool M_verbose;
    callbacks_ptrtype M_callbacks;
    fmi_import_context_t* M_context;

    fmumodel_ptrtype M_model;
    solver_ptrtype M_solver;

    var_list_ptrtype M_export_list;
    double M_tinit, M_tfinal;
};



} //namespace Feel
#endif
#endif
