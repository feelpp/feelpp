/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Jean-Baptiste Wahl
            Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 01 Dec 2019

 Copyright (C) 2018-2019 Feel++ Consortium

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
#ifndef FEELPP_FMUMODELBASE_HPP
#define FEELPP_FMUMODELBASE_HPP 1

#include <fmilib.h>

namespace Feel
{
class FmuModelBase
{
public :
    typedef jm_callbacks callbacks_type;
    typedef std::shared_ptr<jm_callbacks> callbacks_ptrtype;

    typedef std::vector<std::string> var_list_type;
    typedef std::shared_ptr<var_list_type> var_list_ptrtype;

    FmuModelBase( callbacks_ptrtype const& callbacks )
        :
        M_callbacks( callbacks ),
        M_allocated_xml( false ),
        M_allocated_dll( false ),
        M_allocated_fmu( false ),
        M_setup( false )
    {}

    virtual ~FmuModelBase()
    {}

    void setExportList( var_list_ptrtype const& exp_list ) { M_export_list = exp_list; }
    void setExportDirectory( std::string const& path ) { M_export_directory=path; }

    int version() const { return M_version; }
    std::string const& name() const { return M_name; }
    std::string const& guid() const  { return M_guid; }
    std::string const& kind() const { return M_kind; }
    bool isSetup() const { return M_setup; }

    virtual void reset()=0;
    virtual void setupExperiment( double const& t_init, double const& t_final, double const& tol )=0;
    virtual void initialize( double t_init )=0;
    virtual void terminate()=0;
    virtual void doStep( double t_cur, double step, bool newStep )=0;

    virtual void printInfo() const =0;
    virtual void printVariablesInfo() const =0;
    virtual void exportValues() const =0;

    virtual double defaultStartTime() const =0;
    virtual double defaultFinalTime() const =0;
    virtual double defaultTolerance() const =0;

    virtual void setValue( std::string name, double value )=0;
    virtual void setValue( std::string name, int value )=0;
    virtual void setValue( std::string name, std::string value )=0;
    virtual void setValue( std::string name, bool value )=0;

    virtual void getValue( std::string name, double& value ) const =0;
    virtual void getValue( std::string name, int& value ) const =0;
    virtual void getValue( std::string name, std::string& value ) const =0;
    virtual void getValue( std::string name, bool& value ) const =0;

protected :
    callbacks_ptrtype M_callbacks;
    bool M_allocated_xml, M_allocated_dll, M_allocated_fmu, M_setup;
    std::string M_name, M_guid, M_id, M_kind, M_export_directory;
    int M_version;

    var_list_ptrtype M_export_list;
};

}

#endif
