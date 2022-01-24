/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 5 September 2019

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
#ifndef FEELPP_MODELS_INITIALCONDITIONS_HPP
#define FEELPP_MODELS_INITIALCONDITIONS_HPP 1

#include <feel/feelcore/commobject.hpp>
#include <feel/feelmodels/modelmarkers.hpp>
#include <feel/feelmodels/modelexpression.hpp>

namespace Feel {

class FEELPP_EXPORT ModelInitialCondition : public CommObject
{
    using super = CommObject;
  public :
    ModelInitialCondition( worldcomm_ptr_t const& world, std::string const& directoryLibExpr )
        :
        super( world ),
        M_isExpression( false ),
        M_isFile( false ),
        M_directoryLibExpr( directoryLibExpr )
        {}
    ModelInitialCondition( ModelInitialCondition const& ) = default;
    ModelInitialCondition( ModelInitialCondition && ) = default;
    ModelInitialCondition& operator=( ModelInitialCondition const& ) = default;
    ModelInitialCondition& operator=( ModelInitialCondition && ) = default;

    std::string const& name() const { return M_name; }
    void setName( std::string const& name ) { M_name = name; }

    bool isExpression() const { return M_isExpression; }
    bool isFile() const { return M_isFile; }

    ModelExpression const& expression() const { return M_modelExpr; }
    ModelMarkers const& markers() const { return M_markers; }
    std::string const& fileName() const { return M_fileName; }
    std::string const& fileType() const { return M_fileType; }
    std::string const& fileDirectory() const { return M_fileDirectory; }

    nl::json const& jsonSetup() const { return M_jsonSetup; }

    void setParameterValues( std::map<std::string,double> const& mp );
    void setup( nl::json const& jarg, std::string const& typeIC );
  private :
    std::string M_name;
    bool M_isExpression, M_isFile;

    std::string M_directoryLibExpr;
    ModelMarkers M_markers;
    ModelExpression M_modelExpr;
    std::string M_fileName, M_fileType, M_fileDirectory;
    nl::json M_jsonSetup;
};


class FEELPP_EXPORT ModelInitialConditionTimeSet : public CommObject, public std::map<double,std::map<std::string,std::vector<ModelInitialCondition>>>
{
    using super = CommObject;
  public :
    ModelInitialConditionTimeSet( worldcomm_ptr_t const& world = Environment::worldCommPtr(), std::string const& directoryLibExpr = "" )
        :
        super( world ),
        M_directoryLibExpr( directoryLibExpr )
        {}
    ModelInitialConditionTimeSet( ModelInitialConditionTimeSet const& ) = default;
    ModelInitialConditionTimeSet( ModelInitialConditionTimeSet && ) = default;
    ModelInitialConditionTimeSet& operator=( ModelInitialConditionTimeSet const& ) = default;
    ModelInitialConditionTimeSet& operator=( ModelInitialConditionTimeSet && ) = default;

    void setParameterValues( std::map<std::string,double> const& mp );
    void setup( nl::json const& jarg );
  private :
    std::string M_directoryLibExpr;

};

class FEELPP_EXPORT ModelInitialConditions : public CommObject, public std::map<std::string,std::map<std::string,ModelInitialConditionTimeSet>> // toolbox name -> (field name -> ic data)
{
    using super = CommObject;
  public :
    explicit ModelInitialConditions( worldcomm_ptr_t const& world = Environment::worldCommPtr() )
        :
        super( world )
    {}
    ModelInitialConditions( ModelInitialConditions const& ) = default;
    ModelInitialConditions( ModelInitialConditions && ) = default;
    ModelInitialConditions& operator=( ModelInitialConditions const& ) = default;
    ModelInitialConditions& operator=( ModelInitialConditions && ) = default;

    void setPTree( nl::json const& jarg );
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );

    bool has( std::string const& tbname, std::string const& fieldname ) const;
    ModelInitialConditionTimeSet const& get( std::string const& tbname, std::string const& fieldname ) const;
  private :
    void setupInternal( std::string const& tbname, nl::json const& jarg );
  private :
    std::string M_directoryLibExpr;
    ModelInitialConditionTimeSet M_emptylInitialConditionTimeSet;
};


} // namespace Feel

#endif
