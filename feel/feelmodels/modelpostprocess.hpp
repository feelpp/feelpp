/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 11 Apr 2015
 
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
#ifndef FEELPP_MODELPOSTPROCESS_HPP
#define FEELPP_MODELPOSTPROCESS_HPP 1

#include <map>

#include <boost/property_tree/ptree.hpp>

#include <feel/feelvf/ginac.hpp>

namespace Feel
{

namespace pt = boost::property_tree;

class ModelPointPosition
{
  public:
    typedef Eigen::MatrixXd coord_value_type;
    typedef vector_field_expression<3> coord_expr_type;

    ModelPointPosition()
        : M_value( coord_value_type::Zero( 3, 1 ) )
    {
    }

    ModelPointPosition( ModelPointPosition const& ) = default;
    ModelPointPosition( ModelPointPosition&& ) = default;
    ModelPointPosition& operator=( ModelPointPosition const& ) = default;
    ModelPointPosition& operator=( ModelPointPosition&& ) = default;

    std::string const& name() const { return M_name; }
    coord_value_type const& value() const { return M_value; }
    std::string const& meshMarker() const { return M_meshMarker; }
    coord_expr_type const& expression() const
    {
        CHECK( this->hasExpression() ) << "no expression defined";
        return *M_expr;
    }
    bool hasExpression() const { return M_expr.get_ptr() != 0; }

    void setName( std::string const& s ) { M_name = s; }
    void setValue( coord_value_type const& v ) { M_value = v; }
    void setMeshMarker( std::string const& s ) { M_meshMarker = s; }
    void setExpression( std::string const& expression, std::string const& dirLibExpr = "",
                        WorldComm const& world = Environment::worldComm() )
    {
        M_expr = expr<3, 1>( expression, "", world, dirLibExpr );
        M_value = M_expr->evaluate();
    }

    void setParameterValues( std::map<std::string, double> const& mp )
    {
        if ( !this->hasExpression() )
            return;
        M_expr->setParameterValues( mp );
        M_value = M_expr->evaluate();
    }

  private:
    std::string M_name;
    coord_value_type M_value;
    boost::optional<coord_expr_type> M_expr;
    std::string M_meshMarker;
};

class ModelPostprocessPointPosition : public std::pair<ModelPointPosition, std::set<std::string>>
{
    typedef std::pair<ModelPointPosition, std::set<std::string>> super_type;

  public:
    ModelPostprocessPointPosition( WorldComm const& world = Environment::worldComm() )
        : super_type(),
          M_worldComm( world )
    {
    }
    ModelPostprocessPointPosition( ModelPointPosition const& ptPos, WorldComm const& world = Environment::worldComm() )
        : super_type( ptPos, std::set<std::string>() ),
          M_worldComm( world )
    {
    }
    ModelPostprocessPointPosition( ModelPointPosition const& ptPos, std::set<std::string> const& fields,
                                   WorldComm const& world = Environment::worldComm() )
        : super_type( ptPos, fields ),
          M_worldComm( world )
    {
    }
    ModelPostprocessPointPosition( ModelPostprocessPointPosition const& ) = default;
    ModelPostprocessPointPosition( ModelPostprocessPointPosition&& ) = default;
    ModelPostprocessPointPosition& operator=( ModelPostprocessPointPosition const& ) = default;
    ModelPostprocessPointPosition& operator=( ModelPostprocessPointPosition&& ) = default;

    ModelPointPosition const& pointPosition() const { return this->first; }
    ModelPointPosition& pointPosition() { return this->first; }
    std::set<std::string> const& fields() const { return this->second; }
    std::set<std::string>& fields() { return this->second; }

    void setPTree( pt::ptree const& _p, std::string const& name )
    {
        M_p = _p;
        this->setup( name );
    }
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setFields( std::set<std::string> const& fields ) { this->second = fields; }
    void addFields( std::string const& field ) { this->second.insert( field ); }
    void setParameterValues( std::map<std::string, double> const& mp ) { this->pointPosition().setParameterValues( mp ); }
  private:
    void setup( std::string const& name );

  private:
    WorldComm const& M_worldComm;
    pt::ptree M_p;
    std::string M_directoryLibExpr;
};

class ModelExtremum
{
  public:
    ModelExtremum() = default;
    ModelExtremum( ModelExtremum const& ) = default;
    ModelExtremum( ModelExtremum&& ) = default;
    ModelExtremum& operator=( ModelExtremum const& ) = default;
    ModelExtremum& operator=( ModelExtremum&& ) = default;

    std::string const& name() const { return M_name; }
    std::string const& type() const { return M_type; }
    std::list<std::string> const& meshMarkers() const { return M_meshMarkers; }

    void setName( std::string const& s ) { M_name = s; }
    void setType( std::string const& s )
    {
        CHECK( s == "max" || s == "min" ) << "invalid type " << s;
        M_type = s;
    }
    void addMarker( std::string const& mark )
    {
        if ( std::find( M_meshMarkers.begin(), M_meshMarkers.end(), mark ) == M_meshMarkers.end() )
            M_meshMarkers.push_back( mark );
    }

  private:
    std::string M_name;
    std::string M_type;
    std::list<std::string> M_meshMarkers;
};
class ModelPostprocessExtremum : public std::pair<ModelExtremum, std::set<std::string>>
{
  public:
    ModelPostprocessExtremum( WorldComm const& world = Environment::worldComm() )
        : M_worldComm( world )
    {
    }
    ModelPostprocessExtremum( ModelPostprocessExtremum const& ) = default;
    ModelPostprocessExtremum( ModelPostprocessExtremum&& ) = default;
    ModelPostprocessExtremum& operator=( ModelPostprocessExtremum const& ) = default;
    ModelPostprocessExtremum& operator=( ModelPostprocessExtremum&& ) = default;

    ModelExtremum const& extremum() const { return this->first; }
    ModelExtremum& extremum() { return this->first; }
    std::set<std::string> const& fields() const { return this->second; }
    std::set<std::string>& fields() { return this->second; }

    void setPTree( pt::ptree const& _p, std::string const& name )
    {
        M_p = _p;
        this->setup( name );
    }
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setFields( std::set<std::string> const& fields ) { this->second = fields; }
    void addFields( std::string const& field ) { this->second.insert( field ); }
    //void setParameterValues( std::map<std::string,double> const& mp ) { this->extremum().setParameterValues( mp ); }

  private:
    void setup( std::string const& name );

  private:
    WorldComm const& M_worldComm;
    pt::ptree M_p;
    std::string M_directoryLibExpr;
};

class ModelPostprocess : public std::map<std::string, std::vector<std::string>>
{
  public:
    ModelPostprocess( WorldComm const& world = Environment::worldComm() );
    ModelPostprocess( pt::ptree const& p, WorldComm const& world = Environment::worldComm() );
    virtual ~ModelPostprocess();
    pt::ptree const& pTree() const { return M_p; }
    pt::ptree& pTree() { return M_p; }
    std::vector<ModelPostprocessPointPosition> const& measuresPoint() const { return M_measuresPoint; }
    std::vector<ModelPostprocessExtremum> const& measuresExtremum() const { return M_measuresExtremum; }

    void setPTree( pt::ptree const& _p );
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

    void setParameterValues( std::map<std::string, double> const& mp );

    std::map<std::string, double> toPostprocessValues() const
    {
        std::map<std::string, double> pv;
#if 0
            for( auto const& p : *this )
                pv[p.first]=p.second.value();
#endif
        return pv;
    }
    void saveMD( std::ostream& os );

  private:
    void setup();

  private:
    WorldComm const& M_worldComm;
    pt::ptree M_p;
    std::vector<ModelPostprocessPointPosition> M_measuresPoint;
    std::vector<ModelPostprocessExtremum> M_measuresExtremum;
    std::string M_directoryLibExpr;
};
}
#endif
