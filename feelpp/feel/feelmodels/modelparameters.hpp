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
#ifndef FEELPP_MODELPARAMETERS_HPP
#define FEELPP_MODELPARAMETERS_HPP 1

#include <map>

#include <boost/property_tree/ptree.hpp>

#include <feel/feelcore/commobject.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelfit/interpolator.hpp>

namespace Feel {

namespace pt =  boost::property_tree;

struct FEELPP_EXPORT ModelParameter 
{
    ModelParameter() = default;
    ModelParameter( ModelParameter const& ) = default;
    ModelParameter( ModelParameter&& ) = default;
    ModelParameter& operator=( ModelParameter const& ) = default;
    ModelParameter& operator=( ModelParameter && ) = default;
    ModelParameter( std::string const& name, double value, double min = 0., double max = 0. )
        :
        M_name( name ),
        M_type( "value" ),
        M_value( value ),
        M_min( min ),
        M_max( max )
        {}
    ModelParameter( std::string const& name, std::string const& expression,
                    std::string const& dirLibExpr = "",
                    WorldComm const& world = Environment::worldComm(),
                    double min = 0., double max = 0. )
        :
        M_name( name ),
        M_type( "expression" ),
        M_value( 0. ),
        M_min( min ),
        M_max( max ),
        M_expr( expr<2>( expression,"",world,dirLibExpr ) )
        {
            M_value = M_expr->evaluate();
        }
    ModelParameter( std::string const& name, std::shared_ptr<Interpolator> interpolator, std::string const& expression,
                    std::string const& dirLibExpr = "",
                    WorldComm const& world = Environment::worldComm() )
    :
    M_name( name ),
    M_type( "fit" ),
    M_value( 0. ),
    M_min( 0. ),
    M_max( 0. ),
    M_expr( expr<2>( expression,"",world,dirLibExpr ) ),
    M_interpolator( interpolator )
    {}

    std::string const& name() const { return M_name; }
    void setName( std::string const& name ) { M_name = name; }

    std::string const& type() const { return M_type; }

    double value() const { return M_value; }
    void setValue( double v ) { M_value = v; }
    double min() const { return M_min; }
    void setMin( double v ) { M_min = v; }
    double max() const { return M_max; }
    void setMax( double v ) { M_max = v; }
    bool hasMinMax() const { return M_min != 0 || M_max != 0; }

    bool hasExpression() const { return M_expr.get_ptr() != 0; } //M_expr != boost::none; }
    void setParameterValues( std::map<std::string,double> const& mp )
        {
            if ( !this->hasExpression() )
                return;
            M_expr->setParameterValues( mp );
            M_value = M_expr->evaluate();
        }
    scalar_field_expression<2> const& expression() const { CHECK( this->hasExpression() ) << "no expression defined"; return *M_expr; }

    bool hasFitInterpolator() const { return ( M_interpolator )? true : false; }
    std::shared_ptr<Interpolator> fitInterpolator() const { return M_interpolator; }

private:
    std::string M_name, M_type;
    double M_value, M_min, M_max;
    boost::optional<scalar_field_expression<2>> M_expr;
    std::shared_ptr<Interpolator> M_interpolator;

};

//!
//! class for Model Parameters
//!
class ModelParameters: public std::map<std::string,ModelParameter>, public CommObject
{
public:
    using super=CommObject;
    ModelParameters( worldcomm_ptr_t const& world = Environment::worldCommPtr() );
    ModelParameters( pt::ptree const& p, worldcomm_ptr_t const& world = Environment::worldCommPtr() );
    ModelParameters( ModelParameters const& ) = default;
    virtual ~ModelParameters();
    void setPTree( pt::ptree const& _p );
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& mp );
    std::map<std::string,double> toParameterValues() const;

   void saveMD(std::ostream &os);
private:
    void setup();
private:
    pt::ptree M_p;
    std::string M_directoryLibExpr;
};


}
#endif
