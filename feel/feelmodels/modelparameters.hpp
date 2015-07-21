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

namespace Feel {

namespace pt =  boost::property_tree;

struct ModelParameter
{
    ModelParameter() = default;
    ModelParameter( ModelParameter const& ) = default;
    ModelParameter( ModelParameter&& ) = default;
    ModelParameter& operator=( ModelParameter const& ) = default;
    ModelParameter& operator=( ModelParameter && ) = default;
    ModelParameter( std::string name, double value, double min = 0., double max = 0. )
        :
        M_name( name ),
        M_value( value ),
        M_min( min ),
        M_max( max )
        {}
    std::string const& name() const { return M_name; }
    void setName( std::string const& name ) { M_name = name; }
    
    double value() const { return M_value; }
    void setValue( double v ) { M_value = v; }
    double min() const { return M_min; }
    void setMin( double v ) { M_min = v; }
    double max() const { return M_max; }
    void setMax( double v ) { M_max = v; }
    
private:
    std::string M_name;
    double M_value, M_min, M_max;
    
};
class ModelParameters: public std::map<std::string,ModelParameter>
{
public:
    ModelParameters();
    ModelParameters( pt::ptree const& p );
    virtual ~ModelParameters();
    void setPTree( pt::ptree const& _p );
    std::map<std::string,double> toParameterValues() const
        {
            std::map<std::string,double> pv;
            for( auto const& p : *this )
                pv[p.first]=p.second.value();
            return pv;
        }
    
   void saveMD(std::ostream &os); 
private:
    void setup();
private:
    pt::ptree M_p;
};


}
#endif
