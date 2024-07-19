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

#include <feel/feelmodels/modelboundaryconditions.hpp>

namespace Feel
{

ModelBoundaryCondition::ModelBoundaryCondition( worldcomm_ptr_t const& worldComm )
    : super( worldComm )
{
}

ModelBoundaryCondition::ModelBoundaryCondition( pt::ptree const& p, std::string const& name,
                                                std::string const& material, std::string const& e1,
                                                std::string const& e2,
                                                worldcomm_ptr_t const& worldComm )
    : super( worldComm ),
      M_pt( p ),
      M_name( name ),
      M_material( material ),
      M_markers( name ),
      M_expr1( e1 ),
      M_expr2( e2 )
{
    if ( auto markers = M_pt.get_child_optional( "markers" ) )
        M_markers.setPTree( *markers );

    if ( this->hasExpression() )
        M_modelExpr1.setExpr( M_expr1 );
    if ( this->hasExpression2() )
        M_modelExpr2.setExpr( M_expr2 );
}

void
ModelBoundaryCondition::setParameterValues( std::map<std::string,double> const& mp )
{
    if ( this->hasExpression() )
        M_modelExpr1.setParameterValues( mp );
    if ( this->hasExpression2() )
        M_modelExpr2.setParameterValues( mp );
}

ModelBoundaryConditions::ModelBoundaryConditions( worldcomm_ptr_t const& worldComm )
    : super( worldComm )
{
}

void ModelBoundaryConditions::setPTree( pt::ptree const& p )
{
    M_pt = p;
    setup();
}

void ModelBoundaryConditions::setup()
{
    for ( auto const& v : M_pt )
    {
        std::string field = v.first; // field name
        for ( auto const& f : v.second )
        {
            std::string cond = f.first;
            for ( auto const& c : f.second ) // condition
            {
                std::string name = c.first;
                std::string mat = c.second.get( "material", "" );
                std::string e1 = "";
                std::string e2 = "";
                try
                {
                    e1 = c.second.get<std::string>( "expr" );
                    LOG( INFO ) << "adding boundary " << name << " with expression "
                                << e1 << " to condition of type " << cond << " on " << field;
                }
                catch ( ... )
                {
                    try
                    {
                        e1 = c.second.get<std::string>( "expr1" );
                        e2 = c.second.get<std::string>( "expr2" );
                        LOG( INFO ) << "adding boundary " << name
                                    << " with expressions " << e1 << " and " << e2
                                    << " to condition of type " << cond << " on " << field;
                    }
                    catch ( ... )
                    {
                        LOG( INFO ) << "adding boundary " << name << " without expression "
                                    << "to condition of type " << cond << " on " << field;
                    }
                }
                auto bc = ModelBoundaryCondition( c.second, name, mat, e1, e2 );
                this->operator[]( field )[cond][name] = bc;
            }
        }
    }
}

void
ModelBoundaryConditions::setParameterValues( std::map<std::string,double> const& mp )
{
    for( auto& [field, mapf] : *this )
    {
        for( auto& [type, mapt] : mapf )
        {
            for( auto& [name, bc] : mapt )
            {
                bc.setParameterValues( mp );
            }
        }
    }
}

std::map<ModelBoundaryId, ModelBoundaryCondition>
ModelBoundaryConditions::flatten() const
{
    // cout << " - flatten" << std::endl;
    std::map<ModelBoundaryId, ModelBoundaryCondition> f_;
    for ( auto const& [bcfield, bc1] : *this )
    {
        // cout << "   field : " << bcfield << std::endl;
        for ( auto const& [bctype, bc2] : bc1 )
        {
            // cout << "   type : " << bctype << std::endl;
            for ( auto const& [bcname, bc] : bc2 )
            {
                // cout << "   name : " << bctype << std::endl;
                ModelBoundaryId bcid{ std::tuple{ bcfield, bctype, bcname } };
                // cout << "   * " << bcid << std::endl;
                f_[bcid] = bc;
            }
        }
    }
    return f_;
}

std::map<std::string,std::map<std::string,ModelBoundaryCondition> > const&
ModelBoundaryConditions::byField( std::string const& field ) const
{
    auto it = this->find(field);
    if( it != this->end() )
        return it->second;
    else
        return M_emptyField;
}

std::map<std::string,ModelBoundaryCondition> const&
ModelBoundaryConditions::byFieldType( std::string const& field, std::string const& type ) const
{
    auto it = this->find(field);
    if( it == this->end() )
        return M_emptyFieldType;
    auto it2 = it->second.find(type);
    if( it2 != it->second.end() )
        return it2->second;
    else
        return M_emptyFieldType;
}

} // namespace Feel
