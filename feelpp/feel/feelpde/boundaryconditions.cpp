/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 25 Jan 2015
 
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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelpde/boundaryconditions.hpp>

namespace Feel
{

BoundaryConditions::BoundaryConditions( WorldComm const& world, bool tryLoadBcFile )
    :
    BoundaryConditions( "", world, tryLoadBcFile )
{}

BoundaryConditions::BoundaryConditions( std::string const& p, WorldComm const& world, bool tryLoadBcFile )
    :
    super(),
    M_worldComm( world ),
    M_prefix( p )
{
    if ( tryLoadBcFile )
    {
        fs::path bc( Environment::expand( soption("bc-file") ) );

        if ( fs::exists( bc ) )
        {
            LOG(INFO) << "Loading Boundary Condition file " << bc.string();
            load( bc.string() );
        }
        else
        {
            LOG(WARNING) << "Boundary condition file " << bc.string() << " does not exist";
        }
    }
}

void
BoundaryConditions::load(const std::string &filename)
{
    // Create an empty property tree object
    using boost::property_tree::ptree;

    read_json(filename, M_pt);
    setup();
}
void
BoundaryConditions::setPTree( pt::ptree const& p )
{
    M_pt = p;
    setup();
}


void
BoundaryConditions::setup()
{
    for( auto const& v : M_pt )
    {
        
        //std::cout << "v.first:" << v.first  << "\n";
        std::string t = v.first; // field name
        for( auto const& f : v.second )
        {
            std::string k = t+"."+f.first; // condition type
            for( auto const& c : f.second ) // condition
            {
                auto bcdatatype  = c.second.get("type","expression");
                std::set<std::string> markers;
                if ( auto meshMarkers = c.second.get_child_optional("markers") )
                {
                    for( auto const& item : c.second.get_child("markers") )
                        markers.insert(item.second.template get_value<std::string>());
                    if( markers.empty() )
                        markers.insert(c.second.template get<std::string>("markers") );
                }
                if( markers.empty() )
                    if( !c.first.empty() ) // for volumic forces
                        markers.insert(c.first);
                int number = c.second.get("number", 1);
                //std::cout << "bcdatatype = " << bcdatatype << std::endl;
                if ( bcdatatype == "file" )
                {
                    try
                    {
                        auto e= c.second.get<std::string>("filename");
                        auto abscissa= c.second.get<std::string>("abscissa");
                        auto ordinate= c.second.get<std::string>("ordinate");
                        LOG(INFO) << "adding boundary " << c.first << " with filename " << e << " to " << k;
                        this->operator[](t)[f.first].push_back( ExpressionStringAtMarker(std::make_tuple( bcdatatype, c.first, e, abscissa, ordinate ), markers, number) );
                    }
                    catch( ... )
                    {
                        LOG(INFO) << "adding boundary " << c.first << " without filename" << " to " << k;
                        throw std::logic_error("invalid boundary conditioner");
                    }
                }
                if ( bcdatatype == "expression" )
                {
                    std::string mat = c.second.get("material", "");
                    try
                    {
                        auto e= c.second.get<std::string>("expr");
                        LOG(INFO) << "adding boundary " << c.first << " with expression " << e << " to " << k;
                        this->operator[](t)[f.first].push_back( ExpressionStringAtMarker(std::make_tuple( bcdatatype, c.first, e, std::string(""), mat ), markers, number) );
                    }
                    catch( ... )
                    {
                        try
                        {
                            auto e1= c.second.get<std::string>("expr1");
                            auto e2= c.second.get<std::string>("expr2");
                            LOG(INFO) << "adding boundary " << c.first << " with expressions " << e1 << " and " << e2 << " to " << k;
                            this->operator[](t)[f.first].push_back( ExpressionStringAtMarker(std::make_tuple( bcdatatype, c.first, e1, e2, mat ), markers, number) );
                        }
                        catch( ... )
                        {
                            LOG(INFO) << "adding boundary " << c.first << " without expression" << " to " << k;
                            this->operator[]( t )[f.first].push_back( ExpressionStringAtMarker(std::make_tuple( bcdatatype, c.first, std::string(""), std::string(""), mat), markers, number ) );
                        }
                    }
                }
            }
        }
        
    }
    if ( Environment::isMasterRank() )
    {
        for( auto const& s : *this )
        {
            LOG(INFO) << "field " << s.first << "\n";
            for( auto const& t : s.second )
            {
                LOG(INFO) << " - type " << t.first << "\n";
                for( auto const& c : t.second )
                {
                    if ( c.hasExpression2() )
                        LOG(INFO) << "  . boundary  " << c.name() << " expr : " << c.expression1() << " expr2:" << c.expression2() << "\n";
                    else
                        LOG(INFO) << "  . boundary  " << c.name() << " expr : " << c.expression() << "\n";
                }
                
            }
        }
    }
}
    

void
BoundaryConditions::saveMD(std::ostream &os)
{
  os << "### Boundary Conditions\n";
  os << "|Name|Type|Expressions|\n";
  os << "|---|---|---|\n";
  for (auto it = this->begin(); it!= this->end(); it++)
  {
    os << "|**" << it->first << "**"; // Var name
    for(auto iit = it->second.begin(); iit !=  it->second.end(); iit++)
    {
     os << "|" << iit->first; // Type
     os << "|<ul>";
     for(auto iiit = iit->second.begin(); iiit !=  iit->second.end(); iiit++)
     {
       os << "<li>**" << iiit->name()      << "**</li>";
       for( auto const& m : iiit->markers() )
           os << "<li>" << m << "</li>";
       os << "<li>" << iiit->expression()  << "</li>";
       os << "<li>" << iiit->expression1() << "</li>";
       os << "<li>" << iiit->expression2() << "</li>";
     }
    }
    os << "</ul>|\n";
  }
  os << "\n";
}


std::pair<bool,int>
BoundaryConditions::iparam( std::string const& field, std::string const& bc, std::string const& marker, std::string const& param ) const
{
    return this->param<int>( field,bc,marker,param, int(0) );
}
std::pair<bool,double>
BoundaryConditions::dparam( std::string const& field, std::string const& bc, std::string const& marker, std::string const& param ) const
{
    return this->param<double>( field,bc,marker,param, double(0) );
}
std::pair<bool,std::string>
BoundaryConditions::sparam( std::string const& field, std::string const& bc, std::string const& marker, std::string const& param ) const
{
    return this->param<std::string>( field,bc,marker,param, std::string("") );
}
template <typename CastType>
std::pair<bool,CastType>
BoundaryConditions::param( std::string const& field,std::string const& bc, std::string const& marker, std::string const& param, CastType const& defaultValue ) const
{
    for( auto const& v : M_pt )
    {
        std::string t = v.first; // field name
        if ( t != field ) continue;
        for( auto const& f : v.second )
        {
            std::string k = t+"."+f.first; // condition type
            if ( f.first != bc ) continue;
            for( auto const& c : f.second ) // marker
            {
                if ( c.first != marker ) continue;
                try
                {
                    CastType e= c.second.get<CastType>( param );
                    return std::make_pair(true,e);
                }
                catch( ... )
                {
                    continue;
                }
            }
        }
    }
    return std::make_pair(false,defaultValue);
}

std::set<std::string>
BoundaryConditions::markers( std::string const& field, std::string const& type ) const
{
    return this->markers( { { field,type } } );
}

std::set<std::string>
BoundaryConditions::markers( std::initializer_list< std::pair<std::string,std::string > > const& listKeys ) const
{
    std::set<std::string> res;
    for ( auto const& key : listKeys )
    {
        std::string const& field = key.first;
        std::string const& type = key.second;
        auto const& itFindField = this->find(field);
        if ( itFindField == this->end() )
            continue;
        auto const& itFindType = itFindField->second.find(type);
        if ( itFindType == itFindField->second.end() )
            continue;
        for ( auto const& f : itFindType->second )
            for( auto const& m : f.markers() )
                res.insert(m);
    }
    return res;
}

}//Feel
