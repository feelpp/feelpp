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
//! This file provides the implement of the model post-processing desciption
//!
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @author Vincent Chabannes <vincent.chabannes@feelpp.org>
//! @date 11 Apr 2015
//! @copyright 2015-2017 Feel++ Consortium
//!
#include <iostream>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>

#include <feel/feelmodels/modelpostprocess.hpp>

namespace Feel {

template <typename T>
std::vector<T> as_vector(pt::ptree const& pt, pt::ptree::key_type const& key)
{
    std::vector<T> r;
    for (auto& item : pt.get_child(key))
        r.push_back(item.second.template get_value<T>());
    return r;
}

void
ModelPostprocessPointPosition::setup( std::string const& name )
{

    // fields is necessary
    if ( !M_p.get_child_optional("fields") )
        return;

    bool hasCoord = (M_p.get_child_optional("coord"))?true:false;
    bool hasMarker = (M_p.get_child_optional("marker"))?true:false;
    // coord or marker is necessary
    if ( !hasCoord && !hasMarker )
        return;


    this->pointPosition().setName( name );
     
    if ( hasMarker )
    {
        std::string marker = M_p.get<std::string>( "marker" );
        this->pointPosition().setMeshMarker( marker );
    }
    if ( hasCoord )
    {
        std::string coordExpr = M_p.get<std::string>( "coord" );
        //std::cout << "coordExpr : "<< coordExpr << "\n";

        auto parseExpr = GiNaC::parse(coordExpr);
        auto const& coordMatExpr = parseExpr.first.evalm();
        auto const& coordExprSymbol = parseExpr.second;
        int nComp = coordMatExpr.nops();
        ModelPointPosition::coord_value_type coordData = ModelPointPosition::coord_value_type::Zero(3,1);
        //for ( auto symbb : coordExprSymbol )
        //    std::cout << "symbb " << symbb << "\n";
        if ( coordExprSymbol.empty() || ( coordExprSymbol.size() == 1 && coordExprSymbol[0].get_name() == "0" ) )
        {
            int nComp = GiNaC::parse(coordExpr).first.evalm().nops();
            //std::cout << "ncomp " << nComp << "\n";
            for (int comp=0;comp<nComp;++comp )
            {
                std::string compExpr = str( coordMatExpr.op(comp) );
                try
                {
                    coordData(comp) = std::stod( compExpr );
                }
                catch (std::invalid_argument& err)
                {
                    LOG(WARNING) << "cast fail from expr to double\n";
                    coordData(comp) = 0;
                }
            }
            LOG(INFO) << "point coord is a cst expr : " << coordData(0) << "," << coordData(1) << "," << coordData(2);
            this->pointPosition().setValue( coordData );
        }
        else
        {
            LOG(INFO) << "point coord is a symbolic expr : " << coordExpr;
            this->pointPosition().setExpression( coordExpr, M_directoryLibExpr, M_worldComm );
        }
    } // hasCoord

    // store fields name
    std::vector<std::string> fieldList = as_vector<std::string>( M_p, "fields" );
    if ( fieldList.empty() )
    {
        std::string fieldUnique = M_p.get<std::string>( "fields" );
        if ( !fieldUnique.empty() )
            fieldList = { fieldUnique };
    }
    for( std::string const& field : fieldList )
    {
        // std::cout << "add field = " << field << "\n";
        this->addFields( field );
    }

}

void
ModelPostprocessExtremum::setup( std::string const& name )
{
    // fields is necessary
    if ( !M_p.get_child_optional("fields") )
        return;
    // markers is necessary
    if ( !M_p.get_child_optional("markers") )
        return;

    this->extremum().setName( name );

    // store fields name
    std::vector<std::string> fieldList = as_vector<std::string>( M_p, "fields" );
    if ( fieldList.empty() )
    {
        std::string fieldUnique = M_p.get<std::string>( "fields" );
        if ( !fieldUnique.empty() )
            fieldList = { fieldUnique };
    }
    for( std::string const& field : fieldList )
    {
        //std::cout << "add extremum field = " << field << " (with name " << name << ")\n";
        this->addFields( field );
    }

    // store markers
    std::vector<std::string> markerList = as_vector<std::string>( M_p, "markers" );
    if ( markerList.empty() )
    {
        std::string markerUnique = M_p.get<std::string>( "markers" );
        if ( !markerUnique.empty() )
            markerList = { markerUnique };
    }
    for( std::string const& marker : markerList )
    {
        //std::cout << "add extremum marker = " << marker << " (with name " << name << ")\n";
        this->extremum().addMarker( marker );
    }


}


ModelPostprocess::ModelPostprocess( WorldComm const& world )
    :
    M_worldComm( world )
{}

ModelPostprocess::ModelPostprocess(pt::ptree const& p, WorldComm const& world )
    :
    M_worldComm( world )
{}

ModelPostprocess::~ModelPostprocess()
{}

void
ModelPostprocess::setPTree( pt::ptree const& p )
{
    M_p = p;
    setup();
}

void
ModelPostprocess::setup()
{
    auto fields = M_p.get_child_optional("Fields");
    if ( fields )
    {
        for (auto i : as_vector<std::string>(M_p, "Fields"))
        {
            this->operator[]("Fields").push_back( i );
            LOG(INFO) << "add to postprocess field  " << i;
        }
    }
    auto forces = M_p.get_child_optional("Force");
    if ( forces )
    {
        for (auto i : as_vector<std::string>(M_p, "Force"))
        {
            this->operator[]("Force").push_back( i );
            LOG(INFO) << "add to postprocess force  " << i;
        }
    }
    auto stresses = M_p.get_child_optional("Stresses");
    if ( stresses )
    {
        for (auto i : as_vector<std::string>(M_p, "Stresses"))
        {
            this->operator[]("Stresses").push_back( i );
            LOG(INFO) << "add to postprocess stresses  " << i;
        }
    }

    auto measures = M_p.get_child_optional("Measures");
    if ( measures )
    {
        auto evalPointsFile = measures->get_child_optional("File");
        if ( evalPointsFile )
        {
            // store fields name
            std::vector<std::string> fieldList = as_vector<std::string>( M_p, "Measures.Fields" );
            if ( fieldList.empty() )
            {
                std::string fieldUnique = M_p.get<std::string>( "fields" );
                if ( !fieldUnique.empty() )
                    fieldList = { fieldUnique };
            }

            std::string file = M_p.get<std::string>("Measures.File"); 
            std::ifstream in(Environment::expand(file));
            if(in.good())
            {
                double d = 0.;
                int comp = 0;
                std::string lineData;
                /*
                 * for each point
                 * - new ModelPostProcessPointPosition
                 * - add coord & name to it
                 */
                while(getline(in, lineData)) 
                {
                    ModelPostprocessPointPosition myPpPtPos( M_worldComm );
                    ModelPointPosition::coord_value_type coordData = ModelPointPosition::coord_value_type::Zero(3,1);
                    d = 0.;
                    comp = 0;
                    std::stringstream lineStream(lineData);
                    while(lineStream >> d) // Read Point
                    {
                        coordData(comp) = d;
                        comp++;
                    }
                    myPpPtPos.pointPosition().setValue( coordData );
                    myPpPtPos.pointPosition().setName( "Dummy" );
                    for( std::string const& field : fieldList )
                        myPpPtPos.addFields( field );
                    M_measuresPoint.push_back( myPpPtPos );
                }
                if ( !M_measuresPoint.empty() )
                    this->operator[]("Measures").push_back( "Points" );
            } // good
            else
            {
                LOG(ERROR) << "Unable to open " << Environment::expand(file) << "\n";
            }
        }
        auto evalPoints = measures->get_child_optional("Points");
        if ( evalPoints )
        {
            for( auto const& evalPoint : *evalPoints )
            {
                ModelPostprocessPointPosition myPpPtPos( M_worldComm );
                myPpPtPos.setDirectoryLibExpr( M_directoryLibExpr );
                myPpPtPos.setPTree( evalPoint.second, evalPoint.first );
                if ( !myPpPtPos.fields().empty() )
                    M_measuresPoint.push_back( myPpPtPos );
            }
            if ( !M_measuresPoint.empty() )
                this->operator[]("Measures").push_back( "Points" );
        }

        for ( std::string const& extremumType : std::vector<std::string>( { "Maximum","Minimum" } ) )
        {
            auto measuresExtremum = measures->get_child_optional( extremumType );
            if ( measuresExtremum )
            {
                for( auto const& measureExtremum : *measuresExtremum )
                {
                    ModelPostprocessExtremum myPpExtremum( M_worldComm );
                    if ( extremumType == "Maximum" )
                        myPpExtremum.extremum().setType( "max" );
                    else
                        myPpExtremum.extremum().setType( "min" );
                    myPpExtremum.setDirectoryLibExpr( M_directoryLibExpr );
                    myPpExtremum.setPTree( measureExtremum.second, measureExtremum.first );
                    if ( !myPpExtremum.fields().empty() )
                        M_measuresExtremum.push_back( myPpExtremum );
                }
                if ( !M_measuresPoint.empty() )
                    this->operator[]("Measures").push_back( "Maximum" );
            }
        }
    }


}

void
ModelPostprocess::saveMD(std::ostream &os)
{
  os << "### PostProcess\n";
  os << "|Kind | data |\n";
  os << "|---|---|\n";
  for(auto it = this->begin(); it != this->end(); it++)
  {
    os << "|" << it->first;
    os << "|<ul>";
    for(auto iit = it->second.begin(); iit != it->second.end(); iit++)
      os << "<li>" << *iit << "</li>";
    os << "</ul>|\n";
  }
}

void
ModelPostprocess::setParameterValues( std::map<std::string,double> const& mp )
{
    for( auto & p : M_measuresPoint )
        p.setParameterValues( mp );
}


}
