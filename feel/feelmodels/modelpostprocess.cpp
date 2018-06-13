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
ModelPostprocessExports::setup( pt::ptree const& p )
{
    auto fields = p.get_child_optional("fields");
    if ( fields )
    {
        for ( std::string const& fieldName : as_vector<std::string>(p, "fields"))
        {
            M_fields.insert( fieldName );
            LOG(INFO) << "add to postprocess field  " << fieldName;
        }
    }
    if ( auto formatOpt = p.get_optional<std::string>( "format" ) )
        M_format = *formatOpt;
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

void
ModelPostprocessNorm::setup( std::string const& name )
{
    M_name = name;

    if ( auto itField = M_p.get_optional<std::string>("field") )
        M_field = *itField;
    else if ( auto ptexpr = M_p.get_child_optional("expr") )
    {
        M_expr.setExpr( "expr", M_p, M_worldComm, M_directoryLibExpr );
        if ( auto ptgradexpr = M_p.get_child_optional("grad_expr") )
            M_gradExpr.setExpr( "grad_expr", M_p, M_worldComm, M_directoryLibExpr );
    }

    if ( auto ptmarkers = M_p.get_child_optional("markers") )
    {
        for( auto const& item : M_p.get_child("markers") )
            M_markers.insert(item.second.template get_value<std::string>());
        if( M_markers.empty() )
            M_markers.insert(M_p.get<std::string>("markers") );
    }

    if ( auto pttype = M_p.get_child_optional("type") )
    {
        for( auto const& item : M_p.get_child("type") )
            M_types.insert(item.second.template get_value<std::string>());
        if( M_types.empty() )
            M_types.insert(M_p.get<std::string>("type") );
    }

    if ( auto itSol = M_p.get_optional<std::string>("solution") )
        M_solution.setExpr( "solution", M_p, M_worldComm, M_directoryLibExpr );

    if ( auto itSol = M_p.get_optional<std::string>("grad_solution") )
        M_gradSolution.setExpr( "grad_solution", M_p, M_worldComm, M_directoryLibExpr );

    if ( auto itQuad = M_p.get_optional<int>("quad") )
    {
        M_quadOrder = *itQuad;
        if ( auto itQuad1 = M_p.get_optional<int>("quad1") )
            M_quad1Order = *itQuad1;
        else
            M_quad1Order = M_quadOrder;
    }
    else if ( auto itQuad1 = M_p.get_optional<int>("quad1") )
        M_quad1Order = *itQuad1;
}

void
ModelPostprocessNorm::setParameterValues( std::map<std::string,double> const& mp )
{
    M_solution.setParameterValues( mp );
    M_gradSolution.setParameterValues( mp );
}


ModelPostprocess::ModelPostprocess( WorldComm const& world )
    :
    M_worldComm( world ),
    M_useModelName( false )
{}

ModelPostprocess::ModelPostprocess(pt::ptree const& p, WorldComm const& world )
    :
    M_worldComm( world ),
    M_useModelName( false )
{}

ModelPostprocess::~ModelPostprocess()
{}

pt::ptree
ModelPostprocess::pTree( std::string const& name ) const
{
    if ( !M_useModelName )
        return M_p;
    else if ( auto ptreeWithName = M_p.get_child_optional( name ) )
        return *ptreeWithName;
    pt::ptree ptree;
    return ptree;
}

void
ModelPostprocess::setPTree( pt::ptree const& p )
{
    M_p = p;
    setup();
}


void
ModelPostprocess::setup()
{
    if ( auto useModelName = M_p.get_optional<bool>("use-model-name") )
        M_useModelName = *useModelName;

    if ( M_useModelName )
    {
        for( auto const& p1 : M_p )
        {
            this->setup( p1.first,p1.second );
        }
    }
    else
        this->setup( "", M_p );
}
void
ModelPostprocess::setup( std::string const& name, pt::ptree const& p  )
{
    if ( auto exports = p.get_child_optional("Exports") )
    {
        ModelPostprocessExports ppexports;
        ppexports.setup( *exports );
        if ( !ppexports.fields().empty() )
            M_exports[name] = ppexports;
    }

    if ( auto measures = p.get_child_optional("Measures") )
    {
        auto evalPoints = measures->get_child_optional("Points");
        if ( evalPoints )
        {
            for( auto const& evalPoint : *evalPoints )
            {
                ModelPostprocessPointPosition myPpPtPos( M_worldComm );
                myPpPtPos.setDirectoryLibExpr( M_directoryLibExpr );
                myPpPtPos.setPTree( evalPoint.second, evalPoint.first );
                if ( !myPpPtPos.fields().empty() )
                    M_measuresPoint[name].push_back( myPpPtPos );
            }
        }

        for ( std::string const& extremumType : std::vector<std::string>( { "Maximum","Minimum","Mean" } ) )
        {
            auto measuresExtremum = measures->get_child_optional( extremumType );
            if ( measuresExtremum )
            {
                for( auto const& measureExtremum : *measuresExtremum )
                {
                    ModelPostprocessExtremum myPpExtremum( M_worldComm );
                    if ( extremumType == "Maximum" )
                        myPpExtremum.extremum().setType( "max" );
                    else if ( extremumType == "Minimum" )
                        myPpExtremum.extremum().setType( "min" );
                    else if ( extremumType == "Mean" )
                        myPpExtremum.extremum().setType( "mean" );
                    myPpExtremum.setDirectoryLibExpr( M_directoryLibExpr );
                    myPpExtremum.setPTree( measureExtremum.second, measureExtremum.first );
                    if ( !myPpExtremum.fields().empty() )
                        M_measuresExtremum[name].push_back( myPpExtremum );
                }
            }
        }

        auto ptreeNorms = measures->get_child_optional("Norm");
        if ( ptreeNorms )
        {
            for( auto const& ptreeNorm : *ptreeNorms )
            {
                ModelPostprocessNorm ppNorm( M_worldComm );
                ppNorm.setDirectoryLibExpr( M_directoryLibExpr );
                ppNorm.setPTree( ptreeNorm.second, ptreeNorm.first );
                if ( ppNorm.hasField() || ppNorm.hasExpr() )
                    M_measuresNorm[name].push_back( ppNorm );
            }
        }

    }
}



void
ModelPostprocess::saveMD(std::ostream &os)
{
#if 0
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
#endif
}

void
ModelPostprocess::setParameterValues( std::map<std::string,double> const& mp )
{
    for( auto & p : M_measuresPoint )
        for( auto & p2 : p.second )
            p2.setParameterValues( mp );

    for( auto & p : M_measuresNorm )
        for( auto & p2 : p.second )
            p2.setParameterValues( mp );
}

bool
ModelPostprocess::hasExports( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return M_exports.find( nameUsed ) != M_exports.end();
}
bool
ModelPostprocess::hasMeasuresPoint( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return M_measuresPoint.find( nameUsed ) != M_measuresPoint.end();
}
bool
ModelPostprocess::hasMeasuresExtremum( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return M_measuresExtremum.find( nameUsed ) != M_measuresExtremum.end();
}
bool
ModelPostprocess::hasMeasuresNorm( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return M_measuresNorm.find( nameUsed ) != M_measuresNorm.end();
}
ModelPostprocessExports const&
ModelPostprocess::exports( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    //CHECK( this->hasExports( nameUsed ) ) << "no exports with name:"<<name;
    if ( this->hasExports( nameUsed ) )
        return M_exports.find( nameUsed )->second;
    else
        return M_emptyExports;
}
std::vector<ModelPostprocessPointPosition> const&
ModelPostprocess::measuresPoint( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    //CHECK( this->hasMeasuresPoint( nameUsed ) ) << "no measures point with name:"<<name;
    if ( this->hasMeasuresPoint( nameUsed ) )
        return M_measuresPoint.find( nameUsed )->second;
    else
        return M_emptyMeasuresPoint;
}
std::vector<ModelPostprocessExtremum> const&
ModelPostprocess::measuresExtremum( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    //CHECK( this->hasMeasuresExtremum( nameUsed ) ) << "no measures extremum with name:"<<name;
    if ( this->hasMeasuresExtremum( nameUsed ) )
        return M_measuresExtremum.find( nameUsed )->second;
    else
        return M_emptyMeasuresExtremum;
}
std::vector<ModelPostprocessNorm> const&
ModelPostprocess::measuresNorm( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    //CHECK( this->hasMeasuresNorm( nameUsed ) ) << "no measures norm with name:"<<name;
    if ( this->hasMeasuresNorm( nameUsed ) )
        return M_measuresNorm.find( nameUsed )->second;
    else
        return M_emptyMeasuresNorm;
}



}
