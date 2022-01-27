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
#include <boost/property_tree/json_parser.hpp>

namespace Feel {

void
ModelPostprocessExports::setup( nl::json const& jarg  )
{
    if ( jarg.contains( "fields" ) )
    {
        auto const& j_fields = jarg.at("fields");
        if ( j_fields.is_string() )
            M_fields.insert( j_fields.template get<std::string>() );
        else if ( j_fields.is_array() )
        {
            for ( auto const& el : j_fields.items() )
                if ( el.value().is_string() )
                    M_fields.insert( el.value().get<std::string>() );
        }
    }

    if ( jarg.contains( "format" ) )
    {
        auto const& j_format = jarg.at("format");
        if ( j_format.is_string() )
            M_format = j_format.get<std::string>();
    }

    if ( jarg.contains( "expr" ) )
    {
        auto const& j_expr = jarg.at("expr");
        for ( auto const& [key,jvalue] : j_expr.items() )
        {
            std::string exprName = key;//item.first;
            std::set<std::string> representations, tags;

            if ( jvalue.is_string() || jvalue.is_number() ) // name:expr
            {
                ModelExpression modelexpr;
                ModelMarkers markers;
                modelexpr.setExpr( jvalue, this->worldComm(), M_directoryLibExpr/*,indexes*/ );
                CHECK( modelexpr.hasAtLeastOneExpr() ) << "expr not given correctly";
                M_exprs.push_back( std::make_tuple(exprName,modelexpr,markers,representations,tags) );
            }
            else if ( jvalue.is_object() )
            {
                if ( jvalue.contains("representation") )
                {
                    auto const& j_representation = jvalue.at("representation");
                    if ( j_representation.is_string() )
                        representations.insert( j_representation.get<std::string>() );
                    else if ( j_representation.is_array() )
                    {
                        for ( auto const& [jrepkey,jrepval] : j_representation.items() )
                        {
                            if ( jrepval.is_string() )
                                representations.insert( jrepval.get<std::string>() );
                        }
                    }
                }

                if ( jvalue.contains("tag") )
                {
                    auto const& j_tag = jvalue.at("tag");
                    if ( j_tag.is_string() )
                        tags.insert( j_tag.get<std::string>() );
                    else if ( j_tag.is_array() )
                        for ( auto const& [tagkey,tagval] : j_tag.items() )
                            if ( tagval.is_string() )
                                tags.insert( tagval.get<std::string>() );
                }

                if ( jvalue.contains("expr") )
                {
                    ModelExpression modelexpr;
                    ModelMarkers markers;
                    modelexpr.setExpr( jvalue.at("expr"), this->worldComm(), M_directoryLibExpr/*,indexes*/ );
                    if ( jvalue.contains("markers") )
                        markers.setup( jvalue.at("markers") /*, indexes*/ );
                    CHECK( modelexpr.hasAtLeastOneExpr() ) << "expr not given correctly";
                    M_exprs.push_back( std::make_tuple(exprName,modelexpr,markers,representations,tags) );
                }

                if ( jvalue.contains("parts") )
                {
                    auto const& j_parts = jvalue.at("parts");
                    CHECK ( j_parts.is_array() ) << "parts should an array";
                    for ( auto const& [j_partskey,j_partsval] : j_parts.items() )
                    {
                        ModelExpression modelexpr;
                        ModelMarkers markers;
                        CHECK( j_partsval.contains("expr") ) << "expr is missing";
                        modelexpr.setExpr( j_partsval.at("expr"), this->worldComm(), M_directoryLibExpr/*,indexes*/ );
                        if ( j_partsval.contains("markers") )
                            markers.setup( j_partsval.at("markers") /*, indexes*/ );
                        CHECK( modelexpr.hasAtLeastOneExpr() ) << "expr not given correctly";
                        M_exprs.push_back( std::make_tuple(exprName,modelexpr,markers,representations,tags) );
                    }
                }
            } // if ( jvalue.is_object() )
        }
    }

}

void
ModelPostprocessExports::setParameterValues( std::map<std::string,double> const& mp )
{
    for ( auto & edata : M_exprs )
        std::get<1>( edata ).setParameterValues( mp );
}

void
ModelPostprocessQuantities::setup( nl::json const& jarg )
{
    if ( jarg.is_string() )
         M_quantities.insert( jarg.get<std::string>() );
    else if ( jarg.is_array() )
    {
        for ( auto const& [jargkey,jargval] : jarg.items() )
            if ( jargval.is_string() )
                M_quantities.insert( jargval.get<std::string>() );
    }
    else if ( jarg.is_object() )
    {
        if ( jarg.contains("names") )
        {
            auto const& j_names = jarg.at("names");
            if ( j_names.is_string() )
                M_quantities.insert( j_names.get<std::string>() );
            else if ( j_names.is_array() )
                for ( auto const& [j_nameskey,j_namesval] : j_names.items() )
                    if ( j_namesval.is_string() )
                        M_quantities.insert( j_namesval.get<std::string>() );
        }

        if ( jarg.contains("expr") )
        {
            auto const& j_expr = jarg.at("expr");
            for ( auto const& [j_exprkey,j_exprval] : j_expr.items() )
            {
                std::string const& exprName = j_exprkey;
                CHECK( !exprName.empty() ) << "expr name is empty or expr is an array";
                ModelExpression modelexpr;
                modelexpr.setExpr( j_exprval, this->worldComm(), M_directoryLibExpr/*,indexes*/ );
                CHECK( modelexpr.hasAtLeastOneExpr() ) << "expr not given correctly";
                M_exprs.emplace( std::make_pair( exprName, std::move(modelexpr) ) );
            }
        }
    }
}

void
ModelPostprocessQuantities::setParameterValues( std::map<std::string,double> const& mp )
{
    for ( auto & [_name,_expr] : M_exprs )
        _expr.setParameterValues( mp );
}

void
ModelPostprocessSave::setup( nl::json const& jarg )
{
    if ( jarg.contains( "Fields" ) )
    {
        auto const& j_fields = jarg.at("Fields");
        if ( j_fields.contains("names") )
        {
            auto const& j_fields_names = j_fields.at("names");
            if ( j_fields_names.is_string() )
            {
                M_fieldsNames.insert( j_fields_names.get<std::string>() );
            }
            else if ( j_fields_names.is_array() )
            {
                for ( auto const& [fieldnameskey,fieldnamesval] : j_fields_names.items() )
                    if ( fieldnamesval.is_string() )
                        M_fieldsNames.insert( fieldnamesval.get<std::string>() );
            }
        }
        if ( j_fields.contains("format") )
        {
            auto const& j_fields_format = j_fields.at("format");
            CHECK( j_fields_format.is_string() ) << "format should be a string";
            M_fieldsFormat = j_fields_format.get<std::string>();
        }
    }
}

void
ModelPostprocessPointPosition::PointPosition::setup( nl::json const& jData, worldcomm_t const& world, std::string const& directoryLibExpr, ModelIndexes const& indexes )
{
    CHECK( jData.is_string() ) << "PointPosition json should be a string";
    this->setup( jData.get<std::string>(), world, directoryLibExpr, indexes );
}
void
ModelPostprocessPointPosition::PointPosition::setup( std::string const& coordExprStr, worldcomm_t const& world, std::string const& directoryLibExpr, ModelIndexes const& indexes )
{
    std::string coordExpr = indexes.replace( coordExprStr );

    ModelExpression mexpr;
    mexpr.setExpr( coordExpr, world, directoryLibExpr/*,indexes*/ );
    CHECK( mexpr.isScalar() || mexpr.isVector() ) << "wrong kind of expression with coord";

    this->setExpression( std::move( mexpr ) );
}

void
ModelPostprocessPointPosition::PointsOverCoordinates::setup( nl::json const& jarg, worldcomm_t const& world, std::string const& directoryLibExpr, ModelIndexes const& indexes )
{
    std::vector<std::string> coordExprs;
    if ( jarg.is_string() )
    {
        coordExprs.push_back( jarg.get<std::string>() );
    }
    else if ( jarg.is_array() )
    {
        for ( auto const& el : jarg.items() )
        {
            CHECK( el.value().is_string() ) << "coord item should be a string";
            coordExprs.push_back( el.value() );
        }
    }

    for ( std::string coordExprBase : coordExprs )
    {
        PointPosition ptPos;
        ptPos.setup( coordExprBase, world, directoryLibExpr, indexes );
        M_pointPositions.push_back( std::move( ptPos ) );
    }
}

void
ModelPostprocessPointPosition::PointsOverSegment::setup( nl::json const& jarg, worldcomm_t const& world, std::string const& directoryLibExpr, ModelIndexes const& indexes )
{
    CHECK( jarg.contains( "point1") ) << "PointsOverSegment json should be have entry : point1";
    M_point1.setup( jarg.at( "point1"), world, directoryLibExpr, indexes );
    CHECK( jarg.contains( "point2") ) << "PointsOverSegment json should be have entry : point2";
    M_point2.setup( jarg.at( "point2"), world, directoryLibExpr, indexes );

    if ( jarg.contains( "n_points" ) )
    {
        auto const& jNPoints = jarg.at( "n_points" );
        if ( jNPoints.is_string() )
            M_nPoints = std::stoi( jNPoints.get<std::string>() );
        else if ( jNPoints.is_number_integer() )
            M_nPoints = jNPoints.get<int>();
    }
}

void
ModelPostprocessPointPosition::PointsOverSegment::updatePointsSampling( coord_value_type const& pt1, coord_value_type const& pt2 )
{
    this->M_coordinates.resize( M_nPoints );

    coord_value_type vec12 = pt2-pt1;

    for (int k=0;k<M_nPoints;++k)
    {
        this->M_coordinates[k] = pt1 + (((double)k)/(M_nPoints-1))*vec12;
        //std::cout << "pt["<<k<<"] = "<<M_pointsSampling[k] << std::endl;
    }
}

void
ModelPostprocessPointPosition::MeasuresOutput::setup( nl::json const& jarg )
{
    if ( jarg.contains( "name" ) )
        M_name = jarg.at("name").get<std::string>();
    if ( jarg.contains( "type" ) )
        M_type = jarg.at("type").get<std::string>();
}

void
ModelPostprocessPointPosition::setup( nl::json const& jarg, std::string const& name, ModelIndexes const& indexes )
{
    bool hasCoord = jarg.contains( "coord" );
    bool hasMarker = jarg.contains( "marker" );

#if 0
    // coord or marker is necessary
    if ( !hasCoord && !hasMarker )
        return;
#endif
    M_name = name;
    if ( hasMarker )
    {
#if 0
        // TODO
        std::string marker = M_p.get<std::string>( "marker" );
        this->pointPosition().setMeshMarker( marker );
#endif
    }
    if ( hasCoord )
    {
        nl::json const& jCoord = jarg.at( "coord" );

        auto poc = std::make_shared<PointsOverCoordinates>();
        poc->setup( jCoord, this->worldComm(), M_directoryLibExpr, indexes );
        if ( true ) // TODO : check if poc is good
            M_pointsOverAllGeometry.push_back( std::move( poc ) );

    } // hasCoord

    bool hasOverGeo = jarg.contains( "over_geometry" );
    if ( hasOverGeo )
    {
        auto const& jOverGeo = jarg.at( "over_geometry" );
        for ( auto const& el : jOverGeo.items() )
        {
            //std::cout << "el.key() " << el.key() << std::endl;
            if ( el.key() == "segment" )
            {
                CHECK( el.value().is_object() ) << "segment item should be an object";
                auto pol = std::make_shared<PointsOverSegment>();
                pol->setup( el.value(), this->worldComm(), M_directoryLibExpr, indexes );
                // TODO check if pol is operational
                if ( true )
                    M_pointsOverAllGeometry.push_back( std::move( pol ) );
            }
        }
    }

    if ( jarg.contains( "fields" ) )
    {
        auto const& jFields = jarg.at( "fields" );

        std::vector<std::string> fieldList;
        if ( jFields.is_string() )
        {
            fieldList.push_back( jFields.get<std::string>() );
        }
        else if ( jFields.is_array() )
        {
            for ( auto const& el : jFields.items() )
            {
                CHECK( el.value().is_string() ) << "fields item should be a string";
                fieldList.push_back( el.value() );
            }
        }

        for ( std::string field : fieldList )
            this->addFields( indexes.replace( field ) );
    }

    if ( jarg.contains( "expressions" ) )
    {
        auto const& jExprs = jarg.at( "expressions" );

        for ( auto const& el : jExprs.items() )
        {
            std::string const& exprName = el.key();
            auto const& exprJson = el.value();
            CHECK( exprJson.is_string() ) << "expr should be a string";
            std::string const& exprStr = exprJson.get<std::string>();

            ModelExpression mexpr;
            mexpr.setExpr( exprStr, this->worldComm(), M_directoryLibExpr/*,indexes*/ );
            CHECK( mexpr.hasAtLeastOneExpr() ) << "expr not given correctly";
            //std::cout << "MYEXPR " << exprName << " : " << mexpr.exprToString() << std::endl;
            M_exprs.emplace( exprName, std::move( mexpr ) );
        }
    }
    if ( jarg.contains( "include_coordinates" ) )
    {
        auto const& j_include_coordinates = jarg.at("include_coordinates");
        if ( j_include_coordinates.is_boolean() )
            M_includeCoordinates = j_include_coordinates.get<bool>();
        else if ( j_include_coordinates.is_string() )
            M_includeCoordinates = boost::lexical_cast<bool>( j_include_coordinates.get<std::string>() );
    }

    if ( jarg.contains( "output" ) )
        M_measuresOutput.setup( jarg.at("output") );
    if ( M_measuresOutput.name().empty() && M_measuresOutput.type() != "values" )
        M_measuresOutput.setName( this->name() );
}

void
ModelPostprocessPointPosition::setParameterValues( std::map<std::string,double> const& mp )
{
    for ( auto & [exprName,mexpr] : M_exprs )
        mexpr.setParameterValues( mp );

    for ( auto & pointsOverGeometry : M_pointsOverAllGeometry )
    {
        if ( !pointsOverGeometry )
            continue;
        pointsOverGeometry->setParameterValues( mp );
    }
}

void
ModelPostprocessNorm::setup( nl::json const& jarg, std::string const& name, ModelIndexes const& indexes )
{
    M_name = name;

    if ( jarg.contains( "field" ) )
    {
        auto const& j_field = jarg.at("field");
        if ( j_field.is_string() )
            M_field = indexes.replace( j_field.get<std::string>() );
    }
    else if ( jarg.contains( "expr" ) )
    {
        M_expr.setExpr( jarg.at("expr"), this->worldComm(), M_directoryLibExpr, indexes );
        if ( jarg.contains( "grad_expr" ) )
            M_gradExpr.setExpr( jarg.at( "grad_expr" ), this->worldComm(), M_directoryLibExpr, indexes );
    }

    if ( jarg.contains( "markers" ) )
        M_markers.setup( jarg.at( "markers" ), indexes );

    if ( jarg.contains( "type" ) )
    {
        auto const& j_type = jarg.at( "type" );
        if ( j_type.is_string() )
            M_types.insert( indexes.replace( j_type.get<std::string>() ) );
        else if ( j_type.is_array() )
            for ( auto const& [j_typekey,j_typeval] : j_type.items() )
                if ( j_typeval.is_string() )
                    M_types.insert( indexes.replace( j_typeval.get<std::string>() ) );
    }

    if ( jarg.contains( "solution" ) )
        M_solution.setExpr( jarg.at( "solution" ), this->worldComm(), M_directoryLibExpr,indexes );

    if ( jarg.contains( "grad_solution" ) )
        M_gradSolution.setExpr( jarg.at( "grad_solution" ), this->worldComm(), M_directoryLibExpr,indexes );

    if ( jarg.contains( "quad") )
    {
        auto const& j_quad = jarg.at( "quad" );
        if ( j_quad.is_number_integer() )
            M_quadOrder = j_quad.get<int>();
    }

    if ( jarg.contains( "quad1") )
    {
        auto const& j_quad1 = jarg.at( "quad1" );
        if ( j_quad1.is_number_integer() )
            M_quad1Order = j_quad1.get<int>();
    }
}

void
ModelPostprocessNorm::setParameterValues( std::map<std::string,double> const& mp )
{
    M_expr.setParameterValues( mp );
    M_gradExpr.setParameterValues( mp );
    M_solution.setParameterValues( mp );
    M_gradSolution.setParameterValues( mp );
}

void
ModelPostprocessStatistics::setup( nl::json const& jarg, std::string const& name, ModelIndexes const& indexes )
{
    M_name = name;

    if ( jarg.contains( "field" ) )
    {
        auto const& j_field = jarg.at("field");
        if ( j_field.is_string() )
            M_field = indexes.replace( j_field.get<std::string>() );
    }
    else if ( jarg.contains( "expr" ) )
    {
        M_expr.setExpr( jarg.at("expr"), this->worldComm(), M_directoryLibExpr, indexes );
    }

    if ( jarg.contains( "markers" ) )
        M_markers.setup( jarg.at( "markers" ), indexes );

    if ( jarg.contains( "type" ) )
    {
        auto const& j_type = jarg.at( "type" );
        if ( j_type.is_string() )
            M_types.insert( indexes.replace( j_type.get<std::string>() ) );
        else if ( j_type.is_array() )
            for ( auto const& [j_typekey,j_typeval] : j_type.items() )
                if ( j_typeval.is_string() )
                    M_types.insert( indexes.replace( j_typeval.get<std::string>() ) );
    }

    if ( jarg.contains( "quad") )
    {
        auto const& j_quad = jarg.at( "quad" );
        if ( j_quad.is_number_integer() )
            M_quadOrder = j_quad.get<int>();
    }

    if ( jarg.contains( "quad1") )
    {
        auto const& j_quad1 = jarg.at( "quad1" );
        if ( j_quad1.is_number_integer() )
            M_quad1Order = j_quad1.get<int>();
    }
}

void
ModelPostprocessStatistics::setParameterValues( std::map<std::string,double> const& mp )
{
    M_expr.setParameterValues( mp );
}


void
ModelPostprocessCheckerMeasure::setup( nl::json const& jarg, std::string const& name, ModelIndexes const& indexes )
{
    M_name = name;

    CHECK( jarg.contains("value") ) << "no value entry";
    M_valueExpr.setExpr( jarg.at( "value" ), this->worldComm(), M_directoryLibExpr, indexes );
    CHECK( M_valueExpr.hasExprScalar() ) << "require value entry and the value should be a scalar expression";

    if ( jarg.contains("tolerance" ) )
    {
        M_toleranceExpr.setExpr( jarg.at( "tolerance" ), this->worldComm(), M_directoryLibExpr, indexes );
        CHECK( M_toleranceExpr.hasExprScalar() ) << "tolerance entry should be a scalar expression";
    }

}
// std::tuple<bool,double>
// ModelPostprocessCheckerMeasure::run( double val ) const
// {
//     M_value = M_valueExpr.exprScalar().evaluate()(0,0);
//     return this->runImpl( val );
// }
std::tuple<bool,double>
ModelPostprocessCheckerMeasure::runImpl( double val ) const
{
    if ( std::abs(val) < 1e-10 || std::abs(M_value) < 1e-10 )
    {
        double diff = std::abs(val-M_value);
        return std::make_tuple( std::abs(val-M_value) <= M_tolerance, diff );
    }
    else
    {
        double diffAbs = std::abs(val-M_value);
        double maxDiffRel = std::max(diffAbs/std::abs(val),diffAbs/std::abs(M_value));
        return std::make_tuple( diffAbs/std::abs(val) <= M_tolerance && diffAbs/std::abs(M_value) <= M_tolerance, maxDiffRel );
    }
}

void
ModelPostprocessCheckerMeasure::setParameterValues( std::map<std::string,double> const& mp )
{
    M_valueExpr.setParameterValues( mp );
    M_toleranceExpr.setParameterValues( mp );
}

ModelPostprocess::ModelPostprocess( worldcomm_ptr_t const& world )
    :
    super( world ),
    M_useModelName( true )
{}

ModelPostprocess::~ModelPostprocess()
{}

bool
ModelPostprocess::hasJsonProperties( std::string const& name ) const
{
    if ( !M_useModelName )
        return true;
    else
        return M_p.contains(name);
}

nl::json const&
ModelPostprocess::jsonProperties( std::string const& name ) const
{
    if ( !M_useModelName )
        return M_p;
    CHECK( M_p.contains(name) ) << "name " << name << "not in data structure";
    return M_p.at( name );
}

void
ModelPostprocess::setPTree( nl::json const& jarg )
{
    M_p = jarg;
    setup();
}


void
ModelPostprocess::setup()
{
    auto const& jarg = M_p;

    if ( jarg.contains("use-model-name") )
    {
        auto const& j_useModelName = jarg.at("use-model-name");
        if ( j_useModelName.is_boolean() )
            M_useModelName = j_useModelName.template get<bool>();
        else if ( j_useModelName.is_string() )
            M_useModelName = boost::lexical_cast<bool>( j_useModelName.template get<std::string>() );
    }

    if ( M_useModelName )
    {
        for (auto const& el : jarg.items())
        {
            if ( !el.value().is_object() )
                continue;
            this->setup( el.key(),el.value() );
        }
    }
    else
        this->setup( "", jarg );
}
void
ModelPostprocess::setup( std::string const& name, nl::json const& jarg )
{
    auto itFindExports = jarg.find( "Exports" );
    if ( itFindExports != jarg.end() )
    {
        ModelPostprocessExports ppexports( this->worldCommPtr() );
        ppexports.setDirectoryLibExpr( M_directoryLibExpr );
        ppexports.setup( itFindExports.value() );
        if ( !ppexports.fields().empty() || !ppexports.expressions().empty() )
            M_exports.emplace( std::make_pair(name, std::move( ppexports ) ) );
    }

    auto itFindSave = jarg.find( "Save" );
    if ( itFindSave != jarg.end() )
    {
        ModelPostprocessSave ppsave;
        ppsave.setup( itFindSave.value() );
        if ( !ppsave.fieldsNames().empty() )
            M_save.emplace( std::make_pair(name, std::move( ppsave ) ) );
    }

    auto itFindMeasures = jarg.find( "Measures" );
    if ( itFindMeasures != jarg.end() )
    {
        auto const& j_measures = itFindMeasures.value();
        for ( std::string const& quantitiesSectionName : { "Quantities", "quantities" } )
            if ( j_measures.contains( quantitiesSectionName ) )
            {
                ModelPostprocessQuantities ppquantities( this->worldCommPtr() );
                ppquantities.setDirectoryLibExpr( M_directoryLibExpr );
                ppquantities.setup( j_measures.at( quantitiesSectionName ) );
                if( !ppquantities.quantities().empty() || !ppquantities.expressions().empty() )
                    M_measuresQuantities.emplace( std::make_pair( name, std::move( ppquantities ) ) );
                break;
            }

        if ( j_measures.contains( "Points" ) )
        {
            auto const& j_measures_points = j_measures.at( "Points" );
            for ( auto const& [j_measures_pointskey,j_measures_pointsval] : j_measures_points.items() )
            {
                auto indexesAllCases = ModelIndexes::generateAllCases( j_measures_pointsval );
                for ( auto const& indexes : indexesAllCases )
                {
                    ModelPostprocessPointPosition myPpPtPos( this->worldCommPtr() );
                    myPpPtPos.setDirectoryLibExpr( M_directoryLibExpr );
                    myPpPtPos.setup( j_measures_pointsval, indexes.replace( j_measures_pointskey ), indexes );
                    if ( !myPpPtPos.fields().empty() || !myPpPtPos.expressions().empty() )
                        M_measuresPoint[name].push_back( std::move(myPpPtPos) );
                }
            }
        }


        for ( std::string const& normsSectionName : { "Norm","Norms" } )
        {
            if ( j_measures.contains( normsSectionName ) )
            {
                auto const& j_measures_norms = j_measures.at( normsSectionName );
                for ( auto const& [j_measures_normskey,j_measures_normsval] : j_measures_norms.items() )
                {
                    auto indexesAllCases = ModelIndexes::generateAllCases( j_measures_normsval );
                    for ( auto const& indexes : indexesAllCases )
                    {
                        ModelPostprocessNorm ppNorm( this->worldCommPtr() );
                        ppNorm.setDirectoryLibExpr( M_directoryLibExpr );
                        ppNorm.setup( j_measures_normsval, indexes.replace( j_measures_normskey ), indexes );
                        if ( ppNorm.hasField() || ppNorm.hasExpr() )
                            M_measuresNorm[name].push_back( std::move( ppNorm ) );
                    }
                }
            }
        }


        for ( std::string const& statisticsSectionName : { "Statistic","Statistics" } )
        {
            if ( j_measures.contains( statisticsSectionName ) )
            {
                auto const& j_measures_statistics = j_measures.at( statisticsSectionName );
                for ( auto const& [j_measures_statisticskey,j_measures_statisticsval] : j_measures_statistics.items() )
                {
                    auto indexesAllCases = ModelIndexes::generateAllCases( j_measures_statisticsval );
                    for ( auto const& indexes : indexesAllCases )
                    {
                        ModelPostprocessStatistics ppStatistics( this->worldCommPtr() );
                        ppStatistics.setDirectoryLibExpr( M_directoryLibExpr );
                        ppStatistics.setup( j_measures_statisticsval, indexes.replace( j_measures_statisticskey ), indexes );
                        if ( ppStatistics.hasField() || ppStatistics.hasExpr() )
                            M_measuresStatistics[name].push_back( std::move( ppStatistics ) );
                    }
                }
            }
        }
    }


    auto itFindCheckers = jarg.find( "Checkers" );
    if ( itFindCheckers != jarg.end() )
    {
        auto const& j_checkers = itFindCheckers.value();
        if ( j_checkers.contains( "Measures" ) )
        {
            auto const& j_checkers_measures = j_checkers.at( "Measures" );
            for ( auto const& [j_checkers_measureskey,j_checkers_measuresval] : j_checkers_measures.items() )
            {
                auto indexesAllCases = ModelIndexes::generateAllCases( j_checkers_measuresval );
                for ( auto const& indexes : indexesAllCases )
                {
                    ModelPostprocessCheckerMeasure ppCheckerMeasure( this->worldCommPtr() );
                    ppCheckerMeasure.setDirectoryLibExpr( M_directoryLibExpr );
                    ppCheckerMeasure.setup( j_checkers_measuresval, indexes.replace( j_checkers_measureskey ), indexes );
                    M_checkersMeasure[name].push_back( std::move(ppCheckerMeasure) );
                }
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

    for (auto & [name,p] : M_exports )
        p.setParameterValues( mp );

    for( auto & p : M_measuresPoint )
        for( auto & p2 : p.second )
            p2.setParameterValues( mp );

    for( auto & p : M_measuresNorm )
        for( auto & p2 : p.second )
            p2.setParameterValues( mp );

    for( auto & p : M_measuresStatistics )
        for( auto & p2 : p.second )
            p2.setParameterValues( mp );

    for (auto & [name,p] : M_measuresQuantities )
        p.setParameterValues( mp );

    for( auto & p : M_checkersMeasure )
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
ModelPostprocess::hasSave( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return M_save.find( nameUsed ) != M_save.end();
}
bool
ModelPostprocess::hasMeasuresQuantities( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return M_measuresQuantities.find( nameUsed ) != M_measuresQuantities.end();
}
bool
ModelPostprocess::hasMeasuresPoint( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return M_measuresPoint.find( nameUsed ) != M_measuresPoint.end();
}
bool
ModelPostprocess::hasMeasuresNorm( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return M_measuresNorm.find( nameUsed ) != M_measuresNorm.end();
}
bool
ModelPostprocess::hasMeasuresStatistics( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return M_measuresStatistics.find( nameUsed ) != M_measuresStatistics.end();
}
bool
ModelPostprocess::hasCheckersMeasure( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return M_checkersMeasure.find( nameUsed ) != M_checkersMeasure.end();
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
ModelPostprocessSave const&
ModelPostprocess::save( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    //CHECK( this->hasSave( nameUsed ) ) << "no save with name:"<<name;
    if ( this->hasSave( nameUsed ) )
        return M_save.find( nameUsed )->second;
    else
        return M_emptySave;
}
ModelPostprocessQuantities const&
ModelPostprocess::measuresQuantities( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    //CHECK( this->hasExports( nameUsed ) ) << "no measures quantities with name:"<<name;
    if ( this->hasMeasuresQuantities( nameUsed ) )
        return M_measuresQuantities.find( nameUsed )->second;
    else
        return M_emptyMeasuresQuantities;
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
std::vector<ModelPostprocessPointPosition> &
ModelPostprocess::measuresPoint( std::string const& name )
{
    std::string nameUsed = (M_useModelName)? name : "";
    //CHECK( this->hasMeasuresPoint( nameUsed ) ) << "no measures point with name:"<<name;
    if ( this->hasMeasuresPoint( nameUsed ) )
        return M_measuresPoint.find( nameUsed )->second;
    else
        return M_emptyMeasuresPoint;
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
std::vector<ModelPostprocessStatistics> const&
ModelPostprocess::measuresStatistics( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    //CHECK( this->hasMeasuresNorm( nameUsed ) ) << "no measures norm with name:"<<name;
    if ( this->hasMeasuresStatistics( nameUsed ) )
        return M_measuresStatistics.find( nameUsed )->second;
    else
        return M_emptyMeasuresStatistics;
}
std::vector<ModelPostprocessCheckerMeasure> const&
ModelPostprocess::checkersMeasure( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    if ( this->hasCheckersMeasure( nameUsed ) )
        return M_checkersMeasure.find( nameUsed )->second;
    else
        return M_emptyCheckersMeasure;
}



}
