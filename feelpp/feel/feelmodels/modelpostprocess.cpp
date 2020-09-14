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
    if ( auto fields = p.get_child_optional("fields") )
    {
        if ( fields->empty() ) // value case
            M_fields.insert( fields->get_value<std::string>() );
        else // array case
        {
            for ( auto const& item : *fields )
            {
                CHECK( item.first.empty() ) << "should be an array, not a subtree";
                std::string const& fieldName = item.second.template get_value<std::string>();
                M_fields.insert( fieldName );
                LOG(INFO) << "add to postprocess field  " << fieldName;
            }
        }
    }
    if ( auto formatOpt = p.get_optional<std::string>( "format" ) )
        M_format = *formatOpt;

    if ( auto exprTree = p.get_child_optional("expr") )
    {
        for ( auto const& item : *exprTree )
        {
            std::string exprName = item.first;
            std::set<std::string> representations, tags;

            if ( item.second.empty() ) // name:expr
            {
                ModelExpression modelexpr;
                ModelMarkers markers;
                modelexpr.setExpr( exprName,  *exprTree, this->worldComm(), M_directoryLibExpr/*,indexes*/ );
                CHECK( modelexpr.hasAtLeastOneExpr() ) << "expr not given correctly";
                M_exprs.push_back( std::make_tuple(exprName,modelexpr,markers,representations,tags) );
            }
            else
            {
                if ( auto repOpt = item.second.get_child_optional("representation") )
                {
                    if ( repOpt->empty() )
                        representations.insert( repOpt->get_value<std::string>() );
                    else
                    {
                        for ( auto const& therep : *repOpt )
                        {
                            CHECK( therep.first.empty() ) << "should be an array, not a subtree";
                            representations.insert( therep.second.template get_value<std::string>() );
                        }
                    }
                }

                if ( auto tagOpt = item.second.get_child_optional("tag") )
                {
                    if ( tagOpt->empty() )
                        tags.insert( tagOpt->get_value<std::string>() );
                    else
                    {
                        for ( auto const& thetag : *tagOpt )
                            tags.insert( thetag.second.template get_value<std::string>() );
                    }
                }

                if ( item.second.get_child_optional("expr") )
                {
                    ModelExpression modelexpr;
                    ModelMarkers markers;
                    modelexpr.setExpr( "expr",  item.second, this->worldComm(), M_directoryLibExpr/*,indexes*/ );
                    if ( auto ptmarkers = item.second.get_child_optional("markers") )
                        markers.setPTree(*ptmarkers/*, indexes*/);
                    CHECK( modelexpr.hasAtLeastOneExpr() ) << "expr not given correctly";
                    M_exprs.push_back( std::make_tuple(exprName,modelexpr,markers,representations,tags) );
                }

                if ( auto ptparts = item.second.get_child_optional("parts") )
                {
                    for ( auto const& itemPart : *ptparts )
                    {
                        CHECK( itemPart.first.empty() ) << "should be an array, not a subtree";
                        ModelExpression modelexpr;
                        ModelMarkers markers;
                        modelexpr.setExpr( "expr",  itemPart.second, this->worldComm(), M_directoryLibExpr/*,indexes*/ );
                        if ( auto ptmarkers = itemPart.second.get_child_optional("markers") )
                            markers.setPTree(*ptmarkers/*, indexes*/);
                        CHECK( modelexpr.hasAtLeastOneExpr() ) << "expr not given correctly";
                        M_exprs.push_back( std::make_tuple(exprName,modelexpr,markers,representations,tags) );
                    }
                }
            }

        } // for ( auto const& item : *exprTree )
    }
}

void
ModelPostprocessExports::setParameterValues( std::map<std::string,double> const& mp )
{
    for ( auto & edata : M_exprs )
        std::get<1>( edata ).setParameterValues( mp );
}

void
ModelPostprocessQuantities::setup( pt::ptree const& p )
{
    if ( auto fields = p.get_child_optional("quantities") )
    {
        if ( fields->empty() ) // value case
            M_quantities.insert( fields->get_value<std::string>() );
        else // array case
        {
            for ( auto const& item : *fields )
            {
                CHECK( item.first.empty() ) << "should be an array, not a subtree";
                std::string const& fieldName = item.second.template get_value<std::string>();
                M_quantities.insert( fieldName );
                LOG(INFO) << "add to postprocess quantity  " << fieldName;
            }
        }
    }
}

void
ModelPostprocessSave::setup( pt::ptree const& p )
{
    if ( auto fieldsPtree = p.get_child_optional("Fields") )
    {
        if ( auto fieldsNamesPtree = fieldsPtree->get_child_optional("names") )
        {
            if ( fieldsNamesPtree->empty() ) // value case
                M_fieldsNames.insert( fieldsNamesPtree->get_value<std::string>() );
            else // array case
            {
                for ( auto const& item : *fieldsNamesPtree )
                {
                    CHECK( item.first.empty() ) << "should be an array, not a subtree";
                    std::string const& fieldName = item.second.template get_value<std::string>();
                    M_fieldsNames.insert( fieldName );
                }
            }
        }
        if ( auto formatOpt = fieldsPtree->get_optional<std::string>( "format" ) )
            M_fieldsFormat = *formatOpt;
    }
}


void
ModelPostprocessPointPosition::setup( std::string const& name, ModelIndexes const& indexes )
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
        std::string coordExpr = indexes.replace( M_p.get<std::string>( "coord" ) );
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
            this->pointPosition().setExpression( coordExpr, M_directoryLibExpr, this->worldComm() );
        }
    } // hasCoord

    // store fields name
    std::vector<std::string> fieldList = as_vector<std::string>( M_p, "fields" );
    if ( fieldList.empty() )
    {
        std::string fieldUnique = indexes.replace( M_p.get<std::string>( "fields" ) );
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
ModelPostprocessNorm::setup( std::string const& name, ModelIndexes const& indexes )
{
    M_name = name;

    if ( auto itField = M_p.get_optional<std::string>("field") )
        M_field = indexes.replace( *itField );
    else if ( auto ptexpr = M_p.get_child_optional("expr") )
    {
        M_expr.setExpr( "expr", M_p, this->worldComm(), M_directoryLibExpr, indexes );
        if ( auto ptgradexpr = M_p.get_child_optional("grad_expr") )
            M_gradExpr.setExpr( "grad_expr", M_p, this->worldComm(), M_directoryLibExpr, indexes );
    }

    if ( auto ptmarkers = M_p.get_child_optional("markers") )
        M_markers.setPTree(*ptmarkers,indexes);

    if ( auto pttype = M_p.get_child_optional("type") )
    {
        for( auto const& item : M_p.get_child("type") )
            M_types.insert( indexes.replace( item.second.template get_value<std::string>() ) );
        if( M_types.empty() )
            M_types.insert( indexes.replace( M_p.get<std::string>("type") ) );
    }

    if ( auto itSol = M_p.get_optional<std::string>("solution") )
        M_solution.setExpr( "solution", M_p, this->worldComm(), M_directoryLibExpr,indexes );

    if ( auto itSol = M_p.get_optional<std::string>("grad_solution") )
        M_gradSolution.setExpr( "grad_solution", M_p, this->worldComm(), M_directoryLibExpr,indexes );

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
    M_expr.setParameterValues( mp );
    M_gradExpr.setParameterValues( mp );
    M_solution.setParameterValues( mp );
    M_gradSolution.setParameterValues( mp );
}

void
ModelPostprocessStatistics::setup( std::string const& name, ModelIndexes const& indexes )
{
    M_name = name;

    if ( auto itField = M_p.get_optional<std::string>("field") )
        M_field = indexes.replace( *itField );
    else if ( auto ptexpr = M_p.get_child_optional("expr") )
    {
        M_expr.setExpr( "expr", M_p, this->worldComm(), M_directoryLibExpr, indexes );
    }

    if ( auto ptmarkers = M_p.get_child_optional("markers") )
        M_markers.setPTree(*ptmarkers, indexes);

    if ( auto pttype = M_p.get_child_optional("type") )
    {
        for( auto const& item : M_p.get_child("type") )
            M_types.insert( indexes.replace( item.second.template get_value<std::string>() ) );
        if( M_types.empty() )
            M_types.insert( indexes.replace( M_p.get<std::string>("type") ) );
    }

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
ModelPostprocessStatistics::setParameterValues( std::map<std::string,double> const& mp )
{
    M_expr.setParameterValues( mp );
}


void
ModelPostprocessCheckerMeasure::setup( std::string const& name, ModelIndexes const& indexes )
{
    M_name = name;

    M_valueExpr.setExpr( "value", M_p, this->worldComm(), M_directoryLibExpr, indexes );
    CHECK( M_valueExpr.hasExprScalar() ) << "require value entry and the value should be a scalar expression";
    M_value = M_valueExpr.exprScalar().evaluate()(0,0);

    if ( auto itTol = M_p.get_optional<double>("tolerance") )
        M_tolerance = *itTol;
}
std::tuple<bool,double>
ModelPostprocessCheckerMeasure::run( double val ) const
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
    M_value = M_valueExpr.exprScalar().evaluate()(0,0);
}

ModelPostprocess::ModelPostprocess( worldcomm_ptr_t const& world )
    :
    super( world ),
    M_useModelName( false )
{}

ModelPostprocess::ModelPostprocess(pt::ptree const& p, worldcomm_ptr_t const& world )
    :
    super( world ),
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
        ModelPostprocessExports ppexports( this->worldCommPtr() );
        ppexports.setDirectoryLibExpr( M_directoryLibExpr );
        ppexports.setup( *exports );
        if ( !ppexports.fields().empty() || !ppexports.expressions().empty() )
            M_exports[name] = ppexports;
    }
    if ( auto save = p.get_child_optional("Save") )
    {
        ModelPostprocessSave ppsave;
        ppsave.setup( *save );
        if ( !ppsave.fieldsNames().empty() )
            M_save[name] = ppsave;
    }

    if ( auto measures = p.get_child_optional("Measures") )
    {
        ModelPostprocessQuantities ppquantities;
        ppquantities.setup( *measures );
        if( !ppquantities.quantities().empty() )
            M_measuresQuantities[name] = ppquantities;

        auto evalPoints = measures->get_child_optional("Points");
        if ( evalPoints )
        {
            for( auto const& evalPoint : *evalPoints )
            {
                auto indexesAllCases = ModelIndexes::generateAllCases( evalPoint.second );
                for ( auto const& indexes : indexesAllCases )
                {
                    ModelPostprocessPointPosition myPpPtPos( this->worldCommPtr() );
                    myPpPtPos.setDirectoryLibExpr( M_directoryLibExpr );
                    myPpPtPos.setPTree( evalPoint.second, indexes.replace( evalPoint.first ), indexes );
                    if ( !myPpPtPos.fields().empty() )
                        M_measuresPoint[name].push_back( myPpPtPos );
                }
            }
        }

        auto ptreeNorms = measures->get_child_optional("Norm");
        if ( ptreeNorms )
        {
            for( auto const& ptreeNorm : *ptreeNorms )
            {
                auto indexesAllCases = ModelIndexes::generateAllCases(  ptreeNorm.second );
                for ( auto const& indexes : indexesAllCases )
                {
                    ModelPostprocessNorm ppNorm( this->worldCommPtr() );
                    ppNorm.setDirectoryLibExpr( M_directoryLibExpr );
                    ppNorm.setPTree( ptreeNorm.second, indexes.replace( ptreeNorm.first ), indexes );
                    if ( ppNorm.hasField() || ppNorm.hasExpr() )
                        M_measuresNorm[name].push_back( ppNorm );
                }
            }
        }
        auto ptreeStatistics = measures->get_child_optional("Statistics");
        if ( ptreeStatistics )
        {
            for( auto const& ptreeStatistic : *ptreeStatistics )
            {
                auto indexesAllCases = ModelIndexes::generateAllCases( ptreeStatistic.second );
                for ( auto const& indexes : indexesAllCases )
                {
                    ModelPostprocessStatistics ppStatistics( this->worldCommPtr() );
                    ppStatistics.setDirectoryLibExpr( M_directoryLibExpr );
                    ppStatistics.setPTree( ptreeStatistic.second, indexes.replace( ptreeStatistic.first ), indexes );
                    if ( ppStatistics.hasField() || ppStatistics.hasExpr() )
                        M_measuresStatistics[name].push_back( ppStatistics );
                }
            }
        }
    }

    if ( auto checkers = p.get_child_optional("Checkers") )
    {
        if ( auto measures = checkers->get_child_optional("Measures") )
        {
            for( auto const& ptreeCheckerMeasure : *measures )
            {
                auto indexesAllCases = ModelIndexes::generateAllCases( ptreeCheckerMeasure.second );
                for ( auto const& indexes : indexesAllCases )
                {
                    ModelPostprocessCheckerMeasure ppCheckerMeasure( this->worldCommPtr() );
                    ppCheckerMeasure.setDirectoryLibExpr( M_directoryLibExpr );
                    ppCheckerMeasure.setPTree( ptreeCheckerMeasure.second, indexes.replace( ptreeCheckerMeasure.first ), indexes );
                    M_checkersMeasure[name].push_back( ppCheckerMeasure );
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
