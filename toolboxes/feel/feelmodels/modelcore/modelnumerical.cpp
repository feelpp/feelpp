/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2012-01-19

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3.0 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
/**
 \file modelnumerical.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelcore/remotedata.hpp>

namespace Feel
{
namespace FeelModels
{


ModelNumerical::ModelNumerical( std::string const& _theprefix, std::string const& keyword,
                                worldcomm_ptr_t const& _worldComm, std::string const& subPrefix,
                                ModelBaseRepository const& modelRep,
                                ModelBaseCommandLineOptions const& modelCmdLineOpt )
        :
        ModelBase( _theprefix, keyword, _worldComm, subPrefix, modelRep, modelCmdLineOpt ),
        super_type( _theprefix, keyword, _worldComm, subPrefix, modelRep, modelCmdLineOpt ),
        M_isStationary( boption(_name="ts.steady") ),
        M_doRestart( boption(_name="ts.restart") ),
        M_restartPath( soption(_name="ts.restart.path") ),
        M_restartAtLastSave( boption(_name="ts.restart.at-last-save") ),
        M_timeInitial( doption(_name="ts.time-initial") ),
        M_timeFinal( doption(_name="ts.time-final") ),
        M_timeStep( doption(_name="ts.time-step") ),
        M_timeOrder( ioption(_prefix=_theprefix, _name="ts.order",_vm=this->clovm()) ),
        M_tsSaveInFile( boption(_name="ts.save") ),
        M_tsSaveFreq( ioption(_name="ts.save.freq") ),
        M_timeCurrent(M_timeInitial),
        M_manageParameterValues( true ),
        M_manageParameterValuesOfModelProperties( true ),
        M_exporterPath( (fs::path(this->rootRepository())/prefixvm(this->prefix(), prefixvm(this->subPrefix(),"exports"))).string() ),
        M_postProcessSaveRepository( fs::path(this->rootRepository())/prefixvm(this->prefix(), prefixvm(this->subPrefix(),"save")) ),
        M_postProcessMeasures( (fs::path(this->rootRepository())/prefixvm(this->prefix(), prefixvm(this->subPrefix(),"measures"))).string(), this->worldCommPtr() ),
        //M_PsLogger( new PsLogger(prefixvm(this->prefix(),"PsLogger"),this->worldComm() ) )
        M_useChecker( boption(_prefix=this->prefix(),_name="checker",_vm=this->clovm()) )
    {
        //-----------------------------------------------------------------------//
        // move in stationary mode if we have this relation
        if ( M_timeInitial + M_timeStep == M_timeFinal)
            M_isStationary=true;
        //-----------------------------------------------------------------------//
        if (soption(_prefix=this->prefix(),_name="geomap",_vm=this->clovm())=="opt")
            M_geomap=GeomapStrategyType::GEOMAP_OPT;
        else
            M_geomap=GeomapStrategyType::GEOMAP_HO;
        //-----------------------------------------------------------------------//
        std::string modelPropFilename = Environment::expand( soption( _name="filename",_prefix=this->prefix(),_vm=this->clovm()) );
        if ( !modelPropFilename.empty() )
        {
            this->setModelProperties( modelPropFilename );
            if ( auto ptMeshes = M_modelProps->pTree().get_child_optional("Meshes"))
                super_model_meshes_type::setup( *ptMeshes );
        }
    }

   void
   ModelNumerical::setModelProperties( std::string const& filename )
   {
        M_modelProps = std::make_shared<ModelProperties>( filename, this->repository().expr(),
                                                          this->worldCommPtr(), this->prefix() );
   }

   void
   ModelNumerical::setModelProperties( nl::json const& json )
   {
        M_modelProps = std::make_shared<ModelProperties>( json, this->repository().expr(),
                                                          this->worldCommPtr(), this->prefix() );
   }

   void
   ModelNumerical::addParameterInModelProperties( std::string const& symbolName,double value)
   {
        M_modelProps->parameters()[symbolName] = ModelParameter(symbolName,value);
   }

    void
    ModelNumerical::setStationary(bool b)
    {
        if ( M_isStationary != b)
        {
            M_isStationary=b;
            this->setNeedToRebuildCstPart(true);
        }
    }

    void
    ModelNumerical::updateTime(double t)
    {
        M_timeCurrent=t;
        if ( M_modelProps )
            this->addParameterInModelProperties( "t", M_timeCurrent );
        super_model_meshes_type::updateTime( t );
    }

void
ModelNumerical::initPostProcess()
{
    if ( this->hasModelProperties() )
    {
        for ( auto & [tag,fields] : M_postProcessExportsFields )
            std::get<0>( fields ) = this->postProcessExportsFields( tag, this->modelProperties().postProcess().exports( this->keyword() ).fields() );

        M_postProcessSaveFields = this->postProcessSaveFields( this->modelProperties().postProcess().save( this->keyword() ).fieldsNames() );
        M_postProcessSaveFieldsFormat = this->modelProperties().postProcess().save( this->keyword() ).fieldsFormat();
        M_postProcessMeasuresQuantitiesNames = this->postProcessMeasuresQuantitiesNames( this->modelProperties().postProcess().measuresQuantities( this->keyword() ).quantities() );
    }
    if ( M_postProcessSaveFieldsFormat.empty() )
        M_postProcessSaveFieldsFormat = "default";
}

bool
ModelNumerical::hasPostProcessExportsExpr( std::string const& exportTag ) const
{
    if ( this->hasModelProperties() )
    {
        for ( auto const& [_name,_expr,_markers,_rep,_tag] : this->modelProperties().postProcess().exports( this->keyword() ).expressions() )
        {
            if ( !_tag.empty() && _tag.find( exportTag ) == _tag.end() )
                continue;
            if ( _tag.empty() && exportTag != "" )
                continue;
            return true;
        }
    }
    return false;
}

std::set<std::string>
ModelNumerical::postProcessExportsFields( std::string const& tag, std::set<std::string> const& ifields, std::string const& prefix ) const
{
    std::set<std::string> res;
    auto itFindTag = M_postProcessExportsFields.find( tag );
    if ( itFindTag == M_postProcessExportsFields.end() )
        return res;
    for ( auto const& o : ifields )
    {
        for ( auto const& fieldAvailable : std::get<1>( itFindTag->second ) )
            if ( o == prefixvm(prefix,fieldAvailable) || o == prefixvm(prefix,"all") || o == "all" )
                res.insert( fieldAvailable );
    }
    std::string const& pidName = std::get<2>( itFindTag->second );
    if ( !pidName.empty() )
    {
        if ( ifields.find( pidName ) != ifields.end() || ifields.find( "all" ) != ifields.end() )
            res.insert( pidName );
    }
    return res;
}
std::set<std::string>
ModelNumerical::postProcessSaveFields( std::set<std::string> const& ifields, std::string const& prefix ) const
{
    std::set<std::string> res;
    for ( auto const& o : ifields )
    {
        for ( auto const& fieldAvailable : M_postProcessSaveAllFieldsAvailable )
            if ( o == prefixvm(prefix,fieldAvailable) || o == prefixvm(prefix,"all") || o == "all" )
                res.insert( fieldAvailable );
    }
    return res;
}


std::set<std::string>
ModelNumerical::postProcessMeasuresQuantitiesNames( std::set<std::string> const& inames ) const
{
    std::set<std::string> res;
    for ( auto const& o : inames )
    {
        for ( auto const& nameAvailable : M_postProcessMeasuresQuantitiesAllNamesAvailable )
            if ( o == nameAvailable || o == "all" )
                res.insert( nameAvailable );
    }
    return res;
}

} // namespace FeelModels

} // namespace Feel

