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

namespace Feel
{
namespace FeelModels
{

ModelMeasuresIO::ModelMeasuresIO( std::string const& pathFile, WorldComm const& worldComm )
    :
    M_worldComm( worldComm ),
    M_pathFile( pathFile )
{}

void
ModelMeasuresIO::clear()
{
    M_mapParameterData.clear();
    M_mapMeasureData.clear();
}

void
ModelMeasuresIO::start()
{
    if ( M_worldComm.isMasterRank() && !M_mapMeasureData.empty() )
    {
        bool hasAlreadyWrited = false;
        std::ofstream fileWrited(M_pathFile, std::ios::out | std::ios::trunc);
        for ( auto const& data : M_mapParameterData )
        {
            int spacing = std::max(20, int(data.first.size()+2) );
            if ( hasAlreadyWrited )
                fileWrited << ",";
            fileWrited << std::setw( spacing ) << std::left << data.first;
            hasAlreadyWrited = true;
        }
        for ( auto const& data : M_mapMeasureData )
        {
            int spacing = std::max(28, int(data.first.size()+2) );
            if ( hasAlreadyWrited )
                fileWrited << ",";
            fileWrited << std::setw( spacing ) << std::right << data.first;
            hasAlreadyWrited = true;
        }
        fileWrited << std::endl;
        fileWrited.close();
    }
}
void
ModelMeasuresIO::restart( std::string const& paramKey, double val )
{
    if ( M_worldComm.isMasterRank() && !M_mapMeasureData.empty() )
    {
        double ti = val;//this->timeInitial();
        std::ifstream fileI(M_pathFile, std::ios::in);
        double timeLoaded=0;
        double valueLoaded = 0.;

        bool find=false;
        std::ostringstream buffer;
        bool hasAlreadyWrited = false;

            std::string measureTag;
            for ( auto const& data : M_mapParameterData )
            {
                fileI >> measureTag; // e.g. load time
                if ( measureTag == "," )
                    fileI >> measureTag; // e.g. load time
                measureTag.erase( std::remove(measureTag.begin(), measureTag.end(), ','), measureTag.end() );

                int spacing = std::max(20, int(data.first.size()+2) );
                if ( hasAlreadyWrited )
                    buffer << ",";
                buffer << std::setw(spacing) << std::left << measureTag;
                hasAlreadyWrited=true;
            }
            for ( auto const& data : M_mapMeasureData )
            {
                fileI >> measureTag;
                if ( measureTag == "," )
                    fileI >> measureTag; // e.g. load time
                measureTag.erase( std::remove(measureTag.begin(), measureTag.end(), ','), measureTag.end() );

                int spacing = std::max(28, int(data.first.size()+2) );
                if ( hasAlreadyWrited )
                    buffer << ",";
                buffer << std::setw(spacing) << std::right << measureTag;
                hasAlreadyWrited=true;
            }
            buffer << std::endl;

            while ( !fileI.eof() && !find )
            {
                hasAlreadyWrited = false;

#if 0
                for ( auto const& data : M_mapParameterData )
                    if ( paramKey == M_mapParameterData.first )
                        if ( (timeLoaded-1e-7) > ti ) { find=true;break; }
                if ( find ) break;
#endif
                for ( auto const& data : M_mapParameterData )
                {
#if 1
                    fileI >> measureTag;
                    if ( measureTag == "," )
                        fileI >> measureTag; // e.g. load time
                    measureTag.erase( std::remove(measureTag.begin(), measureTag.end(), ','), measureTag.end() );
                    valueLoaded = std::stod(measureTag);
#else
                    fileI >> valueLoaded;
#endif
                    //std::cout << "timeLoaded " << timeLoaded << " ti " << ti << "\n";
                    int spacing = std::max(20, int(data.first.size()+2) );
                    if ( hasAlreadyWrited )
                        buffer << ",";
                    buffer << std::setw(spacing) << std::left << std::setprecision( 9 ) << std::scientific << valueLoaded;
                    hasAlreadyWrited = true;
                    // check if last writing (e.g. time equality)
                    if ( paramKey == data.first && !find )
                        if ( std::abs(valueLoaded-ti) < 1e-9 )
                            find=true;
                }
                for ( auto const& data : M_mapMeasureData )
                {
#if 1
                    fileI >> measureTag;
                    if ( measureTag == "," )
                        fileI >> measureTag; // e.g. load time
                    measureTag.erase( std::remove(measureTag.begin(), measureTag.end(), ','), measureTag.end() );
                    valueLoaded = std::stod(measureTag);
#else
                    fileI >> valueLoaded;
#endif

                    int spacing = std::max(28, int(data.first.size()+2) );
                    if ( hasAlreadyWrited )
                        buffer << ",";
                    buffer << std::setw(spacing) << std::right << std::setprecision( 16 ) << std::scientific << valueLoaded;
                    hasAlreadyWrited = true;

                }
                buffer << std::endl;
            }
            fileI.close();
            std::ofstream fileW(M_pathFile/*.c_str()*/, std::ios::out | std::ios::trunc);
            fileW << buffer.str();
            fileW.close();
    }

}

void
ModelMeasuresIO::exportMeasures()
{
    if ( M_worldComm.isMasterRank() && !M_mapMeasureData.empty() )
    {
        bool hasAlreadyWrited = false;
        std::ofstream fileWrited(M_pathFile, std::ios::out | std::ios::app);
        for ( auto const& data : M_mapParameterData )
        {
            int spacing = std::max(20, int(data.first.size()+2) );
            if ( hasAlreadyWrited )
                fileWrited << ",";
            fileWrited << std::setw(spacing) << std::left << std::setprecision( 9 ) << std::scientific << data.second;
            hasAlreadyWrited = true;
        }
        for ( auto const& data : M_mapMeasureData )
        {
            int spacing = std::max(28, int(data.first.size()+2) );
            if ( hasAlreadyWrited )
                fileWrited << ",";
            fileWrited << std::setw(spacing) << std::right << std::setprecision( 16 ) << std::scientific << data.second;
            hasAlreadyWrited = true;
        }
        fileWrited << std::endl;
        fileWrited.close();
    }
}
void
ModelMeasuresIO::setParameter(std::string const& key,double val)
{
    M_mapParameterData[key] = val;
}
void
ModelMeasuresIO::setMeasure(std::string const& key,double val)
{
    M_mapMeasureData[key] = val;

}
void
ModelMeasuresIO::setMeasureComp( std::string const& key,std::vector<double> const& values )
{
    if ( values.empty() )
        return;
    int nValue = values.size();
    if ( nValue == 1 )
    {
        this->setMeasure( key,values[0]);
        return;
    }
    if ( nValue > 3 )
        return;

    this->setMeasure( key+"_x",values[0] );
    if ( nValue > 1 )
        this->setMeasure( key+"_y",values[1] );
    if ( nValue > 2 )
        this->setMeasure( key+"_z",values[2] );
}


void
ModelMeasuresEvaluatorContext::add( std::string const& field, int ctxId, std::string const& name )
{
    M_mapFieldToMapCtxIdToName[field][ctxId] = name;
}
bool
ModelMeasuresEvaluatorContext::has( std::string const& field ) const
{
    auto itFindField = M_mapFieldToMapCtxIdToName.find(field);
    if ( itFindField == M_mapFieldToMapCtxIdToName.end() )
        return false;
    return true;
}
bool
ModelMeasuresEvaluatorContext::has( std::string const& field, int ctxId ) const
{
    auto itFindField = M_mapFieldToMapCtxIdToName.find(field);
    if ( itFindField == M_mapFieldToMapCtxIdToName.end() )
        return false;
    auto itFindCtxId = itFindField->second.find(ctxId);
    if ( itFindCtxId == itFindField->second.end() )
        return false;
    return true;
}
std::string const&
ModelMeasuresEvaluatorContext::name( std::string const& field, int ctxId ) const
{
    if ( !this->has( field,ctxId ) )
         return M_emptyString;
    else
        return M_mapFieldToMapCtxIdToName.find(field)->second.find(ctxId)->second;
}
int
ModelMeasuresEvaluatorContext::ctxId( std::string const& field, std::string const& name ) const
{
    int ctxIdNull = -1;
    auto itFindField = M_mapFieldToMapCtxIdToName.find(field);
    if ( itFindField == M_mapFieldToMapCtxIdToName.end() )
        return ctxIdNull;

    for ( auto const& ctxIdAndName : itFindField->second )
    {
        int curCtxId = ctxIdAndName.first;
        std::string const& curName = ctxIdAndName.second;
        if ( name == curName )
            return curCtxId;
    }
    return ctxIdNull;
}


ModelNumerical::ModelNumerical( std::string _theprefix, WorldComm const& _worldComm, std::string subPrefix,
                                std::string appliShortRepository )
        :
        super_type( _theprefix, _worldComm, subPrefix, appliShortRepository ),
        M_rebuildMeshPartitions( boption(_name="rebuild_mesh_partitions",_prefix=this->prefix()) ),
        M_isStationary( false /*_isStationary*/),
        M_doRestart( boption(_name="ts.restart") ),
        M_restartPath( soption(_name="ts.restart.path") ),
        M_restartAtLastSave( boption(_name="ts.restart.at-last-save") ),
        M_timeInitial( doption(_name="ts.time-initial") ),
        M_timeFinal( doption(_name="ts.time-final") ),
        M_timeStep( doption(_name="ts.time-step") ),
        M_tsSaveInFile( boption(_name="ts.save") ),
        M_tsSaveFreq( ioption(_name="ts.save.freq") ),
        M_timeCurrent(M_timeInitial),
        M_directoryLibSymbExpr( "" ),
        M_row_startInMatrix(0),
        M_col_startInMatrix(0),
        M_row_startInVector(0),
        M_mshFileStr("FEELMODELS_WARNING_NODEFINE"),
        M_geoFileStr("FEELMODELS_WARNING_NODEFINE"),
        M_exporterPath( this->appliRepository()+"/"+prefixvm(this->prefix(), prefixvm(this->subPrefix(),"exports")) ),
        M_postProcessMeasuresIO( this->appliRepository()+"/"+prefixvm(this->prefix(), prefixvm(this->subPrefix(),"measures.csv")),this->worldComm() )
        //M_PsLogger( new PsLogger(prefixvm(this->prefix(),"PsLogger"),this->worldComm() ) )
    {
        //-----------------------------------------------------------------------//
        // move in stationary mode if we have this relation
        if ( M_timeInitial + M_timeStep == M_timeFinal)
            M_isStationary=true;
        //-----------------------------------------------------------------------//
        if ( Environment::vm().count(prefixvm(this->prefix(),"ginac-expr-directory").c_str()) )
            M_directoryLibSymbExpr = Environment::rootRepository()+"/"+soption(_name="ginac-expr-directory",_prefix=this->prefix());
        else
            M_directoryLibSymbExpr = (fs::path(this->appliRepositoryWithoutNumProc() )/fs::path("symbolic_expr")).string();
        //-----------------------------------------------------------------------//
        // mesh file : .msh
        if (Environment::vm().count(prefixvm(this->prefix(),"mshfile").c_str()))
            M_mshFileStr = Environment::expand( soption(_prefix=this->prefix(),_name="mshfile") );
        // mesh file : .geo
        if (Environment::vm().count(prefixvm(this->prefix(),"geofile").c_str()))
            M_geoFileStr = Environment::expand( soption(_prefix=this->prefix(),_name="geofile") );
        //-----------------------------------------------------------------------//
        if (soption(_prefix=this->prefix(),_name="geomap")=="opt")
            M_geomap=GeomapStrategyType::GEOMAP_OPT;
        else
            M_geomap=GeomapStrategyType::GEOMAP_HO;
        //-----------------------------------------------------------------------//
        M_modelProps = std::make_shared<ModelProperties>( Environment::expand( soption( _name=prefixvm(this->prefix(),"filename")) ),
                                                          M_directoryLibSymbExpr, this->worldComm() );
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
        M_modelProps->parameters()["t"] = ModelParameter("current_time",M_timeCurrent);
    }



    void
    ModelNumerical::saveMSHfilePath( std::string const& fileSavePath, std::string const& meshPath ) const
    {
        std::string meshPathUsed = (meshPath.empty())? this->mshfileStr() : meshPath;
        if (this->verbose()) FeelModels::Log(this->prefix()+".ModelNumerical","saveMSHfilePath",
                                             "fileSavePath :"+ fileSavePath + "\nwrite :\n" + meshPathUsed,
                                             this->worldComm(),this->verboseAllProc());

        if ( this->worldComm().isMasterRank() )
        {
            fs::path thedir = fs::path( fileSavePath ).parent_path();
            if ( !fs::exists(thedir))
                fs::create_directories(thedir);

            std::ofstream file(fileSavePath.c_str(), std::ios::out);
            file << meshPathUsed;
            file.close();
        }
    }


} // namespace FeelModels

} // namespace Feel

