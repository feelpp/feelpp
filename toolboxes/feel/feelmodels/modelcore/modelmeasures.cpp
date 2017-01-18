/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2015-12-02

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
 \file modelpostprocessextra.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2015-12-02
 */

#include <feel/feelmodels/modelcore/modelmeasures.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename T>
std::vector<T> as_vector(pt::ptree const& pt, pt::ptree::key_type const& key)
{
    std::vector<T> r;
    for (auto& item : pt.get_child(key))
        r.push_back(item.second.template get_value<T>());
    return r;
}



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




ModelMeasuresFlowRate::ModelMeasuresFlowRate()
{
    this->setDirection( "interior_normal" );
}

void
ModelMeasuresFlowRate::addMarker( std::string const& mark )
{
    if ( std::find( M_meshMarkers.begin(),M_meshMarkers.end(), mark ) == M_meshMarkers.end() )
        M_meshMarkers.push_back( mark );
}
void
ModelMeasuresFlowRate::setDirection( std::string const& dir )
{
    CHECK( dir == "interior_normal" || dir == "exterior_normal" ) << "invalid dir " << dir;
    M_direction = dir;
}

void
ModelMeasuresFlowRate::setup( pt::ptree const& ptree, std::string const& name )
{
    std::vector<std::string> markerList = as_vector<std::string>( ptree, "markers" );
    if ( markerList.empty() )
    {
        std::string markerUnique = ptree.get<std::string>( "markers" );
        if ( !markerUnique.empty() )
            markerList = { markerUnique };
    }
    std::string direction = ptree.get<std::string>( "direction" );
    //std::cout << "markerList " << markerList.front() << " direction " << direction << "\n";
    this->setName( name );
    for ( auto const& marker : markerList )
        this->addMarker( marker );
    this->setDirection( direction );
}


} // namespace FeelModels
} // namespace Feel
