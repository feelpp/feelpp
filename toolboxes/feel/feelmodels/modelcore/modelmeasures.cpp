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

#include <boost/algorithm/string.hpp>

namespace Feel
{
namespace FeelModels
{


bool
ModelMeasuresStorage::hasValue( std::string const& name, std::string const& key ) const
{
    auto itFindName = M_values.find( name );
    if ( itFindName == M_values.end() )
        return false;
    return itFindName->second.hasValue( key );
}

typename ModelMeasuresStorage::value_type
ModelMeasuresStorage::value( std::string const& name, std::string const& key ) const
{
    return M_values.at( name ).value( key );
}


void
ModelMeasuresStorage::setTable( std::string const& name, Feel::Table && table )
{
    if ( M_tables.find( name ) == M_tables.end() )
        M_tables.emplace( std::make_pair( name, ModelMeasuresStorageTable(name) ) );
    M_tables.at( name ).setTable( std::move( table ) );
}

void
ModelMeasuresStorage::save( double time )
{
    if ( !this->isUpdated() )
        return;

    // TODO check the time
    M_times.push_back( time );

    if ( M_worldComm->isMasterRank() )
    {
        fs::path pdir(M_directory);
        if ( !fs::exists( pdir ) )
            fs::create_directories( pdir );

        for ( auto & [name, values] : M_values )
            if ( values.isUpdated() )
                values.saveCSV( M_directory );
        for ( auto & [name, table] : M_tables )
            if ( table.isUpdated() )
                table.saveCSV( M_directory,M_times.size()-1 );

        nl::json jmeta;
        jmeta["times"] = M_times;
        nl::json j_measures_values, j_measures_tables;
        for ( auto const& [name, values] : M_values )
            j_measures_values[name] = { { "format", "csv" } };
        for ( auto const& [name, table] : M_tables )
            j_measures_tables[name] = { { "format", "csv" } };
        nl::json j_measures;
        if ( !j_measures_values.is_null() )
            j_measures["values"] = j_measures_values;
        if ( !j_measures_tables.is_null() )
            j_measures["tables"] = j_measures_tables;
        if ( !j_measures.is_null() )
            jmeta["measures"] = j_measures;

        std::ofstream ofile( (pdir/"metadata.json").string(), std::ios::trunc);
        if ( ofile )
        {
            ofile << jmeta.dump(1);
            ofile.close();
        }
    }
    M_worldComm->barrier();
}

bool
ModelMeasuresStorage::isUpdated() const
{
    for ( auto const& [name, values] : M_values )
        if ( values.isUpdated() )
            return true;
    for ( auto const& [name, table] : M_tables )
        if ( table.isUpdated() )
            return true;
    return false;
}

void
ModelMeasuresStorage::resetState()
{
    for ( auto & [name, values] : M_values )
        values.resetState();
    for ( auto & [name, table] : M_tables )
        table.resetState();
}

void
ModelMeasuresStorage::updateParameterValues( std::map<std::string,double> & mp, std::string const& prefix_symbol ) const
{
    for ( auto const& [name, values] : M_values )
        values.updateParameterValues( mp, prefix_symbol );
}

void
ModelMeasuresStorage::restart( double time )
{
    std::vector<std::string> storageValuesNames;
    if ( M_worldComm->isMasterRank() )
    {
        this->restartSeqImpl( time );
        for ( auto & [name, values] : M_values )
            storageValuesNames.push_back( name );
    }

    if ( M_worldComm->globalComm().size() > 1 )
    {
        auto dataBroadcasted = std::forward_as_tuple( M_times, storageValuesNames );
        boost::mpi::broadcast( M_worldComm->globalComm(), dataBroadcasted, M_worldComm->masterRank() );

        if ( !M_worldComm->isMasterRank() )
        {
            for ( std::string const& name : storageValuesNames )
            {
                ModelMeasuresStorageValues mmsv( name );
                M_values.erase( name );
                M_values.emplace( std::make_pair( name, std::move( mmsv ) ) );
            }
        }

        for ( auto & [name, values] : M_values )
            boost::mpi::broadcast( M_worldComm->globalComm(), values, M_worldComm->masterRank() );
    }

    this->resetState();
}

void
ModelMeasuresStorage::restartSeqImpl( double time )
{
    fs::path pdir(M_directory);
    std::ifstream ifile( (pdir/"metadata.json").string() );
    if ( !ifile )
        return;

    nl::json jmeta;
    ifile >> jmeta;
    if ( jmeta.contains("times") )
        M_times = jmeta.at( "times" ).get<std::vector<double>>();

    int restartIndex=0;
    for ( int k=0; k<M_times.size(); ++k,++restartIndex )
        if ( std::abs(M_times[k]-time ) < 1e-8 )
            break;

    if ( M_times.empty() )
        return;

    CHECK( restartIndex < M_times.size() ) << "wrong time for restart";
    bool isLastIndex = restartIndex == (M_times.size() - 1);

    if ( jmeta.contains("measures") )
    {
        auto const& j_measures = jmeta.at("measures");
        if ( j_measures.contains( "values" ) )
        {
            for ( auto const& el : j_measures.at("values").items() )
            {
                std::string const& name = el.key();
                ModelMeasuresStorageValues mmsv( name );
                mmsv.restart( M_directory, restartIndex, isLastIndex );
                M_values.erase( name );
                M_values.emplace( std::make_pair( name, std::move( mmsv ) ) );
            }
        }
    }
}


void
ModelMeasuresStorageValues::saveCSV( std::string const& directory, bool append )
{
    fs::path pdir(directory);
    if ( !fs::exists( pdir ) )
        fs::create_directories( pdir );
    std::ofstream ofile( this->filename_CSV( pdir ), append && !M_keys.empty()? std::ios::app : std::ios::trunc);
    if ( ofile )
        this->exportCSV( ofile );
}
void
ModelMeasuresStorageValues::exportCSV( std::ostream &o )
{
    if ( M_data.empty() )
        return;

    Feel::Table table;
    table.format().setFloatingPointPrecision( 16 );
    if ( M_keys.empty() )
    {
        table.format().setFirstRowIsHeader( true );
        table.resize( 2, M_data.size() );
        M_keys.reserve( M_data.size() );
        for ( auto const& [key,val] : M_data )
        {
            int i = M_keys.size();
            M_keys.push_back( key );
            table(0,i) = key;
            table(1,i) = val;
        }
    }
    else
    {
        CHECK( M_data.size() == M_keys.size() ) << "number of key are different from previous save (" << M_data.size() << " != " << M_keys.size() << ") : not compatible with csv format";
        table.resize( 1, M_data.size() );
        for ( int i=0;i<M_keys.size();++i )
        {
            auto const& key = M_keys[i];
            auto itFindVal = M_data.find( key );
            CHECK( itFindVal != M_data.end() ) << "not find key "<<key << "in data values";
            table(0,i) = itFindVal->second;
        }
    }
    table.exportCSV( o );
}

void
ModelMeasuresStorageValues::updateParameterValues( std::map<std::string,double> & mp, std::string const& prefix_symbol ) const
{
    for ( auto const& [key,val] : M_data )
        mp.emplace( std::make_pair( prefixvm(prefix_symbol,key,"_"), val ) );
}


void
ModelMeasuresStorageValues::restart( std::string const& directory, int restartIndex, bool isLastIndex )
{
    std::string filename = this->filename_CSV( directory );

    std::unique_ptr<std::ifstream> ifile;
    std::unique_ptr<std::ofstream> ofile;

    fs::path pfilename = fs::path(filename);

    if ( isLastIndex )
    {
        ifile = std::make_unique<std::ifstream>( pfilename );
    }
    else
    {
        fs::path pfilename_tmp = fs::path(filename);
        pfilename_tmp.replace_extension( ".tmp.csv" );
        fs::rename( pfilename, pfilename_tmp );

        ifile = std::make_unique<std::ifstream>( pfilename_tmp );
        ofile = std::make_unique<std::ofstream>( pfilename, std::ios::trunc );
    }

    std::string lineRead;
    std::vector<std::string> lineReadSplitted;

    // first line
    std::getline( *ifile, lineRead );
    if ( ofile )
        *ofile << lineRead << "\n";

    boost::split( M_keys, lineRead, boost::is_any_of(","), boost::token_compress_on );
    for ( auto & key : M_keys )
        boost::trim( key );
    std::vector<std::string> idToKey( M_keys.size() );
    for (int k=0;k<M_keys.size();++k)
        idToKey[k] = M_keys[k];

    int currentIndex=0;
    while ( !ifile->eof() )
    {
        if ( currentIndex > restartIndex )
        {
            break;
        }
        std::getline(*ifile,lineRead);
        if ( ofile )
            *ofile << lineRead << "\n";
        ++currentIndex;
    }

    boost::split( lineReadSplitted, lineRead, boost::is_any_of(","), boost::token_compress_on );
    CHECK( lineReadSplitted.size() == idToKey.size() ) << "something wrong";
    for (int k=0;k<idToKey.size();++k)
        this->setValue( idToKey[k], std::stod( lineReadSplitted[k] ) );
}


void
ModelMeasuresStorageTable::saveCSV( std::string const& directory, size_type index )
{
    fs::path pdir(directory);
    if ( !fs::exists( pdir ) )
        fs::create_directories( pdir );

    std::ofstream ofile( this->filename_CSV( pdir, index ), std::ios::trunc);
    M_table.exportCSV( ofile );
}




template <typename T>
std::vector<T> as_vector(pt::ptree const& pt, pt::ptree::key_type const& key)
{
    std::vector<T> r;
    for (auto& item : pt.get_child(key))
        r.push_back(item.second.template get_value<T>());
    return r;
}


#if 0
ModelMeasuresIO::ModelMeasuresIO( std::string const& pathFile, worldcomm_ptr_t const& worldComm )
    :
    M_worldComm( worldComm ),
    M_pathFile( pathFile ),
    M_addNewDataIsLocked( false )
{}

void
ModelMeasuresIO::clear()
{
    M_dataNameToIndex.clear();
    M_dataIndexToName.clear();
    M_data.clear();
    M_addNewDataIsLocked = false;
}

void
ModelMeasuresIO::writeHeader()
{
    if ( M_dataNameToIndex.empty() )
        return;
    if ( M_worldComm->isMasterRank() )
    {
        bool hasAlreadyWrited = false;
        std::ofstream fileWrited(M_pathFile, std::ios::out | std::ios::trunc);
        for ( auto const& dataName : M_dataIndexToName )
        {
            int spacing = std::max(28, int(dataName.size()+2) );
            if ( hasAlreadyWrited )
                fileWrited << ",";
            fileWrited << std::setw( spacing ) << std::right << dataName;
            hasAlreadyWrited = true;
        }
        fileWrited << std::endl;
        fileWrited.close();
    }
    M_addNewDataIsLocked = true;
}
void
ModelMeasuresIO::restart( std::string const& paramKey, double val )
{
    if ( M_worldComm->isMasterRank() && fs::exists( M_pathFile ) )
    {
        double ti = val;//this->timeInitial();
        std::ifstream fileI(M_pathFile, std::ios::in);
        std::ostringstream buffer;

        std::string lineRead;
        std::getline(fileI,lineRead);
        std::vector<std::string> lineReadSplitted;
        boost::split( lineReadSplitted, lineRead, boost::is_any_of(","), boost::token_compress_on );
        uint16_type indexKey = invalid_uint16_type_value;
        M_dataNameToIndex.clear();
        M_data.resize( lineReadSplitted.size(), 0. );
        M_dataIndexToName.resize( lineReadSplitted.size(), "" );
        for ( uint16_type k=0;k<lineReadSplitted.size();++k )
        {
            std::string headerName = lineReadSplitted[k];
            boost::trim(headerName);
            if ( paramKey == headerName )
                indexKey = k;
            M_dataNameToIndex[headerName] = k;
            M_dataIndexToName[k] = headerName;
            //std::cout << "headerName="<<headerName<<"\n";
        }
        buffer << lineRead << std::endl;

        while ( !fileI.eof() )
        {
            std::getline(fileI,lineRead);

            buffer << lineRead << std::endl;
            if ( indexKey == invalid_uint16_type_value )
                continue;

            lineReadSplitted.clear();
            boost::split( lineReadSplitted, lineRead, boost::is_any_of(","), boost::token_compress_on );
            std::string valueReadStr = lineReadSplitted[indexKey];
            boost::trim( valueReadStr );
            double valueRead = std::stod( valueReadStr );
            if ( std::abs(valueRead-ti) < 1e-9 )
                break;
        }

        fileI.close();
        std::ofstream fileW(M_pathFile/*.c_str()*/, std::ios::out | std::ios::trunc);
        fileW << buffer.str();
        fileW.close();
    }

    mpi::broadcast( M_worldComm->globalComm(), M_dataIndexToName, M_worldComm->masterRank() );

    if ( !M_worldComm->isMasterRank() )
    {
        M_dataNameToIndex.clear();
        M_data.resize( M_dataIndexToName.size() );
        for (int k=0;k<M_dataIndexToName.size();++k )
            M_dataNameToIndex[ M_dataIndexToName[k] ] = k;
    }
    M_addNewDataIsLocked = true;
}

void
ModelMeasuresIO::exportMeasures()
{
    if ( M_dataNameToIndex.empty() )
        return;
    if ( !M_addNewDataIsLocked )
        this->writeHeader();
    if ( M_worldComm->isMasterRank() )
    {
        bool hasAlreadyWrited = false;
        std::ofstream fileWrited(M_pathFile, std::ios::out | std::ios::app);
        //for ( auto const& dataValue : M_data )
        for ( int k=0;k<M_data.size();++k )
        {
            double dataValue = M_data[k];
            std::string const& dataName = M_dataIndexToName[k];
            int spacing = std::max(28, int(dataName.size()+2) );
            if ( hasAlreadyWrited )
                fileWrited << ",";
            fileWrited << std::setw(spacing) << std::right << std::setprecision( 16 ) << std::scientific << dataValue;
            hasAlreadyWrited = true;
        }
        fileWrited << std::endl;
        fileWrited.close();
    }
}
void
ModelMeasuresIO::setParameter(std::string const& key,double val)
{
    this->setMeasure( key,val );
}
void
ModelMeasuresIO::setMeasure(std::string const& key,double val)
{
    auto itFindKey = M_dataNameToIndex.find( key );
    if ( itFindKey != M_dataNameToIndex.end() )
        M_data[ itFindKey->second ] = val;
    else if ( !M_addNewDataIsLocked )
    {
        M_dataNameToIndex[ key ] =  M_dataIndexToName.size();
        M_dataIndexToName.push_back( key );
        M_data.push_back( val );
    }
    else CHECK( false ) << "add new data is locked";

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
ModelMeasuresIO::setMeasures( std::map<std::string,double> const& m )
{
    for ( auto const& [key,val] : m )
        this->setMeasure(key,val);
}

double
ModelMeasuresIO::measure( std::string const& key ) const
{
    CHECK( this->hasMeasure( key ) ) << "no measure with key " << key;
    uint16_type k = M_dataNameToIndex.find( key )->second;
    return M_data[k];
}

std::map<std::string,double>
ModelMeasuresIO::currentMeasures() const
{
    std::map<std::string,double> res;
    for ( int k=0;k<M_data.size();++k )
    {
        double dataValue = M_data[k];
        std::string const& dataName = M_dataIndexToName[k];
        res[dataName] = dataValue;
    }
    return res;
}

void
ModelMeasuresIO::updateParameterValues( std::map<std::string,double> & mp, std::string const& prefix_symbol ) const
{
    for ( auto const& [name,val] : this->currentMeasures() )
        mp.emplace( std::make_pair( prefixvm(prefix_symbol,name,"_"), val ) );
}

#endif

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


void
ModelMeasuresNormalFluxGeneric::setup( pt::ptree const& _pt, std::string const& name, ModelIndexes const& indexes )
{
    M_name = name;

    if ( auto ptmarkers = _pt.get_child_optional("markers") )
        M_markers.setPTree(*ptmarkers, indexes);

    if ( auto itDir = _pt.get_optional<std::string>("direction") )
        M_direction = indexes.replace( *itDir );
    else
        M_direction = "outward";

    CHECK( M_direction == "inward" || M_direction == "outward" ) << "invalid dir " << M_direction;
}


} // namespace FeelModels
} // namespace Feel
