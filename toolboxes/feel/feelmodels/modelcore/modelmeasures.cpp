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

nl::json
ModelMeasuresStorage::metadata() const
{
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
    return jmeta;
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

        std::ofstream ofile( (pdir/"metadata.json").string(), std::ios::trunc);
        if ( ofile )
        {
            auto jmeta = this->metadata();
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

    std::unique_ptr<fs::ifstream> ifile;
    std::unique_ptr<fs::ofstream> ofile;

    fs::path pfilename = fs::path(filename);

    if ( isLastIndex )
    {
        ifile = std::make_unique<fs::ifstream>( pfilename );
    }
    else
    {
        fs::path pfilename_tmp = fs::path(filename);
        pfilename_tmp.replace_extension( ".tmp.csv" );
        fs::rename( pfilename, pfilename_tmp );

        ifile = std::make_unique<fs::ifstream>( pfilename_tmp );
        ofile = std::make_unique<fs::ofstream>( pfilename, std::ios::trunc );
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




ModelMeasuresFlowRate::ModelMeasuresFlowRate()
{
    this->setDirection( "interior_normal" );
}

void
ModelMeasuresFlowRate::setDirection( std::string const& dir )
{
    CHECK( dir == "interior_normal" || dir == "exterior_normal" ) << "invalid dir " << dir;
    M_direction = dir;
}

void
ModelMeasuresFlowRate::setup( nl::json const& jarg, std::string const& name, ModelIndexes const& indexes )
{
    M_name = name;

    if ( jarg.contains("markers") )
        M_markers.setup( jarg.at("markers"), indexes);

    if ( jarg.contains("direction" ) )
    {
        auto const& j_dir = jarg.at("direction");
        if ( j_dir.is_string() )
            this->setDirection( indexes.replace( j_dir.get<std::string>() ) );
    }
}


void
ModelMeasuresNormalFluxGeneric::setup( nl::json const& jarg, std::string const& name, ModelIndexes const& indexes )
{
    M_name = name;

    if ( jarg.contains("markers") )
        M_markers.setup( jarg.at("markers"), indexes);

    if ( jarg.contains("direction") )
    {
        auto const& j_direction = jarg.at("direction");
        if ( j_direction.is_string() )
            M_direction = indexes.replace( j_direction.get<std::string>() );
    }
    else
        M_direction = "outward";

    CHECK( M_direction == "inward" || M_direction == "outward" ) << "invalid dir " << M_direction;

    if ( jarg.contains( "requires_markers_connection") )
        M_requiresMarkersConnection.setup( jarg.at( "requires_markers_connection"), indexes );

    if ( jarg.contains( "internalfaces_evaluation") )
    {
          auto const& j_ifevaltype = jarg.at( "internalfaces_evaluation");
          if ( j_ifevaltype.is_string() )
              M_internalFacesEvalutationType = indexes.replace( j_ifevaltype.get<std::string>() );
    }
}


} // namespace FeelModels
} // namespace Feel
