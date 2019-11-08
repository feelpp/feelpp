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
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 18 Jun 2019
//! @copyright 2019 Feel++ Consortium
//!

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/importercsv.hpp>

namespace Feel {


void
SensorDescriptionMap::read()
{
    std::string fileSensorPos = soption(_name="sensor-position");
    ImporterCSV readerCSV_position( fileSensorPos );

    std::set<std::pair<size_type,std::string>> sensorUsedInPosFile;
    for ( size_type k=0;k< readerCSV_position.numberOfLineInMemory();++k )
    {
        std::string isInBatBStr = readerCSV_position.value<std::string>( k, "isInBatB" );
        if ( isInBatBStr.empty() )
            continue;
        if ( !readerCSV_position.value<bool>( k, "isInBatB" ) )
            continue;

        std::string sensorArchi = readerCSV_position.value<std::string>( k, "archi" );
        std::vector<std::string> sensorArchiSplitted;
        boost::split( sensorArchiSplitted, sensorArchi, boost::is_any_of(":"), boost::token_compress_on );
        if ( sensorArchiSplitted.empty() )
            continue;

        //std::string sensorName = "zigduino-" + readerCSV_position.value<std::string>( k, "node_id" );
        std::string sensorName = sensorArchiSplitted[0] + "-" + readerCSV_position.value<std::string>( k, "custom_id" );
        if ( readerCSV_measures.hasName( sensorName ) )
            sensorUsedInPosFile.insert( std::make_pair( k, sensorName ) );
    }
    std::vector<double> center(3);
    for ( auto const& [k, sensorName] : sensorUsedInPosFile )
    {


        center[0] = readerCSV_position.value<double>( k, "X_m"/*"x"*/ );
        center[1] = readerCSV_position.value<double>( k, "Y_m"/*"y"*/ );
        center[2] = readerCSV_position.value<double>( k, "Z_m"/*"z"*/ );

    }


    if ( Environment::isMasterRank() )
    {
        std::cout << "---------------------------------------------------\n";
        std::cout << "sensor name used in the visualization : \n";
        for ( auto [id,sensorName] : sensorUsedInPosFile )
            std::cout << sensorName << "\n";
        std::cout << "---------------------------------------------------\n";
    }



}



}
