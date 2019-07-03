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


#include <array>
#include <string>
#include <map>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/importercsv.hpp>


namespace Feel
{

template<int Dim = 3>
class SensorDescription
{
public:

    static constexpr std::array<double,Dim> zero()
        {
            if constexpr ( Dim == 1 )
                return {0};
            if constexpr ( Dim == 2 )
                return {0,0};
            if constexpr ( Dim == 3 )
                return {0,0,0};
        }
    SensorDescription(std::string const& name,
                      std::string const& type="gaussian",
                      double raduis=1.0,
                      std::array<double,Dim> const& pos = zero() ):
        M_name(name),
        M_type(type),
        M_radius(raduis),
        M_pos(pos)
    {}

    void setName(std::string const& name)
    {
        M_name = name;
    }

    void setType(std::string const& type)
    {
        M_type = type;
    }

    void setRadius(double radius)
    {
        M_radius = radius;
    }

    void setPosition(std::array<double,Dim> const& pos)
    {
        M_pos = pos;
    }
    std::string const& type()
    {
        return M_type;
    }
    std::string const& name()
    {
        return M_name;
    }
    std::array<double,Dim> const& position()
    {
        return M_pos;
    }

    double radius()
    {
        return M_radius;
    }
    virtual ~SensorDescription(){}
private:
    std::string M_name;
    std::string M_type;
    std::array<double,Dim> M_pos;
    double M_radius;
};

template<int Dim=3>
class SensorDescriptionMap: public std::map<std::string,boost::shared_ptr<SensorDescription<Dim>>>
{

public:
   //!
    //! \p file provides the name of the csv file describing the sensor network
    //! if \p N is set to -1 then we read and store all the sensor data dsc
    //!
    SensorDescriptionMap(std::string  file, int N=-1):
        M_file(file),
        M_numberOfSensors(N)
    {
      read();
    }

    void setFile(std::string file)
    {
        M_file = file;
    }

    void setNumberOfSensors(int numberOfSensors)
    {
        M_numberOfSensors = numberOfSensors;
    }
    std::string file() const
    {
        return M_file;
    }

    int numberOfSensors() const
    {
        return M_numberOfSensors;
    }

    virtual ~SensorDescriptionMap() = default;

    std::set<std::pair<size_type,std::string>> read()
    {
        //this->reserve( M_numberOfSensors );
        std::array<double,Dim> center;
        ImporterCSV readerCSV_position( M_file );
        std::set<std::pair<size_type,std::string>> sensorUsedInPosFile;
        double r;
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
            if ( readerCSV_position.hasName( sensorName ) )
            {
               sensorUsedInPosFile.insert( std::make_pair( k, sensorName ) );
               center[0] = readerCSV_position.value<double>( k, "X_m"/*"x"*/ );
               center[1] = readerCSV_position.value<double>( k, "Y_m"/*"y"*/ );
               center[2] = readerCSV_position.value<double>( k, "Z_m"/*"z"*/ );
               r = readerCSV_position.value<double>( k, "radius"/*"r_m"*/ );
               std::string type = readerCSV_position.value<std::string>( k, "type" );
               SensorDescription<Dim> newElement( sensorName, type, r, center );
              this->insert(newElement);
            }

        }

    }

private:
    std::string M_file;
    int M_numberOfSensors;

};


}//namespace feel
