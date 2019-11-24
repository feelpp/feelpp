/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2019-04-08

  Copyright (C) 2019 Feel++ Consortium

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

#include <boost/algorithm/string.hpp>
#include <feel/feelfilters/importercsv.hpp>

namespace Feel
{

ImporterCSV::ImporterCSV( std::string const& filename, bool readOnlyHeader )
{
    M_ifs.open( filename );
    CHECK( M_ifs.is_open() ) << "fail to open file : " << filename << "\n";

    std::string lineHeader;
    std::getline( M_ifs, lineHeader );

    boost::trim(lineHeader);
    tokenizer_type tok( lineHeader );
    M_dataIndexToName.assign( tok.begin(), tok.end() );
    for ( int k=0;k<M_dataIndexToName.size();++k )
    {
        M_dataNameToIndex[ M_dataIndexToName[k] ] = k;
    }
    if ( !readOnlyHeader )
    {
        while ( !this->readIsFinished() )
        {
            this->readLine();
        }
    }
}

void
ImporterCSV::readLine( bool clearData )
{
    if ( this->readIsFinished() )
        return;

    std::string line;
    std::getline( M_ifs, line );
    boost::trim(line);
    tokenizer_type tok( line );
    
    std::vector<std::string> vec;
    vec.assign( tok.begin(), tok.end() );
    if ( vec.empty() && this->readIsFinished() )
        return;

    CHECK( M_dataIndexToName.size() == vec.size() ) << "invalid number of line entry" << M_dataIndexToName.size() << " vs " << vec.size();
    if ( clearData )
        M_data.clear();
    M_data.push_back( vec );
}

bool
ImporterCSV::readIsFinished() const
{
    return !M_ifs.is_open() || !M_ifs.good();
}

std::string const&
ImporterCSV::valueToString( size_type lineId, std::string const& key ) const
{
    if ( lineId == invalid_size_type_value && !M_data.empty() )
        lineId = M_data.size() - 1;
    auto itFindKey = M_dataNameToIndex.find( key );
    CHECK( itFindKey != M_dataNameToIndex.end() ) << "key " << key << "is not present";
    CHECK( lineId < M_data.size() ) << "lineId " << lineId << " too large, should be less than " << M_data.size();
    return M_data[lineId][ itFindKey->second ];
}

} // namespace Feel
