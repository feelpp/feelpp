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

#ifndef FEELPP_FILTERS_IMPORTERCSV_HPP
#define FEELPP_FILTERS_IMPORTERCSV_HPP 1

#include <feel/feelcore/feel.hpp>

namespace Feel
{

class ImporterCSV
{
    typedef boost::tokenizer< boost::escaped_list_separator<char> > tokenizer_type;
public :
    ImporterCSV( std::string const& filename, bool readOnlyHeader = false );
    ImporterCSV( ImporterCSV const& ) = default;
    ImporterCSV( ImporterCSV && ) = default;

    //! read the next line in csv file, arg clearData is used to clean previous data read
    void readLine( bool clearData = false );

    //! return true if csv reading is finished
    bool readIsFinished() const;

    //! return names in header
    std::vector<std::string> const& names() const { return M_dataIndexToName; }

    //! return true if csv contains a column with the name key
    bool hasName( std::string const& key ) const { return M_dataNameToIndex.find( key ) != M_dataNameToIndex.end(); }

    //! return the value at line read lineId in the column key
    template<typename T>
    T
    value( size_type lineId, std::string const& key ) const
        {
            return boost::lexical_cast<T>( this->valueToString( lineId, key ) );
        }
    //! return the value at the last line read in the column key
    template<typename T>
    T
    value( std::string const& key ) const
        {
            return this->value<T>( invalid_size_type_value, key );
        }

    //! return the number of line in memory
    size_type numberOfLineInMemory() const { return M_data.size(); }

private :

    std::string const& valueToString( size_type lineId, std::string const& key ) const;

private :
    std::ifstream M_ifs;
    std::map<std::string,uint16_type> M_dataNameToIndex;
    std::vector<std::string> M_dataIndexToName;
    std::vector<std::vector<std::string>> M_data;
};


} // namespace Feel

#endif
