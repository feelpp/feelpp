/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 17 September 2019

 Copyright (C) 2019 Feel++ Consortium

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

#include <feel/feelmodels/modelindexes.hpp>

namespace Feel {

std::vector<std::string>
ModelIndexes::generateIndexFromString( std::string const& input )
{
    std::vector<std::string> res;
    if ( input.empty() )
        return res;
    boost::char_separator<char> sep(":");
    boost::tokenizer<boost::char_separator<char> > kvlist( input, sep);
    int sizeOfSplit = std::distance(kvlist.begin(),kvlist.end());
    if ( sizeOfSplit == 1 )
        res.push_back( input );
    else
    {
        try {
        std::vector<int> rangeDef;
        for( auto const& r : kvlist )
            rangeDef.push_back(std::stoi(r));
        if( sizeOfSplit == 2 )
            rangeDef.push_back(1);
        auto rangeOfInt = boost::irange(rangeDef[0],rangeDef[1],rangeDef[2]);
        for( auto const& i : rangeOfInt )
            res.push_back( std::to_string(i) );
        }
        catch ( std::invalid_argument const& )
        {
            // the char : come from an expression
            res.push_back( input );
        }
    }
    return res;
}

std::vector<ModelIndexes>
ModelIndexes::generateAllCases( nl::json const& jarg, int startIndex )
{
    // evaluate each index
    int currentIndex = startIndex;
    std::map<int,std::vector<std::vector<std::string>>> M_indexes;
    while ( 1 )
    {
        std::string keyIndex = (boost::format( "index%1%")%currentIndex).str();
        if ( !jarg.contains( keyIndex ) )
             break;
        auto const& jKeyIndex = jarg.at( keyIndex );

        M_indexes[currentIndex].clear();
        if ( jKeyIndex.is_string() ) // string value case
        {
            std::vector<std::string> indexReaded = generateIndexFromString( jKeyIndex.get<std::string>() );
            for ( std::string const& s : indexReaded )
                M_indexes[currentIndex].push_back( { s } );
        }
        else if ( jKeyIndex.is_array() ) // array case
        {
            int nSubKey = -1;
            for ( auto const& el : jKeyIndex.items() )
            {
                auto const& el_val = el.value();
                if ( el_val.is_array() ) // array of array
                {
                    int nItem = std::distance(el_val.items().begin(),el_val.items().end());
                    if ( nSubKey > 0 )
                        CHECK( nSubKey == nItem ) << "number of subkey is different (probably a subarray has not the good size)";
                    else
                        nSubKey = nItem;
                    std::vector<std::string> currentIndexValues;
                    currentIndexValues.reserve( nItem );
                    for ( auto const& el2 : el_val.items() )
                    {
                        auto const& el2_val = el2.value();
                        CHECK( el2_val.is_string() ) << "TODO : take into account other type";
                        std::vector<std::string> indexReaded = generateIndexFromString( el2_val );
                        for ( std::string const& s : indexReaded )
                            currentIndexValues.push_back( s );
                    }
                    M_indexes[currentIndex].push_back( currentIndexValues );
                }
                else
                {
                    if ( nSubKey > 0 )
                        CHECK( nSubKey == 1 ) << "number of subkey is different (probably a subarray has not the good size)";
                    else
                        nSubKey = 1;
                    CHECK( el_val.is_string() ) << "TODO : take into account other type";
                    std::vector<std::string> indexReaded = generateIndexFromString( el_val );
                    for ( std::string const& s : indexReaded )
                        M_indexes[currentIndex].push_back( { s } );
                }
            }
        }
        else
            CHECK( false ) << "TODO : an other type " << jKeyIndex;

        if ( M_indexes[currentIndex].empty() )
            break;
        ++currentIndex;
    }

    // create mapping tag to replacement strings for each index
    std::map< int, std::map<std::string, std::vector<std::string> > > indexTagToReplacementStrings;
    for ( auto const& [ indexId, indexes ] : M_indexes )
    {
        CHECK( !indexes.empty() ) << "something wrong";
        int sizeOfSubArray = indexes[0].size();
        if ( sizeOfSubArray == 1 )
        {
            std::string indexTag = "%"+std::to_string(indexId)+"%";
            for ( int k=0;k<indexes.size();++k )
            {
                CHECK( indexes[k].size() == sizeOfSubArray ) << "wrong subarray size";
                indexTagToReplacementStrings[indexId][indexTag].push_back( indexes[k][0] );
            }
        }
        else
        {
            for ( int k=0;k<indexes.size();++k )
            {
                for ( int l=0;l<indexes[k].size();++l )
                {
                    std::string indexTag = "%"+std::to_string(indexId)+"_"+std::to_string(l+1)+"%";
                    indexTagToReplacementStrings[indexId][indexTag].push_back( indexes[k][l] );
                }
            }
        }
    }

    // generate all cases
    std::vector<ModelIndexes> crsfiAllCases;
    if ( indexTagToReplacementStrings.size() == 0 )
    {
        crsfiAllCases.push_back( ModelIndexes() );
    }
    else
    {
        for ( auto const& [ indexId, tagToReplaceStrings ] : indexTagToReplacementStrings )
        {
            std::vector<ModelIndexes> crsfiOnIndex;

            std::vector<std::string> listTags;
            int sizeOfSubArray{0};
            for ( auto const& [ indexTag, indexes ] : tagToReplaceStrings )
            {
                listTags.push_back( indexTag );
                sizeOfSubArray = indexes.size();
            }

            for (int k= 0;k<sizeOfSubArray;++k)
            {
                ModelIndexes _tmp;
                for ( std::string const& indexTag : listTags )
                {
                    auto const& indexes = tagToReplaceStrings.find( indexTag )->second;
                    _tmp[indexTag] = indexes[k];
                }
                crsfiOnIndex.push_back( _tmp );
            }

            if ( crsfiAllCases.empty() )
            {
                crsfiAllCases = crsfiOnIndex;
            }
            else
            {
                // cartesian product
                std::vector<ModelIndexes> crsfiNew;
                for ( auto const& crsfi1 : crsfiAllCases )
                {
                    for ( auto const& crsfi2 : crsfiOnIndex )
                    {
                        ModelIndexes _tmp;
                        for ( auto const& [key,value] : crsfi1 )
                            _tmp[key]=value;
                        for ( auto const& [key,value] : crsfi2 )
                            _tmp[key]=value;
                        crsfiNew.push_back( _tmp );
                    }
                }
                crsfiAllCases = crsfiNew;
            }
        }
    }

    for ( auto & indexToUp : crsfiAllCases )
        indexToUp.setNextFreeIndex( currentIndex );

    return crsfiAllCases;
}


 std::vector<ModelIndexes>
 ModelIndexes::generateAllCases( nl::json const& jarg, ModelIndexes const& indexes )
 {
     std::vector<ModelIndexes> res = ModelIndexes::generateAllCases( jarg, indexes.nextFreeIndex() );
     for ( auto & indexesToUp : res )
         for ( auto const& [key,value] : indexes )
             indexesToUp[key] = value;
     return res;
 }


std::vector<ModelIndexes>
ModelIndexes::generateAllCases( pt::ptree const& pt, int startIndex )
{
    // evaluate each index
    int currentIndex = startIndex;
    std::map<int,std::vector<std::vector<std::string>>> M_indexes;
    while ( 1 )
    {
        std::string keyIndex = (boost::format( "index%1%")%currentIndex).str();
        auto indexPTree = pt.get_child_optional(keyIndex);
        if ( !indexPTree )
            break;

        M_indexes[currentIndex].clear();
        if ( indexPTree->empty() ) // value case
        {
            std::vector<std::string> indexReaded = generateIndexFromString( indexPTree->get_value<std::string>() );
            for ( std::string const& s : indexReaded )
                M_indexes[currentIndex].push_back( { s } );
        }
        else // array case
        {
            bool isAnArrayOfArray = !indexPTree->begin()->second.empty();
            //int sizeOfSubArray = indexPTree->begin()->second.size();
            for ( auto const& item : *indexPTree )
            {
                CHECK( item.first.empty() ) << "should be an array, not a subtree";
                if ( isAnArrayOfArray )
                {
                    int sizeOfSubArray = -1;
                    std::vector<std::string> currentIndexValues;
                    for ( auto const& item2 : item.second )
                    {
                        CHECK( item2.first.empty() ) << "should be an array, not a subtree";
                        std::vector<std::string> indexReaded = generateIndexFromString( item2.second.template get_value<std::string>() );
                        for ( std::string const& s : indexReaded )
                            currentIndexValues.push_back( s );
                        //currentIndexValues.push_back( item2.second.template get_value<std::string>() );
                    }
                    if ( sizeOfSubArray < 0 )
                        sizeOfSubArray = currentIndexValues.size();
                    CHECK( currentIndexValues.size() == sizeOfSubArray ) << "all subarray should be have the same size : " << currentIndexValues.size() << " vs " <<  sizeOfSubArray;
                    M_indexes[currentIndex].push_back( currentIndexValues );
                }
                else
                {
                    std::vector<std::string> indexReaded = generateIndexFromString( item.second.template get_value<std::string>() );
                    for ( std::string const& s : indexReaded )
                        M_indexes[currentIndex].push_back( { s } );
                    //M_indexes[currentIndex].push_back( { item.second.template get_value<std::string>() } );
                }
            }
        }

        if ( M_indexes[currentIndex].empty() )
            break;
        ++currentIndex;
    }

    // create mapping tag to replacement strings for each index
    std::map< int, std::map<std::string, std::vector<std::string> > > indexTagToReplacementStrings;
    for ( auto const& [ indexId, indexes ] : M_indexes )
    {
        CHECK( !indexes.empty() ) << "something wrong";
        int sizeOfSubArray = indexes[0].size();
        if ( sizeOfSubArray == 1 )
        {
            std::string indexTag = "%"+std::to_string(indexId)+"%";
            for ( int k=0;k<indexes.size();++k )
            {
                CHECK( indexes[k].size() == sizeOfSubArray ) << "wrong subarray size";
                indexTagToReplacementStrings[indexId][indexTag].push_back( indexes[k][0] );
            }
        }
        else
        {
            for ( int k=0;k<indexes.size();++k )
            {
                for ( int l=0;l<indexes[k].size();++l )
                {
                    std::string indexTag = "%"+std::to_string(indexId)+"_"+std::to_string(l+1)+"%";
                    indexTagToReplacementStrings[indexId][indexTag].push_back( indexes[k][l] );
                }
            }
        }
    }

    // generate all cases
    std::vector<ModelIndexes> crsfiAllCases;
    if ( indexTagToReplacementStrings.size() == 0 )
    {
        crsfiAllCases.push_back( ModelIndexes() );
    }
    else
    {
        for ( auto const& [ indexId, tagToReplaceStrings ] : indexTagToReplacementStrings )
        {
            std::vector<ModelIndexes> crsfiOnIndex;

            std::vector<std::string> listTags;
            int sizeOfSubArray{0};
            for ( auto const& [ indexTag, indexes ] : tagToReplaceStrings )
            {
                listTags.push_back( indexTag );
                sizeOfSubArray = indexes.size();
            }

            for (int k= 0;k<sizeOfSubArray;++k)
            {
                ModelIndexes _tmp;
                for ( std::string const& indexTag : listTags )
                {
                    auto const& indexes = tagToReplaceStrings.find( indexTag )->second;
                    _tmp[indexTag] = indexes[k];
                }
                crsfiOnIndex.push_back( _tmp );
            }

            if ( crsfiAllCases.empty() )
            {
                crsfiAllCases = crsfiOnIndex;
            }
            else
            {
                // cartesian product
                std::vector<ModelIndexes> crsfiNew;
                for ( auto const& crsfi1 : crsfiAllCases )
                {
                    for ( auto const& crsfi2 : crsfiOnIndex )
                    {
                        ModelIndexes _tmp;
                        for ( auto const& [key,value] : crsfi1 )
                            _tmp[key]=value;
                        for ( auto const& [key,value] : crsfi2 )
                            _tmp[key]=value;
                        crsfiNew.push_back( _tmp );
                    }
                }
                crsfiAllCases = crsfiNew;
            }
        }
    }

    for ( auto & indexToUp : crsfiAllCases )
        indexToUp.setNextFreeIndex( currentIndex );

    return crsfiAllCases;
}


 std::vector<ModelIndexes>
 ModelIndexes::generateAllCases( pt::ptree const& pt, ModelIndexes const& indexes )
 {
     std::vector<ModelIndexes> res = ModelIndexes::generateAllCases( pt, indexes.nextFreeIndex() );
     for ( auto & indexesToUp : res )
         for ( auto const& [key,value] : indexes )
             indexesToUp[key] = value;
     return res;
 }


} // namespace Feel
