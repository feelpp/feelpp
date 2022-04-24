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
#ifndef FEELPP_MODELS_MODELINDEXES_HPP
#define FEELPP_MODELS_MODELINDEXES_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/json.hpp>

namespace Feel {

class FEELPP_EXPORT ModelIndexes : public std::map<std::string,std::string>
{
  public :
    ModelIndexes() : M_nextFreeIndex( 1 ) {}
    ModelIndexes( ModelIndexes const& ) = default;
    ModelIndexes( ModelIndexes&& ) = default;
    ModelIndexes& operator=( ModelIndexes const& ) = default;
    ModelIndexes& operator=( ModelIndexes && ) = default;

    std::string replace( std::string const& input ) const
        {
            std::string res = input;
            for ( auto const& [sin,sout] : *this )
                boost::replace_all( res,  sin, sout );
            return res;
        }

    int nextFreeIndex() const { return M_nextFreeIndex; }
    void setNextFreeIndex( int i ) { M_nextFreeIndex = i; }

    static std::vector<ModelIndexes> generateAllCases( nl::json const& jarg, int startIndex = 1 );
    static std::vector<ModelIndexes> generateAllCases( nl::json const& jarg, ModelIndexes const& indexes );
    static std::vector<ModelIndexes> generateAllCases( pt::ptree const& pt, int startIndex = 1 );
    static std::vector<ModelIndexes> generateAllCases( pt::ptree const& pt, ModelIndexes const& indexes );

  private :

    static std::vector<std::string> generateIndex( nl::json const& input );

  private :
    int M_nextFreeIndex;
};


} // namespace Feel

#endif
