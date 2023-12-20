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

#include <feel/feelmesh/filters.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename MeshType>
class RangeDistributionByMaterialName
{
public :
    RangeDistributionByMaterialName() = default;
    RangeDistributionByMaterialName( RangeDistributionByMaterialName const& ) = default;
    RangeDistributionByMaterialName( RangeDistributionByMaterialName && ) = default;

    //! init material range of elements
    void init( std::map<std::string, std::tuple<Range<MeshType,MESH_ELEMENTS>,Range<MeshType,MESH_ELEMENTS>> > const& rangeMeshElementsByMaterial )
        {
            for ( auto const& rangeMat : rangeMeshElementsByMaterial )
            {
                std::string const& matName = rangeMat.first;
                M_matNames.insert( matName );
                for ( auto const& eltWrap : std::get<0>( rangeMat.second ) )
                {
                    auto const& elt = unwrap_ref( eltWrap );
                    M_mapEltIdToMatName[elt.id()] = matName;
                }
            }
        }


    //! return map of range of faces by material for a given type
    std::map<std::string, Range<MeshType,MESH_FACES>> const& rangeMeshFacesByMaterial( std::string const& type ) const
        {
            auto itFindType = M_rangeMeshFacesByMaterial.find( type );
            CHECK( itFindType != M_rangeMeshFacesByMaterial.end() ) << "type does not find " << type;
            return itFindType->second;
        }

    //! update map of range of faces by material for a given type
    void update( std::string const& type, Range<MeshType,MESH_FACES> const& rangeFaces )
        {
            auto & rangeMatFaces = M_rangeMeshFacesByMaterial[type];
            for ( std::string const& matName : M_matNames )
            {
                if ( rangeMatFaces.find( matName ) == rangeMatFaces.end() )
                {
                    rangeMatFaces[matName] = Range<MeshType,MESH_FACES>(rangeFaces.mesh());
                }
            }

            for ( auto const& faceWrap : rangeFaces )
            {
                auto const& face = unwrap_ref( faceWrap );
                std::vector<size_type> eltIdsConnected;
                if ( face.isConnectedTo0() )
                    eltIdsConnected.push_back( face.element0().id() );
                if ( face.isConnectedTo1() )
                    eltIdsConnected.push_back( face.element1().id() );
                for ( size_type eltId : eltIdsConnected )
                {
                    auto itFindEltId = M_mapEltIdToMatName.find( eltId );
                    if ( itFindEltId == M_mapEltIdToMatName.end() )
                        continue;
                    std::string const& matName = itFindEltId->second;
                    CHECK( rangeMatFaces.find( matName ) != rangeMatFaces.end() ) << "invalid matName " << matName;
                    rangeMatFaces[matName].push_back( faceWrap );
                }
            }
        }
private :

    std::set<std::string> M_matNames;
    std::unordered_map<size_type,std::string> M_mapEltIdToMatName;
    std::map<std::string,std::map<std::string,Range<MeshType,MESH_FACES>>> M_rangeMeshFacesByMaterial;

};

} // namespace FeelModels
} // namespace Feel
