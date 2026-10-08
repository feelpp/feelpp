/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2026-10-02

 Copyright (C) 2026 Feel++ Consortium

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

#ifndef FEELPP_FEELDISCR_PRECICE_HPP
#define FEELPP_FEELDISCR_PRECICE_HPP 1

#if defined( FEELPP_HAS_PRECICE )

#include <boost/mp11/algorithm.hpp>
#include <precice/precice.hpp>
#include <feel/feeldiscr/pdhv.hpp>


namespace Feel
{

namespace concepts
{
template <typename T>
concept RangeOnElements = requires {
    { T::isOnElements() } -> std::convertible_to<bool>;
    requires T::isOnElements();
};

template <typename T>
concept RangeOnFaces = requires {
    { T::isOnFaces() } -> std::convertible_to<bool>;
    requires T::isOnFaces();
};
}

class PreciceData{
public:
    enum class ShapeType { scalar=0, vector };
    enum class ModeType { read=0, write };
    PreciceData( std::string const& name )
        :
        M_name( name )
        {}
    PreciceData( PreciceData const& ) = default;
    PreciceData( PreciceData&& ) = default;
    PreciceData& operator=( PreciceData const& ) = default;
    PreciceData& operator=( PreciceData&& ) = default;

    std::string const& name() const { return M_name; }
    virtual ShapeType shape() const = 0;
    virtual ModeType mode() const = 0;

private:
    std::string M_name;
};

template <typename SpaceType>
class PreciceReadWriteData : public PreciceData
{
public:
    using space_type = SpaceType;
    constexpr static ShapeType shape_c = space_type::is_scalar ? ShapeType::scalar : ShapeType::vector;

    PreciceReadWriteData( std::string const& name )
        :
        PreciceData( name )
        {}
    ShapeType shape() const override { return shape_c; }

    typename space_type::element_type const& field() const { return *M_field; }
    typename space_type::element_type & field() { return *M_field; }
    void setField( std::shared_ptr<typename space_type::element_type> field ) { M_field = field; }
    bool hasField() const { return M_field != nullptr; }

private:
    std::shared_ptr<typename space_type::element_type> M_field;
};
template <typename SpaceType>
class PreciceReadData : public PreciceReadWriteData<SpaceType>
{
    using super_type = PreciceReadWriteData<SpaceType>;
public:
    using self_type = PreciceReadData<SpaceType>;
    using space_type = typename super_type::space_type;
    using callback_dataupdated_type = std::function<void( self_type const& )>;

    PreciceReadData( std::string const& name, callback_dataupdated_type func )
        :
        super_type( name ),
        M_func( func )
        {}
    typename PreciceData::ModeType mode() const override { return PreciceData::ModeType::read; }

    callback_dataupdated_type const& callbackDataUpdated() const { return M_func; }
    template <typename FuncType>
    void setCallbackDataUpdated( FuncType && func ) { M_func = std::forward<FuncType>( func ); }

    void invokeDataUpdated()
        {
            if ( M_func )
                std::invoke( M_func, *this );
        }

private:
    std::shared_ptr<typename space_type::element_type> M_field;
    callback_dataupdated_type M_func;
};

template <typename SpaceType>
class PreciceWriteData : public PreciceReadWriteData<SpaceType>
{
    using super_type = PreciceReadWriteData<SpaceType>;
public:
    using self_type = PreciceWriteData<SpaceType>;
    using space_type = typename super_type::space_type;
    using callback_datarequested_type = std::function<void( self_type & )>;

    PreciceWriteData( std::string const& name, callback_datarequested_type func )
        :
        super_type( name ),
        M_func( func )
        {}
    typename PreciceData::ModeType mode() const override { return PreciceData::ModeType::write; }

    callback_datarequested_type const& callbackDataRequested() const { return M_func; }
    template <typename FuncType>
    void setCallbackDataRequested( FuncType && func ) { M_func = std::forward<FuncType>( func ); }

    void invokeDataRequested()
        {
            if ( M_func )
                std::invoke( M_func, *this );
        }

private:
    std::shared_ptr<typename space_type::element_type> M_field;
    callback_datarequested_type M_func;
};

enum class PreciceCouplingMeshLocation { vertices=0, barycenter };

template <typename MeshType,PreciceCouplingMeshLocation Location = PreciceCouplingMeshLocation::vertices>
class PreciceCouplingMesh
{
    enum class SpaceIndex { lagrange_scalar=0, lagrange_vector };

    template<class T>
    using get_functionspace_element_t = typename T::element_type;
public:
    using mesh_type = MeshType;
    static constexpr PreciceCouplingMeshLocation location_c = Location;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    using space_scalar_type = std::conditional_t<location_c==PreciceCouplingMeshLocation::vertices,
        Pch_type<mesh_type,1>, Pdh_type<mesh_type,0>>;
    using space_scalar_element_type = typename space_scalar_type::element_type;
    using space_vector_type = std::conditional_t<location_c==PreciceCouplingMeshLocation::vertices,
        Pchv_type<mesh_type,1>, Pdhv_type<mesh_type,0>>;
    using space_vector_element_type = typename space_vector_type::element_type;
    using spaces_type = std::tuple<space_scalar_type, space_vector_type>;
    using spaces_element_type = boost::mp11::mp_transform<get_functionspace_element_t,spaces_type>;

    template <PreciceData::ShapeType Shape>
    static constexpr SpaceIndex space_index_from_shape_v = Shape==PreciceData::ShapeType::scalar ? SpaceIndex::lagrange_scalar : SpaceIndex::lagrange_vector;
    template <SpaceIndex S>
    using space_type_t = std::tuple_element_t<std::to_underlying(S), spaces_type>;
    template <PreciceData::ShapeType Shape>
    using space_from_shape_t = space_type_t< space_index_from_shape_v<Shape> >;


    PreciceCouplingMesh( precice::Participant& precice, mesh_ptrtype mesh, std::string const& preciceMeshName )
        :
        M_precice( precice ),
        M_mesh( mesh ),
        M_meshName( preciceMeshName )
        {
            this->initMapping();
        }

    PreciceCouplingMesh( precice::Participant& precice, mesh_ptrtype mesh, Range<mesh_type,MESH_ELEMENTS> rangeElt, std::string const& preciceMeshName )
        :
        M_precice( precice ),
        M_mesh( mesh ),
        M_meshName( preciceMeshName )
        {
            this->initMapping();
        }

    template <SpaceIndex S>
    auto initSpace()
        {
            static constexpr int space_index = std::to_underlying(S);
            using _space_type = std::tuple_element_t<space_index, spaces_type>;
            auto space = std::get<space_index>( M_spaces );
            if ( !space )
                space = _space_type::New(_mesh=M_mesh );
            return space;
        }
    std::shared_ptr<space_scalar_type> initSpaceScalar() { return this->initSpace<SpaceIndex::lagrange_scalar>(); }
    std::shared_ptr<space_vector_type> initSpaceVector() { return this->initSpace<SpaceIndex::lagrange_vector>(); }




    template <PreciceData::ShapeType Shape,typename FuncType = typename PreciceWriteData<space_from_shape_t<Shape>>::callback_datarequested_type>
    PreciceWriteData< space_from_shape_t<Shape> > & initWriteData( std::string const& preciceDataName, FuncType && func = {} )
        {
            return this->initWriteData< space_index_from_shape_v<Shape>,FuncType >( preciceDataName, std::forward<FuncType>( func ) );
        }

    template <SpaceIndex S,typename FuncType = typename PreciceWriteData<space_type_t<S>>::callback_datarequested_type>
    PreciceWriteData<space_type_t<S>> & initWriteData( std::string const& preciceDataName, FuncType && func = {} )
        {
            static constexpr int space_index = std::to_underlying(S);
            using _space_type = std::tuple_element_t<space_index, spaces_type>;
            auto space = this->initSpace<S>();
            auto & writeDataMap = std::get<space_index>( M_writeData );

            auto itFindData = writeDataMap.find( preciceDataName );
            if ( itFindData == writeDataMap.end() )
                writeDataMap.emplace( preciceDataName, PreciceWriteData<_space_type>( preciceDataName, func ) );
            else if constexpr (requires { static_cast<bool>(func); }) {
                if ( func )
                    itFindData->second.setCallbackDataRequested( std::forward<FuncType>( func ) );
            }
            else
                itFindData->second.setCallbackDataRequested( std::forward<FuncType>( func ) );

            auto & writeData = writeDataMap.at( preciceDataName );
            // warning, we use shared internal fields, so if the user calls initWriteData
            // with different data names, the internal field will be shared between them.
            if ( !writeData.hasField() )
            {
                auto & internalField = std::get<space_index>( M_internalFields );
                if ( !internalField )
                    internalField = space->elementPtr();
                writeData.setField( internalField );
            }
            return writeData;
        }
    PreciceWriteData<space_type_t<SpaceIndex::lagrange_scalar>> & initWriteScalarData( std::string const& preciceDataName ) { return this->initWriteData<SpaceIndex::lagrange_scalar>( preciceDataName ); }
    PreciceWriteData<space_type_t<SpaceIndex::lagrange_scalar>> & initWriteVectorData( std::string const& preciceDataName ) { return this->initWriteData<SpaceIndex::lagrange_vector>( preciceDataName ); }

    template <SpaceIndex S>
    void hasWriteData( std::string const& preciceDataName )
        {
            static constexpr int space_index = std::to_underlying(S);
            auto & writeDataMap = std::get<space_index>( M_writeData );
            auto itFindData = writeDataMap.find( preciceDataName );
            return itFindData != writeDataMap.end();
        }

    template <SpaceIndex S,typename FuncType = typename PreciceWriteData<space_type_t<S>>::callback_datarequested_type>
    PreciceWriteData<space_type_t<S>> & getWriteData( std::string const& preciceDataName )
        {
            if ( !this->hasWriteData<S>( preciceDataName ) )
                throw std::runtime_error( "PreciceCouplingMesh::getWriteData : no write data with name " + preciceDataName + " has been initialized" );

            static constexpr int space_index = std::to_underlying(S);
            auto & writeDataMap = std::get<space_index>( M_writeData );
            return writeDataMap.at( preciceDataName );
        }

    template <SpaceIndex S>
    void writeData( std::string const& preciceDataName )
        {
            auto & _writeData = this->getWriteData<S>( preciceDataName );
            this->write( _writeData );
        }


    template <typename PreciceWriteDataType>
    void write( PreciceWriteDataType & _writeData )
        {
            static constexpr SpaceIndex spaceIndex = space_index_from_shape_v< std::decay_t<PreciceWriteDataType>::shape_c>;
            using _space_type = typename std::decay_t<PreciceWriteDataType>::space_type;
            std::string const& preciceDataName = _writeData.name();
            std::vector<double> preciceDataValues( M_preciceIDs.size() * _space_type::nComponents, 0.0 );

            // get write data and init the field if not already done
            auto & writeData = this->initWriteData<spaceIndex>( preciceDataName );
            auto & feelDataField = writeData.field();
            // invoke the callback to update the field values before writing to preCICE
            writeData.invokeDataRequested();

            for ( auto const& p_id : M_preciceIDs )
            {
                index_type f_scalardofid = M_preciceIdToScalarDofId[p_id];
                for ( uint16_type k = 0; k < _space_type::nComponents; ++k )
                {
                    index_type f_dofid = f_scalardofid*_space_type::nComponents + k;
                    preciceDataValues[ _space_type::nComponents * p_id + k ] = feelDataField( f_dofid );
                }
            }
            M_precice.writeData( M_meshName, preciceDataName, M_preciceIDs, preciceDataValues );
        }

    template <PreciceData::ShapeType Shape,typename FuncType = typename PreciceReadData<space_from_shape_t<Shape>>::callback_dataupdated_type>
    PreciceReadData< space_from_shape_t<Shape> > & initReadData( std::string const& preciceDataName, FuncType && func = {} )
        {
            return this->initReadData< space_index_from_shape_v<Shape>,FuncType >( preciceDataName, std::forward<FuncType>( func ) );
        }

    template <SpaceIndex S,typename FuncType = typename PreciceReadData<space_type_t<S>>::callback_dataupdated_type>
    PreciceReadData<space_type_t<S>> & initReadData( std::string const& preciceDataName,
                                                     FuncType && func = {} )
        {
            static constexpr int space_index = std::to_underlying(S);
            using _space_type = std::tuple_element_t<space_index, spaces_type>;
            auto space = this->initSpace<S>();
            auto & readDataMap = std::get<space_index>( M_readData );

            auto itFindData = readDataMap.find( preciceDataName );
            if ( itFindData == readDataMap.end() )
                readDataMap.emplace( preciceDataName, PreciceReadData<_space_type>( preciceDataName, func ) );
            else if constexpr (requires { static_cast<bool>(func); }) {
                if ( func )
                    itFindData->second.setCallbackDataUpdated( std::forward<FuncType>( func ) );
            }
            else
                itFindData->second.setCallbackDataUpdated( std::forward<FuncType>( func ) );

            auto & readData = readDataMap.at( preciceDataName );
            if ( !readData.hasField() )
                readData.setField( space->elementPtr() );
            return readData;
        }
    PreciceReadData<space_type_t<SpaceIndex::lagrange_scalar>> & initReadScalarData( std::string const& preciceDataName ) { return this->initReadData<SpaceIndex::lagrange_scalar>( preciceDataName ); }
    PreciceReadData<space_type_t<SpaceIndex::lagrange_scalar>> & initReadVectorData( std::string const& preciceDataName ) { return this->initReadData<SpaceIndex::lagrange_vector>( preciceDataName ); }

    template <SpaceIndex S>
    void hasReadData( std::string const& preciceDataName )
        {
            static constexpr int space_index = std::to_underlying(S);
            auto & readDataMap = std::get<space_index>( M_readData );
            auto itFindData = readDataMap.find( preciceDataName );
            return itFindData != readDataMap.end();
        }

    template <SpaceIndex S,typename FuncType = typename PreciceReadData<space_type_t<S>>::callback_dataupdated_type>
    PreciceReadData<space_type_t<S>> & getReadData( std::string const& preciceDataName )
        {
            if ( !this->hasReadData<S>( preciceDataName ) )
                throw std::runtime_error( "PreciceCouplingMesh::getReadData : no read data with name " + preciceDataName + " has been initialized" );

            static constexpr int space_index = std::to_underlying(S);
            auto & readDataMap = std::get<space_index>( M_readData );
            return readDataMap.at( preciceDataName );
        }

    template <SpaceIndex S>
    void readData( std::string const& preciceDataName )
        {
            auto & _readData = this->getReadData<S>( preciceDataName );
            this->read( _readData );
        }

    template <typename PreciceReadDataType>
    void read( PreciceReadDataType & _readData )
        {
            static constexpr SpaceIndex spaceIndex = space_index_from_shape_v< std::decay_t<PreciceReadDataType>::shape_c>;
            using _space_type = typename std::decay_t<PreciceReadDataType>::space_type;
            std::string const& preciceDataName = _readData.name();
            std::vector<double> preciceDataValues( M_preciceIDs.size() * _space_type::nComponents );
            M_precice.readData( M_meshName, preciceDataName, M_preciceIDs, 0, preciceDataValues );

            // get read data and init the field if not already done
            auto & readData = this->initReadData<spaceIndex>( preciceDataName );
            auto & feelDataField = readData.field();
            for ( auto const& p_id : M_preciceIDs )
            {
                for ( uint16_type k = 0; k < _space_type::nComponents; ++k )
                {
                    index_type f_dofid = M_preciceIdToScalarDofId[p_id]*_space_type::nComponents + k;
                    feelDataField.set( f_dofid, preciceDataValues[ _space_type::nComponents * p_id + k ] );
                }
            }

            readData.invokeDataUpdated();
        }
    void readScalarData( std::string const& preciceDataName ) { this->readData<SpaceIndex::lagrange_scalar>( preciceDataName ); }
    void readVectorData( std::string const& preciceDataName ) { this->readData<SpaceIndex::lagrange_vector>( preciceDataName ); }


    void readAllData()
        {
            auto applyReadData = [this]( auto & readDataMap )
                                     {
                                         for ( auto & [preciceDataName,_readData] : readDataMap )
                                             this->read( _readData );
                                     };
            std::apply([applyReadData]<typename... T>(T&&... args) {
                    ((applyReadData( args )), ...);
                }, M_readData);
        }

    void writeAllData()
        {
            auto applyWriteData = [this]( auto & writeDataMap )
                                      {
                                          for ( auto & [preciceDataName,_writeData] : writeDataMap )
                                              this->write( _writeData );
                                      };
            std::apply([applyWriteData]<typename... T>(T&&... args) {
                    ((applyWriteData( args )), ...);
                }, M_writeData);
        }
private:
    void initMapping()
        {
            if constexpr (location_c==PreciceCouplingMeshLocation::vertices)
                this->initMappingVertices();
            else if constexpr (location_c==PreciceCouplingMeshLocation::barycenter)
                this->initMappingBarycenter();
        }
    void initMappingVertices()
        {
            auto space = initSpaceScalar();
            auto mesh = space->mesh();
            auto dof = space->dof();
            int dim = mesh->realDimension();
            std::vector<double> coordinates;
            std::map<index_type,index_type> visitedVertices;// dof id to container rank of feelVertexIDs

            // tmp container to store the feel dof ids of the vertices that are sent to precice
            std::vector<index_type> feelVertexIDs;

            bool requiresMeshConnectivity = M_precice.requiresMeshConnectivityFor( M_meshName );

            std::vector<index_type> feelEdgeVertexIDs;
            std::vector<index_type> feelTriangleVertexIDs;
            std::vector<index_type> feelQuadrilateralVertexIDs;
            std::vector<index_type> feelTetrahedronVertexIDs;
            for ( auto const& eltWrap : elements(mesh) )
            {
                auto const& elt = unwrap_ref(eltWrap);
                std::vector<index_type> feelElementVertexIDs( elt.nPoints(), invalid_v<index_type> );
                for( auto const& ldof : dof->localDof( elt.id() ) )
                {
                    index_type thedof = ldof.second.index();
                    uint16_type localDof = ldof.first.localDof();
                    // with Pch1, the localDof can be greater than the number of points in the face
                    if ( localDof >= elt.nPoints() )
                        continue;

                    auto const& point = elt.point( localDof );

                    auto itFind = visitedVertices.find( thedof );
                    if ( itFind == visitedVertices.end() )
                    {
                        itFind = visitedVertices.insert( std::make_pair( thedof, feelVertexIDs.size() ) ).first;
                        feelVertexIDs.push_back( thedof );
                        for ( int d = 0; d < dim; ++d )
                            coordinates.push_back( point[d] );
                    }

                    if ( requiresMeshConnectivity )
                        feelElementVertexIDs[localDof] = itFind->second;
                }

                if ( requiresMeshConnectivity )
                {
                    if constexpr ( mesh_type::nDim == 1 )
                    {
                        feelEdgeVertexIDs.insert( feelEdgeVertexIDs.end(), feelElementVertexIDs.begin(), feelElementVertexIDs.end() );
                    }
                    else if constexpr ( mesh_type::nRealDim == 2 )
                    {
                        if ( elt.nPoints() == 3 )
                            feelTriangleVertexIDs.insert( feelTriangleVertexIDs.end(), feelElementVertexIDs.begin(), feelElementVertexIDs.end() );
                        else if ( elt.nPoints() == 4 )
                            feelQuadrilateralVertexIDs.insert( feelQuadrilateralVertexIDs.end(), feelElementVertexIDs.begin(), feelElementVertexIDs.end() );
                    }
                    else if constexpr ( mesh_type::nRealDim == 3 )
                    {
                        feelTetrahedronVertexIDs.insert( feelTetrahedronVertexIDs.end(), feelElementVertexIDs.begin(), feelElementVertexIDs.end() );
                    }
                }
            }

            std::vector<int> preciceVertexIDs( feelVertexIDs.size() );
            M_precice.setMeshVertices( M_meshName, coordinates, preciceVertexIDs );
            M_preciceIDs = preciceVertexIDs;
            for ( size_t i = 0; i < preciceVertexIDs.size(); ++i )
            {
                int p_id = preciceVertexIDs[i];
                index_type f_dofid = feelVertexIDs[i];
                M_preciceIdToScalarDofId[p_id] = f_dofid;
            }


            for ( std::size_t i = 0; i < feelEdgeVertexIDs.size(); i += 2 )
            {
                index_type k0 = feelEdgeVertexIDs[i];
                index_type k1 = feelEdgeVertexIDs[i+1];
                M_precice.setMeshEdge( M_meshName, preciceVertexIDs[k0], preciceVertexIDs[k1] );
            }

            for ( std::size_t i = 0; i < feelTriangleVertexIDs.size(); i += 3 )
            {
                index_type k0 = feelTriangleVertexIDs[i];
                index_type k1 = feelTriangleVertexIDs[i+1];
                index_type k2 = feelTriangleVertexIDs[i+2];
                M_precice.setMeshTriangle( M_meshName, preciceVertexIDs[k0], preciceVertexIDs[k1], preciceVertexIDs[k2] );
            }
            for ( std::size_t i = 0; i < feelQuadrilateralVertexIDs.size(); i += 4 )
            {
                index_type k0 = feelQuadrilateralVertexIDs[i];
                index_type k1 = feelQuadrilateralVertexIDs[i+1];
                index_type k2 = feelQuadrilateralVertexIDs[i+2];
                index_type k3 = feelQuadrilateralVertexIDs[i+3];
                M_precice.setMeshQuad( M_meshName, preciceVertexIDs[k0], preciceVertexIDs[k1], preciceVertexIDs[k2], preciceVertexIDs[k3] );
            }
            for ( std::size_t i = 0; i < feelTetrahedronVertexIDs.size(); i += 4 )
            {
                index_type k0 = feelTetrahedronVertexIDs[i];
                index_type k1 = feelTetrahedronVertexIDs[i+1];
                index_type k2 = feelTetrahedronVertexIDs[i+2];
                index_type k3 = feelTetrahedronVertexIDs[i+3];
                M_precice.setMeshTetrahedron( M_meshName, preciceVertexIDs[k0], preciceVertexIDs[k1], preciceVertexIDs[k2], preciceVertexIDs[k3] );
            }
        }


        void initMappingBarycenter()
        {
            auto space = initSpaceScalar();
            auto mesh = space->mesh();
            auto dof = space->dof();
            int dim = mesh->realDimension();
            std::vector<double> coordinates;
            // tmp container to store the feel dof ids
            std::vector<index_type> feelDofIds;

            for ( auto const& eltWrap : elements(mesh) )
            {
                auto const& elt = unwrap_ref(eltWrap);

                auto bary = elt.barycenter();

                for ( int d = 0; d < dim; ++d )
                    coordinates.push_back( bary[d] );
                // with Pdh0, the number of local dofs is 1, so we can just take the first one
                for( auto const& ldof : dof->localDof( elt.id() ) )
                {
                    index_type thedof = ldof.second.index();
                    uint16_type localDof = ldof.first.localDof();
                    feelDofIds.push_back( thedof );
                    break;
                }
            }
            std::vector<int> preciceVertexIDs( feelDofIds.size() );
            M_precice.setMeshVertices( M_meshName, coordinates, preciceVertexIDs );
            M_preciceIDs = preciceVertexIDs;
            for ( size_t i = 0; i < preciceVertexIDs.size(); ++i )
            {
                int p_id = preciceVertexIDs[i];
                index_type f_dofid = feelDofIds[i];
                M_preciceIdToScalarDofId[p_id] = f_dofid;
            }
        }

private:
    precice::Participant& M_precice;
    mesh_ptrtype M_mesh;
    std::string M_meshName;
    std::vector<int> M_preciceIDs;
    std::map<int, index_type> M_preciceIdToScalarDofId;

    template<class T>
    using add_shared_ptr_t = std::shared_ptr<T>;
    boost::mp11::mp_transform<add_shared_ptr_t,spaces_type> M_spaces;

    template<class T>
    using add_map_precice_read_data_t = std::map<std::string, PreciceReadData<T> >;
    boost::mp11::mp_transform<add_map_precice_read_data_t,spaces_type> M_readData;

    template<class T>
    using add_map_precice_write_data_t = std::map<std::string, PreciceWriteData<T> >;
    boost::mp11::mp_transform<add_map_precice_write_data_t,spaces_type> M_writeData;

    boost::mp11::mp_transform<add_shared_ptr_t,spaces_element_type> M_internalFields;
};

template <typename MeshType>
class PreciceAdapter
{
public:
    using mesh_type = MeshType;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using trace_mesh_type = trace_mesh_t<mesh_type>;

    PreciceAdapter() = default;

    ~PreciceAdapter() {
        if ( M_precice )
            M_precice->finalize();
    };

    bool hasParticipant() const { return M_precice != nullptr; }
    precice::Participant& participant() { return *M_precice; }

    void initParticipant( std::string const& solverName, std::string const& configFileName,
                          WorldComm const& worldComm = Environment::worldComm() )
        {
            M_precice = std::make_unique< precice::Participant >( solverName, configFileName,
                                                                  /*comm_rank*/ worldComm.globalRank(),
                                                                  /*comm_size*/ worldComm.globalSize() );
        }

    void inspectPreciceConfiguration()
        {
            // TODO : inspect the precice configuration and print some information about the coupling
        }


    template <PreciceCouplingMeshLocation Location = PreciceCouplingMeshLocation::vertices, concepts::RangeOnElements RangeElts>
    auto initCouplingMesh( std::string const& preciceMeshName, RangeElts const& rangeElts )
        {
            using precice_coupling_mesh_type = PreciceCouplingMesh<mesh_type,Location>;
            M_couplingMesh[preciceMeshName] = std::make_unique<precice_coupling_mesh_type>( *M_precice, shared_from_this(rangeElts.mesh()), rangeElts, preciceMeshName );
            return std::get<std::unique_ptr<precice_coupling_mesh_type>>( M_couplingMesh.at( preciceMeshName ) ).get();
        }


    template <PreciceCouplingMeshLocation Location = PreciceCouplingMeshLocation::vertices, concepts::RangeOnFaces RangeFaces>
    auto initCouplingMesh( std::string const& preciceMeshName, RangeFaces const& rangeFaces )
        {
            using precice_coupling_mesh_type = PreciceCouplingMesh<trace_mesh_type,Location>;
            auto submesh = createSubmesh( _range=rangeFaces,_view=true );
            M_couplingMesh[preciceMeshName] = std::make_unique<precice_coupling_mesh_type>( *M_precice, submesh, preciceMeshName );
            return std::get<std::unique_ptr<precice_coupling_mesh_type>>( M_couplingMesh.at( preciceMeshName ) ).get();
        }

    void updateForUse()
        {
            M_precice->initialize();
            this->inspectPreciceConfiguration();
        }

    void readData()
        {
            for ( auto & [meshName,couplingMesh] : M_couplingMesh )
            {
                std::visit( [this]( auto & couplingMeshPtr ){
                    couplingMeshPtr->readAllData();
                }, couplingMesh );
            }
        }

    void writeData()
        {
            for ( auto & [meshName,couplingMesh] : M_couplingMesh )
            {
                std::visit( [this]( auto & couplingMeshPtr ){
                    couplingMeshPtr->writeAllData();
                }, couplingMesh );
            }
        }

    void advance( double dt )
        {
            //double dt = this->M_precice->getMaxTimeStepSize();
            this->M_precice->advance( dt );
        }

protected:
    std::unique_ptr<precice::Participant> M_precice;
    std::map<std::string,std::variant<std::unique_ptr<PreciceCouplingMesh<mesh_type,PreciceCouplingMeshLocation::vertices>>,
                                      std::unique_ptr<PreciceCouplingMesh<mesh_type,PreciceCouplingMeshLocation::barycenter>>,
                                      std::unique_ptr<PreciceCouplingMesh<trace_mesh_type,PreciceCouplingMeshLocation::vertices>>,
                                      std::unique_ptr<PreciceCouplingMesh<trace_mesh_type,PreciceCouplingMeshLocation::barycenter>>
                                      >> M_couplingMesh;
};

} // namespace Feel

#endif // FEELPP_HAS_PRECISE
#endif
