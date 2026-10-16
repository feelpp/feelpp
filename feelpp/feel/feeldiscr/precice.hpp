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

/**
 * @file precice.hpp
 * @brief Integration of the preCICE coupling framework within Feel++.
 * @author Vincent Chabannes
 * @date 2026-10-02
 */

#ifndef FEELPP_FEELDISCR_PRECICE_HPP
#define FEELPP_FEELDISCR_PRECICE_HPP 1

#if defined( FEELPP_HAS_PRECICE )

#include <boost/mp11/algorithm.hpp>
#include <precice/precice.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>


namespace Feel
{

namespace concepts
{

/**
 * @brief Concept checking if a range is defined on elements.
 * @tparam T Range type.
 */
template <typename T>
concept RangeOnElements = requires {
    { T::isOnElements() } -> std::convertible_to<bool>;
    requires T::isOnElements();
};

/**
 * @brief Concept checking if a range is defined on faces/trace.
 * @tparam T Range type.
 */
template <typename T>
concept RangeOnFaces = requires {
    { T::isOnFaces() } -> std::convertible_to<bool>;
    requires T::isOnFaces();
};

} // namespace concepts


/**
 * @enum PreciceCouplingMeshLocation
 * @brief Discretization location for preCICE coupling mesh vertices.
 */
enum class PreciceCouplingMeshLocation { vertices = 0, barycenter };

/**
 * @enum PreciceDataShapeType
 * @brief Shape type for preCICE coupling data (scalar or vector).
 */
enum class PreciceDataShapeType { scalar = 0, vector };


/**
 * @class PreciceSetupData
 * @brief Represents a preCICE coupling data setup.
 */
class PreciceSetupData
{
public:
    PreciceSetupData() = default;
    PreciceSetupData( PreciceSetupData const& ) = default;
    PreciceSetupData( PreciceSetupData&& ) = default;
    PreciceSetupData& operator=( PreciceSetupData const& )= default;
    PreciceSetupData& operator=( PreciceSetupData&& )= default;

    std::string const& name() const { return M_name; }
    PreciceDataShapeType shape() const { return M_shape; }
    nl::json const& properties() const { return M_properties; }

    bool is( PreciceDataShapeType shapeType ) const { return M_shape == shapeType; }
    bool isScalar() const { return M_shape == PreciceDataShapeType::scalar; }
    bool isVector() const { return M_shape == PreciceDataShapeType::vector; }

    static PreciceSetupData fromJson( nl::json const& j )
    {
        std::string dataName = j.at( "name" ).get<std::string>();
        std::string dataType = boost::algorithm::to_lower_copy( j.at( "shape" ).get<std::string>() );
        if ( dataType != "scalar" && dataType != "vector" )
            throw std::runtime_error( fmt::format( "PreciceSetupData: invalid data type: {}. Must be 'scalar' or 'vector'", dataType ) );
        nl::json dataProperties = j.at( "properties" );
        PreciceSetupData res;
        res.M_name = dataName;
        res.M_shape = ( dataType == "scalar" ) ? PreciceDataShapeType::scalar : PreciceDataShapeType::vector;
        res.M_properties = dataProperties;
        return res;
    }
private:
    std::string M_name;
    PreciceDataShapeType M_shape;
    nl::json M_properties;
};

/**
 * @class PreciceSetupCouplingMesh
 * @brief Represents a preCICE coupling mesh setup.
 */
class PreciceSetupCouplingMesh
{
public:
    PreciceSetupCouplingMesh() = default;
    PreciceSetupCouplingMesh( PreciceSetupCouplingMesh const& ) = default;
    PreciceSetupCouplingMesh( PreciceSetupCouplingMesh&& ) = default;
    PreciceSetupCouplingMesh& operator=( PreciceSetupCouplingMesh const& )= default;
    PreciceSetupCouplingMesh& operator=( PreciceSetupCouplingMesh&& )= default;

    std::string const& name() const { return M_name; }
    std::string const& markers() const { return M_markers; }
    PreciceCouplingMeshLocation locations() const { return M_locations; }
    bool isLocatedOnVertices() const { return M_locations == PreciceCouplingMeshLocation::vertices; }
    bool isLocatedOnBarycenter() const { return M_locations == PreciceCouplingMeshLocation::barycenter; }
    std::vector<std::string> const& readDataNames() const { return M_readDataNames; }
    std::vector<std::string> const& writeDataNames() const { return M_writeDataNames; }

    static PreciceSetupCouplingMesh fromJson( nl::json const& j )
    {
        std::string meshName = j.at( "name" ).get<std::string>();
        std::string markers = j.at( "markers" ).get<std::string>();
        std::string locations = boost::algorithm::to_lower_copy( j.at( "locations" ).get<std::string>() );
        if ( locations != "vertices" && locations != "barycenter" )
            throw std::runtime_error( fmt::format( "PreciceSetupCouplingMesh: invalid locations value: {}. Must be 'vertices' or 'barycenter'", locations ) );
        std::vector<std::string> readData = j.at( "readData" ).get<std::vector<std::string>>();
        std::vector<std::string> writeData = j.at( "writeData" ).get<std::vector<std::string>>();
        PreciceSetupCouplingMesh res;
        res.M_name = meshName;
        res.M_markers = markers;
        res.M_locations = ( locations == "vertices" ) ? PreciceCouplingMeshLocation::vertices : PreciceCouplingMeshLocation::barycenter;
        res.M_readDataNames = readData;
        res.M_writeDataNames = writeData;
        return res;
    }
private:
    std::string M_name;
    std::string M_markers;
    PreciceCouplingMeshLocation M_locations;
    std::vector<std::string> M_readDataNames;
    std::vector<std::string> M_writeDataNames;
};

/**
 * @class PreciceParticipantSetup
 * @brief Represents a preCICE participant setup, including meshes and data.
 */
class PreciceParticipantSetup
{
public:
    PreciceParticipantSetup( std::string const& participantName, std::string const& configFileName )
        :
        M_participantName( participantName ),
        M_configFileName( Environment::expand( configFileName ) )
    {
    }

    static PreciceParticipantSetup fromJson( nl::json const& j )
    {
        std::string participantName = j.at( "participant" ).get<std::string>();
        std::string configFileName = j.at( "configFileName" ).get<std::string>();
        PreciceParticipantSetup res( participantName, configFileName );
        for ( auto const& [meshName, meshData] : j.at( "couplingMeshes" ).items() )
        {
            PreciceSetupCouplingMesh couplingMesh = PreciceSetupCouplingMesh::fromJson( meshData );
            res.M_couplingMeshes[couplingMesh.name()] = couplingMesh;
        }
        for ( auto const& [dataName, dataData] : j.at( "data" ).items() )
        {
            PreciceSetupData couplingData = PreciceSetupData::fromJson( dataData );
            res.M_couplingData[couplingData.name()] = couplingData;
        }
        return res;
    }

    /**
     * @brief Get the participant name.
     * @return Reference to the participant name string.
     */
    std::string const& name() const { return M_participantName; }

    /**
     * @brief Get the configuration file name.
     * @return Reference to the configuration file name string.
     */
    std::string const& configFileName() const { return M_configFileName; }

    /**
     * @brief Get the coupling data map.
     * @return Const reference to the map of coupling data.
     */
    std::map<std::string, PreciceSetupData> const& couplingData() const { return M_couplingData; }

    /**
     * @brief Get a specific coupling data by name.
     * @param dataName Name of the coupling data.
     * @return Const reference to the corresponding PreciceSetupData.
     * @throws std::runtime_error if the data name is not found.
     */
    PreciceSetupData const& couplingData( std::string const& dataName ) const
    {
        auto it = M_couplingData.find( dataName );
        if ( it == M_couplingData.end() )
            throw std::runtime_error( "PreciceParticipantSetup: coupling data not found: " + dataName );
        return it->second;
    }

    /**
     * @brief Get the coupling meshes map.
     * @return Const reference to the map of coupling meshes.
     */
    std::map<std::string, PreciceSetupCouplingMesh> const& couplingMeshes() const { return M_couplingMeshes; }

    /**
     * @brief Get a specific coupling mesh by name.
     * @param meshName Name of the coupling mesh.
     * @return Const reference to the corresponding PreciceSetupCouplingMesh.
     * @throws std::runtime_error if the mesh name is not found.
     */
    PreciceSetupCouplingMesh const& couplingMesh( std::string const& meshName ) const
        {
            auto it = M_couplingMeshes.find( meshName );
            if ( it == M_couplingMeshes.end() )
                throw std::runtime_error( "PreciceParticipantSetup: coupling mesh not found: " + meshName );
            return it->second;
        }
private:
    std::string M_participantName;
    std::string M_configFileName;
    std::map<std::string, PreciceSetupData> M_couplingData;
    std::map<std::string, PreciceSetupCouplingMesh> M_couplingMeshes;
};

/**
 * @class PreciceData
 * @brief Abstract base class representing data exchanged via preCICE.
 */
class PreciceData
{
public:
    /**
     * @enum ShapeType
     * @brief Data shape classification (scalar or vector).
     */
    //enum class ShapeType { scalar = 0, vector };
    using ShapeType = PreciceDataShapeType;

    /**
     * @enum ModeType
     * @brief Data transfer mode (read from preCICE or write to preCICE).
     */
    enum class ModeType { read = 0, write };

    /**
     * @brief Constructor for PreciceData.
     * @param name Name of the preCICE data set.
     */
    explicit PreciceData( std::string const& name, nl::json const& properties = {} )
        : M_name( name ),
          M_properties( properties )
        {}

    PreciceData( PreciceData const& ) = default;
    PreciceData( PreciceData&& ) = default;
    PreciceData& operator=( PreciceData const& ) = default;
    PreciceData& operator=( PreciceData&& ) = default;
    virtual ~PreciceData() = default;

    /**
     * @brief Get the preCICE data name.
     * @return Reference to the data name string.
     */
    std::string const& name() const { return M_name; }

    /**
     * @brief Get the shape of the data.
     * @return ShapeType (scalar or vector).
     */
    virtual ShapeType shape() const = 0;

    /**
     * @brief Get the transfer mode of the data.
     * @return ModeType (read or write).
     */
    virtual ModeType mode() const = 0;

    /**
     * @brief Get the properties associated with the data.
     * @return Const reference to a JSON object containing properties.
     */
    nl::json const& properties() const { return M_properties; }

private:
    std::string M_name;
    nl::json M_properties;
};

/**
 * @class PreciceReadWriteData
 * @brief Base wrapper linking preCICE data with a Feel++ FunctionSpace element.
 * @tparam SpaceType Feel++ FunctionSpace type.
 */
template <typename SpaceType>
class PreciceReadWriteData : public PreciceData
{
public:
    using space_type = SpaceType;
    constexpr static ShapeType shape_c = space_type::is_scalar ? ShapeType::scalar : ShapeType::vector;

    /**
     * @brief Constructor for PreciceReadWriteData.
     * @param name Name of the data set.
     */
    explicit PreciceReadWriteData( std::string const& name, nl::json const& properties = {} )
        : PreciceData( name, properties )
        {}

    /**
     * @brief Get the data shape type.
     * @return ShapeType corresponding to the space dimension.
     */
    ShapeType shape() const override { return shape_c; }

    /**
     * @brief Get const reference to the underlying Feel++ field element.
     * @return Const reference to the field element.
     */
    typename space_type::element_type const& field() const { return *M_field; }

    /**
     * @brief Get mutable reference to the underlying Feel++ field element.
     * @return Mutable reference to the field element.
     */
    typename space_type::element_type & field() { return *M_field; }

    /**
     * @brief Set the Feel++ field element pointer.
     * @param field Shared pointer to the function space element.
     */
    void setField( std::shared_ptr<typename space_type::element_type> field ) { M_field = field; }

    /**
     * @brief Check if a field element is attached.
     * @return True if M_field is non-null, false otherwise.
     */
    bool hasField() const { return M_field != nullptr; }

private:
    std::shared_ptr<typename space_type::element_type> M_field;
};

/**
 * @class PreciceReadData
 * @brief Class handling data read from preCICE into Feel++.
 * @tparam SpaceType Feel++ FunctionSpace type.
 */
template <typename SpaceType>
class PreciceReadData : public PreciceReadWriteData<SpaceType>
{
    using super_type = PreciceReadWriteData<SpaceType>;
public:
    using self_type = PreciceReadData<SpaceType>;
    using space_type = typename super_type::space_type;
    using callback_dataupdated_type = std::function<void( self_type const& )>;

    /**
     * @brief Constructor for PreciceReadData.
     * @param name Name of the preCICE data.
     * @param properties Optional JSON object containing additional properties for the data.
     * @param func Optional callback invoked after data is read and updated.
     */
    PreciceReadData( std::string const& name, nl::json const& properties = {}, callback_dataupdated_type func = {} )
        : super_type( name, properties ),
          M_func( func )
        {}

    /**
     * @brief Get the transfer mode type.
     * @return ModeType::read.
     */
    typename PreciceData::ModeType mode() const override { return PreciceData::ModeType::read; }

    /**
     * @brief Get the callback function for data updates.
     * @return Const reference to callback_dataupdated_type.
     */
    callback_dataupdated_type const& callbackDataUpdated() const { return M_func; }

    /**
     * @brief Set or update the data-updated callback.
     * @tparam FuncType Callable type.
     * @param func Callback function.
     */
    template <typename FuncType>
    void setCallbackDataUpdated( FuncType && func ) { M_func = std::forward<FuncType>( func ); }

    /**
     * @brief Trigger the data-updated callback if defined.
     */
    void invokeDataUpdated()
        {
            if ( M_func )
                std::invoke( M_func, *this );
        }

private:
    callback_dataupdated_type M_func;
};

/**
 * @class PreciceWriteData
 * @brief Class handling data written from Feel++ to preCICE.
 * @tparam SpaceType Feel++ FunctionSpace type.
 */
template <typename SpaceType>
class PreciceWriteData : public PreciceReadWriteData<SpaceType>
{
    using super_type = PreciceReadWriteData<SpaceType>;
public:
    using self_type = PreciceWriteData<SpaceType>;
    using space_type = typename super_type::space_type;
    using callback_datarequested_type = std::function<void( self_type & )>;

    /**
     * @brief Constructor for PreciceWriteData.
     * @param name Name of the preCICE data.
     * @param properties Optional JSON object containing additional properties for the data.
     * @param func Optional callback invoked when data is requested before writing.
     */
    PreciceWriteData( std::string const& name, nl::json const& properties = {}, callback_datarequested_type func = {} )
        : super_type( name, properties ),
          M_func( func )
        {}

    /**
     * @brief Get the transfer mode type.
     * @return ModeType::write.
     */
    typename PreciceData::ModeType mode() const override { return PreciceData::ModeType::write; }

    /**
     * @brief Get the callback function for data requests.
     * @return Const reference to callback_datarequested_type.
     */
    callback_datarequested_type const& callbackDataRequested() const { return M_func; }

    /**
     * @brief Set or update the data-requested callback.
     * @tparam FuncType Callable type.
     * @param func Callback function.
     */
    template <typename FuncType>
    void setCallbackDataRequested( FuncType && func ) { M_func = std::forward<FuncType>( func ); }

    /**
     * @brief Trigger the data-requested callback if defined.
     */
    void invokeDataRequested()
        {
            if ( M_func )
                std::invoke( M_func, *this );
        }

private:
    callback_datarequested_type M_func;
};


/**
 * @class PreciceCouplingMesh
 * @brief Manages a preCICE mesh representation mapped to Feel++ meshes and function spaces.
 * @tparam MeshType Feel++ mesh type.
 * @tparam Location Coupling mesh location strategy (vertices P1/P1v or barycenter P0d/P0dv).
 */
template <typename RangeMeshEntitiesType, PreciceCouplingMeshLocation Location = PreciceCouplingMeshLocation::vertices>
class PreciceCouplingMesh
{
    enum class SpaceIndex { lagrange_scalar = 0, lagrange_vector };

    template<class T>
    using get_functionspace_element_t = typename T::element_type;
public:
    using range_mesh_entities_type = RangeMeshEntitiesType;
    using parent_mesh_type = typename range_mesh_entities_type::mesh_t;
    using mesh_type = std::conditional_t< range_mesh_entities_type::mesh_entities == MESH_ELEMENTS, parent_mesh_type, trace_mesh_t<parent_mesh_type> >;

    static constexpr PreciceCouplingMeshLocation location_c = Location;
    static constexpr bool isOnVertices_c = location_c == PreciceCouplingMeshLocation::vertices;
    static constexpr bool isOnBarycenter_c = location_c == PreciceCouplingMeshLocation::barycenter;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    using space_scalar_type = std::conditional_t<location_c == PreciceCouplingMeshLocation::vertices,
                                                 Pch_type<mesh_type, 1>, Pdh_type<mesh_type, 0>>;
    using space_scalar_element_type = typename space_scalar_type::element_type;
    using space_vector_type = std::conditional_t<location_c == PreciceCouplingMeshLocation::vertices,
                                                 Pchv_type<mesh_type, 1>, Pdhv_type<mesh_type, 0>>;
    using space_vector_element_type = typename space_vector_type::element_type;
    using spaces_type = std::tuple<space_scalar_type, space_vector_type>;
    using spaces_element_type = boost::mp11::mp_transform<get_functionspace_element_t, spaces_type>;

    template <PreciceData::ShapeType Shape>
    static constexpr SpaceIndex space_index_from_shape_v = Shape == PreciceData::ShapeType::scalar ? SpaceIndex::lagrange_scalar : SpaceIndex::lagrange_vector;
    template <SpaceIndex S>
    using space_type_t = std::tuple_element_t<std::to_underlying(S), spaces_type>;
    template <PreciceData::ShapeType Shape>
    using space_from_shape_t = space_type_t< space_index_from_shape_v<Shape> >;

    /**
     * @brief Construct coupling mesh on a range of mesh entities (elements, faces, edges, or points).
     * @param precice Reference to the preCICE participant.
     * @param rangeEntities Range of mesh entities (elements, faces, edges, or points).
     * @param preciceMeshName Name of preCICE mesh as specified in precice-config.xml.
     */
    PreciceCouplingMesh( precice::Participant& precice, range_mesh_entities_type const& rangeEntities, std::string const& preciceMeshName )
        : M_precice( precice ),
          M_meshName( preciceMeshName ),
          M_rangeMeshEntities( rangeEntities )
        {
            if constexpr ( range_mesh_entities_type::isOnElements() )
            {
                M_mesh = const_cast<mesh_type*>( rangeEntities.mesh() )->shared_from_this();
                this->initMapping( M_rangeMeshEntities );
            }
            else
            {
                M_mesh = createSubmesh( _range = rangeEntities, _view = true );
                this->initMapping( elements( M_mesh ) );
            }
        }

    /**
     * @brief Get the preCICE participant reference.
     * @return Reference to the preCICE participant.
     */
    precice::Participant& precice() { return M_precice; }

    /**
     * @brief Get the name of the preCICE mesh.
     * @return Const reference to the mesh name string.
     */
    std::string const& name() const { return M_meshName; }

    /**
     * @brief Get the shared pointer to the Feel++ mesh associated with this coupling mesh.
     * @return Const reference to the shared pointer of the mesh.
     */
    mesh_ptrtype const& couplingMesh() const { return M_mesh; }

    /**
     * @brief Get the range of mesh entities associated with this coupling mesh.
     * @return Const reference to the range of mesh entities.
     */
    range_mesh_entities_type const& rangeMeshEntities() const { return M_rangeMeshEntities; }


    /**
     * @brief Initialize or retrieve a function space by SpaceIndex.
     * @tparam S SpaceIndex (lagrange_scalar or lagrange_vector).
     * @return Shared pointer to the requested FunctionSpace.
     */
    template <SpaceIndex S>
    auto initSpace()
        {
            static constexpr int space_index = std::to_underlying(S);
            using _space_type = std::tuple_element_t<space_index, spaces_type>;
            auto & space = std::get<space_index>( M_spaces );
            if ( !space )
            {
                if constexpr ( range_mesh_entities_type::isOnElements() )
                    space = _space_type::New( _mesh = M_mesh, _range = M_rangeMeshEntities );
                else
                    space = _space_type::New( _mesh = M_mesh );
            }
            return space;
        }

    /** @brief Initialize scalar function space. */
    std::shared_ptr<space_scalar_type> initSpaceScalar() { return this->initSpace<SpaceIndex::lagrange_scalar>(); }
    /** @brief Initialize vector function space. */
    std::shared_ptr<space_vector_type> initSpaceVector() { return this->initSpace<SpaceIndex::lagrange_vector>(); }

    /**
     * @brief Retrieve a function space by SpaceIndex.
     * @tparam S SpaceIndex (lagrange_scalar or lagrange_vector).
     * @return Shared pointer to the requested FunctionSpace.
     * @throws std::runtime_error if the space has not been initialized yet.
     */
    template <SpaceIndex S>
    auto getSpace() const
        {
            static constexpr int space_index = std::to_underlying(S);
            auto space = std::get<space_index>( M_spaces );
            if ( !space )
                throw std::runtime_error( "PreciceCouplingMesh::getSpace : space has not been initialized yet" );
            return space;
        }

    /**
     * @brief Retrieve a function space by ShapeType (scalar or vector).
     * @tparam Shape Data shape type.
     * @return Shared pointer to the requested FunctionSpace.
     * @throws std::runtime_error if the space has not been initialized yet.
     */
    template <PreciceData::ShapeType Shape>
    auto getSpace() const { return this->getSpace< space_index_from_shape_v<Shape> >(); }

    /** @brief Retrieve the scalar function space. */
    auto getSpaceScalar() const { return this->getSpace<SpaceIndex::lagrange_scalar>(); }
    /** @brief Retrieve the vector function space. */
    auto getSpaceVector() const { return this->getSpace<SpaceIndex::lagrange_vector>(); }

    /**
     * @brief Initialize write data set by ShapeType.
     * @tparam Shape Data shape (scalar or vector).
     * @tparam FuncType Callback type.
     * @param preciceDataName Name of data in preCICE.
     * @param properties Optional JSON object containing additional properties for the data.
     * @param func Optional callback triggered when write data is requested.
     * @return Reference to initialized PreciceWriteData.
     */
    template <PreciceData::ShapeType Shape, typename FuncType = typename PreciceWriteData<space_from_shape_t<Shape>>::callback_datarequested_type>
    PreciceWriteData< space_from_shape_t<Shape> > & initWriteData( std::string const& preciceDataName, nl::json const& properties = {}, FuncType && func = {} )
        {
            return this->initWriteData< space_index_from_shape_v<Shape>, FuncType >( preciceDataName, properties, std::forward<FuncType>( func ) );
        }

    /**
     * @brief Initialize write data set by SpaceIndex.
     * @tparam S SpaceIndex.
     * @tparam FuncType Callback type.
     * @param preciceDataName Name of data in preCICE.
     * @param properties Optional JSON object containing additional properties for the data.
     * @param func Optional callback triggered when write data is requested.
     * @return Reference to initialized PreciceWriteData.
     */
    template <SpaceIndex S, typename FuncType = typename PreciceWriteData<space_type_t<S>>::callback_datarequested_type>
    PreciceWriteData<space_type_t<S>> & initWriteData( std::string const& preciceDataName, nl::json const& properties = {}, FuncType && func = {} )
        {
            static constexpr int space_index = std::to_underlying(S);
            using _space_type = std::tuple_element_t<space_index, spaces_type>;
            auto space = this->initSpace<S>();
            auto & writeDataMap = std::get<space_index>( M_writeData );

            auto itFindData = writeDataMap.find( preciceDataName );
            if ( itFindData == writeDataMap.end() )
                writeDataMap.emplace( preciceDataName, PreciceWriteData<_space_type>( preciceDataName, properties, func ) );
            else if constexpr ( requires { static_cast<bool>(func); } ) {
                if ( func )
                    itFindData->second.setCallbackDataRequested( std::forward<FuncType>( func ) );
            }
            else
                itFindData->second.setCallbackDataRequested( std::forward<FuncType>( func ) );

            auto & writeData = writeDataMap.at( preciceDataName );
            if ( !writeData.hasField() )
            {
                auto & internalField = std::get<space_index>( M_internalFields );
                if ( !internalField )
                    internalField = space->elementPtr();
                writeData.setField( internalField );
            }
            return writeData;
        }

    /**
     * @brief Check if a write data set has been initialized.
     * @tparam S SpaceIndex.
     * @param preciceDataName Name of data in preCICE.
     * @return True if initialized, false otherwise.
     */
    template <SpaceIndex S>
    bool hasWriteData( std::string const& preciceDataName )
        {
            static constexpr int space_index = std::to_underlying(S);
            auto & writeDataMap = std::get<space_index>( M_writeData );
            auto itFindData = writeDataMap.find( preciceDataName );
            return itFindData != writeDataMap.end();
        }

    /**
     * @brief Retrieve initialized write data object.
     * @tparam S SpaceIndex.
     * @param preciceDataName Name of data in preCICE.
     * @return Reference to PreciceWriteData.
     */
    template <SpaceIndex S>
    PreciceWriteData<space_type_t<S>> & getWriteData( std::string const& preciceDataName )
        {
            if ( !this->hasWriteData<S>( preciceDataName ) )
                throw std::runtime_error( "PreciceCouplingMesh::getWriteData : no write data with name " + preciceDataName + " has been initialized" );

            static constexpr int space_index = std::to_underlying(S);
            auto & writeDataMap = std::get<space_index>( M_writeData );
            return writeDataMap.at( preciceDataName );
        }

    /**
     * @brief Write data to preCICE by name.
     * @tparam S SpaceIndex.
     * @param preciceDataName Data name in preCICE.
     */
    template <SpaceIndex S>
    void writeData( std::string const& preciceDataName )
        {
            auto & _writeData = this->getWriteData<S>( preciceDataName );
            this->write( _writeData );
        }

    /**
     * @brief Perform actual write buffer formatting and dispatch to preCICE API.
     * @tparam PreciceWriteDataType PreciceWriteData specialization.
     * @param _writeData Write data instance.
     */
    template <typename PreciceWriteDataType>
    void write( PreciceWriteDataType & _writeData )
        {
            static constexpr SpaceIndex spaceIndex = space_index_from_shape_v< std::decay_t<PreciceWriteDataType>::shape_c >;
            using _space_type = typename std::decay_t<PreciceWriteDataType>::space_type;
            std::string const& preciceDataName = _writeData.name();
            std::vector<double> preciceDataValues( M_preciceIDs.size() * _space_type::nComponents, 0.0 );

            auto & writeData = this->initWriteData<spaceIndex>( preciceDataName );
            auto & feelDataField = writeData.field();
            writeData.invokeDataRequested();

            for ( auto const& p_id : M_preciceIDs )
            {
                index_type f_scalardofid = M_preciceIdToScalarDofId[p_id];
                for ( uint16_type k = 0; k < _space_type::nComponents; ++k )
                {
                    index_type f_dofid = f_scalardofid * _space_type::nComponents + k;
                    preciceDataValues[ _space_type::nComponents * p_id + k ] = feelDataField( f_dofid );
                }
            }
            M_precice.writeData( M_meshName, preciceDataName, M_preciceIDs, preciceDataValues );
        }

    /**
     * @brief Initialize read data set by ShapeType.
     * @tparam Shape Data shape.
     * @tparam FuncType Callback type.
     * @param preciceDataName Name of data in preCICE.
     * @param properties Optional JSON object containing additional properties for the data.
     * @param func Optional callback triggered after data update.
     * @return Reference to initialized PreciceReadData.
     */
    template <PreciceData::ShapeType Shape, typename FuncType = typename PreciceReadData<space_from_shape_t<Shape>>::callback_dataupdated_type>
    PreciceReadData< space_from_shape_t<Shape> > & initReadData( std::string const& preciceDataName, nl::json const& properties = {}, FuncType && func = {} )
        {
            return this->initReadData< space_index_from_shape_v<Shape>, FuncType >( preciceDataName, properties, std::forward<FuncType>( func ) );
        }

    /**
     * @brief Initialize read data set by SpaceIndex.
     * @tparam S SpaceIndex.
     * @tparam FuncType Callback type.
     * @param preciceDataName Name of data in preCICE.
     * @param properties Optional JSON object containing additional properties for the data.
     * @param func Optional callback triggered after data update.
     * @return Reference to initialized PreciceReadData.
     */
    template <SpaceIndex S, typename FuncType = typename PreciceReadData<space_type_t<S>>::callback_dataupdated_type>
    PreciceReadData<space_type_t<S>> & initReadData( std::string const& preciceDataName, nl::json const& properties = {}, FuncType && func = {} )
        {
            static constexpr int space_index = std::to_underlying(S);
            using _space_type = std::tuple_element_t<space_index, spaces_type>;
            auto space = this->initSpace<S>();
            auto & readDataMap = std::get<space_index>( M_readData );

            auto itFindData = readDataMap.find( preciceDataName );
            if ( itFindData == readDataMap.end() )
                readDataMap.emplace( preciceDataName, PreciceReadData<_space_type>( preciceDataName, properties, func ) );
            else if constexpr ( requires { static_cast<bool>(func); } ) {
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

    /**
     * @brief Check if a read data set has been initialized.
     * @tparam S SpaceIndex.
     * @param preciceDataName Name of data in preCICE.
     * @return True if initialized, false otherwise.
     */
    template <SpaceIndex S>
    bool hasReadData( std::string const& preciceDataName )
        {
            static constexpr int space_index = std::to_underlying(S);
            auto & readDataMap = std::get<space_index>( M_readData );
            auto itFindData = readDataMap.find( preciceDataName );
            return itFindData != readDataMap.end();
        }

    /**
     * @brief Retrieve initialized read data object.
     * @tparam S SpaceIndex.
     * @param preciceDataName Name of data in preCICE.
     * @return Reference to PreciceReadData.
     */
    template <SpaceIndex S>
    PreciceReadData<space_type_t<S>> & getReadData( std::string const& preciceDataName )
        {
            if ( !this->hasReadData<S>( preciceDataName ) )
                throw std::runtime_error( "PreciceCouplingMesh::getReadData : no read data with name " + preciceDataName + " has been initialized" );

            static constexpr int space_index = std::to_underlying(S);
            auto & readDataMap = std::get<space_index>( M_readData );
            return readDataMap.at( preciceDataName );
        }

    /**
     * @brief Read data from preCICE by name.
     * @tparam S SpaceIndex.
     * @param preciceDataName Data name in preCICE.
     */
    template <SpaceIndex S>
    void readData( std::string const& preciceDataName )
        {
            auto & _readData = this->getReadData<S>( preciceDataName );
            this->read( _readData );
        }

    /**
     * @brief Perform actual read operation from preCICE API and map into Feel++ field.
     * @tparam PreciceReadDataType PreciceReadData specialization.
     * @param _readData Read data instance.
     */
    template <typename PreciceReadDataType>
    void read( PreciceReadDataType & _readData )
        {
            static constexpr SpaceIndex spaceIndex = space_index_from_shape_v< std::decay_t<PreciceReadDataType>::shape_c >;
            using _space_type = typename std::decay_t<PreciceReadDataType>::space_type;
            std::string const& preciceDataName = _readData.name();
            std::vector<double> preciceDataValues( M_preciceIDs.size() * _space_type::nComponents );
            M_precice.readData( M_meshName, preciceDataName, M_preciceIDs, 0, preciceDataValues );

            auto & readData = this->initReadData<spaceIndex>( preciceDataName );
            auto & feelDataField = readData.field();
            for ( auto const& p_id : M_preciceIDs )
            {
                for ( uint16_type k = 0; k < _space_type::nComponents; ++k )
                {
                    index_type f_dofid = M_preciceIdToScalarDofId[p_id] * _space_type::nComponents + k;
                    feelDataField.set( f_dofid, preciceDataValues[ _space_type::nComponents * p_id + k ] );
                }
            }
            sync( feelDataField );

            readData.invokeDataUpdated();
        }

    /** @brief Read scalar data set from preCICE. */
    void readScalarData( std::string const& preciceDataName ) { this->readData<SpaceIndex::lagrange_scalar>( preciceDataName ); }
    /** @brief Read vector data set from preCICE. */
    void readVectorData( std::string const& preciceDataName ) { this->readData<SpaceIndex::lagrange_vector>( preciceDataName ); }

    /** @brief Read all registered data sets from preCICE for this coupling mesh. */
    void readAllData()
        {
            auto applyReadData = [this]( auto & readDataMap )
                                     {
                                         for ( auto & [preciceDataName, _readData] : readDataMap )
                                             this->read( _readData );
                                     };
            std::apply([applyReadData]<typename... T>(T&&... args) {
                    ((applyReadData( args )), ...);
                }, M_readData);
        }

    /** @brief Write all registered data sets to preCICE for this coupling mesh. */
    void writeAllData()
        {
            auto applyWriteData = [this]( auto & writeDataMap )
                                      {
                                          for ( auto & [preciceDataName, _writeData] : writeDataMap )
                                              this->write( _writeData );
                                      };
            std::apply([applyWriteData]<typename... T>(T&&... args) {
                    ((applyWriteData( args )), ...);
                }, M_writeData);
        }



    /** @brief Update Feel++ field data from an expression by using either nodal or L2 projection. */
    template <typename DataType,typename ExprType>
    void updateDataField( DataType & data, ExprType const& expr, bool useNodalProjection = false ) const
        {
            auto rangeEntities = this->rangeMeshEntities();
            if ( useNodalProjection )
            {
                data.field().on( _range=rangeEntities, _expr=expr,_close=true );
            }
            else // L2 projection
            {
                if constexpr ( DataType::space_type::is_continuous )
                {
                    throw std::runtime_error( "PreciceCouplingMesh::updateDataField : TODO: L2 projection is not supported for continuous function spaces" );
                }
                else
                {
                    auto P0d = this->getSpace<DataType::shape_c>();
                    auto F = Feel::backend()->newVector( _test=P0d );
                    auto uProj = P0d->element( F );
                    if constexpr ( range_mesh_entities_type::isOnElements() )
                    {
                        form1( _test=P0d, _vector=F )
                            = integrate( _range=rangeEntities,
                                         _expr=inner(expr,id(uProj))/meas() );
                    }
                    else if constexpr ( range_mesh_entities_type::isOnFaces() )
                    {
                        form1( _test=P0d, _vector=F )
                            = integrate( _range=rangeEntities,
                                         _expr=inner(expr,id(uProj))/measFace() );
                    }
                    else
                    {
                        throw std::runtime_error( "PreciceCouplingMesh::updateDataField : unsupported range of mesh entities" );
                    }
                    F->close();
                    data.field() = uProj;
                }
            }
        }

private:
    /** @brief Initialize preCICE mesh mapping based on Location template parameter and provided range of mesh elements. */
    void initMapping( Range<mesh_type, MESH_ELEMENTS> const& rangeElts )
        {
            if constexpr ( location_c == PreciceCouplingMeshLocation::vertices )
                this->initMappingVertices( rangeElts );
            else if constexpr ( location_c == PreciceCouplingMeshLocation::barycenter )
                this->initMappingBarycenter( rangeElts );
        }

    /** @brief Register nodal mesh vertices and element connectivities in preCICE. */
    void initMappingVertices( Range<mesh_type, MESH_ELEMENTS> const& rangeElts )
        {
            auto space = initSpaceScalar();
            auto mesh = space->mesh();
            auto dof = space->dof();
            int dim = mesh->realDimension();
            std::vector<double> coordinates;
            std::map<index_type, index_type> visitedVertices;

            std::vector<index_type> feelVertexIDs;
            bool requiresMeshConnectivity = M_precice.requiresMeshConnectivityFor( M_meshName );

            std::vector<index_type> feelEdgeVertexIDs;
            std::vector<index_type> feelTriangleVertexIDs;
            std::vector<index_type> feelQuadrilateralVertexIDs;
            std::vector<index_type> feelTetrahedronVertexIDs;

            for ( auto const& eltWrap : rangeElts )
            {
                auto const& elt = unwrap_ref(eltWrap);
                std::vector<index_type> feelElementVertexIDs( elt.nPoints(), invalid_v<index_type> );
                for ( auto const& ldof : dof->localDof( elt.id() ) )
                {
                    index_type thedof = ldof.second.index();
                    uint16_type localDof = ldof.first.localDof();
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

    /** @brief Register element barycenters (P0d) in preCICE. */
    void initMappingBarycenter( Range<mesh_type, MESH_ELEMENTS> const& rangeElts )
        {
            auto space = initSpaceScalar();
            auto mesh = space->mesh();
            auto dof = space->dof();
            int dim = mesh->realDimension();
            std::vector<double> coordinates;
            std::vector<index_type> feelDofIds;

            for ( auto const& eltWrap : rangeElts )
            {
                auto const& elt = unwrap_ref(eltWrap);
                auto bary = elt.barycenter();

                for ( int d = 0; d < dim; ++d )
                    coordinates.push_back( bary[d] );

                for ( auto const& ldof : dof->localDof( elt.id() ) )
                {
                    index_type thedof = ldof.second.index();
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
    std::string M_meshName;
    mesh_ptrtype M_mesh;
    range_mesh_entities_type M_rangeMeshEntities;
    std::vector<int> M_preciceIDs;
    std::map<int, index_type> M_preciceIdToScalarDofId;

    template<class T>
    using add_shared_ptr_t = std::shared_ptr<T>;
    boost::mp11::mp_transform<add_shared_ptr_t, spaces_type> M_spaces;

    template<class T>
    using add_map_precice_read_data_t = std::map<std::string, PreciceReadData<T> >;
    boost::mp11::mp_transform<add_map_precice_read_data_t, spaces_type> M_readData;

    template<class T>
    using add_map_precice_write_data_t = std::map<std::string, PreciceWriteData<T> >;
    boost::mp11::mp_transform<add_map_precice_write_data_t, spaces_type> M_writeData;

    boost::mp11::mp_transform<add_shared_ptr_t, spaces_element_type> M_internalFields;
};

/**
 * @class PreciceAdapter
 * @brief High-level adapter managing preCICE participant lifecycle and coupling meshes.
 * @tparam MeshType Primary Feel++ mesh type.
 */
template <typename MeshType>
class PreciceAdapter
{
public:
    using mesh_type = MeshType;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using trace_mesh_type = trace_mesh_t<mesh_type>;
    using variant_mesh_range_type = std::variant<Range<mesh_type,MESH_ELEMENTS>, Range<mesh_type,MESH_FACES> >;

    PreciceAdapter() = default;

    /**
     * @brief Destructor. Calls preCICE finalize if the participant exists.
     */
    ~PreciceAdapter()
        {
            if ( M_precice )
                M_precice->finalize();
        }

    /**
     * @brief Check if preCICE participant is initialized.
     * @return True if participant exists.
     */
    bool hasParticipant() const { return M_precice != nullptr; }

    /**
     * @brief Get reference to the underlying preCICE participant instance.
     * @return Reference to precice::Participant.
     */
    precice::Participant& participant() { return *M_precice; }

    /**
     * @brief Initialize preCICE participant instance.
     * @param participantName Participant name in precice-config.xml.
     * @param configFileName Path to precice-config.xml.
     * @param worldComm MPI communicator context (defaults to Feel++ WorldComm).
     */
    void initParticipant( std::string const& participantName, std::string const& configFileName,
                          WorldComm const& worldComm = Environment::worldComm() )
        {
            M_precice = std::make_unique< precice::Participant >( participantName, configFileName,
                                                                  /*comm_rank*/ worldComm.globalRank(),
                                                                  /*comm_size*/ worldComm.globalSize() );
        }


    void initParticipant( PreciceParticipantSetup const& setup,
                          WorldComm const& worldComm = Environment::worldComm() )
        {
            this->initParticipant( setup.name(), setup.configFileName(), worldComm );
        }

    /**
     * @brief Inspect and log preCICE configuration information.
     */
    void inspectPreciceConfiguration()
        {
            // TODO : inspect the precice configuration and print some information about the coupling
        }

    /**
     * @brief Initialize coupling mesh on an element range.
     * @tparam Location Discretization location (vertices or barycenter).
     * @tparam RangeElts Element range type.
     * @param preciceMeshName Name of preCICE mesh.
     * @param rangeElts Range of elements.
     * @return Raw pointer to created PreciceCouplingMesh instance.
     */
    template <PreciceCouplingMeshLocation Location = PreciceCouplingMeshLocation::vertices, concepts::RangeOnElements RangeElts>
    auto initCouplingMesh( std::string const& preciceMeshName, RangeElts const& rangeElts )
        {
            using precice_coupling_mesh_type = PreciceCouplingMesh<Range<mesh_type, MESH_ELEMENTS>, Location>;
            M_couplingMesh[preciceMeshName] = std::make_unique<precice_coupling_mesh_type>( *M_precice, rangeElts, preciceMeshName );
            return std::get<std::unique_ptr<precice_coupling_mesh_type>>( M_couplingMesh.at( preciceMeshName ) ).get();
        }

    /**
     * @brief Initialize coupling mesh on a face range (trace mesh).
     * @tparam Location Discretization location (vertices or barycenter).
     * @tparam RangeFaces Face range type.
     * @param preciceMeshName Name of preCICE mesh.
     * @param rangeFaces Range of faces.
     * @return Raw pointer to created PreciceCouplingMesh instance.
     */
    template <PreciceCouplingMeshLocation Location = PreciceCouplingMeshLocation::vertices, concepts::RangeOnFaces RangeFaces>
    auto initCouplingMesh( std::string const& preciceMeshName, RangeFaces const& rangeFaces )
        {
            using precice_coupling_mesh_type = PreciceCouplingMesh<Range<mesh_type, MESH_FACES>, Location>;
            M_couplingMesh[preciceMeshName] = std::make_unique<precice_coupling_mesh_type>( *M_precice, rangeFaces, preciceMeshName );
            return std::get<std::unique_ptr<precice_coupling_mesh_type>>( M_couplingMesh.at( preciceMeshName ) ).get();
        }

    /**
     * @brief Initialize all coupling meshes based on provided setup and callback function.
     * @tparam F Callback type.
     * @param setup PreciceParticipantSetup instance containing coupling mesh configurations.
     * @param func Callback function invoked for each initialized coupling mesh and should return a range (element or face) for the mesh.
     */
    template <typename F>
    void initAllCouplingMeshes( PreciceParticipantSetup const& setup, F && func )
        {
            for ( auto const& [name, meshSetup] : setup.couplingMeshes() )
            {
                auto range = std::invoke( func, meshSetup );
                variant_mesh_range_type rangeVariant = std::invoke( func, meshSetup );
                std::visit( [this,&name,&meshSetup]( auto const& range ){
                                if ( meshSetup.isLocatedOnVertices() )
                                    this->initCouplingMesh<PreciceCouplingMeshLocation::vertices>( name, range );
                                else if ( meshSetup.isLocatedOnBarycenter() )
                                    this->initCouplingMesh<PreciceCouplingMeshLocation::barycenter>( name, range );
                            }, rangeVariant );
            }
        }

    /* @brief Initialize all write data sets across registered coupling meshes.
     * @tparam Shape Data shape (scalar or vector).
     * @tparam F Callback type.
     * @param setup PreciceParticipantSetup instance containing coupling mesh and data configurations.
     * @param func Callback function invoked for each initialized write data set.
     */
    template <PreciceData::ShapeType Shape, typename F>
    void initAllWriteData( PreciceParticipantSetup const& setup, F && func )
        {
            for ( auto const& [name, meshSetup] : setup.couplingMeshes() )
            {
                auto itFindMesh = M_couplingMesh.find( name );
                if ( itFindMesh == M_couplingMesh.end() )
                    throw std::runtime_error( fmt::format( "PreciceAdapter::initAllWriteData : no coupling mesh with name {} has been initialized", name ) );
                std::visit( [this,&setup,&meshSetup,&func]( auto & couplingMeshPtr ){
                                for ( std::string const& dataName : meshSetup.writeDataNames() )
                                {
                                    auto dataSetup = setup.couplingData( dataName );
                                    if ( dataSetup.is( Shape ) )
                                    {
                                        couplingMeshPtr->template initWriteData<Shape>
                                            ( dataName,
                                              dataSetup.properties(),
                                              [func,couplingMeshPtr=couplingMeshPtr.get()]( auto & writeData ){
                                                  std::invoke( func, writeData, *couplingMeshPtr );
                                              } );
                                    }
                                }
                            }, itFindMesh->second );
            }
        }

    /* @brief Initialize all read data sets across registered coupling meshes.
     * @tparam Shape Data shape (scalar or vector).
     * @tparam F Callback type.
     * @param setup PreciceParticipantSetup instance containing coupling mesh and data configurations.
     * @param func Callback function invoked for each initialized read data set.
     */
    template <PreciceData::ShapeType Shape,typename F>
    void initAllReadData( PreciceParticipantSetup const& setup, F && func )
        {
            for ( auto const& [name, meshSetup] : setup.couplingMeshes() )
            {
                auto itFindMesh = M_couplingMesh.find( name );
                if ( itFindMesh == M_couplingMesh.end() )
                    throw std::runtime_error( fmt::format( "PreciceAdapter::initAllReadData : no coupling mesh with name {} has been initialized", name ) );
                std::visit( [this,&setup,&meshSetup,&func]( auto & couplingMeshPtr ){
                                for ( std::string const& dataName : meshSetup.readDataNames() )
                                {
                                    auto dataSetup = setup.couplingData( dataName );
                                    if ( dataSetup.is( Shape ) )
                                    {
                                        couplingMeshPtr->template initReadData<Shape>
                                            ( dataName,
                                              dataSetup.properties(),
                                              [func,couplingMeshPtr=couplingMeshPtr.get()]( auto & readData ){
                                                  std::invoke( func, readData, *couplingMeshPtr );
                                              } );
                                    }
                                }
                            }, itFindMesh->second );
            }
        }

    /**
     * @brief Trigger preCICE initialization and configuration checks.
     */
    void updateForUse()
        {
            M_precice->initialize();
            this->inspectPreciceConfiguration();
        }

    /**
     * @brief Read all data across registered coupling meshes from preCICE.
     */
    void readData()
        {
            for ( auto & [meshName, couplingMesh] : M_couplingMesh )
            {
                std::visit( [this]( auto & couplingMeshPtr ){
                                couplingMeshPtr->readAllData();
                            }, couplingMesh );
            }
        }

    /**
     * @brief Write all data across registered coupling meshes to preCICE.
     */
    void writeData()
        {
            for ( auto & [meshName, couplingMesh] : M_couplingMesh )
            {
                std::visit( [this]( auto & couplingMeshPtr ){
                                couplingMeshPtr->writeAllData();
                            }, couplingMesh );
            }
        }

    /**
     * @brief Advance coupling time step in preCICE.
     * @param dt Time step size.
     */
    void advance( double dt )
        {
            this->M_precice->advance( dt );
        }

protected:
    std::unique_ptr<precice::Participant> M_precice;
    std::map<std::string, std::variant<std::unique_ptr<PreciceCouplingMesh<Range<mesh_type, MESH_ELEMENTS>, PreciceCouplingMeshLocation::vertices>>,
                                       std::unique_ptr<PreciceCouplingMesh<Range<mesh_type, MESH_ELEMENTS>, PreciceCouplingMeshLocation::barycenter>>,
                                       std::unique_ptr<PreciceCouplingMesh<Range<mesh_type, MESH_FACES>, PreciceCouplingMeshLocation::vertices>>,
                                       std::unique_ptr<PreciceCouplingMesh<Range<mesh_type, MESH_FACES>, PreciceCouplingMeshLocation::barycenter>>
                                       >> M_couplingMesh;
};

} // namespace Feel

#endif // FEELPP_HAS_PRECICE
#endif // FEELPP_FEELDISCR_PRECICE_HPP
