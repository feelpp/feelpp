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
//! @date 25 Jul 2018
//! @copyright 2018 Feel++ Consortium
//!
#include <feel/feelpython/pybind11/pybind11.h>
#include <feel/feelpython/pybind11/stl.h>
#include <feel/feelpython/pybind11/json.h>

#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feelmodels/modelcore/modelalgebraic.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/options.hpp>


namespace py = pybind11;
using namespace Feel;


template <typename IndexType>
void bind_ModelMeshes( py::module& m )
{
    namespace py = pybind11;
    using MeshBase = Feel::MeshBase<IndexType>;
    using ModelMeshCommon = Feel::FeelModels::ModelMeshCommon<IndexType>;
    using ModelMesh = Feel::FeelModels::ModelMesh<IndexType>;
    using ModelMeshes = Feel::FeelModels::ModelMeshes<IndexType>;

    py::class_<ModelMeshCommon>( m, "ModelMeshCommon",
                                            R"pbdoc(
    A class representing a common model mesh.

    This class provides functionality for importing, manipulating, and accessing model meshes.

    )pbdoc" )
        .def( py::init<>(), R"pbdoc(
    Default constructor.

    Creates an empty ModelMeshCommon instance.

    )pbdoc" )
        .def( py::init<ModelMeshes const&>(), R"pbdoc(
    Constructor with ModelMeshes parameter.

    Creates a ModelMeshCommon instance with the specified ModelMeshes.

    Parameters:
        mMeshes (ModelMeshes): The ModelMeshes instance.

    )pbdoc" )
        .def_property( "importConfig",
                        py::overload_cast<>( &ModelMeshCommon::importConfig, py::const_ ),
                        py::overload_cast<>( &ModelMeshCommon::importConfig ),
                        R"pbdoc(
                The import configuration.

                This property gets or sets the import configuration for the model mesh.

                Returns:
                    ImportConfig: The import configuration.

                )pbdoc" )
        .def( "hasMesh", &ModelMeshCommon::hasMesh, R"pbdoc(
    Check if the model has a mesh.

    Returns:
        bool: True if the model has a mesh, False otherwise.

    )pbdoc" )
        .def(
            "mesh", []( ModelMeshCommon& self )
                    {
                        return self.mesh();
                    },
                    R"pbdoc(
    Get the mesh.

    Returns:
        MeshType: The mesh object.

    )pbdoc" )
        .def( "meshFilename", &ModelMeshCommon::meshFilename, py::return_value_policy::reference, R"pbdoc(
    Get the mesh filename.

    Returns:
        str: The mesh filename.

    )pbdoc" )
        .def( "setMesh", &ModelMeshCommon::setMesh, R"pbdoc(
    Set the mesh.

    Parameters:
        m (MeshBase): The mesh object.
        filename (str, optional): The filename associated with the mesh. Defaults to "".

    )pbdoc" )
        .def( "clearFunctionSpaces", &ModelMeshCommon::clearFunctionSpaces, R"pbdoc(
    Clear all function spaces.

    Removes all function spaces associated with the model.

    Parameters:
        None

    Returns:
        None

    )pbdoc" );

    // ModelMeshCommon::ImportConfig
    py::class_<typename ModelMeshCommon::ImportConfig>( m, "ImportConfig",
                                                        R"pbdoc(
    Configuration settings for importing a model mesh.

    This class provides various options for importing a model mesh.

    )pbdoc" )
        .def( py::init<>(), R"pbdoc(
    Default constructor.

    Creates an empty ImportConfig instance.

    )pbdoc" )
#if 0    
        .def_readonly( "inputFilename", &ModelMeshCommon::ImportConfig::inputFilename, R"pbdoc(
    The input filename.

    )pbdoc" )
        .def_readonly( "meshFilename", &ModelMeshCommon::ImportConfig::meshFilename, R"pbdoc(
    The mesh filename.

    )pbdoc" )
        .def_readonly( "geoFilename", &ModelMeshCommon::ImportConfig::geoFilename, R"pbdoc(
    The geometry filename.

    )pbdoc" )
        .def_readonly( "generatePartitioning", &ModelMeshCommon::ImportConfig::generatePartitioning, R"pbdoc(
    Flag indicating whether to generate partitioning.

    )pbdoc" )
        .def_readonly( "numberOfPartition", &ModelMeshCommon::ImportConfig::numberOfPartition, R"pbdoc(
    The number of partitions.

    )pbdoc" )
        .def_readonly( "meshSize", &ModelMeshCommon::ImportConfig::meshSize, R"pbdoc(
    The mesh size.

    )pbdoc" )
        .def_readonly( "straightenMesh", &ModelMeshCommon::ImportConfig::straightenMesh, R"pbdoc(
    Flag indicating whether to straighten the mesh.

    )pbdoc" )
        .def_readonly( "meshComponents", &ModelMeshCommon::ImportConfig::meshComponents, R"pbdoc(
    The mesh components.

    )pbdoc" )
  
        .def_readonly( "loadByMasterRankOnly", &ModelMeshCommon::ImportConfig::loadByMasterRankOnly, R"pbdoc(
    Flag indicating whether to load by master rank only.

    )pbdoc" )
#endif
        .def( "setStraightenMesh", &ModelMeshCommon::ImportConfig::setStraightenMesh, R"pbdoc(
    Set the straighten mesh flag.

    Parameters:
        b (bool): The value of the flag.

    )pbdoc" )
        .def( "setMeshComponents", &ModelMeshCommon::ImportConfig::setMeshComponents, R"pbdoc(
    Set the mesh components.

    Parameters:
        c (int): The mesh components.

    )pbdoc" )
        .def( "hasMeshFilename", &ModelMeshCommon::ImportConfig::hasMeshFilename, R"pbdoc(
    Check if the import configuration has a mesh filename.

    Returns:
        bool: True if the mesh filename is set, False otherwise.

    )pbdoc" )
        .def( "hasGeoFilename", &ModelMeshCommon::ImportConfig::hasGeoFilename, R"pbdoc(
    Check if the import configuration has a geometry filename.

    Returns:
        bool: True if the geometry filename is set, False otherwise.

    )pbdoc" )
        .def( "setupInputMeshFilenameWithoutApplyPartitioning", &ModelMeshCommon::ImportConfig::setupInputMeshFilenameWithoutApplyPartitioning, R"pbdoc(
    Setup the input mesh filename without applying partitioning.

    Parameters:
        filename (str): The input mesh filename.

    )pbdoc" )
        .def( "setupSequentialAndLoadByMasterRankOnly", &ModelMeshCommon::ImportConfig::setupSequentialAndLoadByMasterRankOnly, R"pbdoc(
    Setup for sequential execution and load by master rank only.

    Parameters:
        None

    )pbdoc" )
        .def( "updateInformationObject", &ModelMeshCommon::ImportConfig::updateInformationObject, R"pbdoc(
    Update the information object.

    Parameters:
        p (dict): The information object.

    )pbdoc" )
        .def_static( "tabulateInformations", &ModelMeshCommon::ImportConfig::tabulateInformations, R"pbdoc(
    Tabulate informations.

    Parameters:
        jsonInfo (str): The JSON information.
        tabInfoProp (TabulateInformationProperties): The tabulate information properties.

    Returns:
        TabulateInformations: The tabulated informations.

    )pbdoc" );

    py::class_<ModelMesh>( m, "ModelMesh", py::module_local() )
        .def( py::init<std::string const&>(), py::arg( "name" ) )
        .def( py::init<std::string const&, ModelMeshes const&>(), py::arg( "name" ), py::arg( "mMeshes" ) )
        .def( "metadata", &ModelMesh::metadata, py::return_value_policy::reference, R"pbdoc(
        Return the JSON metadata.

        Returns:
            nl::json: The JSON metadata.

        )pbdoc" )
        .def( "setup", &ModelMesh::setup, py::arg( "jarg" ), py::arg( "mMeshes" ), R"pbdoc(
        Setup the model mesh.

        Parameters:
            jarg (nl::json): The JSON configuration.
            mMeshes (ModelMeshes): The model meshes.

        Returns:
            None

        )pbdoc" )
        .def( "setupRestart", &ModelMesh::setupRestart, py::arg( "mMeshes" ), R"pbdoc(
        Setup the model mesh for restart.

        Parameters:
            mMeshes (ModelMeshes): The model meshes.

        Returns:
            None

        )pbdoc" )
        .def( "setMesh", &ModelMesh::setMesh, py::arg( "m" ), py::arg( "meshFilename" ) = "", R"pbdoc(
        Set the mesh.

        Parameters:
            m (MeshBase): The mesh.
            meshFilename (str): The mesh filename.

        Returns:
            None

        )pbdoc" )
        .def( "setAsShared", &ModelMesh::setAsShared, py::arg( "m" ), R"pbdoc(
        Set the model mesh as shared.

        Parameters:
            m (ModelMesh): The model mesh.

        Returns:
            None

        )pbdoc" )
        .def( "importConfig", py::overload_cast<>( &ModelMesh::importConfig ), py::return_value_policy::reference_internal, R"pbdoc(
        Get the import configuration.

        Returns:
            ImportConfig: The import configuration.

        )pbdoc" )
        .def( "importConfig", py::overload_cast<>( &ModelMesh::importConfig, py::const_ ), py::return_value_policy::reference_internal, R"pbdoc(
        Get the import configuration (const version).

        Returns:
            ImportConfig: The import configuration.

        )pbdoc" )
#if 0        
        .def( "updateForUse", &ModelMesh::template updateForUse<MeshBase>, py::arg( "mMeshes" ), R"pbdoc(
        Update the model mesh for use.

        Parameters:
            mMeshes (ModelMeshes): The model meshes.

        Returns:
            None

        )pbdoc" )
        .def( "initMeasurePointsEvaluationTool", &ModelMesh::template initMeasurePointsEvaluationTool<MeshBase>, R"pbdoc(
        Initialize the measure points evaluation tool.

        Parameters:
            None

        Returns:
            None

        )pbdoc" )
        .def( "measurePointsEvaluationTool", &ModelMesh::template measurePointsEvaluationTool<MeshBase>, R"pbdoc(
        Get the measure points evaluation tool.

        Parameters:
            None

        Returns:
            MeasurePointsEvaluationTool: The measure points evaluation tool.

        )pbdoc" )

        .def( "mesh", py::overload_cast<>( &ModelMesh::template mesh<MeshBase>, py::const_ ), py::return_value_policy::reference_internal, R"pbdoc(
        Get the mesh.

        Returns:
            MeshBase: The mesh.

        )pbdoc" )
#endif        
        .def( "hasMeshMotion", &ModelMesh::hasMeshMotion, R"pbdoc(
        Check if the model mesh has mesh motion.

        Returns:
            bool: True if the model mesh has mesh motion, False otherwise.

        )pbdoc" )
#if 0        
        .def( "meshMotionTool", &ModelMesh::template meshMotionTool<MeshBase>, py::return_value_policy::reference_internal, R"pbdoc(
        Get the mesh motion tool.

        Returns:
            MeshMotionTool: The mesh motion tool.

        )pbdoc" )

        .def( "collectionOfDataByMeshEntity", &ModelMesh::collectionOfDataByMeshEntity, py::return_value_policy::reference_internal, R"pbdoc(
        Get the collection of data by mesh entity.

        Returns:
            dict: The collection of data by mesh entity.

        )pbdoc" )

        .def( "basisFieldTypeSupported", &ModelMesh::template basisFieldTypeSupported<MeshBase>, py::return_value_policy::copy, R"pbdoc(
        Return the supported basis field types.

        Returns:
            tuple: A tuple of tuples containing the names and types of the basis fields supported by the mesh.

        )pbdoc" )
       
        .def( "modelFields", &ModelMesh::template modelFields<MeshBase>, py::arg( "prefix_field" ) = "", py::arg( "prefix_symbol" ) = "", R"pbdoc(
        Return the model fields for a given mesh.

        Parameters:
            prefix_field (str): A prefix to use for the field names.
            prefix_symbol (str): A prefix to use for the symbol names.

        Returns:
            tuple: A tuple of model fields for the given mesh.

        )pbdoc" )
        .def( "symbolsExpr", &ModelMesh::template symbolsExpr<MeshBase>, py::arg( "prefix_symbol" ) = "", py::return_value_policy::copy, R"pbdoc(
        Return the symbols expressions for the model mesh.

        Parameters:
            prefix_symbol (str): A prefix to use for the symbol names.

        Returns:
            dict: The symbols expressions for the model mesh.

        )pbdoc" )
#endif         
        .def( "updateTime", &ModelMesh::updateTime, py::arg( "time" ), R"pbdoc(
        Update the model mesh time.

        Parameters:
            time (float): The time value.

        Returns:
            None

        )pbdoc" )
        .def( "updateInformationObject", &ModelMesh::updateInformationObject, py::arg( "p" ), py::arg( "prefix_symbol" ), R"pbdoc(
        Update the information object.

        Parameters:
            p (nl::json): The information object.
            prefix_symbol (str): A prefix to use for the symbol names.

        Returns:
            None

        )pbdoc" )
        .def_static( "tabulateInformations", &ModelMesh::tabulateInformations, py::arg( "p" ), py::arg( "tabInfoProp" ), R"pbdoc(
        Tabulate the informations.

        Parameters:
            p (nl::json): The JSON configuration.
            tabInfoProp (TabulateInformationProperties): The tabulate information properties.

        Returns:
            TabulateInformationsPtr: The tabulate informations.

        )pbdoc" )
        .def( "setParameterValues", &ModelMesh::setParameterValues, py::arg( "paramValues" ), R"pbdoc(
        Set the parameter values.

        Parameters:
            paramValues (dict): The parameter values.

        Returns:
            None

        )pbdoc" )
#if 0        
        .def( "updateDistanceToRange", &ModelMesh::template updateDistanceToRange<MeshBase>, py::arg( "MeshType" ), R"pbdoc(
        Update the distance to range.

        Parameters:
            MeshType: The mesh type.

        Returns:
            None

        )pbdoc" )
        .def( "updateMeshMotion", &ModelMesh::template updateMeshMotion<MeshBase>, py::arg( "MeshType" ), py::arg( "se" ), R"pbdoc(
        Update the mesh motion.

        Parameters:
            MeshType: The mesh type.
            se: The symbols expression.

        Returns:
            None

        )pbdoc" )
        .def( "updateMeshAdaptation", &ModelMesh::template updateMeshAdaptation<MeshBase>, py::arg( "event" ), py::arg( "se" ), R"pbdoc(
        Update the mesh adaptation.

        Parameters:
            event: The mesh adaptation event.
            se: The symbols expression.

        Returns:
            None

        )pbdoc" )
#endif        
        .def( "setFunctionApplyRemesh", &ModelMesh::setFunctionApplyRemesh, py::arg( "f" ), R"pbdoc(
        Set the function for applying remesh.

        Parameters:
            f (function): The function.

        Returns:
            None

        )pbdoc" )
#if 0        
        .def( "applyRemesh", &ModelMesh::template applyRemesh<MeshBase>, py::arg( "newMesh" ), R"pbdoc(
        Apply remesh.

        Parameters:
            newMesh (MeshType): The new mesh.

        Returns:
            None

        )pbdoc" )
#endif        
#if 0        
        .def( "updateField",
              [](, ModelMesh & self, std::string const& name, std::string const& e, std::string const& basis, bool initZero )
                {
                    self.updateField( name, expr(e), basis, initZero );
                },
              py::arg( "name" ), py::arg( "expr" ), py::arg( "range" ), py::arg( "basis" ) = "", py::arg( "initZero" ) = false, R"pbdoc(
        Update a field based on an expression and a range.

        Parameters:
            name (str): The name of the field.
            expr (Expr): The expression to update the field.
            range (Range): The range of elements.
            basis (str): The basis function space to be used for the field.
            initZero (bool): Flag indicating whether to initialize the field to zero before updating.

        Returns:
            None

        )pbdoc" )
        .def( "updateField", 
              [](, ModelMesh & self, std::string const& name, std::string const& e, std::string const& basis )
                {
                    self.updateField( name, expr(e), basis );
                },
               py::arg( "name" ), py::arg( "expr" ), py::arg( "basis" ) = "", R"pbdoc(
        Update a field based on an expression.

        Parameters:
            name (str): The name of the field.
            expr (Expr): The expression to update the field.
            basis (str): The basis function space to be used for the field.

        Returns:
            None

        )pbdoc" )
        .def( "updateField", py::overload_cast<std::string const&, FieldType const&, std::string const&>( &ModelMesh::template updateField<MeshBase, FieldType> ), py::arg( "name" ), py::arg( "u" ), py::arg( "basis" ) = "", R"pbdoc(
        Update a field with data from another field.

        Parameters:
            name (str): The name of the field to be updated.
            u (FieldType): The field to be used for the update.
            basis (str): The basis function space to be used for the field.

        Returns:
            None

        )pbdoc" )
#endif        
        .def( "removeField", &ModelMesh::removeField, py::arg( "name" ), R"pbdoc(
        Remove a field from the model mesh.

        Parameters:
            name (str): The name of the field to remove.

        Returns:
            bool: True if the field was successfully removed, False otherwise.

        )pbdoc" );

    // ModelMeshes
    py::class_<ModelMeshes, std::shared_ptr<ModelMeshes>>( m, "ModelMeshes" )
        .def( py::init<>() )
        .def( "hasModelMesh", &ModelMeshes::hasModelMesh, py::arg( "meshName" ), R"pbdoc(
        Check if a model mesh exists.

        Parameters:
            meshName (str): The name of the model mesh.

        Returns:
            bool: True if the model mesh exists, False otherwise.

        )pbdoc" )
        .def( "modelMesh", py::overload_cast<const std::string&>( &ModelMeshes::modelMesh ), py::arg( "meshName" ), R"pbdoc(
        Get a reference to a model mesh.

        Parameters:
            meshName (str): The name of the model mesh.

        Returns:
            ModelMesh: A reference to the model mesh.

        )pbdoc" )
        .def( "modelMesh", py::overload_cast<>( &ModelMeshes::modelMesh ), R"pbdoc(
        Get a reference to the default model mesh.

        Returns:
            ModelMesh: A reference to the default model mesh.

        )pbdoc" )
        .def( "modelMesh", py::overload_cast<const std::string&>( &ModelMeshes::modelMesh, py::const_ ), py::arg( "meshName" ), R"pbdoc(
        Get a const reference to a model mesh.

        Parameters:
            meshName (str): The name of the model mesh.

        Returns:
            ModelMesh: A const reference to the model mesh.

        )pbdoc" )
        .def( "modelMesh", py::overload_cast<>( &ModelMeshes::modelMesh, py::const_ ), R"pbdoc(
        Get a const reference to the default model mesh.

        Returns:
            ModelMesh: A const reference to the default model mesh.

        )pbdoc" )
        .def( "setup", py::overload_cast<const nl::json&, const std::set<std::string>&>( &ModelMeshes::setup ), py::arg( "jarg" ), py::arg( "keywordsToSetup" ), R"pbdoc(
        Setup the model meshes.

        Parameters:
            jarg (dict): The JSON configuration.
            keywordsToSetup (set): A set of keywords to setup.

        Returns:
            None

        )pbdoc" )
        .def( "setup", py::overload_cast<const nl::json&>( &ModelMeshes::setup ), py::arg( "jarg" ), R"pbdoc(
        Setup the model meshes.

        Parameters:
            jarg (dict): The JSON configuration.

        Returns:
            None

        )pbdoc" )
        .def( "setupRestart", &ModelMeshes::setupRestart, py::arg( "meshName" ), R"pbdoc(
        Setup the model mesh for restart.

        Parameters:
            meshName (str): The name of the model mesh.

        Returns:
            None

        )pbdoc" )
        .def( "setMesh", &ModelMeshes::setMesh, py::arg( "meshName" ), py::arg( "m" ), R"pbdoc(
        Set the mesh for a model mesh.

        Parameters:
            meshName (str): The name of the model mesh.
            m (MeshBase): The mesh to set.

        Returns:
            None

        )pbdoc" )
        .def( "setModelMeshAsShared", py::overload_cast<const std::string&, const ModelMesh&>( &ModelMeshes::setModelMeshAsShared ), py::arg( "meshName" ), py::arg( "m" ), R"pbdoc(
        Set a model mesh as shared.

        Parameters:
            meshName (str): The name of the model mesh.
            m (ModelMesh): The model mesh to set as shared.

        Returns:
            None

        )pbdoc" )
        .def( "setModelMeshAsShared", py::overload_cast<const ModelMesh&>( &ModelMeshes::setModelMeshAsShared ), py::arg( "m" ), R"pbdoc(
        Set the default model mesh as shared.

        Parameters:
            m (ModelMesh): The model mesh to set as shared.

        Returns:
            None

        )pbdoc" )
#if 0        
        .def( "updateForUse", &ModelMeshes::template updateForUse<typename ModelMesh::mesh_base_type>, py::arg( "meshName" ), R"pbdoc(
        Update the model mesh for use.

        Parameters:
            meshName (str): The name of the model mesh.

        Returns:
            None

        )pbdoc" )
        .def( "initMeasurePointsEvaluationTool", &ModelMeshes::template initMeasurePointsEvaluationTool<typename ModelMesh::mesh_base_type>, py::arg( "meshName" ), R"pbdoc(
        Initialize the measure points evaluation tool for a model mesh.

        Parameters:
            meshName (str): The name of the model mesh.

        Returns:
            None

        )pbdoc" )
        .def( "measurePointsEvaluationTool", &ModelMeshes::template measurePointsEvaluationTool<typename ModelMesh::mesh_base_type>, py::arg( "meshName" ), R"pbdoc(
        Get the measure points evaluation tool for a model mesh.

        Parameters:
            meshName (str): The name of the model mesh.

        Returns:
            MeasurePointsEvaluationTool: The measure points evaluation tool.

        )pbdoc" )
        .def( "mesh", py::overload_cast<const std::string&>( &ModelMeshes::template mesh<typename ModelMesh::mesh_base_type> ), py::arg( "meshName" ), R"pbdoc(
        Get the mesh of a model mesh.

        Parameters:
            meshName (str): The name of the model mesh.

        Returns:
            MeshBase: The mesh.

        )pbdoc" )
        .def( "modelFields", &ModelMeshes::template modelFields<typename ModelMesh::mesh_base_type>, py::arg( "prefix_field" ), py::arg( "prefix_symbol" ) = "", R"pbdoc(
        Get the model fields of the model meshes.

        Parameters:
            prefix_field (str): The prefix for the field names.
            prefix_symbol (str): The prefix for the symbol names.

        Returns:
            list: The model fields.

        )pbdoc" )
        .def( "symbolsExpr", &ModelMeshes::template symbolsExpr<typename ModelMesh::mesh_base_type, true>, py::arg( "prefix_symbol" ) = "meshes", R"pbdoc(
        Get the symbolic expressions of the model meshes.

        Parameters:
            prefix_symbol (str): The prefix for the symbol names.

        Returns:
            list: The symbolic expressions.

        )pbdoc" )
#endif        
        .def( "updateTime", &ModelMeshes::updateTime, py::arg( "time" ), R"pbdoc(
        Update the time for all model meshes.

        Parameters:
            time (float): The time value.

        Returns:
            None

        )pbdoc" )
#if 0        
        .def( "updateInformationObject", &ModelMeshes::updateInformationObject, py::arg( "p" ), py::arg( "prefix_symbol" ) = "meshes", R"pbdoc(
        Update the information object.

        Parameters:
            p (dict): The information object.
            prefix_symbol (str): The prefix for the symbol names.

        Returns:
            None

        )pbdoc" )
#endif        
        .def( "tabulateInformations", &ModelMeshes::tabulateInformations, py::arg( "jsonInfo" ), py::arg( "tabInfoProp" ), R"pbdoc(
        Tabulate the information object.

        Parameters:
            jsonInfo (dict): The JSON configuration for tabulation.
            tabInfoProp (TabulateInformationProperties): The tabulation information properties.

        Returns:
            TabulateInformations: The tabulated information object.

        )pbdoc" )
        .def( "setParameterValues", &ModelMeshes::setParameterValues, py::arg( "paramValues" ), R"pbdoc(
        Set the parameter values for all model meshes.

        Parameters:
            paramValues (dict): The parameter values.

        Returns:
            None

        )pbdoc" )
        .def( "hasMeshMotion", &ModelMeshes::hasMeshMotion, py::arg( "meshName" ), R"pbdoc(
        Check if a model mesh has mesh motion.

        Parameters:
            meshName (str): The name of the model mesh.

        Returns:
            bool: True if the model mesh has mesh motion, False otherwise.

        )pbdoc" )
#if 0        
        .def( "meshMotionTool", &ModelMeshes::template meshMotionTool<typename ModelMesh::mesh_base_type>, py::arg( "meshName" ), R"pbdoc(
        Get the mesh motion tool for a model mesh.

        Parameters:
            meshName (str): The name of the model mesh.

        Returns:
            MeshMotionTool: The mesh motion tool.

        )pbdoc" )
        .def( "updateMeshMotion", &ModelMeshes::template updateMeshMotion<typename ModelMesh::mesh_base_type>, py::arg( "meshName" ), py::arg( "se" ), R"pbdoc(
        Update the mesh motion for a model mesh.

        Parameters:
            meshName (str): The name of the model mesh.
            se (SymbolsExprType): The symbolic expressions.

        Returns:
            None

        )pbdoc" )
        .def( "updateMeshAdaptation", &ModelMeshes::template updateMeshAdaptation<typename ModelMesh::mesh_base_type>, py::arg( "meshName" ), py::arg( "event" ), py::arg( "se" ), R"pbdoc(
        Update the mesh adaptation for a model mesh.

        Parameters:
            meshName (str): The name of the model mesh.
            event (MeshAdaptation.Event): The mesh adaptation event.
            se (SymbolsExprType): The symbolic expressions.

        Returns:
            None

        )pbdoc" )
        .def( "applyRemesh", &ModelMeshes::template applyRemesh<typename ModelMesh::mesh_base_type>, py::arg( "meshName" ), py::arg( "newMesh" ), R"pbdoc(
        Apply remeshing to a model mesh.

        Parameters:
            meshName (str): The name of the model mesh.
            newMesh (MeshType): The new mesh.

        Returns:
            None

        )pbdoc" )
#endif        
        .def( "repository_meshes", &ModelMeshes::repository_meshes, R"pbdoc(
        Get the repository path for the model meshes.

        Returns:
            str: The repository path.

        )pbdoc" )
        .def( "saveMetadata", &ModelMeshes::saveMetadata, R"pbdoc(
        Save the metadata for the model meshes.

        Returns:
            None

        )pbdoc" );

}

template void bind_ModelMeshes<uint32_t>( py::module& m );
//template void bind_ModelMeshes<int64_t>( py::module& m );