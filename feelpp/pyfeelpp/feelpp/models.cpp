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
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <mpi4py/mpi4py.h>
#include<feel/feelmodels/modelproperties.hpp>

namespace py = pybind11;

using namespace Feel;

PYBIND11_MAKE_OPAQUE(ModelMaterial);
PYBIND11_MAKE_OPAQUE(ModelMaterials);
PYBIND11_MAKE_OPAQUE(ModelParameters);

PYBIND11_MODULE(_models, m )
{
    using namespace Feel;

    auto parameter_to_str = [](const ModelParameter &mat )
    {
        std::ostringstream s;
        s << mat.name() << " : " << mat.value() << std::endl;
        return s.str();
    };
    std::string pyclass_name;
    pyclass_name = "ModelParameter";
    py::class_<ModelParameter>(m,pyclass_name.c_str())
        .def(py::init<>())
        .def("name",&ModelParameter::name, "name of the parameter")
        .def("type",&ModelParameter::type, "type of the parameter: value, expression, fit")
        .def("value",&ModelParameter::value, "value of the parameter")
        .def("setValue",&ModelParameter::setValue, "set value of the parameter")
        //.def("hasExpression",&ModelParameter::hasExpression, "return true if the parameter has an expression, false otherwise")
        .def("hasFitInterpolator",&ModelParameter::hasFitInterpolator, "return true if the parameter has a fit interpolator, false otherwise")
        .def("setParameterValues",&ModelParameter::setParameterValues, "set parameter values from a map of string/double pairs")
        .def("__str__", parameter_to_str, "");



    pyclass_name = "ModelParameters";
    py::class_<ModelParameters>(m,pyclass_name.c_str())
        .def(py::init<worldcomm_ptr_t const&>(),
             py::arg("worldComm"))
        .def("clear", &ModelParameters::clear)
        .def("__len__", [](const ModelParameters &v) { return v.size(); })
        .def("__iter__", [](ModelParameters &v) {
                return py::make_iterator(v.begin(), v.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def("__str__", [](ModelParameters &mp) {
                std::ostringstream s;
                s << "{";
                for( auto const& [k,v] : mp )
                {
                    s << "{" << v.name() << "," << v.value() << "},";
                }
                s << "}";
                return s.str();
            },"")
        .def("at",static_cast<ModelParameter& (ModelParameters::*)(std::string const&)>(&ModelParameters::at),"",
             py::return_value_policy::reference, py::keep_alive<0, 1>())
        .def("setParameterValues",&ModelParameters::setParameterValues, "set parameter values from a map of string/double pairs")
        .def("toParameterValues",&ModelParameters::toParameterValues, "get a dictionary from the map of parameter values");

    pyclass_name = "ModelMaterial";
    py::class_<ModelMaterial>(m,pyclass_name.c_str())
        .def(py::init<>())
        .def("__getitem__", [](const ModelMaterial &map, std::string key) {
                try { return map.property(key).exprScalar().expression().exprDesc(); }
                catch (const std::out_of_range&) {
                    throw py::key_error("key '" + key + "' does not exist");
                }
            })

        .def("setProperty",[](ModelMaterial& m, const std::string& prop, const std::string& e )
             {
                 m.setProperty( prop, e );
             }, "returns true of the property exists, false otherwise")
        .def("hasProperty",&ModelMaterial::hasProperty, "returns true of the property exists, false otherwise")
        .def("hasPropertyConstant",&ModelMaterial::hasPropertyConstant, "returns true of the property exists and is constant, false otherwise")
        .def("hasPropertyScalar",&ModelMaterial::hasPropertyExprScalar, "returns true of the property exists and is a scalar expression, false otherwise")
        .def("propertyConstant",&ModelMaterial::propertyConstant, "return the value of the constant property")
        .def("setParameterValues",&ModelMaterial::setParameterValues, "set parameter values from a map of string/double pairs")
        .def("__str__", [](const ModelMaterial &mat )
             {
                 std::ostringstream s;
                 s << mat.name() << std::endl;
                 s << " . markers: ";
                 for( auto const& p : mat.meshMarkers() )
                     s << p << " ";
                 s << std::endl;
                 s << " . physics: ";
                 for( auto const& p : mat.physics() )
                     s << p << " ";
                 s << std::endl;
                 s << " . properties: " << std::endl;
                 for( auto const& [p,value] : mat.properties() )
                 {
                     if ( value.isConstant() )
                         s << "   {" << p << ", " << value.value() << "}" << std::endl;
                     if ( value.isScalar() )
                         s << "   {" << p << ", " << value.exprScalar().expression().exprDesc() << "}" << std::endl;
                     if ( value.hasExprVectorial2() )
                         s << "   {" << p << ", " << value.exprVectorial2().expression().exprDesc() << "}" << std::endl;
                     if ( value.hasExprVectorial3() )
                         s << "   {" << p << ", " << value.exprVectorial3().expression().exprDesc() << "}" << std::endl;
                 }
                 s << std::endl;
                 return s.str();
             }, "")

        ;

    pyclass_name = "ModelMaterials";
    py::class_<ModelMaterials>(m,pyclass_name.c_str())
        .def(py::init<worldcomm_ptr_t const&>(),
             py::arg("worldComm"))
        .def("clear", &ModelMaterials::clear)
        .def("__len__", [](const ModelMaterials &v) { return v.size(); })
        .def("__str__", [](ModelMaterials &v) {
                std::ostringstream s;
                for (auto const& [key,mat] : v )
                    s << mat.name() << " " << std::endl;
                return s.str();
            })
#if 0
        .def("__getitem__", [](const ModelMaterials &map, std::string key) {
                try { return map.at(key); }
                catch (const std::out_of_range&) {
                    throw py::key_error("key '" + key + "' does not exist");
                }
            })
#endif
        .def("__getitem__", [](ModelMaterials &map, std::string key) {
                try { return map.at(key); }
                catch (const std::out_of_range&) {
                    throw py::key_error("key '" + key + "' does not exist");
                }
            })
        .def("at", static_cast<ModelMaterial& (ModelMaterials::*)(std::string const&)>(&ModelMaterials::at),"")
        .def("__iter__", [](ModelMaterials &v) {
                return py::make_iterator(v.begin(), v.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def("setParameterValues",&ModelMaterials::setParameterValues, "set parameter values from a map of string/double pairs")
        .def("items", [](ModelMaterials &map) { return py::make_iterator(map.begin(), map.end()); },
             py::keep_alive<0, 1>());

    pyclass_name = "ModelProperties";
    py::class_<ModelProperties,std::shared_ptr<ModelProperties>>(m,pyclass_name.c_str())
        .def(py::init<std::string const&, std::string const&, worldcomm_ptr_t const&, std::string const&>(),"initialize ModelProperties",py::arg("filename")="",py::arg("directoryLibExpr")="",py::arg("worldComm"),py::arg("prefix")="")
        .def("parameters",static_cast<ModelParameters& (ModelProperties::*)()>(&ModelProperties::parameters), "get parameters of the model",py::return_value_policy::reference)
        .def("materials",static_cast<ModelMaterials& (ModelProperties::*)()>(&ModelProperties::materials), "get the materials of the model",py::return_value_policy::reference);

}
