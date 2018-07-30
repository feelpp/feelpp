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
#include <pybind11/pybind11.h>

#include <feel/feeldiscr/pch.hpp>
#include <mpi4py/mpi4py.h>

namespace py = pybind11;
using namespace Feel;

template<typename MeshT, int Order = 1>
class MyElement: public Pch_type<MeshT,Order>::element_type
{
public:
    using base = typename Pch_type<MeshT,Order>::element_type;
    MyElement() : base() {}
};
template<typename MeshT, int Order = 1>
void defDiscr(py::module &m)
{
    using namespace Feel;
    
    using pch_t = Pch_type<MeshT,Order>;
    using pch_ptr_t = Pch_ptrtype<MeshT,Order>;
    using mesh_support_vector_t = typename pch_t::mesh_support_vector_type;
    using periodicity_t = typename pch_t::periodicity_type;
    using mesh_ptr_t = std::shared_ptr<MeshT>;
    using element_t = typename pch_t::element_type;
    std::string pyclass_name;

    std::string suffix = std::to_string(MeshT::nDim)+std::string("D_P") + std::to_string(Order);
    pyclass_name = std::string("mesh_support_vector_")+suffix;
    py::class_<mesh_support_vector_t>(m,pyclass_name.c_str()).def(py::init<>());
    pyclass_name = std::string("periodicity_")+suffix;
    py::class_<periodicity_t>(m,pyclass_name.c_str()).def(py::init<>());
    
    
    pyclass_name = std::string("Pch_") + suffix;
    py::class_<pch_t,std::shared_ptr<pch_t>>(m,pyclass_name.c_str())
        .def(py::init<mesh_ptr_t const&,mesh_support_vector_t const&, size_type, periodicity_t, std::vector<WorldComm> const&, std::vector<bool>>(),
             py::arg("mesh"),
             py::arg("support")=mesh_support_vector_t(),
             py::arg("components")=MESH_RENUMBER | MESH_CHECK,
             py::arg("periodicity")=periodicity_t(),
             py::arg("worldsComm") = Environment::worldsComm(1),
             py::arg("extendedDofTable") = std::vector<bool>(1,false) )
        .def("nDof",static_cast<size_type(pch_t::*)() const>(&pch_t::nDof), "get the number of degrees of freedom over the whole domain")
        .def("nLocalDof",static_cast<size_type(pch_t::*)() const>(&pch_t::nLocalDof), "get the number of degrees of freedom over the current subdomain")
        .def("nLocalDofWithGhost",static_cast<size_type(pch_t::*)() const>(&pch_t::nLocalDofWithGhost), "get the number of degrees of freedom over the current subdomain withthe ghost")
        .def("nLocalDofWithoutGhost",static_cast<size_type(pch_t::*)() const>(&pch_t::nLocalDofWithoutGhost), "get the number of degrees of freedom over the current subdomain without the ghost")
        .def("basisName",static_cast<std::string (pch_t::*)() const>(&pch_t::basisName), "get the basis function name")

        .def("mesh",static_cast<mesh_ptr_t const&(pch_t::*)() const>(&pch_t::mesh), "get the mesh of the function space")
        .def("element",static_cast<element_t (pch_t::*)(std::string const&, std::string const&)>(&pch_t::element), "get an element of the function space", py::arg("name")="u", py::arg("desc")="u")
        ;

    // Element
    pyclass_name = std::string("Element_Pch_") + suffix;
    py::class_<element_t,std::shared_ptr<element_t>>(m,pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<std::shared_ptr<pch_t> const&, std::string const&, std::string const&, size_type, ComponentType>(),py::arg("space"), py::arg("name"), py::arg("desc"), py::arg("start")=0, py::arg("ct")= ComponentType::NO_COMPONENT)
        .def("functionSpace",static_cast<pch_ptr_t const&(element_t::*)() const>(&element_t::functionSpace), "Get funtion space from element")
        .def("size", static_cast< size_type (element_t::*)() const>(&element_t::size), "Get size of element")
        ;
}
    
PYBIND11_MODULE(discr, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;

    std::string pyclass_name = std::string("ComponentType");
    py::class_<ComponentType>(m,pyclass_name.c_str()).def(py::init<>());
    
    //defDiscr<Mesh<Simplex<1>>,1>( m );
    defDiscr<Mesh<Simplex<2>>,1>( m );
    //defDiscr<Mesh<Simplex<2>>,2>( m );
}

