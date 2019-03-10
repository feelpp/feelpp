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
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelvf/ginac.hpp>
#include <mpi4py/mpi4py.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace Feel;

//PYBIND11_MAKE_OPAQUE(Feel::worldscomm_ptr_t);

template<typename MeshT, int Order = 1>
class MyElement: public Pch_type<MeshT,Order>::element_type
{
public:
    using base = typename Pch_type<MeshT,Order>::element_type;
    MyElement() : base() {}
};
template<typename SpaceT>
void defDiscr(py::module &m)
{
    using namespace Feel;
    
    using space_t = SpaceT;
    using space_ptr_t = std::shared_ptr<space_t>;
    using mesh_support_vector_t = typename space_t::mesh_support_vector_type;
    using periodicity_t = typename space_t::periodicity_type;
    using mesh_t = typename space_t::mesh_type;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    using element_t = typename space_t::element_type;
    std::string pyclass_name;
    int Order = space_t::basis_0_type::nOrder;
    std::string suffix = std::to_string(mesh_t::nDim)+std::string("D_P") + std::to_string(Order);
    //std::cout << "suffix discr: " << suffix << std::endl;
    //py::bind_vector<worldscomm_ptr_t>(m, "WorldsComm");
    if ( space_t::is_continuous )
        pyclass_name = std::string("Pch_") + suffix;
    else
        pyclass_name = std::string("Pdh_") + suffix;
    py::class_<space_t,std::shared_ptr<space_t>>(m,pyclass_name.c_str())
        .def(py::init<mesh_ptr_t const&,mesh_support_vector_t const&, size_type, periodicity_t, worldscomm_ptr_t const&, std::vector<bool>>(),
             py::arg("mesh"),
             py::arg("support")=mesh_support_vector_t(),
             py::arg("components")=MESH_RENUMBER | MESH_CHECK,
             py::arg("periodicity")=periodicity_t(),
             py::arg("worldsComm") = Environment::worldsComm(1),
             py::arg("extendedDofTable") = std::vector<bool>(1,false) )
        .def("nDof",static_cast<size_type(space_t::*)() const>(&space_t::nDof), "get the number of degrees of freedom over the whole domain")
        .def("nLocalDof",static_cast<size_type(space_t::*)() const>(&space_t::nLocalDof), "get the number of degrees of freedom over the current subdomain")
        .def("nLocalDofWithGhost",static_cast<size_type(space_t::*)() const>(&space_t::nLocalDofWithGhost), "get the number of degrees of freedom over the current subdomain withthe ghost")
        .def("nLocalDofWithoutGhost",static_cast<size_type(space_t::*)() const>(&space_t::nLocalDofWithoutGhost), "get the number of degrees of freedom over the current subdomain without the ghost")
        .def("basisName",static_cast<std::string (space_t::*)() const>(&space_t::basisName), "get the basis function name")

        .def("mesh",static_cast<mesh_ptr_t const&(space_t::*)() const>(&space_t::mesh), "get the mesh of the function space")
        .def("element",static_cast<element_t (space_t::*)(std::string const&, std::string const&)>(&space_t::element), "get an element of the function space", py::arg("name")="u", py::arg("desc")="u")
        .def("elementFromExpr",static_cast<element_t (space_t::*)(std::string const&, std::string const&, std::string const& )>(&space_t::elementFromExpr), "get an element of the function space interpolating the expression", py::arg("expr"),py::arg("name")="u", py::arg("desc")="u")
        ;

    // Element
    if ( space_t::is_continuous )
        pyclass_name = std::string("Element_Pch_") + suffix;
    else
        pyclass_name = std::string("Element_Pdh_") + suffix;
    py::class_<element_t,std::shared_ptr<element_t>>(m,pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<std::shared_ptr<space_t> const&, std::string const&, std::string const&, size_type, ComponentType>(),py::arg("space"), py::arg("name"), py::arg("desc"), py::arg("start")=0, py::arg("ct")= ComponentType::NO_COMPONENT)
        .def("functionSpace",static_cast<space_ptr_t const&(element_t::*)() const>(&element_t::functionSpace), "Get funtion space from element")
        .def("size", static_cast< size_type (element_t::*)() const>(&element_t::size), "Get size of element")

        .def("save", &element_t::saveImpl, py::arg("path"), py::arg("type")="binary", py::arg("suffix")="", py::arg("sep")="", "save functionspace element in file ")
        .def("load", &element_t::loadImpl, py::arg("path"), py::arg("type")="binary", py::arg("suffix")="", py::arg("sep")="", "load functionspace element from file ")

        .def("on", static_cast<void (element_t::*)( elements_reference_wrapper_t<mesh_ptr_t> const&, Expr<GinacEx<2>> const&, std::string const&, GeomapStrategyType, bool, bool)>(&element_t::template onImpl<elements_reference_wrapper_t<mesh_ptr_t>,Expr<GinacEx<2>>>),
             py::arg("range"), py::arg("expr"), py::arg("prefix")="",
             py::arg("geomap")=GeomapStrategyType::GEOMAP_OPT, py::arg("accumulate")=false, py::arg("verbose")=false, "build the interpolant of the expression expr on a range of elements")
        ;

    if ( space_t::is_continuous )
        pyclass_name = std::string("I_Pch_") + suffix;
    else
        pyclass_name = std::string("I_Pdh_") + suffix;
    //using  Feel::detail::opinterprangetype<IteratorRange>::type
    py::class_<I_t<space_t,space_t>,std::shared_ptr<I_t<space_t,space_t>>>(m,pyclass_name.c_str())
        .def(py::init<>())
        //.def(py::init<std::shared_ptr<space_t> const&, std::shared_ptr<space_t> const&>())
        ;
    
}
template<typename space_t>
void defDiscrDiscontinuous(py::module &m )
{
    if ( (space_t::is_continuous == false) && ( space_t::basis_0_type::nOrder == 0 ) )
    {
        m.def( "pid", &regionProcess<space_t>, "get an piecewise constant function storing the process ids", py::arg("space") );
        //m.def( "marker", &regionMarker<space_t>, "get an piecewise constant function storing the marker of the mesh", py::arg("space") );
    }
}
    
PYBIND11_MODULE(_discr, m )
{
    using namespace Feel;
    using namespace Feel::vf;

    if (import_mpi4py()<0) return ;

    std::string pyclass_name = std::string("ComponentType");
    py::enum_<ComponentType>(m,pyclass_name.c_str())
        .value("NO_COMPONENT", ComponentType::NO_COMPONENT )
        .value("X", ComponentType::X )
        .value("Y", ComponentType::Y )
        .value("Z", ComponentType::Z )
        .value("NX", ComponentType::NX )
        .value("NY", ComponentType::NY )
        .value("NZ", ComponentType::NZ )
        .value("TX", ComponentType::TX )
        .value("TY", ComponentType::TY )
        .value("TZ", ComponentType::TZ )
        .export_values();

    pyclass_name = std::string("Periodic");
    py::class_<Periodic<double>>(m,pyclass_name.c_str()).def(py::init<>());
    pyclass_name = std::string("PeriodicityPeriodic");
    py::class_<Periodicity<Periodic<double>>>(m,pyclass_name.c_str()).def(py::init<>());
    pyclass_name = std::string("NoPeriodicity");
    py::class_<NoPeriodicity>(m,pyclass_name.c_str()).def(py::init<>());
    pyclass_name = std::string("PeriodicityNoPeriodicity");
    py::class_<Periodicity<NoPeriodicity>>(m,pyclass_name.c_str()).def(py::init<>());

    
    defDiscr<Pch_type<Mesh<Simplex<1>>,1>>( m );
    defDiscr<Pch_type<Mesh<Simplex<1>>,2>>( m );
    defDiscr<Pch_type<Mesh<Simplex<1>>,3>>( m );
    
    defDiscr<Pch_type<Mesh<Simplex<2>>,1>>( m );
    defDiscr<Pch_type<Mesh<Simplex<2>>,2>>( m );
    defDiscr<Pch_type<Mesh<Simplex<2>>,3>>( m );

    defDiscr<Pch_type<Mesh<Simplex<3>>,1>>( m );
    defDiscr<Pch_type<Mesh<Simplex<3>>,2>>( m );
    defDiscr<Pch_type<Mesh<Simplex<3>>,3>>( m );
    
    defDiscr<Pdh_type<Mesh<Simplex<2>>,0>>( m );
    defDiscrDiscontinuous<Pdh_type<Mesh<Simplex<2>>,0>>( m );
    defDiscr<Pdh_type<Mesh<Simplex<2>>,1>>( m );
    defDiscr<Pdh_type<Mesh<Simplex<2>>,2>>( m );
    defDiscr<Pdh_type<Mesh<Simplex<2>>,3>>( m );

    defDiscr<Pdh_type<Mesh<Simplex<3>>,0>>( m );
    defDiscrDiscontinuous<Pdh_type<Mesh<Simplex<3>>,0>>( m );
    defDiscr<Pdh_type<Mesh<Simplex<3>>,1>>( m );
    defDiscr<Pdh_type<Mesh<Simplex<3>>,2>>( m );
    defDiscr<Pdh_type<Mesh<Simplex<3>>,3>>( m );
    //defDiscr<Mesh<Simplex<2>>,2>( m );
}

