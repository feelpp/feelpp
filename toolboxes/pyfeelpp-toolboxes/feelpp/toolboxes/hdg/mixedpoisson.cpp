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

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/hdg/mixedpoisson.hpp>

namespace py = pybind11;
using namespace Feel;

template<int nDim, int OrderPotential>
void defMP(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using toolbox_t = MixedPoisson< Simplex<nDim, 1>, OrderPotential>;
    using space_flux_ptrtype = typename toolbox_t::space_flux_ptrtype;
    using element_flux_type = typename toolbox_t::element_flux_type;
    using element_flux_ptrtype = typename toolbox_t::element_flux_ptrtype;
    using space_potential_ptrtype = typename toolbox_t::space_potential_ptrtype;
    using element_potential_type = typename toolbox_t::element_potential_type;
    using element_potential_ptrtype = typename toolbox_t::element_potential_ptrtype;
    using space_postpotential_ptrtype = typename toolbox_t::space_postpotential_ptrtype;
    using element_postpotential_type = typename toolbox_t::element_postpotential_type;

    std::string pyclass_name = std::string("mixedpoisson_") + std::to_string(nDim) + std::string("DP") + std::to_string(OrderPotential);
    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,MixedPoissonPhysics const&,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("physic")=MixedPoissonPhysics::None,
             py::arg("worldComm")=Environment::worldCommPtr(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the electric mechanics toolbox"
             )
        .def("init",&toolbox_t::init, "initialize the electric mechanics toolbox",py::arg("buildModelAlgebraicFactory")= true)

        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        .def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )

        // elements
        .def( "spaceFlux", &toolbox_t::spaceFlux, "returns the flux function space")
        .def( "fieldFlux", &toolbox_t::fieldFlux, "returns the flux field")
        .def( "fieldFluxPtr", &toolbox_t::fieldFluxPtr, "returns the flux field shared_ptr")
        .def( "spacePotential", &toolbox_t::spacePotential, "returns the potential function space")
        .def( "fieldPotential", &toolbox_t::fieldPotential, "returns the potential field")
        .def( "fieldPotentialPtr", &toolbox_t::fieldPotentialPtr, "returns the potential field shared_tr")
        .def( "spacePostPotential", &toolbox_t::spacePostPotential, "returns the post processed potential function space")
        .def( "fieldPostPotential", &toolbox_t::fieldPostPotential, "returns the post processed potential field")
        .def( "fieldPostPotentialPtr", &toolbox_t::fieldPostPotentialPtr, "returns the post processed potential field shared_tr")

        // solve
        .def("solve",&toolbox_t::solve, "solve the electric mechanics problem, set boolean to true to update velocity and acceleration")
        .def("exportResults",static_cast<void (toolbox_t::*)()>(&toolbox_t::exportResults), "export the results of the electric mechanics problem")
        .def("exportResults",static_cast<void (toolbox_t::*)( double )>(&toolbox_t::exportResults), "export the results of the electric mechanics problem", py::arg("time"))
        ;
}


PYBIND11_MODULE(_hdgpoisson, m )
{
    using namespace Feel;
    using namespace Feel::FeelModels;

    py::enum_<MixedPoissonPhysics>(m,"MixedPoissonPhysics")
        .value("none", MixedPoissonPhysics::None)
        .value("electric", MixedPoissonPhysics::Electric)
        .value("heat", MixedPoissonPhysics::Heat);

    defMP<2,1>(m);
    defMP<2,2>(m);
    defMP<3,1>(m);
    defMP<3,2>(m);

}

