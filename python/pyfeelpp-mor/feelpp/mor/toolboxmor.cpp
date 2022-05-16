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
//!
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
// #include <pybind11/eigen.h>

#include <feel/feelmor/crb.hpp>
#include <feel/feelmor/toolboxmor.hpp>

namespace py = pybind11;
using namespace Feel;

template<typename SpaceType, int Options>
void defToolboxMor(py::module &m)
{
    using namespace Feel;
    using mor_t = ToolboxMor<SpaceType, Options>;
    static constexpr uint16_type nDim = SpaceType::nDim;
    static const bool is_time_dependent = ((Options&TimeDependent)==TimeDependent);
    static const bool is_linear = !((Options&NonLinear)==NonLinear);
    static const bool by_block = (Options&UseBlock)==UseBlock;
    using affine_decomposition_type = typename mor_t::affine_decomposition_type;

    std::string opt = "";
    if( is_time_dependent )
        opt += "_dt";
    if( !is_linear )
        opt += "_nl";
    if( by_block )
        opt += "_block";
    std::string pyclass_name = std::string("ToolboxMor_") + std::to_string(nDim) + std::string("D") + opt;
    py::class_<mor_t,std::shared_ptr<mor_t>>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,std::string const&>(),
             py::arg("name")=std::string("toolboxmor").c_str(),
             py::arg("prefix")=std::string("").c_str(),
             "Initialize the toolboxmor"
             )
        .def("initModel", &mor_t::initModel, "initialized the model" )
        .def("postInitModel", &mor_t::postInitModel, "finalize the init of the model")
        .def("setInitialized", &mor_t::setInitialized, "set the model to initialized", py::arg("b") )
        .def("functionSpace", &mor_t::functionSpace, "returns the function space" )
        .def("setFunctionSpaces",&mor_t::setFunctionSpaces, "set the function spaces", py::arg("Vh"))
        .def("setAssembleDEIM", static_cast<void (mor_t::*)(std::function<vector_ptrtype(ParameterSpaceX::Element const&)> const&) >(&mor_t::setAssembleDEIM), "set the function to assemble DEIM", py::arg("fct") )
        .def("setAssembleMDEIM", &mor_t::setAssembleMDEIM, "set the function to assemble MDEIM", py::arg("fct"))
        .def("setOnlineAssembleDEIM", &mor_t::setOnlineAssembleDEIM, "set the function to assemble DEIM for the online model", py::arg("fct"))
        .def("setOnlineAssembleMDEIM", &mor_t::setOnlineAssembleMDEIM, "set the function to assemble MDEIM for the online model", py::arg("fct"))
        .def("getDEIMReducedMesh", &mor_t::getDEIMReducedMesh, "get the reduced mesh of DEIM" )
        .def("getMDEIMReducedMesh", &mor_t::getMDEIMReducedMesh, "get the reduced mesh of MDEIM" )
        .def("parameterSpace", &mor_t::parameterSpace, "get the parameter space" )
        .def("getAffineDecomposition",
             [](mor_t& self) {
                 auto AF = self.computeAffineDecomposition();
                 if constexpr( is_time_dependent ) {
                     auto Mqm = AF.template get<0>();
                     auto Aqm = AF.template get<1>();
                     auto Fqm = AF.template get<2>();
                     return std::make_tuple(Aqm, Fqm, Mqm);
                 } else {
                     auto Aqm = AF.template get<0>();
                     auto Fqm = AF.template get<1>();
                     return std::make_tuple(Aqm, Fqm);
                 }
             })
        .def("computeBetaQm", [](mor_t& self, ParameterSpaceX::Element const& mu) {
                                 auto betaB = self.computeBetaQm(mu);
                                 if constexpr( is_time_dependent ) {
                                     auto betaMqm = betaB.template get<0>();
                                     auto betaAqm = betaB.template get<1>();
                                     auto betaFqm = betaB.template get<2>();
                                     return std::make_tuple(betaAqm, betaFqm, betaMqm);
                                 } else {
                                     auto betaAqm = betaB.template get<0>();
                                     auto betaFqm = betaB.template get<1>();
                                     return std::make_tuple(betaAqm, betaFqm);
                                 }
                              }, "compute the coefficients for parameter mu" )
        .def("computeBetaQm", [](mor_t& self, ParameterSpaceX::Element const& mu, double time) {
                                 auto betaB = self.computeBetaQm(mu, time);
                                 if constexpr( is_time_dependent ){
                                     auto betaMqm = betaB.template get<0>();
                                     auto betaAqm = betaB.template get<1>();
                                     auto betaFqm = betaB.template get<2>();
                                     return std::make_tuple(betaAqm, betaFqm, betaMqm);
                                 } else {
                                     auto betaAqm = betaB.template get<0>();
                                     auto betaFqm = betaB.template get<1>();
                                     return std::make_tuple(betaAqm, betaFqm);
                                 }
                              }, "compute the coefficients for parameter mu" )
        .def("modelProperties", &mor_t::modelProperties, "return model properties")
        ;
    std::string modelnew_name = std::string("toolboxmor_") + std::to_string(nDim) +std::string("d");
    m.def(modelnew_name.c_str(), [](std::string const& name="toolboxmor") { return std::make_shared<mor_t>(name); }," return a pointer on model");

    std::string crbmodelclass_name = std::string("CRBModel_") + pyclass_name;
    py::class_<CRBModel<mor_t>,std::shared_ptr<CRBModel<mor_t> > >(m, crbmodelclass_name.c_str())
        // .def(py::init<>(), "init")
        ;
    std::string crbmodelnew_name = std::string("crbmodel_toolboxmor_") + std::to_string(nDim) +std::string("d");
    m.def(crbmodelnew_name.c_str(), [](std::shared_ptr<mor_t>& m) { return std::make_shared<CRBModel<mor_t> >(m, crb::stage::offline); }," return a pointer on crbmodel");

    std::string crbclass_name = std::string("CRB_") + pyclass_name;
    py::class_<CRB<CRBModel<mor_t> >,std::shared_ptr<CRB<CRBModel<mor_t> > > >(m, crbclass_name.c_str())
        .def(py::init<std::string const&,
             std::shared_ptr<CRBModel<mor_t>> const&,
             crb::stage,
             std::string const&>(),
             py::arg("name"),
             py::arg("model"),
             py::arg("stage")=crb::stage::online,
             py::arg("prefixExt")=std::string(""),
             "init")
        // get rid of the return
        .def("offline", [](CRB<CRBModel<mor_t> >& c) { c.offline(); }, "run the offline step")
        .def("online", [](CRB<CRBModel<mor_t>>& c, ParameterSpaceX::Element const& mu) {
                           int timeSteps = 1, N = c.dimension();
                           std::vector<Feel::vectorN_type> uNs(timeSteps), uNolds(timeSteps), uNdus(timeSteps), uNduolds(timeSteps);
                           std::vector<double> outputs(timeSteps, 0);
                           c.fixedPoint(N, mu, uNs, uNdus, uNolds, uNduolds, outputs, 0, false);
                           return uNs[0];
                       })
        ;
    std::string crbnew_name = std::string("crb_toolboxmor_") + std::to_string(nDim) +std::string("d");
    m.def(crbnew_name.c_str(), [](std::shared_ptr<CRBModel<mor_t> >& m, std::string const& name = "toolboxmor"/*, crb::stage stage = crb::stage::online*/) { return CRB<CRBModel<mor_t> >::New(name, m, crb::stage::offline); }," return a pointer on crb");
}


PYBIND11_MODULE(_toolboxmor, m )
{
    using namespace Feel;
    using space_2d_type = FunctionSpace<Mesh<Simplex<2> >, bases<Lagrange<1, Scalar,Continuous,PointSetFekete> > >;
    using space_3d_type = FunctionSpace<Mesh<Simplex<3> >, bases<Lagrange<1, Scalar,Continuous,PointSetFekete> > >;
    

    defToolboxMor<space_2d_type, 0>(m);
    defToolboxMor<space_2d_type, 0x1>(m);
    defToolboxMor<space_3d_type, 0>(m);
    defToolboxMor<space_3d_type, 0x1>(m);

    m.def("makeToolboxMorOptions", &makeToolboxMorOptions, "get options for the model" );

}

