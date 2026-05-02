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
#pragma once

#include <fmt/core.h>
#include <feel/feelpython/pybind11/pybind11.h>
#include <feel/feelpython/pybind11/eigen.h>
#include <pybind11/operators.h>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feelvf/mean.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/normh1.hpp>
#include <feel/feelvf/evaluator.hpp>
#include <feel/feelvf/ginac.hpp>

namespace py = pybind11;
using namespace Feel;

PYBIND11_MAKE_OPAQUE(std::vector<DofTableExtendedType>);

template<typename RangeT, typename FunctionT>
double
f_norml2( RangeT const& elts, FunctionT const& f )
{
    return normL2( _range=elts, _expr=idv(f) );
}
template<typename RangeT, typename FunctionT>
double
f_normh1( RangeT const& elts, FunctionT const& f )
{
    return normH1( _range=elts, _expr=idv(f), _grad_expr=gradv(f) );
}
using eigen_v_t = Eigen::Matrix<double,Eigen::Dynamic,1>;
template<typename RangeT, typename FunctionT>
eigen_v_t
f_mean( RangeT const& elts, FunctionT const& f )
{
    return mean( _range=elts, _expr=idv(f) );
}
using eigen_v2_t = Eigen::Matrix<double,Eigen::Dynamic,2>;
template<typename RangeT, typename FunctionT>
std::tuple<double,double,eigen_v2_t>
f_minmax( RangeT const& elts, FunctionT const& f )
{
    auto e = Feel::vf::minmax( _range=elts, _pset=_Q<3>(), _expr=idv(f) );
    return std::tuple{e.min(),e.max(),e.coords()};
}

template<typename SpaceT, int SpaceOrder>
void
defDiscr( py::module& m, std::string const& family_prefix = "" )
{
    using namespace Feel;

    using space_t = SpaceT;
    using space_ptr_t = std::shared_ptr<space_t>;
    using mesh_support_vector_t = typename space_t::mesh_support_vector_type;
    using mesh_t = typename space_t::mesh_type;
    using size_type = typename mesh_t::size_type;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    using element_t = typename space_t::element_type;
    constexpr int BasisOrder = space_t::basis_0_type::nOrder;
    std::string order_label = ( SpaceOrder == Dynamic ) ? "Dynamic" : std::to_string( BasisOrder );
    std::string family = family_prefix;
    if ( family.empty() )
    {
        if ( space_t::is_continuous && space_t::is_scalar )
            family = "Pch";
        else if ( space_t::is_continuous && space_t::is_vectorial )
            family = "Pchv";
        else if ( !space_t::is_continuous && space_t::is_scalar )
            family = "Pdh";
        else if ( !space_t::is_continuous && space_t::is_vectorial )
            family = "Pdhv";
        else
            family = "Space";
    }
    std::string pyclass_name = fmt::format( "{}_{}D_P{}_G{}", family, mesh_t::nDim, order_label, mesh_t::nOrder );
    if ( py::hasattr( m, pyclass_name.c_str() ) )
    {
        // Some template combinations can map to the same exported Python class name.
        return;
    }

    py::class_<space_t,std::shared_ptr<space_t>>(m,pyclass_name.c_str())
        .def(py::init([]( mesh_ptr_t const& mesh,
                          mesh_support_vector_t const& support,
                          size_type components,
                          worldscomm_ptr_t const& worldsComm,
                          std::vector<DofTableExtendedType> extendedDofTable,
                          int runtimeOrder )
            {
                return std::make_shared<space_t>( mesh,
                                                  support,
                                                  components,
                                                  worldsComm,
                                                  std::move( extendedDofTable ),
                                                  RuntimeOrder::checked( runtimeOrder ) );
            }),
             py::arg("mesh"),
             py::arg("support")=mesh_support_vector_t(),
             py::arg("components")=MESH_RENUMBER | MESH_CHECK,
             py::arg("worldsComm"),
             py::arg("extendedDofTable") = std::vector<DofTableExtendedType>(space_t::nSpaces,DofTableExtendedType::DEFAULT),
             py::arg("runtimeOrder") = 1
             )
        .def("nDof",static_cast<size_type(space_t::*)() const>(&space_t::nDof), "get the number of degrees of freedom over the whole domain")
        .def("nLocalDof",static_cast<size_type(space_t::*)() const>(&space_t::nLocalDof), "get the number of degrees of freedom over the current subdomain")
        .def("nLocalDofWithGhost",static_cast<size_type(space_t::*)() const>(&space_t::nLocalDofWithGhost), "get the number of degrees of freedom over the current subdomain withthe ghost")
        .def("nLocalDofWithoutGhost",static_cast<size_type(space_t::*)() const>(&space_t::nLocalDofWithoutGhost), "get the number of degrees of freedom over the current subdomain without the ghost")
        .def("basisName",static_cast<std::string (space_t::*)() const>(&space_t::basisName), "get the basis function name")
        .def("order",&space_t::order, "get the polynomial order (runtime order for dynamic spaces)")
        .def("mapPtr",&space_t::mapPtr, "return the datamap")

        .def("mesh",static_cast<mesh_ptr_t const&(space_t::*)() const>(&space_t::mesh), "get the mesh of the function space")
        .def("element",static_cast<element_t (space_t::*)(std::string const&, std::string const&)>(&space_t::element), "get an element of the function space", py::arg("name")="u", py::arg("desc")="u")
        .def("elementFromExpr",static_cast<element_t (space_t::*)(std::string const&, std::string const&, std::string const& )>(&space_t::elementFromExpr), "get an element of the function space interpolating the expression", py::arg("expr"),py::arg("name")="u", py::arg("desc")="u")
        .def("element", []( std::shared_ptr<space_t> & Xh, Vector<double> const& v, int blockIdStart ) { return Xh->element( v );
            }, py::arg("vec"), py::arg("start") = 0, "get an element from a vector")
        .def("element", []( std::shared_ptr<space_t> & Xh, VectorPetsc<double> const& v, int blockIdStart ) { return Xh->element( v, blockIdStart );
            }, py::arg("vec"), py::arg("start") = 0, "get an element from a vector")
        ;

    std::string e_pyclass_name = std::string("Element_") + pyclass_name;
    py::class_<element_t,std::shared_ptr<element_t>,VectorUblas<double>> elt(m,e_pyclass_name.c_str());
    elt.def( py::init<>() )
        .def( py::init<std::shared_ptr<space_t> const&, std::string const&, std::string const&, size_type, ComponentType>(), py::arg( "space" ), py::arg( "name" ), py::arg( "desc" ), py::arg( "start" ) = 0, py::arg( "ct" ) = ComponentType::NO_COMPONENT )
        .def( "functionSpace", static_cast<space_ptr_t const& (element_t::*)() const>( &element_t::functionSpace ), "Get funtion space from element" )
        .def( "size", static_cast<size_type ( element_t::* )() const>( &element_t::size ), "Get size of element" )
        .def( "min", static_cast<double ( element_t::* )() const>( &element_t::min ), "get the minimum of the element vector representation" )
        .def( "max", static_cast<double ( element_t::* )() const>( &element_t::max ), "get the maximum of the element vector representation" )
        .def( "save", &element_t::saveImpl, py::arg( "path" ), py::arg( "name" ), py::arg( "type" ) = "default", py::arg( "suffix" ) = "", py::arg( "sep" ) = "", "save functionspace element in file " )
        .def( "load", &element_t::loadImpl, py::arg( "path" ), py::arg( "name" ), py::arg( "type" ) = "default", py::arg( "suffix" ) = "", py::arg( "sep" ) = "", py::arg("space_path") = "", "load functionspace element from file " )
        .def( pybind11::detail::self + pybind11::detail::self )
        .def( pybind11::detail::self - pybind11::detail::self )
        .def( double() + pybind11::detail::self )
        .def( double() - pybind11::detail::self )
        .def( pybind11::detail::self + double() )
        .def( pybind11::detail::self - double() )
        .def( pybind11::detail::self * double() )
        .def( double() * pybind11::detail::self )
        .def( pybind11::detail::self += pybind11::detail::self )
        .def( pybind11::detail::self -= pybind11::detail::self )
        .def( pybind11::detail::self *= double() )
        .def( -pybind11::detail::self );

    elt.def("printMatlab", []( element_t& element, std::string const& fname ) {
        element.printMatlab(fname);
    }, py::arg("filename"), "print element to matlab format");

    if constexpr ( space_t::is_scalar )
    {
        elt.def( "on", static_cast<void ( element_t::* )( Range<mesh_ptr_t,MESH_ELEMENTS> const&, Expr<GinacEx<2>> const&, std::string const&, GeomapStrategyType, bool, bool )>( &element_t::template onImpl<Range<mesh_ptr_t,MESH_ELEMENTS>, Expr<GinacEx<2>>> ),
                 py::arg( "range" ), py::arg( "expr" ), py::arg( "prefix" ) = "",
                 py::arg( "geomap" ) = GeomapStrategyType::GEOMAP_OPT, py::arg( "accumulate" ) = false, py::arg( "verbose" ) = false, "build the interpolant of the expression expr on a range of elements" );
    }
    elt.def( "on", []( element_t& element, Range<mesh_ptr_t,MESH_ELEMENTS> const& r,
                        Expr<GinacMatrix<element_t::nComponents1,element_t::nComponents2,2>> const& e, std::string const& p, GeomapStrategyType g, bool a, bool v ){
                            element.on( _range=r, _expr=e );
                    },
            py::arg( "range" ), py::arg( "expr" ), py::arg( "prefix" ) = "",
            py::arg( "geomap" ) = GeomapStrategyType::GEOMAP_OPT, py::arg( "accumulate" ) = false, py::arg( "verbose" ) = false, "build the interpolant of the expression expr on a range of elements" );

    m.def( "normL2", static_cast<double (*)( Range<mesh_ptr_t,MESH_ELEMENTS> const&,element_t const&)>( &f_norml2<Range<mesh_ptr_t,MESH_ELEMENTS> const&,element_t const&> ), "compute L2 norm of function over a range of elements", py::arg("range"), py::arg("expr") );
    m.def( "normH1", static_cast<double (*)( Range<mesh_ptr_t,MESH_ELEMENTS> const&,element_t const&)>( &f_normh1<Range<mesh_ptr_t,MESH_ELEMENTS> const&,element_t const&> ), "compute H1 norm of function over a range of elements", py::arg("range"), py::arg("expr") );
    m.def( "mean", static_cast<eigen_v_t (*)( Range<mesh_ptr_t,MESH_ELEMENTS> const&,element_t const&)>( &f_mean<Range<mesh_ptr_t,MESH_ELEMENTS> const&,element_t const&> ), "compute mean of function over a range of elements", py::arg("range"), py::arg("expr") );
    m.def( "mean", static_cast<eigen_v_t ( * )( Range<mesh_ptr_t,MESH_FACES> const&, element_t const& )>( &f_mean<Range<mesh_ptr_t,MESH_FACES> const&, element_t const&> ), "compute mean of function over a range of facets", py::arg( "range" ), py::arg( "expr" ) );
    m.def( "minmax", static_cast<std::tuple<double,double,eigen_v2_t> (*)( Range<mesh_ptr_t,MESH_ELEMENTS> const&,element_t const&)>( &f_minmax<Range<mesh_ptr_t,MESH_ELEMENTS> const&,element_t const&> ), "compute min max argmin argmax of a function over a range of elements", py::arg("range"), py::arg("expr") );
}

template<typename space_t>
void
defDiscrDiscontinuous( py::module& m )
{
    if ( (space_t::is_continuous == false) && ( space_t::basis_0_type::nOrder == 0 ) )
    {
        m.def( "pid", &regionProcess<space_t>, "get an piecewise constant function storing the process ids", py::arg("space") );
    }
}

void bindDiscrCommon( py::module& m );
void bindDiscrPch( py::module& m );
void bindDiscrPchv( py::module& m );
void bindDiscrPdh( py::module& m );
void bindDiscrPdhv( py::module& m );
void bindDiscrDh( py::module& m );
