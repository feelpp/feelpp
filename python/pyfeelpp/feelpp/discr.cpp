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
#include <fmt/core.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelvf/mean.hpp>
#include <feel/feelvf/evaluator.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/normh1.hpp>
#include <feel/feelvf/ginac.hpp>
#include <mpi4py/mpi4py.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace Feel;

//PYBIND11_MAKE_OPAQUE(Feel::worldscomm_ptr_t);

template<typename MeshT, int Order = 1>
class MyElement: public Pch_type<MeshT,Order,double,PointSetFekete>::element_type
{
public:
    using base = typename Pch_type<MeshT,Order,double,PointSetFekete>::element_type;
    MyElement() : base() {}
};

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
    auto e = minmax( _range=elts, _pset=_Q<3>(), _expr=idv(f) );
    return std::tuple{e.min(),e.max(),e.coords()};
}

template<typename SpaceT>
void defDiscr(py::module &m, std::string const& suffix = "")
{
    using namespace Feel;
    
    using space_t = SpaceT;
    using space_ptr_t = std::shared_ptr<space_t>;
    using mesh_support_vector_t = typename space_t::mesh_support_vector_type;
    using periodicity_t = typename space_t::periodicity_type;
    using mesh_t = typename space_t::mesh_type;
    using size_type = typename mesh_t::size_type;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    using element_t = typename space_t::element_type;
    std::string pyclass_name;
    int Order = space_t::basis_0_type::nOrder;
    std::string suffix2 = fmt::format("{}D_P{}_G{}",mesh_t::nDim,Order,mesh_t::nOrder);
    //std::cout << "suffix discr: " << suffix2 << std::endl;
    //py::bind_vector<worldscomm_ptr_t>(m, "WorldsComm");
    if ( space_t::is_continuous && space_t::is_scalar )
        pyclass_name = std::string("Pch_") + suffix2;
    if ( space_t::is_continuous && space_t::is_vectorial )
        pyclass_name = std::string("Pchv_") + suffix2;
    if ( !space_t::is_continuous && space_t::is_scalar )
        pyclass_name = std::string("Pdh_") + suffix2;
    if ( !space_t::is_continuous && space_t::is_vectorial )
        pyclass_name = std::string("Pdhv_") + suffix2;
    
    py::class_<space_t,std::shared_ptr<space_t>>(m,pyclass_name.c_str())
        .def(py::init<mesh_ptr_t const&,mesh_support_vector_t const&, size_type, periodicity_t, worldscomm_ptr_t const&, std::vector<bool>>(),
             py::arg("mesh"),
             py::arg("support")=mesh_support_vector_t(),
             py::arg("components")=MESH_RENUMBER | MESH_CHECK,
             py::arg("periodicity")=periodicity_t(),
             py::arg("worldsComm"),
             py::arg("extendedDofTable") = std::vector<bool>(1,false) )
        .def("nDof",static_cast<size_type(space_t::*)() const>(&space_t::nDof), "get the number of degrees of freedom over the whole domain")
        .def("nLocalDof",static_cast<size_type(space_t::*)() const>(&space_t::nLocalDof), "get the number of degrees of freedom over the current subdomain")
        .def("nLocalDofWithGhost",static_cast<size_type(space_t::*)() const>(&space_t::nLocalDofWithGhost), "get the number of degrees of freedom over the current subdomain withthe ghost")
        .def("nLocalDofWithoutGhost",static_cast<size_type(space_t::*)() const>(&space_t::nLocalDofWithoutGhost), "get the number of degrees of freedom over the current subdomain without the ghost")
        .def("basisName",static_cast<std::string (space_t::*)() const>(&space_t::basisName), "get the basis function name")
        .def("mapPtr",&space_t::mapPtr, "return the datamap")

        .def("mesh",static_cast<mesh_ptr_t const&(space_t::*)() const>(&space_t::mesh), "get the mesh of the function space")
        .def("element",static_cast<element_t (space_t::*)(std::string const&, std::string const&)>(&space_t::element), "get an element of the function space", py::arg("name")="u", py::arg("desc")="u")
        .def("elementFromExpr",static_cast<element_t (space_t::*)(std::string const&, std::string const&, std::string const& )>(&space_t::elementFromExpr), "get an element of the function space interpolating the expression", py::arg("expr"),py::arg("name")="u", py::arg("desc")="u")
        .def("elementFromVec",static_cast<element_t (space_t::*)(std::shared_ptr<Vector<double>> const&, int)>(&space_t::element), "get an element from a vector")
        ;

    // Element
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
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( double() + py::self )
        .def( double() - py::self )
        .def( py::self + double() )
        .def( py::self - double() )
        .def( py::self * double() )
        .def( double() * py::self )
        //.def(double() / py::self)
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self *= double() )
        .def( -py::self )
        //        .def(py::self /= double())
        //        .def(py::self *= py::self)
        //        .def(py::self /= py::self)
        ;
    elt.def("printMatlab", []( element_t& element, std::string const& fname ) {
        element.printMatlab(fname);
    }, py::arg("filename"), "print element to matlab format");
    if constexpr ( space_t::is_scalar )
    {
        elt.def( "on", static_cast<void ( element_t::* )( elements_reference_wrapper_t<mesh_ptr_t> const&, Expr<GinacEx<2>> const&, std::string const&, GeomapStrategyType, bool, bool )>( &element_t::template onImpl<elements_reference_wrapper_t<mesh_ptr_t>, Expr<GinacEx<2>>> ),
                 py::arg( "range" ), py::arg( "expr" ), py::arg( "prefix" ) = "",
                 py::arg( "geomap" ) = GeomapStrategyType::GEOMAP_OPT, py::arg( "accumulate" ) = false, py::arg( "verbose" ) = false, "build the interpolant of the expression expr on a range of elements" );
    }
    elt.def( "on", []( element_t& element, elements_reference_wrapper_t<mesh_ptr_t> const& r, 
                        Expr<GinacMatrix<element_t::nComponents1,element_t::nComponents2,2>> const& e, std::string const& p, GeomapStrategyType g, bool a, bool v ){
                            element.on( _range=r, _expr=e );
                    },
            py::arg( "range" ), py::arg( "expr" ), py::arg( "prefix" ) = "",
            py::arg( "geomap" ) = GeomapStrategyType::GEOMAP_OPT, py::arg( "accumulate" ) = false, py::arg( "verbose" ) = false, "build the interpolant of the expression expr on a range of elements" );
#if 0
    std::string I_pyclass_name = std::string("I_") + pyclass_name;
    py::class_<I_t<space_t,space_t>,std::shared_ptr<I_t<space_t,space_t>>>(m,I_pyclass_name.c_str())
        .def(py::init<>())
        //.def(py::init<std::shared_ptr<space_t> const&, std::shared_ptr<space_t> const&>())
        ;
#endif        
    m.def( "normL2", static_cast<double (*)( elements_reference_wrapper_t<mesh_ptr_t> const&,element_t const&)>( &f_norml2<elements_reference_wrapper_t<mesh_ptr_t> const&,element_t const&> ), "compute L2 norm of function over a range of elements", py::arg("range"), py::arg("expr") );
    m.def( "normH1", static_cast<double (*)( elements_reference_wrapper_t<mesh_ptr_t> const&,element_t const&)>( &f_normh1<elements_reference_wrapper_t<mesh_ptr_t> const&,element_t const&> ), "compute H1 norm of function over a range of elements", py::arg("range"), py::arg("expr") );
    m.def( "mean", static_cast<eigen_v_t (*)( elements_reference_wrapper_t<mesh_ptr_t> const&,element_t const&)>( &f_mean<elements_reference_wrapper_t<mesh_ptr_t> const&,element_t const&> ), "compute mean of function over a range of elements", py::arg("range"), py::arg("expr") );
    m.def( "mean", static_cast<eigen_v_t ( * )( faces_reference_wrapper_t<mesh_ptr_t> const&, element_t const& )>( &f_mean<faces_reference_wrapper_t<mesh_ptr_t> const&, element_t const&> ), "compute mean of function over a range of facets", py::arg( "range" ), py::arg( "expr" ) );
    m.def( "minmax", static_cast<std::tuple<double,double,eigen_v2_t> (*)( elements_reference_wrapper_t<mesh_ptr_t> const&,element_t const&)>( &f_minmax<elements_reference_wrapper_t<mesh_ptr_t> const&,element_t const&> ), "compute min max argmin argmax of a function over a range of elements", py::arg("range"), py::arg("expr") );

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

    auto dimt = hana::make_tuple( hana::int_c<1>, hana::int_c<2>, hana::int_c<3>);
    auto ordert = hana::make_tuple( hana::int_c<0>, hana::int_c<1>, hana::int_c<2>, hana::int_c<3> );
    auto geot = hana::make_tuple( hana::int_c<1>, hana::int_c<2> );
    hana::for_each(hana::cartesian_product(hana::make_tuple(dimt, ordert, geot)),
                   [&m]( auto const& d )
                       {
                           constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                           constexpr int _order = std::decay_t<decltype(hana::at_c<1>(d))>::value;
                           constexpr int _geo = std::decay_t<decltype( hana::at_c<2>( d ) )>::value;
                           defDiscr<Pch_type<Mesh<Simplex<_dim, _geo>>, _order>>( m );
                           defDiscr<Pdh_type<Mesh<Simplex<_dim, _geo>>, _order>>( m );
                           defDiscr<Pchv_type<Mesh<Simplex<_dim, _geo>>, _order>>( m );
                           defDiscr<Pdhv_type<Mesh<Simplex<_dim, _geo>>, _order>>( m );
                       });
    defDiscrDiscontinuous<Pdh_type<Mesh<Simplex<2>>,0>>( m );
    defDiscrDiscontinuous<Pdh_type<Mesh<Simplex<3>>,0>>( m );
}

