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

#include <feel/feelmesh/filters.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/meshstructured.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <mpi4py/mpi4py.h>

#include <boost/shared_ptr.hpp>
#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/binding.hpp>
//#include <boost/parameter/python.hpp>
#include <boost/mpl/vector.hpp>

#include<feel/feelcore/environment.hpp>
#include <feel/feelmesh/remesh.hpp>
namespace py = pybind11;

using namespace Feel;

template<typename MeshT>
std::shared_ptr<MeshT>
loadmesh( std::shared_ptr<MeshT> const& m, std::string const& n, double h )
{
    
    return loadMesh( _mesh=new MeshT, _filename=n, _h=h );
}

template<typename MeshT>
elements_pid_t<MeshT>
elementsByPid( MeshT m ) 
{
    return elements( m );
}

template<typename MeshT>
markedelements_t<MeshT>
elementsByMarker( MeshT m, std::string const& tag ) 
{
    return markedelements( m, tag );
}

template<typename MeshT>
boundaryfaces_t<MeshT>
boundaryfacesByPid( MeshT m ) 
{
    return boundaryfaces( m );
}

template<typename MeshT>
auto
nElementsTuple( elements_pid_t<MeshT> const& r, bool global ) 
{
    return nelements( r, global );
}

template<typename MeshT>
auto
nFacesTuple( faces_reference_wrapper_t<MeshT> const& r, bool global ) 
{
    return nelements( r, global );
}

template<typename MeshT>
void defMesh(py::module &m)
{
    using namespace Feel;
    using mesh_t = MeshT;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    int Dim = MeshT::nDim;
    int Order = MeshT::nOrder;
    int RealDim = MeshT::nRealDim;
    std::string pyclass_name = std::string("Mesh_");
    std::string suffix = "S";
    if constexpr ( is_hypercube_v<typename MeshT::element_type> )
        suffix = std::string("H");
    suffix += std::to_string(Dim)+"DG"+ std::to_string(Order) + std::string("R") + std::to_string(RealDim);
    pyclass_name += suffix;
    py::class_<mesh_t,std::shared_ptr<mesh_t>>(m,pyclass_name.c_str())
        .def(py::init<worldcomm_ptr_t const&>(),py::arg("worldComm"),"Construct a new mesh")
        .def_static("create",&mesh_t::New,"Construct a new shared_ptr mesh")
        .def("dimension",&mesh_t::dimension,"get topological dimension")
        .def("realDimension",&mesh_t::realDimension,"get real dimension")


        // elements counts
        .def("numGlobalElements",&mesh_t::numGlobalElements,"get the number of elements over the whole mesh, requires communication if the mesh is parallel")
        .def("numGlobalFaces",&mesh_t::numGlobalFaces,"get the number of faces over the whole mesh, requires communication if the mesh is parallel")
        .def("numGlobalEdges",&mesh_t::numGlobalEdges,"get the number of edges over the whole mesh, requires communication if the mesh is parallel") 
        .def("numGlobalPoints",&mesh_t::numGlobalPoints,"get the number of points over the whole mesh, requires communication if the mesh is parallel") 
        // measures
        .def("hAverage",&mesh_t::hAverage,"get the average edge length of the mesh")
        .def("hMin",&mesh_t::hMin,"get the minimum edge length of the mesh")
        .def("hMax",&mesh_t::hMax,"get the maximum edge length of the mesh")
        .def("updateMeasures",&mesh_t::updateMeasures,"update the measures of the mesh")
        .def("measure",&mesh_t::measure,py::arg("parallel") = true,"get the measure of the mesh")
        .def("measureBoundary",&mesh_t::measureBoundary,"get the measure of the boundary of the mesh");

    pyclass_name = std::string("simplex_elements_reference_wrapper_")+suffix;
    py::class_<elements_reference_wrapper_t<mesh_ptr_t>>(m,pyclass_name.c_str())
        .def(py::init<>());
    pyclass_name = std::string("faces_reference_wrapper_")+suffix;
    py::class_<faces_reference_wrapper_t<mesh_ptr_t>>(m,pyclass_name.c_str())
        .def(py::init<>());

    pyclass_name = std::string("mesh_support_")+suffix;
    py::class_<MeshSupport<mesh_t>,std::shared_ptr<MeshSupport<mesh_t>>>(m,pyclass_name.c_str())
        .def(py::init<mesh_ptr_t const&>(),py::arg("mesh"))
        .def(py::init<mesh_ptr_t const&, elements_reference_wrapper_t<mesh_t> const&>(),py::arg("mesh"), py::arg("range"))
        ;

    pyclass_name = std::string("mesh_support_vector_")+suffix;
    py::class_<fusion::vector<std::shared_ptr<MeshSupport<mesh_t>>>>(m,pyclass_name.c_str())
        .def(py::init<>());
        

    // load mesh
    m.def("load",&loadmesh<mesh_t>,"load a mesh from a file");

    m.def("elements", &elementsByPid<mesh_ptr_t>,"get iterator over the elements of the mesh", py::arg("mesh"));
    m.def("markedelements", &elementsByMarker<mesh_ptr_t>,"get iterator over the marked elements of the mesh", py::arg("mesh"),py::arg("tag"));
    m.def("boundaryfaces", &boundaryfacesByPid<mesh_ptr_t>,"get iterator over the boundary faces of the mesh", py::arg("mesh"));
    m.def("nelements", &nElementsTuple<mesh_ptr_t>,"get the number of elements in range, the local one if global is false", py::arg("range"),py::arg("global") = true );
    m.def("nfaces", &nFacesTuple<mesh_ptr_t>,"get the number of faces in range, the local one if global is false", py::arg("range"),py::arg("global") = true );
    //m.def("markedfaces", &markedfaces<mesh_t>,"get iterator over the marked faces of the mesh");
    m.def( "markedelements", []( mesh_ptr_t const& r, std::string const& m ) {
            return markedelements(r,m);
        },py::arg("range"), py::arg("marker"), "return the range of elements of the mesh with marker" );
    m.def( "markedelements", []( mesh_ptr_t const& r, std::vector<std::string> const& m ) {
            return markedelements(r,m);
        },py::arg("range"), py::arg("markers"), "return the range of elements of the mesh with markers" );


    m.def(
        "markedfaces", []( mesh_ptr_t const& r, std::string const& m ) {
            return markedfaces( r, m );
        },
        py::arg( "range" ), py::arg( "marker" ), "return the range of facets of the mesh with marker" );
    m.def(
        "markedfaces", []( mesh_ptr_t const& r, std::vector<std::string> const& m ) {
            return markedfaces( r, m );
        },
        py::arg( "range" ), py::arg( "markers" ), "return the range of facets of the mesh with marker" );
    if constexpr ( mesh_t::nDim == 3 )
    {
        m.def(
            "markededges", []( mesh_ptr_t const& r, std::string const& m ) {
                return markededges( r, m );
            },
            py::arg( "range" ), py::arg( "marker" ), "return the range of facets of the mesh with marker" );

    }
    if constexpr ( mesh_t::nDim >= 2 )
    {
        pyclass_name = std::string("remesher_")+suffix;
        py::class_<Remesh<mesh_t>, std::shared_ptr<Remesh<mesh_t>>>( m, pyclass_name.c_str() )
            .def( py::init<mesh_ptr_t const&, std::vector<std::string> const&, std::vector<std::string> const&>(),
                  py::arg( "mesh" ), py::arg( "required_elts" ) = std::vector<std::string>{}, py::arg( "required_facets" ) = std::vector<std::string>{}, "Initialize a remesher" )
            .def( "execute", &Remesh<mesh_t>::execute, py::arg( "run" ) = true, "execute remesh task" )
            .def( "setMetric", &Remesh<mesh_t>::setMetric, py::arg( "metric" ), "set the metric" );
        m.def(
            "remesher", []( mesh_ptr_t const& r ) {
                return std::make_shared<Remesh<mesh_t>>( r );
            },
            py::arg( "mesh" ), "create a Remesher data structure", py::return_value_policy::copy );
        m.def(
            "remesher", []( mesh_ptr_t const& r, std::vector<std::string> const& req_elts ) {
                return std::make_shared<Remesh<mesh_t>>( r, req_elts );
            },
            py::return_value_policy::copy,py::arg( "mesh" ), py::arg( "required_elts" ), "create a Remesher data structure" );
        m.def(
            "remesher", []( mesh_ptr_t const& r, std::vector<std::string> const& req_elts, std::vector<std::string> const& req_facets ) {
                return std::make_shared<Remesh<mesh_t>>(  r, req_elts, req_facets );
            },
            py::return_value_policy::copy,py::arg( "mesh" ), py::arg( "required_elts" )=std::vector<std::string>{},  py::arg( "required_facets" )=std::vector<std::string>{}, 
            "create a Remesher data structure" );
    }
}
    

PYBIND11_MODULE(_mesh, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;

    py::enum_<GeomapStrategyType>(m,"GeomapStrategyType")
        .value("GEOMAP_OPT", GeomapStrategyType::GEOMAP_OPT)
        .value("GEOMAP_O1", GeomapStrategyType::GEOMAP_O1)
        .value("GEOMAP_HO", GeomapStrategyType::GEOMAP_HO)
        .export_values();

    // 1D
    defMesh<Mesh<Simplex<1,1,1>>>(m);
    defMesh<Mesh<Simplex<1,2,1>>>(m);
    // 2D
    defMesh<Mesh<Simplex<2,1,2>>>(m);
    defMesh<Mesh<Simplex<2,2,2>>>(m);
    //defMesh<Simplex<1,1,2>>(m);
    //defMesh<Simplex<1,2,2>>(m);
    // 3D
    defMesh<Mesh<Simplex<3,1,3>>>(m);
    defMesh<Mesh<Simplex<3,2,3>>>(m);
    //defMesh<Simplex<2,1,3>>(m);
    //defMesh<Simplex<2,2,3>>(m);
    //defMesh<Simplex<1,1,3>>(m);
    //defMesh<Simplex<1,2,3>>(m);
    defMesh<MeshStructured<Hypercube<2>>>(m);
}
