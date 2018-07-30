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

#include <feel/feel.hpp>
#include <feel/feelmesh/filters.hpp>
#include <mpi4py/mpi4py.h>

#include <boost/shared_ptr.hpp>
#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/binding.hpp>
//#include <boost/parameter/python.hpp>
#include <boost/mpl/vector.hpp>

#include<feel/feelcore/environment.hpp>
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
size_type
nElementsTuple( elements_pid_t<MeshT> const& r, bool global ) 
{
    return nelements( r, global );
}

template<typename MeshT>
size_type
nFacesTuple( faces_reference_wrapper_t<MeshT> const& r, bool global ) 
{
    return nelements( r, global );
}

template<typename ConvexT>
void defMesh(py::module &m)
{
    using namespace Feel;
    using mesh_t = Mesh<ConvexT>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    int Dim = ConvexT::nDim;
    int Order = ConvexT::nOrder;
    int RealDim = ConvexT::nRealDim;
    std::string pyclass_name = std::string("Mesh_");
    std::string suffix = std::to_string(Dim);
    if ( Order > 1 )
        suffix += std::string("_") + std::to_string(Order);
    if ( RealDim != Dim && Order == 1 )
        suffix += std::string("_") + std::to_string(Order) + std::string("_") + std::to_string(RealDim);
    else if ( RealDim != Dim && Order == 1 )
        suffix += std::string("_") + std::to_string(RealDim);
    pyclass_name += suffix;
    py::class_<mesh_t,std::shared_ptr<mesh_t>>(m,pyclass_name.c_str())
        .def(py::init<>())
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

    pyclass_name = std::string("elements_reference_wrapper_")+suffix;
    py::class_<elements_reference_wrapper_t<mesh_ptr_t>>(m,pyclass_name.c_str())
        .def(py::init<>());
    pyclass_name = std::string("faces_reference_wrapper_")+suffix;
    py::class_<faces_reference_wrapper_t<mesh_ptr_t>>(m,pyclass_name.c_str())
        .def(py::init<>());

    // load mesh
    m.def("load",&loadmesh<mesh_t>,"load a mesh from a file");

    m.def("elements", &elementsByPid<mesh_ptr_t>,"get iterator over the elements of the mesh", py::arg("mesh"));
    m.def("markedelements", &elementsByMarker<mesh_ptr_t>,"get iterator over the marked elements of the mesh", py::arg("mesh"),py::arg("tag"));
    m.def("boundaryfaces", &boundaryfacesByPid<mesh_ptr_t>,"get iterator over the boundary faces of the mesh", py::arg("mesh"));
    m.def("nelements", &nElementsTuple<mesh_ptr_t>,"get the number of elements in range, the local one if global is false", py::arg("range"),py::arg("global") = true );
    m.def("nfaces", &nFacesTuple<mesh_ptr_t>,"get the number of faces in range, the local one if global is false", py::arg("range"),py::arg("global") = true );
    //m.def("markedfaces", &markedfaces<mesh_t>,"get iterator over the marked faces of the mesh");
        
}
    

PYBIND11_MODULE(mesh, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;

    defMesh<Simplex<1>>(m);
    defMesh<Simplex<2>>(m);
    //defMesh<3>();
}
