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

#include <feel/feelfilters/exporter.hpp>
#include <mpi4py/mpi4py.h>

namespace py = pybind11;
using namespace Feel;

template<typename MeshT, int Order = 1>
std::shared_ptr<Exporter<MeshT,Order>> exporter_impl( std::shared_ptr<MeshT> const& m,
                                                      std::string const& n,
                                                      std::string const& geo )
{
    using namespace Feel;
    
    using exporter_t = Exporter<MeshT,Order>;
    using exporter_ptr_t = std::shared_ptr<exporter_t>;
    using mesh_ptr_t = std::shared_ptr<MeshT>;
    return exporter( _mesh=m, _name=n, _geo=geo );
}
template<typename MeshT, int Order = 1>
class PyExporter : public Exporter<MeshT, Order>
{
public:
    using base_t = Exporter<MeshT,Order>;
    using base_t::Exporter;
    void save() const override
        {
            PYBIND11_OVERLOAD_PURE(void,base_t,save);
        }
    void visit( MeshT* m ) override
        {
            PYBIND11_OVERLOAD_PURE(void,base_t,visit,m);
        }
};
template<typename MeshT, int Order = 1>
void defExporter(py::module &m)
{
    using namespace Feel;
    
    using exporter_t = Exporter<MeshT,Order>;
    using exporter_ptr_t = std::shared_ptr<exporter_t>;
    using mesh_ptr_t = std::shared_ptr<MeshT>;
    std::string pyclass_name;

    std::string suffix = std::to_string(MeshT::nDim)+"D";
    pyclass_name = std::string("Exporter") + suffix;
    py::class_<exporter_t,PyExporter<MeshT,1>,exporter_ptr_t>(m,pyclass_name.c_str())
        .def(py::init<>() )
        .def("save", &exporter_t::save,"save Exporter data set");
        //.def(py::init(static_cast<exporter_ptr_t (exporter_t::*)( std::string, WorldComm const&)>(&exporter_t::New)), py::arg("prefix"), py::arg("worldcomm")=Environment::worldComm() );;

    m.def( "exporter",
           &exporter_impl<MeshT,Order>,
           py::arg("mesh"),
           py::arg("name") = "Exporter",
           py::arg("geo") = "change_coords_only",
           "create an Exporter object" );
        
}
    
PYBIND11_MODULE(exporter, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;
    
    defExporter<Mesh<Simplex<2>>,1>( m );
}

