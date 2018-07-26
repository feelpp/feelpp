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

#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <mpi4py/mpi4py.h>

namespace py = pybind11;
using namespace Feel;

template<typename MeshT>
std::shared_ptr<MeshT>
loadmesh( std::string const& n, double h )
{
    return loadMesh( _mesh=new MeshT, _filename=n, _h=h );
}
template<int nDim, int OrderVelocity=2, int OrderPressure = 1, int OrderGeo = 1>
void defFM(py::module &m)
{
    using namespace Feel;
    using fm_t = FeelModels::FluidMechanics< Simplex<nDim,OrderGeo>,
                                             Lagrange<OrderVelocity, Vectorial,Continuous,PointSetFekete>,
                                             Lagrange<OrderPressure, Scalar,Continuous,PointSetFekete> > ;
    std::string pyclass_name = std::string("Fluid_P") + std::to_string(OrderVelocity) + std::string("P") + std::to_string(OrderPressure) + std::string("G") + std::to_string(OrderGeo);
    py::class_<fm_t,std::shared_ptr<fm_t>>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,bool,WorldComm const&,std::string const&, ModelBaseRepository const&>(),
             py::arg(_buildmesh)=true,
             py::arg(_worldComm)=Environment::worldComm(),
             py::arg(_subprefix)="",
             py::arg(_modelRep) = ModelBaseRepository() 
             );
        
}
    

PYBIND11_MODULE(fluid, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;

    
    defFM<2,2,1,1>(m);

}

