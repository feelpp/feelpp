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

#include <feel/feel.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <mpi4py/mpi4py.h>

#include <boost/shared_ptr.hpp>
#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/binding.hpp>
#include <boost/parameter/python.hpp>
#include <boost/python.hpp>
#include <boost/mpl/vector.hpp>

#include<feel/feelcore/environment.hpp>

namespace Feel {


using namespace Feel;

namespace py = boost::parameter::python;

FEELPP_NO_EXPORT void export_core()
{
    if (import_mpi4py()<0) return ;

    class_<Feel::Environment,boost::noncopyable>("Environment", init<boost::python::list>()) 
        .def("worldComm",&Feel::Environment::worldComm,return_value_policy<copy_non_const_reference>())
        .staticmethod("worldComm");

    class_<WorldComm>("WorldComm",init<>());
}

}
