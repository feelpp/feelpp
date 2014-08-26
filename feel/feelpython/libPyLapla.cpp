/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Thomas Lantz <lantz.thomas0@gmail.com>

  Date: 2014-08-26

  Copyright (C)  2014 Université de Strasborg

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file libPyLapla.cpp
   \author Thomas Lantz <lantz.thomas0@gmail.com>
   \date 2014-08-26

*/

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
#include<feel/feelfilters/loadmesh.hpp>
#include<feel/feelfilters/exporter.hpp>
#include<feel/feelfilters/detail/mesh.hpp>
#include<feel/feeldiscr/functionspace.hpp>

using namespace boost::python;
using namespace Feel;

namespace py = boost::parameter::python;

// retrieval of complicated return types with the "give bad type to have the good one" trick 

typedef Feel::vf::detail::LinearForm<Feel::FunctionSpace<Feel::Mesh<Feel::Simplex<2, 1, 2>, double, 0>,Feel::bases<Feel::Lagrange<1, Scalar, Feel::Continuous, PointSetEquiSpaced, 0> >, double,Feel::Periodicity<Feel::NoPeriodicity>, Feel::mortars<Feel::NoMortar> >, Feel::Vector<double>,Feel::Vector<double> > form1_return_type;

typedef Feel::vf::detail::BilinearForm<Feel::FunctionSpace<Feel::Mesh<Feel::Simplex<2, 1, 2>, double, 0>,Feel::bases<Feel::Lagrange<1, Scalar, Feel::Continuous, PointSetEquiSpaced, 0> >, double,Feel::Periodicity<Feel::NoPeriodicity>, Feel::mortars<Feel::NoMortar> >,Feel::FunctionSpace<Feel::Mesh<Feel::Simplex<2, 1, 2>, double, 0>, Feel::bases<Feel::Lagrange<1, Scalar,Feel::Continuous, PointSetEquiSpaced, 0> >, double, Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar> >, Feel::VectorUblas<double, ublas::vector<double> > > form2_return_type;

typedef FunctionSpace<Mesh<Simplex<2>>,Feel::bases<Feel::Lagrange<1,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>> function_space_type;

//definition of all methods we need for the wrapping and wich are create from a BOOST_PARAMETER_FUNCTION or a method with default value arguments  

boost::shared_ptr<Mesh<Simplex<2>>> unitSquare_w ()
{
    return unitSquare();
}


function_space_type::element_type element_w(boost::shared_ptr<Feel::meta::Pch<Mesh<Simplex<2>>,1>::type> f)
    {
        return f->element("u");

            }
//form 

form1_return_type form1_w (boost::shared_ptr<Feel::meta::Pch<Mesh<Simplex<2>>,1>::type> f)
{
    return form1(_test=f);
}

form2_return_type form2_w (boost::shared_ptr<Feel::meta::Pch<Mesh<Simplex<2>>,1>::type> f,boost::shared_ptr<Feel::meta::Pch<Mesh<Simplex<2>>,1>::type> f2)
{
    return form2(_trial=f,_test=f2);
}

//integrate

form1_return_type integrate_form1 (form1_return_type f,boost::shared_ptr<Mesh<Feel::Simplex<2>>> mesh, function_space_type::element_type v)
{
  f=integrate(_range=elements(mesh),_expr=id(v));
  return f;
}


form2_return_type integrate_form2 (form2_return_type f,boost::shared_ptr<Mesh<Feel::Simplex<2>>> mesh,function_space_type::element_type u,function_space_type::element_type v)
{
 f=integrate(_range=elements(mesh),_expr=gradt(u)*trans(grad(v)));
 return f;
}

form2_return_type on_form2 (form2_return_type f,form1_return_type l, boost::shared_ptr<Mesh<Simplex<2>>> mesh,function_space_type::element_type u)
{
 f+=on(_range=boundaryfaces(mesh),_rhs=l,_element=u,_expr=constant(0.));
 return f;
}

//solve

 function_space_type::element_type solve_w(form2_return_type a,form1_return_type l,function_space_type::element_type u)
{
 a.solve(_rhs=l,_solution=u);
 return u;
}

//exporter

template<typename MeshType,int N>
void expo_w ( boost::shared_ptr<MeshType> m,function_space_type::element_type u)
{
    auto x=Exporter<MeshType,N>::New();
    x->setMesh(m);
    x->add("u",u);
    x->save();
}

//creation of the libPyLapla library 

BOOST_PYTHON_MODULE(libPyLapla)
{

    if (import_mpi4py()<0) return ;


    
    def("unitSquare",unitSquare_w);
       
       
  //definition of form object and methods link to them 
  class_<form1_return_type>("Form1",no_init); 
  class_<form2_return_type>("Form2",no_init); 

    def("form1",form1_w);
    def("form2",form2_w);

    def("integrate",integrate_form1);
    def("integrate",integrate_form2);


  class_<function_space_type::element_type>("element",no_init);
    
     def("element",element_w);
    
    def("on",on_form2);
    def("solve",solve_w);
    def("exporter",expo_w<Mesh<Simplex<2>>,1>);
}

