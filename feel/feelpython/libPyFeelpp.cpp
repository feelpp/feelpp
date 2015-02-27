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
   \file libPyFeelpp.cpp
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
#include<feel/feelvf/operators.hpp>


// definiton of Macros that will be used in the library creation with BOOST_PP_REPEAT

#define SIMPLEX(_,n,type) type<n+1,1>("Simplex");
#define HYPERCUBE(_,n,type) type<n+1,1>("Hypercube");

#define PCH1(_,n,type) type<1,n+1>();
#define PCH2(_,n,type) type<2,n+1>();
#define PCH3(_,n,type) type<3,n+1>();

// try to implemente all the Feelpp operators 
/*
#define VF_CLASS_DEF(_,OT)\
    VF_CLASS_DEF2 OT;

#define VF_CLASS_DEF2(O,T) \
     class_<Expr<BOOST_PP_CAT(expr_t,BOOST_PP_CAT(VF_OPERATOR_SYMBOL(O),VF_OP_TYPE_SUFFIX(T)))>>(BOOST_PP_STRINGIZE(BOOST_PP_CAT(expr_t,BOOST_PP_CAT(VF_OPERATOR_SYMBOL(O),VF_OP_TYPE_SUFFIX(T)))),no_init)

#define VF_METHODS_DEF(_,OT)\
    VF_METHODS_DEF2 OT;
    
#define VF_METHODS_DEF2(O,T)\
    def(BOOST_PP_STRINGIZE(BOOST_PP_CAT(VF_OPERATOR_SYMBOL(O),VF_OP_TYPE_SUFFIX(T))),Feel::vf::BOOST_PP_CAT(VF_OPERATOR_SYMBOL(O),VF_OP_TYPE_SUFFIX(T)))       
*/

using namespace boost::python;
using namespace Feel;

namespace py = boost::parameter::python;

//definition of all methods we need for the wrapping and wich are create from a BOOST_PARAMETER_FUNCTION or a method with default value arguments
 
    template<typename MeshType, int N>
void expo_w ( boost::shared_ptr<MeshType> m)
{
    auto x=Exporter<MeshType,N>::New();
    x->setMesh(m);
    x->addRegions();
    x->save();
}

template<typename MeshType>
boost::shared_ptr<MeshType> loadMesh_w (MeshType* mesh)
{
    return loadMesh(_mesh=mesh);
}

template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced,
         typename MeshType,
         int Tag = 0>

typename meta::Pch<MeshType,Order,T,Pts,Tag>::ptrtype
Pch_w( boost::shared_ptr<MeshType> mesh)
{
    return Pch<Order,T,Pts,MeshType,Tag>(mesh);
}



// method that will define all the elements link to Mesh object into the Python library

template <int n,int N>
void def_wrapper (std::string s)
{    
    std::ostringstream f;
    std::ostringstream g;

    if(s.compare("Simplex")==0)
    {
        f<<"Simplex"<<n;
        g<<"MeshS"<<n;
        class_<Feel::Simplex<n>>(f.str().c_str(),init<>());
        class_<Feel::Mesh<Simplex<n>>,boost::shared_ptr<Feel::Mesh<Simplex<n>>>,boost::noncopyable>(g.str().c_str(),init<>())
            .def("new",&Feel::Mesh<Simplex<n>>::New)
            .staticmethod("new");

        def("loadMesh",loadMesh_w<Mesh<Simplex<n>>>);
        def("export",expo_w<Mesh<Simplex<n>>,N>);

    }

    else if(s.compare("Hypercube")==0)
    {
        f<<"Hypercube"<<n;
        g<<"MeshH"<<n;   
        class_<Feel::Hypercube<n>>(f.str().c_str(),init<>());
        class_<Feel::Mesh<Hypercube<n>>,boost::shared_ptr<Feel::Mesh<Hypercube<n>>>,boost::noncopyable>(g.str().c_str(),init<>())
            .def("new",&Feel::Mesh<Hypercube<n>>::New)
            .staticmethod("new");

        def("loadMesh",loadMesh_w<Mesh<Hypercube<n>>>);
        def("export",expo_w<Mesh<Hypercube<n>>,N>);
    }

}

// method that will define all the elements link to Pch object into the Python library

    template <int n,int k>
void def_wrapper_Pch ()
{
    std::ostringstream h;
    std::ostringstream i;
    std::ostringstream j;

    j<<"newPch"<<k;

    h<<"PchS"<<n<<k;
    i<<"FunctSpaceS"<<n<<k;
    
    class_<Feel::meta::Pch<Mesh<Simplex<n>>,k>,boost::shared_ptr<Feel::meta::Pch<Mesh<Simplex<n>>,k>>>(h.str().c_str(),no_init);

    class_<Feel::FunctionSpace<Mesh<Simplex<n>>,Feel::bases<Feel::Lagrange<k,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>>,boost::shared_ptr<Feel::FunctionSpace<Mesh<Simplex<n>>,Feel::bases<Feel::Lagrange<k,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>>>,boost::python::bases<Feel::FunctionSpaceBase>>(i.str().c_str(),no_init);

    def(j.str().c_str(),Pch_w<k,double,PointSetEquiSpaced,Mesh<Simplex<n>>,0>);

    h.str("");
    i.str(""); 

    h<<"PchH"<<n<<k;
    i<<"FunctSpaceH"<<n<<k;

    class_<Feel::meta::Pch<Mesh<Hypercube<n>>,k>,boost::shared_ptr<Feel::meta::Pch<Mesh<Hypercube<n>>,k>>>(h.str().c_str(),no_init);

    class_<Feel::FunctionSpace<Mesh<Hypercube<n>>,Feel::bases<Feel::Lagrange<k,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>>,boost::shared_ptr<Feel::FunctionSpace<Mesh<Hypercube<n>>,Feel::bases<Feel::Lagrange<k,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>>>,boost::python::bases<Feel::FunctionSpaceBase>>(i.str().c_str(),no_init);

    def(j.str().c_str(),Pch_w<k,double,PointSetEquiSpaced,Mesh<Hypercube<n>>,0>);

    h.str("");
    i.str(""); 
}

/*
#define VF_DEF(_,OT) \
    typedef VF_OPERATOR_NAME(O)<ELEM,VF_OP_TYPEOBJECT(T)> BOOST_PP_CAT(expr_t,BOOST_PP_CAT(VF_OPERATOR_SYMBOL(0) , VF_OP_TYPE_SUFFIX(T)));

    BOOST_PP_LIST_FOR_EACH_PRODUCT(VF_DEF,2,(VF_OPERATORS,VF_OPERATORS_TYPE))

*/

//creation of the libPyFeelpp library 

BOOST_PYTHON_MODULE(libPyFeelpp)
{

    if (import_mpi4py()<0) return ;

    class_<Feel::Environment,boost::noncopyable>("Environment", init<boost::python::list>()) 
        .def("worldComm",&Feel::Environment::worldComm,return_value_policy<copy_non_const_reference>())
        .staticmethod("worldComm");

    class_<WorldComm>("WorldComm",init<>());

    class_<Feel::FunctionSpaceBase>("FunctionSpaceBase",no_init);

    /*
    def_wrapper<1,1>("Simplex");
    def_wrapper<2,1>("Simplex");
    def_wrapper<3,1>("Simplex");
    def_wrapper<1,1>("Hypercube");
    def_wrapper<2,1>("Hypercube");
    def_wrapper<3,1>("Hypercube");
    
    
    def_wrapper_Pch<1,1>();
    def_wrapper_Pch<1,2>();
    def_wrapper_Pch<1,3>();

    def_wrapper_Pch<2,1>();
    def_wrapper_Pch<2,2>();
    def_wrapper_Pch<2,3>();

    def_wrapper_Pch<3,1>();
    def_wrapper_Pch<3,2>();
    def_wrapper_Pch<3,3>();
    */


    // use of BOOST_PP_REPEAT with macro and method define before to define objects we need for all dimension

    BOOST_PP_REPEAT(3,SIMPLEX,def_wrapper)
    BOOST_PP_REPEAT(3,HYPERCUBE,def_wrapper)


    BOOST_PP_REPEAT(3,PCH1,def_wrapper_Pch)
    BOOST_PP_REPEAT(3,PCH2,def_wrapper_Pch)
    BOOST_PP_REPEAT(3,PCH3,def_wrapper_Pch)

    /*
    BOOST_PP_LIST_FOR_EACH_PRODUCT(VF_CLASS_DEF,2,(VF_OPERATORS,VF_OPERATORS_TYPE))

    BOOST_PP_LIST_FOR_EACH_PRODUCT(VF_METHODS_DEF,2,(VF_OPERATORS,VF_OPERATORS_TYPE))
    */
}

