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


    // definition of the Environment object and methods and classes link to it
    class_<Feel::detail::Environment,boost::noncopyable>("Environment", init<boost::python::list>()) 
        .def("worldComm",&Feel::detail::Environment::worldComm,return_value_policy<copy_non_const_reference>())
        .staticmethod("worldComm");

    // definition of the geometrical object (Simplex and Hypercube) and of the Mesh class 
    class_<Feel::Simplex<2>>("Simplex",init<>())
        .def("dim",&Feel::Simplex<2>::dimension);

    class_<Feel::Mesh<Feel::Simplex<2>>,boost::shared_ptr<Feel::Mesh<Feel::Simplex<2>>>,boost::noncopyable>("Mesh",init<>())
        .def("new",&Feel::Mesh<Simplex<2>>::New)
        .staticmethod("new")
        .def("clear",&Feel::Mesh<Simplex<2>>::clear);

    def("unitSquare",unitSquare_w);
       
      //definition of Pch class with others classes link to her and his constructor 
     class_<Feel::meta::Pch<Mesh<Simplex<2>>,1>,boost::shared_ptr<Feel::meta::Pch<Mesh<Simplex<2>>,1>>>("Pch",no_init);

    class_<Feel::FunctionSpaceBase>("FunctionSpaceBase",no_init);

    class_<Feel::FunctionSpace<Mesh<Simplex<2>>,Feel::bases<Feel::Lagrange<1,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>>,boost::shared_ptr<Feel::FunctionSpace<Mesh<Simplex<2>>,Feel::bases<Feel::Lagrange<1,Feel::Scalar,Feel::Continuous,Feel::PointSetEquiSpaced,0>>,double,Feel::Periodicity<Feel::NoPeriodicity>,Feel::mortars<Feel::NoMortar>>>,boost::python::bases<Feel::FunctionSpaceBase>>("FunctSpace",no_init);

def("newPch",Feel::Pch<1,PointSetEquiSpaced,Mesh<Simplex<2>>,0>); 
  
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

