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
   \file libPyInteg.cpp
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

// retrieval of complicated return types with the "give bad type to have the good one" trick  

typedef boost::tuples::tuple<mpl_::size_t<0ul>, boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double>, std::allocator<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double> > > > > > > > > > >, boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double>, std::allocator<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double> > > > > > > > > > >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> elements_return_type; 

typedef boost::tuples::tuple<mpl_::size_t<0ul>, boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double>, std::allocator<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double> > > > > > >, boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double>, std::allocator<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double> > > > > > >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> boundaryelements_return_type;

typedef boost::tuples::tuple<mpl_::size_t<1ul>, boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement2D<(unsigned short)3, Feel::Simplex<(unsigned short)2, (unsigned short)1, (unsigned short)3>, Feel::SubFaceOf<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double> >, double>, std::allocator<Feel::GeoElement2D<(unsigned short)3, Feel::Simplex<(unsigned short)2, (unsigned short)1, (unsigned short)3>, Feel::SubFaceOf<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double> >, double> > > > >, boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement2D<(unsigned short)3, Feel::Simplex<(unsigned short)2, (unsigned short)1, (unsigned short)3>, Feel::SubFaceOf<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double> >, double>, std::allocator<Feel::GeoElement2D<(unsigned short)3, Feel::Simplex<(unsigned short)2, (unsigned short)1, (unsigned short)3>, Feel::SubFaceOf<Feel::GeoElement3D<(unsigned short)3, Feel::Simplex<(unsigned short)3, (unsigned short)1, (unsigned short)3>, double> >, double> > > > >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type> boundaryfaces_return_type;

typedef Feel::vf::Expr<Feel::vf::Integrator<boost::tuples::tuple<mpl_::size_t<0>,boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement3D<3, Feel::Simplex<3, 1, 3>, double>, std::allocator<Feel::GeoElement3D<3, Feel::Simplex<3, 1, 3>, double> > > > > > > > > > >,boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement3D<3,Feel::Simplex<3, 1, 3>, double>, std::allocator<Feel::GeoElement3D<3, Feel::Simplex<3, 1, 3>, double> > > > > > > > > > >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type>,Feel::_Q<2>, Feel::vf::Expr<Feel::vf::GinacEx<2> >, Feel::_Q<2> > >  integrate_return_type;


typedef Feel::vf::Expr<Feel::vf::Integrator<boost::tuples::tuple<mpl_::size_t<0>,boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement3D<3,Feel::Simplex<3, 1, 3>, double>, std::allocator<Feel::GeoElement3D<3, Feel::Simplex<3, 1, 3>, double> > > > > > >,boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement3D<3,Feel::Simplex<3, 1, 3>, double>, std::allocator<Feel::GeoElement3D<3, Feel::Simplex<3, 1, 3>, double> > > > > > >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type,boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type>,Feel::_Q<2>, Feel::vf::Expr<Feel::vf::GinacEx<2> >, Feel::_Q<2> > > integratebound_return_type;


typedef Feel::vf::Expr<Feel::vf::Integrator<boost::tuples::tuple<mpl_::size_t<0>,boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement3D<3,Feel::Simplex<3, 1, 3>, double>, std::allocator<Feel::GeoElement3D<3, Feel::Simplex<3, 1, 3>, double> > > > > > > > > > >,boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement3D<3,Feel::Simplex<3, 1, 3>, double>, std::allocator<Feel::GeoElement3D<3, Feel::Simplex<3, 1, 3>, double> > > > > > > > > > >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type,boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type>,Feel::_Q<2>, Feel::vf::Expr<Feel::vf::GinacMatrix<1, 2, 2> >, Feel::_Q<2> > > integrategrad_return_type;


typedef Feel::vf::Expr<Feel::vf::Integrator<boost::tuples::tuple<mpl_::size_t<1>,boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement2D<3,Feel::Simplex<2, 1, 3>, Feel::SubFaceOf<Feel::GeoElement3D<3, Feel::Simplex<3, 1, 3>, double> >, double>,std::allocator<Feel::GeoElement2D<3, Feel::Simplex<2, 1, 3>, Feel::SubFaceOf<Feel::GeoElement3D<3,Feel::Simplex<3, 1, 3>, double> >, double> > > > >,boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<boost::multi_index::detail::index_node_base<Feel::GeoElement2D<3,Feel::Simplex<2, 1, 3>, Feel::SubFaceOf<Feel::GeoElement3D<3, Feel::Simplex<3, 1, 3>, double> >, double>,std::allocator<Feel::GeoElement2D<3, Feel::Simplex<2, 1, 3>, Feel::SubFaceOf<Feel::GeoElement3D<3,Feel::Simplex<3, 1, 3>, double> >, double> > > > >, boost::tuples::null_type, boost::tuples::null_type,boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type,boost::tuples::null_type>, Feel::_Q<2>, Feel::vf::Expr<Feel::vf::GinacEx<2> >, Feel::_Q<2> > > integrateboundfaces_return_type;


using namespace boost::python;
using namespace Feel;

namespace py = boost::parameter::python;


//definition of all methods we need for the wrapping and wich are create from a BOOST_PARAMETER_FUNCTION or a method with default value arguments  

boost::shared_ptr<Exporter<Mesh<Simplex<3>>,1>> New1 (po::variables_map const& x,std::string y,WorldComm const& z) 
{
    return Exporter<Mesh<Simplex<3>>,1>::New(x,y,z);
}

boost::shared_ptr<Exporter<Mesh<Simplex<3>>,1>> New2 () 
{
    return Exporter<Mesh<Simplex<3>>,1>::New();
}



    template<typename MeshType,int N>
void expo_w ( boost::shared_ptr<MeshType> m)
{
    auto x=Exporter<MeshType,N>::New();
    x->setMesh(m);
    x->addRegions();
    x->save();
}


boost::shared_ptr<Mesh<Simplex<3>>> loadMesh_w (Mesh<Simplex<3>>* mesh)
{
    return loadMesh(_mesh=mesh);
}

std::string soption_w (std::string name)
{
    return soption(_name=name);
}

Expr< GinacEx<2> > expr_w( std::string const& s )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    return Expr< GinacEx<2> >(  GinacEx<2>( g.first, g.second,"") );
}

template<int Order=2>
    Expr<GinacMatrix<1,1,Order> >
laplacian_w( Expr<GinacEx<Order>> const& s)
{
    return expr<1,1,Order>( GiNaC::laplacian(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), "");
}

template<int M,int Order=2>
    Expr<GinacMatrix<1,M,Order> >
grad_w( Expr<GinacEx<Order>> const& s)
{
    return expr<1,M,Order>( GiNaC::grad(s.expression().expression(),s.expression().symbols()), s.expression().symbols(),"" );
}


// elements on the mesh 

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
    typename MeshTraits<MeshType>::element_const_iterator,
    typename MeshTraits<MeshType>::element_const_iterator>
elements_w( MeshType const& mesh )
{
    return elements( mesh);
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
    typename MeshTraits<MeshType>::location_element_const_iterator,
    typename MeshTraits<MeshType>::location_element_const_iterator>
boundaryelements_w( MeshType const& mesh)
{
    return boundaryelements(mesh);
}



template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
    typename MeshTraits<MeshType>::location_face_const_iterator,
    typename MeshTraits<MeshType>::location_face_const_iterator>
boundaryfaces_w( MeshType const& mesh  )
{
    return boundaryfaces(mesh);
}


// integrate

integrate_return_type integrate_w ( elements_return_type e , Expr<GinacEx<2>> g)
{
    return integrate(_range=e,_expr=g);
}

integratebound_return_type integrate_w2 ( boundaryelements_return_type e , Expr<GinacEx<2>> g)
{
    return integrate(_range=e,_expr=g);
}

integrategrad_return_type integrate_w3 ( elements_return_type e,Expr<GinacMatrix<1,2,2>> grad)
{
    return integrate(_range=e,_expr=grad);
}

integrateboundfaces_return_type integrate_w4 ( boundaryfaces_return_type e , Expr<GinacEx<2>> g)
{
    return integrate(_range=e,_expr=g);
}

//evaluate

integrate_return_type::value_type evaluate_w (integrate_return_type i)
{
    return i.evaluate();
}

integratebound_return_type::value_type evaluate_w2 (integratebound_return_type i)
{
    return i.evaluate();
}

integrategrad_return_type::value_type evaluate_w3 (integrategrad_return_type i)
{
    return i.evaluate();
}

integrateboundfaces_return_type::value_type evaluate_w4 (integrateboundfaces_return_type i)
{
    return i.evaluate();
}

// definition of print methods for some objects

void printMa1 (integrate_return_type::value_type m)
{
    std::cout<< m << std::endl;
}

void printMa2 (integratebound_return_type::value_type m)
{
    std::cout<< m << std::endl;
}

void printMa3 (integrategrad_return_type::value_type m)
{
    std::cout<< m << std::endl;
}

void printMa4 (integrateboundfaces_return_type::value_type m)
{
    std::cout<< m << std::endl;
}

void printExpr1(Expr<GinacEx<2>> e)
{
    std::cout<< e << std::endl;
}

void printExpr2(Expr<GinacMatrix<1,2,2>> e)
{
    std::cout<< e << std::endl;
}

//creation of the libPyInteg library  

BOOST_PYTHON_MODULE(libPyInteg)
{

    if (import_mpi4py()<0) return ;


// definition of the Expr object 

    class_<Expr<GinacEx<2>>>("Expr",no_init);
    def("soption",soption_w);
    def("expr",expr_w);

    class_<Expr<GinacMatrix<1,1,2>>>("ExprLapla",no_init);
    class_<Expr<GinacMatrix<1,2,2>>>("ExprLapla",no_init);


    def("laplacian",laplacian_w<2>);

// definition of methods link to elements on the mesh 

    class_<elements_return_type>("Elements_return_type",no_init);
    def("elements",elements_w<Mesh<Simplex<3>>>);

    class_<boundaryelements_return_type>("Boundaryelements_return_type",no_init);
    def("boundaryelements",boundaryelements_w<Mesh<Simplex<3>>>);

    class_<boundaryfaces_return_type>("Boundaryfaces_return_type",no_init);
    def("boundaryfaces",boundaryfaces_w<Mesh<Simplex<3>>>);

//definition of all the integrate method we need for this example
    
    def("grad",grad_w<2,2>);

    class_<integrate_return_type>("Integrate_return_type",no_init);
    class_<integratebound_return_type>("Integratebound_return_type",no_init);
    class_<integrategrad_return_type>("Integrategrad_return_type",no_init);
    class_<integrateboundfaces_return_type>("Integrateboundfaces_return_type",no_init);

    def("integrate",integrate_w);
    def("integrate",integrate_w2);
    def("integrate",integrate_w3);
    def("integrate",integrate_w4);

// definition of all the evaluate method we need for this example

    class_<Eigen::Matrix<double,1,1,0,1,1>>("Matrix1",no_init);
    class_<Eigen::Matrix<double,1,3,1,1,3>>("Matrix1",no_init);

    def("evaluate",evaluate_w);
    def("evaluate",evaluate_w2);
    def("evaluate",evaluate_w3);
    def("evaluate",evaluate_w4);

//definition of the print methods 

    def("printSol",printMa1);
    def("printSol",printMa2);
    def("printSol",printMa3);
    def("printSol",printMa4);

    def("printExpr",printExpr1);
    def("printExpr",printExpr2);

}

