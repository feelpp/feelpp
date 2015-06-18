/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
   -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
   \file dist2wallsoptimized.cpp
   \author Guillaume Dolle <gdolle at unistra.fr>
   \date 2014-01-21
 */

//#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE precAFP
#include <testsuite/testsuite.hpp>
#endif

#include <feel/feel.hpp>
#include <feel/feelpde/preconditionerblockms.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feeldiscr/ned1h.hpp>

using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description opts( "test_precAFP" );
    opts.add_options()
    ( "myModel", po::value<std::string>()->default_value( "model.mod" ), "name of the model" )
    ;
    return opts.add( Feel::feel_options() )
        .add(Feel::backend_options("ms"));
}

inline
AboutData
makeAbout()
{
    AboutData about( "precAFP" ,
                     "precAFP" ,
                     "0.1",
                     "test precAFP",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Vincent HUBER", "developer", "vincent.huber@cemosis.fr", "" );

    return about;
}

///     \tparam DIM         Topological dimension.
template<int DIM>
class TestPrecAFP : public Application
{
    private:
    typedef Application super;
    //! Numerical type is double
    typedef double value_type;

    //! Simplexes of order ORDER
    typedef Simplex<DIM> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! Hcurl space
    typedef Nedelec<0,NedelecKind::NED1 > curl_basis_type;
    typedef FunctionSpace<mesh_type, bases<curl_basis_type>> curl_space_type;
    typedef boost::shared_ptr<curl_space_type> curl_space_ptrtype;
    typedef typename curl_space_type::element_type curl_element_type;

    //! Pch space
    typedef Lagrange<1, Scalar> lag_basis_type; 
    typedef FunctionSpace<mesh_type, bases<lag_basis_type>> lag_space_type;
    typedef boost::shared_ptr<lag_space_type> lag_space_ptrtype;
    typedef typename lag_space_type::element_type lag_element_type;

    //! Pch 0 space
    typedef Lagrange<2, Scalar, Discontinuous> lag_0_basis_type; 
    typedef FunctionSpace<mesh_type, bases<lag_0_basis_type>, Continuous> lag_0_space_type;
    typedef boost::shared_ptr<lag_0_space_type> lag_0_space_ptrtype;
    typedef typename lag_0_space_type::element_type lag_0_element_type;

    //! Pchv space
    typedef Lagrange<1, Vectorial> lag_v_basis_type;
    typedef FunctionSpace<mesh_type, bases<lag_v_basis_type>> lag_v_space_type;
    typedef boost::shared_ptr<lag_v_space_type> lag_v_space_ptrtype;
    typedef typename lag_v_space_type::element_type lag_v_element_type;

#ifndef STAB
    typedef FunctionSpace<mesh_type, bases<curl_basis_type,lag_basis_type>> comp_space_type;
    typedef boost::shared_ptr<comp_space_type> comp_space_ptrtype;
    typedef typename comp_space_type::element_type comp_element_type;
#endif
    
    //! The exporter factory
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    //! Backends factory
    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    public:

    /// Init the geometry with a circle/sphere from radius and characteristic length
    ///     \param radius   Circle or sphere radius.
    ///     \param h        Mesh size.
    TestPrecAFP( ) 
    {
        auto M_mesh = loadMesh(_mesh=new mesh_type);
        auto Xh = comp_space_type::New(M_mesh); // curl x lag
        auto Mh = lag_0_space_type::New(M_mesh); // lag_0
        auto Jh = lag_v_space_type::New(M_mesh); // lag_v

        auto model = ModelProperties(Environment::expand(soption("myModel")));
        auto J = vf::project(_space=Jh,
                             _range=elements(M_mesh),
                             _expr=expr<DIM,1>(soption("functions.j")));
        auto M_mu = vf::project(_space=Mh,
                               _range=elements(M_mesh),
                               _expr=expr(soption("functions.m")));
        
        auto U = Xh->element();
        auto u = U.template element<0>();
        auto v = U.template element<0>();
        auto phi = U.template element<1>();
        auto psi = U.template element<1>();
        
        auto f2 = form2(_test=Xh,_trial=Xh);
        auto f1 = form1(_test=Xh);
       
       auto M_prec = blockms(
          _space = Xh,
          _space2 = Mh, 
          _matrix = f2.matrixPtr(),
          _bc = model.boundaryConditions());
        map_vector_field<DIM,1,2> m_dirichlet {model.boundaryConditions().getVectorFields<DIM>("u","Dirichlet")};
        map_scalar_field<2> m_dirichlet_phi {model.boundaryConditions().getScalarFields<DIM>("phi","Dirichlet")};
        
        f1 = integrate(_range=elements(M_mesh),
                       _expr = inner(idv(J),id(u)));    // rhs
        f2 = integrate(_range=elements(M_mesh),
                       _expr = 
                       - inner(trans(id(v)),gradt(phi)) // -grad(phi)
                       + inner(trans(idt(u)),grad(psi)) // div(u) = 0
                       + (1./idv(M_mu))*(trans(curlxt(u))*curlx(u)) // curl curl :: 2D VERSION
                       );
        for(auto const & it : m_dirichlet)
            f2 += on(_range=markedfaces(M_mesh,it.first),
                     _rhs=f1,
                     _element=u,
                     _expr = it.second);
        for(auto const & it : m_dirichlet_phi)
            f2 += on(_range=markedfaces(M_mesh,it.first),
                     _rhs=f1,
                     _element=phi,
                     _expr = it.second);
        if(soption("ms.pc-type") == "AFP" ){
        M_prec->update(f2.matrixPtr(),M_mu);
        f2.solveb(_rhs=f1,
                  _solution=U,
                  _backend=backend(_name="ms"),
                  _prec=M_prec);
        }else{
        f2.solveb(_rhs=f1,
                  _solution=U,
                  _backend=backend(_name="ms"));
        }
    }

private:
    mesh_ptrtype M_mesh;

};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( levelset )

// Test 2D
BOOST_AUTO_TEST_CASE( test_2d )
{
    TestPrecAFP<2> test;
}

//// Test 3D
//BOOST_AUTO_TEST_CASE( test_3d )
//{
//    TestPrecAFP<3> test;
//}

BOOST_AUTO_TEST_SUITE_END()

#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _about=makeAbout(),
                           _desc=makeOptions() );
    
    TestPrecAFP<2> tl2;
    //TestPrecAFP<3> tl3;

    return 0;
}
#endif
