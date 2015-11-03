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
#define BOOST_TEST_MODULE Regul
#include <testsuite/testsuite.hpp>
#endif

#include <feel/feel.hpp>
#include <feel/feelpde/preconditionerblockms.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feeldiscr/ned1h.hpp>

#if FEELPP_DIM == 2
#define curl_op curlx
#define curlt_op curlxt
#define curlv_op curlxv
#else
#define curl_op curl
#define curlt_op curlt
#define curlv_op curlv
#endif

using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description opts( "test_Regul" );
    opts.add_options()
    ( "penaldir", po::value<double>()->default_value( 0 ), "Use penaldir > 0 for weak BC" )
    ( "saveTimers", po::value<bool>()->default_value( true ), "print timers" )
    ;
    return opts.add( Feel::feel_options() )
        .add(Feel::backend_options("ms"));
}

inline
AboutData
makeAbout()
{
#if FEELPP_DIM==2
    AboutData about( "Regul2D" ,
                     "Regul2D" ,
                     "0.1",
                     "test Regul2D",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );
#else
    AboutData about( "Regul3D" ,
                     "Regul3D" ,
                     "0.1",
                     "test Regul3D",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );
#endif

    about.addAuthor( "Vincent HUBER", "developer", "vincent.huber@cemosis.fr", "" );

    return about;
}

///     \tparam DIM         Topological dimension.
template<int DIM>
class TestRegul : public Application
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
    typedef Lagrange<0, Scalar, Discontinuous> lag_0_basis_type; 
    typedef FunctionSpace<mesh_type, bases<lag_0_basis_type>> lag_0_space_type;
    typedef boost::shared_ptr<lag_0_space_type> lag_0_space_ptrtype;
    typedef typename lag_0_space_type::element_type lag_0_element_type;

    //! Pchv space
    typedef Lagrange<1, Vectorial> lag_v_basis_type;
    typedef FunctionSpace<mesh_type, bases<lag_v_basis_type>> lag_v_space_type;
    typedef boost::shared_ptr<lag_v_space_type> lag_v_space_ptrtype;
    typedef typename lag_v_space_type::element_type lag_v_element_type;

    typedef FunctionSpace<mesh_type, bases<curl_basis_type,lag_basis_type>> comp_space_type;
    typedef boost::shared_ptr<comp_space_type> comp_space_ptrtype;
    typedef typename comp_space_type::element_type comp_element_type;

    //! The exporter factory
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    //! Backends factory
    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef backend_type::solve_return_type solve_ret_type;

    public:

    /// Init the geometry with a circle/sphere from radius and characteristic length
    ///     \param radius   Circle or sphere radius.
    ///     \param h        Mesh size.
    TestRegul( ) 
    {
        auto M_mesh = loadMesh(_mesh=new mesh_type);
        auto Xh = curl_space_type::New(M_mesh); 
       
        // Exact Solution 
        auto f_M_a = expr<DIM,1>(soption("functions.a"));
        // Rhs - projected on marker COIL
        auto rhs = expr<DIM,1>(soption("functions.j"));

        auto u = Xh->element();
        auto v = Xh->element();
        
        auto f2 = form2(_test=Xh,_trial=Xh);
        auto f1 = form1(_test=Xh);

        std::cout << "[a;b;c] = [ " << doption("parameters.a") << ";" <<  doption("parameters.b") << ";" <<  doption("parameters.c") << "]" <<  std::endl;

        f1 = integrate(_range=markedelements(M_mesh,"COIL"),
                       _expr = doption("parameters.c")*inner(rhs,id(v)));    // rhs
        f2 = integrate(_range=elements(M_mesh),
                       _expr = 
                         doption("parameters.a")*(trans(curlt_op(u))*curl_op(v)) // (curl, curl)
                       + doption("parameters.b")*inner(idt(u),id(v))             // regul
                       );

        if(doption("penaldir")>0.)
        {
            std::cout << "Using weak BC\n";
#if DIM==2
          f1 += integrate(_range=boundaryfaces(M_mesh), _expr=
              - doption("parameters.a")*trans(curl_op(v))*cross(N(),vec(cst(0),cst(0))) 
              + doption("parameters.a")*doption("penaldir")/(hFace())*inner(cross(vec(cst(0),cst(0)),N()),cross(id(v),N())) );
#else
          f1 += integrate(_range=boundaryfaces(M_mesh), _expr=
              - doption("parameters.a")*trans(curl_op(v))*cross(N(),vec(cst(0),cst(0),cst(0))) 
              + doption("parameters.a")*doption("penaldir")/(hFace())*inner(cross(vec(cst(0),cst(0),cst(0)),N()),cross(id(v),N())) );
#endif
          f2 += integrate(_range=boundaryfaces(M_mesh), 
              _expr=- doption("parameters.a")*trans(curlt_op(u))*(cross(N(),id(v)) )
                    - doption("parameters.a")*trans(curl_op(v))*(cross(N(),idt(u)) )
                    + doption("parameters.a")*doption("penaldir")/(hFace())*inner(cross(idt(u),N()),cross(id(v),N())) );
        }
        else
        {
            std::cout << "Using strong BC\n";
        f2 += on(_range=boundaryfaces(M_mesh),
                 _rhs=f1,
                 _element=u,
#if DIM==2
                 _expr=vec(cst(0),cst(0))
#else
                 _expr=vec(cst(0),cst(0),cst(0))
#endif
                );
        }
        tic();
        f2.solveb(_rhs=f1,
                  _solution=u,
                  _backend=backend(_name="ms"));
        toc("Inverse",FLAGS_v>0);
        Environment::saveTimers(boption("saveTimers")); 
        auto e21 = normL2(_range=elements(M_mesh), _expr=(f_M_a-idv(u)));
        auto e22 = normL2(_range=elements(M_mesh), _expr=f_M_a);
        
        // export
        if(boption("exporter.export")){
            auto ex = exporter(_mesh=M_mesh);
            ex->add("potential"          ,u  );
            ex->save();
        }
    }

private:
    mesh_ptrtype M_mesh;
};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( Regul )

BOOST_AUTO_TEST_CASE( test )
{
    TestRegul<FEELPP_DIM> test;
}

//// Test 3D
//BOOST_AUTO_TEST_CASE( test_3d )
//{
//    TestRegul<3> test;
//}

BOOST_AUTO_TEST_SUITE_END()

#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _about=makeAbout(),
                           _desc=makeOptions() );
    
    TestRegul<FEELPP_DIM> t_regul;

    return 0;
}
#endif
