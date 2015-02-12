/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Author(s): Cecile Daversin  <cecile.daversin@lncmi.cnrs.fr>
Date: 2011-12-07

Copyright (C) 2011 UJF
Copyright (C) 2011 CNRS

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3.0 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
  \file test_hcurl_lag.cpp
  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
  \date 2011-12-07
  */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE H_curl approximation
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/pchv.hpp>

#include <feel/feelfilters/loadmesh.hpp>

namespace Feel
{

  /**
   * This routine returns the list of options using the
   * boost::program_options library. The data returned is typically used
   * as an argument of a Feel::Application subclass.
   *
   * \return the list of options
   */
  inline
    po::options_description
    makeOptions()
    {
      po::options_description testhcurloptions( "test h_curl options" );
      testhcurloptions.add_options()
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ;
      return testhcurloptions.add( Feel::feel_options() );
    }

  inline
    AboutData
    makeAbout()
    {
      AboutData about( "test_hcurl_lag" ,
          "test_hcurl_lag" ,
          "0.1",
          "Test for h_curl space (Dim=2 Order=1)",
          AboutData::License_GPL,
          "Copyright (c) 2009 Universite Joseph Fourier" );
      about.addAuthor( "Vincent HUBER", "research engineer", "vincent.huber@cemosis.fr", "" );
      return about;

    }

  using namespace Feel;

  class TestHCurlLag
    :
      public Application
  {
    typedef Application super;

    public:

    //! numerical type is double
    typedef double value_type;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<HL_DIM,HL_ORDER> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    //typedef bases<Nedelec<0> > basis_type;
    typedef Nedelec<0,NedelecKind::NED1 > basis_type;
    typedef Lagrange<1,Scalar> lagrange_basis_s_type; //P1 scalar space
    typedef Lagrange<1,Vectorial> lagrange_basis_v_type; //P1 vectorial space

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, bases<basis_type> > space_type;
    typedef FunctionSpace<mesh_type, bases<lagrange_basis_s_type> > lagrange_space_s_type;
    typedef FunctionSpace<mesh_type, bases<lagrange_basis_v_type> > lagrange_space_v_type;
    typedef FunctionSpace<mesh_type, bases<lagrange_basis_v_type, lagrange_basis_s_type> > lagrange_space_v_s_type;
    typedef FunctionSpace<mesh_type, bases<basis_type, lagrange_basis_s_type> > mixed_space_s_type;
    typedef FunctionSpace<mesh_type, bases<basis_type, lagrange_basis_v_type> > mixed_space_v_type;

    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<lagrange_space_s_type> lagrange_space_s_ptrtype;
    typedef boost::shared_ptr<lagrange_space_v_type> lagrange_space_v_ptrtype;
    typedef boost::shared_ptr<lagrange_space_v_s_type> lagrange_space_v_s_ptrtype;
    typedef boost::shared_ptr<mixed_space_s_type> mixed_space_s_ptrtype;
    typedef boost::shared_ptr<mixed_space_v_type> mixed_space_v_ptrtype;

    /**
     * Constructor
     */
    TestHCurlLag()
      :
        super()
    {
      std::cout << "[TestHCurlLag]\n";
    }

    /**
     * run the application
     */
    void exampleProblem1();

    private:
  }; //TestHCurlLag

  // Resolve problem curl(curl(u)) + u = f with cross_prod(u,n) = 0 on boundary
  void
    TestHCurlLag::exampleProblem1()
    {
      auto mesh = loadMesh(_mesh = new mesh_type );

      space_ptrtype               Ah = space_type::New(mesh)           ;
      lagrange_space_s_ptrtype    Bh = lagrange_space_s_type::New(mesh);
      lagrange_space_v_ptrtype    Ch = lagrange_space_v_type::New(mesh);
      lagrange_space_v_s_ptrtype  Dh = lagrange_space_v_s_type::New(mesh);
      mixed_space_s_ptrtype       Eh = mixed_space_s_type::New(mesh)   ;
      mixed_space_v_ptrtype       Fh = mixed_space_v_type::New(mesh)   ;

      // Curl
      auto A1 = Ah->element();
      auto A2 = Ah->element();

      // Lagrange scalaire
      auto B1 = Bh->element();
      auto B2 = Bh->element();

      // Lagrange vectoriel
      auto C1 = Ch->element();
      auto C2 = Ch->element();

      // Lagrange vectoriel - scalaire
      auto D1  = Dh->element();
      auto D2  = Dh->element();
      auto d11 = D1.element<0>();
      auto d12 = D1.element<1>();
      auto d21 = D2.element<0>();
      auto d22 = D2.element<1>();

      // curl x lagrange scalaire
      auto E1  = Eh->element();
      auto E2  = Eh->element();
      auto e11 = E1.element<0>(); //curl
      auto e12 = E1.element<1>(); //scal
      auto e21 = E2.element<0>();
      auto e22 = E2.element<1>();

      // curl x lagrange vectoriel
      auto F1  = Fh->element();
      auto F2  = Fh->element();
      auto f11 = F1.element<0>();
      auto f12 = F1.element<1>();
      auto f21 = F2.element<0>();
      auto f22 = F2.element<1>();
      A1 = vf::project(_space=Ah, _range=elements(mesh), _expr=idv(A2));   // Curl
      B1 = vf::project(_space=Bh, _range=elements(mesh), _expr=idv(B2));   // Lag s
      C1 = vf::project(_space=Ch, _range=elements(mesh), _expr=idv(C2));   // Lag v
      d11 = vf::project(_space=Dh->functionSpace<0>(), _range=elements(mesh), _expr=idv(d21));
      d12 = vf::project(_space=Dh->functionSpace<1>(), _range=elements(mesh), _expr=idv(d22));
      e11 = vf::project(_space=Eh->functionSpace<0>(), _range=elements(mesh), _expr=idv(e21));
      e12 = vf::project(_space=Eh->functionSpace<1>(), _range=elements(mesh), _expr=idv(e22));
      f11 = vf::project(_space=Fh->functionSpace<0>(), _range=elements(mesh), _expr=idv(f21));
      f12 = vf::project(_space=Fh->functionSpace<1>(), _range=elements(mesh), _expr=idv(f22));

      auto a = form2(_test=Eh,_trial=Eh);
      a = integrate(_range=elements(mesh), _expr=inner(id(e21),trans(gradt(e12))));

    }
}

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( space )

BOOST_AUTO_TEST_CASE( test_hcurl_lag_example_1 )
{
  BOOST_TEST_MESSAGE( "test_hcurl_lag on example 1" );
  Feel::TestHCurlLag t;
  t.exampleProblem1();
  BOOST_TEST_MESSAGE( "test_hcurl_lag_N0 on example 1 done" );
}

BOOST_AUTO_TEST_SUITE_END()
#else

  int
main( int argc, char* argv[] )
{
  Feel::Environment env( argc,argv,
      makeAbout(), makeOptions() );

  Feel::TestHCurlLag app_hcurl_lag;

  app_hcurl_lag.exampleProblem1();
}

#endif
