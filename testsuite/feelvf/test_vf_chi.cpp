/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#define BOOST_TEST_MODULE test_vf_chi

#include <testsuite/testsuite.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( test_vf_chi )

BOOST_AUTO_TEST_CASE( test_vf_chi_1 )
{
    using namespace Feel;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();
    auto myVec = vec( cst(1.),cst(1.) );

    double eval = integrate(_range=elements(mesh),
                            _expr=chi( (myVec(0,0)) > cst(0.) ) ).evaluate()(0,0);

    form1(_test=Vh) =
        integrate(_range=elements(mesh),
                  _expr=chi( (myVec(0,0)) > cst(0.) )*id(u) );
    form2(_test=Vh,_trial=Vh) =
        integrate(_range=elements(mesh),
                  _expr=chi( (myVec(0,0)) > cst(0.) )*id(u)*idt(u) );

}
BOOST_AUTO_TEST_SUITE_END()
