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
#include <feel/feelalg/matrixeigensparse.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/dh.hpp>

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
        ( "solveIt", po::value<bool>()->default_value( true ), "Solve the pb ?" )
        ( "model", po::value<std::string>()->default_value( "torus_quart.mod" ), "Name of the model to use" )
        ( "penaldir", po::value<double>()->default_value( 0 ), "Use penaldir > 0 for weak BC" )
        ("vec_curl" , po::value<std::string>()->default_value( "{1,1,1}:x:y:z" ), "function st curl(vec_curl) is exported" )
        ("scal_grad", po::value<std::string>()->default_value( "1:x:y:z" ), "function st grad(scal_grad) is exported" )
        ("vec_div"  , po::value<std::string>()->default_value( "{1,1,1}:x:y:z" ), "function st div(vec_div) is exported" )
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

public:

    /// Init the geometry with a circle/sphere from radius and characteristic length
    ///     \param radius   Circle or sphere radius.
    ///     \param h        Mesh size.
    TestRegul( ) 
    {
        LOG(INFO) << "DIM = " << DIM << std::endl;
        auto M_mesh = loadMesh(_mesh=new Mesh<Simplex<DIM>>);
        /*
         * P1 - grad -> Hcurl - curl -> Hdiv - div -> L2
         */
        auto curl_h = Ned1h<0>(M_mesh);    // Hcurl
        auto div_h  = Dh<0>(M_mesh);       // HDiv - RT
        auto h1_h   = Pch<1>(M_mesh);      // Lagrange
        auto l2_h   = Pdh<0>(M_mesh);      // Lagrange P_0 disc

        auto Idiv  = Div( _domainSpace=div_h, _imageSpace=l2_h);
        auto Icurl = Curl( _domainSpace=curl_h, _imageSpace=div_h);
        auto Igrad = Grad( _domainSpace=h1_h, _imageSpace=curl_h);

        auto ozz = curl_h->element("ozz");
        auto zoz = curl_h->element("zoz");
        auto zzo = curl_h->element("zzo");
        ozz.on(_range=elements(M_mesh), _expr=vec(cst(1.),cst(0.),cst(0.)));
        zoz.on(_range=elements(M_mesh), _expr=vec(cst(0.),cst(1.),cst(0.)));
        zzo.on(_range=elements(M_mesh), _expr=vec(cst(0.),cst(0.),cst(1.)));

        auto ozzl = h1_h->element("ozz");
        ozzl.on(_range=elements(M_mesh), _expr=cst(1.));

        std::cout << "Idiv : " << Idiv.matPtr()->size1() << " : " << Idiv.matPtr()->size2() << std::endl;
        std::cout << "Icurl : " << Icurl.matPtr()->size1() << " : " << Icurl.matPtr()->size2() << std::endl;
        std::cout << "Igrad : " << Igrad.matPtr()->size1() << " : " << Igrad.matPtr()->size2() << std::endl;

        Idiv.matPtr()->printMatlab("idiv.m");
        Icurl.matPtr()->printMatlab("icurl.m");
        Igrad.matPtr()->printMatlab("igrad.m");
        ozz.printMatlab("cstNed1.m");
        zoz.printMatlab("cstNed2.m");
        zzo.printMatlab("cstNed3.m");
        ozzl.printMatlab("cstLag1.m");

        // Exporting element to graphically find out where on the mesh bug can
        // be
        auto curl_cst = div_h->element();
        auto curl_cst_exp = curl_h->element();
        curl_cst_exp.on(_range=elements(M_mesh), _expr=expr<DIM,1>(soption("vec_curl")));
        Icurl.matPtr()->multVector(curl_cst_exp, curl_cst);

        auto grad_cst = curl_h->element();
        auto grad_cst_exp = h1_h->element();
        grad_cst_exp.on(_range=elements(M_mesh), _expr=expr(soption("scal_grad")));
        Igrad.matPtr()->multVector(grad_cst_exp, grad_cst);

        auto div_cst = l2_h->element();
        auto div_cst_exp = div_h->element();
        div_cst_exp.on(_range=elements(M_mesh), _expr=expr<DIM,1>(soption("vec_div")));
        Idiv.matPtr()->multVector(div_cst_exp, div_cst);

        // Now, we solve the pb
        auto u = curl_h->element(); // Potenial - unknown
        auto e = curl_h->element(); // exact
        auto v = curl_h->element(); // test 
        auto w = div_h->element();  // curl(potentential)
        auto ec= div_h->element();  // curl(exact)

        auto f2 = form2(_test=curl_h,_trial=curl_h,_properties=MatrixProperties::SPD);
        auto f1 = form1(_test=curl_h);
        if(boption("solveIt"))
        { 

            ModelProperties model;

            double a = model.parameters()["a"].value();
            double b = model.parameters()["b"].value();
            double c = model.parameters()["c"].value();

            std::cout << "[a;b;c] = [ " << a << ";" <<  b << ";" <<  c << "]" <<  std::endl;

            for(auto it:model.materials())
            {

                std::string entry = "Materials.";
                entry += marker(it)+".ex";
                LOG(INFO) << "reading " << entry << "...\n";
                // Exact 
                e.on(_range=markedelements(M_mesh,marker(it)),
                     _expr=expr<DIM,1>(model.getEntry(entry)));
                // Rhs
                entry = "Materials.";
                entry += marker(it)+".j";
                LOG(INFO) << "reading " << entry << "...\n";
                f1 += integrate(_range=markedelements(M_mesh,marker(it)),
                                _expr = c*inner(expr<DIM,1>(model.getEntry(entry)),id(v)));    // rhs
            }
            // Lhs
            f2 = integrate(_range=elements(M_mesh),
                           _expr = 
                           a*(trans(curlt_op(u))*curl_op(v)) // (curl, curl)
                           + b*inner(idt(u),id(v))             // regul
                          );

            if(doption("penaldir")>0.)
            {
                std::cout << "Using weak BC\n";
                for(auto it: model.boundaryConditions().getVectorFields<DIM> ( "u", "Dirichlet" ) )
                {
                    f1 += integrate(_range=markedfaces(M_mesh,it.first), 
                                    _expr=
                                    - a*trans(curl_op(v))*cross(N(),it.second) 
                                    + a*doption("penaldir")/(hFace())*inner(cross(it.second,N()),cross(id(v),N())) );
                    f2 += integrate(_range=boundaryfaces(M_mesh), 
                                    _expr=
                                    - a*trans(curlt_op(u))*(cross(N(),id(v)) )
                                    - a*trans(curl_op(v))*(cross(N(),idt(u)) )
                                    + a*doption("penaldir")/(hFace())*inner(cross(idt(u),N()),cross(id(v),N())) );
                }
            }
            else
            {
                std::cout << "Using strong BC\n";
                for(auto it: model.boundaryConditions().getVectorFields<DIM> ( "u", "Dirichlet" ) )
                {
                    f2 += on(_range=markedfaces(M_mesh,it.first),
                             _rhs=f1,
                             _element=u,
                             _expr=it.second);
                }
            }
            CHECK( f2.matrixPtr()->isSymmetric(true) ) << "Matrix is not symmetric!";
            CHECK( f2.matrixPtr()->isPositiveDefinite() ) << "Matrix is not Definite Positive!";
            CHECK( f2.matrixPtr()->isSPD() ) << "Matrix is not Symmetric Definite Positive!";
            tic();
            f2.solveb(_rhs=f1,
                      _solution=u,
                      _backend=backend(_name="ms"));
            toc("Inverse",FLAGS_v>0);

#if DIM==2 
            w.on( _range=elements(M_mesh),  _expr=vec(curlv_op(u),cst(0.)) );
            ec.on( _range=elements(M_mesh), _expr=vec(curlv_op(e),cst(0.)) );
#else
            w.on( _range=elements(M_mesh), _expr=curlv_op(u) );
            ec.on( _range=elements(M_mesh), _expr=curlv_op(e) );
#endif
            auto e21 = normL2(_range=elements(M_mesh), _expr=(idv(e)-idv(u)));
            auto e22 = normL2(_range=elements(M_mesh), _expr=idv(e));
            auto ec21 = normL2(_range=elements(M_mesh), _expr=(curlv_op(e)-curlv_op(u)));
            auto ec22 = normL2(_range=elements(M_mesh), _expr=curlv_op(e));
            std::cout << "Erreur " << e21 << "\t" << e21/e22  << "\t" << ec21 << "\t" << ec21/ec22 << std::endl;
        } 
        // export
        if(boption("exporter.export"))
        {
            auto ex = exporter(_mesh=M_mesh);
            for(int i = 0; i < 2; i++)
            {
                v.zero();
                v = f1.vector();
                ex->step(i)->add("curl_cst_exp",curl_cst_exp);
                ex->step(i)->add("grad_cst_exp",grad_cst_exp);
                ex->step(i)->add("div_cst_exp",div_cst_exp);
                ex->step(i)->add("curl_cst",curl_cst);
                ex->step(i)->add("grad_cst",grad_cst);
                ex->step(i)->add("div_cst",div_cst);
                if(boption("solveIt"))
                {
                    ex->step(i)->add("exactP", e);
                    ex->step(i)->add("potential", u);
                    ex->step(i)->add("j", v);
                    ex->step(i)->add("B", w);
                    ex->step(i)->add("Be", ec);
                }
                ex->save();
            }
        }
    }
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
