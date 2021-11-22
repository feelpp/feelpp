/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
   \file test_ginac.cpp
   \author Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2013-04-04
*/

#define BOOST_TEST_MODULE test_ginac
#include <feel/feelcore/testsuite.hpp>

#include <iostream>
#include <string>
#include <list>

#include <boost/algorithm/string/split.hpp>
#include <boost/foreach.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feelvf/vf.hpp>



using namespace Feel;

class MyVisitor :
    public GiNaC::visitor,
    public GiNaC::symbol::visitor
{
public:
    const std::list<symbol> & symbols()
    {
        return syms;
    }

    const std::list<std::string> & symbol_names()
    {
        sym_names.sort();
        sym_names.unique();
        return sym_names;
    }

    bool hassymbol(std::string const& s, symbol &found_symbol)
    {
        bool result = false;
        boost::for_each(syms, [&s, &result, &found_symbol](symbol const& sym)
                        {
                            if (sym.get_name() == s)
                                {
                                    result = true;
                                    found_symbol = sym;
                                }
                        });
        return result;
    }

private:
    std::list<symbol> syms;
    std::list<std::string> sym_names;

    void visit(const symbol & s)
    {
        syms.push_back(s);
        sym_names.push_back(s.get_name());
    }

};
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description ginacoptions("Ginac options");
    ginacoptions.add_options()
        ("dim", Feel::po::value<int>()->default_value( 1 ), "geometric dimension")
        ("params", Feel::po::value<std::string>()->default_value( "a;b" ), "name of parameters")
        ("exact", Feel::po::value<std::string>()->default_value( "a*x+b" ), "name of the input")
        ;
    return ginacoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_ginac" ,
                           "test_ginac" ,
                           "0.1",
                           "test ginac integration with Feelpp",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2012-2013 Laboratoire national des Champs magnetiques Intenses");

    about.addAuthor("Christophe Trophime", "developer", "christophe.trophime@lncmi.cnrs.fr", "");
    return about;

}


void runTest0()
{
    using GiNaC::symtab;
    using GiNaC::symbol;

    ex exact_parsed;
    std::vector<symbol> vars, parameters;
    std::vector<std::string> lst_params;
    std::string exact;
    std::string params;
    int dim = 0;

    if ( boption(_name="ginac.strict-parser") )
        std::cout << "Strict Ginac Parser enabled\n";
    if ( !ioption(_name="dim") )
        {
            std::cout << "dim=" << std::flush;
            std::cin >> dim;
            std::cout << std::flush;
        }
        else
            dim = ioption(_name="dim");

    // Load param list
    if ( !soption(_name="params").empty() )
        {
            boost::split(lst_params, params, boost::is_any_of(";"));
            parameters = symbols(lst_params);
        }

    // load exact
    if ( soption(_name="exact").empty() )
        {
            // Load param list
            std::cout << "params (eg params=\"a;b\")=" << std::flush; std::cin >> params;  std::cout << std::flush;
            if ( !params.empty() )
                {
                    boost::split(lst_params, params, boost::is_any_of(";"));
                    parameters = symbols(lst_params);
                }
            std::cout << parameters.size() << " params read\n";

            std::cout << "exact=" << std::flush;
            std::cin >> exact;
            std::cout << std::flush;
        }
    else
        exact = soption(_name="exact");

    switch (dim) {
    case(1) : {
        vars = symbols<1>();
        break;
    }
    case(2) : {
        vars = symbols<2>();
        break;
    }
    case(3) : {
        vars = symbols<3>();
        break;
    }
    default: {
        CHECK( false ) << "wrong dimension - should be lesser or egal to 3 - is equatl to " << dim << "\n";
    }
    }

#if 0
    if ( !parameters.size() )
        exact_parsed = parse(exact, vars);
    else
#endif
        exact_parsed = parse(exact, vars, parameters);

#if 0
    // retrieve symbols - not working because nops is not "recursive"
    std::cout << "Loading symbols from : " << exact_parsed << std::endl << std::flush;
    std::cout << "Contains " << exact_parsed.GiNaC::ex::nops() << " Ginac:ex objects\n" << std::flush;

    for (GiNaC::const_iterator i=exact_parsed.begin(); i!=exact_parsed.end(); ++i)
        {
            if ( GiNaC::is_a<symbol>(*i) )
                std::cout << "Found Symbol : " << GiNaC::ex_to<symbol>(*i).get_name() << std::endl;
        }
    std::cout << std::flush;

    //not working because ex is not mutable
    //boost::for_each(exact_parsed, [](GiNaC::ex const& e) 
    //{
    //if (GiNaC::is_a<symbol>(e)) std::cout << "Found Symbol : " <<  GiNaC::ex_to<symbol>(e).get_name() << "\n";
    //});
#endif

    // Retrieve each symbols using visitor
    std::cout << "Loading symbols from : " << exact_parsed << " (visitor)" << std::endl << std::flush;
    MyVisitor v;
    exact_parsed.traverse(v);
    std::list<std::string> symbols = v.symbol_names();
    boost::for_each(symbols, [](std::string const& s) {std::cout << "Found Symbol :" << s << std::endl;});

    // Substitute expressions
    std::cout << "Replace vars\n";
    boost::for_each(vars, [&exact_parsed](symbol const& sym)
                    {
                        ex f = exact_parsed;
                        std::cout << "Replace " << sym.get_name() << " in " << f << " :";
                        std::cout << substitute(f, sym, 1.0) << std::endl;
                    });

    std::cout << "Replace params\n";
    boost::for_each(parameters, [&exact_parsed](symbol const& sym)
                    {
                        ex f = exact_parsed;
                        if (f.has(sym))
                            {
                                std::cout << "Replace " << sym.get_name() << " in " << f << " :";
                                std::cout << substitute(f, sym, 1.0) << std::endl;
                            }
                    });

    std::cout << "Replace params using strings\n";
    boost::for_each(lst_params, [&exact_parsed](std::string const& s)
                    {
                        ex f = exact_parsed;
                        symbol sym;
                        //use MyVisitor class to find if symbol is present in f
                        MyVisitor f_v;
                        f.traverse(f_v);
                        if ( f_v.hassymbol(s, sym) )
                            {
                                std::cout << "Replace " << s << " in " << f << " :";
                                std::cout << substitute(f, sym, 1.0) << std::endl;
                            }
                    });

    std::cout << "Replace params using strings by function g()\n";
    boost::for_each(lst_params, [&exact_parsed, &exact, &vars, &parameters](std::string const& s)
                    {
                        ex f = exact_parsed;
                        ex f1;
                        symbol sym;
                        //use MyVisitor class to find if symbol is present in f
                        MyVisitor f_v;
                        f.traverse(f_v);
                        if ( f_v.hassymbol(s, sym) )
                            {
                                std::cout << "Replace " << s << " in " << f << " by :";
                                std::cin >> exact;
                                ex f1 = parse(exact, vars, parameters);
                                std::cout << std::flush;
                                std::cout << substitute(f, sym, f1) << std::endl;
                            }
                    });

    exact="sin(x)";
    auto g = expr(exact, vars);
}

template <typename ElementType>
void checkEqualElements( std::string const& info, ElementType const& u1, ElementType const& u2, double tol = 1e-10  )
{
    BOOST_TEST_MESSAGE( "checkEqualElements start  : "<< info );
    auto Xh = u1.functionSpace();
    for ( size_type k=0;k<Xh->nLocalDof();++k )
    {
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( u1(k), u2(k), tol );
#endif
        if ( std::abs(u1(k)-u2(k) ) > tol )
            break;
    }
    BOOST_TEST_MESSAGE( "checkEqualElements finish : "<< info );
}
void runTest1()
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto XhScalar = Pch<3>( mesh );
    auto XhVectorial = Pchv<3>( mesh );
    auto XhTensor2 = Pchm<3>( mesh );

    double a=2, b=4.5;
    std::map<std::string,double> params = { { "a", a }, { "b",b } };
    std::map<std::string,double> paramsV2 = { { "a", a } };

    auto exprScalarGinacV1 = expr( "a*x*x*cos(pi*y)+b*2*y*sin(x)*exp(x*x):x:y:a:b"/*, "scalarExpr"*/ );
    exprScalarGinacV1.setParameterValues( params );
    auto exprScalarGinacV2 = expr( exprScalarGinacV1, symbolExpr("b",cst(b)) );
    exprScalarGinacV2.setParameterValues( paramsV2 );
    auto exprScalarGinacV3 = expr( exprScalarGinacV1, symbolExpr("a",cst(a)), symbolExpr("b",cst(b)) );

    auto testExprScalar = hana::make_tuple( std::make_pair("V1 scalar",exprScalarGinacV1),
                                            std::make_pair("V2 scalar",exprScalarGinacV2),
                                            std::make_pair("V3 scalar",exprScalarGinacV3) );
    hana::for_each( testExprScalar, [&](auto const& t )
                    {
                        std::string testTag = t.first;
                        auto const& exprScalarGinac = t.second;
                        // id scalar
                        auto exprScalarFeel = a*Px()*Px()*cos(M_PI*Py()) + b*2*Py()*sin(Px())*exp(Px()*Px());
                        auto uFeel = XhScalar->element( exprScalarFeel );
                        auto uGinac = XhScalar->element( exprScalarGinac );
                        checkEqualElements( testTag + " id", uFeel, uGinac );
                        // grad scalar
                        auto exprScalarFeelGrad = vec( 2*a*Px()*cos(M_PI*Py()) + 4*b*Px()*Py()*exp(Px()*Px())*sin(Px()) + 2*b*Py()*exp(Px()*Px())*cos(Px()),
                                                       -M_PI*a*Px()*Px()*sin(M_PI*Py()) + 2*b*exp(Px()*Px())*sin(Px()) );
                        auto exprScalarGinacGrad = trans( grad<2>( exprScalarGinac/*, "scalarExprGrad"*/ ) );
                        auto uFeelGrad = XhVectorial->element( exprScalarFeelGrad );
                        auto uGinacGrad = XhVectorial->element( exprScalarGinacGrad );
                        checkEqualElements( testTag + " grad", uFeelGrad, uGinacGrad );
                        // laplacian scalar
                        auto exprScalarFeelLaplacian = 2*(a*cos(M_PI*Py())+4*b*Px()*Px()*Py()*exp(Px()*Px())*sin(Px()) + 4*b*Px()*Py()*exp(Px()*Px())*cos(Px()) + b*Py()*exp(Px()*Px())*sin(Px()))  - M_PI*M_PI*a*Px()*Px()*cos(M_PI*Py());
                        auto exprScalarGinacLaplacian = laplacian( exprScalarGinac/*, "scalarExprLaplacian"*/ );
                        auto uFeelLaplacian = XhScalar->element( exprScalarFeelLaplacian );
                        auto uGinacLaplacian = XhScalar->element( exprScalarGinacLaplacian );
                        checkEqualElements( testTag + " laplacian", uFeelLaplacian, uGinacLaplacian );
                        // diff scalar
                        auto exprScalarFeelDiff = vec( Px()*Px()*cos(M_PI*Py()), 2*Py()*exp(Px()*Px())*sin(Px()) );
                        auto exprScalarGinacDiffA = diff( exprScalarGinac, "a" );
                        auto exprScalarGinacDiffB = diff( exprScalarGinac, "b" );
                        auto uFeelDiff = XhVectorial->element( exprScalarFeelDiff );
                        auto uGinacDiff = XhVectorial->element( vec( exprScalarGinacDiffA, exprScalarGinacDiffB ) );
                        checkEqualElements( testTag + " diff", uFeelDiff, uGinacDiff );
                        auto exprScalarGinacDiffX = diff( exprScalarGinac, "x" );
                        auto exprScalarGinacDiffY = diff( exprScalarGinac, "y" );
                        auto uGinacDiffX_Y = XhVectorial->element( vec( exprScalarGinacDiffX, exprScalarGinacDiffY ) );
                        checkEqualElements( testTag + " diff (grad)", uFeelGrad, uGinacDiffX_Y );
                        auto exprScalarGinacDiffXX = diff( exprScalarGinac, "x", 2 );
                        auto exprScalarGinacDiffYY = diff( exprScalarGinac, "y", 2 );
                        auto uGinacDiffXX_YY = XhScalar->element( exprScalarGinacDiffXX + exprScalarGinacDiffYY );
                        checkEqualElements( testTag + " diff (laplacian)", uFeelLaplacian, uGinacDiffXX_YY );
                    });



    auto exprVectorialGinacV1 = expr<2,1>( "{a*x*x*cos(pi*y),b*2*y*sin(x)*exp(x*x)}:x:y:a:b"/*, "vectorialExpr"*/ );
    exprVectorialGinacV1.setParameterValues( params );
    auto exprVectorialGinacV2 = expr( exprVectorialGinacV1, symbolExpr("b",cst(b)) );
    exprVectorialGinacV2.setParameterValues( paramsV2 );
    auto exprVectorialGinacV3 = expr( exprVectorialGinacV1, symbolExpr("a",cst(a)), symbolExpr("b",cst(b)) );

    auto testExprVectorial = hana::make_tuple( std::make_pair("V1 vectorial",exprVectorialGinacV1),
                                               std::make_pair("V2 vectorial",exprVectorialGinacV2),
                                               std::make_pair("V3 vectorial",exprVectorialGinacV3) );
    hana::for_each( testExprVectorial, [&](auto const& t )
                    {
                        std::string testTag = t.first;
                        auto const& exprVectorialGinac = t.second;
                        // id vectorial
                        auto exprVectorialFeel = vec( a*Px()*Px()*cos(M_PI*Py()),
                                                      b*2*Py()*sin(Px())*exp(Px()*Px()) );
                        auto uVectorialFeel = XhVectorial->element( exprVectorialFeel );
                        auto uVectorialGinac = XhVectorial->element( exprVectorialGinac );
                        checkEqualElements( testTag+" id", uVectorialFeel, uVectorialGinac );
                        // div vectorial
                        auto exprVectorialFeelDiv = 2*a*Px()*cos(M_PI*Py()) + 2*b*sin(Px())*exp(Px()*Px());
                        auto exprVectorialGinacDiv = div/*<2>*/( exprVectorialGinac/*, "vectorialExprDiv"*/ );
                        auto uVectorialFeelDiv = XhScalar->element( exprVectorialFeelDiv );
                        auto uVectorialGinacDiv = XhScalar->element( exprVectorialGinacDiv );
                        checkEqualElements( testTag+" div", uVectorialFeelDiv, uVectorialGinacDiv );
                        // grad vectorial
                        auto exprVectorialFeelGrad = mat<2,2>( 2*a*Px()*cos(M_PI*Py()), -M_PI*a*Px()*Px()*sin(M_PI*Py()),
                                                               4*b*Px()*Py()*exp(Px()*Px())*sin(Px()) + 2*b*Py()*exp(Px()*Px())*cos(Px()), 2*b*sin(Px())*exp(Px()*Px()) );
                        auto exprVectorialGinacGrad = grad<2>( exprVectorialGinac );
                        auto uVectorialFeelGrad = XhTensor2->element( exprVectorialFeelGrad );
                        auto uVectorialGinacGrad = XhTensor2->element( exprVectorialGinacGrad );
                        checkEqualElements( testTag+" grad", uVectorialFeelGrad, uVectorialGinacGrad );
                        // laplacian vectorial
                        auto exprVectorialFeelLaplacian = vec( -M_PI*M_PI*a*Px()*Px()*cos(M_PI*Py()) + 2*a*cos(M_PI*Py()),
                                                               2*b*Py()*(4*Px()*Px()*sin(Px())+4*Px()*cos(Px())+sin(Px()))*exp(Px()*Px()) );
                        auto exprVectorialGinacLaplacian = laplacian( exprVectorialGinac/*, "toto"*/ );
                        auto uVectorialFeelLaplacian = XhVectorial->element( exprVectorialFeelLaplacian );
                        auto uVectorialGinacLaplacian = XhVectorial->element( exprVectorialGinacLaplacian );
                        checkEqualElements( testTag+" laplacian", uVectorialFeelLaplacian, uVectorialGinacLaplacian, 1e-8 );
                        // diff vectorial
                        auto exprVectorialFeelDiffA = vec(Px()*Px()*cos(M_PI*Py()),cst(0.) );
                        auto exprVectorialFeelDiffB = vec(cst(0.),2*Py()*sin(Px())*exp(Px()*Px()) );
                        auto exprVectorialGinacDiffA = diff( exprVectorialGinac, "a" );
                        auto exprVectorialGinacDiffB = diff( exprVectorialGinac, "b" );
                        auto uVectorialFeelDiffA = XhVectorial->element( exprVectorialFeelDiffA );
                        auto uVectorialGinacDiffA = XhVectorial->element( exprVectorialGinacDiffA );
                        auto uVectorialFeelDiffB = XhVectorial->element( exprVectorialFeelDiffB );
                        auto uVectorialGinacDiffB = XhVectorial->element( exprVectorialGinacDiffB );
                        checkEqualElements( testTag+" diff a", uVectorialFeelDiffA, uVectorialGinacDiffA, 1e-8 );
                        checkEqualElements( testTag+" diff b", uVectorialFeelDiffB, uVectorialGinacDiffB, 1e-8 );
                    });

    // id vf
    auto exprScalarField = Px()*Px()*cos(M_PI*Py()) + 2*Py()*sin(Px())*exp(Px()*Px());
    auto scalarField = XhScalar->element( exprScalarField );
    auto exprScalarFeelVF = a*2*Px()*sin(Py()*b)*exp(pow(idv(scalarField),2));
    auto exprScalarGinacVF = expr( "a*2*x*sin(y*b)*exp(u^2):x:y:a:b:u", "u", idv(scalarField) );
    exprScalarGinacVF.setParameterValues( params );
    auto uScalarFeelVF = XhScalar->element( exprScalarFeelVF );
    auto uScalarFeelGinacVF = XhScalar->element( exprScalarGinacVF );
    checkEqualElements( "id vf", uScalarFeelVF, uScalarFeelGinacVF );
    // diff vf
    auto exprScalarFeelVFDiff = a*2*Px()*sin(Py()*b)*2*idv(scalarField)*exp(pow(idv(scalarField),2));
    auto exprScalarGinacVFDiffStr = diff( "a*2*x*sin(y*b)*exp(u^2):x:y:a:b:u", "u", "u", idv(scalarField) );
    exprScalarGinacVFDiffStr.setParameterValues( params );
    auto uScalarFeelVFDiff = XhScalar->element( exprScalarFeelVFDiff );
    auto uScalarGinacVFDiffStr = XhScalar->element( exprScalarGinacVFDiffStr );
    checkEqualElements( "diff vf str", uScalarFeelVFDiff, uScalarGinacVFDiffStr );
#if 1
    // diff vf str
    auto exprScalarGinacVFDiff = diff( exprScalarGinacVF , "u" );
    auto uScalarGinacVFDiff = XhScalar->element( exprScalarGinacVFDiff );
    checkEqualElements( "diff vf", uScalarFeelVFDiff, uScalarGinacVFDiff );
#endif
}

void runTest2()
{
    using GiNaC::symbol;
    using GiNaC::ex;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    double area = integrate( _range=elements(mesh), _expr=cst(1.) ).evaluate()(0,0);
    std::vector<symbol> vars = symbols<2>();

    std::string ex1_s = "1";
    std::string ex2_s = "2";

    ex ex1 = parse(ex1_s, vars);
    ex ex2 = parse(ex2_s, vars);

    double int_expr1 = area;
    double int_expr2 = 2*area;

    // Expressions computed from a string
    auto expr1_1 = expr(ex1_s, vars, "ex1_s_file");
    auto expr2_1 = expr(ex2_s, vars, "ex2_s_file");

    double int_expr1_string = integrate( _range=elements(mesh), _expr=expr1_1 ).evaluate()(0,0);
#if defined(USE_BOOST_TEST)
    BOOST_CHECK_CLOSE( int_expr1_string, int_expr1, 1e-10 );
#endif
    double int_expr2_string = integrate( _range=elements(mesh), _expr=expr2_1 ).evaluate()(0,0);
#if defined(USE_BOOST_TEST)
    BOOST_CHECK_CLOSE( int_expr2_string, int_expr2, 1e-10 );
#endif
    // Expressions computed from a string
    auto expr1_2 = expr(ex1, vars, "ex1_file");
    auto expr2_2 = expr(ex2, vars, "ex2_file");

    double int_expr1_ex = integrate( _range=elements(mesh), _expr=expr1_2 ).evaluate()(0,0);
#if defined(USE_BOOST_TEST)
    BOOST_CHECK_CLOSE( int_expr1_ex, int_expr1, 1e-10 );
#endif
    double int_expr2_ex = integrate( _range=elements(mesh), _expr=expr2_2 ).evaluate()(0,0);
#if defined(USE_BOOST_TEST)
    BOOST_CHECK_CLOSE( int_expr2_ex, int_expr2, 1e-10 );
#endif
}

void runTest3()
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2,1>>);
    auto Vh1 = Pch<1>( mesh );
    auto Vh2 = Pch<2>( mesh );
    auto Vh3 = Pch<3>( mesh );
    auto u1 = Vh1->element(Px());
    auto u1b = Vh1->element(Py());
    auto u2 = Vh2->element(Py());
    auto u3 = Vh3->element(Px()*Px()*Py());

    // scalar expressions
    auto exprA = expr("2*x*y+x*x*y:x:y");
    auto exprB = expr("2*u*y+x*x*y:u:x:y");
    auto exprBvf = expr(exprB, symbolExpr("u",idv(u1)) );

    auto exprC = expr("2*u*v+x*x*y:u:v:x:y");
    auto exprCvf = expr(exprC,symbolExpr("u",idv(u1)), symbolExpr("v",idv(u2)) );
    auto exprD = expr("2*u*v+w:u:v:w");
    auto exprDvf = expr(exprD,symbolExpr("u",idv(u1)), symbolExpr("v",idv(u2) ), symbolExpr("w",idv(u3) ) );
    auto exprE = expr("2*u*v+w:u:v:w");
    auto exprEvf = expr(exprE, symbolsExpr( symbolExpr("u",idv(u1)), symbolExpr("v",idv(u2) ), symbolExpr("w",idv(u1)*idv(u1)*idv(u2) ) ) );

    auto exprF = expr("2*u*v+x*u*v:u:v:x:y");
    //auto exprFvf = expr(exprC,symbolExpr( { { "u",idv(u1) }, { "v",idv(u1b) } }) );
    //auto exprFvf = expr(exprC,symbolExpr( { std::make_pair("u",idv(u1)) , std::make_pair("v",idv(u1b)) }) );
    auto exprFvf = expr(exprC,symbolExpr( { std::make_pair( std::string("u"),idv(u1)) , std::make_pair( std::string("v"),idv(u1b)) }) );

    double intA = integrate(_range=elements(mesh),_expr=exprA,_quad=5).evaluate()(0,0);
    double intB = integrate(_range=elements(mesh),_expr=exprBvf,_quad=5).evaluate()(0,0);
    double intC = integrate(_range=elements(mesh),_expr=exprCvf,_quad=5).evaluate()(0,0);
    double intD = integrate(_range=elements(mesh),_expr=exprDvf,_quad=5).evaluate()(0,0);
    double intE = integrate(_range=elements(mesh),_expr=exprEvf,_quad=5).evaluate()(0,0);
    double intF = integrate(_range=elements(mesh),_expr=exprFvf,_quad=5).evaluate()(0,0);

    BOOST_CHECK_CLOSE( intA, intB, 1e-8 );
    BOOST_CHECK_CLOSE( intA, intC, 1e-8 );
    BOOST_CHECK_CLOSE( intA, intD, 1e-8 );
    BOOST_CHECK_CLOSE( intA, intE, 1e-8 );
    BOOST_CHECK_CLOSE( intA, intF, 1e-8 );

    // vectorial expressions
    auto exprVecA = expr<2,1>("{2*x*y+x*x*y,4*x-y}:x:y");
    auto exprVecB = expr<2,1>("{2*u*y+u*x*y,4*u-y}:u:x:y");
    auto exprVecBvf = expr(exprVecB, symbolExpr( "u",idv(u1) ) );
    auto exprVecC = expr<2,1>("{2*u*v+x*x*y,4*u-v}:u:v:x:y");
    auto exprVecCvf = expr(exprVecC, symbolExpr( "u",idv(u1) ), symbolExpr("v",idv(u2) ) );
    auto exprVecD = expr<2,1>("{2*u*v+w,4*u-v}:u:v:w");
    auto exprVecDvf = expr(exprVecD, symbolExpr("u",idv(u1)), symbolExpr("v",idv(u2)), symbolExpr("w",idv(u3)) );
    auto exprVecE = expr<2,1>("{2*u*v+w,4*u-v}:u:v:w");
    auto exprVecEvf = expr(exprVecE,symbolsExpr( symbolExpr( "u",idv(u1) ), symbolExpr("v",idv(u2) ), symbolExpr("w",idv(u1)*idv(u1)*idv(u2) ) ) );

    double intVecA = integrate(_range=elements(mesh),_expr=inner(exprVecA),_quad=5).evaluate()(0,0);
    double intVecB = integrate(_range=elements(mesh),_expr=inner(exprVecBvf),_quad=5).evaluate()(0,0);
    double intVecC = integrate(_range=elements(mesh),_expr=inner(exprVecCvf),_quad=5).evaluate()(0,0);
    double intVecD = integrate(_range=elements(mesh),_expr=inner(exprVecDvf),_quad=5).evaluate()(0,0);
    double intVecE = integrate(_range=elements(mesh),_expr=inner(exprVecEvf),_quad=5).evaluate()(0,0);

    BOOST_CHECK_CLOSE( intVecA, intVecB, 1e-8 );
    BOOST_CHECK_CLOSE( intVecA, intVecC, 1e-8 );
    BOOST_CHECK_CLOSE( intVecA, intVecD, 1e-8 );
    BOOST_CHECK_CLOSE( intVecA, intVecE, 1e-8 );
}

void runTest4()
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2,1>>);
    auto Vh = Pch<3>( mesh );
    auto u = Vh->element( Px()*Px()*Py() );

    std::map<std::string,double> params;
    params["a"]=2;
    params["b"]=4.5;

    // constant
    double cstValue = 3.14;
    auto exprCst = expr( cstValue );
    BOOST_CHECK( exprCst.expression().isConstant() );
    BOOST_CHECK( exprCst.isPolynomial() );
    BOOST_CHECK( exprCst.polynomialOrder() == 0 );
    BOOST_CHECK_CLOSE( exprCst.evaluate()(0,0), cstValue, 1e-12 );
    double intCst1 = integrate(_range=elements(mesh),_expr=cst(cstValue)).evaluate()(0,0);
    double intCst2 = integrate(_range=elements(mesh),_expr=exprCst).evaluate()(0,0);
    BOOST_CHECK_CLOSE( intCst1, intCst2, 1e-12 );

    // scalar
    auto exprA = expr("2*a*b:a:b");
    exprA.setParameterValues( params );
    BOOST_CHECK( exprA.expression().isConstant() );
    BOOST_CHECK( exprA.isPolynomial() );
    BOOST_CHECK( exprA.polynomialOrder() == 0 );
    BOOST_CHECK_CLOSE( exprA.evaluate()(0,0), 2*2*4.5, 1e-12 );
    double intA1 = integrate(_range=elements(mesh),_expr=cst(2*2*4.5)).evaluate()(0,0);
    double intA2 = integrate(_range=elements(mesh),_expr=exprA).evaluate()(0,0);
    BOOST_CHECK_CLOSE( intA1, intA2, 1e-12 );

    auto exprAu = expr(exprA, symbolsExpr( symbolExpr( "u",idv(u) ) ) );
    BOOST_CHECK( exprAu.expression().isConstant() );
    BOOST_CHECK( exprAu.isPolynomial() );
    BOOST_CHECK( exprAu.polynomialOrder() == 0 );
    BOOST_CHECK_CLOSE( exprAu.evaluate()(0,0), 2*2*4.5, 1e-12 );
    double intAu2 = integrate(_range=elements(mesh),_expr=exprAu).evaluate()(0,0);
    BOOST_CHECK_CLOSE( intA1, intAu2, 1e-12 );

    auto exprB = expr("2*a*b+x^2*y+(x+y)^4:a:b:x:y");
    exprB.setParameterValues( params );
    BOOST_CHECK( !exprB.expression().isConstant() );
    BOOST_CHECK( exprB.isPolynomial() );
    BOOST_CHECK( exprB.polynomialOrder() == 4 );

    auto exprC = expr<3>("2*a+x^2*y+(x+y)^4+cos(y):a:x:y");
    exprC.setParameterValues( params );
    BOOST_CHECK( !exprC.expression().isConstant() );
    BOOST_CHECK( !exprC.isPolynomial() );
    BOOST_CHECK( exprC.polynomialOrder() == 3 );


    auto exprDtmp = expr("2*a+x^2*y+u*(x+y)^4:a:x:y:u");
    auto exprD = expr( exprDtmp,symbolsExpr( symbolExpr( "u",idv(u) ) ) );
    exprD.setParameterValues( params );
    BOOST_CHECK( !exprD.expression().isConstant() );
    BOOST_CHECK( exprD.isPolynomial() );
    BOOST_CHECK( exprD.polynomialOrder() == 7 );

    auto exprEtmp = expr<3>("2*a+x^2*y+u*(x+y)^4:a:x:y:u");
    auto exprE = expr( exprEtmp,symbolsExpr( symbolExpr( "u",cos(Px()) ) ) );
    BOOST_CHECK( !exprE.expression().isConstant() );
    BOOST_CHECK( !exprE.isPolynomial() );
    BOOST_CHECK( exprE.polynomialOrder() == 3 );

    auto exprF = expr("4*sin(3*Pi/2)*exp(3):x:a:u");
    BOOST_CHECK( exprF.expression().isConstant() );
    BOOST_CHECK( exprF.isPolynomial() );
    BOOST_CHECK( exprF.polynomialOrder() == 0 );
    double evalFvalue = 4*std::sin(3*M_PI/2.)*std::exp(3.);
    BOOST_CHECK_CLOSE( exprF.evaluate()(0,0), evalFvalue, 1e-12 );
    double intF1 = integrate(_range=elements(mesh),_expr=cst(evalFvalue)).evaluate()(0,0);
    double intF2 = integrate(_range=elements(mesh),_expr=exprF).evaluate()(0,0);
    BOOST_CHECK_CLOSE( intF1, intF2, 1e-12 );
    auto exprFu = expr( exprF,symbolExpr( "u",idv(u) ) );
    BOOST_CHECK( exprFu.expression().isConstant() );
    BOOST_CHECK( exprFu.isPolynomial() );
    BOOST_CHECK( exprFu.polynomialOrder() == 0 );
    BOOST_CHECK_CLOSE( exprFu.evaluate()(0,0), evalFvalue, 1e-12 );

    // vectorial
    auto exprVA = expr<2,1>("{2*a*b,a+b}:a:b");
    exprVA.setParameterValues( params );
    BOOST_CHECK( exprVA.expression().isConstant() );
    BOOST_CHECK( exprVA.isPolynomial() );
    BOOST_CHECK( exprVA.polynomialOrder() == 0 );
    auto evalVA = exprVA.evaluate();
    BOOST_CHECK_CLOSE( evalVA(0,0), 2*2*4.5, 1e-12 );
    BOOST_CHECK_CLOSE( evalVA(1,0), 2+4.5, 1e-12 );
    auto intVA1 = integrate(_range=elements(mesh),_expr=vec(cst(2*2*4.5),cst(2+4.5))).evaluate();
    auto intVA2 = integrate(_range=elements(mesh),_expr=exprVA).evaluate();
    BOOST_CHECK_CLOSE( intVA1(0,0), intVA2(0,0), 1e-12 );
    BOOST_CHECK_CLOSE( intVA1(1,0), intVA2(1,0), 1e-12 );

    auto exprVB = expr<2,1>("{2*a*b+x^2*y,(x+y)^4}:a:b:x:y");
    exprVB.setParameterValues( params );
    BOOST_CHECK( !exprVB.expression().isConstant() );
    BOOST_CHECK( exprVB.isPolynomial() );
    BOOST_CHECK( exprVB.polynomialOrder() == 4 );

    auto exprVC = expr<2,1,3>("{2*a*b+x^2*cos(y),(x+y)^4}:a:b:x:y");
    exprVC.setParameterValues( params );
    BOOST_CHECK( !exprVC.expression().isConstant() );
    BOOST_CHECK( !exprVC.isPolynomial() );
    BOOST_CHECK( exprVC.polynomialOrder() == 3 );

    auto exprVDtmp = expr<2,1>("{2*a*b+x^2*y*u,(x+y)^4}:a:b:x:y:u");
    auto exprVD = expr( exprVDtmp,symbolsExpr( symbolExpr( "u",idv(u) ) ) );
    exprVD.setParameterValues( params );
    BOOST_CHECK( !exprVD.expression().isConstant() );
    BOOST_CHECK( exprVD.isPolynomial() );
    BOOST_CHECK( exprVD.polynomialOrder() == 6 );

    auto exprVEtmp = expr<2,1,3>("{2*a*b+x^2*y*u,(x+y)^4}:a:b:x:y:u");
    auto exprVE = expr( exprVEtmp,symbolsExpr( symbolExpr( "u",cos(Px()) ) ) );
    exprVE.setParameterValues( params );
    BOOST_CHECK( !exprVE.expression().isConstant() );
    BOOST_CHECK( !exprVE.isPolynomial() );
    BOOST_CHECK( exprVE.polynomialOrder() == 3 );

    // tensor2
    auto exprMA = expr<2,2>("{2*a*b,a+b,a-b,a/b}:a:b");
    exprMA.setParameterValues( params );
    BOOST_CHECK( exprMA.expression().isConstant() );
    BOOST_CHECK( exprMA.isPolynomial() );
    BOOST_CHECK( exprMA.polynomialOrder() == 0 );
    auto evalMA = exprMA.evaluate();
    BOOST_CHECK_CLOSE( evalMA(0,0), 2*2*4.5, 1e-12 );
    BOOST_CHECK_CLOSE( evalMA(0,1), 2+4.5, 1e-12 );
    BOOST_CHECK_CLOSE( evalMA(1,0), 2-4.5, 1e-12 );
    BOOST_CHECK_CLOSE( evalMA(1,1), 2/4.5, 1e-12 );
    auto intMA1 = integrate(_range=elements(mesh),_expr=mat<2,2>(cst(2*2*4.5),cst(2+4.5),cst(2-4.5),cst(2/4.5))).evaluate();
    auto intMA2 = integrate(_range=elements(mesh),_expr=exprMA).evaluate();
    for ( int i=0;i<2;++i )
        for ( int j=0;j<2;++j )
            BOOST_CHECK_CLOSE( intMA1(i,j), intMA2(i,j), 1e-12 );


    auto exprMB = expr<2,2>("{cos(Pi/3),exp(2),sin(2),4*16}:a:x");
    BOOST_CHECK( exprMB.expression().isConstant() );
    BOOST_CHECK( exprMB.isPolynomial() );
    BOOST_CHECK( exprMB.polynomialOrder() == 0 );
    auto evalMB = exprMB.evaluate();
    BOOST_CHECK_CLOSE( evalMB(0,0), std::cos(M_PI/3.), 1e-12 );
    BOOST_CHECK_CLOSE( evalMB(0,1), std::exp(2.), 1e-12 );
    BOOST_CHECK_CLOSE( evalMB(1,0), std::sin(2.), 1e-12 );
    BOOST_CHECK_CLOSE( evalMB(1,1), 4*16, 1e-12 );
    auto evalMBvalue = mat<2,2>( cst(std::cos(M_PI/3.)), cst(std::exp(2.)),
                                 cst(std::sin(2.)), cst(4.*16) );
    auto intMB1 = integrate(_range=elements(mesh),_expr=evalMBvalue).evaluate();
    auto intMB2 = integrate(_range=elements(mesh),_expr=exprMB).evaluate();
    for ( int i=0;i<2;++i )
        for ( int j=0;j<2;++j )
            BOOST_CHECK_CLOSE( intMB1(i,j), intMB2(i,j), 1e-12 );
}

void runTest5()
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2,1>>);
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();

    auto exprScalar1 = expr("2*u*v_2+v_0+v_1:u:v_0:v_1:v_2");
    auto exprScalar1vf = expr(exprScalar1,symbolExpr("u",cst(3.14)), symbolExpr("v",vec(cst(1.),cst(2.),cst(3.)), SymbolExprComponentSuffix( 3,1 ) ) );
    BOOST_CHECK_CLOSE( exprScalar1vf.evaluate()(0,0), 21.84, 1e-12 );
    u.on(_range=elements(mesh),_expr=exprScalar1vf);
    BOOST_CHECK_CLOSE( u.max(), 21.84, 1e-12 );

    auto exprScalar2 = expr("2*u*v_z+v_x+v_y:u:v_x:v_y:v_z");
    auto exprScalar2vf = expr(exprScalar2,symbolExpr("u",cst(3.14)), symbolExpr("v",vec(cst(1.),cst(2.),cst(3.)), SymbolExprComponentSuffix( 3,1,true ) ) );
    BOOST_CHECK_CLOSE( exprScalar2vf.evaluate()(0,0), 21.84, 1e-12 );
    u.on(_range=elements(mesh),_expr=exprScalar2vf);
    BOOST_CHECK_CLOSE( u.max(), 21.84, 1e-12 );


    auto exprScalar3 = expr("2*u*v_xz+v_xx+v_yz:u:v_xx:v_xz:v_yz");
    auto exprScalar3vf = expr(exprScalar3,symbolExpr("u",cst(3.14)), symbolExpr("v",mat<2,3>(cst(1.),cst(2.),cst(3.),
                                                                                             cst(4.),cst(5.),cst(6.)), SymbolExprComponentSuffix( 2,3, true ) ) );
    BOOST_CHECK_CLOSE( exprScalar3vf.evaluate()(0,0), 25.84, 1e-12 );
    u.on(_range=elements(mesh),_expr=exprScalar3vf);
    BOOST_CHECK_CLOSE( u.max(), 25.84, 1e-12 );

    auto exprVec1 = expr<2,1>("{2*u*v_z+v_x, 1+u+v_y}:u:v_x:v_y:v_z");
    auto exprVec1vf = expr(exprVec1,symbolExpr("u",cst(3.14)), symbolExpr("v",vec(cst(1.),cst(2.),cst(3.)), SymbolExprComponentSuffix( 3,1,true ) ) );
    BOOST_CHECK_CLOSE( exprVec1vf.evaluate()(0,0), 19.84, 1e-12 );
    BOOST_CHECK_CLOSE( exprVec1vf.evaluate()(1,0), 6.14, 1e-12 );
    u.on(_range=elements(mesh),_expr=exprVec1vf(0,0));
    BOOST_CHECK_CLOSE( u.max(), 19.84, 1e-12 );
    u.on(_range=elements(mesh),_expr=exprVec1vf(1,0));
    BOOST_CHECK_CLOSE( u.max(), 6.14, 1e-12 );
}

void runTest6()
{
    auto exprScalar1 = expr("2+u:u");
    auto exprScalar2 = expr("3*u+v:u:v");
    auto exprScalar3 = expr("v+w:v:w");
    auto exprScalar4 = expr("xx+v+w:xx:v:w");
    //v=2+3=5 w=3*3+5=14 xx=5+14=19 yy=19+5+14=38
    auto thesymbolExprA = symbolsExpr( symbolExpr("u",cst(3.0)),
                                       symbolExpr("v",exprScalar1),
                                       symbolExpr("w",exprScalar2),
                                       symbolExpr("xx",exprScalar3),
                                       symbolExpr("yy",exprScalar4)
                                       );
    auto exprScalar3vf = expr(exprScalar3,thesymbolExprA);
    BOOST_CHECK_CLOSE( exprScalar3vf.evaluate()(0,0), 19, 1e-12 );
    auto exprScalar4vf = expr(exprScalar4,thesymbolExprA);
    BOOST_CHECK_CLOSE( exprScalar4vf.evaluate()(0,0), 38, 1e-12 );

    auto exprVectorial1 = expr<2,1>("{v+w,v-w}:v:w");
    auto exprVectorial2 = expr<2,1>("{xx_x+v+w,xx_y+3*(v-w)}:xx_x:xx_y:v:w");
    //v=2+3=5 w=3*3+5=14 xx_x=5+14=19 xx_y=5-14=-9   yy_x=19+5+14=38   yy_y=19+3*(5-14)=-36
    auto thesymbolExprB = symbolsExpr( symbolExpr("u",cst(3.0)),
                                       symbolExpr("v",exprScalar1),
                                       symbolExpr("w",exprScalar2),
                                       symbolExpr("xx",exprVectorial1, SymbolExprComponentSuffix( 2,1,true )),
                                       symbolExpr("yy",exprVectorial2, SymbolExprComponentSuffix( 2,1,true ) )
                                       );

    auto exprVectorial1vf = expr(exprVectorial1,thesymbolExprB);
    BOOST_CHECK_CLOSE( exprVectorial1vf.evaluate()(0,0), 19, 1e-12 );
    BOOST_CHECK_CLOSE( exprVectorial1vf.evaluate()(1,0), -9, 1e-12 );
    auto exprVectorial2vf = expr(exprVectorial2,thesymbolExprB);
    BOOST_CHECK_CLOSE( exprVectorial2vf.evaluate()(0,0), 38, 1e-12 );
    BOOST_CHECK_CLOSE( exprVectorial2vf.evaluate()(1,0), -36, 1e-12 );

    using mesh_t = Mesh<Simplex<2>>;
    auto mesh = loadMesh( _mesh = new mesh_t );
    double evalIntegrate4 = integrate(_range=elements(mesh),_expr=exprScalar4vf ).evaluate()(0,0);
    double evalIntegrate4_check = integrate(_range=elements(mesh),_expr=cst(38.) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE( evalIntegrate4, evalIntegrate4_check, 1e-12 );
    auto evalIntegrateVectorial2 = integrate(_range=elements(mesh),_expr=exprVectorial2vf ).evaluate();
    auto evalIntegrateVectorial2_check = integrate(_range=elements(mesh),_expr=vec(cst(38.),cst(-36.)) ).evaluate();
    BOOST_CHECK_CLOSE( evalIntegrateVectorial2(0,0), evalIntegrateVectorial2_check(0,0), 1e-12 );
    BOOST_CHECK_CLOSE( evalIntegrateVectorial2(1,0), evalIntegrateVectorial2_check(1,0), 1e-12 );
}

template <typename ElementType>
struct MyUpdateFunctor
{
    MyUpdateFunctor( ElementType & u ) : M_u( u ), M_t( 0 ) {}

    void setValue( double t ) { M_t = t; }

    void update() { M_u.setConstant( M_t ); }

private :
    ElementType & M_u;
    double M_t;
};

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( inner_suite )
BOOST_AUTO_TEST_CASE( test_0 )
{
    runTest0();
}
BOOST_AUTO_TEST_CASE( test_1 )
{
    runTest1();
}
BOOST_AUTO_TEST_CASE( test_2 )
{
    runTest2();
}
BOOST_AUTO_TEST_CASE( test_3 )
{
    runTest3();
}
BOOST_AUTO_TEST_CASE( test_4 )
{
    runTest4();
}
BOOST_AUTO_TEST_CASE( test_5 )
{
    runTest5();
}
BOOST_AUTO_TEST_CASE( test_6 )
{
    runTest6();
}
BOOST_AUTO_TEST_CASE( test_7 )
{
    auto a1 = expr( "3*u+v:u:v");
    a1.setParameterValues( { { "u", 2 }, { "v",5 } } );
    BOOST_CHECK_CLOSE( a1.evaluate()(0,0), 11, 1e-12 );
    auto a2 = expr( "{3*u+v}:u:v");
    a2.setParameterValues( { { "u", 2 }, { "v",5 } } );
    BOOST_CHECK_CLOSE( a2.evaluate()(0,0), 11, 1e-12 );
    auto a3 = expr<1,1>( "3*u+v:u:v");
    a3.setParameterValues( { { "u", 2 }, { "v",5 } } );
    BOOST_CHECK_CLOSE( a3.evaluate()(0,0), 11, 1e-12 );
    auto a4 = expr<1,1>( "{3*u+v}:u:v");
    a4.setParameterValues( { { "u", 2 }, { "v",5 } } );
    BOOST_CHECK_CLOSE( a4.evaluate()(0,0), 11, 1e-12 );
    auto a5 = expr<1,1>( "{3*u*x^2+2*v*y^2}:u:v:x:y");
    a5.setParameterValues( { { "u", 2 }, { "v",5 } } );
    auto lap_a5 = laplacian( a5 );
    BOOST_CHECK_CLOSE( lap_a5.evaluate()(0,0), 32, 1e-12 );
}

BOOST_AUTO_TEST_CASE( test_8 )
{
    auto se = symbolsExpr( symbolExpr( "u", _e1 ), symbolExpr( "v", _e2 ) );
    auto a1b = expr("2*u+v*w:u:v:w");
    a1b.setParameterValues( { { "w", 2 } } );
    auto a1l = expr( a1b, se );
    auto a1 = a1l( cst(3.),cst(5.) );
    BOOST_CHECK_CLOSE( a1.evaluate()(0,0), 16, 1e-12 );
}
BOOST_AUTO_TEST_CASE( test_symbolsexpr_update_function )
{
    using mesh_t = Mesh<Simplex<2, 1>>;
    auto mesh = loadMesh( _mesh = new mesh_t );
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    double t=0;
    auto up_lambda = [&u,&t]() {
        u.setConstant( t );
    };

    using element_type = std::decay_t<decltype(v)>;
    MyUpdateFunctor<element_type> up_functor( v );

    auto se = symbolsExpr( symbolExpr( "k", cst(3.) ),
                           symbolExpr( "u", idv(u),SymbolExprComponentSuffix(1,1), up_lambda ),
                           symbolExpr( "v", idv(v),SymbolExprComponentSuffix(1,1), std::bind( &MyUpdateFunctor<element_type>::update, &up_functor ) )
                           );

    for ( ; t<1 ; t+=0.1 )
    {
        double int_exact = integrate(_range=elements(mesh),_expr=cst(t) ).evaluate()(0,0);
        double int_u = integrate(_range=elements(mesh),_expr=expr(expr("u:u"),se) ).evaluate()(0,0);
        BOOST_CHECK_CLOSE( int_u,int_exact, 1e-10 );

        up_functor.setValue( t );
        double int_v = integrate(_range=elements(mesh),_expr=expr(expr("v:v"),se) ).evaluate()(0,0);
        BOOST_CHECK_CLOSE( int_v,int_exact, 1e-10 );
    }
}

BOOST_AUTO_TEST_CASE( test_mod )
{
    auto se = symbolsExpr( symbolExpr( "u", _e1 ), symbolExpr( "v", _e2 ) );
    auto a1b = expr("mod(u,v):u:v");
    a1b.setParameterValues( { { "u", 2 }, { "v", 1 }} );
    BOOST_CHECK_SMALL( a1b.evaluate()(0,0), 1e-12 );
    a1b.setParameterValues( { { "u", 3 }, { "v", 6 }} );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 3, 1e-12 );
    a1b.setParameterValues( { { "u", 6.1 }, { "v", 3 }} );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 0.1, 1e-12 );
}
BOOST_AUTO_TEST_CASE( test_step1 )
{
    auto se = symbolsExpr( symbolExpr( "u", _e1 ), symbolExpr( "v", _e2 ) );
    auto a1b = expr("step1(u,v):u:v");
    a1b.setParameterValues( { { "u", 0 }, { "v", 0.5 }} );
    BOOST_CHECK_SMALL( a1b.evaluate()(0,0), 1e-12 );
    a1b.setParameterValues( { { "u", 1 }, { "v", 0.5 }} );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 1, 1e-12 );
    a1b.setParameterValues( { { "u", 0.49999 }, { "v", 0.5 }} );
    BOOST_CHECK_SMALL( a1b.evaluate()(0,0), 1e-12 );
    a1b.setParameterValues( { { "u", 0.5 }, { "v", 0.5 }} );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 1, 1e-12 );
}
BOOST_AUTO_TEST_CASE( test_rectangle )
{
    auto se = symbolsExpr( symbolExpr( "u", _e1 ), symbolExpr( "v", _e2 ) );
    auto a1b = expr("rectangle(t,1,2):t");
    a1b.setParameterValues( { { "t", 0 } } );
    BOOST_CHECK_SMALL( a1b.evaluate()(0,0), 1e-12 );
    a1b.setParameterValues( { { "t", 1.5 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 1, 1e-12 );
    a1b.setParameterValues( { { "t", 1 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 1, 1e-12 );
    a1b.setParameterValues( { { "t", 2 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 1, 1e-12 );
    a1b.setParameterValues( { { "t", 3} } );
    BOOST_CHECK_SMALL( a1b.evaluate()(0,0), 1e-12 );
}
BOOST_AUTO_TEST_CASE( test_triangle )
{
    auto a1b = expr("triangle(t,1,2):t");
    a1b.setParameterValues( { { "t", -2 } } );
    BOOST_CHECK_SMALL( a1b.evaluate()(0,0), 1e-12 );
    a1b.setParameterValues( { { "t", 3 } } );
    BOOST_CHECK_SMALL( a1b.evaluate()(0,0), 1e-12 );
    a1b.setParameterValues( { { "t", 0.5 } } );
    BOOST_CHECK_SMALL( a1b.evaluate()(0,0), 1e-12 );
    a1b.setParameterValues( { { "t", 1.5 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 1, 1e-12 );
    a1b.setParameterValues( { { "t", 1 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 0, 1e-12 );
    a1b.setParameterValues( { { "t", 2 } } );
    BOOST_CHECK_SMALL( a1b.evaluate()(0,0), 1e-12 );
    a1b.setParameterValues( { { "t", 1.75 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 0.5, 1e-12 );
    a1b.setParameterValues( { { "t", 1.25 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 0.5, 1e-12 );
}

BOOST_AUTO_TEST_CASE( test_mapabcd )
{
    auto a1b = expr("mapabcd(t,1,2,-1,1):t");
    a1b.setParameterValues( { { "t", -2 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), -1, 1e-12 );
    a1b.setParameterValues( { { "t", 3 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 1, 1e-12 );
    a1b.setParameterValues( { { "t", 1 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), -1, 1e-12 );
    a1b.setParameterValues( { { "t", 1.5 } } );
    BOOST_CHECK_SMALL( a1b.evaluate()(0,0), 1e-12 );
    a1b.setParameterValues( { { "t", 2 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 1, 1e-12 );
    a1b.setParameterValues( { { "t", 1.25 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), -0.5, 1e-12 );
    a1b.setParameterValues( { { "t", 1.75 } } );
    BOOST_CHECK_CLOSE( a1b.evaluate()(0,0), 0.5, 1e-12 );
}
BOOST_AUTO_TEST_SUITE_END()
