/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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

#define USE_BOOST_TEST 1
//#undef USE_BOOST_TEST
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_ginac
#include <testsuite/testsuite.hpp>
#endif

#include <iostream>
#include <string>
#include <list>

#include <boost/algorithm/string/split.hpp>
#include <boost/foreach.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
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

void runTest1()
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto XhScalar = Pch<3>( mesh );
    auto XhVectorial = Pchv<3>( mesh );

    // id scalar
    auto exprScalarFeel = Px()*Px()*cos(M_PI*Py()) + 2*Py()*sin(Px())*exp(Px()*Px());
    auto exprScalarGinac = expr( "x*x*cos(pi*y)+2*y*sin(x)*exp(x*x):x:y"/*, "scalarExpr"*/ );
    auto uFeel = XhScalar->element( exprScalarFeel );
    auto uGinac = XhScalar->element( exprScalarGinac );
    for ( size_type k=0;k<XhScalar->nLocalDof();++k )
    {
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( uFeel(k), uGinac(k), 1e-10 );
#endif
        if ( std::abs(uFeel(k)-uGinac(k) ) > 1e-10 )
            break;
    }

    // grad scalar
    auto exprScalarFeelGrad = vec( 2*Px()*cos(M_PI*Py()) + 2*Py()*( cos(Px())*exp(Px()*Px()) + sin(Px())*2*Px()*exp(Px()*Px()) ),
                                   -M_PI*Px()*Px()*sin(M_PI*Py()) + 2*sin(Px())*exp(Px()*Px()) );
    auto exprScalarGinacGrad = trans( grad<2>( exprScalarGinac/*, "scalarExprGrad"*/ ) );
    auto uFeelGrad = XhVectorial->element( exprScalarFeelGrad );
    auto uGinacGrad = XhVectorial->element( exprScalarGinacGrad );
    for ( size_type k=0;k<XhVectorial->nLocalDof();++k )
    {
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( uFeelGrad(k), uGinacGrad(k), 1e-10 );
#endif
        if ( std::abs(uFeelGrad(k)-uGinacGrad(k) ) > 1e-10 )
            break;
    }

    // laplacian scalar
    auto exprScalarFeelLaplacian = 2*exp(Px()*Px())*(4*Px()*Px()+1)*Py()*sin(Px())+8*exp(Px()*Px())*Px()*Py()*cos(Px())+(2-M_PI*M_PI*Px()*Px())*cos(M_PI*Py());
    auto exprScalarGinacLaplacian = laplacian( exprScalarGinac/*, "scalarExprLaplacian"*/ );
    auto uFeelLaplacian = XhScalar->element( exprScalarFeelLaplacian );
    auto uGinacLaplacian = XhScalar->element( exprScalarGinacLaplacian );
    for ( size_type k=0;k<XhScalar->nLocalDof();++k )
    {
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( uFeelLaplacian(k), uGinacLaplacian(k), 1e-10 );
#endif
        if ( std::abs(uFeelLaplacian(k)-uGinacLaplacian(k) ) > 1e-10 )
            break;
    }

    // id vectorial
    auto exprVectorialFeel = vec( Px()*Px()*cos(M_PI*Py()),
                                  2*Py()*sin(Px())*exp(Px()*Px()) );
    auto exprVectorialGinac = expr<2,1>( "{x*x*cos(pi*y),2*y*sin(x)*exp(x*x)}:x:y"/*, "vectorialExpr"*/ );
    auto uVectorialFeel = XhVectorial->element( exprVectorialFeel );
    auto uVectorialGinac = XhVectorial->element( exprVectorialGinac );
    for ( size_type k=0;k<XhVectorial->nLocalDof();++k )
    {
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( uVectorialFeel(k), uVectorialGinac(k), 1e-10 );
#endif
        if ( std::abs(uVectorialFeel(k)-uVectorialGinac(k) ) > 1e-10 )
            break;
    }

    // div vectorial
    auto exprVectorialFeelDiv = 2*Px()*cos(M_PI*Py()) + 2*sin(Px())*exp(Px()*Px());
    auto exprVectorialGinacDiv = div<2>( exprVectorialGinac/*, "vectorialExprDiv"*/ );
    auto uVectorialFeelDiv = XhScalar->element( exprVectorialFeelDiv );
    auto uVectorialGinacDiv = XhScalar->element( exprVectorialGinacDiv );
    for ( size_type k=0;k<XhScalar->nLocalDof();++k )
    {
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( uVectorialFeelDiv(k), uVectorialGinacDiv(k), 1e-10 );
#endif
        if ( std::abs(uVectorialFeelDiv(k)-uVectorialGinacDiv(k) ) > 1e-10 )
            break;
    }

    // id vf
    auto exprScalarField = Px()*Px()*cos(M_PI*Py()) + 2*Py()*sin(Px())*exp(Px()*Px());
    auto scalarField = XhScalar->element( exprScalarField );
    auto exprScalarFeelVF = 2*Px()*sin(Py())*exp(pow(idv(scalarField),2));
    auto exprScalarGinacVF = expr( "2*x*sin(y)*exp(u^2):x:y:u", "u", idv(scalarField) );
    auto uScalarFeelVF = XhScalar->element( exprScalarFeelVF );
    auto uScalarFeelGinacVF = XhScalar->element( exprScalarGinacVF );
    for ( size_type k=0;k<XhScalar->nLocalDof();++k )
    {
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( uScalarFeelVF(k), uScalarFeelGinacVF(k), 1e-10 );
#endif
        if ( std::abs(uScalarFeelVF(k)-uScalarFeelGinacVF(k) ) > 1e-10 )
            break;
    }

    // diff vf
    auto exprScalarFeelVFDiff = 2*Px()*sin(Py())*2*idv(scalarField)*exp(pow(idv(scalarField),2));
    auto exprScalarGinacVFDiff = diff( "2*x*sin(y)*exp(u^2):x:u:y", "u", "u", idv(scalarField) );
    auto uScalarFeelVFDiff = XhScalar->element( exprScalarFeelVFDiff );
    auto uScalarGinacVFDiff = XhScalar->element( exprScalarGinacVFDiff );
    for ( size_type k=0;k<XhScalar->nLocalDof();++k )
    {
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( uScalarFeelVFDiff(k), uScalarGinacVFDiff(k), 1e-10 );
#endif
        if ( std::abs(uScalarFeelVFDiff(k)-uScalarGinacVFDiff(k) ) > 1e-10 )
            break;
    }

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

#if defined(USE_BOOST_TEST)

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
BOOST_AUTO_TEST_SUITE_END()

#else

int main( int argc, char* argv[] )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );
    runTest0();
    runTest1();
    runTest2();
}
#endif
