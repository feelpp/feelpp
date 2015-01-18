/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-07-05

  Copyright (C) 2013 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_integration_ginac.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-07-05
 */
#include <sstream>
#include <boost/timer.hpp>

//#define BOOST_TEST_MODULE ginac integration testsuite
//#include <testsuite/testsuite.hpp>

#include <boost/mpl/list.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

#if 0
FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( ginacsuite )
//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
typedef boost::mpl::list<boost::mpl::int_<2> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( ginacint, T, dim_types )
{
    BOOST_TEST_MESSAGE( "check 1 ginac int  = " << T::value << "\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "hypercube-%1%" )  % T::value ).str() ,
                                                _usenames=true,
                                                _addmidpoint=false,
                                                _shape="hypercube",
                                                _h=0.2 ),
                                        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    auto P1h = Pch<1>( mesh );
    //auto vars = symbols<T::value>();
    symbol x("x"), y ("y");

    //auto f1g = vars[0];auto f2g = vars[1];
    ex f1g = sin(x)/y ;
    ex f2g = x*cos(y);
    std::vector<symbol> s = {x,y};
    auto f1= expr( f1g, s, "a" );
    auto f2 = expr(f2g, s, "b");
    auto xg = integrate(_range=elements(mesh), _expr=cst(2.)*f1/2. ).evaluate()(0,0);
    auto yg = integrate(_range=elements(mesh), _expr=cst(2.)*f2/2. ).evaluate()(0,0);

    BOOST_CHECK_CLOSE( xg, 0.5, 1e-10 );
    BOOST_CHECK_CLOSE( yg, 0.5, 1e-10 );
    BOOST_TEST_MESSAGE( "test gravity center " << xg << "," << yg );

}

BOOST_AUTO_TEST_SUITE_END()

#if 1
int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}

#endif
#else
#if 0
void ginacint()
{
    using namespace Feel;
    typedef Mesh<Simplex<2,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "hypercube-%1%" )  % 2 ).str() ,
                                                      _usenames=true,
                                                      _addmidpoint=false,
                                                      _shape="hypercube",
                                                      _h=0.2,
                                                      _xmin=0., _ymin=0., _zmin=0. ),
                                        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    auto P1h = Pch<1>( mesh );
    auto vars = symbols<2>();
    //symbol x("x"), y ("y");

    //auto f1g = vars[0];auto f2g = vars[1];
    //ex f1g = x;
    //ex f2g = y;
    //std::vector<symbol> s = {x,y};
    auto f1g = parse( "x", vars );
    auto f2g = parse( "y", vars );
    auto f1= expr( f1g, vars, "a" );
    auto f2 = expr("y", vars, "b");
    auto xg = integrate(_range=elements(mesh), _expr=cst(2.)*f1/2. ).evaluate()(0,0);
    auto yg = integrate(_range=elements(mesh), _expr=cst(2.)*f2/2. ).evaluate()(0,0);

    CHECK( math::abs(xg- 0.5)< 1e-10 ) << "check failed : xg = " << xg;
    CHECK( math::abs(yg- 0.5)< 1e-10 ) << "check failed : yg = " << yg;
    LOG(INFO) << "test done\n";

    matrix u_exact_g = matrix(2,1);
    u_exact_g = f1g,f2g;
    auto u_exact = expr<2,1,2>( u_exact_g, vars, "x_exact" );
    auto Xg = integrate(_range=elements(mesh), _expr=u_exact ).evaluate();
    LOG(INFO) << "Xg = " << Xg;
    CHECK( math::abs(Xg(0,0) - 0.5)< 1e-10 ) << "check failed : xg = " << Xg(0,0);
    CHECK( math::abs(Xg(1,0) - 0.5)< 1e-10 ) << "check failed : yg = " << Xg(1,0);
}
#endif
void ginacinthypercube()
{
    using namespace Feel;
    typedef Mesh<Hypercube<2,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "hypercube-structured-%1%" )  % 2 ).str() ,
                                                      _usenames=true,
                                                      _addmidpoint=false,
                                                      _structured=1,
                                                      _convex="hypercube",
                                                      _shape="hypercube",
                                                      _h=option(_name="gmsh.hsize").as<double>(),
                                                      _xmin=-1, _ymin=-1, _zmin=-1. ),
                                        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    auto P1h = Pch<1>( mesh );
    auto vars = symbols<2>();
    auto f1= expr( "x*y", vars, "a" );
    auto f2 = expr( "x*(y+1)", vars, "b");
    auto u = P1h->element();
    auto v = P1h->element();
    auto w = P1h->element();
    auto us = P1h->element();
    u = project( _space=P1h, _range=elements(mesh), _expr=Px()*Py() );
    v = project( _space=P1h, _range=elements(mesh), _expr=Px()*(Py()+1) );
    w = project( _space=P1h, _range=elements(mesh), _expr=Px()+Py() );

    double v1 = integrate(_range=elements(mesh), _expr=f1 ).evaluate()(0,0);
    double v11 = integrate(_range=elements(mesh), _expr=Px()*Py() ).evaluate()(0,0);
    double v111 = integrate(_range=elements(mesh), _expr=idv(u) ).evaluate()(0,0);
    double v2 = integrate(_range=elements(mesh), _expr=f2 ).evaluate()(0,0);
    double v22 = integrate(_range=elements(mesh), _expr=Px()*(Py()+1) ).evaluate()(0,0);
    double v222 = integrate(_range=elements(mesh), _expr=idv(v) ).evaluate()(0,0);
    std::cout << "int x*y=" << v1 << ", " << v11 << ", " << v111 << " error: " << math::abs(v1-v11) << "\n";
    std::cout << "int x*(y+1=" << v2 << ", " << v22 << " error: " << math::abs(v2-v22) << "\n";

    auto gf1= expr<2,1,1>( GiNaC::grad("x*y",vars), vars, "ga" );
    auto gf2 = expr<2,1,1>( GiNaC::grad("x*(y+1)",vars), vars, "gb");
    auto gv1 = integrate(_range=elements(mesh), _expr=gf1 ).evaluate();
    auto gv11 = integrate(_range=elements(mesh), _expr=vec(Py(),Px() ) ).evaluate();
    auto gv111 = integrate(_range=elements(mesh), _expr=gradv(u) ).evaluate();
    auto gv2 = integrate(_range=elements(mesh), _expr=gf2 ).evaluate();
    auto gv22 = integrate(_range=elements(mesh), _expr=vec(Py()+1,Px()) ).evaluate();
    auto gv222 = integrate(_range=elements(mesh), _expr=gradv(v) ).evaluate();
    std::cout << "int nabla x*y=" << gv1 << ", " << gv11 << ",  " << gv111 << " error: " << (gv1-gv11).norm() << "\n";
    std::cout << "int nabla x*(y+1=" << gv2 << ", " << gv22 << ", " << gv222 << " error: " << (gv2-gv22).norm() << "\n";



    auto l1 = form1(_test=P1h);
    l1 = integrate(_range=elements(mesh), _expr=grad(v)*ones<2,1>() );
    std::cout << "l1 u = " << l1(u) << "\n";
    std::cout << "l1 v = " << l1(v) << "\n";

    auto l3 = form1(_test=P1h);
    l3 = integrate(_range=elements(mesh), _expr=idv(u)*id(u) + gradv(u)*trans(grad(v)));
    auto l2 = form2(_trial=P1h, _test=P1h);
    l2 = integrate(_range=elements(mesh), _expr=
#if 0
                   option(_name="parameters.alpha").as<double>()*print(print(dxt(u),"dxt=")*print(dx(v),"dx(v)"),"dxudxv=")+
                   + 0.*print(Px(),"x=")+ 0.*print(Py(),"y=")+
                   + 0.*print(JinvT()(0,0),"B(0,0)=")+ 0.*print(JinvT()(0,1),"B(0,1)=")
                   + 0.*print(JinvT()(1,0),"B(1,0)=")+ 0.*print(JinvT()(1,1),"B(1,1)=")+
#else
#if 1
                   option(_name="parameters.alpha").as<double>()*dxt(u)*dx(v)+0.*Py()+
#else
                   option(_name="parameters.alpha").as<double>()*dxt(u)*dx(v)+ 0*JinvT()(0,0)+0.*Px()+
#endif
#endif
                   option(_name="parameters.beta").as<double>()*dyt(u)*dy(v) );
    l2 += integrate(_range=elements(mesh), _expr=idt(u)*id(v));
    l2.solve( _solution=us, _rhs=l3 );

    l2.matrixPtr()->printMatlab("A.m");
    auto e1= dxt(u)*dx(u);
    auto e2= dxt(u)*dx(u)+0.*Px();
    std::cout << "e1.context = " << e1.context << "\n";
    std::cout << "e2.context = " << e2.context << "\n";
    std::cout << "||u-us||_2 = " << normL2( _range=elements(mesh), _expr=idv(u)-idv(us) ) << "\n";

    std::cout << "l2 u 1 = " << l2(u,w) << "\n";
    std::cout << "l2 1 u  = " << l2(w,u) << "\n";
    std::cout << "l2 v 1 = " << l2(v,w) << "\n";
    std::cout << "l2 1 v  = " << l2(w,v) << "\n";


    if ( math::abs( option(_name="gmsh.hsize" ).as<double>()-2 ) < 1e-12 )
    {
        auto l4 = form2(_trial=P1h, _test=P1h);
        l4 = integrate(_range=elements(mesh), _expr=dxt(u)*dx(v));
        l4.matrixPtr()->printMatlab("dxdx.m");
        CHECK(P1h->nLocalDof() == 4 ) << "Invalid number of dof :  " << P1h->nLocalDof() << " should be 4";
        // (1-xi)*(1-eta)/4
        u.zero();
        u[0] = 1;
        double nl2 = normL2( _range=elements(mesh), _expr=idv(u)-(1-Px())*(1-Py())/4.);
        std::cout << "(1-xi)*(1-eta)/4 = " << nl2 << "\n";
        double nl2x = normL2( _range=elements(mesh), _expr=dxv(u)+(1-Py())/4);
        double nl2y = normL2( _range=elements(mesh), _expr=dyv(u)+(1-Px())/4);
        std::cout << "dx (1-xi)*(1-eta)/4 = " << nl2x << "\n";
        std::cout << "dy (1-xi)*(1-eta)/4 = " << nl2y << "\n";
        l3 = integrate(_range=elements(mesh), _expr=id(u)*(1-Px())*(1-Py())/4. + gradv(u)*trans(grad(v)));
        l2.solve( _solution=us, _rhs=l3 );
        std::cout << "us=" << us << "\n";
        std::cout << "||.||_2 = " << normL2( _range=elements(mesh), _expr=(1-Px())*(1-Py())/4.-idv(us) ) << "\n";
        std::cout << "dxdx  = " << integrate(_range=elements(mesh), _expr=dxv(u)*dxv(u)).evaluate()(0,0) << "\n";
        // (1+xi)*(1-eta)/4
        u.zero();
        u[1] = 1;
        nl2 = normL2( _range=elements(mesh), _expr=idv(u)-(1+Px())*(1-Py())/4.);
        std::cout << "(1+xi)*(1-eta)/4 = " << nl2 << "\n";
        nl2x = normL2( _range=elements(mesh), _expr=dxv(u)-(1-Py())/4);
        nl2y = normL2( _range=elements(mesh), _expr=dyv(u)+(1+Px())/4);
        std::cout << "dx (1+xi)*(1-eta)/4 = " << nl2x << "\n";
        std::cout << "dy (1+xi)*(1-eta)/4 = " << nl2y << "\n";
        l3 = integrate(_range=elements(mesh), _expr=id(u)*(1+Px())*(1-Py())/4. + gradv(u)*trans(grad(v)));
        l2.solve( _solution=us, _rhs=l3 );
        std::cout << "us=" << us << "\n";
        std::cout << "||.||_2 = " << normL2( _range=elements(mesh), _expr=(1+Px())*(1-Py())/4.-idv(us) ) << "\n";
        std::cout << "dxdx  = " << integrate(_range=elements(mesh), _expr=dxv(u)*dxv(u)).evaluate()(0,0) << "\n";
        // (1+xi)*(1+eta)/4
        u.zero();
        u[2] = 1;
        nl2 = normL2( _range=elements(mesh), _expr=idv(u)-(1+Px())*(1+Py())/4.);
        std::cout << "(1+xi)*(1+eta)/4 = " << nl2 << "\n";
        nl2x = normL2( _range=elements(mesh), _expr=dxv(u)-(1+Py())/4);
        nl2y = normL2( _range=elements(mesh), _expr=dyv(u)-(1+Px())/4);
        std::cout << "dx (1+xi)*(1+eta)/4 = " << nl2x << "\n";
        std::cout << "dy (1+xi)*(1+eta)/4 = " << nl2y << "\n";
        l3 = integrate(_range=elements(mesh), _expr=id(u)*(1+Px())*(1+Py())/4. + gradv(u)*trans(grad(v)));
        l2.solve( _solution=us, _rhs=l3 );
        std::cout << "us=" << us << "\n";
        std::cout << "||.||_2 = " << normL2( _range=elements(mesh), _expr=(1+Px())*(1+Py())/4.-idv(us) ) << "\n";
        std::cout << "dxdx  = " << integrate(_range=elements(mesh), _expr=dxv(u)*dxv(u)).evaluate()(0,0) << "\n";
        // (1-xi)*(1+eta)/4
        u.zero();
        u[3] = 1;
        nl2 = normL2( _range=elements(mesh), _expr=idv(u)-(1-Px())*(1+Py())/4.);
        std::cout << "(1-xi)*(1+eta)/4 = " << nl2 << "\n";
        nl2x = normL2( _range=elements(mesh), _expr=dxv(u)+(1+Py())/4);
        nl2y = normL2( _range=elements(mesh), _expr=dyv(u)-(1-Px())/4);
        std::cout << "dx (1-xi)*(1+eta)/4 = " << nl2x << "\n";
        std::cout << "dy (1-xi)*(1+eta)/4 = " << nl2y << "\n";
        l3 = integrate(_range=elements(mesh), _expr=id(u)*(1-Px())*(1+Py())/4. + gradv(u)*trans(grad(v)));
        l2.solve( _solution=us, _rhs=l3 );
        std::cout << "us=" << us << "\n";
        std::cout << "||.||_2 = " << normL2( _range=elements(mesh), _expr=(1-Px())*(1+Py())/4.-idv(us) ) << "\n";
        std::cout << "dxdx  = " << integrate(_range=elements(mesh), _expr=dxv(u)*dxv(u)).evaluate()(0,0) << "\n";

        Eigen::MatrixXd A(4,4);
        for(int i = 0; i < 4; i ++ )
        {
            u.zero();
            u[i]=1;
            for(int j = 0; j < 4; j ++ )
            {
                v.zero();
                v[j]=1;
                A(i,j) = integrate(_range=elements(mesh), _expr=dxv(u)*dxv(v)).evaluate()(0,0);
            }
        }
        std::cout << "A=" << A << "\n";
    }
}


#if 0

void poiseuille()
{
    using namespace Feel;
    typedef Mesh<Simplex<2,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "hypercube-%1%" )  % 2 ).str() ,
                                                      _usenames=true,
                                                      _addmidpoint=false,
                                                      _shape="hypercube",
                                                      _h=0.2,
                                                      _xmin=0., _ymin=0., _zmin=0. ),
                                        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    auto P1h = Pch<1>( mesh );
    auto vars = symbols<2>();
    symbol x("x"), y ("y");
    ex u_exact_x=(1-y*y);
    ex u_exact_y=0;
    ex p_exact_exp=-2*x+5;

#if 0
    std::string dim_str =  boost::str( boost::format( "2D" ) );
    std::string u1_str = soption(_name="u_exact_x",_prefix=dim_str);
    std::string u2_str = soption(_name="u_exact_y",_prefix=dim_str);
    std::string p_str = soption(_name="p_exact_ex",_prefix=dim_str);
    LOG(INFO) << "ux = " << u1_str;
    LOG(INFO) << "uy = " << u2_str;
    LOG(INFO) << "p = " << p_str;
#endif

    auto u1 = parse( "1-y*y" , vars );
    auto u2 = parse( "0", vars );
    auto p_exact_g = parse( "-2*x+5", vars );
    matrix u_exact_g = matrix(2,1);
    u_exact_g = u1,u2;
    auto gradu_exact_g = grad( u_exact_g, vars );
    auto divu_exact_g = div( u_exact_g, vars );
    std::cout << "u_exact = " << u_exact_g;
    LOG(INFO) << "p_exact = " << p_exact_g;
    LOG(INFO) << "gradu_exact_g = " << gradu_exact_g;
    LOG(INFO) << "divu_exact_g = " << divu_exact_g;

    auto u_exact = expr<2,1,7>( u_exact_g, vars, "u_exact" );
    auto u1_exact = expr<7>( u1, vars, "u1_exact" );
    auto u2_exact = expr<7>( u2, vars, "u2_exact" );
    auto p_exact = expr<7>( p_exact_g, vars, "p_exact" );
    auto gradu_exact = expr<2,2,7>( gradu_exact_g, vars, "gradu_exact" );
    auto divu_exact = expr<1,1,7>( divu_exact_g, vars, "divu_exact" );

    auto u11=1-Py()*Py();
    auto u22= cst(0.);
    //matrix u_exact_strong = matrix(2,1);
    //u_exact_strong = u11,u22;
    auto u_exact_strong = vec( u11,u22 ) ;
    auto p_exact_strong = -2*Px()+5;

    auto du_dx = cst(0.);
    auto du_dy = -2*Py();
    auto dv_dx =  cst(0.);
    auto dv_dy =  cst(0.);

    auto gradu_exact_strong =  mat<2,2>( du_dx, du_dy, dv_dx, dv_dy) ;
    auto divu_exact_strong = du_dx + dv_dy;

    double u_errorL2 = normL2( _range=elements( mesh ), _expr=( print(print(u_exact,"gu") - print(u_exact_strong,"fu"),"error") ) );
    LOG(INFO) << "u mean g = " << mean( _range=elements( mesh ), _expr=u_exact );
    LOG(INFO) << "u mean f = " << mean( _range=elements( mesh ), _expr=u_exact_strong );
    double u1_errorL2 = normL2( _range=elements( mesh ), _expr=( u11 - u1_exact ) );
    double u2_errorL2 = normL2( _range=elements( mesh ), _expr=( u22 - u2_exact ) );
    double p_errorL2 = normL2( _range=elements( mesh ), _expr=( p_exact-p_exact_strong ) );
    double gradu_errorL2 = normL2( _range=elements( mesh ), _expr=( gradu_exact-gradu_exact_strong ) );
    double divu_errorL2 = normL2( _range=elements( mesh ), _expr=( divu_exact-divu_exact_strong ) );


    LOG(INFO) <<"u-error-L2 = "<<u_errorL2<< "\n" ;
    LOG(INFO) <<"u1-error-L2 = "<<u1_errorL2<< "\n" ;
    LOG(INFO) <<"u2-error-L2 = "<<u2_errorL2<< "\n" ;
    LOG(INFO) <<"p-error-L2 = "<<p_errorL2<< "\n" ;
    LOG(INFO) <<"gradu-error-L2 = "<<gradu_errorL2<< "\n" ;
    LOG(INFO) <<"divu-error-L2 = "<<divu_errorL2<< "\n" ;

}
#endif

int main(int argc, char**argv )
{
    using namespace Feel;
    Environment env( Feel::_argc=argc,
                     Feel::_argv=argv );
    ginacinthypercube();
    //poiseuille();
}

#endif
