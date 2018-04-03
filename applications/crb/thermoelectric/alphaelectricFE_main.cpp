//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>

#include <boost/math/special_functions/binomial.hpp>

int main( int argc, char** argv)
{
    using boost::math::binomial_coefficient;
    using namespace Feel;

    po::options_description opt("options");
    opt.add_options()
        ( "zb", po::value<double>()->default_value(0.5), "z of the base" )
        ( "zt", po::value<double>()->default_value(4.5), "z of the top" )
        ( "p", po::value<std::vector<double> >()->multitoken(), "control points" )
        ( "sigma", po::value<double>()->default_value(58000), "sigma" )
        ( "dif", po::value<double>()->default_value(3333), "dif" )
        ;
    Environment env( _argc=argc, _argv=argv,
                     _desc=opt.add( backend_options("V")));

    auto mesh = loadMesh( new Mesh<Simplex<3> >, _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    auto meshCond = createSubmesh( mesh, markedelements(mesh, "cond"));
    auto meshHelix = createSubmesh( mesh, markedelements(mesh, {"helix1","helix2"}));

    /*********************** ALPHA ***********************/
    double zb = doption("zb");
    double zt = doption("zt");
    std::vector<double> p = vdoption("p");
    p.insert(p.begin(), 0);
    p.push_back(0);
    int n = p.size() - 1;

    std::stringstream alphaStream;
    alphaStream << "0 + (z<" << zt << ")*(z>" << zb << ")*(";
    for( int k = 1; k < n; ++k )
    {
        double ckn = binomial_coefficient<double>(n,k);
        alphaStream << " + " << ckn
                    << "*((z-" << zb << ")/" << zt-zb << ")^" << k
                    << "*(1-(z-" << zb << ")/" << zt-zb << ")^" << n-k
                    << "*" << p[k];
    }
    alphaStream << "):x:y:z";
    std::string alphaStr = alphaStream.str();
    auto alphaExpr = expr(alphaStr);

    std::stringstream alphaPrimeStream;
    alphaPrimeStream << "0 + (z<" << zt << ")*(z>" << zb << ")*(";
    alphaPrimeStream << n << "/" << zt-zb << "*(";
    for( int k = 0; k < n; ++k )
    {
        double ckn1 = binomial_coefficient<double>(n-1,k);
        alphaPrimeStream << " + " << (p[k+1]-p[k]) << "*" << ckn1
                         << "*((z-" << zb << ")/" << zt-zb << ")^" << k
                         << "*(1-(z-" << zb << ")/" << zt-zb << ")^" << n-1-k;
    }
    alphaPrimeStream << ")):x:y:z";
    std::string alphaPrimeStr = alphaPrimeStream.str();
    auto alphaPrimeExpr = expr(alphaPrimeStr);

    Feel::cout << "alpha(z) : " << alphaStr << std::endl
               << "alpha\'(z): " << alphaPrimeStr << std::endl;

    auto Xh = Pch<1>(meshCond);
    auto alpha = Xh->element(alphaExpr);

    auto e = exporter(_mesh=meshCond, _name="conductor");
    e->add("alpha", alpha);

    auto Vh = Pch<1>(meshHelix);
    auto alphaH = Vh->element(alphaExpr);

    auto eH = exporter(_mesh=meshHelix, _name="helix");
    eH->add("alpha", alphaH);
    eH->save();

    /********************** ELECTRO **********************/
    auto V = Xh->element();
    auto phiV = Xh->element();
    double gamma = doption("parameters.gamma");
    double sigma = doption("sigma");
    double current = doption("dif");

    // auto psi = vec( cos(alphaExpr)*Px() + sin(alphaExpr)*Py(),
    //                 -sin(alphaExpr)*Px() + cos(alphaExpr)*Py(),
    //                 Pz() );
    // auto J = mat<3,3>( cos(alphaExpr), sin(alphaExpr), alphaPrimeExpr*(-sin(alphaExpr)*Px()+cos(alphaExpr)*Py()),
    //                    -sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*(-cos(alphaExpr)*Px()-sin(alphaExpr)*Py()),
    //                    cst(0.), cst(0.), cst(1.) );
    // // auto J = grad<3,3>(psi);
    // auto Jac = abs(det(J));
    // Feel::cout << "Jacobien : " << Jac << std::endl;
    auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                          sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                          cst(0.), cst(0.), cst(1.) );
    // auto Jinv = inv(J);
    auto gradVJ = gradt(V)*Jinv;
    auto gradPhiVJ = grad(phiV)*Jinv;

    tic();
    auto aV = form2(_test=Xh, _trial=Xh);
    aV = integrate( elements(meshCond), inner(gradVJ,gradPhiVJ) );
    aV += integrate( markedfaces(meshCond, "base"),
                     gamma/hFace()*inner(idt(V),id(phiV))
                     -inner(gradVJ*N(),id(phiV))
                     -inner(gradPhiVJ*N(),idt(V)) );
    aV += integrate( markedfaces(meshCond, "top"),
                     gamma/hFace()*inner(idt(V),id(phiV))
                     -inner(gradVJ*N(),id(phiV))
                     -inner(gradPhiVJ*N(),idt(V)) );

    auto fV = form1(_test=Xh);
    fV = integrate( markedfaces(meshCond, "top"),
                    current*(gamma/hFace()*id(phiV) - gradPhiVJ*N()) );
    aV.solve(_rhs=fV, _solution=V, _name="V");
    toc("solve V");

    auto V0 = Xh->element();
    auto a0 = form2(_test=Xh, _trial=Xh);
    a0 = integrate( elements(meshCond), inner(gradt(V0),grad(phiV)) );
    a0 += integrate( markedfaces(meshCond, "base"),
                     gamma/hFace()*inner(idt(V),id(phiV))
                     -inner(gradt(V0)*N(),id(phiV))
                     -inner(grad(phiV)*N(),idt(V0)) );
    auto f0 = form1(_test=Xh);
    f0 = integrate( markedfaces(meshCond, "top"),
                    current*(gamma/hFace()*id(phiV) - grad(phiV)*N()) );
    a0.solve(_rhs=f0, _solution=V0, _name="V");

    e->add("V",V);
    e->add("V0",V0);
    e->save();

    auto am = aV.matrixPtr();
    auto fv = fV.vectorPtr();
    auto a0m = a0.matrixPtr();
    auto f0v = f0.vectorPtr();
    double amN = am->l1Norm();
    double fvN = fv->l2Norm();
    double a0N = a0m->l1Norm();
    double f0N = f0v->l2Norm();
    am->addMatrix(-1., a0m);
    fv->add(-1., f0v);
    double aN = am->l1Norm();
    double fN = fv->l2Norm();
    Feel:: cout << "||A|| = " << amN << "\t||F|| = " << fvN << std::endl;
    Feel:: cout << "||A0|| = " << a0N << "\t||F0|| = " << f0N << std::endl;
    Feel:: cout << "||errA|| = " << aN << "\t||errF|| = " << fN << std::endl;

    return 0;
}
