/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 This file is part of the Feel++ library
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date     : Tue Feb 25 12:13:15 2014
 Copyright (C) 2014-2015 Feel++ Consortium
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
#include <feel/feel.hpp>
#include <feel/feelpde/boundaryconditions.hpp>
#include <feel/feelpde/preconditionerblockns.hpp>
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/product.hpp>

template<typename ElementType>
struct f_evaluate
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef value_type evaluate_type;
    typedef Feel::uint16_type uint16_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;
    f_evaluate( ElementType& u1, ElementType& u2, double R ) : M_u1(u1) M_u2(u2); M_R(R)
    {}
    double operator()( uint16_type c1, uint16_type c2, ublas::vector<double> const& x, ublas::vector<double> const& /*n*/ ) const
    {
        return M_u1(x)[0][0][0]*((c1==0)*x[0]/M_R+(c1==1)*x[1]/M_R)+M_u2(x)[0][0][0]*((c1==2)*(1-(x[0]*x[0]+x[1]*x[1])/(M_R*M_R)));
    }
    ElementType& M_u1;
    ElementType& M_u2;
    double M_R;
};


int main(int argc, char**argv )
{
    constexpr int dim = FEELPP_DIM;
    constexpr int order_p= FEELPP_ORDER_P;
    
    using namespace Feel;
    po::options_description stokes1Doptions( "Steady NS options" );
    stokes1Doptions.add_options()
    ( "rho", po::value<double>()->default_value( 1.0 ), "coeff" )
    ( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
    ( "R", po::value<double>()->default_value( 1.0 ), "Cylinder radius" )
    ( "sym", po::value<bool>()->default_value( 0 ), "symmetric formulation of the stress tensor" )
    ( "penaldir", po::value<double>()->default_value( 100 ), "coeff" )
    ( "stokes.preconditioner", po::value<std::string>()->default_value( "petsc" ), "Stokes preconditioner: petsc, PM, Blockns" )
    ;
    stokes1Doptions.add( backend_options( "stokes" ) );
    
    Environment env( _argc=argc, _argv=argv,
                    _desc=stokes1Doptions,
                    _about=about(_name=(boost::format("stokes_1D")).str(),
                                 _author="Feel++ Consortium",
                                 _email="feelpp-devel@feelpp.org"));
    
    auto mesh3d = loadMesh(_mesh=new Mesh<Simplex<3>>);
    std::cout<<"loading mesh 3D: DONE \n";
    auto Xh3 = THch<order_p>( mesh3d );
    std::cout<<"creating Xh3: DONE \n";
    auto U3 = Xh3->element();
    auto u3 = U3.element<0>();
    auto p3 = U3.element<1>();
    std::cout<<"Creating U 3D: DONE \n";
    
    auto mesh1d = createSubmesh(mesh3d, markededges(mesh3d,"centerline"));
    //auto mesh1d = loadMesh(_mesh=new Mesh<Simplex<1,1,3>>);
    std::cout<<"subtrating  mesh 1D: DONE \n";
    
    
    //////////////////////////////
    
    /*typedef boost::shared_ptr<Mesh<Simplex<1>>> mesh_ptrtype;
    typedef Lagrange<2, Scalar> basis_u1_type;
    typedef Lagrange<2, Scalar> basis_u2_type;
    typedef Lagrange<1, Scalar> basis_p_type;
    typedef bases<basis_u1_type,basis_u2_type,basis_p_type> basis_type;
    typedef FunctionSpace<Mesh<Simplex<1>>, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    mesh_ptrtype mesh1d;
    space_ptrtype Xh1 ;*/
    //////////////////////////////
    typedef FunctionSpace<Mesh<Simplex<1,1,3>>, bases<Lagrange<2, Scalar>, Lagrange<2, Scalar>, Lagrange<1, Scalar> > > space_type;
    
    auto Xh1 = space_type::New( mesh1d );
    std::cout<<"creating Xh1: DONE \n";
    auto U = Xh1->element();
    auto V = Xh1->element();
    
    auto u1 = U.element<0>();
    auto u2 = U.element<1>();
    auto p = U.element<2>();
    
    auto v1 = V.element<0>();
    auto v2 = V.element<1>();
    auto q = V.element<2>();
    std::cout<<"creating u1, v1, u2, v2, p and q: DONE \n";
    double mu = doption(_name="mu");
    double rho = doption(_name="rho");
    double R = doption(_name="R");

    auto e1 = exporter( _mesh=mesh1d );
    auto e3 = exporter( _mesh=mesh3d );
    //e3->save();
    
    std::cout<<"Measure of Xh1 = "<<nelements(_range=elements( mesh1d ))<<"\n";
    
    auto l = form1( _test=Xh1 );
    std::cout<<"creating the linear form l: DONE \n";
    auto a = form2( _trial=Xh1, _test=Xh1);
    std::cout<<"creating the bilinear form a: DONE \n";

    a += integrate( _range=elements( mesh1d ), _expr=4*Pi*mu*inner(idt(u1),id(v1))  + Pi*mu*inner(idt(u2)*id(v2)) );
        
    a +=integrate( _range=elements( mesh1d ), _expr=-2*Pi*R*idt(p)*id(v1) - (Pi*R*R/2)*idt(p)*div(v2) );
    a +=integrate( _range=elements( mesh1d ), _expr=-2*Pi*R*id(q)*idt(u1) - (Pi*R*R/2)*id(q)*divt(u2) );
    std::cout<<"Assembling the bilinear terms: DONE \n";
    
    
    tic();
    a+=on(_range=boundaryfaces(mesh1d), _rhs=l, _element=u1,
                   _expr=cst(0.) ) ;
    a+=on(_range=boundaryfaces(mesh1d), _rhs=l, _element=u2,
          _expr=cst(1.) ) ;
    std::cout<<"setting boundary condition: DONE \n";
    
    toc("Setting boundary conditions...");
    
    auto b = backend(_prefix="stokes",_name="stokes");
    
    auto precPetscStokes = preconditioner( _prefix="stokes",
                                    _matrix=a.matrixPtr(),_pc=b->pcEnumType(),
                                    _pcfactormatsolverpackage=b->matSolverPackageEnumType(),
                                    _backend=b->shared_from_this(),
                                    _worldcomm=b->comm() );
    
    tic();
    b->solve(_matrix=a.matrixPtr(),_solution=U,_rhs=l.vectorPtr(),_prec=precPetscStokes );
    
    toc(" - Solving Stokes...");
    
    
#if 0
    u3.on(_range=elements(mesh3d), _expr=idf(f_evaluate(u1,u2,R)));
#endif
    /*tic();
    e1->step(0)->add( "u1", u1 );
    e1->step(0)->add( "u2", u2 );
    e1->step(0)->add( "p", p );
    e1->save();
    toc(" - Exporting Stokes results...");*/
    Environment::saveTimers( true );



    return 0;
}
