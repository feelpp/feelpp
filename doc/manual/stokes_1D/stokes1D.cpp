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
//#include <feel/feelpde/preconditionerblockns.hpp>
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/product.hpp>

constexpr int dim = FEELPP_DIM;
constexpr int order_p= FEELPP_ORDER_P;

using namespace Feel;

template<typename ElementType1, typename ElementType2>
struct f_evaluate
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef value_type evaluate_type;
    typedef Feel::uint16_type uint16_type;
    static const uint16_type rank = 1;
    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;
    using element1_type = decay_type<ElementType1>;
    using element2_type = decay_type<ElementType2>;
    using node_t = boost::numeric::ublas::vector<double>;
    f_evaluate( element1_type& u1, element2_type& u2, double R ) : M_u1(u1), M_u2(u2), M_R(R), M_z(3)
    {
        M_z[0]=0;
        M_z[1]=0;
        //M_z[1]=0;
        //M_z[2]=0;
        
    }
    double operator()( uint16_type c1, uint16_type c2, boost::numeric::ublas::vector<double> const& x, boost::numeric::ublas::vector<double> const& ) const
    {
        M_z[2]=x[2];//M_z[0]=x[0];
        return M_u1(M_z)(0,0,0)*((c1==0)*x[0]/M_R+(c1==1)*x[1]/M_R)+M_u2(M_z)(0,0,0)*((c1==2)*(1-(x[0]*x[0]+x[1]*x[1])/(M_R*M_R)));
        //return M_u1(M_z)(0,0,0)*((c1==1)*x[1]/M_R+(c1==2)*x[2]/M_R)+M_u2(M_z)(0,0,0)*((c1==0)*(1-(x[1]*x[1]+x[2]*x[2])/(M_R*M_R)));
    }
    element1_type& M_u1;
    element2_type& M_u2;
    mutable node_t M_z;
    double M_R;
};


int main(int argc, char**argv )
{
    po::options_description stokes1Doptions( "Steady NS options" );
    stokes1Doptions.add_options()
    ( "rho", po::value<double>()->default_value( 1.0 ), "coeff" )
    ( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
    ( "R", po::value<double>()->default_value( 1.0 ), "Cylinder radius" )
    ( "sym", po::value<bool>()->default_value( 0 ), "symmetric formulation of the stress tensor" )
    ( "penaldir", po::value<double>()->default_value( 100 ), "coeff" )
    ( "stokes.preconditioner", po::value<std::string>()->default_value( "petsc" ), "Stokes preconditioner: petsc, PM, Blockns" )
    ( "markername", po::value<std::string>()->default_value( "centerline" ), "marker name" )
    ( "filename1d", po::value<std::string>()->default_value( "tige1D.geo" ), "1D mesh name" )
    ( "filename3d", po::value<std::string>()->default_value( "cylinder3D.geo" ), "3D mesh name" )
    ;
    stokes1Doptions.add( backend_options( "stokes" ) );
    
    Environment env( _argc=argc, _argv=argv,
                    _desc=stokes1Doptions,
                    _about=about(_name=(boost::format("stokes_1D")).str(),
                                 _author="Feel++ Consortium",
                                 _email="feelpp-devel@feelpp.org"));
    
    std::string markername = soption(_name="markername");
    std::string filename1d = soption(_name="filename1d");
    std::string filename3d = soption(_name="filename3d");
    
    auto mesh3d = loadMesh(_mesh=new Mesh<Simplex<3,1,3>>,_filename=filename3d);
    std::cout<<"loading mesh 3D: DONE \n";
    auto Xh3 = THch<order_p>( mesh3d );
    std::cout<<"creating Xh3: DONE \n";
    auto U3 = Xh3->element();
    auto u3 = U3.element<0>();
    auto p3 = U3.element<1>();
    auto V3 = Xh3->element();
    auto v3 = V3.element<0>();
    auto q3 = V3.element<1>();

    std::cout<<"Creating U 3D: DONE \n";
    if ( Environment::isMasterRank() )
    {
        std::cout << "=========================================================== " << std::endl;
        std::cout << "========================= 3D mesh ========================= " << std::endl;
        std::cout << " - mesh entities" << std::endl;
        std::cout << "      number of elements : " << mesh3d->numGlobalElements() << std::endl;
        std::cout << "         number of faces : " << mesh3d->numGlobalFaces() << std::endl;
        std::cout << "      number of points : " << mesh3d->numGlobalPoints() << std::endl;
        std::cout << "    number of vertices : " << mesh3d->numGlobalVertices() << std::endl;
        std::cout << " - mesh sizes" << std::endl;
        std::cout << "                h max : " << mesh3d->hMax() << std::endl;
        std::cout << "                h min : " << mesh3d->hMin() << std::endl;
        std::cout << "                h avg : " << mesh3d->hAverage() << std::endl;
        std::cout << "              measure : " << mesh3d->measure() << std::endl;
        std::cout << "------------------------------------------------------------ " << std::endl;
        std::cout << "FunctionSpace\tLocalDOF\tu3\tp\n";
        std::cout.width(16);
        std::cout << std::left << Xh3->nDof();
        std::cout.width(16);
        std::cout << std::left << Xh3->nLocalDof();
        std::cout.width(16);
        std::cout << std::left << Xh3->functionSpace<0>()->nDof();
        std::cout.width(16);
        std::cout << std::left << Xh3->functionSpace<1>()->nDof()<< "\n";
        std::cout<<"Number of Xh3 elements = "<< nelements(elements( mesh3d ))<<"\n";
        std::cout << "=========================================================== " << std::endl;
    }
    
    //auto mesh1d = createSubmesh(mesh3d, markededges(mesh3d,"markername"));
    auto mesh1d = loadMesh(_mesh=new Mesh<Simplex<1,1,3>>,_filename=filename1d);
    std::cout<<"creating  mesh 1D: DONE \n";
    
    
    //////////////////////////////
    
    /*typedef std::shared_ptr<Mesh<Simplex<1>>> mesh_ptrtype;
    typedef Lagrange<2, Scalar> basis_u1_type;
    typedef Lagrange<2, Scalar> basis_u2_type;
    typedef Lagrange<1, Scalar> basis_p_type;
    typedef bases<basis_u1_type,basis_u2_type,basis_p_type> basis_type;
    typedef FunctionSpace<Mesh<Simplex<1>>, basis_type> space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    mesh_ptrtype mesh1d;
    space_ptrtype Xh1 ;*/
    //////////////////////////////
    typedef FunctionSpace<Mesh<Simplex<1,1,3>>, bases<Lagrange<2, Scalar,Continuous,PointSetEquiSpaced,0>, Lagrange<2, Scalar,Continuous,PointSetEquiSpaced,1>, Lagrange<1, Scalar> > > space_type;
    
    auto Xh1 = space_type::New( mesh1d );
    std::cout<<"creating Xh1: DONE \n";
    if ( Environment::isMasterRank() )
    {
        std::cout << "=========================================================== " << std::endl;
        std::cout << "========================= 1D mesh ========================= " << std::endl;
        std::cout << " - mesh entities" << std::endl;
        std::cout << "      number of elements : " << mesh1d->numGlobalElements() << std::endl;
        std::cout << "         number of faces : " << mesh1d->numGlobalFaces() << std::endl;
        std::cout << "      number of points : " << mesh1d->numGlobalPoints() << std::endl;
        std::cout << "    number of vertices : " << mesh1d->numGlobalVertices() << std::endl;
        std::cout << " - mesh sizes" << std::endl;
        std::cout << "                h max : " << mesh1d->hMax() << std::endl;
        std::cout << "                h min : " << mesh1d->hMin() << std::endl;
        std::cout << "                h avg : " << mesh1d->hAverage() << std::endl;
        std::cout << "              measure : " << mesh1d->measure() << std::endl;
        std::cout << "------------------------------------------------------------ " << std::endl;
    }
    auto f=expr(soption("functions.f"));
    Feel::cout<< "Measure 1D: "<<integrate( _range=elements( mesh1d ), _expr=cst(1.)).evaluate()<< std::endl;
    Feel::cout<< "Int z^2 =  "<<integrate( _range=elements( mesh1d ), _expr=f).evaluate()<< std::endl;
    auto U = Xh1->element();
    auto V = Xh1->element();
    
    auto u1 = U.element<0>();
    auto u2 = U.element<1>();
    auto p = U.element<2>();
    
    u1.on(_range=elements(mesh1d), _expr=cst(0.));
    u2.on(_range=elements(mesh1d), _expr=cst(1.));
    
    auto v1 = V.element<0>();
    auto v2 = V.element<1>();
    auto q = V.element<2>();
    double mu = doption(_name="mu");
    double rho = doption(_name="rho");
    double R = doption(_name="R");

    auto e1 = exporter( _mesh=mesh1d );
    auto e3 = exporter( _mesh=mesh3d );
    //e3->save();
    
    if ( Environment::isMasterRank() )
    {
        std::cout << "FunctionSpace\tLocalDOF\tu1\t\tu2\t\tp\n";
        std::cout.width(16);
        std::cout << std::left << Xh1->nDof();
        std::cout.width(16);
        std::cout << std::left << Xh1->nLocalDof();
        std::cout.width(16);
        std::cout << std::left << Xh1->functionSpace<0>()->nDof();
        std::cout.width(16);
        std::cout << std::left << Xh1->functionSpace<1>()->nDof();
        std::cout.width(16);
        std::cout << std::left << Xh1->functionSpace<2>()->nDof() << "\n";
        std::cout<<"Number of Xh1 elements = "<< nelements(elements( mesh1d ))<<"\n";
        std::cout << "=========================================================== " << std::endl;
    }
    
    
    auto l = form1( _test=Xh1 );
    std::cout<<"creating the linear form l: DONE \n";
    auto a = form2( _trial=Xh1, _test=Xh1);
    std::cout<<"creating the bilinear form a: DONE \n";
#if 0
    a += integrate( _range=elements( mesh1d ), _expr=4*Pi*mu*inner(idt(u1),id(v1))  + 2*Pi*mu*inner(idt(u2)*id(v2)) );
    auto b=a;
    a += integrate( _range=elements( mesh1d ), _expr=-Pi*R*mu*inner(idt(u2),dz(v1)) -Pi*R*mu*inner(id(v2),dzt(u1)) );
    auto c=a;
    a += integrate( _range=elements( mesh1d ), _expr=0.5*Pi*mu*R*R*inner(dzt(u1),dz(v1)) + (2/3)*Pi*mu*R*R*inner(dzt(u2),dz(v2)));
    auto d=a;
    a +=integrate( _range=elements( mesh1d ), _expr=-2*Pi*R*idt(p)*id(v1) - (Pi*R*R/2)*idt(p)*dz(v2) );
    auto e=a;
    a +=integrate( _range=elements( mesh1d ), _expr=2*Pi*R*id(q)*idt(u1) + (Pi*R*R/2)*id(q)*dzt(u2) );
    auto g=a;
    std::cout<<"Assembling the bilinear terms: DONE \n";
#else
    a += integrate( _range=elements( mesh1d ), _expr=4*Pi*mu*inner(idt(u1),id(v1))  + 2*Pi*mu*inner(idt(u2)*id(v2)) );
    auto b=a;
    a += integrate( _range=elements( mesh1d ), _expr=-Pi*R*mu*inner(idt(u2),dx(v1)) -Pi*R*mu*inner(id(v2),dxt(u1)) );
    auto c=a;
    a += integrate( _range=elements( mesh1d ), _expr=0.5*Pi*mu*R*R*inner(dxt(u1),dx(v1)) + (2/3)*Pi*mu*R*R*inner(dxt(u2),dx(v2)));
    auto d=a;
    a +=integrate( _range=elements( mesh1d ), _expr=-2*Pi*R*idt(p)*id(v1) - (Pi*R*R/2)*idt(p)*dx(v2) );
    auto e=a;
    a +=integrate( _range=elements( mesh1d ), _expr=2*Pi*R*id(q)*idt(u1) + (Pi*R*R/2)*id(q)*dxt(u2) );
    auto g=a;
    std::cout<<"Assembling the bilinear terms: DONE \n";
#endif
    
    tic();
    
#if 0    //Dirichlet-Neumann BC
    a+=on(_range=markedfaces(mesh1d,"inlet"), _rhs=l, _element=u1, _expr=cst(0.) ) ;
    a+=on(_range=markedfaces(mesh1d,"inlet"), _rhs=l, _element=u2, _expr=cst(1.) ) ;
    std::cout<<"setting boundary condition: DONE \n";
#else    //Dirichlet-Dirichlet BC
    a+=on(_range=boundaryfaces(mesh1d), _rhs=l, _element=u1, _expr=cst(0.) ) ;
    a+=on(_range=boundaryfaces(mesh1d), _rhs=l, _element=u2, _expr=cst(1.) ) ;
#endif
    
    
    toc("Setting boundary conditions...");
    
    auto bck = backend(_prefix="stokes",_name="stokes");
    
    auto precPetscStokes = preconditioner( _prefix="stokes",
                                    _matrix=a.matrixPtr(),_pc=bck->pcEnumType(),
                                    _pcfactormatsolverpackage=bck->matSolverPackageEnumType(),
                                    _backend=bck->shared_from_this(),
                                    _worldcomm=bck->comm() );
    
    tic();
    bck->solve(_matrix=a.matrixPtr(),_solution=U,_rhs=l.vectorPtr(),_prec=precPetscStokes );
    
    toc(" - Solving Stokes...");
    
    u1.printMatlab("u1.m");
    u2.printMatlab("u2.m");
    b.matrixPtr()->printMatlab("b.m");
    c.matrixPtr()->printMatlab("c.m");
    d.matrixPtr()->printMatlab("d.m");
    e.matrixPtr()->printMatlab("e.m");
    g.matrixPtr()->printMatlab("g.m");


#if 1
#if 0
    u3.on(_range=elements(mesh3d), _expr=idf(f_evaluate<decltype(u1),decltype(u2)>( u1, u2,R)));
    e3->step(0)->add( "u3", u3 );
    u1.on(_range=elements(mesh1d), _expr=cst(0.));
    u2.on(_range=elements(mesh1d), _expr=cst(1.));
    v3.on(_range=elements(mesh3d), _expr=idf(f_evaluate<decltype(u1),decltype(u2)>( u1, u2,R)));
#endif
    e3->step(0)->add( "u3_exact", v3 );

        //e3->step(0)->add( "p", p );
    e3->save();
#else
    tic();
    e1->step(0)->add( "u1", u1 );
    e1->step(0)->add( "u2", u2 );
    e1->step(0)->add( "p", p );
    e1->save();
    toc(" - Exporting Stokes results...");
#endif
    
    Environment::saveTimers( true );



    return 0;
}
