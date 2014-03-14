/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
             Abdoulaye Samake <abdoulaye.samake@imag.fr>
       Date: 2012-04-25

  Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

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
   \file mortar.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \author Abdoulaye Samake <abdoulaye.samake@imag.fr>
   \date 2012-04-25
 */
#include <feel/feel.hpp>
#include <feel/feeltiming/tic.hpp>

/** use Feel namespace */
using namespace Feel;

/**
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description mortaroptions("Mortar options");
    mortaroptions.add_options()
        ("hsize1", po::value<double>()->default_value( 0.1 ), "mesh size for first domain")
        ("hsize2", po::value<double>()->default_value( 0.1 ), "mesh size for second domain")
        ("shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
        ("coeff", po::value<double>()->default_value( 1 ), "grad.grad coefficient")
        ("weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ("penaldir", Feel::po::value<double>()->default_value( 10 ),
         "penalisation parameter for the weak boundary Dirichlet formulation")
        ;
    return mortaroptions.add( Feel::feel_options() );
}

/**
 * \class MortarProd
 *
 * MortarProd Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=2 or 3)
 */
template<int Dim, int Order1, int Order2>
class MortarProd
    :
    public Simget
{
    typedef Simget super;
public:

    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef Mesh< Simplex<Dim,1,Dim>,value_type,0 > sub1_mesh_type;
    typedef Mesh< Simplex<Dim,1,Dim>,value_type,1 > sub2_mesh_type;
    typedef typename sub1_mesh_type::trace_mesh_type trace_mesh_type;
    typedef meshes<sub1_mesh_type,sub2_mesh_type,trace_mesh_type> mesh_type;
    typedef bases<Lagrange<Order1,Scalar>,Lagrange<Order2,Scalar>,Lagrange<Order1,Scalar> > basis_type;
    typedef FunctionSpace< mesh_type, basis_type, mortars<NoMortar,NoMortar,Mortar> > space_type;
    typedef typename space_type::element_type element_type;
    typedef Exporter<sub1_mesh_type> export1_type;
    typedef Exporter<sub2_mesh_type> export2_type;
    typedef boost::shared_ptr<export1_type> export1_ptrtype;
    typedef boost::shared_ptr<export2_type> export2_ptrtype;
    typedef Exporter<trace_mesh_type> trace_export_type;
    typedef boost::shared_ptr<trace_export_type> trace_export_ptrtype;

    /**
     * Constructor
     */
    MortarProd()
        :
        super(),
        M_backend( backend_type::build( this->vm() ) ),
        mesh1Size( this->vm()["hsize1"].template as<double>() ),
        mesh2Size( this->vm()["hsize2"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() ),
        exporter1( export1_type::New( this->vm(),
                                      (boost::format( "%1%-%2%-%3%" )
                                       % this->about().appName()
                                       % Dim % int(1) ).str() ) ),
        exporter2( export2_type::New( this->vm(),
                                      (boost::format( "%1%-%2%-%3%" )
                                       % this->about().appName()
                                       % Dim % int(2) ).str() ) ),
        trace_exporter( trace_export_type::New( this->vm(),
                                                (boost::format( "%1%-%2%-%3%" )
                                                 % this->about().appName()
                                                 % Dim % int(3) ).str() ) )
    {}

    void exportResults();
    void computeErrors();
    void run();

private:

    backend_ptrtype M_backend;
    double mesh1Size;
    double mesh2Size;
    std::string shape;
    export1_ptrtype exporter1;
    export2_ptrtype exporter2;
    trace_export_ptrtype trace_exporter;
    std::vector<int> outside1;
    std::vector<int> outside2;
    int gamma1;
    int gamma2;
    element_type u;

}; // MortarProd


template<int Dim, int Order1, int Order2>
void
MortarProd<Dim, Order1, Order2>::exportResults()
{
    auto Xh=u.functionSpace();

    double pi = M_PI;
    auto g = sin(pi*Px())*cos(pi*Py())*cos(pi*Pz());

    auto e1 = vf::project( Xh->template functionSpace<0>(), elements(Xh->template mesh<0>()), g );
    auto e2 = vf::project( Xh->template functionSpace<1>(), elements(Xh->template mesh<1>()), g );

    exporter1->step(0)->setMesh( Xh->template mesh<0>() );
    exporter1->step(0)->add( "solution", (boost::format( "solution-%1%" ) % int(1) ).str(), u.template element<0>() );
    exporter1->step(0)->add( "exact", (boost::format( "exact-%1%" ) % int(1) ).str(), e1 );
    exporter1->save();

    exporter2->step(0)->setMesh( Xh->template mesh<1>() );
    exporter2->step(0)->add( "solution",(boost::format( "solution-%1%" ) % int(2) ).str(), u.template element<1>() );
    exporter2->step(0)->add( "exact",(boost::format( "exact-%1%" ) % int(2) ).str(), e2 );
    exporter2->save();

    trace_exporter->step(0)->setMesh( Xh->template mesh<2>() );
    trace_exporter->step(0)->add( "lambda",(boost::format( "lambda-%1%" ) % int(3) ).str(), u.template element<2>() );
    trace_exporter->save();

    std::ofstream ofs( (boost::format( "%1%.sos" ) % this->about().appName() ).str().c_str() );

    if ( ofs )
    {
        ofs << "FORMAT:\n"
            << "type: master_server gold\n"
            << "SERVERS\n"
            << "number of servers: " << int(2) << "\n";
        for( int j = 1; j <= 2; ++ j )
        {
            ofs << "#Server " << j << "\n";
            ofs << "machine id: " << mpi::environment::processor_name()  << "\n";
            ofs << "executable:\n";
            ofs << "data_path: .\n";
            ofs << "casefile: space_prod-" << Dim << "-" << j << "-1_0.case\n";
        }
    }

} // MortarProd::export

template<int Dim, int Order1, int Order2>
void
MortarProd<Dim, Order1, Order2>::computeErrors()
{
    auto Xh = u.functionSpace();
    value_type pi = M_PI;
    auto g = sin(pi*Px())*cos(pi*Py())*cos(pi*Pz());
    auto gradg = trans( +pi*cos(pi*Px())*cos(pi*Py())*cos(pi*Pz())*unitX()
                        -pi*sin(pi*Px())*sin(pi*Py())*cos(pi*Pz())*unitY()
                        -pi*sin(pi*Px())*cos(pi*Py())*sin(pi*Pz())*unitZ() );

    double L2error12 =integrate(elements(Xh->template mesh<0>()),(idv(u.template element<0>())-g)*(idv(u.template element<0>())-g) ).evaluate()(0,0);
    double L2error22 =integrate(elements(Xh->template mesh<1>()),(idv(u.template element<1>())-g)*(idv(u.template element<1>())-g) ).evaluate()(0,0);

    double semi_H1error1 =integrate(elements(Xh->template mesh<0>()),
                                    ( gradv(u.template element<0>())-gradg )*trans( (gradv(u.template element<0>())-gradg) ) ).evaluate()(0,0);

    double semi_H1error2 =integrate(elements(Xh->template mesh<1>()),
                                    ( gradv(u.template element<1>())-gradg )*trans( (gradv(u.template element<1>())-gradg) ) ).evaluate()(0,0);

    double error =integrate(elements(Xh->template mesh<2>()), (idv(u.template element<0>())-idv(u.template element<1>()))*
                            (idv(u.template element<0>())-idv(u.template element<1>())) ).evaluate()(0,0);

    double global_error = math::sqrt(L2error12 + L2error22 + semi_H1error1 + semi_H1error2);

    std::cout << "----------L2 errors---------- \n" ;
    std::cout << "||u1_error||_L2=" << math::sqrt(L2error12) << "\n";
    std::cout << "||u2_error||_L2=" << math::sqrt(L2error22) << "\n";
    std::cout << "----------H1 errors---------- \n" ;
    std::cout << "||u1_error||_H1=" << math::sqrt( L2error12 + semi_H1error1 ) << "\n";
    std::cout << "||u2_error||_H1=" << math::sqrt( L2error22 + semi_H1error2 ) << "\n";
    std::cout << "||u_error||_H1=" << global_error << "\n";
    std::cout << "L2 norm of jump at interface  \n" ;
    std::cout << "||u1-u2||_L2=" << math::sqrt(error) << "\n";
    std::cout << "----------------------------- \n" ;
}

template<int Dim, int Order1, int Order2>
void
MortarProd<Dim, Order1, Order2>::run()
{

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "research/hamm/%1%/%2%-%3%/P%4%-P%5%/h_%6%-%7%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order1
                                       % Order2
                                       % mesh1Size
                                       % mesh2Size );

    std::cout<<"create meshes starts\n";
    tic();

    auto mesh1 = createGMSHMesh( _mesh=new sub1_mesh_type,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                 _desc=domain( _name=(boost::format( "%1%-%2%-%3%" ) % shape % Dim % int(1)).str() ,
                                               _addmidpoint=false, _usenames=false, _shape=this->shape, _dim=Dim,
                                               _h=mesh1Size, _xmin=0., _xmax=0.5, _ymin=0., _ymax=1.,
                                               _zmin=0., _zmax=1.) );

    auto mesh2 = createGMSHMesh( _mesh=new sub2_mesh_type,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                 _desc=domain( _name=(boost::format( "%1%-%2%-%3%" ) % shape % Dim % int(2)).str() ,
                                               _addmidpoint=false, _usenames=false, _shape=this->shape, _dim=Dim,
                                               _h=mesh2Size, _xmin=0.5, _xmax=1., _ymin=0., _ymax=1.,
                                               _zmin=0., _zmax=1.) );

    std::cout<<"create meshes done in " << toc() <<"s\n";

    if ( Dim == 2 )
    {
        using namespace boost::assign;
        outside1 += 1,2,4;
        outside2 += 2,3,4;
        gamma1 = 3;
        gamma2 = 1;
    }
    else if( Dim == 3 )
    {
        using namespace boost::assign;
        outside1 += 6,15,19,23,28;
        outside2 += 6,15,23,27,28;
        gamma1 = 27;
        gamma2 = 19;
    }

    auto trace_mesh = mesh1->trace( markedfaces(mesh1,gamma1) );
    auto mesh = fusion::make_vector(mesh1, mesh2, trace_mesh );

    /**
     * The function space and some associated elements(functions) are then defined
     */
    std::cout<<"function space construct starts\n";
    tic();
    auto Xh = space_type::New( _mesh=mesh );
    std::cout<<"function space construct done in "<< toc() <<"s\n";

    value_type pi = M_PI;
    auto g = sin(pi*Px())*cos(pi*Py())*cos(pi*Pz());

    auto gradg = trans( +pi*cos(pi*Px())*cos(pi*Py())*cos(pi*Pz())*unitX()
                        -pi*sin(pi*Px())*sin(pi*Py())*cos(pi*Pz())*unitY()
                        -pi*sin(pi*Px())*cos(pi*Py())*sin(pi*Pz())*unitZ() );

    auto f = pi*pi*Dim*g;

    bool weakdir = this->vm()["weakdir"].template as<int>();
    value_type penaldir = this->vm()["penaldir"].template as<double>();
    value_type coeff = this->vm()["coeff"].template as<double>();

    u = Xh->element();
    auto u1 = u.template element<0>();
    auto u2 = u.template element<1>();
    auto mu = u.template element<2>();

    auto v = Xh->element();
    auto v1 = v.template element<0>();
    auto v2 = v.template element<1>();
    auto nu = v.template element<2>();

    auto F = M_backend->newVector( Xh );

    std::cout<<"assembly_F starts\n";
    tic();
    form1( _test=Xh, _vector=F, _init=true );

    form1( _test=Xh, _vector=F ) +=
        integrate(elements(Xh->template mesh<0>()),f*id(v1) );

    form1( _test=Xh, _vector=F ) +=
        integrate(elements(Xh->template mesh<1>()),f*id(v2) );


    for ( int marker : outside1 )
    {
        form1( _test=Xh, _vector=F ) +=
            integrate( markedfaces(Xh->template mesh<0>(),marker),
                       g*(-grad(v1)*vf::N()+penaldir*id(v1)/hFace()) );
    }

    for( int marker : outside2 )
    {
        form1( _test=Xh, _vector=F ) +=
            integrate( markedfaces(Xh->template mesh<1>(),marker),
                       g*(-grad(v2)*vf::N()+penaldir*id(v2)/hFace()) );
    }
    std::cout<<"assembly_F done in " << toc() <<"s\n";
    F->close();
    F->printMatlab("F.m");

    auto A = M_backend->newMatrix( _trial=Xh, _test=Xh );
    std::cout<<"assembly_A starts\n";
    tic();
    form2( _trial=Xh, _test=Xh, _matrix=A, _init=true );

    LOG(INFO) << "A...";
    form2( _trial=Xh, _test=Xh, _matrix=A ) +=
        integrate( elements(Xh->template mesh<0>()), coeff*gradt(u1)*trans(grad(v1)) );

    form2( _trial=Xh, _test=Xh, _matrix=A ) +=
        integrate( elements(Xh->template mesh<1>()), coeff*gradt(u2)*trans(grad(v2)) );

    LOG(INFO) << "outside1...";
    for( int marker : outside1 )
    {
        form2( _trial=Xh, _test=Xh, _matrix=A ) +=
            integrate( markedfaces(Xh->template mesh<0>(),marker),
                       -(gradt(u1)*vf::N())*id(v1)
                       -(grad(v1)*vf::N())*idt(u1)
                       +penaldir*id(v1)*idt(u1)/hFace());
    }

    LOG(INFO) << "outside2...";
    for( int marker : outside2 )
    {
        form2( _trial=Xh, _test=Xh, _matrix=A ) +=
            integrate( markedfaces(Xh->template mesh<1>(),marker),
                       -(gradt(u2)*vf::N())*id(v2)
                       -(grad(v2)*vf::N())*idt(u2)
                       +penaldir*id(v2)*idt(u2)/hFace());
    }
    LOG(INFO) << "B, B^T..." ;
    form2( _trial=Xh, _test=Xh, _matrix=A ) +=
        integrate( markedfaces(mesh1,gamma1), idt(u1)*id(nu) + idt(mu)*id(v1)
                                              -idt(u2)*id(nu)-idt(mu)*id(v2) );

    std::cout<<"assembly_A done in " << toc() <<"s\n";

    A->close();
    A->printMatlab("A.m");

    std::cout<<"solve starts\n";
    tic();
    M_backend->solve(_matrix=A, _solution=u, _rhs=F, _pcfactormatsolverpackage="mumps");
    std::cout<<"solve done in " << toc() <<"s\n";

    std::cout<<"exportResults starts\n";
    tic();
    this->exportResults();
    std::cout<<"exportResults done in " << toc() <<"s\n";

    std::cout<<"computeErrors starts\n";
    tic();
    this->computeErrors();
    std::cout<<"computeErrors done in " << toc() <<"s\n";

} // MortarProd::run

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="mortar_prod",
                                  _author="Abdoulaye Samake",
                                  _email="abdoulaye.samake@imag.fr") );
    Application app;

    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    app.add( new MortarProd<2,1,1>() );
    app.run();
}
