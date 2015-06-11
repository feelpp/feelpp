
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake@imag.fr>
             Date: 2011-08-24

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
   \file mortar_lag.cpp
   \author Abdoulaye Samake <abdoulaye.samake@imag.fr>
   \date 2011-08-24
 */

#include <feel/feel.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;

/**
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description mortarlagoptions("MortarLag options");
    mortarlagoptions.add_options()
        ("hsize1", po::value<double>()->default_value( 0.015 ), "mesh size for first domain")
        ("hsize2", po::value<double>()->default_value( 0.02 ), "mesh size for second domain")
        ("shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
        ("coeff", po::value<double>()->default_value( 1 ), "grad.grad coefficient")
        ("weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ("penaldir", Feel::po::value<double>()->default_value( 10 ),
         "penalisation parameter for the weak boundary Dirichlet formulation")
        ;
    return mortarlagoptions.add( Feel::feel_options() );
}


namespace Feel
{
template<int Dim,int Order1,int Order2>
class MortarLag
    :
    public Simget
{
    typedef Simget super;
public:

    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef Mesh< Simplex<Dim,1,Dim> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef typename mesh_type::trace_mesh_type trace_mesh_type;
    typedef typename mesh_type::trace_mesh_ptrtype trace_mesh_ptrtype;
    typedef bases<Lagrange<Order1,Scalar> > basis1_type;
    typedef bases<Lagrange<Order2,Scalar> > basis2_type;
    typedef FunctionSpace<mesh_type, basis1_type> space1_type;
    typedef FunctionSpace<mesh_type, basis2_type> space2_type;
    typedef FunctionSpace<trace_mesh_type, basis1_type, Periodicity<NoPeriodicity>, mortars<Mortar> > lagmult_space_type;
    //typedef typename space1_type::trace_functionspace_type lagmult_space_type;
    typedef typename space1_type::element_type element1_type;
    typedef typename space2_type::element_type element2_type;
    typedef typename lagmult_space_type::element_type trace_element_type;
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef Exporter<trace_mesh_type> trace_export_type;
    typedef boost::shared_ptr<trace_export_type> trace_export_ptrtype;

    /**
     * Constructor
     */
    MortarLag()
        :
        super(),
        M_backend( backend_type::build( this->vm() ) ),
        mesh1Size( this->vm()["hsize1"].template as<double>() ),
        mesh2Size( this->vm()["hsize2"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() ),
        M_firstExporter( export_type::New( this->vm(),
                                           (boost::format( "%1%-%2%-%3%" )
                                            % this->about().appName()
                                            % Dim
                                            % int(1)).str() ) ),
        M_secondExporter( export_type::New( this->vm(),
                                            (boost::format( "%1%-%2%-%3%" )
                                             % this->about().appName()
                                             % Dim
                                             % int(2)).str() ) ),
        M_trace_exporter( trace_export_type::New( this->vm(),
                                            (boost::format( "%1%-%2%-%3%" )
                                             % this->about().appName()
                                             % Dim
                                             % int(3)).str() ) ),
        timers()
    {}

    mesh_ptrtype createMesh(  double xmin, double xmax, double meshsize, int id );

    void exportResults( element1_type& u,element2_type& v, trace_element_type& t );

    void run();

private:

    backend_ptrtype M_backend;
    double mesh1Size;
    double mesh2Size;
    std::string shape;
    export_ptrtype M_firstExporter;
    export_ptrtype M_secondExporter;
    trace_export_ptrtype M_trace_exporter;
    std::map<std::string, std::pair<boost::timer, double> > timers;
    std::vector<int> outside1;
    std::vector<int> outside2;
    int gamma1;
    int gamma2;

}; // MortarLag

template<int Dim,int Order1,int Order2>
typename MortarLag<Dim,Order1,Order2>::mesh_ptrtype
MortarLag<Dim,Order1,Order2>::createMesh(  double xmin, double xmax, double meshsize, int id )
{

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                        _desc=domain( _name=(boost::format( "%1%-%2%-%3%" ) % shape % Dim % id).str() ,
                                                      _addmidpoint=false,
                                                      _usenames=false,
                                                      _shape=shape,
                                                      _dim=Dim,
                                                      _h=meshsize,
                                                      _xmin=xmin,
                                                      _xmax=xmax,
                                                      _ymin=0.,
                                                      _ymax=1.,
                                                      _zmin=0.,
                                                      _zmax=1.) );

        return mesh;
}

template<int Dim,int Order1,int Order2>
void
MortarLag<Dim,Order1,Order2>::exportResults( element1_type& u, element2_type& v, trace_element_type& t )
{

    auto Xh1=u.functionSpace();
    auto mesh1=Xh1->mesh();
    auto Xh2=v.functionSpace();
    auto mesh2=Xh2->mesh();
    auto trace_mesh = mesh1->trace( markedfaces(mesh1,gamma1) );

    double pi = M_PI;
    using namespace vf;
    auto g = sin(pi*Px())*cos(pi*Py())*cos(pi*Pz());

    auto e1 = Xh1->element();
    e1 = vf::project( Xh1, elements(mesh1), g );

    auto e2 = Xh2->element();
    e2 = vf::project( Xh2, elements(mesh2), g );

    std::cout << "exportResults starts\n";
    timers["export"].first.restart();

    M_firstExporter->step(0)->setMesh( mesh1 );
    M_firstExporter->step(0)->add( "solution", (boost::format( "solution-%1%" ) % int(1) ).str(), u );
    M_firstExporter->step(0)->add( "exact", (boost::format( "exact-%1%" ) % int(1) ).str(), e1 );
    M_firstExporter->save();

    M_secondExporter->step(0)->setMesh( mesh2 );
    M_secondExporter->step(0)->add( "solution",(boost::format( "solution-%1%" ) % int(2) ).str(), v );
    M_secondExporter->step(0)->add( "exact",(boost::format( "exact-%1%" ) % int(2) ).str(), e2 );
    M_secondExporter->save();

    M_trace_exporter->step(0)->setMesh( trace_mesh );
    M_trace_exporter->step(0)->add( "lambda",(boost::format( "lambda-%1%" ) % int(3) ).str(), t );
    M_trace_exporter->save();

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
            ofs << "casefile: mortar-lag-" << Dim << "-" << j << ".case\n";//"-1_0.case\n";
        }
    }
    std::cout << "exportResults done\n";
    timers["export"].second = timers["export"].first.elapsed();
    std::cout << "[timer] exportResults(): " << timers["export"].second << "s\n";
} // MortarLag::export


template<int Dim,int Order1,int Order2>
void
MortarLag<Dim,Order1,Order2>::run()
{

    //if ( !this->vm().count( "nochdir" ) )
    Environment::changeRepository( boost::format( "doc/manual/%1%/%2%-%3%/P%4%-P%5%/h_%6%-%7%/" )
                                   % this->about().appName()
                                   % shape
                                   % Dim
                                   % Order1
                                   % Order2
                                   % mesh1Size
                                   % mesh2Size );

    mesh_ptrtype mesh1 = createMesh(0.,0.5,mesh1Size,1);

    mesh_ptrtype mesh2 = createMesh(0.5,1.,mesh2Size,2);


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

    /**
     * The function space and some associated elements(functions) are then defined
     */
    auto Xh1 = space1_type::New( mesh1 );
    auto u1 = Xh1->element();
    auto v1 = Xh1->element();

    auto trace_mesh = mesh1->trace( markedfaces(mesh1,gamma1) );
    auto Lh1 = lagmult_space_type::New( trace_mesh );
    auto mu = Lh1->element();
    auto nu = Lh1->element();

    auto Xh2 = space2_type::New( mesh2 );
    auto u2 = Xh2->element();
    auto v2 = Xh2->element();

    value_type pi = M_PI;

    auto g = sin(pi*Px())*cos(pi*Py())*cos(pi*Pz());

    auto gradg = trans( +pi*cos(pi*Px())*cos(pi*Py())*cos(pi*Pz())*unitX()
                        -pi*sin(pi*Px())*sin(pi*Py())*cos(pi*Pz())*unitY()
                        -pi*sin(pi*Px())*cos(pi*Py())*sin(pi*Pz())*unitZ() );

    auto f = pi*pi*Dim*g;

    bool weakdir = ioption("weakdir");
    value_type penaldir = doption("penaldir");
    value_type coeff = doption("coeff");

    auto F1 = M_backend->newVector( Xh1 );
    form1( _test=Xh1, _vector=F1, _init=true ) =
        integrate( elements(mesh1), 10000*f*id(v1) );

    BOOST_FOREACH( int marker, outside1 )
    {
        form1( _test=Xh1, _vector=F1 ) +=
            integrate( markedfaces(mesh1,marker),
                       10000*g*(-grad(v1)*vf::N()+penaldir*id(v1)/hFace()) );
    }

    F1->close();

    auto D1 = M_backend->newMatrix( _trial=Xh1, _test=Xh1 );

    form2( _trial=Xh1, _test=Xh1, _matrix=D1, _init=true ) =
        integrate( elements(mesh1), 10000*coeff*gradt(u1)*trans(grad(v1)) );

    BOOST_FOREACH( int marker, outside1 )
    {
        form2( _trial=Xh1, _test=Xh1, _matrix=D1 ) +=
            integrate( markedfaces(mesh1,marker),
                       10000*(-(gradt(u1)*vf::N())*id(v1)
                           -(grad(v1)*vf::N())*idt(u1)
                           +penaldir*id(v1)*idt(u1)/hFace()));
    }

    D1->close();

    auto B1 = M_backend->newMatrix( _trial=Xh1, _test=Lh1 );

    std::cout << "assembly_B1 starts\n";
    timers["assemby_B1"].first.restart();

    form2( _trial=Xh1, _test=Lh1, _matrix=B1, _init=true) =
        integrate( markedfaces(mesh1,gamma1), idt(u1)*id(nu) );

    timers["assemby_B1"].second = timers["assemby_B1"].first.elapsed();
    std::cout << "assemby_B1 done in " << timers["assemby_B1"].second << "s\n";

    B1->close();

    auto F2 = M_backend->newVector( Xh2 );
    form1( _test=Xh2, _vector=F2, _init=true ) =
        integrate( elements(mesh2), f*id(v2) );

    BOOST_FOREACH( int marker, outside2 )
    {
        form1( _test=Xh2, _vector=F2 ) +=
            integrate( markedfaces(mesh2,marker),
                       g*(-grad(v2)*vf::N()+penaldir*id(v2)/hFace()) );
    }

    F2->close();

    auto D2 = M_backend->newMatrix( _trial=Xh2, _test=Xh2 );

    form2( _trial=Xh2, _test=Xh2, _matrix=D2, _init=true ) =
        integrate( elements(mesh2), coeff*gradt(u2)*trans(grad(v2)) );

    BOOST_FOREACH( int marker, outside2 )
    {
        form2( _trial=Xh2, _test=Xh2, _matrix=D2 ) +=
            integrate( markedfaces(mesh2,marker),
                       -(gradt(u2)*vf::N())*id(v2)
                       -(grad(v2)*vf::N())*idt(u2)
                       +penaldir*id(v2)*idt(u2)/hFace());
    }

    D2->close();

    auto B2 = M_backend->newMatrix( _trial=Xh2, _test=Lh1, _buildGraphWithTranspose=true );

    std::cout << "assembly_B2 starts\n";
    timers["assembly_B2"].first.restart();

    form2( _trial=Xh2, _test=Lh1, _matrix=B2, _init=true ) =
        integrate( markedfaces(mesh1,gamma1), -idt(u2)*id(nu) );

    timers["assembly_B2"].second = timers["assembly_B2"].first.elapsed();
    std::cout << "assembly_B2 done in " << timers["assembly_B2"].second << "s\n";

    B2->close();

    auto B12 = M_backend->newZeroMatrix( _trial=Xh2, _test=Xh1 );

    auto B21 = M_backend->newZeroMatrix( _trial=Xh1, _test=Xh2 );

    auto BLL = M_backend->newZeroMatrix( _trial=Lh1, _test=Lh1 );

    auto B1t = M_backend->newMatrix( _trial=Lh1, _test=Xh1, _buildGraphWithTranspose=true );

    B1->transpose(B1t);

    auto B2t = M_backend->newMatrix( _trial=Lh1, _test=Xh2, _buildGraphWithTranspose=true );

    B2->transpose(B2t);

    auto myb = BlocksSparseMatrix<3,3>()<< D1 << B12 << B1t
                                        << B21 << D2 << B2t
                                        << B1 << B2  << BLL ;

    auto AbB = M_backend->newBlockMatrix(myb);
    AbB->close();

    auto FbB = M_backend->newVector( u1.size()+u2.size()+mu.size(),u1.size()+u2.size()+mu.size() );
    auto UbB = M_backend->newVector( u1.size()+u2.size()+mu.size(),u1.size()+u2.size()+mu.size() );

    for (size_type i = 0 ; i < F1->size(); ++ i)
        FbB->set(i, (*F1)(i) );

    for (size_type i = 0 ; i < F2->size(); ++ i)
        FbB->set(F1->size()+i, (*F2)(i) );

    std::cout << "number of dof(u1): " << Xh1->nDof() << "\n";
    std::cout << "number of dof(u2): " << Xh2->nDof() << "\n";
    std::cout << "number of dof(lambda): " << Lh1->nDof() << "\n";
    std::cout << "size of linear system: " << FbB->size() << "\n";

    std::cout << "solve starts\n";
    timers["solve"].first.restart();

    M_backend->solve(_matrix=AbB,
                     _solution=UbB,
                     _rhs=FbB,
                     _pcfactormatsolverpackage="umfpack");

    timers["solve"].second = timers["solve"].first.elapsed();
    std::cout << "solve done in " << timers["solve"].second << "s\n";

    for (size_type i = 0 ; i < u1.size(); ++ i)
        u1.set(i, (*UbB)(i) );

    for (size_type i = 0 ; i < u2.size(); ++ i)
        u2.set(i, (*UbB)(u1.size()+i) );

    for (size_type i = 0 ; i < mu.size(); ++ i)
    mu.set(i, (*UbB)(u1.size()+u2.size()+i) );

    double L2error12 =integrate(elements(mesh1),(idv(u1)-g)*(idv(u1)-g) ).evaluate()(0,0);
    double L2error1 =   math::sqrt( L2error12 );

    double L2error22 =integrate(elements(mesh2),(idv(u2)-g)*(idv(u2)-g) ).evaluate()(0,0);
    double L2error2 =   math::sqrt( L2error22 );

    double semi_H1error1 =integrate(elements(mesh1),
                                    ( gradv(u1)-gradg )*trans( (gradv(u1)-gradg) ) ).evaluate()(0,0);

    double semi_H1error2 =integrate(elements(mesh2),
                                    ( gradv(u2)-gradg )*trans( (gradv(u2)-gradg) ) ).evaluate()(0,0);

    double H1error1 = math::sqrt( L2error12 + semi_H1error1 );

    double H1error2 = math::sqrt( L2error22 + semi_H1error2 );

    double error =integrate(elements(trace_mesh), (idv(u1)-idv(u2))*(idv(u1)-idv(u2)) ).evaluate()(0,0);

    double global_error = math::sqrt(L2error12 + L2error22 + semi_H1error1 + semi_H1error2);

    std::cout << "----------L2 errors---------- \n" ;
    std::cout << "||u1_error||_L2=" << L2error1 << "\n";
    std::cout << "||u2_error||_L2=" << L2error2 << "\n";
    std::cout << "----------H1 errors---------- \n" ;
    std::cout << "||u1_error||_H1=" << H1error1 << "\n";
    std::cout << "||u2_error||_H1=" << H1error2 << "\n";
    std::cout << "||u_error||_H1=" << global_error << "\n";
    std::cout << "L2 norm of jump at interface  \n" ;
    std::cout << "||u1-u2||_L2=" << math::sqrt(error) << "\n";

    this->exportResults(u1,u2,mu);

} // MortarLag::run

} // Feel
/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="mortar-lag",
                                  _author="Abdoulaye Samake",
                                  _email="abdoulaye.samake@imag.fr") );
    Application app;

    app.add( new MortarLag<2,2,2>() ) ;
    app.add( new MortarLag<2,2,3>() ) ;
    app.run();
}
