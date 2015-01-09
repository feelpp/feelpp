/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <Abdoulaye.Samake@imag.fr>
  Date: 2011-12-15

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file threefields.cpp
   \author Abdoulaye Samake <Abdoulaye.Samake@imag.fr>
   \date 2011-12-15
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelalg/matrixblock.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>
#include <feel/feelpoly/polynomialset.hpp>

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
    po::options_description threefieldsoptions( "ThreeFields options" );
    threefieldsoptions.add_options()
    ( "hsize1", po::value<double>()->default_value( 0.03 ), "mesh size for first domain" )
    ( "hsize2", po::value<double>()->default_value( 0.02 ), "mesh size for second domain" )
    ( "hsize3", po::value<double>()->default_value( 0.01 ), "mesh size for interface" )
    ( "split", po::value<double>()->default_value( 0.5 ), "interface between subdomains" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "coeff", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
    ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
    ( "penaldir", Feel::po::value<double>()->default_value( 10 ),
      "penalisation parameter for the weak boundary Dirichlet formulation" ),
      ( "export-matlab", "export matrix and vectors in matlab" )
      ;
    return threefieldsoptions.add( Feel::feel_options() );
}

/**
 * \return some data about the application.
 */
inline
AboutData
makeAbout()
{
    AboutData about( "threefields" ,
                     "threefields" ,
                     "0.2",
                     "nD(n=2,3) ThreeFields using threefields",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Abdoulaye Samake", "developer", "Abdoulaye.Samake@imag.fr", "" );
    return about;

}

/**
 * \class ThreeFieldsLaplacian
 *
 * ThreeFieldsLaplacian Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=2 or 3)
 */
template<int Dim, int Order1, int Order2, int Order3>
class ThreeFieldsLaplacian
    :
public Simget
{
    typedef Simget super;
public:

    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef Simplex<Dim,1,Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef typename mesh_type::trace_mesh_type trace_mesh_type;
    typedef typename mesh_type::trace_mesh_ptrtype trace_mesh_ptrtype;
    typedef bases<Lagrange<Order1,Scalar> > basis1_type;
    typedef bases<Lagrange<Order2,Scalar> > basis2_type;
    typedef FunctionSpace<mesh_type, basis1_type> space1_type;
    typedef FunctionSpace<mesh_type, basis2_type> space2_type;
    typedef typename space1_type::trace_functionspace_type trace1_space_type;
    typedef typename space2_type::trace_functionspace_type trace2_space_type;
    typedef typename space1_type::element_type element1_type;
    typedef typename space2_type::element_type element2_type;
    typedef typename trace1_space_type::element_type trace1_element_type;
    typedef typename trace2_space_type::element_type trace2_element_type;
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef Exporter<trace_mesh_type> trace_export_type;
    typedef boost::shared_ptr<trace_export_type> trace_export_ptrtype;

    // @name typedef for interfaces

    typedef bases<Lagrange<Order3,Scalar> > interfaces_basis_type;
    typedef FunctionSpace<trace_mesh_type, interfaces_basis_type> interfaces_space_type;
    typedef boost::shared_ptr<interfaces_space_type> interfaces_space_ptrtype;
    typedef typename interfaces_space_type::element_type interfaces_element_type;
    /**
     * Constructor
     */
    ThreeFieldsLaplacian()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        split( doption("split") ),
        mesh1Size( doption("hsize1") ),
        mesh2Size( doption("hsize2") ),
        mesh3Size( doption("hsize3") ),
        shape( soption("shape") ),
        timers(),
        M_firstExporter( export_type::New( this->vm(),
                                           ( boost::format( "%1%-%2%-%3%" )
                                             % this->about().appName()
                                             % Dim
                                             % int( 1 ) ).str() ) ),
        M_secondExporter( export_type::New( this->vm(),
                                            ( boost::format( "%1%-%2%-%3%" )
                                              % this->about().appName()
                                              % Dim
                                              % int( 2 ) ).str() ) ),
        M_trace1_exporter( trace_export_type::New( this->vm(),
                           ( boost::format( "%1%-%2%-%3%" )
                             % this->about().appName()
                             % Dim
                             % int( 3 ) ).str() ) ),
        M_trace2_exporter( trace_export_type::New( this->vm(),
                           ( boost::format( "%1%-%2%-%3%" )
                             % this->about().appName()
                             % Dim
                             % int( 4 ) ).str() ) )
    {}

    mesh_ptrtype createMesh(  double xmin, double xmax, double meshsize, int id );

    trace_mesh_ptrtype createMesh( double meshSize, double interface );

    void exportResults( element1_type& u,element2_type& v, trace1_element_type& t1, trace2_element_type& t2 );

    void run();

private:

    backend_ptrtype M_backend;
    double split;
    double mesh1Size;
    double mesh2Size;
    double mesh3Size;
    std::string shape;
    std::map<std::string, std::pair<boost::timer, double> > timers;
    export_ptrtype M_firstExporter;
    export_ptrtype M_secondExporter;
    trace_export_ptrtype M_trace1_exporter,M_trace2_exporter;
    std::vector<int> outside1;
    std::vector<int> outside2;
    int gamma1;
    int gamma2;

}; // ThreeFieldsLaplacian

template<int Dim, int Order1, int Order2, int Order3>
typename ThreeFieldsLaplacian<Dim, Order1, Order2, Order3>::mesh_ptrtype
ThreeFieldsLaplacian<Dim, Order1, Order2, Order3>::createMesh(  double xmin, double xmax, double meshsize, int id )
{

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                        _desc=domain( _name=( boost::format( "%1%-%2%-%3%" ) % shape % Dim % id ).str() ,
                                                      _addmidpoint=false,
                                                      _usenames=false,
                                                      _shape=this->shape,
                                                      _dim=Dim,
                                                      _h=meshsize,
                                                      _xmin=xmin,
                                                      _xmax=xmax,
                                                      _ymin=0.,
                                                      _ymax=1.,
                                                      _zmin=0.,
                                                      _zmax=1. ) );

    return mesh;

} // ThreeFieldsLaplacian::createMesh

template<int Dim, int Order1, int Order2, int Order3>
typename ThreeFieldsLaplacian<Dim, Order1, Order2, Order3>::trace_mesh_ptrtype
ThreeFieldsLaplacian<Dim, Order1, Order2, Order3>::createMesh( double meshSize, double interface )
{
    trace_mesh_ptrtype mesh( new trace_mesh_type );

    Gmsh __gmsh( Dim-1 );
    __gmsh.setOrder( GMSH_ORDER_ONE );
    std::ostringstream ostr;
    std::ostringstream nameStr;
    __gmsh.setRecombine( false );
    __gmsh.setCharacteristicLength( meshSize );
    ostr << __gmsh.preamble() << "\n";

    switch ( Dim )
    {
    case 2:
        ostr << "a=" << interface << ";\n"
             << "Point(1) = {a,0,0,h};\n"
             << "Point(2) = {a,1,0,h};\n"
             << "Line(1) = {1,2};\n"
             << "Physical Point(1) = {1};\n"
             << "Physical Point(2) = {2};\n"
             << "Physical Line(\"Mat1\") = {1};\n";
        nameStr << "interface2D";
        break;

    case 3:
        ostr << "a=" << interface << ";\n"
             << "Point(1) = {a,0,0,h};\n"
             << "Point(2) = {a,1,0,h};\n"
             << "Point(3) = {a,1,1,h};\n"
             << "Point(4) = {a,0,1,h};\n"
             << "Line(1) = {1,2};\n"
             << "Line(2) = {2,3};\n"
             << "Line(3) = {3,4};\n"
             << "Line(4) = {4,1};\n"
             << "Line Loop(5) = {1,2,3,4};\n"
             << "Plane Surface(6) = {5};\n"
             << "Physical Surface(\"Mat1\") = {6};\n";
        nameStr << "interface3D";
        break;

    default:
        std::ostringstream os;
        os << "invalid dimension: " << Dim;
        throw std::logic_error( os.str() );
    }

    std::string fname = __gmsh.generate( nameStr.str(), ostr.str() ).template get<0>();
    ImporterGmsh<trace_mesh_type> import( fname );
    mesh->accept( import );

    mesh->components().set( MESH_RENUMBER | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
    mesh->updateForUse();
    return mesh;

} // ThreeFieldsLaplacian::createMesh

template<int Dim, int Order1, int Order2, int Order3>
void
ThreeFieldsLaplacian<Dim, Order1, Order2, Order3>::exportResults( element1_type& u, element2_type& v, trace1_element_type& t1, trace2_element_type& t2 )
{
    auto Xh1=u.functionSpace();
    auto mesh1=Xh1->mesh();
    auto Xh2=v.functionSpace();
    auto mesh2=Xh2->mesh();
    auto trace_mesh1 = mesh1->trace( markedfaces( mesh1,gamma1 ) );
    auto trace_mesh2 = mesh2->trace( markedfaces( mesh2,gamma2 ) );

    double pi = M_PI;
    using namespace vf;
    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );

    auto e1 = Xh1->element();
    e1 = vf::project( Xh1, elements( mesh1 ), g );

    auto e2 = Xh2->element();
    e2 = vf::project( Xh2, elements( mesh2 ), g );

    LOG(INFO) << "exportResults starts\n";
    timers["export"].first.restart();

    M_firstExporter->step( 0 )->setMesh( mesh1 );
    M_firstExporter->step( 0 )->add( "solution", ( boost::format( "solution-%1%" ) % int( 1 ) ).str(), u );
    M_firstExporter->step( 0 )->add( "exact", ( boost::format( "exact-%1%" ) % int( 1 ) ).str(), e1 );
    M_firstExporter->save();

    M_secondExporter->step( 0 )->setMesh( mesh2 );
    M_secondExporter->step( 0 )->add( "solution",( boost::format( "solution-%1%" ) % int( 2 ) ).str(), v );
    M_secondExporter->step( 0 )->add( "exact",( boost::format( "exact-%1%" ) % int( 2 ) ).str(), e2 );
    M_secondExporter->save();

    M_trace1_exporter->step( 0 )->setMesh( trace_mesh1 );
    M_trace1_exporter->step( 0 )->add( "trace1",( boost::format( "solution-%1%" ) % int( 3 ) ).str(), t1 );
    M_trace1_exporter->save();

    M_trace2_exporter->step( 0 )->setMesh( trace_mesh2 );
    M_trace2_exporter->step( 0 )->add( "trace2",( boost::format( "solution-%1%" ) % int( 4 ) ).str(), t2 );
    M_trace2_exporter->save();

    std::ofstream ofs( ( boost::format( "%1%.sos" ) % this->about().appName() ).str().c_str() );

    if ( ofs )
    {
        ofs << "FORMAT:\n"
            << "type: master_server gold\n"
            << "SERVERS\n"
            << "number of servers: " << int( 2 ) << "\n";

        for ( int j = 1; j <= 2; ++ j )
        {
            ofs << "#Server " << j << "\n";
            ofs << "machine id: " << mpi::environment::processor_name()  << "\n";
            ofs << "executable:\n";
            ofs << "data_path: .\n";
            ofs << "casefile: threefields-" << Dim << "-" << j << "-1_0.case\n";
        }
    }

    LOG(INFO) << "exportResults done\n";
    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] exportResults(): " << timers["export"].second << "\n";

} // ThreeFieldsLaplacian::export


template<int Dim, int Order1, int Order2, int Order3>
void
ThreeFieldsLaplacian<Dim, Order1, Order2, Order3>::run()
{
    LOG(INFO) << "Execute ThreeFieldsLaplacian<" << Dim << "," << Order1 << "," << Order2 << "," << Order3 << ">\n";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "doc/manual/%1%/%2%-%3%/P%4%-P%5%-P%6%/h_%7%-%8%-%9%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order1
                                       % Order2
                                       % Order3
                                       % mesh1Size
                                       % mesh2Size
                                       % mesh3Size );

    LOG(INFO) << "create meshes starts\n";
    mesh_ptrtype mesh1 = createMesh( 0.,split,mesh1Size,1 );
    mesh_ptrtype mesh2 = createMesh( split,1.,mesh2Size,2 );

    auto trace_mesh1 = mesh1->trace( markedfaces( mesh1,gamma1 ) );
    auto trace_mesh2 = mesh2->trace( markedfaces( mesh2,gamma2 ) );

    trace_mesh_ptrtype interface_mesh = createMesh( mesh3Size,split );

    LOG(INFO) << "create meshes done\n";

    if ( Dim == 2 )
    {
        using namespace boost::assign;
        outside1 += 1,2,4;
        outside2 += 2,3,4;
        gamma1 = 3;
        gamma2 = 1;
    }

    else if ( Dim == 3 )
    {
        using namespace boost::assign;
        outside1 += 6,15,19,23,28;
        outside2 += 6,15,23,27,28;
        gamma1 = 27;
        gamma2 = 19;
    }

    auto Xh1 = space1_type::New( mesh1 );
    auto u1 = Xh1->element();
    auto v1 = Xh1->element();

    auto Xh2 = space2_type::New( mesh2 );
    auto u2 = Xh2->element();
    auto v2 = Xh2->element();

    auto Lh1 = trace1_space_type::New( mesh1->trace( markedfaces( mesh1,gamma1 ) ) );
    auto mu1 = Lh1->element();
    auto nu1 = Lh1->element();

    auto Lh2 = trace2_space_type::New( mesh2->trace( markedfaces( mesh2,gamma2 ) ) );
    auto mu2 = Lh2->element();
    auto nu2 = Lh2->element();

    auto Lh = interfaces_space_type::New( interface_mesh );
    auto mu = Lh->element();
    auto nu = Lh->element();

    // buildGraphWithTranspose ask
    bool buildGraphWithTrans1=false, buildGraphWithTrans2=false;

    if ( mu.size() > mu1.size() ) buildGraphWithTrans1=true;

    if ( mu.size() > mu2.size() ) buildGraphWithTrans2=true;

    value_type pi = M_PI;
    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );

    auto gradg = trans( +pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() )*unitX()
                        -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() )*unitY()
                        -pi*sin( pi*Px() )*cos( pi*Py() )*sin( pi*Pz() )*unitZ() );

    auto f = pi*pi*Dim*g;

    bool weakdir = ioption("weakdir");
    value_type penaldir = doption("penaldir");
    value_type coeff = doption("coeff");

    LOG(INFO) << "assembly_F1 starts\n";
    timers["assemby_F1"].first.restart();

    auto F1 = M_backend->newVector( Xh1 );
    form1( _test=Xh1, _vector=F1, _init=true ) =
        integrate( elements( mesh1 ), f*id( v1 ) );

    BOOST_FOREACH( int marker, outside1 )
    {
        form1( _test=Xh1, _vector=F1 ) +=
            integrate( markedfaces( mesh1,marker ),
                       g*( -grad( v1 )*vf::N()+penaldir*id( v1 )/hFace() ) );
    }

    timers["assemby_F1"].second = timers["assemby_F1"].first.elapsed();
    LOG(INFO) << "assemby_F1 done in " << timers["assemby_F1"].second << "s\n";

    F1->close();

    LOG(INFO) << "assembly_D1 starts\n";
    timers["assemby_D1"].first.restart();

    auto D1 = M_backend->newMatrix( Xh1, Xh1 );

    form2( _trial=Xh1, _test=Xh1, _matrix=D1, _init=true ) =
        integrate( elements( mesh1 ), coeff*gradt( u1 )*trans( grad( v1 ) ) );

    BOOST_FOREACH( int marker, outside1 )
    {
        form2( _trial=Xh1, _test=Xh1, _matrix=D1 ) +=
            integrate( markedfaces( mesh1,marker ),
                       -( gradt( u1 )*vf::N() )*id( v1 )
                       -( grad( v1 )*vf::N() )*idt( u1 )
                       +penaldir*id( v1 )*idt( u1 )/hFace() );
    }

    timers["assemby_D1"].second = timers["assemby_D1"].first.elapsed();
    LOG(INFO) << "assemby_D1 done in " << timers["assemby_D1"].second << "s\n";

    D1->close();

    if ( this->vm().count( "export-matlab" ) )
        D1->printMatlab( "D1.m" );

    LOG(INFO) << "assembly_B1 starts\n";
    timers["assemby_B1"].first.restart();

    auto B1 = M_backend->newMatrix( Xh1, Lh1 );
    form2( _trial=Xh1, _test=Lh1, _matrix=B1, _init=true ) +=
        integrate( elements( Lh->mesh() ), -idt( u1 )*id( nu1 ) );

    timers["assemby_B1"].second = timers["assemby_B1"].first.elapsed();
    LOG(INFO) << "assemby_B1 done in " << timers["assemby_B1"].second << "s\n";

    B1->close();

    if ( this->vm().count( "export-matlab" ) )
        B1->printMatlab( "B1.m" );

    LOG(INFO) << "assembly_C1 starts\n";
    timers["assemby_C1"].first.restart();

    auto C1 = M_backend->newMatrix( Lh, Lh1, _buildGraphWithTranspose=!buildGraphWithTrans1 );
    form2( _trial=Lh, _test=Lh1, _matrix=C1, _init=true ) +=
        integrate( elements( Lh->mesh() ), idt( mu )*id( nu1 ) );

    timers["assemby_C1"].second = timers["assemby_C1"].first.elapsed();
    LOG(INFO) << "assemby_C1 done in " << timers["assemby_C1"].second << "s\n";

    C1->close();

    if ( this->vm().count( "export-matlab" ) )
        C1->printMatlab( "C1.m" );

    LOG(INFO) << "assembly_F2 starts\n";
    timers["assemby_F2"].first.restart();

    auto F2 = M_backend->newVector( Xh2 );
    form1( _test=Xh2, _vector=F2, _init=true ) =
        integrate( elements( mesh2 ), f*id( v2 ) );

    BOOST_FOREACH( int marker, outside2 )
    {
        form1( _test=Xh2, _vector=F2 ) +=
            integrate( markedfaces( mesh2,marker ),
                       g*( -grad( v2 )*vf::N()+penaldir*id( v2 )/hFace() ) );
    }

    timers["assemby_F2"].second = timers["assemby_F2"].first.elapsed();
    LOG(INFO) << "assemby_F2 done in " << timers["assemby_F2"].second << "s\n";

    F2->close();

    LOG(INFO) << "assembly_D2 starts\n";
    timers["assemby_D2"].first.restart();

    auto D2 = M_backend->newMatrix( Xh2, Xh2 );

    form2( _trial=Xh2, _test=Xh2, _matrix=D2, _init=true ) =
        integrate( elements( mesh2 ), coeff*gradt( u2 )*trans( grad( v2 ) ) );

    BOOST_FOREACH( int marker, outside2 )
    {
        form2( _trial=Xh2, _test=Xh2, _matrix=D2 ) +=
            integrate( markedfaces( mesh2,marker ),
                       -( gradt( u2 )*vf::N() )*id( v2 )
                       -( grad( v2 )*vf::N() )*idt( u2 )
                       +penaldir*id( v2 )*idt( u2 )/hFace() );
    }

    timers["assemby_D2"].second = timers["assemby_D2"].first.elapsed();
    LOG(INFO) << "assemby_D2 done in " << timers["assemby_D2"].second << "s\n";

    D2->close();
    LOG(INFO) << "matrix D2 assembly done\n";

    if ( this->vm().count( "export-matlab" ) )
        D2->printMatlab( "D2.m" );

    LOG(INFO) << "assembly_B2 starts\n";
    timers["assemby_B2"].first.restart();

    auto B2 = M_backend->newMatrix( Xh2, Lh2 );
    form2( _trial=Xh2, _test=Lh2, _matrix=B2, _init=true ) +=
        integrate( elements( Lh->mesh() ), -idt( u2 )*id( nu2 ) );

    timers["assemby_B2"].second = timers["assemby_B2"].first.elapsed();
    LOG(INFO) << "assemby_B2 done in " << timers["assemby_B2"].second << "s\n";

    B2->close();


    if ( this->vm().count( "export-matlab" ) )
        B2->printMatlab( "B2.m" );

    LOG(INFO) << "assembly_C2 starts\n";
    timers["assemby_C2"].first.restart();

    auto C2 = M_backend->newMatrix( Lh, Lh2,_buildGraphWithTranspose=!buildGraphWithTrans2 );
    form2( _trial=Lh, _test=Lh2, _matrix=C2, _init=true ) +=
        integrate( elements( Lh->mesh() ), idt( mu )*id( nu2 ) );

    timers["assemby_C2"].second = timers["assemby_C2"].first.elapsed();
    LOG(INFO) << "assemby_C2 done in " << timers["assemby_C2"].second << "s\n";

    C2->close();

    if ( this->vm().count( "export-matlab" ) )
        C2->printMatlab( "C2.m" );

    // transposes
    LOG(INFO) << "assembly_B1t starts\n";
    timers["assemby_B1t"].first.restart();

    auto B1t = M_backend->newMatrix( Lh1, Xh1, _buildGraphWithTranspose=true );

    B1->transpose( B1t );

    timers["assemby_B1t"].second = timers["assemby_B1t"].first.elapsed();
    LOG(INFO) << "assemby_B1t done in " << timers["assemby_B1t"].second << "s\n";

    if ( this->vm().count( "export-matlab" ) )
        B1t->printMatlab( "B1t.m" );

    LOG(INFO) << "assembly_B2t starts\n";
    timers["assemby_B2t"].first.restart();

    auto B2t = M_backend->newMatrix( Lh2, Xh2, _buildGraphWithTranspose=true );

    B2->transpose( B2t );

    timers["assemby_B2t"].second = timers["assemby_B2t"].first.elapsed();
    LOG(INFO) << "assemby_B2t done in " << timers["assemby_B2t"].second << "s\n";

    if ( this->vm().count( "export-matlab" ) )
        B2t->printMatlab( "B2t.m" );

    LOG(INFO) << "matrix C1t assembly starts\n";

    LOG(INFO) << "assemply_C1t starts\n";
    timers["assemply_C1t"].first.restart();

    auto C1t = M_backend->newMatrix( Lh1, Lh, _buildGraphWithTranspose=buildGraphWithTrans1 );

    C1->transpose( C1t );

    timers["assemply_C1t"].second = timers["assemply_C1t"].first.elapsed();
    LOG(INFO) << "[timer] assemply_C1t: " << timers["assemply_C1t"].second << "\n";

    if ( this->vm().count( "export-matlab" ) )
        C1t->printMatlab( "C1t.m" );

    LOG(INFO) << "assemply_C2t starts\n";
    timers["assemply_C2t"].first.restart();

    auto C2t = M_backend->newMatrix( Lh2, Lh, _buildGraphWithTranspose=buildGraphWithTrans2 );

    C2->transpose( C2t );

    timers["assemply_C2t"].second = timers["assemply_C2t"].first.elapsed();
    LOG(INFO) << "[timer] assemply_C2t: " << timers["assemply_C2t"].second << "\n";

    if ( this->vm().count( "export-matlab" ) )
        C2t->printMatlab( "C2t.m" );

    // zero matrices

    auto zero31 = M_backend->newZeroMatrix( _trial=Lh, _test=Xh1 );
    auto zero41 = M_backend->newZeroMatrix( _trial=Lh2, _test=Xh1 );
    auto zero51 = M_backend->newZeroMatrix( _trial=Xh2, _test=Xh1 );
    auto zero22 = M_backend->newZeroMatrix( _trial=Lh1, _test=Lh1 );
    auto zero42 = M_backend->newZeroMatrix( _trial=Lh2, _test=Lh1 );
    auto zero52 = M_backend->newZeroMatrix( _trial=Xh2, _test=Lh1 );
    auto zero13 = M_backend->newZeroMatrix( _trial=Xh1, _test=Lh );
    auto zero33 = M_backend->newZeroMatrix( _trial=Lh, _test=Lh  );
    auto zero53 = M_backend->newZeroMatrix( _trial=Xh2, _test=Lh );
    auto zero14 = M_backend->newZeroMatrix( _trial=Xh1, _test=Lh2 );
    auto zero24 = M_backend->newZeroMatrix( _trial=Lh1, _test=Lh2 );
    auto zero44 = M_backend->newZeroMatrix( _trial=Lh2, _test=Lh2 );
    auto zero15 = M_backend->newZeroMatrix( _trial=Xh1, _test=Xh2 );
    auto zero25 = M_backend->newZeroMatrix( _trial=Lh1, _test=Xh2 );
    auto zero35 = M_backend->newZeroMatrix( _trial=Lh, _test=Xh2 );

    LOG(INFO) << "assemply_Block starts\n";
    timers["assemply_Block"].first.restart();

    auto myb = BlocksSparseMatrix<5,5>()<< D1 << B1t << zero31 << zero41 << zero51
                                        << B1 << zero22 << C1 << zero42 << zero52
                                        << zero13 << C1t << zero33 << C2t << zero53
                                        << zero14 << zero24 << C2 << zero44 << B2
                                        << zero15 << zero25 << zero35 << B2t << D2 ;

    auto AbB = M_backend->newBlockMatrix( _block=myb );

    timers["assemply_Block"].second = timers["assemply_Block"].first.elapsed();
    LOG(INFO) << "[timer] assemply_Block: " << timers["assemply_Block"].second << "\n";

    AbB->close();

    LOG(INFO) <<"***********************************************\n";
    LOG(INFO) << "full matrix size1= " << AbB->size1() <<"\n";
    LOG(INFO) << "full matrix size2= " << AbB->size2() <<"\n";
    LOG(INFO) <<"***********************************************\n";

    auto FbB = M_backend->newVector( u1.size()+u2.size()+mu.size()+mu1.size()+mu2.size(),
                                     u1.size()+u2.size()+mu.size()+mu1.size()+mu2.size() );

    auto UbB = M_backend->newVector( u1.size()+u2.size()+mu.size()+mu1.size()+mu2.size(),
                                     u1.size()+u2.size()+mu.size()+mu1.size()+mu2.size() );


    for ( size_type i = 0 ; i < F1->size(); ++ i )
        FbB->set( i, ( *F1 )( i ) );

    for ( size_type i = 0 ; i < F2->size(); ++ i )
        FbB->set( u1.size()+mu1.size()+mu.size()+mu2.size()+i, ( *F2 )( i ) );


    LOG(INFO) <<"***********************************************\n";
    LOG(INFO) << "rhs size=      " << FbB->size() <<"\n";
    LOG(INFO) << "solution size= " << UbB->size() <<"\n";
    LOG(INFO) <<"***********************************************\n";

    if ( this->vm().count( "export-matlab" ) )
    {
        AbB->printMatlab( "AbB.m" );
        FbB->printMatlab( "FbB.m" );
        F1->printMatlab( "F1.m" );
        F2->printMatlab( "F2.m" );
    }

    LOG(INFO) << "solve starts\n";
    timers["solve"].first.restart();

    M_backend->solve( _matrix=AbB,
                      _solution=UbB,
                      _rhs=FbB,
                      _pcfactormatsolverpackage="umfpack" );

    timers["solve"].second = timers["solve"].first.elapsed();
    LOG(INFO) << "[timer] solve: " << timers["solve"].second << "\n";


    for ( size_type i = 0 ; i < u1.size(); ++ i )
        u1.set( i, ( *UbB )( i ) );

    for ( size_type i = 0 ; i < mu1.size(); ++ i )
        mu1.set( i, ( *UbB )( u1.size()+i ) );

    for ( size_type i = 0 ; i < mu.size(); ++ i )
        mu.set( i, ( *UbB )( u1.size()+mu1.size()+i ) );

    for ( size_type i = 0 ; i < mu2.size(); ++ i )
        mu2.set( i, ( *UbB )( u1.size()+mu1.size()+mu.size()+i ) );

    for ( size_type i = 0 ; i < u2.size(); ++ i )
        u2.set( i, ( *UbB )( u1.size()+mu1.size()+mu.size()+mu2.size()+i ) );

    // compute errors
    double L2error12 =integrate( elements( mesh1 ),( idv( u1 )-g )*( idv( u1 )-g ) ).evaluate()( 0,0 );
    double L2error1 =   math::sqrt( L2error12 );

    double L2error22 =integrate( elements( mesh2 ),( idv( u2 )-g )*( idv( u2 )-g ) ).evaluate()( 0,0 );
    double L2error2 =   math::sqrt( L2error22 );

    double semi_H1error1 =integrate( elements( mesh1 ),
                                     ( gradv( u1 )-gradg )*trans( ( gradv( u1 )-gradg ) ) ).evaluate()( 0,0 );

    double semi_H1error2 =integrate( elements( mesh2 ),
                                     ( gradv( u2 )-gradg )*trans( ( gradv( u2 )-gradg ) ) ).evaluate()( 0,0 );

    double H1error1 = math::sqrt( L2error12 + semi_H1error1 );

    double H1error2 = math::sqrt( L2error22 + semi_H1error2 );

    double global_error = math::sqrt( L2error12 + L2error22 + semi_H1error1 + semi_H1error2 );

    std::cout << "----------L2 errors---------- \n" ;
    std::cout << "||u1_error||_L2=" << L2error1 << "\n";
    std::cout << "||u2_error||_L2=" << L2error2 << "\n";
    std::cout << "----------H1 errors---------- \n" ;
    std::cout << "||u1_error||_H1=" << H1error1 << "\n";
    std::cout << "||u2_error||_H1=" << H1error2 << "\n";
    std::cout << "||u_error||_H1=" << global_error << "\n";

    LOG(INFO) << "----------L2 errors---------- \n" ;
    LOG(INFO) << "||u1_error||_L2=" << L2error1 << "\n";
    LOG(INFO) << "||u2_error||_L2=" << L2error2 << "\n";
    LOG(INFO) << "----------H1 errors---------- \n" ;
    LOG(INFO) << "||u1_error||_H1=" << H1error1 << "\n";
    LOG(INFO) << "||u2_error||_H1=" << H1error2 << "\n";
    LOG(INFO) << "||u_error||_H1=" << global_error << "\n";

    this->exportResults( u1,u2,mu1,mu2 );

    auto interface_exporter = trace_export_type::New( this->vm(),
                                                      ( boost::format( "interface_%1%-%2%" )
                                                        % this->about().appName()
                                                        % Dim ).str() ) ;

    interface_exporter->step( 0 )->setMesh( interface_mesh );
    interface_exporter->step( 0 )->add( "mu", mu );
    interface_exporter->save();

} // ThreeFieldsLaplacian::run

int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );
    Application app;
    app.add( new ThreeFieldsLaplacian<2,2,2,3>() );

    app.run();
}
