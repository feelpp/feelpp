/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake@e.ujf-grenoble.fr>
       Date: 2011-08-10

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

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
   \file test.cpp
   \author Abdoulaye Samake<abdoulaye.samake@e.ujf-grenoble.fr>
   \date 2011-08-10
*/

#define BOOST_TEST_MODULE test_lift
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;


inline
po::options_description
makeOptions()
{
    po::options_description testliftoptions( "TestLift options" );
    testliftoptions.add_options()
        ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
        ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
        ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
        ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ( "proj", po::value<int>()->default_value( 0 ), "use weak Dirichlet condition" )
        ( "penaldir", Feel::po::value<double>()->default_value( 20 ),
          "penalisation parameter for the weak boundary Dirichlet formulation" )
        ;
    return testliftoptions.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_lift" ,
                     "test_lift" ,
                     "0.2",
                     "nD(n=1,2,3) TestLift on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier" );

    about.addAuthor( "Abdoulaye Samake", "developer", "abdoulaye.samake@imag.fr", "" );
    return about;

}


template<int Dim>
class TestLift
    :
public Simget
{
    typedef Simget super;
public:

    static const uint16_type Order = 2;
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef Simplex<Dim,1,Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef bases<Lagrange<Order,Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef typename space_type::element_type element_type;
    typedef Exporter<mesh_type> export_type;

    /**
     * Constructor
     */
    TestLift()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() )
    {
    }

    void run();
    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:

    backend_ptrtype M_backend;
    double meshSize;
    std::string shape;
}; // TestLift

template<int Dim> const uint16_type TestLift<Dim>::Order;

template<int Dim>
void
TestLift<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute TestLift<" << Dim << ">\n";
    std::vector<double> X( 2 );
    X[0] = meshSize;

    if ( shape == "hypercube" )
        X[1] = 1;

    else // default is simplex
        X[1] = 0;

    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template<int Dim>
void
TestLift<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) shape = "simplex";

    if ( X[1] == 1 ) shape = "hypercube";

    if ( !this->vm().count( "nochdir" ) )

        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/P%4%/h_%5%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name=( boost::format( "%1%-%2%-%3%" ) % shape % Dim % 1 ).str() ,
                                        _usenames=true,
                                        _shape=shape,
                                        _dim=Dim,
                                        _h=X[0],
                                        _xmin=0.,
                                        _xmax=1.,
                                        _ymin=0.,
                                        _ymax=1.,
                                        _zmin=0.,
                                        _zmax=1. ) );


    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();

    value_type pi = M_PI;
    auto g = sin( pi*Px()/4 )*cos( pi*Py()/3 )*cos( pi*Pz() );
    auto f = cst( 0. );

    bool weakdir = this->vm()["weakdir"].template as<int>();
    value_type penaldir = this->vm()["penaldir"].template as<double>();
    value_type nu = this->vm()["nu"].template as<double>();

    auto F =  M_backend->newVector( Xh ) ;

    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements( mesh ), f*id( v ) ) ;

    if ( this->comm().size() != 1 || weakdir )
        {
            form1( _test=Xh, _vector=F ) +=
                integrate( markedfaces( mesh,"Dirichlet" ),
                           g*( -grad( v )*vf::N()+penaldir*id( v )/hFace() ) );
        }
    F->close();

    auto D = M_backend->newMatrix( Xh, Xh ) ;

    form2( _test=Xh, _trial=Xh, _matrix=D, _init=true ) =
        integrate( elements( mesh ), nu*gradt( u )*trans( grad( v ) ) );

    if ( this->comm().size() != 1 || weakdir )
        {

            form2( _test=Xh, _trial=Xh, _matrix=D ) +=
                integrate( markedfaces( mesh,"Dirichlet" ),
                           -( gradt( u )*vf::N() )*id( v )
                           -( grad( v )*vf::N() )*idt( u )
                           +penaldir*id( v )*idt( u )/hFace() );
            D->close();
        }

    else
        {
            D->close();

            form2( _test=Xh, _trial=Xh, _matrix=D ) +=
                on( markedfaces( mesh, "Dirichlet" ), u, F, g );
        }

    M_backend->solve( _matrix=D, _solution=u, _rhs=F );

    DirichletType dir_type = (DirichletType)this->vm()["weakdir"].template as<int>();
    Feel::ProjectorType proj_type = (Feel::ProjectorType)this->vm()["proj"].template as<int>();

    /** create operators using different  interfaces */
    auto op_lift1 = projector( Xh, Xh, M_backend, proj_type, 0.01, penaldir, dir_type );
    auto op_lift2 = opLift( _domainSpace=Xh, _backend=M_backend, _type=dir_type, _penaldir=penaldir );
    auto op_projection = opProjection( _domainSpace=Xh, _imageSpace=Xh, _backend=M_backend, _type=proj_type, _penaldir=penaldir );

    /** apply operators on the function g */
    auto glift1 = op_lift1->project( _expr=trans( g ), _range=markedfaces( mesh,"Dirichlet" ) );
    auto glift2 = op_lift2->project( _expr=trans( g ), _range=markedfaces( mesh,"Dirichlet" ) );
    auto glift2a = ( *op_lift2 )( _expr= trans(g), _range=markedfaces( mesh,"Dirichlet" ) );

    /** project the function g for export */
    auto gproj =  vf::project( _space=Xh, _range=elements( mesh ), _expr=g );

    /** compute errors */
    double L2error2 =integrate( elements( mesh ),
                                ( idv( u )-idv( glift2 ) )*( idv( u )-idv( glift2 ) ) ).evaluate()( 0,0 );

    double semi_H1error2 =integrate( elements( mesh ),
                                     ( gradv( u )-gradv( glift2 ) )*trans( ( gradv( u )-gradv( glift2 ) ) ) ).evaluate()( 0,0 );

    std::cout << " -- ||u-glift||_L2  =" << math::sqrt( L2error2 ) << "\n";
    std::cout << " -- ||u-glift||_H1  =" << math::sqrt( L2error2+semi_H1error2 ) << "\n";

    /** export results */
    auto exporter( export_type::New( this->vm(),
                                     ( boost::format( "%1%-%2%-%3%" )
                                       % this->about().appName()
                                       % shape
                                       % Dim ).str() ) );

    if ( exporter->doExport() )
        {
            LOG(INFO) << "exportResults starts\n";
            exporter->step( 0 )->setMesh( mesh );
            exporter->step( 0 )->add( "u", u );
            exporter->step( 0 )->add( "glift", glift2 );
            exporter->step( 0 )->add( "g", gproj );
            exporter->save();
            LOG(INFO) << "exportResults done\n";
        }
} // TestLift::run


/**
 * main code
 */
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( lift )
BOOST_AUTO_TEST_CASE( MyLiftCase )
{
    Application app;

    //app.add( new TestLift<1>( app.vm(), app.about() ) );
    app.add( new TestLift<2>() );
    //app.add( new TestLift<3>( app.vm(), app.about() ) );
    app.run();
}

BOOST_AUTO_TEST_SUITE_END()
