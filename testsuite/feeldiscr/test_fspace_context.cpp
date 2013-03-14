/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-03-14

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
   \file test_fspace_context.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2012-03-14
*/

#define BOOST_TEST_MODULE test_fspace_context
#include <testsuite/testsuite.hpp>

#include <fstream>

#include <boost/tuple/tuple.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <feel/feelcore/serialization.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/options.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;


inline
po::options_description
makeOptions()
{
    po::options_description testfspacecontext( "TestElementSerialize options" );
    testfspacecontext.add_options()
        ( "hsize", po::value<double>()->default_value( 0.2 ), "mesh size" )
        ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)" )
        ;
    return testfspacecontext.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_fspace_context" ,
                     "test_fspace_context" ,
                     "0.2",
                     "nD(n=1,2,3) test context of functionspace",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}


template<int Dim>
class TestFspaceContext
    :
public Simget
{
    typedef Simget super;
public:

    static const uint16_type Order = 1;

    typedef double value_type;
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef bases<Lagrange<Order,Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef space_type::element_type element_type;

    typedef Eigen::Matrix<double, nDim, 1> node_type;
    typedef space_type::Context ctx_type;

    /**
     * Constructor
     */
    TestFspaceContext( bool rebuild_database )
        :
        super(),
        M_meshSize( this->vm()["hsize"].template as<double>() ),
        M_shape( this->vm()["shape"].template as<std::string>() ),
        M_rebuild_database( rebuild_database )
    {
    }

    void run();
    void run( const double* X, unsigned long P, double* Y, unsigned long N );

    void saveDB();
    void loadDB();
    fs::path dbPath();

    friend class boost::serialization::access;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;
    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;
    bool existDB();
    void setRebuildDatabase( bool b);

private:

    double M_meshSize;
    std::string M_shape;
    bool M_rebuild_database;
    space_type M_Xh ;
    ctx_type M_ctx;
}; // TestFspaceContext


template<int Dim>
void
TestFspaceContext<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute TestFspaceContext<" << Dim << ">\n";
    std::vector<double> X( 2 );
    X[0] = M_meshSize;

    if ( M_shape == "hypercube" )
        X[1] = 1;

    else // default is simplex
        X[1] = 0;

    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}

template<int Dim>
void
TestFspaceContext<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) M_shape = "simplex";

    if ( X[1] == 1 ) M_shape = "hypercube";

#if 1
    if ( !this->vm().count( "nochdir" ) )

        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/h_%4%/" )
                                       % this->about().appName()
                                       % M_shape
                                       % Dim
                                       % M_meshSize );
#endif


    auto mesh = mesh_type::New();

    auto is_mesh_loaded = mesh->load( _name="mymesh",_path=".",_type="text" );

    if ( ! is_mesh_loaded )
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc=domain( _name=( boost::format( "%1%-%2%-%3%" ) % M_shape % Dim % 1 ).str() ,
                                             _usenames=true,
                                             _shape=M_shape,
                                             _dim=Dim,
                                             _h=X[0],
                                             _xmin=0.,
                                             _xmax=1.,
                                             _ymin=0.,
                                             _ymax=1.,
                                             _zmin=0.,
                                             _zmax=1. ) );
        mesh->save( _name="mymesh",_path=".",_type="text" );
    }


    if( ! existDB() || M_rebuild_database )
    {
        std::cout<<"no database..."<<std::endl;
        M_Xh = space_type::New( mesh );
        //std::cout<<"DB will be created "<<std::endl;
        //this->saveDB();
        //std::cout<<"DB successfully created "<<std::endl;
    }
    else
    {
        throw std::logic_error( "[TestFspaceContext] ERROR no database for now " );
        //std::cout<<"DB will be loaded"<<std::endl;
        //this->loadDB();
        //std::cout<<"DB successfully loaded"<<std::endl;
    }

    //expression we want to evaluate
    auto x = Px();
    auto y = Py();
    auto theta = 2*atan(Py()/(Px()+sqrt(Px()*Px()+Py()*Py())));
    auto r= sqrt(Px()*Px()+Py()*Py());

    //projection on the mesh
    auto px = vf::project( M_Xh , elements(mesh), x );
    auto py = vf::project( M_Xh , elements(mesh), y );
    auto ptheta = vf::project( M_Xh , elements(mesh), theta );
    auto pr = vf::project( M_Xh , elements(mesh), r );


    //nodes where we want evaluate expressions x,y,r,theta

#if DIM==1
    node_type t1; /*x*/ t1(0)=0.1;
    M_ctx.add( t1 );
    node_type t2; /*x*/ t2(0)=0.2;
    M_ctx.add( t2 );
    node_type t3; /*x*/ t3(0)=1;
    M_ctx.add( t3 );
#elif DIM==2
    node_type t1; /*x*/ t1(0)=0.1; /*y*/ t1(1)=0.2;
    M_ctx.add( t1 );
    node_type t2; /*x*/ t2(0)=0.1; /*y*/ t2(1)=0.8;
    M_ctx.add( t2 );
    node_type t3; /*x*/ t3(0)=1; /*y*/ t3(1)=1;
    M_ctx.add( t3 );
#elif DIM==3
    node_type t1; /*x*/ t1(0)=0.1; /*y*/ t1(1)=0.2; /*z*/ t1(2)=1;
    M_ctx.add( t1 );
    node_type t2; /*x*/ t2(0)=0.1; /*y*/ t2(1)=0.8; /*z*/ t2(2)=0.4;
    M_ctx.add( t2 );
    node_type t3; /*x*/ t3(0)=1; /*y*/ t3(1)=1; /*z*/ t2(2)=0.6;
    M_ctx.add( t3 );
#else
    throw std::logic_error("ERROR with the dimension ( dim > 3 ) " );
#endif

    //evaluation on all nodes via functionspace and M_ctx
    auto px_evaluate = px->evaluate( M_ctx );
    auto py_evaluate = py->evaluate( M_ctx );
    auto ptheta_evaluate = ptheta->evaluate( M_ctx );
    auto pr_evaluate = pr->evaluate( M_ctx );

    //true expressions (for verification)
    double x_t1=t1(0);
    double x_t2=t2(0);
    double x_t3=t3(0);
    double y_t1=t1(1);
    double y_t2=t2(1);
    double y_t3=t3(1);
    double r_t1 = sqrt( x_t1*x_t1 + y_t1*y_t1 );
    double r_t2 = sqrt( x_t2*x_t2 + y_t2*y_t2 );
    double r_t3 = sqrt( x_t3*x_t3 + y_t3*y_t3 );
    double theta_t1 = 2*atan(y_t1 / (x_t1 + r_t1 ) );
    double theta_t2 = 2*atan(y_t2 / (x_t2 + r_t2 ) );
    double theta_t3 = 2*atan(y_t3 / (x_t3 + r_t3 ) );

    //store true expressions in a vectore
    std::vector<double> solution_x, solution_y, solution_theta, solution_r;
    solution_x.resize(M_ctx.nPoints());
    solution_y.resize(M_ctx.nPoints());
    solution_theta.resize(M_ctx.nPoints());
    solution_r.resize(M_ctx.nPoints());

    //fill vectors solution
    solution_x[0] = x_t1;
    solution_x[1] = x_t2;
    solution_x[2] = x_t3;
    solution_y[0] = y_t1;
    solution_y[1] = y_t2;
    solution_y[2] = y_t3;
    solution_theta[0] = theta_t1;
    solution_theta[1] = theta_t2;
    solution_theta[2] = theta_t3;
    solution_r[0] = r_t1;
    solution_r[1] = r_t2;
    solution_r[2] = r_t3;

    //verification step
    for( int i=0; i<M_ctx.nPoints(); i++)
     {
         //check for expression x
         double evaluation_x = px_evaluate( i );
         double evaluation_x_node = px->evaluate(M_ctx , i);
         BOOST_CHECK_CLOSE( evaluation_x, evaluation_x_node, 1e-13 );
         BOOST_CHECK_CLOSE( evaluation_x, solution_x[i], 1e-13 );

         //check for expression y
         double evaluation_y = py_evaluate( i );
         double evaluation_y_node = py->evaluate(M_ctx , i);
         BOOST_CHECK_CLOSE( evaluation_y, evaluation_y_node, 1e-13 );
         BOOST_CHECK_CLOSE( evaluation_y, solution_y[i], 1e-13 );

         //check for expression theta
         double evaluation_theta = ptheta_evaluate( i );
         double evaluation_theta_node = ptheta->evaluate(M_ctx , i);
         BOOST_CHECK_CLOSE( evaluation_theta, evaluation_theta_node, 1e-13 );
         BOOST_CHECK_CLOSE( evaluation_theta, solution_theta[i], 1e-13 );

         //check for expression r
         double evaluation_r = pr_evaluate( i );
         double evaluation_r_node = pr->evaluate(M_ctx , i);
         BOOST_CHECK_CLOSE( evaluation_r, evaluation_r_node, 1e-13 );
         BOOST_CHECK_CLOSE( evaluation_r, solution_r[i], 1e-13 );
     }


} // TestFspaceContext::run

template< int Dim>
bool
TestFspaceContext<Dim>::existDB()
{
    return fs::exists( this->dbPath() / "database" );
}

template< int Dim>
template<class Archive>
void
TestFspaceContext<Dim>::save( Archive & ar, const unsigned int version ) const
{
}

template< int Dim>
template<class Archive>
void
TestFspaceContext<Dim>::load( Archive & ar, const unsigned int version )
{
}


template< int Dim>
fs::path
TestFspaceContext<Dim>::dbPath()
{
    std::string localpath = ( boost::format( "%1%/testsuite/feeldiscr/%2%/%3%-%4%/h_%5%/db/" )
                              % Feel::Environment::rootRepository()
                              % this->about().appName()
                              % M_shape
                              % Dim
                              % M_meshSize).str();

    fs::path rep_path = localpath;

    return rep_path ;
}


template< int Dim >
void
TestFspaceContext<Dim>::saveDB()
{

    fs::path path = this->dbPath();
    fs::create_directories( path );

    fs::ofstream ofs ( this->dbPath() / "database" );

    if ( ofs )
    {
        boost::archive::text_oarchive oa( ofs );
        oa << *this;
    }
    else
        throw std::logic_error(" can't use ofs in saveDB() " );
}

template< int Dim >
void
TestFspaceContext<Dim>::loadDB()
{

    fs::path db = this->dbPath();

    fs::ifstream ifs( db / "database" );

    if ( ifs )
    {
        boost::archive::text_iarchive ia( ifs );
        ia >> *this;
    }
}


template < int Dim >
void
TestFspaceContext<Dim>::setRebuildDatabase( bool b )
{
    M_rebuild_database = b;
}

/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_SUITE( fspace_context )

BOOST_AUTO_TEST_CASE( MyElementSerializeCase )
{
    Application app;

    //the first one :  create database ( if doesn't exist )
    //app.add( new TestFspaceContext<1>( true ) );
    //the second one : load database
    //app.add( new TestFspaceContext<1>(false ) );

    app.add( new TestFspaceContext<2>(true ) );
    //app.add( new TestFspaceContext<2>(false ) );

    //app.add( new TestFspaceContext<3>(true ) );
    //app.add( new TestFspaceContext<3>(false ) );

    app.run();

}
BOOST_AUTO_TEST_SUITE_END()


