/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2012-05-23

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file test_eim.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-05-23
 */
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE eim testsuite

#if defined(USE_BOOST_TEST)
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#include <boost/test/floating_point_comparison.hpp>
#endif

#include <boost/timer.hpp>
#include <boost/smart_ptr/enable_shared_from_this.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelcrb/eim.hpp>

#define FEELAPP( argc, argv, about, options )                           \
    Feel::Application app( argc, argv, about, options );                \
    if ( app.vm().count( "help" ) )                                     \
    {                                                                   \
        std::cout << app.optionsDescription() << "\n";                  \
    }


namespace Feel
{
inline
po::options_description
makeOptions()
{
    po::options_description simgetoptions( "test_eim options" );
    simgetoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ;
    return simgetoptions.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_simget" ,
                     "test_simget" ,
                     "0.1",
                     "SimGet tests",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "" );
    return about;

}

/**
 *
 */
class model:
    public Simget,
    public boost::enable_shared_from_this<model>
{
public:
    typedef Mesh<Simplex<2> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type,bases<Lagrange<1> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;
    typedef space_type::element_type element_type;


    typedef ParameterSpace<1> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;
    typedef parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef parameterspace_type::sampling_type sampling_type;
    typedef parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef EIMFunctionBase<space_type, parameterspace_type> fun_type;
    typedef boost::shared_ptr<fun_type> fun_ptrtype;
    typedef std::vector<fun_ptrtype> funs_type;

    model( po::variables_map const& vm, AboutData const& about )
        :
        Simget( vm, about ),
        meshSize( vm["hsize"].as<double>() )
        {

            mesh = createGMSHMesh( _mesh=new mesh_type,
                                   _desc=domain( _name=( boost::format( "%1%-%2%" ) % "hypercube" % 2 ).str() ,
                                                 _usenames=true,
                                                 _shape="hypercube",
                                                 _dim=2,
                                                 _h=0.1 ) );

            Xh =  space_type::New( mesh );
            BOOST_CHECK( Xh );
            u = Xh->element();
            Dmu = parameterspace_type::New();
            BOOST_CHECK( Dmu );

            parameter_type mu_min( Dmu );
            mu_min << 0.2;
            Dmu->setMin( mu_min );
            parameter_type mu_max( Dmu );
            mu_max << 50;
            Dmu->setMax( mu_max );
            mu = Dmu->element();
            BOOST_CHECK_EQUAL( mu.parameterSpace(), Dmu );


            BOOST_TEST_MESSAGE( "Allocation done" );
            BOOST_TEST_MESSAGE( "pushing function to be empirically interpolated" );

            using namespace vf;
            //auto p = this->shared_from_this();
            //BOOST_CHECK( p );
            //BOOST_TEST_MESSAGE( "shared from this" );
            auto e = eim( _model=this,
                          _element=u,
                          _parameter=mu,
                          _expr=sin(cst_ref(mu(0)))*idv(u)*idv(u),
                          _name="q_1" );
            BOOST_TEST_MESSAGE( "create eim" );
            BOOST_CHECK( e );

            M_funs.push_back( e );
            BOOST_TEST_MESSAGE( "function to apply eim pushed" );

        }
    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return Dmu;
    }

    space_ptrtype functionSpace() { return Xh; }

    element_type solve( parameter_type const& mu )
        {
            using namespace vf;
            auto M = backend( _vm=this->vm() )->newMatrix( Xh, Xh );
            auto a = form2( _test=Xh, _trial=Xh, _matrix=M );
            auto F = backend( _vm=this->vm() )->newVector( Xh );
            auto rhs = form1( _test=Xh, _vector=F, _init=true );
            auto v = Xh->element();
            auto w = Xh->element();
            rhs = integrate( _range=elements( mesh ), _expr=id(v) );
            a = integrate( _range=elements( mesh ), _expr=mu(0)*gradt( w )*trans( grad( v ) ) + idt(w)*id(v) );
            a+=on( _range=boundaryfaces(mesh),
                   _element=w, _rhs=F, _expr=cst(0.) );
            backend(_vm=this->vm())->solve( _matrix=M, _solution=w, _rhs=F );
            return w;
        }
    void run()
        {
            auto exporter = Exporter<mesh_type>::New( this->vm(), this->about().appName() );
            exporter->step(0)->setMesh( mesh );
            auto S = Dmu->sampling();
            S->logEquidistribute(10);
            BOOST_FOREACH( auto fun, M_funs )
            {
                BOOST_FOREACH( auto p, *S )
                {
                    auto v = fun->operator()( p );
                    exporter->step(0)->add( (boost::format( "u(%1%)" ) % p(0) ).str(), v );

                }
            }
            exporter->save();

        }
    void run( const double*, long unsigned int, double*, long unsigned int ) {}
private:
    double meshSize;
    mesh_ptrtype mesh;
    space_ptrtype Xh;
    element_type u;
    parameterspace_ptrtype Dmu;
    parameter_type mu;

    funs_type M_funs;
};


//model circle
class model_circle:
    public Simget,
    public boost::enable_shared_from_this<model_circle>
{
public:
    typedef Mesh<Simplex<2> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type,bases<Lagrange<1> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;
    typedef space_type::element_type element_type;

    typedef ParameterSpace<4> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;
    typedef parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef parameterspace_type::sampling_type sampling_type;
    typedef parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef EIMFunctionBase<space_type, parameterspace_type> fun_type;
    typedef boost::shared_ptr<fun_type> fun_ptrtype;
    typedef std::vector<fun_ptrtype> funs_type;

    model_circle( po::variables_map const& vm, AboutData const& about )
        :
        Simget( vm, about ),
        meshSize( vm["hsize"].as<double>() )
        {

            mesh = createGMSHMesh( _mesh=new mesh_type,
                                   _desc=domain( _name=( boost::format( "%1%-%2%" ) % "hypercube" % 2 ).str() ,
                                                 _usenames=true,
                                                 _shape="hypercube",
                                                 _dim=2,
                                                 _h=0.1 ) );


            Xh =  space_type::New( mesh );
            BOOST_CHECK( Xh );
            u = Xh->element();
            Dmu = parameterspace_type::New();
            BOOST_CHECK( Dmu );

            //function :
            // \alpha ( x - cx )^2 + \beta ( y - cy )^2
            parameter_type mu_min( Dmu );
            mu_min << /*alpha*/0.01, /*beta*/ 0.01, /*cx*/ 1e-8, /*cy*/ 1e-8;
            Dmu->setMin( mu_min );
            parameter_type mu_max( Dmu );
            mu_max << /*alpha*/2, /*beta*/ 2, /*cx*/ 1, /*cy*/ 1;
            Dmu->setMax( mu_max );
            mu = Dmu->element();
            BOOST_CHECK_EQUAL( mu.parameterSpace(), Dmu );

            BOOST_TEST_MESSAGE( "Allocation done" );
            BOOST_TEST_MESSAGE( "pushing function to be empirically interpolated" );

            using namespace vf;
            //auto p = this->shared_from_this();
            //BOOST_CHECK( p );
            //BOOST_TEST_MESSAGE( "shared from this" );
            auto e = eim( _model=this,
                          _element=u,
                          _parameter=mu,
                          _expr= cst_ref(mu(0)) *( Px() - cst_ref(mu(2)) )*( Px() - cst_ref(mu(2)) )+cst_ref(mu(1)) *( Py() - cst_ref(mu(3)) )*( Py() - cst_ref(mu(3)) ),
                          _name="q_1");
            BOOST_TEST_MESSAGE( "create eim" );
            BOOST_CHECK( e );

            M_funs.push_back( e );
            BOOST_TEST_MESSAGE( "function to apply eim pushed" );

        }
    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return Dmu;
    }

    space_ptrtype functionSpace() { return Xh; }

    element_type solve( parameter_type const& mu )
        {
            using namespace vf;
            auto w = Xh->element();
            return w;
        }
    void run()
        {
            auto exporter = Exporter<mesh_type>::New( this->vm(), "model_circle" );
            exporter->step(0)->setMesh( mesh );
            auto S = Dmu->sampling();
            S->logEquidistribute(10);
            BOOST_FOREACH( auto fun, M_funs )
            {
                BOOST_FOREACH( auto p, *S )
                //auto p = Dmu->element();
                //p(0)=1; p(1)=1; p(2)=0.5; p(3)=0.5;
                {
                    auto v = fun->operator()( p );
                    exporter->step(0)->add( (boost::format( "eim_circle(%1%-%2%-%3%-%4%)" ) %p(0) %p(1) %p(2) %p(3) ).str(), v );

                }
            }
            exporter->save();

        }
    void run( const double*, long unsigned int, double*, long unsigned int ) {}
private:
    double meshSize;
    mesh_ptrtype mesh;
    space_ptrtype Xh;
    element_type u;
    parameterspace_ptrtype Dmu;
    parameter_type mu;

    funs_type M_funs;
};

} // Feel

using namespace Feel;

BOOST_AUTO_TEST_SUITE( eimsuite )

BOOST_AUTO_TEST_CASE( test_eim1 )
{
    Environment env;


    FEELAPP( 0,//boost::unit_test::framework::master_test_suite().argc,
             boost::unit_test::framework::master_test_suite().argv,
             makeAbout(), makeOptions() );
    BOOST_CHECK( mpi::environment::initialized() );
    BOOST_TEST_MESSAGE( "adding simget" );
    app.add( new model( app.vm(), app.about() ) );
    app.add( new model_circle( app.vm(), app.about() ) );
    app.run();

    BOOST_TEST_MESSAGE( "test_eim1 done" );

}

BOOST_AUTO_TEST_SUITE_END()


