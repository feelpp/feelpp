/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Author(s): Stephane Veys <stephane.veys@imag.fr>
  Author(s): Cecile Daversin <daversin@math.unistra.fr>
       Date: 2014-01-06

  Copyright (C) 2011 - 2014 Feel++ Consortium

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
   \file test_eim_grepl.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane veys <stephane.veys@imag.fr>
   \author Cecile Daversin <daversin@math.unistra.fr>
   \date 2014-01-06
 */

/**
  This test is from M.A.Grepl thesis, and appears also in the following paper
  EFFICIENT REDUCED-BASIS TREATMENT OF NONAFFINE
  AND NONLINEAR PARTIAL DIFFERENTIAL EQUATIONS
  authors :
  Martin A. Grepl, Yvon Maday, Ngoc C. Nguyen and Anthony T. Patera
  ESAIM: Mathematical Modelling and Numerical Analysis
 */
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE eim_composite testsuite

#include <testsuite/testsuite.hpp>

#include <boost/timer.hpp>
#include <boost/smart_ptr/enable_shared_from_this.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/unitsquare.hpp>
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
    po::options_description eimCompositeoptions( "test_eim_composite options" );
    eimCompositeoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ( "chrono-online-step" , po::value<bool>()->default_value( false ), "give access to computational time during online step if true" )
    ( "n-eval", po::value<int>()->default_value( 10 ), "number of evaluations" )
    ( "cvg-study" , po::value<bool>()->default_value( false ), "run a convergence study if true" )
    ;
    return eimCompositeoptions.add( eimOptions() ).add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_eim_composite" ,
                     "test_eim_composite" ,
                     "0.1",
                     "SimGet tests",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2012 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    return about;

}

/**
 *
 */
class EimModel:
        public boost::enable_shared_from_this<EimModel>
{
    typedef boost::enable_shared_from_this<EimModel> super;
public:
    typedef Mesh<Simplex<2> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    //typedef FunctionSpace<mesh_type,bases<Lagrange<1> > > space_type;
    typedef Lagrange<1> basis_type;
    typedef bases<basis_type,basis_type> prod_basis_type;
    typedef FunctionSpace<mesh_type, bases< Lagrange<1> > > u1_space_type;
    typedef FunctionSpace<mesh_type, prod_basis_type > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<u1_space_type> u1_space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;
    typedef space_type::element_type element_type;
    typedef space_type::element_ptrtype element_ptrtype;

    typedef u1_space_type::element_type u1_element_type;
    typedef u1_space_type::element_ptrtype u1_element_ptrtype;

    typedef ParameterSpace<2> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;
    typedef parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef parameterspace_type::sampling_type sampling_type;
    typedef parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef EIMFunctionBase<u1_space_type, space_type, parameterspace_type> fun_type;
    typedef boost::shared_ptr<fun_type> fun_ptrtype;
    typedef std::vector<fun_ptrtype> funs_type;

    typedef Eigen::VectorXd vectorN_type;

    EimModel()
        :
        super()
        {
        }
    void init()
        {
            mesh = unitSquare();

            Xh =  space_type::New( mesh );
            Xh1 = u1_space_type::New( mesh );
            LOG(INFO) << " nb dofs (composite): "<<Xh->nDof()<<"\n";
            LOG(INFO) << " nb dofs (no composite): "<<Xh1->nDof()<<"\n";
            BOOST_CHECK( Xh );
            u = Xh->element();
            std::cout << "u size = " << u.size() << std::endl;
            u1 = Xh1->element();
            u2 = Xh1->element();
            u1 = u.element<0>();
            u2 = u.element<1>();
            LOG(INFO) << "u1 size = " << u1.size() << std::endl;
            LOG(INFO) << "u2 size = " << u2.size() << std::endl;
            Dmu = parameterspace_type::New();
            BOOST_CHECK( Dmu );
            parameter_type mu_min( Dmu );
            mu_min << -1,-1;
            Dmu->setMin( mu_min );
            parameter_type mu_max( Dmu );
            mu_max << -0.01,-0.01;
            Dmu->setMax( mu_max );
            mu = Dmu->element();

            auto Pset = Dmu->sampling();
            //specify how many elements we take in each direction
            std::vector<int> N(2);
            //40 elements in each direction
            N[0]=40; N[1]=40;
            Pset->equidistributeProduct( N );
            BOOST_TEST_MESSAGE( "Allocation done" );
            BOOST_TEST_MESSAGE( "pushing function to be empirically interpolated" );

            using namespace vf;

            std::cout << "before eim" << std::endl;
            auto e = eim( _model=eim_no_solve(this->shared_from_this()),
                          _element=u,
                          _space=Xh1,
                          _parameter=mu,
                          _expr=1/sqrt( (Px()-cst_ref(mu(0)))*(Px()-cst_ref(mu(0))) + (Py()-cst_ref(mu(1)))*(Py()-cst_ref(mu(1))) ),
                          //_expr = cst_ref( mu(0) ),
                          _sampling=Pset,
                          _name="q" );
            BOOST_CHECK( e );
            M_funs.push_back( e );
        }
    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
        {
            return Dmu;
        }
    std::string modelName() const { return std::string("test_eim_composite" );}

    space_ptrtype functionSpace() { return Xh; }

    element_type solve( parameter_type const& mu )
        {
            element_type e = Xh->element();
            return e;
        }
    void run()
        {
            auto solution = Xh->elementPtr();
            auto e = exporter( _mesh=mesh, _name=Environment::about().appName() );
            auto S = Dmu->sampling();
            int n = option("n-eval").as<int>();
            LOG(INFO)<<"will compute "<<n<<" evaluations\n";
            bool chrono = option("chrono-online-step").as<bool>();
            bool cvg_study = option("cvg-study").as<bool>();
            S->equidistribute(n);
            int fun_number=0;
            std::vector<vectorN_type> time_vector( M_funs.size() );
            for( auto fun : M_funs )
            {
                time_vector[fun_number].resize( n );
                int mu_number=0;
                for( auto p : *S )
                {
                    *solution = solve(p);
                    int max = fun->mMax();

                    if( ! chrono )
                    {
                        LOG(INFO) << "evaluate model at p = " << p << "\n";
                        auto v = fun->operator()( p );
                        e->add( (boost::format( "%1%(%2%)" ) % fun->name() % p(0) ).str(), v );
                        LOG(INFO) << "evaluate eim interpolant at p = " << p << "\n";
                    }
                    if( cvg_study )
                    {
                        *solution = solve(p);
                        std::vector< std::string > all_file_name;
                        all_file_name.push_back( "EimConvergenceL2.dat");
                        all_file_name.push_back( "EimConvergenceL2estimated.dat");
                        all_file_name.push_back( "EimConvergenceL2ratio.dat");
                        all_file_name.push_back( "EimConvergenceLINF.dat");
                        all_file_name.push_back( "EimConvergenceLINFestimated.dat");
                        all_file_name.push_back( "EimConvergenceLINFratio.dat");
                        fun->studyConvergence( p , *solution , all_file_name );
                    }
                    boost::mpi::timer timer;
                    auto w = fun->interpolant( p );
                    double t=timer.elapsed();
                    e->add( (boost::format( "%1%-eim(%2%)" ) % fun->name() % p(0) ).str(), w );
                    time_vector[fun_number]( mu_number )=t;
                    mu_number++;
                }
                fun_number++;
            }
            e->save();

            //some statistics
            LOG(INFO)<<"Computational time during online step ( "<<n<<" evaluations )\n";
            if( option("eim.use-dimension-max-functions").as<bool>() )
                LOG(INFO)<<option("eim.dimension-max").as<int>()<<" eim basis functions were used\n";
            Eigen::MatrixXf::Index index;
            for(int expression=0; expression<time_vector.size(); expression++)
            {
                LOG(INFO)<<"expression number "<<expression<<"\n";
                double min = time_vector[expression].minCoeff(&index);
                double max = time_vector[expression].maxCoeff(&index);
                double mean = time_vector[expression].mean();
                LOG(INFO)<<"min : "<<min<<" - max : "<<max<<" - mean : "<<mean<<std::endl;
            }

        }
    void run( const double*, long unsigned int, double*, long unsigned int ) {}
private:
    double meshSize;
    mesh_ptrtype mesh;
    space_ptrtype Xh;
    u1_space_ptrtype Xh1;
    element_type u;
    u1_element_type u1,u2;
    parameterspace_ptrtype Dmu;
    parameter_type mu;

    funs_type M_funs;
};

} // Feel

using namespace Feel;

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( eimsuite_composite )

BOOST_AUTO_TEST_CASE( test_eim_composite )
{
    BOOST_TEST_MESSAGE( "test_eim_composite starts..." );

    boost::shared_ptr<EimModel> m( new EimModel);
    m->init();
    m->run();

    BOOST_TEST_MESSAGE( "test_eim1 done" );

}

BOOST_AUTO_TEST_SUITE_END()
