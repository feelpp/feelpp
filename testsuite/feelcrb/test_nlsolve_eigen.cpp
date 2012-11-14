/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-11-14

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-11-14
 */
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE testsuite

#include <testsuite/testsuite.hpp>

#include <boost/timer.hpp>
#include <boost/smart_ptr/enable_shared_from_this.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

#include <feel/feelalg/solvernonlinear.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>


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
    po::options_description simgetoptions( "test_nlsolve_eigen options" );
    simgetoptions.add_options()
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
                     "Copyright (c) 2012 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}

/**
 *
 */
class linear:
    public Simget,
    public boost::enable_shared_from_this<linear>
{
public:

    typedef Eigen::VectorXd vectorN_type;
    typedef Eigen::MatrixXd matrixN_type;
    typedef Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > map_dense_matrix_type;
    typedef Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > map_dense_vector_type;

    typedef linear self_type;

    linear()
        :
        Simget(),
        M_nlsolver( SolverNonLinear<double>::build( SOLVERS_PETSC, Environment::worldComm() ) )
        {
        }

        void updateResidual(const map_dense_vector_type& map_X, map_dense_vector_type& map_R )
        {
            int size=2;

            matrixN_type A(size,size);
            A(0,0) = 2; A(0,1) = 1 ;
            A(1,0) = 0; A(1,1) = 3 ;
            vectorN_type F(size);
            F(0)=10;
            F(1)=6;

            map_R = A*map_X - F;
        }
        void updateJacobian(const map_dense_vector_type& map_X, map_dense_matrix_type& map_J )
        {
            map_J(0,0) = 2; map_J(0,1) = 1 ;
            map_J(1,0) = 0; map_J(1,1) = 3 ;
        }


    void run()
        {
            M_nlsolver->map_dense_residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
            M_nlsolver->map_dense_jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

            //initial guess
            vectorN_type solution(2);
            solution<<1,5;

            matrixN_type J(2,2);
            vectorN_type R(2);

            BOOST_TEST_MESSAGE( "initialization done, now convert R,J,solution via map" );


            double *r_data = R.data();
            double *j_data = J.data();
            double *solution_data = solution.data();

            int size=2;
            Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_R ( r_data, size );
            Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_solution ( solution_data, size );
            Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , Eigen::Dynamic> > map_J ( j_data, size , size);

            updateResidual( map_solution, map_R );
            updateJacobian( map_solution, map_J );

            M_nlsolver->solve( map_J, map_solution, map_R, 1e-10, 1 );
            BOOST_TEST_MESSAGE( "system solved" );

            std::cout<<"solution of linear problem : ("<< solution(0)<<" , "<<solution(1)<<") with initial guess (1 , 5) \n"<<std::endl;
            BOOST_CHECK_CLOSE( solution(0), 4, 1e-10);
            BOOST_CHECK_CLOSE( solution(1), 2, 1e-10);

            BOOST_TEST_MESSAGE( "solution checked" );

        }
    void run( const double*, long unsigned int, double*, long unsigned int ) {}
private:

    boost::shared_ptr<SolverNonLinear<double> > M_nlsolver;
};



class NL22:
    public Simget,
    public boost::enable_shared_from_this<NL22>
{
public:

    typedef Eigen::VectorXd vectorN_type;
    typedef Eigen::MatrixXd matrixN_type;
    typedef Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > map_dense_matrix_type;
    typedef Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > map_dense_vector_type;

    typedef NL22 self_type;

    NL22( const vectorN_type& initial_guess , const vectorN_type & exact_solution)
        :
        Simget(),
        M_nlsolver( SolverNonLinear<double>::build( SOLVERS_PETSC, Environment::worldComm() ) ),
        M_initial_guess ( initial_guess ),
        M_exact_solution( exact_solution )
        {
        }

        void updateResidual(const map_dense_vector_type& map_X, map_dense_vector_type& map_R )
        {
            map_R(0) = map_X(0)*map_X(0) - 2*map_X(0)*map_X(1) - 2 ;
            map_R(1) = map_X(0) + map_X(1)*map_X(1) + 1;
        }
        void updateJacobian(const map_dense_vector_type& map_X, map_dense_matrix_type& map_J )
        {
            map_J(0,0) = 2*map_X(0)-2*map_X(1); map_J(0,1) = -2*map_X(0) ;
            map_J(1,0) = 1;                     map_J(1,1) = 2*map_X(1) ;
        }


    void run()
        {
            M_nlsolver->map_dense_residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
            M_nlsolver->map_dense_jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

            //initial guess
            vectorN_type solution(2);
            solution = M_initial_guess;

            matrixN_type J(2,2);
            vectorN_type R(2);

            BOOST_TEST_MESSAGE( "initialization done, now convert R,J,solution via map" );


            double *r_data = R.data();
            double *j_data = J.data();
            double *solution_data = solution.data();

            int size=2;
            Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_R ( r_data, size );
            Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_solution ( solution_data, size );
            Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , Eigen::Dynamic> > map_J ( j_data, size , size);

            updateResidual( map_solution, map_R );
            updateJacobian( map_solution, map_J );

            M_nlsolver->solve( map_J, map_solution, map_R, 1e-10, 10 );
            BOOST_TEST_MESSAGE( "system solved" );

            std::cout<<"solution of non linear problem : ("<< solution(0)<<" , "<<solution(1)<<") with initial guess ("<<M_initial_guess(0)<<" , "<<M_initial_guess(1)<<")"<<std::endl;
            double norm_solution = solution.norm();
            double norm_exact_solution = M_exact_solution.norm();
            BOOST_CHECK_CLOSE( norm_solution, norm_exact_solution, -10);
            BOOST_TEST_MESSAGE( "solution checked" );

        }
    void run( const double*, long unsigned int, double*, long unsigned int ) {}
private:

    boost::shared_ptr<SolverNonLinear<double> > M_nlsolver;
    vectorN_type M_initial_guess;
    vectorN_type M_exact_solution;
};


} // Feel

using namespace Feel;

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( nlsolveeigensuite )

BOOST_AUTO_TEST_CASE( test_nlsolve_eigen1 )
{

    typedef Eigen::VectorXd vectorN_type;

    BOOST_CHECK( mpi::environment::initialized() );
    Application app;
    BOOST_TEST_MESSAGE( "adding simget" );
    app.add( new linear );

    vectorN_type initial_guess(2);
    vectorN_type exact_solution(2);
    initial_guess<<-5,-2;
    exact_solution<<-3.93432,-1.71298;
    app.add( new NL22( initial_guess , exact_solution ) );

    initial_guess<<-0.9,2;
    exact_solution<<-1.1150879943748373 , 0.33924621563435658;
    app.add( new NL22( initial_guess , exact_solution ) );

    app.run();

    BOOST_TEST_MESSAGE( "test_nlsolve_eigen1 done" );

}

BOOST_AUTO_TEST_SUITE_END()
