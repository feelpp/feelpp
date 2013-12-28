/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-12-28

  Copyright (C) 2011 - 2014 Feel++ Consortium

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
   \file test_bdf2.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-12-28
*/

#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_bdf2
#include <testsuite/testsuite.hpp>
#endif

#include <feel/feel.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description testbdf2( "bdf2 test options" );
    testbdf2.add_options()
        ( "hsize", po::value<double>()->default_value( 0.05 ), "mesh size" )
        ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)" )
        ;
    return testbdf2.add( Feel::feel_options().add( Feel::bdf_options( "test_bdf2" ) ) );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_bdf2" ,
                     "test_bdf2" ,
                     "0.2",
                     "nD(n=2,3) test bdf2",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}

template<int Dim, int Order>
class Test:
    public Simget,
    public boost::enable_shared_from_this< Test<Dim,Order> >
{

public :

    typedef Mesh<Simplex<Dim> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type,bases<Lagrange<Order> >, Periodicity <NoPeriodicity> > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef double value_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;
    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;


    Test()
        :
        Simget(),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() )
    {
    }

    void run()
    {
        auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                    _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                                  _usenames=true,
                                                  _shape=shape,
                                                  _order=Order,
                                                  _dim=Dim,
                                                  _h=meshSize ) );

        auto Xh = Pch<1>( mesh );
        auto u = Xh->element();
        auto v = Xh->element();
        auto solution = Xh->element();
        auto A = backend()->newMatrix( Xh , Xh );
        auto D = backend()->newMatrix( Xh , Xh );
        auto M = backend()->newMatrix( Xh , Xh );
        auto F = backend()->newVector( Xh );
        auto Rhs = backend()->newVector( Xh );
        double mu0=0.5;
        //stifness matrix
        form2( _test=Xh, _trial=Xh, _matrix=A)  = integrate( _range= elements( mesh ), _expr= gradt( u )*trans( grad( v ) ) );
        //we have a robin condition
        form2( _test=Xh, _trial=Xh, _matrix=A) += integrate( _range= markedfaces( mesh, "Dirichlet" ), _expr= mu0 * idt( u )*id( v ) );
        //mass matrix
        form2( _test=Xh, _trial=Xh, _matrix=M)  = integrate( _range= elements( mesh ), _expr= idt( u )* id( v ) );
        //Rhs
        form1( Xh , F )  = integrate( _range=markedfaces( mesh,"Dirichlet" ), _expr= mu0 * id( v ) ) ;

        A->close();
        F->close();
        M->close();

        double bdf_coeff;
        auto vec_bdf_poly = backend()->newVector( Xh );

        bdf_ptrtype mybdf;
        mybdf = bdf( _space=Xh, _vm=this->vm() , _name="mybdf" );
        for ( mybdf->start(solution);  mybdf->isFinished() == false; mybdf->next() )
        {
            bdf_coeff = mybdf->polyDerivCoefficient( 0 );

            auto bdf_poly = mybdf->polyDeriv();

            *D = *A;
            *Rhs = *F;
            D->addMatrix( bdf_coeff, M );
            *vec_bdf_poly = bdf_poly;
            Rhs->addVector( *vec_bdf_poly, *M );

            D->close();
            Rhs->close();

            backend()->solve( _matrix=D, _solution=solution, _rhs=Rhs);

            mybdf->shiftRight(solution);
        }

        //check that we obtain the same result in sequential or in parallel
        double energy = M->energy( solution, solution );
        if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
            std::cout<<"energy : "<<energy<<std::endl;
    }

private :
    double meshSize;
    std::string shape;
};




#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( bdf2 )

BOOST_AUTO_TEST_CASE( test_1 )
{
    //Test<2,1> test ( new Test<2,1>() );
    boost::shared_ptr<Test<2,1> > test ( new Test<2,1>() );
    test->run();
}

BOOST_AUTO_TEST_SUITE_END()
#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=feel_options() );
    boost::shared_ptr<Test<2,1> > test ( new Test<2,1>() );
    //Test<2,1> test ( new Test<2,1>() );
    test->run();
}
#endif





