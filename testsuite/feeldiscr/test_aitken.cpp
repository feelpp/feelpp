/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Abdoulaye Samake <abdoulaye.samake1@ujf-grenoble.fr>
   Date: 2011-04-14

   Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

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

#define BOOST_TEST_MODULE test_element_component
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/aitken.hpp>

using namespace Feel;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description laplacianoptions( "Aitken testsuite options" );
    laplacianoptions.add_options()
    ( "relaxmethod", po::value<int>()->default_value( 0 ), "use DD (=0) or DN (=1) method" )
    ( "additive", po::value<bool>()->default_value( false ), "use relax_aik additive method" )
    ( "x1max", Feel::po::value<double>()->default_value( 1.25 ),  " x1max for first subdomain" )
    ( "x2min", Feel::po::value<double>()->default_value( 1.25 ), " x2min for second subdomain" )
    ( "tol", Feel::po::value<double>()->default_value( 1e-06 ),  " tolerance " )
    ( "imax", Feel::po::value<double>()->default_value( 50 ), " maximum number of iteration" )
    ( "theta", Feel::po::value<double>()->default_value( 1. ), " relaxation parameter" )
    ;
    return laplacianoptions.add( Feel::feel_options() );
}

enum DDMethod
{
    // Dirichlet-Dirichlet
    DD = 0,
    // Dirichlet-Neumann
    DN = 1
};


template<typename SpaceType>
class LocalProblem
{
public :
    using space_ptrtype = std::shared_ptr<SpaceType>;
    LocalProblem( space_ptrtype Xh )
        :
        M_backend(  backend_type::build( soption( _name="backend" ) ) ),
        M_A( M_backend->newMatrix( _test=Xh, _trial=Xh ) ),
        M_B( M_backend->newVector( Xh ) )
        {}
    template<typename ElementType, typename RhsExpr, typename DirichletExpr, typename InterfaceExpr>
    void
    solve( ElementType& u,
           std::vector<int> const& dirichletFlags, DirichletExpr const& gD,
           RhsExpr const& f,
           std::vector<int> const& interfaceFlags, InterfaceExpr const& w,
           DDMethod choice )
        {
            auto Xh = u.functionSpace();
            auto mesh = Xh->mesh();
            auto const& v = u;

            form1( _test=Xh,_vector=M_B ) =
                integrate( _range=elements( mesh ), _expr=f*id( v ) );

            if ( choice == DDMethod::DN )
            {
                form1( _test=Xh,_vector=M_B ) +=
                    integrate( _range=markedfaces( mesh, interfaceFlags ), _expr=w*id( v ) );
            }

            form2( _test=Xh, _trial=Xh, _matrix=M_A ) =
                integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) );

            form2( _test=Xh, _trial=Xh, _matrix=M_A ) +=
                on( _range=markedfaces( mesh, dirichletFlags ), _element=u, _rhs=M_B, _expr=gD );

            if ( choice == DDMethod::DD )
            {
                form2( _test=Xh, _trial=Xh, _matrix=M_A ) +=
                    on( _range=markedfaces( mesh, interfaceFlags ) , _element=u, _rhs=M_B, _expr=w );
            }

            M_backend->solve( _matrix=M_A, _solution=u, _rhs=M_B );

        }
private :
    typedef Backend<double> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
    backend_ptrtype M_backend;
    backend_type::sparse_matrix_ptrtype M_A;
    backend_type::vector_ptrtype M_B;
};

template<typename MeshType>
void
runTestAitken()
{
    using mesh_type = MeshType;
    using mesh_ptrtype = std::shared_ptr<MeshType>;
    constexpr uint16_type Dim = mesh_type::nDim;
    std::string shape = "hypercube";

    uint32_type k0 = 0;
    uint32_type k1 = 1;

    double x1max = doption(_name="x1max");
    double x2min =  doption(_name="x2min");
    double tol =  doption(_name="tol");
    double imax =  doption(_name="imax");
    double theta =  doption(_name="theta");

    mesh_ptrtype mesh1 = createGMSHMesh( _mesh=new mesh_type,
                                         _desc=domain( _name=( boost::format( "%1%-%2%-%3%-%4%" ) % shape % Dim % 1 % k0 ).str() ,
                                                       _addmidpoint=false,
                                                       _usenames=false,
                                                       _shape=shape,
                                                       _dim=Dim,
                                                       //_h=X[0],
                                                       _xmin=0.,
                                                       _xmax=x1max,
                                                       _ymin=0.,
                                                       _ymax=2.,
                                                       _zmin=0.,
                                                       _zmax=2.
                                                       ) );

    mesh_ptrtype mesh2 = createGMSHMesh( _mesh=new mesh_type,
                                         _desc=domain( _name=( boost::format( "%1%-%2%-%3%-%4%" ) % shape % Dim % 1 % k1 ).str() ,
                                                       _addmidpoint=false,
                                                       _usenames=false,
                                                       _shape=shape,
                                                       _dim=Dim,
                                                       //_h=X[0],
                                                       _xmin=x2min,
                                                       _xmax=2.,
                                                       _ymin=0.,
                                                       _ymax=2.,
                                                       _zmin=0.,
                                                       _zmax=2.
                                                       ) );

    // flags for dirichlet boundary conditions
    std::vector<int> dirichletFlags1;
    std::vector<int> dirichletFlags2;

    // flags for interface conditions
    std::vector<int> interfaceFlags1;
    std::vector<int> interfaceFlags2;


    if ( Dim == 1 )
    {
        dirichletFlags1 = { 1 };
        dirichletFlags2 = { 3 };
        interfaceFlags1 = { 3 };
        interfaceFlags2 = { 1 };
    }
    else if ( Dim == 2 )
    {
        dirichletFlags1 = { 1,2,4 };
        dirichletFlags2 = { 2,3,4 };
        interfaceFlags1 = { 3 };
        interfaceFlags2 = { 1 };
    }
    else if ( Dim == 3 )
    {
        dirichletFlags1 = { 6,15,19,23,28 };
        dirichletFlags2 = { 6,15,23,27,28 };
        interfaceFlags1 = { 27 };
        interfaceFlags2 = { 19 };
    }

    // The function space and some associated elements(functions) are then defined
    auto Xh1 = Pch<2>( mesh1 );
    auto Xh2 = Pch<2>( mesh2 );
    auto u1 = Xh1->element();
    auto u2 = Xh2->element();
    auto u1old = Xh1->element();
    auto uu = Xh1->element();
    auto lambda = Xh2->element();
    auto residual = Xh2->element();

    LocalProblem<decay_type<decltype(Xh1)>> lp1( Xh1 );
    LocalProblem<decay_type<decltype(Xh2)>> lp2( Xh2 );

    AitkenType relaxmethod = ( AitkenType ) ioption(_name="relaxmethod");
    auto aitkenRelax =  aitken( _space= Xh2,
                                _type=( relaxmethod == 0 ) ? "standard" : "method1",
                                _initial_theta=theta,
                                _tolerance=tol );
    aitkenRelax.initialize( _residual=residual, _currentElt=lambda );

    double pi = M_PI;

    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    auto f = pi*pi*Dim*g;

    auto gradg = trans( +pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() )*unitX()+
                        -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() )*unitY()+
                        -pi*sin( pi*Px() )*cos( pi*Py() )*sin( pi*Pz() )*unitZ() );


    bool additive = boption(_name="additive");

    double err1 = 1.;
    double err2 = 1.;

    double error1 = 2.;
    double error2 = 2.;

    auto Ih12 = opInterpolation( _domainSpace=Xh1, _imageSpace=Xh2, _range=markedfaces( Xh2->mesh(), interfaceFlags2[0] ) );
    auto Ih21 = opInterpolation( _domainSpace=Xh2, _imageSpace=Xh1, _range=markedfaces( Xh1->mesh(), interfaceFlags1[0] ) );

    aitkenRelax.restart();

    while ( !aitkenRelax.isFinished() &&  aitkenRelax.nIterations() <= imax )
    {
        BOOST_TEST_MESSAGE( "============================================================" );
        BOOST_TEST_MESSAGE( "iteration  : " + std::to_string( aitkenRelax.nIterations() ) );
        BOOST_TEST_MESSAGE( "L2erroru1  : " + std::to_string( err1 ) );
        BOOST_TEST_MESSAGE( "L2erroru2  : " + std::to_string( err2 ) );
        BOOST_TEST_MESSAGE( "H1erroru1  : " + std::to_string( error1 ) );
        BOOST_TEST_MESSAGE( "H1erroru2  : " + std::to_string( error2 ) );

        u1old = u1;
        if ( additive )
        {
            if ( aitkenRelax.nIterations()==1 )
                BOOST_TEST_MESSAGE( "test_aitken additive method" );
        }
        else
        {
            if ( aitkenRelax.nIterations()==1 )
                BOOST_TEST_MESSAGE( "test_aiken multiplicative method" );
        }

        Ih21->apply( _residual=u2, _currentElt=uu );

        lp1.solve( u1,
                   dirichletFlags1, /*dirichlet*/g,
                   /*rhs*/f,
                   /**/interfaceFlags1,idv( uu ),
                   DDMethod::DD );

        if ( !additive )
            u1old = u1;

        lambda = u2;

        lp2.solve( u2,
                   dirichletFlags2, /*dirichlet*/g,
                   /*rhs*/f,
                   /**/ interfaceFlags2,gradv( u1old )*vf::N(),
                   DDMethod::DN );

        residual = u2-lambda;

        u2 = aitkenRelax.apply( _residual=residual, _currentElt=u2  );

        aitkenRelax.printInfo();

        aitkenRelax.saveConvergenceHistory( (boost::format( "history_%1%d.dat" )%Dim).str() );

        ++aitkenRelax;

        double L2error1 = integrate( _range=elements( mesh1 ),
                                    _expr=( idv( u1 )-g )*( idv( u1 )-g ) ).evaluate()( 0,0 );
        err1 = math::sqrt( L2error1 );

        double  L2error2 = integrate( _range=elements( mesh2 ),
                                      _expr=( idv( u2 )-g )*( idv( u2 )-g ) ).evaluate()( 0,0 );
        err2 = math::sqrt( L2error2 );

        double semi_H1error1 = 2.;
        double semi_H1error2 = 2.;


        semi_H1error1 =integrate( _range=elements( mesh1 ),
                                  _expr=( gradv( u1 )-gradg )*trans( ( gradv( u1 )-gradg ) ) ).evaluate()( 0,0 );
        semi_H1error2 =integrate( _range=elements( mesh2 ),
                                  _expr=( gradv( u2 )-gradg )*trans( ( gradv( u2 )-gradg ) ) ).evaluate()( 0,0 );

        error1 = math::sqrt( L2error1 + semi_H1error1 );
        error2 = math::sqrt( L2error2 + semi_H1error2 );

    }; // iteration loop


    BOOST_CHECK( aitkenRelax.nIterations() > 2 );
    BOOST_CHECK( aitkenRelax.nIterations() < imax );
#if 0
    auto exporter1 = exporter( _mesh=mesh1, _name=( boost::format( "%1%-%2%-%3%-%4%" )
                                                   % this->about().appName() % shape
                                                   % k0 % Dim ).str() );
    auto exporter2 = exporter( _mesh=mesh2, _name=( boost::format( "%1%-%2%-%3%-%4%" )
                                                    % this->about().appName() % shape
                                                    % k1 % Dim ).str() );
    if ( exporter1->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";
        exporter1->step( 0 )->add( "u1", u1 );
        exporter1->save();
        exporter2->step( 1 )->add( "u2", u2 );
        exporter2->save();
        LOG(INFO) << "exportResults done\n";
    }
#endif

}


FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault("test_aitken"), makeOptions() )

BOOST_AUTO_TEST_SUITE( test_aitken )

BOOST_AUTO_TEST_CASE( test_1d )
{
    runTestAitken<Mesh<Simplex<1>>>();
 }
BOOST_AUTO_TEST_CASE( test_2d )
{
    runTestAitken<Mesh<Simplex<2>>>();
}

BOOST_AUTO_TEST_SUITE_END()
