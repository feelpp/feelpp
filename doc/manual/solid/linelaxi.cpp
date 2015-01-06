/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-03-20

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
   \file beamaxi2D.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-03-20
 */
#include <feel/feel.hpp>
#include <feel/feeltiming/tic.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description linelaxioptions( "LinElAxi options" );
    linelaxioptions.add_options()
        ( "E", Feel::po::value<double>()->default_value( 21e5 ), "Young modulus" )
        ( "nu", Feel::po::value<double>()->default_value( 0.28 ), "Poisson coefficient" )
        ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
        ( "bctype", Feel::po::value<int>()->default_value( 1 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
        ( "bccoeff", Feel::po::value<double>()->default_value( 1.0e+5 ), "coeff for weak Dirichlet conditions" )
        ( "gr", Feel::po::value<double>()->default_value( 1.0 ), "component in r of the surfacic force" )
        ( "gz", Feel::po::value<double>()->default_value( 1.0 ), "component in z of the surfacic force" )
        ;
    return linelaxioptions.add( Feel::feel_options() ) ;
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "linelaxi" ,
                           "linelaxi" ,
                           "0.1",
                           "Elasticity axisym  on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2012 University Joseph Fourier Grenoble 1" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


namespace Feel
{
/**
   \page LinElasAxi Axisymmetric Linear Elasticity
   \author Christophe Prud'homme
   
   \see beamaxi2D.cpp
 */
template<int Order>
class LinElAxi
    :
public Application
{
    typedef Application super;
public:
    typedef double value_type;

    /*mesh*/
    typedef Simplex<2> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef Lagrange<Order, Scalar> basis_scalar_type;
    typedef bases<basis_scalar_type,basis_scalar_type> basis_type;
    typedef Lagrange<Order, Vectorial> basis_vector_type;
    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef FunctionSpace<mesh_type, bases<basis_vector_type> > disp_space_type;
    typedef boost::shared_ptr<disp_space_type> disp_space_ptrtype;
    typedef typename disp_space_type::element_type disp_element_type;


    /* export */
    typedef Exporter<mesh_type> export_type;


    LinElAxi()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize"  ) ),
        bcCoeff(  doption("bccoeff") ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
    {
        LOG(INFO) << "[LinElAxi] hsize = " << meshSize << "\n";
        LOG(INFO) << "[LinElAxi] bccoeff = " << bcCoeff << "\n";
        LOG(INFO) << "[LinElAxi] export = " << this->vm().count( "export" ) << "\n";

    }

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double ,element_type& u );

private:

    backend_ptrtype M_backend;

    double meshSize;
    double bcCoeff;

    boost::shared_ptr<export_type> exporter;

}; // LinElAxi


template<int Order>
void
LinElAxi<Order>::run()
{
    tic();

    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    this->changeRepository( boost::format( "doc/manual/solid/%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % doption("hsize")
                          );
    /*
     * First we create the mesh
     */
    tic();
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name="beamaxi",
                                                _shape="hypercube",
                                                _usenames=true,
                                                _xmin=0, _xmax=1,
                                                _ymin=0, _ymax=10,
                                                _h=meshSize ),
                                        _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    toc();


    /*
     * The function space and some associate elements are then defined
     */
    tic();
    space_ptrtype Xh = space_type::New( mesh );

    //Xh->dof()->showMe();
    element_type U( Xh, "u" );
    element_type V( Xh, "v" );

    auto ur = U.template element<0>();
    auto vr = V.template element<0>();

    auto uz = U.template element<1>();
    auto vz = V.template element<1>();

    toc();

    /*
     * Data associated with the simulation
     */
    const double tol = 1e-5;
    const double E = Environment::vm(_name="E").template as<double>();
    const double nu = Environment::vm(_name="nu").template as<double>();
    const double mu = E/( 2*( 1+nu ) );
    const double lambda = E*nu/( ( 1+nu )*( 1-2*nu ) );
    const double density = 50;
    //    const double gravity = -density*9.81;
    const double gravity = -1.0;
    LOG(INFO) << "lambda = " << lambda << "\n"
          << "mu     = " << mu << "\n"
          << "gravity= " << gravity << "\n";

    auto gr = doption("gr");
    auto gz = doption("gz");
    /*
     * Construction of the constant right hand side
     *
     * \f$ f = \int_\Omega g * v \f$ where \f$ g \f$ is a vector
     * directed in the \f$ z \f$ direction.
     */

    auto rhs = M_backend->newVector( Xh );
    auto lhs = M_backend->newMatrix( Xh, Xh );

    using namespace Feel::vf;
    form1( _test=Xh, _vector=rhs ) =
        integrate( _range=elements( mesh ), _expr=gravity*id( vz )*Px() );
    form1( _test=Xh, _vector=rhs ) +=
        integrate( _range=markedfaces( mesh,"Dirichlet" ),
                   _expr=gr*id( vr )*Px()+gz*id( vz )*Px() );

    form2( _test=Xh, _trial=Xh, _matrix=lhs ) =
        integrate( _range=elements( mesh ),
                   _expr=Px()*lambda*( dxt( ur )+idt( ur )/Px()+
                                       dyt( uz ) )*( dx( vr )+idt( vr )/Px() ) );
    form2( _test=Xh, _trial=Xh, _matrix=lhs ) +=
        integrate( _range=elements( mesh ),
                   _expr=Px()*mu*( 2*( dxt( ur )*dx( vr )+idt( ur )*id( vr )/( Px()*Px() ) )+
                                   dyt( ur )*dy( vr )+dxt( uz )*dy( vr ) ) );
    form2( _test=Xh, _trial=Xh, _matrix=lhs ) +=
        integrate( _range=elements( mesh ),
                   _expr=( Px()*mu*( dyt( ur )*dx( vz )+dxt( uz )*dx( vz )+2*dyt( uz )*dy( vz ) )+
                           Px()*lambda*( dxt( ur )+idt( ur )/Px()+dyt( uz ) )*dy( vz ) ) );


    form2( _test=Xh, _trial=Xh, _matrix=lhs ) +=
        on( _range=markedfaces( mesh,"Neumann" ), _element=uz, _rhs=rhs, _expr=cst( 0. ) );
    form2( _test=Xh, _trial=Xh, _matrix=lhs ) +=
        on( _range=markedfaces( mesh,"Neumann" ), _element=ur, _rhs=rhs, _expr=cst( 0. ) );
    M_backend->solve( _matrix=lhs, _solution=U, _rhs=rhs );


    this->exportResults( 0, U );
    toc();

} //run


template<int Order>
void
LinElAxi<Order>::exportResults( double time, element_type& U )
{
    tic();
    using namespace Feel::vf;
    auto mesh = U.functionSpace()->mesh();
    auto ur = U.template element<0>();
    auto uz = U.template element<1>();
    disp_space_ptrtype Dh = disp_space_type::New( mesh );
    auto disp = Dh->element();
    disp = vf::project( _range=elements( mesh ), _space=Dh, _expr=vec( idv( ur ),idv( uz ) ) );

    exporter->step( time )->setMesh( mesh );
    exporter->step( time )->add( "ur", ur );
    exporter->step( time )->add( "uz", uz );
    exporter->step( time )->add( "displacement", disp );

    exporter->save();
    LOG(INFO) << "[timer] exportResults(): " << toc( "export results", false ) << "\n";
} // LinElAxi::export


} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );
    tic();

    typedef Feel::LinElAxi<2> linelaxi_type;

    /* assertions handling */
    Feel::Assert::setLog( "linelaxi.assert" );

    /* define and run application */
    linelaxi_type linelaxi;
    linelaxi.run();
    toc();
}
