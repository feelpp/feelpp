/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
#include <feel/feelcrb/crbplugin.hpp>
#include <unsteadyheat1d.hpp>

namespace Feel
{

po::options_description
makeUnsteadyHeat1DOptions()
{
    po::options_description unsteadyheat1doptions( "UnsteadyHeat1D options" );
    unsteadyheat1doptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.01 ), "mesh size" )
    ( "mu1", po::value<double>()->default_value( 0.2 ), "mu1" )
    ( "mu2", po::value<double>()->default_value( 0.2 ), "mu2" )
    ( "mu3", po::value<double>()->default_value( -1 ), "mu3" )
    ( "mu4", po::value<double>()->default_value( 0.1 ), "mu4" )
    ( "alpha" , po::value<double>()->default_value( 1 ), "temporal coefficient" )
    ;
    return unsteadyheat1doptions.add( bdf_options( "unsteadyHeat1d" ) ).add( ts_options( "unsteadyHeat1d" ) );
}
AboutData
makeUnsteadyHeat1DAbout( std::string const& str)
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "unsteady 1D Heat Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010,2011 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;
}



gmsh_ptrtype
createGeo( double meshSize  )
{
    std::ostringstream ostr;
    ostr << "P=" << meshSize << ";\n"
         << "Point(1) = {-0.5,0,0,P};\n"
         << "Point(2) = {-0.1,0,0,P};\n"
         << "Point(3) = {0,0,0,P};\n"
         << "Point(4) = {0.1,0,0.0,P};\n"
         << "Point(5) = {0.5,0,0.0,P};\n"
         << "\n"

         << "Line(1) = {1, 2};\n"
         << "Line(2) = {2, 3};\n"
         << "Line(3) = {3, 4};\n"
         << "Line(4) = {4, 5};\n"
         << "Line Loop(5) = {1,2,3,4};\n"
         << "Physical Line(\"k1_1\") = {1};\n"
         << "Physical Line(\"k1_2\") = {2};\n"
         << "Physical Line(\"k2_1\") = {3};\n"
         << "Physical Line(\"k2_2\") = {4};\n"
         << "Physical Point(\"left\") = {1};\n"
         << "Physical Point(\"right\") = {5};\n"
         << "\n";

    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( "unsteadyheat1d" );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}

UnsteadyHeat1D::UnsteadyHeat1D()
    :
    super_type( "heat1d" ),
    alpha( 1 )
{}

void
UnsteadyHeat1D::initModel()
{
    /*
     * First we create the mesh
     */
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=createGeo( doption("gmsh.hsize") ) );

    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    setFunctionSpaces(Xh);
    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << Xh->nLocalDof() << "\n";
        std::cout << "Number of dof " << Xh->nDof() << "\n";
    }


    M_bdf = bdf( _space=Xh, _name="unsteadyHeat1d" , _prefix="unsteadyHeat1d" );

    Dmu->setDimension( 4 );
    auto mu_min = Dmu->element();
    mu_min << 0.2, 0.2, 0.01, 0.1;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 50, 50, 5, 5;
    Dmu->setMax( mu_max );

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    assemble();
} // UnsteadyHeat1d::initModel


void
UnsteadyHeat1D::assemble()
{
    auto u = Xh->element();
    auto v = Xh->element();

    //lhs
    auto a0 = form2( _trial=Xh, _test=Xh);
    a0 = integrate( _range=elements( mesh ), _expr=0.1*( gradt( u )*trans( grad( v ) ) ) ) +
         integrate( _range=markedfaces( mesh,"right" ), _expr=idt( u )*id( v ) );
    this->addLhs( {a0 , "1"} );

    auto a1 = form2( _trial=Xh, _test=Xh);
    a1 = integrate( _range=markedelements( mesh,"k1_1" ),  _expr=gradt( u )*trans( grad( v ) )  );
    this->addLhs(  {a1 , "mu0"} );

    auto a2 = form2( _trial=Xh, _test=Xh);
    a2 = integrate( _range=markedelements( mesh,"k2_1" ),  _expr=gradt( u )*trans( grad( v ) )  );
    this->addLhs( {a2 , "mu1"} );

    //mass matrix
    auto m = form2( _trial=Xh, _test=Xh);
    m = integrate ( _range=elements( mesh ), _expr=alpha*idt( u )*id( v ) );
    this->addMass( {m , "1"} );

    //rhs
    auto f0 = form1( _test=Xh );
    auto f1 = form1( _test=Xh );
    f0 = integrate( _range=markedfaces( mesh,"left" ), _expr=id( v ) );
    this->addRhs( {f0 , "mu2"} );
    f1 =  integrate( _range=elements( mesh ), _expr=id( v ) );
    this->addRhs( {f1 , "mu3"} );

    //output
    auto out = form1( _test=Xh );
    out = integrate( _range=markedelements( mesh,"k1_2" ), _expr=id( v )/0.2 ) +
          integrate( _range=markedelements( mesh,"k2_1" ), _expr=id( v )/0.2 );
    this->addOutput( {out , "1" } );

    auto energy = form2( _trial=Xh, _test=Xh);
    energy = integrate( _range=elements( mesh ), _expr=0.1*( gradt( u )*trans( grad( v ) ) ) ) +
             integrate( _range=markedfaces( mesh,"right" ), _expr=idt( u )*id( v ) ) +
             integrate( _range=markedelements( mesh,"k1_1" ),  _expr=0.2 * gradt( u )*trans( grad( v ) ) )  +
             integrate( _range=markedelements( mesh,"k2_1" ),  _expr=0.2 * gradt( u )*trans( grad( v ) ) )  ;
    this->addEnergyMatrix( energy.matrixPtr() );

    auto mass = form2( _trial=Xh, _test=Xh);
    mass = integrate( _range=elements( mesh ), _expr=idt( u ) * id( v ) ) ;
    this->addMassMatrix( mass.matrixPtr() );
}


double
UnsteadyHeat1D::output( int output_index, parameter_type const& mu, element_type& u, bool need_to_solve )
{

    CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    auto v = Xh->element();

    double output=0;
    // right hand side (compliant)
    if ( output_index == 0 )
    {
        output  = integrate( _range=markedfaces( mesh,"left" ), _expr=mu(2)*idv( u ) ).evaluate()( 0,0 );
        output += integrate( _range=elements( mesh ), _expr=mu(3)*idv( u ) ).evaluate()( 0,0 );
    }
    // output
    else if ( output_index == 1 )
    {
        output = integrate( _range=elements( mesh ),
                            _expr=chi( ( Px() >= -0.1 ) && ( Px() <= 0.1 ) )*idv( u ) ).evaluate()( 0,0 )/0.2;
    }
    else
        throw std::logic_error( "[unsteadyHeat1d::output] error with output_index : only 0 or 1 " );

    return output;
}

FEELPP_CRB_PLUGIN( UnsteadyHeat1D, unsteadyheat1d )
}

