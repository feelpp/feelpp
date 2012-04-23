/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-07

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)

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
   \file myintegrals.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-07
 */
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>


#include <feel/feelvf/vf.hpp>


using namespace Feel;
//# marker8 #
inline
po::options_description
makeOptions()
{
    po::options_description myintegralsoptions( "MyIntegrals options" );
    myintegralsoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.2 ), "mesh size" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ;
    return myintegralsoptions.add( Feel::feel_options() );
}
//# endmarker8 #

//# marker9 #
inline
AboutData
makeAbout()
{
    AboutData about( "myintegrals" ,
                     "myintegrals" ,
                     "0.3",
                     "nD(n=1,2,3) MyIntegrals on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2010 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "" );
    return about;
}
//# endmarker9 #


/**
 * MyIntegrals: compute integrals over a domain
 * \see the \ref ComputingIntegrals section in the tutorial
 * @author Christophe Prud'homme
 */
template<int Dim>
class MyIntegrals
    :
public Simget
{
    typedef Simget super;
public:
    typedef double value_type;

    /*mesh*/
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    MyIntegrals( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>()  )
    {}

    //# marker10 #
    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );
    //# endmarker10 #

private:

    double meshSize;
    std::string shape;
}; // MyIntegrals


template<int Dim>
void
MyIntegrals<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute MyIntegrals<" << Dim << ">\n";
    std::vector<double> X( 2 );
    X[0] = meshSize;

    if ( shape == "ellipsoid" )
        X[1] = 2;

    else if ( shape == "hypercube" )
        X[1] = 1;

    else // default is simplex
        X[1] = 0;

    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template<int Dim>
void
MyIntegrals<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    using namespace Feel::vf;

    if ( X[1] == 0 ) shape = "simplex";

    if ( X[1] == 1 ) shape = "hypercube";

    if ( X[1] == 2 ) shape = "ellipsoid";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "doc/tutorial/%1%/%2%/h_%3%/" )
                                       % this->about().appName()
                                       % shape
                                       % meshSize );

    /*
     * First we create the mesh
     */
    auto mesh = mesh_type::New();
    bool is_mesh_loaded = mesh->load( _name="mymesh",_path=".",_type="text" );

    if ( !is_mesh_loaded )
    {
        //# marker11 #
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                               _desc=domain( _name= ( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                             _shape=shape,
                                             _order=1,
                                             _dim=Dim,
                                             _h=X[0] ),
                               _partitions=this->comm().size() );
        mesh->save( _name="mymesh",_path=".",_type="text" );
    }

    //# endmarker11 #

    /*
     * Compute domain Area
     */
    //# marker1 #
    double local_domain_area = integrate( _range=elements( mesh ),
                                          _expr=constant( 1.0 ) ).evaluate( false )( 0,0 );

    //# endmarker1 #
    //# marker2 #
    double global_domain_area= integrate( _range=elements( mesh ),
                                          _expr=constant( 1.0 ) ).evaluate()( 0,0 );
    //# endmarker2 #
    //# marker3 #
    Log() << "int_Omega 1 = " << global_domain_area
          << "[ " << local_domain_area << " ]\n";

    //# endmarker3 #
    if ( Dim > 1 )
    {
        Log() << "nb faces on boundary " << nelements( boundaryfaces(mesh) ) << "\n";
        /*
         * Compute domain perimeter
         */
        //# marker4 #
        double local_boundary_length =  integrate( boundaryfaces( mesh ),
                                        constant( 1.0 ) ).evaluate( false )( 0,0 );
        double global_boundary_length = integrate( boundaryfaces( mesh ),
                                        constant( 1.0 ) ).evaluate()( 0,0 );
        Log() << "int_BoundaryOmega (1)= " << global_boundary_length
              << "[ " << local_boundary_length << " ]\n";
        //# endmarker4 #
    }

    /*
     * Compute \int f where f= x^2 + y^2 + z^2
     * \note Py() = Pz() = 0 in 1D
     * \note Pz() = 0 in 2D
     */
    //# marker5 #
    double local_intf = integrate( elements( mesh ),
                                   Px()*Px() + Py()*Py() + Pz()*Pz() // trans(P())*P()
                                 ).evaluate( false )( 0,0 );
    double global_intf = integrate( elements( mesh ),
                                    Px()*Px() + Py()*Py() + Pz()*Pz() // trans(P())*P()
                                  ).evaluate()( 0,0 );
    //# endmarker5 #
    Log() << "int_Omega (x^2+y^2+z^2) = " << global_intf
          << "[ " << local_intf << " ]\n";

    //# marker6 #
    double global_intsin = integrate( elements( mesh ),
                                      sin( Px()*Px() + Py()*Py() + Pz()*Pz() )
                                    ).evaluate()( 0,0 );
    double local_intsin = integrate( elements( mesh ),
                                     sin( Px()*Px() + Py()*Py() + Pz()*Pz() ) ).evaluate( false )( 0,0 );
    //# endmarker6 #
    Log() << "int_Omega (sin(x^2+y^2+z^2)) [with order 4 max exact integration]= " << global_intsin << "\n";
    Log() << "int_Omega[" << this->comm().rank() << "] (sin(x^2+y^2+z^2)) [with order 4 max exact integration]= " << local_intsin << "\n";


    //# marker7 #
    double local_intsin2 = integrate( elements( mesh ),
                                      sin( Px()*Px() + Py()*Py() + Pz()*Pz() ),
                                      _Q<2>()
                                    ).evaluate( false )( 0,0 );
    double global_intsin2 = integrate( elements( mesh ),
                                       sin( Px()*Px() + Py()*Py() + Pz()*Pz() ),
                                       _Q<2>() ).evaluate()( 0,0 );
    //# endmarker7 #
    Log() << "int_Omega (sin(x^2+y^2+z^2)) [with order 2 max exact integration] = " << global_intsin2
          << "[ " << local_intsin2 << " ]\n";


} // MyIntegrals::run

int
main( int argc, char** argv )
{
    Feel::Environment env( argc, argv );

    Application app( argc, argv, makeAbout(), makeOptions() );

    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    //app.add( new MyIntegrals<1>( app.vm(), app.about() ) );
    app.add( new MyIntegrals<2>( app.vm(), app.about() ) );
    //app.add( new MyIntegrals<3>( app.vm(), app.about() ) );

    app.run();
}





