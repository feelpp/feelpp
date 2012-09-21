/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-07-15

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file myfunctionspace.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-07-15
 */
#include <feel/feel.hpp>

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description myintegralsoptions( "MyFunctionSpace options" );
    myintegralsoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.2 ), "mesh size" )
    ( "dim", po::value<int>()->default_value( 0 ), "mesh dimension (0: all dimensions, 1,2 or 3)" )
    ( "order", po::value<int>()->default_value( 0 ), "approximation order (0: all orders, 1,2,3,4 or 5)" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex, hypercube or ellipsoid)" )
    ( "alpha", Feel::po::value<double>()->default_value( 3 ), "Regularity coefficient for function f" )
    ;
    return myintegralsoptions.add( Feel::feel_options() );
}

/**
 * MyFunctionSpace: compute integrals over a domain
 * \see the \ref ComputingIntegrals section in the tutorial
 * @author Christophe Prud'homme
 */
template<int Dim, int Order = 2>
class MyFunctionSpace
    :
public Simget
{
    typedef Simget super;
public:
    typedef double value_type;

    //! mesh
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! function space that holds piecewise constant (\f$P_0\f$) functions (e.g. to store material properties or partitioning
    //# marker1 #
    typedef FunctionSpace<mesh_type,bases<Lagrange<0,Scalar,Discontinuous> > >
    p0_space_type;
    //# endmarker1 #

    //! an element type of the \f$P_0\f$ discontinuous function space
    typedef typename p0_space_type::element_type p0_element_type;

    //# marker2 #
    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar,Continuous> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;
    //# endmarker2 #

    /* export */
    typedef Exporter<mesh_type,1> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    MyFunctionSpace()
        :
        super(),
        dim( this->vm()["dim"].template as<int>() ),
        order( this->vm()["order"].template as<int>() ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>()  ),
        exporter()
    {
    }

    void run();

private:

    int dim;
    int order;
    double meshSize;
    std::string shape;
    export_ptrtype exporter;
}; // MyFunctionSpace


template<int Dim, int Order>
void
MyFunctionSpace<Dim,Order>::run()
{
    Environment::changeRepository( boost::format( "doc/manual/tutorial/%1%/%2%/h_%3%/" )
                                   % this->about().appName()
                                   % shape
                                   % meshSize );

    //# marker31 #
    //! create the mesh
    mesh_ptrtype mesh =
        createGMSHMesh( _mesh=new mesh_type,
                        _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % shape % Dim % Order ).str() ,
                                      _shape=shape,
                                      _dim=Dim,
                                      _order=Order,
                                      _h=meshSize ) );

    //# endmarker31 #
    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    // function space \f$ X_h \f$
    //# marker32 #
    space_ptrtype Xh = space_type::New( mesh );
    //# endmarker32 #

    //# marker33 #
    // an element of the function space X_h
    auto u = Xh->element( "u" );
    // another element of the function space X_h
    element_type v( Xh, "v" );
    auto w = Xh->element( "w" );
    //# endmarker33 #
    /** \endcode */

    value_type alpha = this->vm()["alpha"].template as<double>();
    value_type pi = M_PI;

    //# marker4 #
    auto g = sin( 2*pi*Px() )*cos( 2*pi*Py() )*cos( 2*pi*Pz() );
    auto f =( 1-Px()*Px() )*( 1-Py()*Py() )*( 1-Pz()*Pz() )*pow( trans( vf::P() )*vf::P(),( alpha/2.0 ) );
    //# endmarker4 #

    //# marker5 #
    u = vf::project( _space=Xh, _range=elements( mesh ), _expr=g );
    v = vf::project( _space=Xh, _range=elements( mesh ), _expr=f );
    w = vf::project( _space=Xh, _range=elements( mesh ), _expr=idv( u )-g );
    //# endmarker5 #

    //# marker6 #
    double L2g2 = integrate( elements( mesh ), g*g ).evaluate()( 0,0 );
    double L2uerror2 = integrate( elements( mesh ), ( idv( u )-g )*( idv( u )-g ) ).evaluate()( 0,0 );
    LOG(INFO) << "||u-g||_0=" << math::sqrt( L2uerror2/L2g2 ) << "\n";
    double L2f2 = integrate( elements( mesh ), f*f ).evaluate()( 0,0 );
    double L2verror2 = integrate( elements( mesh ), ( idv( v )-f )*( idv( v )-f ) ).evaluate()( 0,0 );
    LOG(INFO) << "||v-f||_0=" << math::sqrt( L2verror2/L2f2 ) << "\n";
    //# endmarker6 #

    //# marker7 #
    // exporting to paraview or gmsh
    LOG(INFO) << "exporting\n" << std::endl;
    exporter = export_ptrtype( export_type::New( this->vm(), ( boost::format( "%1%-%2%-%3%-%4%" ) % this->about().appName() % shape % Dim % Order ).str() ) );
    LOG(INFO) << "exporting mesh \n" << std::endl;
    exporter->step( 0 )->setMesh( mesh );
    auto P0h = p0_space_type::New( mesh );
    LOG(INFO) << "saving pid\n" << std::endl;
    exporter->step( 0 )->add( "pid", regionProcess( P0h ) );

    exporter->step( 0 )->add( "g", u );
    exporter->step( 0 )->add( "u-g", w );
    exporter->step( 0 )->add( "f", v );

    exporter->save();
    //# endmarker7 #
} // MyFunctionSpace::run

int
main( int argc, char** argv )
{
    /**
     * Initialize Feel++ Environment
     */
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="myfunctionspace",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    Application app;

    app.add( new MyFunctionSpace<2,3>() );
    app.add( new MyFunctionSpace<3,3>() );
    app.run();
}
