/* -*- mode: c++ coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-07-15
 */
#include <life/options.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifepoly/polynomialset.hpp>
#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/region.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmsh.hpp>


#include <life/lifevf/vf.hpp>


using namespace Life;

inline
po::options_description
makeOptions()
{
    po::options_description myintegralsoptions("MyFunctionSpace options");
    myintegralsoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.2 ), "mesh size")
        ("shape", Life::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
        ("alpha", Life::po::value<double>()->default_value( 3 ), "Regularity coefficient for function f")
        ;
    return myintegralsoptions.add( Life::life_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "myfunctionspace" ,
                     "myfunctionspace" ,
                     "0.3",
                     "nD(n=1,2,3) MyFunctionSpace on simplices or simplex products",
                     Life::AboutData::License_GPL,
                     "Copyright (c) 2008-2010 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


/**
 * MyFunctionSpace: compute integrals over a domain
 * \see the \ref ComputingIntegrals section in the tutorial
 * @author Christophe Prud'homme
 */
template<int Dim>
class MyFunctionSpace
    :
    public Simget
{
    typedef Simget super;
public:
    //! Polynomial order \f$P_2\f$
    static const uint16_type Order = 2;

    typedef double value_type;

    //! mesh
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //# marker1 #
    //! function space that holds piecewise constant (\f$P_0\f$) functions (e.g. to store material properties or partitioning
    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar> >, Discontinuous > p0_space_type;
    //! an element type of the \f$P_0\f$ discontinuous function space
    typedef typename p0_space_type::element_type p0_element_type;
    //# endmarker1 #

    //# marker2 #
    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;
    //# endmarker2 #

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    MyFunctionSpace( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>()  ),
        exporter()
    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:

    double meshSize;
    std::string shape;
    export_ptrtype exporter;
}; // MyFunctionSpace

template<int Dim> const uint16_type MyFunctionSpace<Dim>::Order;

template<int Dim>
void
MyFunctionSpace<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute MyFunctionSpace<" << Dim << ">\n";
    std::vector<double> X( 2 );
    X[0] = meshSize;
    if ( shape == "hypercube" )
        X[1] = 1;
    else // default is simplex
        X[1] = 0;
    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template<int Dim>
void
MyFunctionSpace<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    using namespace Life::vf;

    if ( X[1] == 0 ) shape = "simplex";
    if ( X[1] == 1 ) shape = "hypercube";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "doc/tutorial/%1%/%2%/h_%3%/" )
                                       % this->about().appName()
                                       % shape
                                       % meshSize );

    //# marker31 #
    //! create the mesh
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name= (boost::format( "%1%-%2%" ) % shape % Dim).str() ,
                                                      _shape=shape,
                                                      _dim=Dim,
                                                      _h=X[0] ) );
    //# endmarker31 #

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    //# marker32 #
    // function space \f$ X_h \f$
    space_ptrtype Xh = space_type::New( mesh );
    // an element of the function space \f$ X_h \f$
    auto u = Xh->element( "u" );
    // another element of the function space \f$ X_h \f$
    element_type v( Xh, "v" );
    //# endmarker32 #
    /** \endcode */

    value_type alpha = this->vm()["alpha"].template as<double>();
    value_type pi = M_PI;

    //# marker4 #
    auto g = sin(pi*Px()/2)*cos(pi*Py()/2)*cos(pi*Pz()/2);
    auto f = (1-Px()*Px())*(1-Py()*Py())*(1-Pz()*Pz())*pow(trans(vf::P())*vf::P(),(alpha/2.0));
    //# endmarker4 #

    //# marker5 #
    u = vf::project( Xh, elements(mesh), g );
    v = vf::project( Xh, elements(mesh), f );
    //# endmarker5 #

    //# marker6 #
    double L2g2 = integrate( elements(mesh), g*g ).evaluate()(0,0);
    double L2uerror2 = integrate( elements(mesh), (idv(u)-g)*(idv(u)-g) ).evaluate()(0,0);
    Log() << "||u-g||_0=" << math::sqrt( L2uerror2/L2g2 ) << "\n";
    double L2f2 = integrate( elements(mesh), f*f ).evaluate()(0,0);
    double L2verror2 = integrate( elements(mesh), (idv(v)-f)*(idv(v)-f) ).evaluate()(0,0);
    Log() << "||v-f||_0=" << math::sqrt( L2verror2/L2f2 ) << "\n";
    //# endmarker6 #

    //# marker7 #
    exporter = export_ptrtype( Exporter<mesh_type>::New( this->vm(), (boost::format( "%1%-%2%-%3%" ) % this->about().appName() % shape % Dim).str() ) );
    exporter->step(0)->setMesh( mesh );
    exporter->step(0)->add( "g", u );
    exporter->step(0)->add( "f", v );
    exporter->save();
    //# endmarker7 #

} // MyFunctionSpace::run

int
main( int argc, char** argv )
{
    Application app( argc, argv, makeAbout(), makeOptions() );

    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    app.add( new MyFunctionSpace<1>( app.vm(), app.about() ) );
    app.add( new MyFunctionSpace<2>( app.vm(), app.about() ) );
    app.add( new MyFunctionSpace<3>( app.vm(), app.about() ) );

    app.run();
}





