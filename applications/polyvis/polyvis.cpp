/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-04-17

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
/**
   \file polyvis.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-04-17
 */
/** include predefined life command line options */
#include <life/options.hpp>

#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>

/** include function space class */
#include <life/lifediscr/functionspace.hpp>

/** include helper function to define \f$P_0\f$ functions associated with regions  */
#include <life/lifediscr/region.hpp>

/** include integration methods */
#include <life/lifepoly/im.hpp>
/** include gmsh mesh importer */
#include <life/lifefilters/importergmsh.hpp>

/** include exporter factory class */
#include <life/lifefilters/exporter.hpp>

/** include gmsh generator for tensorized domains */
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifefilters/gmshsimplexdomain.hpp>

/** include  polynomialset header */
#include <life/lifepoly/polynomialset.hpp>
#include <life/lifepoly/lagrange.hpp>
//#include <life/lifepoly/raviartthomas.hpp>



/** include  the header for the variational formulation language (vf) aka FEEL++ */
#include <life/lifevf/vf.hpp>

/** use Life namespace */
using namespace Life;
using namespace Life::vf;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Life::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description polyvisoptions("Polyvis options");
    polyvisoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.5 ), "mesh size")
        ;
    return polyvisoptions.add( Life::life_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Life::Application subclass.
 *
 * \return some data about the application.
 */
inline
AboutData
makeAbout()
{
    AboutData about( "polyvis" ,
                     "polyvis" ,
                     "0.2",
                     "nD(n=1,2,3) Polynomial on simplices or simplex products",
                     Life::AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


/**
 * \class Polyvis
 *
 */
template<int Dim, typename Basis>
class Polyvis
    :
    public Application
{
    typedef Application super;
public:

    //! Polynomial order \f$P_2\f$
    static const uint16_type Order = 0;

    //! numerical type is double
    typedef double value_type;

    //! geometry entities type composing the mesh, here Simplex in
    //! Dimension Dim of Order 1
    typedef Simplex<Dim> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    typedef bases<Basis> basis_type;
    //typedef bases<Lagrange<Order,Scalar> > basis_type;
    //typedef bases<Lagrange<Order,Vectorial> > basis_type;
    //typedef bases<RaviartThomas<Order> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;
    //! the time set type  to save data over time
    typedef typename export_type::timeset_type timeset_type;

    struct Factory
    {
        typedef Singleton< Factory< Polyvis<Dim, Basis>, std::string > > type;
    };

    /**
     * Constructor
     */
    Polyvis( int argc, char** argv,
             AboutData const& ad,
             po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm()["exporter"].template as<std::string>() )->setOptions( this->vm() ) ),
        timeSet( new timeset_type( "polyvis" ) )
    {
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "polyvis" );
    }

    static Polyvis<Dim,Basis>* New( std::string const& exporter )
    {
        return Factory::type::instance().createObject( exporter );
    }

    /**
     * create the mesh using mesh size \c meshSize
     *
     * \param meshSize the mesh characteristic size
     * \return a mesh of an hyper cube of dimension Dim
     */
    mesh_ptrtype createMesh( double meshSize );

    /**
     * run the convergence test
     */
    void run();

private:

    //! mesh characteristic size
    double meshSize;

    //! exporter factory
    export_ptrtype exporter;

    //! timeset data structure the holds results over time
    typename export_type::timeset_ptrtype timeSet;

}; // Polyvis

template<int Dim> const uint16_type Polyvis<Dim>::Order;

template<int Dim>
typename Polyvis<Dim>::mesh_ptrtype
Polyvis<Dim>::createMesh( double meshSize )
{
    /** instantiate a new mesh */
    /** \code */
    mesh_ptrtype mesh( new mesh_type );
    /** \endcode */

    //! generate a tensorized domain (hyper cube of dimension Dim)
    /** \code */
    GmshSimplexDomain<convex_type::nDim, convex_type::nOrder> td( Gmsh::GMSH_REFERENCE_DOMAIN );
    td.setCharacteristicLength( meshSize );
    std::string fname = td.generate( convex_type::name().c_str() );
    /** \endcode */

    //! importer the mesh generated by td
    /** \code */
    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    /** \endcode */

    return mesh;
} // Polyvis::createMesh


template<int Dim>
void
Polyvis<Dim>::run()
{
    /**
     * print help if --help is passed to the command line
     */
    /** \code */
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }
    /** \endcode */

    /**
     * we change to the directory where the results and logs will be
     * stored
     */
    /** \code */
    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % convex_type::name()
                            % Order
                            % this->vm()["hsize"].template as<double>()
                            );
    /** \endcode */
    /**
     * First we create the mesh
     */
    /** \code */
    mesh_ptrtype oneelement_mesh = createMesh( 2 );
    mesh_ptrtype mesh = createMesh( meshSize );
    /** \endcode */

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    space_ptrtype Xh = space_type::New( oneelement_mesh );
    element_type U( Xh, "U" );
    element_type V( Xh, "V" );
    space_ptrtype Yh = space_type::New( mesh );
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    /** \endcode */

    typename timeset_type::step_ptrtype timeStep = timeSet->step( 0 );
    timeStep->setMesh( mesh );

    for( int i = 0;i < Xh->nLocalDof(); ++i )
        {
            U.zero();
            U( i ) = 1;
            std::ostringstream ostr;
            ostr << Xh->basis()->familyName() << "-" << i;
            timeStep->add( ostr.str(), U );
        }

    exporter->save();

} // Polyvis::run

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    /**
     * intantiate a Polyvis<Dim> class with Dim=2 (e.g. geometric dimension is 2)
     */
    /** \code */
    Polyvis<2> polyvis( argc, argv, makeAbout(), makeOptions() );
    /** \encode */

    /**
     * run the application
     */
    /** \code */
    polyvis.run();
    /** \endcode */
}






