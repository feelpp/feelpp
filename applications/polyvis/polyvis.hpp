/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-04-17

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file polyvis.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-04-17
 */
#ifndef __Polyvis_H
#define __Polyvis_H 1

#include <boost/parameter.hpp>

/** include predefined life command line options */
#include <life/options.hpp>

#include <life/lifecore/parameter.hpp>

/** include function space class */
#include <life/lifediscr/functionspace.hpp>

/** include helper function to define \f$P_0\f$ functions associated with regions  */
#include <life/lifediscr/region.hpp>

/** include integration methods */
#include <life/lifepoly/im.hpp>
/** include gmsh mesh importer */
#include <life/lifefilters/gmsh.hpp>

/** include exporter factory class */
#include <life/lifefilters/exporter.hpp>

/** include  polynomialset header */
#include <life/lifepoly/polynomialset.hpp>
#include <life/lifepoly/lagrange.hpp>
//#include <life/lifepoly/raviartthomas.hpp>



/** include  the header for the variational formulation language (vf) aka FEEL++ */
#include <life/lifevf/vf.hpp>

#include "polyvisbase.hpp"

namespace Life
{
typedef parameter::parameters<
    parameter::required<tag::convex_type, boost::is_base_and_derived<ConvexBase,_> >,
    parameter::required<tag::basis_type, boost::is_class<_> >
    > polyvis_signature;


/**
 * \class Polyvis
 * \brief class for polynomial visualisation
 *
 * \code
 *  Polyvis<dim=2,basis_type=Lagrange>
 * \endcode
 */
template<
    typename A0,
    typename A1,
    typename A2 = parameter::void_,
    typename A3 = parameter::void_,
    typename A4 = parameter::void_>
class Polyvis
    :
    public PolyvisBase
{
    typedef PolyvisBase super;
public:

    typedef typename polyvis_signature::bind<A0,A1,A2,A3,A4>::type args;




    //! geometry entities type composing the mesh, here Simplex in
    //! Dimension Dim of Order 1
    typedef typename parameter::binding<args, tag::convex_type, Simplex<2> >::type convex_type;

    //! the basis type of our approximation space
    typedef typename parameter::binding<args, tag::basis_type, Lagrange<3> >::type basis_type;

    //! geometric dimension
    static const uint16_type Dim = convex_type::nDim;

    //! Polynomial order \f$P_2\f$
    static const uint16_type Order = basis_type::nOrder;

    //! numerical type is double
    typedef double value_type;

    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    //! the approximation function space type
    typedef FunctionSpace<mesh_type, bases<basis_type> > space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    Polyvis()
        :
        super(),
        meshSize( 0.1 ),
        exporter()
    {
    }

    /**
     * Constructor
     */
    Polyvis( po::variables_map vm )
        :
        super( vm ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
    }
    void init( po::variables_map const& vm )
    {
        super::init( vm );
        meshSize = vm["hsize"].template as<double>();
        exporter = export_ptrtype( Exporter<mesh_type>::New( vm, this->name() ) );
    }

    /**
     * run the convergence test
     */
    void run();

private:

    //! mesh characteristic size
    double meshSize;

    //! exporter factory
    export_ptrtype exporter;

}; // Polyvis

template<
    typename A0,
    typename A1,
    typename A2,
    typename A3,
    typename A4>
const uint16_type Polyvis<A0,A1,A2,A3,A4>::Order;

template<
    typename A0,
    typename A1,
    typename A2,
    typename A3,
    typename A4>
void
Polyvis<A0,A1,A2,A3,A4>::run()
{
    // First we create the mesh with one element
    mesh_ptrtype oneelement_mesh = createGMSHMesh( _mesh=new mesh_type,
                                                   _desc=domain( _name="one-elt",
                                                                 _shape="simplex",
                                                                 _dim=Dim,
                                                                 _h=2.0,
                                                                 _addmidpoint=false,
                                                                 _xmin=this->vm()["xmin"].template as<double>(),
                                                                 _ymin=this->vm()["ymin"].template as<double>(),
                                                                 _zmin=this->vm()["zmin"].template as<double>() ) );
    // then a fine mesh which we use to export the basis function to
    // visualize them
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name="fine",
                                                      _shape="simplex",
                                                      _dim=Dim,
                                                      _h=meshSize,
                                                      _xmin=this->vm()["xmin"].template as<double>(),
                                                      _ymin=this->vm()["ymin"].template as<double>(),
                                                      _zmin=this->vm()["zmin"].template as<double>() ) );
    /** \endcode */

    /**
     * The function space and some associated elements(functions) are
     * then defined
     */
    space_ptrtype Xh = space_type::New( oneelement_mesh );
    std::cout << "Family = " << Xh->basis()->familyName() << "\n"
              << "Dim    = " << Xh->basis()->nDim << "\n"
              << "Order  = " << Xh->basis()->nOrder << "\n"
              << "NDof   = " << Xh->nLocalDof() << "\n";
    element_type U( Xh, "U" );

    // set the mesh of the exporter, we use the fine mesh and the
    // exporter does all the interpolation
    exporter->step(0)->setMesh( mesh );

    for( size_type i = 0;i < Xh->nLocalDof(); ++i )
        {
            U.zero();
            U( i ) = 1;
#if 0 // for rtk
            using namespace vf;
            std::cout << "flux " << i << " = " << integrate( boundaryfaces(oneelement_mesh),_Q<4>(),
                                                             trans(idv(U))*N() ).evaluate()( 0, 0 ) << "\n";;
            std::cout << "div " << i << " = " << integrate( elements(oneelement_mesh), _Q<4>(),
                                                            divv(U) ).evaluate()( 0, 0 ) << "\n";
#endif
            std::ostringstream ostr;
            ostr << Xh->basis()->familyName() << "-" << i;
            exporter->step(0)->add( ostr.str(), U );
        }

    exporter->save();

} // Polyvis::run
template<
    typename A0,
    typename A1,
    typename A2,
    typename A3,
    typename A4>
const uint16_type Polyvis<A0,A1,A2,A3,A4>::Dim;
}

#endif // Polyvis_HPP
