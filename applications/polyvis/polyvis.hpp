/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-17

  Copyright (C) 2009-2011 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-17
 */
#ifndef __Polyvis_H
#define __Polyvis_H 1

#include <boost/parameter.hpp>

/** include predefined feel command line options */
#include <feel/options.hpp>

#include <feel/feelcore/parameter.hpp>

/** include function space class */
#include <feel/feeldiscr/functionspace.hpp>

/** include helper function to define \f$P_0\f$ functions associated with regions  */
#include <feel/feeldiscr/region.hpp>

/** include integration methods */
#include <feel/feelpoly/im.hpp>
/** include gmsh mesh importer */
#include <feel/feelfilters/gmsh.hpp>

/** include exporter factory class */
#include <feel/feelfilters/exporter.hpp>

/** include  polynomialset header */
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/lagrange.hpp>


/** include  the header for the variational formulation language (vf) aka FEEL++ */
#include <feel/feelvf/vf.hpp>

/** include linear algebra backend */
#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/domain.hpp>

#include "polyvisbase.hpp"

using namespace Feel;
namespace Feel
{
typedef parameter::parameters<
parameter::required<tag::convex_type, boost::is_base_and_derived<ConvexBase,_> >,
          parameter::required<tag::basis_type, boost::is_class<_> >
          > polyvis_signature;

//Generates one-element mesh for reference element (works independently from gmsh version)
gmsh_ptrtype oneelement_geometry_ref();

//Generates one-element mesh (!= reference) (works independently from gmsh version)
// This element corresponds to apply an homothetic transformation (ratio=2, center=(0,0) )
gmsh_ptrtype oneelement_geometry_real();


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

    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;


    //! geometry entities type composing the mesh, here Simplex in
    //! Dimension Dim of Order 1
    typedef typename parameter::binding<args, tag::convex_type, Simplex<2> >::type convex_type;

    //! the basis type of our approximation space
    typedef typename parameter::binding<args, tag::basis_type, Lagrange<3> >::type basis_type;

    //! geometric dimension
    static const uint16_type Dim = convex_type::nDim;

    //! Polynomial order \f$P_2\f$
    //static const uint16_type Order = basis_type::nOrder;

    // //! numerical type is double
    // typedef double value_type;

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
        super(),
        meshSize( Environment::vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( Environment::vm() ) )
        {
            this->init( vm );
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

    //! linear algebra backend
    backend_ptrtype M_backend;

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
void
Polyvis<A0,A1,A2,A3,A4>::run()
{
    std::string fname = soption(_name=(boost::format("meshes-%1%d")%Dim).str())+".msh";
    // First we create the mesh with one element
    mesh_ptrtype oneelement_mesh = loadMesh( _mesh=new mesh_type, _filename=fname );

#if 0
    mesh_ptrtype oneelement_mesh_ref = createGMSHMesh( _mesh=new mesh_type,
                                                       _desc = oneelement_geometry_ref() );

    mesh_ptrtype oneelement_mesh_real = createGMSHMesh( _mesh=new mesh_type,
                                                        _desc = oneelement_geometry_real() );
#endif
    // then a fine mesh which we use to export the basis function to
    // visualize them
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name="fine",
                                                      _shape="simplex",
                                                      _dim=Dim,
                                                      _h=meshSize,
                                                      _xmin=Environment::vm()["xmin"].template as<double>(),
                                                      _ymin=Environment::vm()["ymin"].template as<double>(),
                                                      _zmin=Environment::vm()["zmin"].template as<double>() ) );
    /** \endcode */

    using namespace Feel::vf;

    /**
     * The function space and some associated elements(functions) are
     * then defined
     */
    space_ptrtype Xh_ref = space_type::New( oneelement_mesh ); // space associated with reference element
    //space_ptrtype Xh_real = space_type::New( oneelement_mesh_real ); // space associated with real element

    std::cout << "Family = " << Xh_ref->basis()->familyName() << "\n"
              << "Dim    = " << Xh_ref->basis()->nDim << "\n"
              << "Order  = " << Xh_ref->basis()->nOrder << "\n"
              << "NDof   = " << Xh_ref->nLocalDof() << "\n";

    // U = shape function on current dof (on reference element)
    element_type U( Xh_ref, "U" );

    // set the mesh of the exporter, we use the fine mesh and the
    // exporter does all the interpolation
    exporter->step( 0 )->setMesh( mesh );

    for ( size_type i = 0; i < Xh_ref->nLocalDof(); ++i )
    {
        U.zero();
        U( i ) = 1;
#if 0 // for rtk
        using namespace vf;
        std::cout << "flux " << i << " = " << integrate( boundaryfaces( oneelement_mesh_ref ),
                                                         trans( idv( U ) )*N(), _Q<4>() ).evaluate()( 0, 0 ) << "\n";;
        std::cout << "div " << i << " = " << integrate( elements( oneelement_mesh_ref ),
                                                        divv( U ), _Q<4>() ).evaluate()( 0, 0 ) << "\n";
#endif
        std::ostringstream ostr;
        ostr << Xh_ref->basis()->familyName() << "-" << i;
        exporter->step( 0 )->add( ostr.str(), U );
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
