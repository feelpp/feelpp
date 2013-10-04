/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2007-07-12

  Copyright (C) 2007 Unil

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
   \file pbeqspace_implementation.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-07-12
 */
#ifndef _PBEQSPACE_IMPLEMENTATION_HPP
#define _PBEQSPACE_IMPLEMENTATION__HPP

#include <set>
#include <vector>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <feel/feelpoly/operations.hpp>
#include <feel/feelpoly/polynomial.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelvf/vf.hpp>

#include "mymesh.hpp"
#include "pbeqspace.hpp"
#include "heavysidefunction.hpp"

namespace Feel
{

PbeqSpace::PbeqSpace( )
    :
    M_meshSize( 1 ),
    M_farFactor( 1 ),
    M_farBnd( 1 ),
    M_stretch( Dim ),
    M_invStretch( Dim ),
    M_translation( Dim ),
    M_meshSetted( false ),
    M_spacesSetted( false ),
    M_mesh(),
    M_HeavysideSpace(),
    M_Space(),
    timers(),
    stats()
{
    std::fill( M_stretch.begin(), M_stretch.end(), 1.0 );
    std::fill( M_translation.begin(), M_translation.end(), 0.0 );

    setUpInvStretch();

} // PbeqSpace::PbeqSpace

PbeqSpace::PbeqSpace( value_type const meshSize,
                      value_type const farFactor,
                      value_type const farBnd,
                      node_type& stretch,
                      node_type& translation )
    :
    M_meshSize( meshSize ),
    M_farFactor( farFactor ),
    M_farBnd( farBnd ),
    M_stretch( stretch ),
    M_invStretch( stretch ),
    M_translation( translation ),
    M_meshSetted( false ),
    M_spacesSetted( false ),
    M_mesh(),
    M_HeavysideSpace(),
    M_Space(),
    timers(),
    stats()
{
    setUpInvStretch();
} // PbeqSpace::PbeqSpace

PbeqSpace::PbeqSpace( value_type const meshSize,
                      value_type const farFactor,
                      value_type const farBnd )
    :
    M_meshSize( meshSize ),
    M_farFactor( farFactor ),
    M_farBnd( farBnd ),
    M_stretch( Dim ),
    M_invStretch( Dim ),
    M_translation( Dim ),
    M_meshSetted( false ),
    M_spacesSetted( false ),
    M_mesh(),
    M_HeavysideSpace(),
    M_Space(),
    timers(),
    stats()
{
    std::fill( M_stretch.begin(), M_stretch.end(), 1.0 );
    std::fill( M_translation.begin(), M_translation.end(), 0.0 );

    setUpInvStretch();

} // PbeqSpace::PbeqSpace


PbeqSpace::PbeqSpace( PbeqSpace const& tc )
    :
    M_meshSize( tc.M_meshSize ),
    M_farFactor( tc.M_farFactor ),
    M_farBnd( tc.M_farBnd ),
    M_stretch( tc.M_stretch ),
    M_invStretch( tc.M_invStretch ),
    M_translation( tc.M_translation ),
    M_meshSetted( tc.M_meshSetted ),
    M_spacesSetted( tc.M_spacesSetted ),
    M_mesh( tc.M_mesh ),
    M_HeavysideSpace( tc.M_HeavysideSpace ),
    M_Space( tc.M_Space ),
    timers( tc.timers ),
    stats( tc.stats )
{
    LOG(INFO) << "[PbeqSpace] hsize = " << M_meshSize << "\n";
} // PbeqSpace::PbeqSpace


//const uint16_type PbeqSpace::Dim(3);

void  PbeqSpace::setUpInvStretch()
{
    M_invStretch.resize( M_stretch.size() );

    for ( uint i( 0 ); i < M_stretch.size(); ++i )
        M_invStretch[i] = 1/M_stretch[i];

}

void
PbeqSpace::createSpaces()
{

    if ( M_spacesSetted )
    {
        LOG(INFO) << "[PbeqSpace] createSpaces(): trying to rebuilding spaces: returning right away" << "\n";
        return;
    }

    using namespace Feel::vf;

    /*
     * since the mesh is used by a P0 space and P1 space, we update
     * the mesh with the proper components for each space. Thus we do
     * not let the first space(here p0) decide what to update in the
     * mesh for the next space (here p1).
     */
    //M_mesh->setComponents(size_type( MESH_CHECK | MESH_UPDATE_FACES | MESH_RENUMBER | MESH_PARTITION ) );
    M_mesh->setComponents( size_type( MESH_ALL_COMPONENTS ) );
    M_mesh->updateForUse();

    /*
     * The function space and some associate elements are then defined
     */
    timers["init"].first.restart();
    M_HeavysideSpace = heavyside_space_type::New( M_mesh );
    M_Space = space_type::New( M_mesh );

    timers["init"].second = timers["init"].first.elapsed();
    M_HeavysideSpace->dof()->showMe();
    M_Space->dof()->showMe();
    stats["ndof"] = M_HeavysideSpace->nDof();

    /*
     * a quadrature rule for numerical integration
     */

    LOG(INFO) << "[timer] createSpaces() (" << M_mesh->numElements() << " Elems,  "
          << M_HeavysideSpace->dof()->nDof() << " P" << OrderHeavyside << " DOFs, "
          << M_Space->dof()->nDof() << " P" << Order << "DOFs): "
          << timers["init"].second << "\n";

    M_spacesSetted = true;


} // PbeqSpace::createSpaces


void
PbeqSpace::createMesh( bool const geoOnly )
{

    if ( M_meshSetted && !geoOnly )
    {
        stats["nelt"] = M_mesh->elements().size();
        createSpaces();
        return;
    }

    timers["mesh"].first.restart();

    myMesh td;
    td.setCharacteristicLength( M_meshSize );

    td.setFarCharacteristic( M_farFactor, M_farBnd );

    if ( geoOnly )
    {
        td.generateGeo( entity_type::name().c_str() );
        return;
    }

    timers["mesh"].second = timers["mesh"].first.elapsed();
    LOG(INFO) << "[timer] createMesh(): " << timers["mesh"].second << "\n";

    loadMesh( td.generate( entity_type::name().c_str() ) ) ;
}

bool
PbeqSpace::loadMesh( std::string const& meshname )
{

    LOG(INFO) << "[PbeqSpace] hsize = " << M_meshSize << "\n";

    timers["mesh"].first.restart();
    M_mesh.reset( new mesh_type );

    namespace fs = boost::filesystem;
    LOG(INFO) << "mesh file name: " << meshname << "\n";
    LOG(INFO) << "does mesh file name exists ?: " << fs::exists( meshname ) << "\n";
    fs::path meshpath( meshname );

    if ( !fs::exists( meshpath ) ) return false;

    ImporterGmsh<mesh_type> import( meshname );
    M_mesh->accept( import );

    timers["mesh"].second = timers["mesh"].first.elapsed();
    LOG(INFO) << "[timer] loadMesh(): " << timers["mesh"].second << "\n";

    M_meshSetted = true;

    /*
     * Now we create the spaces
     */

    stats["nelt"] = M_mesh->elements().size();
    createSpaces();

    return true;

} // PbeqSpace::createMesh




PbeqSpace::heavyside_element_type
PbeqSpace::heavyside( value_type const& radius,
                      node_type  const& center ) const
{
    assert( M_meshSetted );

    using namespace Feel::vf;

    value_type const R2( radius*radius );

    heavyside_element_type Ha( M_HeavysideSpace, "Ha" );

    node_type centerHat( center );

    for ( int i( 0 ); i< Dim; i++ )
        centerHat[i] -= M_translation[i];

    AUTO( Ellipse,
          ( M_stretch[0]*Px() - centerHat[0] )*( M_stretch[0]*Px() - centerHat[0] )  +
          ( M_stretch[1]*Py() - centerHat[1] )*( M_stretch[1]*Py() - centerHat[1] )  +
          ( M_stretch[2]*Pz() - centerHat[2] )*( M_stretch[2]*Pz() - centerHat[2] )  );

    Ha = project( M_HeavysideSpace, elements( *M_mesh ), chi( ( R2 - ( Ellipse ) ) < 0 ) );

    return Ha;

} // end PbeqSpace::Heavyside

PbeqSpace::heavyside_element_type
PbeqSpace::heavyside( molecule_type const& molecule ) const
{
    using namespace Feel::vf;
    LOG(INFO) << "heavyside call\n";
    heavyside_element_type H( M_HeavysideSpace, "H" );

    std::fill( H.begin(), H.end(), 1.0 );

    molecule_type::atoms_const_iterator_type atom( molecule.begin() );

    if ( atom == molecule.end() ) return H;

    heavyside_element_type H1;

    for ( ; atom != molecule.end(); atom++ )
    {
        H1 = heavyside( atom->radius(), atom->center() );

        //H = ublas::element_prod( H, heavyside(atom->radius(), atom->center())  );
        //H = ublas::element_prod( H, H1  );
        for ( size_type i = 0; i < H.localSize(); ++i )
            H( H.firstLocalIndex()+ i ) *= H1( H1.firstLocalIndex() + i );

    }

    LOG(INFO) << "heavyside call done\n";
    return H;


} // end PbeqSpace::Heavyside


PbeqSpace::heavyside_element_type
PbeqSpace::chargeDensity( value_type const& radius,
                          value_type const& charge,
                          node_type  const& center ) const
{
    assert( M_meshSetted );

    using namespace Feel::vf;

    value_type R2, qa;
    // radius of the dirac distribution
    R2 = radius*radius;
    // qa = charge * area of supp delta
    qa = ( charge ) *  2.35619449019234 * R2*radius; // 3/4*pi*r^3


    heavyside_element_type rhoa( M_HeavysideSpace, "rhoa" );

    node_type centerHat( center );

    for ( int i( 0 ); i< Dim; i++ )
        centerHat[i] -= M_translation[i];

    AUTO( Ellipse,
          ( M_stretch[0]*Px() - centerHat[0] )*( M_stretch[0]*Px() - centerHat[0] )  +
          ( M_stretch[1]*Py() - centerHat[1] )*( M_stretch[1]*Py() - centerHat[1] )  +
          ( M_stretch[2]*Pz() - centerHat[2] )*( M_stretch[2]*Pz() - centerHat[2] )  );

    rhoa = project( M_HeavysideSpace, elements( *M_mesh ), qa*chi( ( R2 - ( Ellipse ) ) > 0 ) );

    return rhoa;

} // end PbeqSpace::chargeDensity



void
PbeqSpace::intvrho( molecule_type const& molecule,
                    vector_ptrtype rhs ) const
{

    std::map<std::string,std::pair<boost::timer,value_type> > timer;
    timer["init intvrho"].first.restart();

    element_type v( M_Space, "v" );
    std::fill( v.begin(), v.end(), .0 );

    molecule_type::atoms_const_iterator_type atom( molecule.begin() );

    if ( atom == molecule.end() ) return;

    mesh_type::Inverse meshinv( M_mesh );
    std::vector<value_type> atomcharges( molecule.size() );

    /* initialisation of the mesh::inverse data structure */
    for ( size_type atomid = 0; atom != molecule.end(); ++atom, ++atomid )
    {
        meshinv.addPointWithId( element_prod( M_invStretch , atom->center() - M_translation ),
                                atomid, 0 );
        atomcharges[atomid] = atom->charge();
    }

    meshinv.distribute();

    std::vector<bool> dof_done( molecule.size() );
    std::fill( dof_done.begin(), dof_done.end(), false );
    std::vector<boost::tuple<size_type,uint16_type > > itab;
    matrix_node<value_type>::type pts( mesh_type::nDim, 1 );
    typedef mesh_type::element_type geoelement_type;
    typedef mesh_type::element_iterator mesh_element_iterator;

    mesh_element_iterator it = M_mesh->beginElementWithProcessId( Application::processId() );
    mesh_element_iterator en = M_mesh->endElementWithProcessId( Application::processId() );

    // geometric mapping context
    typedef mesh_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef gm_type::Context<vm::POINT, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;

    // basis
    typedef space_type::basis_type basis_type;
    basis_type const* __basis = M_Space->basis().get();
    gm_ptrtype __gm = M_Space->gm();
    typedef gm_type::precompute_ptrtype geopc_ptrtype;
    typedef gm_type::precompute_type geopc_type;
    geopc_ptrtype __geopc( new geopc_type( __gm, __basis->dual().points() ) );
    gmc_ptrtype __c( new gmc_type( __gm, *it, __geopc ) );


    /* note: meshmover.hpp has this two lines: */
    if ( !v.areGlobalValuesUpdated() )
        v.updateGlobalValues();

    /*shall we use this ?
    */

    timer["init intvrho"].second = timer["init intvrho"].first.elapsed();
    LOG(INFO) << "[timer] init intvrho(): " << timer["init intvrho"].second << "\n";

    timer["intvrho"].first.restart();

    // size_type first_dof = M_Space->dof()->firstDof();
    for ( ; it != en; ++ it )
    {
        __c->update( *it, __geopc );
        meshinv.pointsInConvex( it->id(), itab );

        if ( itab.size() == 0 )
            continue;

        for ( size_type i = 0; i < itab.size(); ++i )
        {
            // get dof id in target dof table
            size_type dof = boost::get<0>( itab[i] );

            if ( !dof_done[dof] )
            {
                dof_done[dof]=true;
                ublas::column( pts, 0 ) = meshinv.referenceCoords()[dof];

                __geopc->update( pts );
                __c->update( *it, __geopc );
                element_type::pc_type pc( v.functionSpace()->fe(), pts );

                for ( uint16_type loc_ind=0; loc_ind < basis_type::nLocalDof ; loc_ind++ )
                {
                    for ( uint16_type comp = 0; comp < basis_type::nComponents; ++comp )
                    {
                        size_type globaldof = boost::get<0>( M_Space->dof()->localToGlobal( it->id(), loc_ind, comp ) );

                        // update only values on the processor
                        if ( globaldof >= v.firstLocalIndex() &&
                                globaldof < v.lastLocalIndex() )
                        {
                            v.setGlobalValue( globaldof, 1 );

                            element_type::id_type interpfunc( v.id( *__c, pc ) );
                            //std::cout << "interpfunc :  " << interpfunc << "\n";

                            rhs->add( globaldof, atomcharges[dof] * interpfunc( comp, 0, 0 ) );
                            // DVLOG(2) << "rhs( " << globaldof << ")=" << (*rhs)( globaldof )
                            //           << " (just added " << atomcharges[dof] * interpfunc( comp, 0, 0 ) << " )" << "\n";
                            v.setGlobalValue( globaldof, 0 );

                        }
                    }
                }

            }
        }
    } // element

    timer["intvrho"].second = timer["intvrho"].first.elapsed();
    LOG(INFO) << "[timer] intvrho(): " << timer["intvrho"].second << "\n";


    return;
}

PbeqSpace::heavyside_element_type
PbeqSpace::chargeDensity( molecule_type const& molecule ) const
{
    using namespace Feel::vf;

    heavyside_element_type rho( M_HeavysideSpace, "rho" );

    std::fill( rho.begin(), rho.end(), .0 );

    value_type factor( 0 );

    for ( int i( 0 ); i< Dim; i++ )
        factor += M_stretch[i];

    factor = 0.5*Dim/factor;

    molecule_type::atoms_const_iterator_type atom( molecule.begin() );

    if ( atom == molecule.end() ) return rho;

    for ( ; atom != molecule.end(); atom++ )
        rho += chargeDensity( atom->radius() * factor, atom->charge(), atom->center() );

    return rho;


} // end PbeqSpace::chargeDensity


PbeqSpace::heavyside_element_type
PbeqSpace::fastHeavyside( molecule_type const& molecule ) const
{

    heavyside_element_type H( M_HeavysideSpace, "H" );

    //std::fill( H.begin(), H.end(), 1.0 );

    /* note: meshmover.hpp has this two lines: shall we use this ?
       THESE LINES ARE DANGEROUS : if you uncomment, H does not behave well anymore :-( )

    if ( !H.areGlobalValuesUpdated() )
        H.updateGlobalValues();
    */

    heavysideFunction heavyside;

    heavyside.setMolecule( &molecule );
    heavyside.setStretch( &M_stretch );
    heavyside.setTranslation( &M_translation );

    heavyside.setSmoothWindow( M_sW );

    heavyside_space_type::dof_type::dof_points_const_iterator it_dofpt = M_HeavysideSpace->dof()->dofPointBegin();
    heavyside_space_type::dof_type::dof_points_const_iterator en_dofpt = M_HeavysideSpace->dof()->dofPointEnd();

    for ( ; it_dofpt != en_dofpt ; ++it_dofpt )
    {
        size_type globaldof = boost::get<1>( *it_dofpt ) + M_HeavysideSpace->dof()->firstDof();

        if ( globaldof >= H.firstLocalIndex() && globaldof < H.lastLocalIndex() )
            H( globaldof ) = heavyside( boost::get<0>( *it_dofpt ) );

    }

    return H;

} // end fastHeavyside

} // end namespace Feel

#endif
