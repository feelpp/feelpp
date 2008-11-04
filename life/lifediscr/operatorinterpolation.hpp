/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-01

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file operatorinterpolation.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-01
 */
#ifndef __OperatorInterpolation_H
#define __OperatorInterpolation_H 1

#include <life/lifediscr/operatorlinear.hpp>

namespace Life
{
/**
 * \class OperatorInterpolation
 * \brief Global interpolation operator
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename DomainSpaceType, typename ImageSpaceType>
class OperatorInterpolation : public OperatorLinear<DomainSpaceType, ImageSpaceType>
{
    typedef OperatorLinear<DomainSpaceType, ImageSpaceType> super;
public:


    /** @name Typedefs
     */
    //@{

    /*
     * domain
     */
    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename domain_space_type::mesh_type domain_mesh_type;
    typedef typename domain_mesh_type::element_type domain_geoelement_type;
    typedef typename domain_mesh_type::element_iterator domain_mesh_element_iterator;


    typedef typename super::backend_ptrtype backend_ptrtype;



    /*
     * image
     */
    typedef typename super::dual_image_space_type dual_image_space_type;
    typedef typename super::dual_image_space_ptrtype dual_image_space_ptrtype;
    typedef typename dual_image_space_type::value_type value_type;
    typedef typename dual_image_space_type::mesh_type image_mesh_type;
    typedef typename image_mesh_type::element_type image_geoelement_type;
    typedef typename image_mesh_type::element_iterator image_mesh_element_iterator;

    // geometric mapping context
    typedef typename image_mesh_type::gm_type image_gm_type;
    typedef typename image_mesh_type::gm_ptrtype image_gm_ptrtype;
    typedef typename image_mesh_type::template gmc<vm::POINT>::type image_gmc_type;
    typedef typename image_mesh_type::template gmc<vm::POINT>::ptrtype image_gmc_ptrtype;

    // dof
    typedef typename dual_image_space_type::dof_type dof_type;

    // basis
    typedef typename dual_image_space_type::basis_type image_basis_type;
    typedef typename domain_space_type::basis_type domain_basis_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     */
    OperatorInterpolation() : super() {}

    /**
     * Construction the global interpolation operator from \p
     * domainspace to \p imagespace and represent it in matrix form
     * using the backend \p backend
     */
    OperatorInterpolation( domain_space_ptrtype const& domainspace,
                           dual_image_space_ptrtype const& imagespace,
                           backend_ptrtype const& backend );

    /**
     * copy constructor
     */
    OperatorInterpolation( OperatorInterpolation const & oi ) : super( oi ) {}

    ~OperatorInterpolation() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

    virtual void update();

private:



};

template<typename DomainSpaceType, typename ImageSpaceType>
OperatorInterpolation<DomainSpaceType, ImageSpaceType>::OperatorInterpolation( domain_space_ptrtype const& domainspace,
                                                                               dual_image_space_ptrtype const& imagespace,
                                                                               backend_ptrtype const& backend )
    :
    super( domainspace, imagespace, backend )
{
    update();
}
template<typename DomainSpaceType, typename ImageSpaceType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType>::update()
{
    if ( this->dualImageSpace()->mesh()->numElements() == 0 )
        return;
#if 0
    /**
     * first check if the mesh are the same as well as the basis, if
     * the same the matrix is the identity
     */
    const bool same_basis = boost::is_same<image_basis_type, domain_basis_type>::value;
    const bool same_mesh = this->domainSpace()->mesh() == this->dualImageSpace()->mesh();
    Debug( 5034 ) << "[interpolate] are the basis the same " << same_basis << "\n";
    Debug( 5034 ) << "[interpolate] are the meshes the same " << same_mesh << "\n";

    if ( same_basis && same_mesh )
        {
            Debug( 5034 ) << "[interpolate] Same mesh and same space\n";
            // construct identity matrix
            return;
        }
#endif
    dof_type const* imagedof = this->dualImageSpace()->dof().get();
    typedef typename domain_space_type::dof_type domain_dof_type;
    domain_dof_type const* domaindof = this->domainSpace()->dof().get();
    image_basis_type const* imagebasis = this->dualImageSpace()->basis().get();
    domain_basis_type const* domainbasis = this->domainSpace()->basis().get();

    int nComponents = this->dualImageSpace()->dof()->nComponents;
    ublas::matrix<value_type,ublas::row_major> Mloc( this->dualImageSpace()->dof()->nDofPerElement()*nComponents,
                                                     this->domainSpace()->dof()->nDofPerElement()*nComponents );

    // Local assembly: compute the Mloc matrix by evaluating
    // the domain space basis function at the dual image space
    // dof points (nodal basis) since we have only computation
    // in the ref elements and the basis and dof points in ref
    // element are the same, we compute Mloc outside the
    // element loop.
    Mloc = ublas::trans( domainbasis->evaluate( imagebasis->dual().points() ) );

    Debug( 5034 ) << "domain ndof = " <<  this->domainSpace()->nDof() << "\n";
    Debug( 5034 ) << "image ndof = " <<  this->dualImageSpace()->nDof() << "\n";

    this->mat().setInitialized( true );
    this->mat().init( this->dualImageSpace()->nDof(), this->domainSpace()->nDof(), 0, 0, 0, 0 );

    Debug( 5034 ) << "matrix size1 = " <<  this->mat().size1() << "\n";
    Debug( 5034 ) << "matrix size2 = " <<  this->mat().size2() << "\n";

    LIFE_ASSERT( this->mat().isInitialized() ).warn( "[OperatorInterpolation] matrix not initialized" );

    // if same mesh but not same function space (e.g. different polynomial order, different basis)
    if ( this->dualImageSpace()->mesh().get() == (image_mesh_type*)this->domainSpace()->mesh().get() )
        {
            Debug( 5034 ) << "[interpolate] Same mesh but not same space\n";

            image_mesh_element_iterator it = this->dualImageSpace()->mesh()->beginElementWithProcessId( Application::processId() );
            image_mesh_element_iterator en = this->dualImageSpace()->mesh()->endElementWithProcessId( Application::processId() );
            for( ; it != en; ++ it )
                {

                    // Global assembly
                    //std::cout << "interpfunc :  " << interpfunc << "\n";
                    for ( uint16_type iloc = 0; iloc < image_basis_type::nLocalDof; ++iloc )
                        {
                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                {
                                    for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                                        {
                                            size_type i =  boost::get<0>(imagedof->localToGlobal( it->id(), iloc, comp ));
                                            size_type j =  boost::get<0>(domaindof->localToGlobal( it->id(), jloc, comp ));
                                            value_type v = Mloc( image_basis_type::nLocalDof*comp + iloc, domain_basis_type::nLocalDof*comp + jloc );
                                            std::cout << "value ( " << i << ", " << j << ")=" << v << "\n";
                                            this->mat().set( i, j, v );
                                        }
                                }
                        }
                }
            Debug( 5034 ) << "[interpolate] Same mesh but not same space done\n";
        } // same mesh
    else
        {
            Debug( 5034 ) << "[interpolate] different meshes\n";
            typename domain_mesh_type::Inverse meshinv( this->domainSpace()->mesh() );

            /* initialisation of the mesh::inverse data structure */
            typename dual_image_space_type::dof_type::dof_points_const_iterator it_dofpt = this->dualImageSpace()->dof()->dofPointBegin();
            typename dual_image_space_type::dof_type::dof_points_const_iterator en_dofpt = this->dualImageSpace()->dof()->dofPointEnd();
            size_type nbpts = 0;
            for( ; it_dofpt != en_dofpt; ++it_dofpt, ++nbpts  )
                {
                    meshinv.addPointWithId( *it_dofpt );
                }
            LIFE_ASSERT( meshinv.nPoints() == nbpts )( meshinv.nPoints() )( nbpts ).error( "invalid number of points " );
            meshinv.distribute();
            Debug( 5034 ) << "number of points inserted : " << nbpts << "\n";

            std::vector<bool> dof_done( nbpts );
            std::fill( dof_done.begin(), dof_done.end(), false );
            std::vector<boost::tuple<size_type,uint16_type > > itab;
            typename matrix_node<value_type>::type pts( image_mesh_type::nDim, 1 );
            domain_mesh_element_iterator it = this->domainSpace()->mesh()->beginElementWithProcessId( Application::processId() );
            domain_mesh_element_iterator en = this->domainSpace()->mesh()->endElementWithProcessId( Application::processId() );

            //size_type ndofcomp = this->dualImageSpace()->nLocalDof()/image_basis_type::nComponents;
            size_type first_dof = this->dualImageSpace()->dof()->firstDof();
            for( ; it != en; ++ it )
                {
                    meshinv.pointsInConvex( it->id(), itab );
                    if (itab.size() == 0)
                        continue;

                    for (size_type i = 0; i < itab.size(); ++i)
                        {
                            // get dof id in target dof table
                            size_type gdof;
                            uint16_type  comp;
                            boost::tie( gdof, comp ) = itab[i];

#if !defined( NDEBUG )
                            Debug( 5034 ) << "[OperatorInterpolation] element : " << it->id() << " npts: " << itab.size() << " ptid: " << i
                                          << " gdof: " << gdof  << " comp: " << comp << "\n";
#endif
                            if ( !dof_done[gdof-first_dof] )
                                {
                                    dof_done[gdof-first_dof]=true;
                                    ublas::column( pts, 0 ) = meshinv.referenceCoords()[gdof];
                                    Mloc = ublas::trans( domainbasis->evaluate( pts ) );

                                    for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                        {
                                            for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                                                {
                                                    size_type globaldof =  gdof;//first_dof+ndofcomp*comp+gdof;
                                                    size_type i =  globaldof;
                                                    size_type j =  boost::get<0>(domaindof->localToGlobal( it->id(),
                                                                                                           jloc,
                                                                                                           comp ));
                                                    //value_type v = Mloc( 0, domain_basis_type::nLocalDof*comp + jloc );
                                                    value_type v = Mloc( 0, nComponents*domain_basis_type::nLocalDof*comp + nComponents*jloc + comp );
                                                    this->mat().set( i, j, v );
                                                }
                                        }
                                }
                        }
                }

            for( size_type i = 0; i < dof_done.size(); ++i )
                {
                    LIFE_ASSERT( dof_done[i] == true )( i ).warn ( "invalid dof, was not treated" );
                }

        }
    Debug( 5034 ) << "[OperatorInterpolation] closing interpolation matrix\n";
    //this->mat().close();
    //this->mat().printMatlab( "Interp.m" );
}
} // Life
#endif /* __OperatorInterpolation_H */
