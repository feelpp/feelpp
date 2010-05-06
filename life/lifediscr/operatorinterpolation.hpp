/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-01

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
    typedef typename dual_image_space_type::dof_type dual_image_dof_type;
    domain_dof_type const* domaindof = this->domainSpace()->dof().get();
    dual_image_dof_type const* dualimagedof = this->dualImageSpace()->dof().get();
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

#if 0
    this->mat().setInitialized( true );
    this->mat().init( this->dualImageSpace()->nDof(), this->domainSpace()->nDof(), 0, 0, 0, 0 );
#endif

    Debug( 5034 ) << "matrix size1 = " <<  this->mat().size1() << "\n";
    Debug( 5034 ) << "matrix size2 = " <<  this->mat().size2() << "\n";

    LIFE_ASSERT( this->mat().isInitialized() ).warn( "[OperatorInterpolation] matrix not initialized" );

    // if same mesh but not same function space (e.g. different polynomial order, different basis)
    if ( this->dualImageSpace()->mesh().get() == (image_mesh_type*)this->domainSpace()->mesh().get() )
        {
            Debug( 5034 ) << "[interpolate] Same mesh but not same space\n";

            image_mesh_element_iterator it = this->dualImageSpace()->mesh()->beginElementWithProcessId( this->dualImageSpace()->mesh()->comm().rank() );
            image_mesh_element_iterator en = this->dualImageSpace()->mesh()->endElementWithProcessId( this->dualImageSpace()->mesh()->comm().rank() );
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
                                            //std::cout << "value ( " << i << ", " << j << ")=" << v << "\n";
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
#if 1
            typedef typename domain_mesh_type::Localization::localization_ptrtype localization_ptrtype;
            typedef typename domain_mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
            typedef typename domain_mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

            typename dual_image_space_type::dof_type::dof_points_const_iterator it_dofpt = this->dualImageSpace()->dof()->dofPointBegin();
            typename dual_image_space_type::dof_type::dof_points_const_iterator en_dofpt = this->dualImageSpace()->dof()->dofPointEnd();

            //init the localization tool
            this->domainSpace()->mesh()->tool_localization()->updateForUse();

            //image_mesh_element_iterator it = this->dualImageSpace()->mesh()->beginElementWithProcessId( this->dualImageSpace()->mesh()->comm().rank() );
            //image_mesh_element_iterator en = this->dualImageSpace()->mesh()->endElementWithProcessId( this->dualImageSpace()->mesh()->comm().rank() );

            //enregistre les numero locales par rapport à la place dans __ptsReal(colonne)
            //std::vector<boost::tuple<uint,size_type> > __memLocDof(image_basis_type::nLocalDof);

            typename matrix_node<value_type>::type __ptsReal( image_mesh_type::nDim, 1/*image_basis_type::nLocalDof*/);
            typename matrix_node<value_type>::type ptsRef(image_mesh_type::nDim , /*nbPtsElt*/1 );
            typename matrix_node<value_type>::type MlocEval(image_mesh_type::nDim,1);
            analysis_iterator_type itanal,itanal_end;
            analysis_output_iterator_type itL,itL_end;

            // create analysys map : id -> List of pt
            localization_ptrtype __loc = this->domainSpace()->mesh()->tool_localization();
            __loc->kdtree()->nbNearNeighbor(3);

            for( ; it_dofpt != en_dofpt; ++ it_dofpt )
                {
                    size_type gdof = boost::get<1>(*it_dofpt);// boost::get<0>(this->dualImageSpace()->dof()->localToGlobal( it->id(), i, 0 ));// 0 est le numero de composantes
                    ublas::column(__ptsReal,0 )= boost::get<0>(*it_dofpt);

                    __loc->run_analysis(__ptsReal);
                    itanal = __loc->result_analysis_begin();
                    itanal_end = __loc->result_analysis_end();

                    for ( ;itanal!=itanal_end;++itanal)
                        {
                            //iterate in the list pt for a element
                            itL=itanal->second.begin();
                            itL_end=itanal->second.end();

                            //compute a point matrix with the list of point
                            ublas::column( ptsRef, 0 ) = boost::get<1>(*itL);

                            /*typename matrix_node<value_type>::type*/ MlocEval = domainbasis->evaluate( ptsRef );

                            itL=itanal->second.begin();
                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                {
                                    size_type j =  boost::get<0>(domaindof->localToGlobal( itanal->first,jloc, 0/*comp*/ ));
                                    value_type v = MlocEval( jloc, 0 );
                                    this->matPtr()->set( gdof, j, v );
                                }
                        }
                }



#elif 0
            typedef typename domain_mesh_type::Localization::localization_ptrtype localization_ptrtype;
            typedef typename domain_mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
            typedef typename domain_mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

            typename dual_image_space_type::dof_type::dof_points_const_iterator it_dofpt = this->dualImageSpace()->dof()->dofPointBegin();
            typename dual_image_space_type::dof_type::dof_points_const_iterator en_dofpt = this->dualImageSpace()->dof()->dofPointEnd();

            uint nbpts = std::distance(it_dofpt,en_dofpt);
            std::vector<bool> dof_done( nbpts );
            std::fill( dof_done.begin(), dof_done.end(), false );

            //init the localization tool
            this->domainSpace()->mesh()->tool_localization()->updateForUse();

            image_mesh_element_iterator it = this->dualImageSpace()->mesh()->beginElementWithProcessId( this->dualImageSpace()->mesh()->comm().rank() );
            image_mesh_element_iterator en = this->dualImageSpace()->mesh()->endElementWithProcessId( this->dualImageSpace()->mesh()->comm().rank() );

            //enregistre les numero locales par rapport à la place dans __ptsReal(colonne)
            std::vector<boost::tuple<uint,size_type> > __memLocDof(image_basis_type::nLocalDof);

            typename matrix_node<value_type>::type __ptsReal( image_mesh_type::nDim, image_basis_type::nLocalDof);

            // create analysys map : id -> List of pt
            localization_ptrtype __loc = this->domainSpace()->mesh()->tool_localization();
            __loc->kdtree()->nbNearNeighbor(3);

            //__loc->run_analysis(__ptsReal);
            for( ; it != en; ++ it )
                {
                    //contient les points real que l'on souhaite evaluer sur un element
                    //typename matrix_node<value_type>::type __ptsReal( image_mesh_type::nDim, image_basis_type::nLocalDof);
                    uint cpt=0;
                    for (uint i=0;i<image_basis_type::nLocalDof;++i)
                        {
                            size_type gdof = boost::get<0>(this->dualImageSpace()->dof()->localToGlobal( it->id(), i, 0 ));// 0 est le numero de composantes
                            ublas::column(__ptsReal,i ) = (boost::get<0>(this->dualImageSpace()->dof()->dofPoint(gdof)));
                            //if ( true/*!dof_done[gdof]*/ ) {
                            //dof_done[gdof]= true;
                            //__memLocDof[cpt]= boost::make_tuple(i,gdof);//i is the num local
                            //ublas::column(__ptsReal,/*i*/cpt ) = (boost::get<0>(this->dualImageSpace()->dof()->dofPoint(gdof)));
                            // ++cpt;
                            //}
                        }

                    //__ptsReal.resize(image_mesh_type::nDim,cpt);
                    /*typename matrix_node<value_type>::type __ptsReal( image_mesh_type::nDim, cpt);
                      for (uint i=0;i<cpt;++i){
                      ublas::column(__ptsReal,i ) = (boost::get<0>(this->dualImageSpace()->dof()->dofPoint( boost::get<1>(__memLocDof[i]))));
                      }*/

                    __loc->run_analysis(__ptsReal);
                    analysis_iterator_type itanal = __loc->result_analysis_begin();
                    analysis_iterator_type itanal_end = __loc->result_analysis_end();
                    analysis_output_iterator_type itL,itL_end;

                    //alocate a point matrix
                    size_type nbPtsElt=itanal->second.size();
                    uint nbCoord=boost::get<1>(*(itanal->second.begin())).size();
                    typename matrix_node<value_type>::type ptsRef( nbCoord, nbPtsElt );

                    for ( ;itanal!=itanal_end;++itanal)
                        {
                            nbPtsElt = itanal->second.size();

                            //iterate in the list pt for a element
                            itL=itanal->second.begin();
                            itL_end=itanal->second.end();

                            //compute a point matrix with the list of point
                            ptsRef= typename matrix_node<value_type>::type( nbCoord, nbPtsElt );
                            for (size_type q=0;q<nbPtsElt;++q,++itL)
                                ublas::column( ptsRef, q ) = boost::get<1>(*itL);

                            //Mloc = ublas::trans( domainbasis->evaluate( pts ) );
                            typename matrix_node<value_type>::type Mloc = /*ublas::trans*/( domainbasis->evaluate( ptsRef ) );

                            itL=itanal->second.begin();
                            for ( uint16_type k1 = 0; k1 < Mloc.size2(); ++k1,++itL )
                                {
                                    size_type i =  boost::get<0>(dualimagedof->localToGlobal( it->id(),  boost::get<0>(*itL) , 0/*comp*/ ));
                                    //size_type i =  boost::get<0>(dualimagedof->localToGlobal( it->id(), __memLocDof[ boost::get<0>(*itL)] , 0/*comp*/ ));
                                    //size_type i =  boost::get<1>(__memLocDof[ boost::get<0>(*itL)]);
                                    for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                        {
                                            size_type j =  boost::get<0>(domaindof->localToGlobal( itanal->first,jloc, 0/*comp*/ ));
                                            //value_type v = Mloc( nComponents*domain_basis_type::nLocalDof*comp + nComponents*jloc + comp,k1 );
                                            value_type v = Mloc( jloc,k1 );
                                            this->mat().set( i, j, v );
                                        }
                                }
                        }
                }


#elif 0 //old version
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
            domain_mesh_element_iterator it = this->domainSpace()->mesh()->beginElementWithProcessId( this->dualImageSpace()->mesh()->comm().rank() );
            domain_mesh_element_iterator en = this->domainSpace()->mesh()->endElementWithProcessId( this->dualImageSpace()->mesh()->comm().rank() );

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

#endif

        }
    Debug( 5034 ) << "[OperatorInterpolation] closing interpolation matrix\n";
    //this->mat().close();
    //this->mat().printMatlab( "Interp.m" );
}
} // Life
#endif /* __OperatorInterpolation_H */
