/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-01

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
   \author Chabannes Vincent <vincent.chabannes@imag.fr>
   \date 2008-02-01
 */
#ifndef __OperatorInterpolation_H
#define __OperatorInterpolation_H 1

#include <feel/feeldiscr/operatorlinear.hpp>

namespace Feel
{

struct InterpolationNonConforme
{
    static const bool value=false;
};

struct InterpolationConforme
{
    static const bool value=true;
};


namespace detailsup
{

template < typename EltType >
size_type
idElt( EltType & elt,mpl::size_t<MESH_ELEMENTS> )
{
    return elt.id();
}

template < typename EltType >
size_type
idElt( EltType & elt,mpl::size_t<MESH_FACES> )
{
    return elt.element0().id();
}


} //detailsup


/**
 * \class OperatorInterpolation
 * \brief Global interpolation operator
 *
 * @author Christophe Prud'homme
 * @see
 */

template<typename DomainSpaceType,
         typename ImageSpaceType,
         typename IteratorRange = boost::tuple<mpl::size_t<MESH_ELEMENTS>,
         typename MeshTraits<typename ImageSpaceType::mesh_type>::element_const_iterator,
         typename MeshTraits<typename ImageSpaceType::mesh_type>::element_const_iterator>,
         typename InterpType = InterpolationNonConforme >
class OperatorInterpolation : public OperatorLinear<DomainSpaceType, ImageSpaceType >
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

    typedef typename boost::tuples::template element<0, IteratorRange>::type idim_type;
    typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_type;
    typedef IteratorRange range_iterator;


    static const uint16_type nLocalDofInDualImageElt = mpl::if_< mpl::equal_to< idim_type ,mpl::size_t<MESH_ELEMENTS> >,
                             mpl::int_< image_basis_type::nLocalDof > ,
                             mpl::int_< image_mesh_type::face_type::numVertices*dual_image_space_type::fe_type::nDofPerVertex +
                             image_mesh_type::face_type::numEdges*dual_image_space_type::fe_type::nDofPerEdge +
                             image_mesh_type::face_type::numFaces*dual_image_space_type::fe_type::nDofPerFace > >::type::value;

    // type conforme or non conforme
    typedef InterpType interpolation_type;

    // matrix graph
    typedef GraphCSR graph_type;
    typedef boost::shared_ptr<graph_type> graph_ptrtype;


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

    OperatorInterpolation( domain_space_ptrtype const& domainspace,
                           dual_image_space_ptrtype const& imagespace,
                           IteratorRange const& r,
                           backend_ptrtype const& backend );


    /**
     * copy constructor
     */
    OperatorInterpolation( OperatorInterpolation const & oi )
        :
        super( oi ),
        _M_range( oi._M_range )
    {}

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

    void updateSameMesh();
    void updateNoRelationMesh();
    void updateNoRelationMeshMPI();


    range_iterator _M_range;

};

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::OperatorInterpolation( domain_space_ptrtype const& domainspace,
        dual_image_space_ptrtype const& imagespace,
        backend_ptrtype const& backend )
    :
    super( domainspace, imagespace, backend, false ),
    _M_range( elements( imagespace->mesh() ) )
{
    update();
}


template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::OperatorInterpolation( domain_space_ptrtype const& domainspace,
        dual_image_space_ptrtype const& imagespace,
        IteratorRange const& r,
        backend_ptrtype const& backend )
    :
    super( domainspace, imagespace, backend, false ),
    _M_range( r )
{
    update();
}


template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::update()
{
    if ( this->dualImageSpace()->mesh()->numElements() == 0 )
        return;

    // if same mesh but not same function space (e.g. different polynomial order, different basis)
    if ( this->dualImageSpace()->mesh().get() == ( image_mesh_type* )this->domainSpace()->mesh().get() )
    {
        this->updateSameMesh();
    }

    else // no relation between meshes
    {
        if ( this->dualImageSpace()->worldsComm()[0].localSize() > 1 )
            this->updateNoRelationMeshMPI();

        else
            this->updateNoRelationMesh();
    }

    // close matrix after build
    this->mat().close();
}

//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::updateSameMesh()
{
#if !defined(FEELPP_ENABLE_MPI_MODE) // NOT MPI
    const size_type proc_id           = this->dualImageSpace()->mesh()->comm().rank();
    const size_type nrow_dof_on_proc    = this->dualImageSpace()->nLocalDof();
    const size_type firstrow_dof_on_proc = this->dualImageSpace()->dof()->firstDof( proc_id );
    const size_type lastrow_dof_on_proc = this->dualImageSpace()->dof()->lastDof( proc_id );
    const size_type firstcol_dof_on_proc = this->domainSpace()->dof()->firstDof( proc_id );
    const size_type lastcol_dof_on_proc = this->domainSpace()->dof()->lastDof( proc_id );
#else
    const size_type proc_id           = this->dualImageSpace()->worldsComm()[0].localRank();
    const size_type nrow_dof_on_proc    = this->dualImageSpace()->nLocalDof();
    const size_type firstrow_dof_on_proc = this->dualImageSpace()->dof()->firstDofGlobalCluster( proc_id );
    const size_type lastrow_dof_on_proc = this->dualImageSpace()->dof()->lastDofGlobalCluster( proc_id );
    const size_type firstcol_dof_on_proc = this->domainSpace()->dof()->firstDofGlobalCluster( proc_id );
    const size_type lastcol_dof_on_proc = this->domainSpace()->dof()->lastDofGlobalCluster( proc_id );
#endif
    graph_ptrtype sparsity_graph( new graph_type( nrow_dof_on_proc,
                                  firstrow_dof_on_proc, lastrow_dof_on_proc,
                                  firstcol_dof_on_proc, lastcol_dof_on_proc ) );

    auto const* imagedof = this->dualImageSpace()->dof().get();
    auto const* domaindof = this->domainSpace()->dof().get();
    auto const* imagebasis = this->dualImageSpace()->basis().get();
    auto const* domainbasis = this->domainSpace()->basis().get();


    std::vector<bool> dof_done( nrow_dof_on_proc, false );
    std::vector< std::list<std::pair<size_type,double> > > memory_valueInMatrix( nrow_dof_on_proc );

    // Local assembly: compute the Mloc matrix by evaluating
    // the domain space basis function at the dual image space
    // dof points (nodal basis) since we have only computation
    // in the ref elements and the basis and dof points in ref
    // element are the same, we compute Mloc outside the
    // element loop.

    //typename matrix_node<value_type>::type Mloc(domain_basis_type::nLocalDof*domain_basis_type::nComponents1,1);
    auto const& Mloc = domainbasis->evaluate( imagebasis->dual().points() );

    Debug( 5034 ) << "[interpolate] Same mesh but not same space\n";

    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = _M_range;

    for ( ; it != en; ++ it )
    {
        auto idElem = detailsup::idElt( *it,idim_type() );

        // Global assembly
        for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
        {
            for ( uint16_type comp = 0; comp < image_basis_type::nComponents; ++comp )
            {
                size_type i =  boost::get<0>( imagedof->localToGlobal( *it, iloc, comp ) );

                if ( !dof_done[i] )
                {
#if !defined(FEELPP_ENABLE_MPI_MODE) // NOT MPI
                    const auto ig1 = i;
                    const auto theproc = imagedof->worldComm().localRank();
#else // WITH MPI
                    const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[i];
                    const auto theproc = imagedof->procOnGlobalCluster( ig1 );
#endif
                    auto& row = sparsity_graph->row( ig1 );
                    row.get<0>() = theproc;
                    const size_type il1 = ig1 - imagedof->firstDofGlobalCluster( theproc );
                    row.get<1>() = il1;
                    //row.get<1>() = i;

                    uint16_type ilocprime=imagedof->localDofInElement( *it, iloc, comp ) ;

                    for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                    {
                        // get column
                        const size_type j =  boost::get<0>( domaindof->localToGlobal( idElem, jloc, comp ) );
                        //up the pattern graph
#if !defined(FEELPP_ENABLE_MPI_MODE) // NOT MPI
                        row.get<2>().insert( j );
#else // WITH MPI
                        row.get<2>().insert( domaindof->mapGlobalProcessToGlobalCluster()[j] );
#endif
                        // get interpolated value
                        const value_type v = Mloc( domain_basis_type::nComponents1*jloc +
                                                   comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof +
                                                   comp,
                                                   ilocprime );
                        // save in matrux
                        memory_valueInMatrix[i].push_back( std::make_pair( j,v ) );
                    }

                    dof_done[i]=true;
                }
            }
        }
    }

    //-----------------------------------------
    // compute graph
    sparsity_graph->close();
    //-----------------------------------------
    // create matrix
    this->matPtr() = this->backend()->newMatrix( this->domainSpace()->mapOnOff(),
                     this->dualImageSpace()->mapOn(),
                     sparsity_graph  );
    //-----------------------------------------

    // assemble matrix
    for ( size_type idx_i=0 ; idx_i<nrow_dof_on_proc; ++idx_i )
    {
        for ( auto it_j=memory_valueInMatrix[idx_i].begin(),en_j=memory_valueInMatrix[idx_i].end() ; it_j!=en_j ; ++it_j )
        {
            this->matPtr()->set( idx_i,it_j->first,it_j->second );
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::updateNoRelationMesh()
{
    typedef typename matrix_node<typename image_mesh_type::value_type>::type matrix_node_type;

    Debug( 5034 ) << "[interpolate] different meshes\n";

    const size_type proc_id           = this->dualImageSpace()->mesh()->comm().rank();
    const size_type n1_dof_on_proc    = this->dualImageSpace()->nLocalDof();
    const size_type firstrow_dof_on_proc = this->dualImageSpace()->dof()->firstDof( proc_id );
    const size_type lastrow_dof_on_proc = this->dualImageSpace()->dof()->lastDof( proc_id );
    const size_type firstcol_dof_on_proc = this->domainSpace()->dof()->firstDof( proc_id );
    const size_type lastcol_dof_on_proc = this->domainSpace()->dof()->lastDof( proc_id );

    graph_ptrtype sparsity_graph( new graph_type( n1_dof_on_proc,
                                  firstrow_dof_on_proc, lastrow_dof_on_proc,
                                  firstcol_dof_on_proc, lastcol_dof_on_proc ) );

    auto const* imagedof = this->dualImageSpace()->dof().get();
    auto const* domaindof = this->domainSpace()->dof().get();
    auto const* domainbasis = this->domainSpace()->basis().get();

    //-----------------------------------------
    //init the localization tool
    auto locTool = this->domainSpace()->mesh()->tool_localization();
    locTool->updateForUse();
    //locTool->kdtree()->nbNearNeighbor(3);
    //locTool->kdtree()->nbNearNeighbor(this->domainSpace()->mesh()->numElements());
    //locTool->setExtrapolation(false);

    //-----------------------------------------
    // usefull data
    matrix_node_type ptsReal( image_mesh_type::nRealDim, 1 );
    matrix_node_type ptsRef( image_mesh_type::nRealDim , 1 );
    typename domain_mesh_type::Localization::container_search_iterator_type itanal,itanal_end;
    typename domain_mesh_type::Localization::container_output_iterator_type itL,itL_end;
    matrix_node_type MlocEval( domain_basis_type::nLocalDof*domain_basis_type::nComponents1,1 );

    std::vector<bool> dof_done( this->dualImageSpace()->nDof(), false );
    std::vector< std::list<std::pair<size_type,double> > > memory_valueInMatrix( this->dualImageSpace()->nDof() );

    //-----------------------------------------
    // for each element in range
    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = _M_range;
    size_type eltIdLocalised = 0;

    for ( ; it != en; ++ it )
        {
            for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
                {
                    for ( uint16_type comp = 0; comp < image_basis_type::nComponents; ++comp )
                        {
                            const auto gdof =  boost::get<0>(imagedof->localToGlobal( *it, iloc, comp ));
                            if (!dof_done[gdof])
                                {
                                    //------------------------
                                    // get the graph row
#if !defined(FEELPP_ENABLE_MPI_MODE) // NOT MPI
                                    const auto ig1 = gdof;
                                    const auto theproc = imagedof->worldComm().localRank();
#else // WITH MPI
                                    const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[gdof];
                                    const auto theproc = imagedof->procOnGlobalCluster( ig1 );
#endif
                                    auto& row = sparsity_graph->row(ig1);
                                    row.get<0>() = theproc;
                                    row.get<1>() = gdof;
                                    //------------------------
                                    // the dof point
                                    ublas::column(ptsReal,0 ) = boost::get<0>(imagedof->dofPoint(gdof));
                                    //------------------------
                                    // localisation process
                                    eltIdLocalised = locTool->run_analysis(ptsReal,eltIdLocalised,it->vertices()/*it->G()*/,mpl::bool_<interpolation_type::value>()).get<1>();
                                    //------------------------
                                    // for each localised points
                                    itanal = locTool->result_analysis_begin();
                                    itanal_end = locTool->result_analysis_end();
                                    for ( ;itanal!=itanal_end;++itanal)
                                        {
                                            itL=itanal->second.begin();
                                            ublas::column( ptsRef, 0 ) = boost::get<1>( *itL );

                                            MlocEval = domainbasis->evaluate( ptsRef );

                                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                {
                                                    //get global dof
                                                    size_type j =  boost::get<0>( domaindof->localToGlobal( itanal->first,jloc,comp ) );
                                                    value_type v = MlocEval( domain_basis_type::nComponents1*jloc
                                                                             + comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof
                                                                             + comp,
                                                                             0 );
#if !defined(FEELPP_ENABLE_MPI_MODE) // NOT MPI
                                                    row.get<2>().insert( j );
#else // WITH MPI
                                                    row.get<2>().insert( domaindof->mapGlobalProcessToGlobalCluster()[j] );
#endif
                                                    memory_valueInMatrix[gdof].push_back( std::make_pair( j,v ) );
                                                }
                                        }
                                    dof_done[gdof]=true;
                                } // if (!dof_done[gdof])
                        } //  for ( uint16_type comp = 0; comp < image_basis_type::nComponents; ++comp )
                } // for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
        } // for( ; it != en; ++ it )

    //-----------------------------------------
    // compute graph
    sparsity_graph->close();//sparsity_graph->printPython("mygraphpython.py");
    //-----------------------------------------
    // create matrix
    this->matPtr() = this->backend()->newMatrix(this->dualImageSpace()->nDof(), this->domainSpace()->nDof() ,
                                               this->dualImageSpace()->nLocalDof(), this->domainSpace()->nLocalDof(),
                                               sparsity_graph );
    //-----------------------------------------
    // assemble matrix
    for (size_type idx_i=0 ; idx_i<this->dualImageSpace()->nDof() ;++idx_i)
        {
            for (auto it_j=memory_valueInMatrix[idx_i].begin(),en_j=memory_valueInMatrix[idx_i].end() ; it_j!=en_j ; ++it_j)
                {
                    this->matPtr()->set(idx_i,it_j->first,it_j->second);
                }
        }

}

//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//


template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::updateNoRelationMeshMPI()
{
#if defined(FEELPP_ENABLE_MPI_MODE) // WITH MPI

    typedef typename matrix_node<typename image_mesh_type::value_type>::type matrix_node_type;

    std::cout << "updateNoRelationMeshMPI()"<<std::endl;

    const size_type proc_id           = this->dualImageSpace()->worldsComm()[0].localRank();
    const size_type nrow_dof_on_proc    = this->dualImageSpace()->nLocalDof();
    const size_type firstrow_dof_on_proc = this->dualImageSpace()->dof()->firstDofGlobalCluster( proc_id );
    const size_type lastrow_dof_on_proc = this->dualImageSpace()->dof()->lastDofGlobalCluster( proc_id );
    const size_type firstcol_dof_on_proc = this->domainSpace()->dof()->firstDofGlobalCluster( proc_id );
    const size_type lastcol_dof_on_proc = this->domainSpace()->dof()->lastDofGlobalCluster( proc_id );

    graph_ptrtype sparsity_graph( new graph_type( nrow_dof_on_proc,
                                                  firstrow_dof_on_proc, lastrow_dof_on_proc,
                                                  firstcol_dof_on_proc, lastcol_dof_on_proc ) );

    auto const* imagedof = this->dualImageSpace()->dof().get();
    auto const* domaindof = this->domainSpace()->dof().get();
    auto const* imagebasis = this->dualImageSpace()->basis().get();
    auto const* domainbasis = this->domainSpace()->basis().get();

    //-----------------------------------------
    //init the localization tool
    auto locTool = this->domainSpace()->mesh()->tool_localization();
    locTool->updateForUse();
    //locTool->kdtree()->nbNearNeighbor(3);
    //locTool->kdtree()->nbNearNeighbor(this->domainSpace()->mesh()->numElements());
    //locTool->setExtrapolation(false);


    //-----------------------------------------
    // usefull data
    matrix_node_type ptsReal( image_mesh_type::nRealDim, 1);
    matrix_node_type ptsRef(image_mesh_type::nRealDim , 1 );
    typename domain_mesh_type::Localization::container_search_iterator_type itanal,itanal_end;
    typename domain_mesh_type::Localization::container_output_iterator_type itL,itL_end;
    matrix_node_type MlocEval(domain_basis_type::nLocalDof*domain_basis_type::nComponents1,1);

    std::vector<bool> dof_done( this->dualImageSpace()->nLocalDof(), false);
    std::vector<bool> dof_done2( this->dualImageSpace()->nLocalDof(), false);
    std::vector< std::list<std::pair<size_type,double> > > memory_valueInMatrix( this->dualImageSpace()->nLocalDof() );
    //std::vector< std::pair<bool,typename image_mesh_type::node_type > > memory_localisationFail( this->dualImageSpace()->nDof() );
    std::list<boost::tuple<size_type,uint16_type> > memory_localisationFail;// gdof,comp

    size_type eltIdLocalised = this->domainSpace()->mesh()->beginElementWithId(this->domainSpace()->mesh()->worldComm().localRank())->id();

    // WARNING in PARALLELE!!!!!!!!!!!!!!!!
    locTool->setExtrapolation(false);
    //-----------------------------------------
    // for each element in range
    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = _M_range;
    for( ; it != en; ++ it )
        {
            for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
                {
                    for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                        {
                            const auto gdof =  boost::get<0>(imagedof->localToGlobal( *it, iloc, comp ));
                            if (!dof_done[gdof])
                                {
                                    //------------------------
                                    // get the graph row
                                    const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[gdof];
                                    const auto theproc = imagedof->procOnGlobalCluster(ig1);
                                    auto& row = sparsity_graph->row(ig1);
                                    row.get<0>() = theproc;
                                    row.get<1>() = ig1 - imagedof->firstDofGlobalCluster(theproc);//   gdof;
                                    //------------------------
                                    // the dof point
                                    ublas::column(ptsReal,0 ) = imagedof->dofPoint(gdof).get<0>();
                                    //------------------------
                                    // localisation process
                                    auto resLocalisation = locTool->run_analysis(ptsReal,eltIdLocalised,it->vertices()/*it->G()*/,mpl::bool_<interpolation_type::value>());
                                    if (!resLocalisation.get<0>()[0]) // not find
                                        {
                                            memory_localisationFail.push_back(boost::make_tuple(gdof,comp) );
                                        }
                                    else // point found
                                        {
                                            //------------------------
                                            // for each localised points
                                            eltIdLocalised = resLocalisation.get<1>();
                                            itanal = locTool->result_analysis_begin();
                                            itanal_end = locTool->result_analysis_end();
                                            for ( ;itanal!=itanal_end;++itanal)
                                                {
                                                    itL=itanal->second.begin();
                                                    ublas::column( ptsRef, 0 ) = boost::get<1>(*itL);
                                                    // evaluate basis functions for this point
                                                    MlocEval = domainbasis->evaluate( ptsRef );
#if 1
                                                    for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                        {
                                                            //get global dof
                                                            size_type j =  boost::get<0>(domaindof->localToGlobal( itanal->first,jloc,comp ));
                                                            // up graph
                                                            row.get<2>().insert(domaindof->mapGlobalProcessToGlobalCluster()[j]);
                                                            // get value
                                                            value_type v = MlocEval( domain_basis_type::nComponents1*jloc +
                                                                                     comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof +
                                                                                     comp, 0 );
                                                            // save value
                                                            memory_valueInMatrix[gdof].push_back(std::make_pair(j,v));
                                                        }
#endif
                                                }
                                            //dof_done[gdof]=true;
                                        } // else // point found
                                    dof_done[gdof]=true;
                                } // if (!dof_done[gdof])
                        } // for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                }  // for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
        } // for( ; it != en; ++ it )

    //-----------------------------------------------------------------------------------------
    this->dualImageSpace()->mesh()->worldComm().localComm().barrier();
    //-----------------------------------------------------------------------------------------

    std::cout << "opinterp  finish " << this->domainSpace()->worldComm().godRank()<< std::endl;
    std::cout << "GodRank " << this->domainSpace()->worldComm().godRank()
              << " memory_localisationFail.size() " << memory_localisationFail.size()
              << " loc->barycenter() " << locTool->barycenter()
              << std::endl;

    //------------------------------
    // memory map (loc index pt) -> global dofs
    std::vector<size_type> memmapGdof( memory_localisationFail.size() );
    // memory map (loc index pt) -> comp
    std::vector<uint16_type> memmapComp( memory_localisationFail.size() );
    // points to lacalize
    //std::vector<typename image_mesh_type::node_type> pointsSearched( memory_localisationFail.size() );
    std::vector<double> pointsSearched_X( memory_localisationFail.size() );
    std::vector<double> pointsSearched_Y( memory_localisationFail.size() );
    std::vector<double> pointsSearched_Z( memory_localisationFail.size() );

    // build
    size_type cpt=0;
    for(auto it_locFail=memory_localisationFail.begin(), en_locFail=memory_localisationFail.end(); it_locFail!=en_locFail;++it_locFail,++cpt)
        {
            //const auto gdof = it_locFail->get<0>();
            memmapGdof[cpt]=it_locFail->get<0>();
            memmapComp[cpt]=it_locFail->get<1>();
            //pointsSearched[cpt]=imagedof->dofPoint(it_locFail->get<0>()).get<0>();
            pointsSearched_X[cpt]=imagedof->dofPoint(it_locFail->get<0>()).get<0>()(0);
            if (image_mesh_type::nRealDim>1) pointsSearched_Y[cpt]=imagedof->dofPoint(it_locFail->get<0>()).get<0>()(1);
            if (image_mesh_type::nRealDim>2) pointsSearched_Z[cpt]=imagedof->dofPoint(it_locFail->get<0>()).get<0>()(2);
        }

    boost::tie( boost::tuples::ignore, it, en ) = _M_range;

    //std::vector<typename image_mesh_type::node_type> dataToRecv;
    //std::vector<typename image_mesh_type::node_type> dataToSend(2);
    std::vector<size_type> pointsSearchedSizeWorld(this->dualImageSpace()->mesh()->worldComm().localComm().size());

    std::vector<double> dataToRecv_X(1,0);
    std::vector<double> dataToRecv_Y(1,0);
    std::vector<double> dataToRecv_Z(1,0);

    std::vector<double> pointsRefFinded_X(1,0);
    std::vector<double> pointsRefFinded_Y(1,0);
    std::vector<double> pointsRefFinded_Z(1,0);
    std::vector<int> pointsRefIsFinded(1,0);
    std::vector<int> pointsIdEltFinded(1,0);
    std::vector<uint16_type> pointsComp(1,0);
    //------------------------------
    // proc apres proc
    for (int proc=0;proc<this->dualImageSpace()->mesh()->worldComm().localSize();++proc)
        {
            //-----------------------------------------------------------------------------------------
            this->dualImageSpace()->mesh()->worldComm().localComm().barrier();
            //-----------------------------------------------------------------------------------------

            mpi::all_gather( this->dualImageSpace()->mesh()->worldComm().localComm(),
                             pointsSearched_X.size(),
                             pointsSearchedSizeWorld );
            size_type nDataRecv = pointsSearchedSizeWorld[proc];
            std::cout << "GodRank " << this->domainSpace()->worldComm().godRank()
                      << " np pt to send " << pointsSearched_X.size()
                      << " nDataRecv " << nDataRecv
                      << std::endl;

            const int tag_X = 0, tag_Y = 1, tag_Z = 2, tag_IsFind = 3, tag_IdElt = 4;


            /*if (proc_id ==proc)
                for (int ii=0;ii<pointsSearched_X.size();++ii)
                std::cout << "(" << pointsSearched_X[ii] << "," << pointsSearched_Y[ii] << ")";*/

            std::cout << std::endl;
            this->dualImageSpace()->mesh()->worldComm().localComm().barrier();


#if 1
            if (proc_id ==proc)
                {
                    int rankToSend = 0;if (proc==0) rankToSend=1; else rankToSend=0;
                    this->dualImageSpace()->mesh()->worldComm().localComm().send(rankToSend,tag_X,pointsSearched_X);
                    if (image_mesh_type::nRealDim>1) this->dualImageSpace()->mesh()->worldComm().localComm().send(rankToSend,tag_Y,pointsSearched_Y);
                    if (image_mesh_type::nRealDim>2) this->dualImageSpace()->mesh()->worldComm().localComm().send(rankToSend,tag_Z,pointsSearched_Z);
                }
            else
                {
                    int rankToRecv = proc;
                    this->dualImageSpace()->mesh()->worldComm().localComm().recv(rankToRecv,tag_X,dataToRecv_X);
                    if (image_mesh_type::nRealDim>1) this->dualImageSpace()->mesh()->worldComm().localComm().recv(rankToRecv,tag_Y,dataToRecv_Y);
                    if (image_mesh_type::nRealDim>2) this->dualImageSpace()->mesh()->worldComm().localComm().recv(rankToRecv,tag_Z,dataToRecv_Z);
                    /*for (int ii=0;ii<dataToRecv_X.size();++ii)
                      std::cout << "(" << dataToRecv_X[ii] << "," << dataToRecv_Y[ii] << ")";*/
                }
#endif

            //-----------------------------------------------------------------------------------------
            std::cout << std::endl;
            this->dualImageSpace()->mesh()->worldComm().localComm().barrier();
            //-----------------------------------------------------------------------------------------
#if 1
            pointsRefFinded_X.resize(nDataRecv);
            if (image_mesh_type::nRealDim>1) pointsRefFinded_Y.resize(nDataRecv);
            if (image_mesh_type::nRealDim>2) pointsRefFinded_Z.resize(nDataRecv);
            pointsRefIsFinded.resize(nDataRecv,0);std::fill(pointsRefIsFinded.begin(),pointsRefIsFinded.end(),0);
            pointsIdEltFinded.resize(nDataRecv);

            if (proc_id !=proc)
                {
                    for (size_type k=0;k<dataToRecv_X.size();++k)
                        {
                            ublas::column(ptsReal,0)(0) = dataToRecv_X[k];
                            if (image_mesh_type::nRealDim>1) ublas::column(ptsReal,0)(1) = dataToRecv_Y[k];
                            if (image_mesh_type::nRealDim>2) ublas::column(ptsReal,0)(2) = dataToRecv_Z[k];
                            auto resLocalisation = locTool->run_analysis(ptsReal,eltIdLocalised,it->vertices()/*it->G()*/,mpl::bool_<interpolation_type::value>());
                            if (resLocalisation.get<0>()[0])
                                {
                                    eltIdLocalised = resLocalisation.get<1>();
                                    pointsRefIsFinded[k]=1;
                                    pointsIdEltFinded[k]=eltIdLocalised;
                                    if ( this->domainSpace()->mesh()->element(eltIdLocalised).isGhostCell()) std::cout << "BAD";
                                    itanal = locTool->result_analysis_begin();
                                    itanal_end = locTool->result_analysis_end();
                                    for ( ;itanal!=itanal_end;++itanal)
                                        {
                                            //if (eltIdLocalised !=itanal->first) std::cout << "NEW BAD" << std::endl;

                                            itL=itanal->second.begin();
                                            ublas::column( ptsRef, 0 ) = boost::get<1>(*itL);

                                            pointsRefFinded_X[k]=ublas::column( ptsRef, 0 )(0);
                                            if (image_mesh_type::nRealDim>1) pointsRefFinded_Y[k]=ublas::column( ptsRef, 0 )(1);
                                            if (image_mesh_type::nRealDim>2) pointsRefFinded_Z[k]=ublas::column( ptsRef, 0 )(2);
                                        }
                                    std::cout << "F";
                                }
                            else
                                {
                                    std::cout << "NOT FIND"<<std::endl;
                                }
                        }
                }
            std::cout << std::endl;
#endif
            //-----------------------------------------------------------------------------------------
            this->dualImageSpace()->mesh()->worldComm().localComm().barrier();
            //-----------------------------------------------------------------------------------------
#if 1
            if (proc_id !=proc)
                {
                    int rankToSend = proc;
                    this->dualImageSpace()->mesh()->worldComm().localComm().send(rankToSend,tag_X,pointsRefFinded_X);
                    if (image_mesh_type::nRealDim>1) this->dualImageSpace()->mesh()->worldComm().localComm().send(rankToSend,tag_Y,pointsRefFinded_Y);
                    if (image_mesh_type::nRealDim>2) this->dualImageSpace()->mesh()->worldComm().localComm().send(rankToSend,tag_Z,pointsRefFinded_Z);
                    this->dualImageSpace()->mesh()->worldComm().localComm().send(rankToSend,tag_IsFind,pointsRefIsFinded);
                    this->dualImageSpace()->mesh()->worldComm().localComm().send(rankToSend,tag_IdElt,pointsIdEltFinded);
                }
            else
                {
                    int rankToRecv = 0;if (proc==0) rankToRecv=1; else rankToRecv=0;
                    this->dualImageSpace()->mesh()->worldComm().localComm().recv(rankToRecv,tag_X,pointsRefFinded_X);
                    if (image_mesh_type::nRealDim>1) this->dualImageSpace()->mesh()->worldComm().localComm().recv(rankToRecv,tag_Y,pointsRefFinded_Y);
                    if (image_mesh_type::nRealDim>2) this->dualImageSpace()->mesh()->worldComm().localComm().recv(rankToRecv,tag_Z,pointsRefFinded_Z);
                    this->dualImageSpace()->mesh()->worldComm().localComm().recv(rankToRecv,tag_IsFind,pointsRefIsFinded);
                    this->dualImageSpace()->mesh()->worldComm().localComm().recv(rankToRecv,tag_IdElt,pointsIdEltFinded);

                    std::cout  << "GodRank " << this->domainSpace()->worldComm().godRank()
                               << " pointsRefFinded_X.size() " << pointsRefFinded_X.size()
                               << " pointsRefFinded_Y.size() " << pointsRefFinded_Y.size()
                               << " pointsRefFinded_Z.size() " << pointsRefFinded_Z.size()
                               << " pointsRefIsFinded.size() " << pointsRefIsFinded.size()
                               << " pointsIdEltFinded.size() " << pointsIdEltFinded.size()
                               << std::endl;
                }
#endif
            //-----------------------------------------------------------------------------------------
            this->dualImageSpace()->mesh()->worldComm().localComm().barrier();
            //-----------------------------------------------------------------------------------------
#if 1
            if (proc_id==proc)
                {
                    for ( int k=0;k<pointsRefFinded_X.size();++k)
                        {
                            if (pointsRefIsFinded[k]>0)
                                { //std::cout << "T";
                                    ublas::column( ptsRef, 0 )(0) = pointsRefFinded_X[k];
                                    if (image_mesh_type::nRealDim>1) ublas::column( ptsRef, 0 )(1) = pointsRefFinded_Y[k];
                                    if (image_mesh_type::nRealDim>2) ublas::column( ptsRef, 0 )(2) = pointsRefFinded_Z[k];

                                    //std::cout << ptsRef;
                                    MlocEval = domainbasis->evaluate( ptsRef );

                                    const auto i_gdof = memmapGdof[k];
                                    if (!dof_done2[i_gdof])
                                        {
#if 1
                                            const auto comp = memmapComp[k];
                                            const size_type myidElt = pointsIdEltFinded[k];
                                            const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[i_gdof];
                                            const auto theproc = imagedof->procOnGlobalCluster(ig1);
                                            auto& row = sparsity_graph->row(ig1);
                                            row.get<0>() = theproc;
                                            //row.get<1>() = i_gdof;
                                            row.get<1>() = ig1 - imagedof->firstDofGlobalCluster(theproc);

                                            //std::cout << "("<<gdof<<"-"<<myidElt<< "-" << comp << ")";
                                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                {
                                                    //std::cout << "("<<myidElt<< "-" << comp << ")";
                                                    //get global dof
                                                    // L'ERREUR est la! myodElt n'est pas connu !!!!!!
                                                    size_type j_gdof =  boost::get<0>(domaindof->localToGlobal( myidElt,jloc,comp ));
                                                    // up graph
                                                    row.get<2>().insert(domaindof->mapGlobalProcessToGlobalCluster()[j_gdof]);
                                                    // get value
                                                    value_type v = MlocEval( domain_basis_type::nComponents1*jloc +
                                                                             comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof +
                                                                             comp, 0 );
                                                    // save value
                                                    memory_valueInMatrix[i_gdof].push_back(std::make_pair(j_gdof,v));
                                                }
#endif
                                            dof_done2[i_gdof]=true;
                                        }
                                }
                        }
                }

#endif

            //-----------------------------------------------------------------------------------------
            this->dualImageSpace()->mesh()->worldComm().localComm().barrier();
            //-----------------------------------------------------------------------------------------


        } // for (int proc=0;proc<this->dualImageSpace()->mesh()->worldComm().localSize();++proc)



    //-----------------------------------------------------------------------------------------
    this->dualImageSpace()->mesh()->worldComm().localComm().barrier();
    //-----------------------------------------------------------------------------------------
    std::cout << "\nBEFORECOMPUTE GRAPH"<< std::endl;

    //-----------------------------------------
    // compute graph
    //sparsity_graph->zero();
    sparsity_graph->close(); //sparsity_graph->printPython("mygraphpython_2.py");
    //-----------------------------------------
    // create matrix
    this->matPtr() = this->backend()->newMatrix(this->domainSpace()->mapOnOff(),
                                                this->dualImageSpace()->mapOn(),
                                                sparsity_graph  );

   // assemble matrix
    for (size_type idx_i=0 ; idx_i<nrow_dof_on_proc;++idx_i)
        {
            for (auto it_j=memory_valueInMatrix[idx_i].begin(),en_j=memory_valueInMatrix[idx_i].end() ; it_j!=en_j ; ++it_j)
                {
                    this->matPtr()->set(idx_i,it_j->first,it_j->second);
                }
        }

#endif  // WITH MPI
}


//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//


template<typename DomainSpaceType, typename ImageSpaceType, typename IteratorRange, typename InterpType >
boost::shared_ptr<OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType> >
opInterp( boost::shared_ptr<DomainSpaceType> const& domainspace,
          boost::shared_ptr<ImageSpaceType> const& imagespace,
          IteratorRange const& r,
          typename OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::backend_ptrtype const& backend,
          InterpType /**/
        )
{
    typedef OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType> operatorinterpolation_type;

    boost::shared_ptr<operatorinterpolation_type> opI( new operatorinterpolation_type( domainspace,imagespace,r,backend ) );

    return opI;
}



template<typename Args>
struct compute_opInterpolation_return
{
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::domainSpace>::type>::type::element_type domain_space_type;
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::imageSpace>::type>::type::element_type image_space_type;

    typedef typename boost::remove_const<
    typename boost::remove_reference<
    typename parameter::binding<Args,
             tag::range,
             typename OperatorInterpolation<domain_space_type, image_space_type>::range_iterator
             >::type >::type >::type iterator_range_type;

    typedef typename boost::remove_const<
    typename boost::remove_reference<
    typename parameter::binding<Args,
             tag::type,
             InterpolationNonConforme
             >::type >::type >::type interpolation_type;


    typedef boost::shared_ptr<OperatorInterpolation<domain_space_type, image_space_type,iterator_range_type,interpolation_type> > type;
};

BOOST_PARAMETER_FUNCTION(
    ( typename compute_opInterpolation_return<Args>::type ), // 1. return type
    opInterpolation,                        // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required
      ( domainSpace,    *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
      ( imageSpace,     *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
    ) // required
    ( optional
      ( range,          *, elements( imageSpace->mesh() )  )
      ( backend,        *, Backend<typename compute_opInterpolation_return<Args>::domain_space_type::value_type>::build() )
      ( type,           *, InterpolationNonConforme()  )
    ) // optionnal
)
{
    Feel::detail::ignore_unused_variable_warning( args );

    return opInterp( domainSpace,imageSpace,range,backend,type );

} // opInterpolation




} // Feel
#endif /* __OperatorInterpolation_H */
