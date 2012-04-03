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


namespace detailsup {

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
        _M_range(oi._M_range)
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


    range_iterator _M_range;

};

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::OperatorInterpolation( domain_space_ptrtype const& domainspace,
                                                                                                        dual_image_space_ptrtype const& imagespace,
                                                                                                        backend_ptrtype const& backend )
    :
    super( domainspace, imagespace, backend, false ),
    _M_range( elements(imagespace->mesh() ) )
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
    if ( this->dualImageSpace()->mesh().get() == (image_mesh_type*)this->domainSpace()->mesh().get() )
        {
            this->updateSameMesh();
        }
    else // no relation between meshes
        {
            this->updateNoRelationMesh();
        }

    // close matrix after build
    this->mat().close();
}

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::updateSameMesh()
{
    const size_type proc_id           = this->dualImageSpace()->mesh()->comm().rank();
    const size_type n1_dof_on_proc    = this->dualImageSpace()->nLocalDof();
    const size_type firstrow_dof_on_proc = this->dualImageSpace()->dof()->firstDof( proc_id );
    const size_type lastrow_dof_on_proc = this->dualImageSpace()->dof()->lastDof( proc_id );
    const size_type firstcol_dof_on_proc = this->domainSpace()->dof()->firstDof( proc_id );
    const size_type lastcol_dof_on_proc = this->domainSpace()->dof()->lastDof( proc_id );
    graph_ptrtype sparsity_graph( new graph_type( n1_dof_on_proc,
                                                  firstrow_dof_on_proc, lastrow_dof_on_proc,
                                                  firstrow_dof_on_proc, lastrow_dof_on_proc ) );

    auto const* imagedof = this->dualImageSpace()->dof().get();
    auto const* domaindof = this->domainSpace()->dof().get();
    auto const* imagebasis = this->dualImageSpace()->basis().get();
    auto const* domainbasis = this->domainSpace()->basis().get();


    std::vector<bool> dof_done( this->dualImageSpace()->nDof(), false);
    std::vector< std::list<std::pair<size_type,double> > > memory_valueInMatrix( this->dualImageSpace()->nDof() );

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
    for( ; it != en; ++ it )
        {
            auto idElem = detailsup::idElt(*it,idim_type());

            // Global assembly
            for ( uint16_type iloc = 0; iloc < nLocalDofInDualImageElt; ++iloc )
                {
                    for ( uint16_type comp = 0;comp < image_basis_type::nComponents;++comp )
                        {
                            size_type i =  boost::get<0>(imagedof->localToGlobal( *it, iloc, comp ));
                            if (!dof_done[i])
                                {
                                    const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[i];
                                    const auto theproc = imagedof->procOnGlobalCluster(ig1);
                                    auto& row = sparsity_graph->row(ig1);
                                    row.get<0>() = theproc;
                                    row.get<1>() = i;

                                    uint16_type ilocprime=imagedof->localDofInElement(*it, iloc, comp) ;

                                    for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                        {
                                            // get column
                                            const size_type j =  boost::get<0>(domaindof->localToGlobal( idElem, jloc, comp ));
                                            //up the pattern graph
                                            row.get<2>().insert(domaindof->mapGlobalProcessToGlobalCluster()[j]);
                                            // get interpolated value
                                            const value_type v = Mloc( domain_basis_type::nComponents1*jloc +
                                                                       comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof +
                                                                       comp,
                                                                       ilocprime );
                                            // save in matrux
                                            memory_valueInMatrix[i].push_back(std::make_pair(j,v));
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

template<typename DomainSpaceType, typename ImageSpaceType,typename IteratorRange,typename InterpType>
void
OperatorInterpolation<DomainSpaceType, ImageSpaceType,IteratorRange,InterpType>::updateNoRelationMesh()
{

    Debug( 5034 ) << "[interpolate] different meshes\n";

    const size_type proc_id           = this->dualImageSpace()->mesh()->comm().rank();
    const size_type n1_dof_on_proc    = this->dualImageSpace()->nLocalDof();
    const size_type firstrow_dof_on_proc = this->dualImageSpace()->dof()->firstDof( proc_id );
    const size_type lastrow_dof_on_proc = this->dualImageSpace()->dof()->lastDof( proc_id );
    const size_type firstcol_dof_on_proc = this->domainSpace()->dof()->firstDof( proc_id );
    const size_type lastcol_dof_on_proc = this->domainSpace()->dof()->lastDof( proc_id );

    graph_ptrtype sparsity_graph( new graph_type( n1_dof_on_proc,
                                                  firstrow_dof_on_proc, lastrow_dof_on_proc,
                                                  firstrow_dof_on_proc, lastrow_dof_on_proc ) );

    auto const* imagedof = this->dualImageSpace()->dof().get();
    auto const* domaindof = this->domainSpace()->dof().get();
    auto const* domainbasis = this->domainSpace()->basis().get();

    //-----------------------------------------
    //init the localization tool
    auto __loc = this->domainSpace()->mesh()->tool_localization();
    __loc->updateForUse();
    //__loc->kdtree()->nbNearNeighbor(3);
    //__loc->kdtree()->nbNearNeighbor(this->domainSpace()->mesh()->numElements());
    //__loc->setExtrapolation(false);

    //-----------------------------------------
    // usefull data
    typename matrix_node<value_type>::type ptsReal( image_mesh_type::nRealDim, 1);
    typename matrix_node<value_type>::type ptsRef(image_mesh_type::nRealDim , 1 );
    typename domain_mesh_type::Localization::container_search_iterator_type itanal,itanal_end;
    typename domain_mesh_type::Localization::container_output_iterator_type itL,itL_end;
    typename matrix_node<value_type>::type MlocEval(domain_basis_type::nLocalDof*domain_basis_type::nComponents1,1);

    std::vector<bool> dof_done( this->dualImageSpace()->nDof(), false);
    std::vector< std::list<std::pair<size_type,double> > > memory_valueInMatrix( this->dualImageSpace()->nDof() );

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
                            size_type gdof =  boost::get<0>(imagedof->localToGlobal( *it, iloc, comp ));
                            if (!dof_done[gdof])
                                {
                                    //------------------------
                                    // get the graph row
                                    const auto ig1 = imagedof->mapGlobalProcessToGlobalCluster()[gdof];
                                    const auto theproc = imagedof->procOnGlobalCluster(ig1);
                                    auto& row = sparsity_graph->row(ig1);
                                    row.get<0>() = theproc;
                                    row.get<1>() = gdof;
                                    //------------------------
                                    // the dof point
                                    ublas::column(ptsReal,0 ) = boost::get<0>(imagedof->dofPoint(gdof));
                                    //------------------------
                                    // localisation process
                                    __loc->run_analysis(ptsReal,it->vertices()/*it->G()*/,mpl::bool_<interpolation_type::value>());
                                    //------------------------
                                    // for each localised points
                                    itanal = __loc->result_analysis_begin();
                                    itanal_end = __loc->result_analysis_end();
                                    for ( ;itanal!=itanal_end;++itanal)
                                        {
                                            itL=itanal->second.begin();

                                            ublas::column( ptsRef, 0 ) = boost::get<1>(*itL);

                                            MlocEval = domainbasis->evaluate( ptsRef );

                                            for ( uint16_type jloc = 0; jloc < domain_basis_type::nLocalDof; ++jloc )
                                                {
                                                    //get global dof
                                                    size_type j =  boost::get<0>(domaindof->localToGlobal( itanal->first,jloc,comp ));
                                                    value_type v = MlocEval( domain_basis_type::nComponents1*jloc
                                                                             + comp*domain_basis_type::nComponents1*domain_basis_type::nLocalDof
                                                                             + comp,
                                                                             0 );
                                                    row.get<2>().insert(domaindof->mapGlobalProcessToGlobalCluster()[j]);
                                                    memory_valueInMatrix[gdof].push_back(std::make_pair(j,v));
                                                }
                                        }
                                    dof_done[gdof]=true;
                                }
                        }
                }
        } // for( ; it != en; ++ it )

    //-----------------------------------------
   // compute graph
   sparsity_graph->close();
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

    boost::shared_ptr<operatorinterpolation_type> opI( new operatorinterpolation_type(domainspace,imagespace,r,backend));

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
                         (typename compute_opInterpolation_return<Args>::type), // 1. return type
                         opInterpolation,                        // 2. name of the function template
                         tag,                                        // 3. namespace of tag types
                         (required
                          (domainSpace,    *(boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> >))
                          (imageSpace,     *(boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> >))
                          ) // required
                         (optional
                          (range,          *, elements(imageSpace->mesh())  )
                          (backend,        *, Backend<typename compute_opInterpolation_return<Args>::domain_space_type::value_type>::build())
                          (type,           *, InterpolationNonConforme()  )
                          ) // optionnal
                         )
{
    Feel::detail::ignore_unused_variable_warning(args);

    return opInterp(domainSpace,imageSpace,range,backend,type);

} // opInterpolation




} // Feel
#endif /* __OperatorInterpolation_H */
