/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-04-07

  Copyright (C) 2014 Feel++ Consortium

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
   \file element_impl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-04-07
 */
#ifndef FEELPP_ELEMENT_IMPL_HPP
#define FEELPP_ELEMENT_IMPL_HPP 1
namespace Feel{

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<int i>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::template sub_element<i>::type
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::elementImpl( std::string const& name,
                                                                      bool updateOffViews )
{
    size_type nbdof_start =  fusion::accumulate( this->functionSpaces(),
                                                 size_type( 0 ),
                                                 Feel::detail::NLocalDof<mpl::bool_<true> >( this->worldsComm(), false, 0, i ) );

    typename mpl::at_c<functionspace_vector_type,i>::type space( M_functionspace->template functionSpace<i>() );
    DVLOG(2) << "Element <" << i << ">::start :  "<< nbdof_start << "\n";
    DVLOG(2) << "Element <" << i << ">::size :  "<<  space->nDof()<< "\n";
    DVLOG(2) << "Element <" << i << ">::local size :  "<<  space->nLocalDof()<< "\n";
    DVLOG(2) << "Element <" << -1 << ">::size :  "<<  this->size() << "\n";

    if ( this->functionSpace()->template functionSpace<i>()->worldComm().isActive() )
    {
        ct_type ct( *this, ublas::range( nbdof_start, nbdof_start+space->nLocalDof() ),
                    M_functionspace->template functionSpace<i>()->dof() );

        // update M_containersOffProcess<i> : send
        if ( this->worldComm().globalSize()>1 && updateOffViews && !this->functionSpace()->hasEntriesForAllSpaces() )
        {
            std::vector<double> dataToSend( ct.begin(), ct.end() );

            if ( !M_containersOffProcess ) M_containersOffProcess = boost::in_place();

            fusion::for_each( *M_containersOffProcess, Feel::detail::SendContainersOn<i,functionspace_type>( this->functionSpace(), dataToSend ) );
        }

        DVLOG(2) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
        DVLOG(2) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";
        return typename sub_element<i>::type( space, ct, name );
    }

    else
    {
        // initialize if not the case
        if ( !M_containersOffProcess ) M_containersOffProcess = boost::in_place();

        fusion::for_each( *M_containersOffProcess, Feel::detail::InitializeContainersOff<i,functionspace_type>( this->functionSpace() ) );

        // update M_containersOffProcess<i> : recv
        if ( this->worldComm().globalSize()>1 && updateOffViews && !this->functionSpace()->hasEntriesForAllSpaces() )
        {
            fusion::for_each( *M_containersOffProcess, Feel::detail::RecvContainersOff<i,functionspace_type>( this->functionSpace() ) );
        }

        // build a subrange view identical
        ct_type ct( *fusion::at_c<i>( *M_containersOffProcess ),
                    ublas::range( 0, space->nLocalDof() ),
                    M_functionspace->template functionSpace<i>()->dof() );

        DVLOG(2) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
        DVLOG(2) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";

        return typename sub_element<i>::type( space, ct, name );
    }
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<int i,typename ExprT>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::template sub_element<i>::type &
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::element( ExprT e, std::string const& name,
                                                             bool updateOffViews,
                                                             typename std::enable_if<std::is_base_of<ExprBase,ExprT>::value >::type*  )
{
#if 0
    auto u = this->element<i>(name,updateOffViews);
    u.on( _range=elements(this->mesh()), _expr=e );
    return u;
#else
    this->element<i>(name,updateOffViews).on( _range=elements(this->mesh()), _expr=e );
    return this->element<i>(name,updateOffViews);
#endif
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<int i>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::template sub_element<i>::type
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::elementImpl( std::string const& name, bool updateOffViews ) const
{
    size_type nbdof_start =  fusion::accumulate( M_functionspace->functionSpaces(),
                                                 size_type( 0 ),
                                                 Feel::detail::NLocalDof<mpl::bool_<true> >( this->worldsComm(), false, 0, i ) );
    typename mpl::at_c<functionspace_vector_type,i>::type space( M_functionspace->template functionSpace<i>() );

    DVLOG(2) << "Element <" << i << ">::start :  "<< nbdof_start << "\n";
    DVLOG(2) << "Element <" << i << ">::size :  "<<  space->nDof()<< "\n";
    DVLOG(2) << "Element <" << i << ">::local size :  "<<  space->nLocalDof()<< "\n";
    DVLOG(2) << "Element <" << -1 << ">::size :  "<<  this->size() << "\n";

    if ( this->functionSpace()->worldsComm()[i].isActive() )
    {
        ct_type ct( const_cast<VectorUblas<value_type>&>( dynamic_cast<VectorUblas<value_type> const&>( *this ) ),
                    ublas::range( nbdof_start, nbdof_start+space->nLocalDof() ),
                    M_functionspace->template functionSpace<i>()->dof() );

        // update M_containersOffProcess<i> : send
        if ( this->worldComm().globalSize()>1 && updateOffViews && !this->functionSpace()->hasEntriesForAllSpaces() )
        {
            std::vector<double> dataToSend( ct.begin(), ct.end() );

            if ( !M_containersOffProcess ) M_containersOffProcess = boost::in_place();

            fusion::for_each( *M_containersOffProcess, Feel::detail::SendContainersOn<i,functionspace_type>( this->functionSpace(), dataToSend ) );
        }

        DVLOG(2) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
        DVLOG(2) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";
        return typename sub_element<i>::type( space, ct, name );

    }

    else
    {
        // update M_containersOffProcess<i> : recv
        if ( this->worldComm().globalSize()>1 && updateOffViews && !this->functionSpace()->hasEntriesForAllSpaces() )
        {
            fusion::for_each( *M_containersOffProcess, Feel::detail::RecvContainersOff<i,functionspace_type>( this->functionSpace() ) );
        }

        // build a subrange view identical
        ct_type ct( *fusion::at_c<i>( *M_containersOffProcess ),
                    ublas::range( 0, space->nLocalDof() ),
                    M_functionspace->template functionSpace<i>()->dof() );

        DVLOG(2) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
        DVLOG(2) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";
        return typename sub_element<i>::type( space, ct, name );
    }


}


//
// Element implementation
//
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element()
    :
    super(),
    M_start( 0 ),
    M_ct( NO_COMPONENT ),
    M_containersOffProcess( boost::none )
{}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( Element const& __e )
    :
    super( __e ),
    M_functionspace( __e.M_functionspace ),
    M_name( __e.M_name ),
    M_start( __e.M_start ),
    M_ct( __e.M_ct ),
    M_containersOffProcess( __e.M_containersOffProcess )
{
    DVLOG(2) << "Element<copy>::range::start = " << this->start() << "\n";
    DVLOG(2) << "Element<copy>::range::size = " << this->size() << "\n";
    this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
        std::string const& __name,
        size_type __start,
        ComponentType __ct )
    :
    super( __functionspace->dof() ),
    M_functionspace( __functionspace ),
    M_name( __name ),
    M_start( __start ),
    M_ct( __ct ),
    M_containersOffProcess( boost::none )
{
    DVLOG(2) << "Element::start = " << this->start() << "\n";
    DVLOG(2) << "Element::size = " << this->size() << "\n";
    DVLOG(2) << "Element::ndof = " << this->nDof() << "\n";
    DVLOG(2) << "Element::nlocaldof = " << this->nLocalDof() << "\n";
    this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
        container_type const& __c,
        std::string const& __name,
        size_type __start,
        ComponentType __ct )
    :
    super( __c ),
    M_functionspace( __functionspace ),
    M_name( __name ),
    M_start( __start ),
    M_ct( __ct ),
    M_containersOffProcess( boost::none )
{
    DVLOG(2) << "Element<range>::range::start = " << __c.start() << "\n";
    DVLOG(2) << "Element<range>::range::size = " << __c.size() << "\n";
    DVLOG(2) << "Element<range>::start = " << this->start() << "\n";
    DVLOG(2) << "Element<range>::size = " << this->size() << "\n";
    DVLOG(2) << "Element<range>::ndof = " << this->nDof() << "\n";
    DVLOG(2) << "Element<range>::nlocaldof = " << this->nLocalDof() << "\n";
    M_start = __c.start();
    this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::~Element()
{
    VLOG(1) << "Element destructor...";
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::initFromSpace( functionspace_ptrtype const& __functionspace,
        container_type const& __c )
{
    M_functionspace = __functionspace;
    ( container_type )*this = __c;
}

namespace detail
{
template<typename ElementType>
struct InitializeElement
{
    InitializeElement( ElementType * element )
        :
        M_element( element )
        {}
    template <typename T>
    void operator()( T & x ) const
    {
        typedef typename T::first_type key_type;
        typedef typename T::second_type::element_type myelt_type;
        std::string name = (boost::format("%1%_%2%")%M_element->name() %key_type::value).str();
        if( !x.second )
            x = std::make_pair(key_type(), boost::shared_ptr<myelt_type>( new myelt_type( M_element->template elementImpl<key_type::value>( name ) ) ) );
    }
    ElementType * M_element;
};
} // namespace detail
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::initSubElementView( mpl::true_ )
{
    fusion::for_each( M_elements,
                      Feel::detail::InitializeElement<Element<Y,Cont>>(this) );
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator=( Element<Y,Cont> const& __e )
{
    if (  this != &__e )
    {
        M_functionspace = __e.M_functionspace;

        if ( __e.M_name != "unknown" )
            M_name = __e.M_name;

        M_start = __e.M_start;
        M_ct = __e.M_ct;
        M_containersOffProcess = __e.M_containersOffProcess;

        this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );

        this->resize( M_functionspace->nLocalDof() );
        super::operator=( __e );
        this->outdateGlobalValues();
    }

    return *this;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename VectorExpr>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator=( VectorExpr const& __v )
{
    if (  __v.size() != this->size() )
    {
        std::ostringstream __err;
        __err << "Invalid vector size this->size()=" << this->size()
              << " and v.size()=" << __v.size();
        throw std::logic_error( __err.str() );
    }

    this->outdateGlobalValues();
    super::operator=( __v );
    return *this;
}


//
// Interpolation tools
//
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::id_type
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator()( node_type const& __x, bool extrapolate ) const
{
    this->updateGlobalValues();

    node_type __x_ref;
    size_type __cv_id;
    rank_type rank = functionSpace()->mesh()->comm().rank();
    rank_type nprocs = functionSpace()->mesh()->comm().size();
    std::vector<uint8_type/*int*/> found_pt( nprocs, 0 );
    std::vector<uint8_type/*int*/> global_found_pt( nprocs, 0 );

    if ( functionSpace()->findPoint( __x, __cv_id, __x_ref ) || extrapolate )
    {
        DVLOG(2) << "Point " << __x << " is in element " << __cv_id << " pt_ref=" << __x_ref << "\n";

        gm_ptrtype __gm = functionSpace()->gm();
        typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
        typedef typename gm_type::precompute_type geopc_type;
        typename matrix_node<value_type>::type pts( __x_ref.size(), 1 );
        ublas::column( pts, 0 ) = __x_ref;
        geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );


        typedef typename gm_type::template Context<vm::POINT|vm::GRAD|vm::KB|vm::JACOBIAN, geoelement_type> gmc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        gmc_ptrtype __c( new gmc_type( __gm,
                                       functionSpace()->mesh()->element( __cv_id ),
                                       __geopc ) );

        DCHECK( ublas::norm_2( __x-__c->xReal(0) ) < 1e-6 ) << "Point " << __x <<  " was not properly found, got " << __c->xReal(0);
        DVLOG(2) << "Point x=" << __x << " and c->xreal = " << __c->xReal(0);

        pc_ptrtype pc( new pc_type( this->functionSpace()->fe(), pts ) );

        typedef typename mesh_type::element_type geoelement_type;
        typedef typename functionspace_type::fe_type fe_type;
        typedef typename fe_type::template Context<vm::POINT|vm::GRAD|vm::KB|vm::JACOBIAN, fe_type, gm_type, geoelement_type> fectx_type;
        typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
        fectx_ptrtype fectx( new fectx_type( this->functionSpace()->fe(),
                                             __c,
                                             pc ) );
        found_pt[ rank ] = 1;

#if defined(FEELPP_HAS_MPI)

        if ( nprocs > 1 )
        {
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<uint8_type/*int*/>() );
        }

#else
        global_found_pt[ 0 ] = found_pt[ 0 ];
#endif /* FEELPP_HAS_MPI */


        rank_type procIdEval = rank;

        // in case where a point are localised on interprocess face, various process can find this point
        // we take only the process of smaller rank for evaluated the function at point
        if ( nprocs > 1 )
            for ( procIdEval=0 ; procIdEval < global_found_pt.size(); ++procIdEval )
                if ( global_found_pt[procIdEval] != 0 )
                {
                    DVLOG(2) << "processor " << procIdEval << " has the point " << __x << "and evaluated this one\n";
                    break;
                }

        id_type __id;
        if ( procIdEval == rank )
        {
            // this resize is crutial for serialisation!
            __id.resize( this->idExtents(*fectx ) );
            __id = this->id( *fectx );
            DVLOG(2) << "[interpolation]  id = " << __id << "\n";
        }

#if defined(FEELPP_HAS_MPI)
        DVLOG(2) << "sending interpolation context to all processors from " << functionSpace()->mesh()->comm().rank() << "\n";

        if ( nprocs > 1 )
        {
            mpi::broadcast( functionSpace()->mesh()->comm(), __id, procIdEval );
        }

        DVLOG(2) << "[interpolation] after broadcast id = " << __id << "\n";
#endif /* FEELPP_HAS_MPI */
        return __id;
    }

    else
    {
#if defined(FEELPP_HAS_MPI)
        if ( nprocs > 1 )
        {
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<uint8_type/*int*/>() );
        }

#endif /* FEELPP_HAS_MPI */
        bool found = false;
        rank_type i = 0;

        if ( nprocs > 1 )
            for ( ; i < global_found_pt.size(); ++i )
                if ( global_found_pt[i] != 0 )
                {
                    DVLOG(2) << "processor " << i << " has the point " << __x << "\n";
                    found = true;
                    break;
                }

        id_type __id;
        if ( found )
        {
            DVLOG(2) << "receiving interpolation context from processor " << i << "\n";
#if defined(FEELPP_HAS_MPI)

            if ( functionSpace()->mesh()->comm().size() > 1 )
            {
                mpi::broadcast( functionSpace()->mesh()->comm(), __id, i );
            }

#endif /* FEELPP_HAS_MPI */
        }

        else
        {
            LOG(WARNING) << "no processor seems to have the point " << __x << "\n";
        }

        return __id;
    }

} // operator()
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename Context_t>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::id_( Context_t const & context, id_array_type& v ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_size_type_value )
        return;

    const uint16_type nq = context.xRefs().size2();

    //double vsum=0;
    auto const& s = M_functionspace->dof()->localToGlobalSigns( elt_id );
    //array_type v( boost::extents[nComponents1][nComponents2][context.xRefs().size2()] );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( typename array_type::index c1 = 0; c1 < ncdof; ++c1 )
        {
            typename array_type::index ldof = basis_type::nDof*c1+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            //std::cout << "ldof = " << ldof << "\n";
            //std::cout << "gdof = " << gdof << "\n";

            DLOG_IF(WARNING, !( gdof < this->size() ) ) << "FunctionSpace::Element invalid access index "
                                                        << context.eId()
                                                        << "    l=" << l
                                                        << "   c1=" << c1
                                                        << " ldof=" << ldof
                                                        << " gdof=" << gdof
                                                        << " size=" << this->size()
                                                        << " localSize=" << this->localSize();

            value_type v_ = this->globalValue( gdof );

            //std::cout << "v_ =" << v_ << "\n";
            //for( typename array_type::index c2 = 0; c2 < nComponents2; ++c2 )
            for ( uint16_type q = 0; q < nq; ++q )
            {
                for ( typename array_type::index i = 0; i < nComponents1; ++i )
                    //for( typename array_type::index j = 0; j < nComponents2; ++j )
                {
                    v[q]( i,0 ) += s(ldof)*v_*context.id( ldof, i, 0, q );
                    //vsum +=v_*context.id( ldof, i, 0, q );
                    //v[q](i,0) += v_*context.gmc()->J(*)*context.pc()->phi( ldof, i, 0, q );
                }
            }
        }
    }

    //LOG( INFO ) << "vsum : "<<vsum;
    //return v;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::idInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );

    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );

    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );

        //update precompute of basis functions
        __pc->update( pts );

        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        id_type __id( this->id( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
        {
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
            {
                v[boost::get<0>( *itL )]( i,0 ) =  __id( i,0,k );
            }
        }
    }

}

//
// Grad
//
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::grad_type
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::grad( node_type const& __x ) const
{
    this->updateGlobalValues();

    node_type __x_ref;
    size_type __cv_id;
    std::vector<int> found_pt( functionSpace()->mesh()->comm().size(), 0 );
    std::vector<int> global_found_pt( functionSpace()->mesh()->comm().size(), 0 );

    if ( functionSpace()->findPoint( __x, __cv_id, __x_ref ) )
    {
#if !defined( NDEBUG )
        DVLOG(2) << "Point " << __x << " is in element " << __cv_id << " pt_ref=" << __x_ref << "\n";
#endif
        gm_ptrtype __gm = functionSpace()->gm();
        typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
        typedef typename gm_type::precompute_type geopc_type;
        typename matrix_node<value_type>::type pts( __x_ref.size(), 1 );
        ublas::column( pts, 0 ) = __x_ref;
        geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );


        typedef typename gm_type::template Context<vm::POINT|vm::JACOBIAN|vm::KB, geoelement_type> gmc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        gmc_ptrtype __c( new gmc_type( __gm,
                                       functionSpace()->mesh()->element( __cv_id ),
                                       __geopc ) );
        pc_ptrtype pc( new pc_type( this->functionSpace()->fe(), pts ) );
        typedef typename mesh_type::element_type geoelement_type;
        typedef typename functionspace_type::fe_type fe_type;
        typedef typename fe_type::template Context<vm::GRAD, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
        typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
        fectx_ptrtype fectx( new fectx_type( this->functionSpace()->fe(),
                                             __c,
                                             pc ) );

        found_pt[ functionSpace()->mesh()->comm().rank() ] = 1;

#if defined(FEELPP_HAS_MPI)

        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<int>() );
        }

#else
        global_found_pt[ 0 ] = found_pt[ 0 ];
#endif /* FEELPP_HAS_MPI */

        grad_type g_( this->grad( *fectx ) );
        //DVLOG(2) << "[interpolation]  id = " << v << "\n";
#if defined(FEELPP_HAS_MPI)
        DVLOG(2) << "sending interpolation context to all processors from " << functionSpace()->mesh()->comm().rank() << "\n";

        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            mpi::broadcast( functionSpace()->mesh()->comm(), g_, functionSpace()->mesh()->comm().rank() );
        }

        //DVLOG(2) << "[interpolation] after broadcast g_ = " << g_ << "\n";
#endif /* FEELPP_HAS_MPI */
        return g_;
    }

    else
    {
#if defined(FEELPP_HAS_MPI)

        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<int>() );
        }

#endif /* FEELPP_HAS_MPI */
        bool found = false;
        size_type i = 0;

        for ( ; i < global_found_pt.size(); ++i )
            if ( global_found_pt[i] != 0 )
            {
                DVLOG(2) << "processor " << i << " has the point " << __x << "\n";
                found = true;
                break;
            }

        grad_type g_;

        if ( found )
        {
            DVLOG(2) << "receiving interpolation context from processor " << i << "\n";
#if defined(FEELPP_HAS_MPI)

            if ( functionSpace()->mesh()->comm().size() > 1 )
                mpi::broadcast( functionSpace()->mesh()->comm(), g_, i );

#endif /* FEELPP_HAS_MPI */

            //DVLOG(2) << "[interpolation] after broadcast id = " << v << "\n";
        }

        else
        {
            Warning() << "no processor seems to have the point " << __x << "\n";
        }

        return g_;
    }

} // grad
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::grad_( ContextType const & context, grad_array_type& v ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_size_type_value )
        return;

    //std::cout << "coeff=" << coeff << "\n";
    //array_type v( boost::extents[nComponents1][nRealDim][context.xRefs().size2()] );
    //std::fill( v.data(), v.data()+v.num_elements(), value_type( 0 ) );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            //value_type v_ = (*this)( gdof );
            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            {
                for ( int k = 0; k < nComponents1; ++k )
                    for ( int j = 0; j < nRealDim; ++j )
                    {
                        v[q]( k,j ) += v_*context.grad( ldof, k, j, q );
                    }
            }
        }
    }

#if 0
    std::cout << "xref = " << context.xRefs() << "\n";
    std::cout << "xreal = " << context.xReal() << "\n";
    std::cout << "pts = " << context.G() << "\n";

    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        for ( typename array_type::index i = 0; i < nDim; ++i )
        {
            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
                std::cout << "Gref[ " << l << "," << i << "," << q << "] = " << pc.grad( l, 0, i, q ) << "\n";

        }
    }

    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        for ( int j = 0; j < nRealDim; ++j )
            for ( typename array_type::index i = 0; i < nDim; ++i )
            {
                for ( size_type q = 0; q < context.xRefs().size2(); ++q )
                    std::cout << "G[ " << l << "," << j << "," << i << "," << q << "] = " << context.B( q )( j, i )*pc.grad( l, 0, i, q ) << "\n";

            }
    }

    for ( int k = 0; k < nComponents1; ++k )
        for ( int j = 0; j < nRealDim; ++j )
        {
            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            {
                std::cout << "v[ " << k << "," << j << "," << q << "]= " << v[q]( k,j )  << "\n";
                std::cout << "B(" << q << ")=" << context.B( q ) << "\n";
                std::cout << "J(" << q << ")=" << context.J( q ) << "\n";
                std::cout << "K(" << q << ")=" << context.K( q ) << "\n";
            }
        }

#endif
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::gradInterpolate(  matrix_node_type __ptsReal, grad_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );

    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::KB|vm::GRAD|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );


        //update precompute of basis functions
        __pc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        grad_type __grad( this->grad( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
        {
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
                for ( uint j = 0; j < nRealDim; ++j )
                {
                    v[boost::get<0>( *itL )]( i,j ) = __grad( i,j,k );
                }
        }

    }

}


template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::div_( ContextType const & context, div_array_type& v ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_size_type_value )
        return;

    const size_type Q = context.xRefs().size2();

    auto const& s = M_functionspace->dof()->localToGlobalSigns( elt_id );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );
            //value_type v_ = (*this)( gdof );
            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < Q; ++q )
            {
#if 0
                std::cout << "v(" << gdof << ")=" << v_ << "\n";
                std::cout << "context.div(" << ldof << "," << q << ")="
                          << context.div( ldof, 0, 0, q ) << "\n" ;
#endif
                v[q]( 0,0 ) += s(ldof)*v_*context.div( ldof, 0, 0, q );
            }
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::divInterpolate( matrix_node_type __ptsReal, div_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::DIV|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        div_type __div( this->div( *__ctx ) );

        //update the output data
        //for( typename array_type::index i = 0; i < nComponents1; ++i )
        //   {
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            v[boost::get<0>( *itL )]( 0,0 ) =  __div( 0,0,k );

        //   }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curl_( ContextType const & context, curl_array_type& v ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_size_type_value )
        return;

    auto const& s = M_functionspace->dof()->localToGlobalSigns( elt_id );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            //value_type v_ = (*this)( gdof );
            value_type v_ = this->globalValue( gdof );
            const uint16_type nq = context.xRefs().size2();

            for ( uint16_type q = 0; q < nq ; ++q )
            {
                if ( nDim == 3 )
                {
                    for ( typename array_type::index i = 0; i < nDim; ++i )
                    {
                        v[q]( i,0 ) += s(ldof)*v_*context.curl( ldof, i, 0, q );
                    }
                }

                else if ( nDim == 2 )
                {
                    v[q]( 0,0 ) += s(ldof)*v_*context.curl( ldof, 0, 0, q );
                }

            }
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curl_( ContextType const & context, comp_curl_array_type& v, int comp ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_size_type_value )
        return;

    auto const& s = M_functionspace->dof()->localToGlobalSigns( elt_id );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            //value_type v_ = (*this)( gdof );
            value_type v_ = this->globalValue( gdof );
            const uint16_type nq = context.xRefs().size2();

            for ( uint16_type q = 0; q < nq ; ++q )
            {
                if ( nDim == 3 )
                {
                    v[q]( 0,0 ) += s(ldof)*v_*context.curl( ldof, comp, 0, q );
                }

                else if ( nDim == 2 )
                {
                    v[q]( 0,0 ) += s(ldof)*v_*context.curl( ldof, 2, 0, q );
                }

            }
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlInterpolate( matrix_node_type __ptsReal, curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );


    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curl_type __curl( this->curl( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curl( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curl( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curl( 2,0,k );
            }
        }

        else if ( nDim == 2 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
                v[boost::get<0>( *itL )]( 0,0 ) =  __curl( 0,0,k );
        }
    }

}


template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlxInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );


    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curlx_type __curlx( this->curlx( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curlx( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curlx( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curlx( 2,0,k );
            }
        }

        else if ( nDim == 2 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
                v[boost::get<0>( *itL )]( 0,0 ) =  __curlx( 0,0,k );
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlyInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curly_type __curly( this->curly( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curly( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curly( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curly( 2,0,k );
            }
        }

        else if ( nDim == 2 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
                v[boost::get<0>( *itL )]( 0,0 ) =  __curly( 0,0,k );
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlzInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curlz_type __curlz( this->curlz( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curlz( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curlz( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curlz( 2,0,k );
            }
        }

        else if ( nDim == 2 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
                v[boost::get<0>( *itL )]( 0,0 ) =  __curlz( 0,0,k );
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::d_( int N, ContextType const & context, id_array_type& v ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    for ( int i = 0; i < basis_type::nDof; ++i )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            size_type ldof = basis_type::nDof*c1 + i;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( context.eId(), i, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( i )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            {
                for ( typename array_type::index i = 0; i < nComponents1; ++i )
                {
                    v[q]( i,0 ) += v_*context.d( ldof, i, N, q );
                }
            }

        }
    }
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dxInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::KB|vm::GRAD|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        dx_type __dx( this->dx( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
            {
                v[boost::get<0>( *itL )]( i,0 ) =  __dx( i,0,k );
            }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dyInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::KB|vm::GRAD|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        dy_type __dy( this->dy( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
            {
                v[boost::get<0>( *itL )]( i,0 ) =  __dy( i,0,k );
            }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dzInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::KB|vm::GRAD|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );
        //evaluate element for these points
        dz_type __dz( this->dz( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
            {
                v[boost::get<0>( *itL )]( i,0 ) =  __dz( i,0,k );
            }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::hess_( ContextType const & context, hess_array_type& v ) const
{
    hess_( context, v, mpl::int_<rank>() );
}
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::hess_( ContextType const & context, hess_array_type& v, mpl::int_<0> ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    //FEELPP_ASSERT( comp < nRealDim )( comp )( nRealDim ).error( "[FunctionSpace::Element] grad: invalid component" );
    for ( int i = 0; i < basis_type::nDof; ++i )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            size_type ldof = basis_type::nDof*c1 + i;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( context.eId(), i, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( i )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            {
                for ( int i = 0; i < nRealDim; ++i )
                    for ( int j = 0; j < nRealDim; ++j )
                    {
                        v[q]( i,j ) += v_*context.hess( ldof, i, j, q );
                    } // i,j
            } // q

        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::hessInterpolate( matrix_node_type __ptsReal, hess_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::KB|vm::HESSIAN|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        hess_type __hess( this->hess( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            for ( int i = 0; i < nRealDim; ++i )
                for ( int j = 0; j < nRealDim; ++j )
                    v[boost::get<0>( *itL )]( i,j ) =  __hess( i,j,k );
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename IteratorType,  typename ExprType>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::onImpl( std::pair<IteratorType,IteratorType> const& r,
                                                            ExprType const& ex,
                                                            std::string const& prefix,
                                                            GeomapStrategyType geomap_strategy,
                                                            bool accumulate,
                                                            bool verbose,
                                                            mpl::int_<MESH_ELEMENTS> )
{
    const size_type context = ExprType::context|vm::POINT|vm::KB|vm::JACOBIAN;
    typedef ExprType expression_type;
    typedef Element<Y,Cont> element_type;
    // mesh element
    typedef typename functionspace_type::mesh_type::element_type geoelement_type;

    // geometric mapping context
    typedef typename functionspace_type::mesh_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef typename gm_type::template Context<context, geoelement_type> gm_context_type;
    typedef typename geoelement_type::gm1_type gm1_type;
    typedef typename gm1_type::template Context<context, geoelement_type> gm1_context_type;
    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;
    typedef boost::shared_ptr<gm1_context_type> gm1_context_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm1_context_ptrtype> > map_gmc1_type;


    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

    // basis
    typedef typename element_type::functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context< context, fe_type, gm_type, geoelement_type> fecontext_type;
    typedef boost::shared_ptr<fecontext_type> fecontext_ptrtype;
    //typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, fecontext_ptrtype> > map_gmc_type;

    // expression
    //typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;
    typedef typename expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename expression_type::template tensor<map_gmc1_type> t_expr1_type;
    typedef typename t_expr_type::shape shape;

    //
    // start
    //
    boost::timer __timer;

    std::vector<int> dofs;
    std::vector<value_type> values;
    auto it = r.first;
    auto en = r.second;

    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    auto* __fe = this->functionSpace()->fe().get();
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( this->functionSpace()->gm(),
                                                                                         __fe->points() ) );
    typename gm1_type::precompute_ptrtype __geopc1( new typename gm1_type::precompute_type( this->mesh()->gm1(),
                                                                                            __fe->points() ) );



    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    // return if no elements
    if ( it == en )
        return;

    gm_context_ptrtype __c( new gm_context_type( this->functionSpace()->gm(),*it,__geopc ) );
    gm1_context_ptrtype __c1( new gm1_context_type( this->mesh()->gm1(),*it,__geopc1 ) );

    typedef typename t_expr_type::shape shape;
    static const bool is_rank_ok = ( shape::M == nComponents1 &&
                                     shape::N == nComponents2 );

    BOOST_MPL_ASSERT_MSG( is_rank_ok,//mpl::bool_<is_rank_ok>::value,
                          INVALID_TENSOR_RANK,
                          ( mpl::int_<shape::M>, mpl::int_<nComponents>, shape ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    //t_expr_type tensor_expr( basis_type::isomorphism( ex ), mapgmc );
    t_expr_type tensor_expr( ex, mapgmc );

    map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

    t_expr1_type tensor_expr1( ex, mapgmc1 );

    std::vector<bool> points_done( this->functionSpace()->dof()->nLocalDof()/this->nComponents );
    std::fill( points_done.begin(), points_done.end(),false );

    auto IhLoc = __fe->localInterpolant();
    for ( ; it!=en ; ++it )
    {
        geoelement_type const& curElt = boost::unwrap_ref(*it);

        switch ( geomap_strategy )
        {
        case GeomapStrategyType::GEOMAP_HO:
        {
            __c->update( *it );
            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
            tensor_expr.update( mapgmc );
            __fe->interpolate( tensor_expr, IhLoc );
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->update( *it );
            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
            tensor_expr1.update( mapgmc1 );
            __fe->interpolate( tensor_expr1, IhLoc );
        }
        break;

        case GeomapStrategyType::GEOMAP_OPT:
        {
            if ( curElt.isOnBoundary() )
            {
                // HO if on boundary
                __c->update( *it );
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                tensor_expr.update( mapgmc );
                __fe->interpolate( tensor_expr, IhLoc );
            }

            else
            {
                __c1->update( *it );
                map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                tensor_expr1.update( mapgmc1 );
                __fe->interpolate( tensor_expr1, IhLoc );
            }
        }
        break;
        }
        if ( accumulate )
            this->plus_assign( curElt, IhLoc );
        else
            this->assign( curElt, IhLoc );


        //if P0 continuous finish loop here
        if ( isP0Continuous<fe_type>::result )
        {
            break;
        }
    }
} // onImpl (MESH_ELEMENTS)

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename IteratorType,  typename ExprType>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::onImpl( std::pair<IteratorType,IteratorType> const& r,
                                                            ExprType const& ex,
                                                            std::string const& prefix,
                                                            GeomapStrategyType geomap_strategy,
                                                            bool accumulate,
                                                            bool verbose,
                                                            mpl::int_<MESH_FACES> )
{
    typedef ExprType expression_type;
    typedef Element<Y,Cont> element_type;

    // mesh element
    typedef typename element_type::functionspace_type::mesh_type::element_type geoelement_type;
    typedef typename element_type::functionspace_type::mesh_type::face_type face_type;
    //typedef typename geoelement_type::face_type face_type;

    // geometric mapping context
    typedef typename geoelement_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename geoelement_type::gm1_type gm1_type;
    typedef boost::shared_ptr<gm1_type> gm1_ptrtype;

    typedef typename element_type::functionspace_type::fe_type fe_type;
    const size_type context = mpl::if_< mpl::or_<is_hdiv_conforming<fe_type>, is_hcurl_conforming<fe_type> >,
                                        mpl::int_<ExprType::context|vm::POINT|vm::JACOBIAN>,
                                        mpl::int_<ExprType::context|vm::POINT> >::type::value;


    typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef typename gm1_type::template Context<context, geoelement_type> gmc1_type;
    typedef boost::shared_ptr<gmc1_type> gmc1_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;


    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

    // basis
    typedef typename fe_type::template Context< context, fe_type, gm_type, geoelement_type> fecontext_type;
    typedef boost::shared_ptr<fecontext_type> fecontext_ptrtype;
    typedef typename fe_type::template Context< context, fe_type, gm1_type, geoelement_type> fecontext1_type;
    typedef boost::shared_ptr<fecontext1_type> fecontext1_ptrtype;
    //typedef fusion::map<fusion::pair<vf::detail::gmc<0>, fecontext_ptrtype> > map_gmc_type;

    // expression
    //typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;
    //typedef decltype( basis_type::isomorphism( M_expr ) ) the_expression_type;
    typedef expression_type the_expression_type;
    typedef typename boost::remove_reference<typename boost::remove_const<the_expression_type>::type >::type iso_expression_type;
    typedef typename iso_expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename iso_expression_type::template tensor<map_gmc1_type> t_expr1_type;
    typedef typename t_expr_type::shape shape;

    //
    // start
    //
    DVLOG(2)  << "assembling Dirichlet conditions\n";

    dof_type const* __dof = this->functionSpace()->dof().get();

    fe_type const* __fe = this->functionSpace()->fe().get();

    auto __face_it = r.first;
    auto const __face_en = r.second;

    if ( __face_it == __face_en )
        return;

    gm_ptrtype __gm( new gm_type );
    gm1_ptrtype __gm1( new gm1_type );



    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typedef typename geoelement_type::permutation_type permutation_type;
    typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
    typedef typename gm_type::precompute_type geopc_type;
    typedef typename gm1_type::precompute_ptrtype geopc1_ptrtype;
    typedef typename gm1_type::precompute_type geopc1_type;

    DVLOG(2)  << "[integratoron] numTopologicalFaces = " << geoelement_type::numTopologicalFaces << "\n";
    std::vector<std::map<permutation_type, geopc_ptrtype> > __geopc( geoelement_type::numTopologicalFaces );
    std::vector<std::map<permutation_type, geopc1_ptrtype> > __geopc1( geoelement_type::numTopologicalFaces );

    for ( uint16_type __f = 0; __f < geoelement_type::numTopologicalFaces; ++__f )
    {
        for ( permutation_type __p( permutation_type::IDENTITY );
                __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            __geopc[__f][__p] = geopc_ptrtype(  new geopc_type( __gm, __fe->points( __f ) ) );
            __geopc1[__f][__p] = geopc1_ptrtype(  new geopc1_type( __gm1, __fe->points( __f ) ) );
            //DVLOG(2) << "[geopc] FACE_ID = " << __f << " ref pts=" << __fe->dual().points( __f ) << "\n";
            FEELPP_ASSERT( __geopc[__f][__p]->nPoints() ).error( "invalid number of points for geopc" );
            FEELPP_ASSERT( __geopc1[__f][__p]->nPoints() ).error( "invalid number of points for geopc1" );
        }
    }
    face_type const& firstFace = *__face_it;
    uint16_type __face_id = firstFace.pos_first();
    gmc_ptrtype __c( new gmc_type( __gm, firstFace.element( 0 ), __geopc, __face_id ) );
    gmc1_ptrtype __c1( new gmc1_type( __gm1, firstFace.element( 0 ), __geopc1, __face_id ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    t_expr_type expr( ex, mapgmc );
    map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
    t_expr1_type expr1( ex, mapgmc1 );




    size_type nbFaceDof = invalid_size_type_value;

    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );

    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    DVLOG(2)  << "[projector::operator(MESH_FACES)] nbFaceDof = " << nbFaceDof << "\n";

    auto IhLoc = __fe->faceLocalInterpolant();
    for ( ; __face_it != __face_en; ++__face_it )
    {
        face_type const& curFace = boost::unwrap_ref(*__face_it);
        FEELPP_ASSERT( curFace.isOnBoundary() && !curFace.isConnectedTo1() )
        ( curFace.marker() )
        ( curFace.isOnBoundary() )
        ( curFace.ad_first() )
        ( curFace.pos_first() )
        ( curFace.ad_second() )
        ( curFace.pos_second() )
        ( curFace.id() ).warn( "inconsistent data face" );
        DVLOG(2) << "[projector] FACE_ID = " << curFace.id()
                      << " element id= " << curFace.ad_first()
                      << " pos in elt= " << curFace.pos_first()
                      << " marker: " << curFace.marker() << "\n";
        DVLOG(2) << "[projector] FACE_ID = " << curFace.id() << " real pts=" << curFace.G() << "\n";

        uint16_type __face_id = curFace.pos_first();

        std::pair<size_type,size_type> range_dof( std::make_pair( this->start(),
                this->functionSpace()->nDof() ) );
        DVLOG(2)  << "[projector] dof start = " << range_dof.first << "\n";
        DVLOG(2)  << "[projector] dof range = " << range_dof.second << "\n";

        switch ( geomap_strategy )
        {
        default:
        case GeomapStrategyType::GEOMAP_OPT:
        case GeomapStrategyType::GEOMAP_HO:
        {
            __c->update( curFace.element( 0 ), __face_id );
            DVLOG(2) << "[projector] FACE_ID = " << curFace.id() << "  ref pts=" << __c->xRefs() << "\n";
            DVLOG(2) << "[projector] FACE_ID = " << curFace.id() << " real pts=" << __c->xReal() << "\n";

            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

            expr.update( mapgmc );
            __fe->faceInterpolate( expr, IhLoc );
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->update( curFace.element( 0 ), __face_id );
            DVLOG(2) << "[projector] FACE_ID = " << curFace.id() << "  ref pts=" << __c1->xRefs() << "\n";
            DVLOG(2) << "[projector] FACE_ID = " << curFace.id() << " real pts=" << __c1->xReal() << "\n";

            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

            expr1.update( mapgmc1 );
            __fe->faceInterpolate( expr1, IhLoc );
        }
        break;
        }
        if ( accumulate )
            this->plus_assign( curFace, IhLoc );
        else
            this->assign( curFace, IhLoc );
    } // face_it

}


}

#endif
