/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-04-07

  Copyright (C) 2014-2016 Feel++ Consortium

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


#if BOOST_VERSION >= 105900
#include <boost/utility/in_place_factory.hpp>
#endif

#include <feel/feelvf/detail/gmc.hpp>

namespace Feel{

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<int i>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::template sub_element<i>::type
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::elementImpl( std::string const& name,
                                                                      bool updateOffViews )
{
    size_type nbdof_start = this->functionSpace()->nLocalDofWithoutGhostStart( i );
    size_type nbdofWithGhost_start = this->functionSpace()->nLocalDofWithGhostStart( i );
    size_type startDofIndexGhost = nbdofWithGhost_start - nbdof_start;
    if ( !Cont::is_shallow_array_adaptor_vector )
        startDofIndexGhost += this->functionSpace()->dof()->nLocalDofWithoutGhost();


    typename mpl::at_c<functionspace_vector_type,i>::type space( M_functionspace->template functionSpace<i>() );
    DVLOG(2) << "Element <" << i << ">::start :  "<< nbdof_start << "\n";
    DVLOG(2) << "Element <" << i << ">::size :  "<<  space->nDof()<< "\n";
    DVLOG(2) << "Element <" << i << ">::local size :  "<<  space->nLocalDof()<< "\n";
    DVLOG(2) << "Element <" << -1 << ">::size :  "<<  this->size() << "\n";

    if ( this->functionSpace()->template functionSpace<i>()->worldComm().isActive() )
    {
        ct_type ct( *this, ublas::range( nbdof_start, nbdof_start+space->dof()->nLocalDofWithoutGhost() ),
                    ublas::range( startDofIndexGhost, startDofIndexGhost+space->dof()->nLocalGhosts() ),
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
                    ublas::range( 0, 0 ),
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
#if 0
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<int i>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::template sub_element<i>::type
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::elementImpl( std::string const& name, bool updateOffViews ) const
{
    size_type nbdof_start = this->functionSpace()->nLocalDofWithoutGhostStart( i );
    size_type nbdofWithGhost_start = this->functionSpace()->nLocalDofWithGhostStart( i );
    size_type startDofIndexGhost = nbdofWithGhost_start - nbdof_start;
    if ( !Cont::is_shallow_array_adaptor_vector )
        startDofIndexGhost += this->functionSpace()->dof()->nLocalDofWithoutGhost();
    typename mpl::at_c<functionspace_vector_type,i>::type space( M_functionspace->template functionSpace<i>() );

    DVLOG(2) << "Element <" << i << ">::start :  "<< nbdof_start << "\n";
    DVLOG(2) << "Element <" << i << ">::size :  "<<  space->nDof()<< "\n";
    DVLOG(2) << "Element <" << i << ">::local size :  "<<  space->nLocalDof()<< "\n";
    DVLOG(2) << "Element <" << -1 << ">::size :  "<<  this->size() << "\n";

    if ( this->functionSpace()->worldsComm()[i].isActive() )
    {
        ct_type ct( const_cast<VectorUblas<value_type>&>( dynamic_cast<VectorUblas<value_type> const&>( *this ) ),
                    ublas::range( nbdof_start, nbdof_start+space->dof()->nLocalDofWithoutGhost() ),
                    ublas::range( startDofIndexGhost, startDofIndexGhost+space->dof()->nLocalGhosts() ),
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
                    ublas::range( 0, 0 ),
                    M_functionspace->template functionSpace<i>()->dof() );

        DVLOG(2) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
        DVLOG(2) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";
        return typename sub_element<i>::type( space, ct, name );
    }


}
#endif

//
// Element implementation
//
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element()
    :
    super(),
    M_start( 0 ),
    M_ct( ComponentType::NO_COMPONENT ),
    M_ct2( ComponentType::NO_COMPONENT ),
    M_containersOffProcess( boost::none )
{
    this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( Element const& __e )
    :
    super( __e ),
    M_functionspace( __e.M_functionspace ),
    M_name( __e.M_name ),
    M_start( __e.M_start ),
    M_ct( __e.M_ct ),
    M_ct2( __e.M_ct2 ),
    M_containersOffProcess( __e.M_containersOffProcess )
{
    DVLOG(2) << "Element<copy>::range::start = " << this->start() << "\n";
    DVLOG(2) << "Element<copy>::range::size = " << this->size() << "\n";
    this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
}

#if 0
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( Element && __e )
    :
    super( __e ),
    M_functionspace( std::move(__e.M_functionspace) ),
    M_name( std::move(__e.M_name) ),
    M_start(0 ),
    M_ct( ComponentType::NO_COMPONENT ),
    M_ct2( ComponentType::NO_COMPONENT ),
    M_containersOffProcess( boost::none )
{
    M_start = __e.M_start;
    M_ct = __e.M_ct;
    M_ct2 = __e.M_ct2;
    M_containersOffProcess = __e.M_containersOffProcess;

    super::operator=( __e );
    this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );

    __e.M_functionspace = nullptr;
    __e.M_start = 0;
    __e.M_ct = ComponentType::NO_COMPONENT;
    __e.M_ct2 = ComponentType::NO_COMPONENT;
    __e.M_containersOffProcess = boost::none;
}
#endif

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
                                                             std::string const& __name,
                                                             std::string const& __desc,
                                                             size_type __start,
                                                             ComponentType __ct )
    :
    super( __functionspace->dof() ),
    M_functionspace( __functionspace ),
    M_name( __name ),
    M_desc( __desc ),
    M_start( __start ),
    M_ct( __ct ),
    M_ct2( ComponentType::NO_COMPONENT ),
    M_containersOffProcess( boost::none )
{
    LOG(INFO) << "creating element " << name() << " : " << description();
    DVLOG(2) << "Element::start = " << this->start() << "\n";
    DVLOG(2) << "Element::size = " << this->size() << "\n";
    DVLOG(2) << "Element::ndof = " << this->nDof() << "\n";
    DVLOG(2) << "Element::nlocaldof = " << this->nLocalDof() << "\n";
    this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
                                                             std::string const& __name,
                                                             size_type __start,
                                                             ComponentType __ct )
    :
    Element( __functionspace, __name, __name, __start, __ct )
{
}


template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
                                                             container_type const& __c,
                                                             std::string const& __name,
                                                             std::string const& __desc,
                                                             size_type __start,
                                                             ComponentType __ct, ComponentType __ct2 )
    :
    super( __c ),
    M_functionspace( __functionspace ),
    M_name( __name ),
    M_desc( __desc ),
    M_start( __start ),
    M_ct( __ct ),
    M_ct2( __ct2 ),
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
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
                                                             container_type const& __c,
                                                             std::string const& __name,
                                                             size_type __start,
                                                             ComponentType __ct, ComponentType __ct2 )
    :
    Element( __functionspace, __c, __name, __name, __start, __ct, __ct2 )
{}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
                                                             size_type nActiveDof, value_type* arrayActiveDof,
                                                             size_type nGhostDof, value_type* arrayGhostDof )
    :
    super( nActiveDof,arrayActiveDof,nGhostDof,arrayGhostDof, __functionspace->dof() ),
    M_functionspace( __functionspace ),
    //M_name( __name ),
    //M_desc( __desc ),
    M_start( 0 ),
    M_ct( ComponentType::NO_COMPONENT ),
    M_ct2( ComponentType::NO_COMPONENT ),
    M_containersOffProcess( boost::none )
{
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

        if( M_element->functionSpace() )
        {
            // build view if not built or built with an empty space
            if ( !x.second || ( !x.second->functionSpace() ) )
            {
                auto e = M_element->template elementImpl<key_type::value>( name );
                auto sp = std::make_shared<myelt_type>( e );
                x = std::make_pair(key_type(), sp );
            }
        }
        else if ( !x.second )
        {
            x = std::make_pair(key_type(), nullptr );
        }
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
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator=( Element<Y,Cont> const& __e )
{
    if (  this != &__e )
    {
        M_functionspace = __e.M_functionspace;

        if ( __e.M_name != "unknown" )
            M_name = __e.M_name;

        M_start = __e.M_start;
        M_ct = __e.M_ct;
        M_ct2 = __e.M_ct2;
        M_containersOffProcess = __e.M_containersOffProcess;

        super::operator=( __e );
        this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
    }

    return *this;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator=( Element<Y,Cont> && __e )
{
    if (  this != &__e )
    {
        M_functionspace = std::move(__e.M_functionspace);
        M_name = std::move(__e.M_name);
        M_start = __e.M_start;
        M_ct = __e.M_ct;
        M_ct2 = __e.M_ct2;
        M_containersOffProcess = __e.M_containersOffProcess;

        super::operator=( __e );
        this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );

        __e.M_functionspace = nullptr;
        __e.M_start = 0;
        __e.M_ct = ComponentType::NO_COMPONENT;
        __e.M_ct2 = ComponentType::NO_COMPONENT;
        __e.M_containersOffProcess = boost::none;
    }

    return *this;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContOtherType>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator=( Element<Y,ContOtherType> const& v )
{
    if ( !M_functionspace )
    {
        M_functionspace = v.M_functionspace;
        if ( v.M_name != "unknown" )
            M_name = v.M_name;
    }
    super::operator=( v );
    this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
    return *this;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename VectorExpr>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator=( VectorExpr const& __v )
{
    super::operator=( __v );
    this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
    if ( !M_functionspace )
        LOG(WARNING) << "no function space, only algebraic view is operational";

    return *this;
}


//
// Interpolation tools
//
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::id_type
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator()( node_type const& __x, bool extrapolate, bool parallel ) const
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

        const size_type gmc_v = vm::JACOBIAN|vm::KB|vm::POINT;
        const size_type fec_v = gmc_v;
        auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( __cv_id ), __geopc );
        pc_ptrtype pc( new pc_type( this->functionSpace()->fe(), pts ) );
        fec_ptr_t<fec_v,gmc_v> fectx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                             __c,
                                                                             pc );
        DCHECK( ublas::norm_2( __x-__c->xReal(0) ) < 1e-6 ) << "Point " << __x <<  " was not properly found, got " << __c->xReal(0);
        DVLOG(2) << "Point x=" << __x << " and c->xreal = " << __c->xReal(0);
        
        found_pt[ rank ] = 1;

#if defined(FEELPP_HAS_MPI)

        if ( nprocs > 1 && parallel )
        {
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<uint8_type/*int*/>() );
        }

#else
        global_found_pt[ 0 ] = found_pt[ 0 ];
#endif /* FEELPP_HAS_MPI */


        rank_type procIdEval = rank;

        // in case where a point are localised on interprocess face, various process can find this point
        // we take only the process of smaller rank for evaluated the function at point
        if ( nprocs > 1 && parallel )
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

        if ( nprocs > 1 && parallel )
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
        if ( nprocs > 1 && parallel )
        {
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<uint8_type/*int*/>() );
        }

#endif /* FEELPP_HAS_MPI */
        bool found = false;
        rank_type i = 0;

        if ( nprocs > 1 && parallel )
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

            if ( nprocs > 1 && parallel )
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
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::id_( Context_t const & context, id_array_type& v ) const
{
#if 0 // TODO VINCENT
    if ( is_hcurl_conforming && ( ( context.gmContext()->context & vm::KB ) == 0 ) )
        throw std::logic_error("invalid geometric mapping context for hcurl " + std::to_string(context.gmContext()->context) );
    if ( is_hdiv_conforming && ( ( context.gmContext()->context & vm::JACOBIAN ) == 0  ) )
        throw std::logic_error("invalid geometric mapping context for hdiv " + std::to_string(context.gmContext()->context) );
#endif
    index_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_v<index_type> )
        return;

    const uint16_type nq = context.xRefs().size2();
    auto const& s = M_functionspace->dof()->localToGlobalSigns( elt_id );

    for( auto const& ldof : M_functionspace->dof()->localDof( elt_id ) )
    {
        size_type index = ldof.second.index();
        uint16_type local_dof = ldof.first.localDof();
        auto v_ = super::operator[]( index )*s(local_dof);
        for ( uint16_type q = 0; q < nq; ++q )
        {
            v[q] += context.id(local_dof,q)*v_;
        }
    }
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::idInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    // create analysys map : id -> List of pt
    auto __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
    auto it = __loc->result_analysis_begin();
    auto it_end = __loc->result_analysis_end();

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

    auto __c = __gm->template context<DefaultCTX>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<> __ctx = std::make_shared<fec_t<>>( this->functionSpace()->fe(),
                                                   __c,
                                                   __pc );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        auto itL=it->second.begin();
        auto itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );

        //update precompute of basis functions
        __pc->update( pts );

        __c->template update<DefaultCTX>( this->functionSpace()->mesh()->element( it->first ), __geopc );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        id_type __id( this->id( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
        {
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
                for ( typename array_type::index j = 0; j < nComponents2; ++j )
            {
                v[boost::get<0>( *itL )]( i,j ) =  __id( i,j,k );
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

        const size_type gmc_v = vm::JACOBIAN|vm::KB|vm::POINT;
        const size_type fec_v = gmc_v|vm::GRAD;
        auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( __cv_id ), __geopc );
        pc_ptrtype pc( new pc_type( this->functionSpace()->fe(), pts ) );
        fec_ptr_t<fec_v,gmc_v> fectx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                             __c,
                                                                             pc );

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
            LOG( WARNING ) << "no processor seems to have the point " << __x << "\n";
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

    index_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_v<index_type> )
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
            size_type gdof = M_functionspace->dof()->localToGlobal( elt_id, l, c1 ).index();
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
#if 0
                em_fixed_size_matrix_t<nComponents1,nRealDim> mv( v[q].data() );
                em_fixed_size_cmatrix_t<nComponents1,nRealDim> mg( context.grad( ldof, q ).data() );
                mv.noalias()= v_*mg;
#else
                //v[q] = v_*context.grad( ldof, q );
                for ( int k = 0; k < nComponents1; ++k )
                    for ( int j = 0; j < nRealDim; ++j )
                    {
                        v[q]( k,j ) += v_*context.grad( ldof, k, j, q );
                    }
#endif
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
template<typename ContextType,typename EType>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::symmetricGradient( ContextType const & context,
                                                                       grad_array_type& v,
                                                                       std::enable_if_t<EType::is_vectorial>* ) const
{
    index_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_v<index_type> )
        return;

    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = M_functionspace->dof()->localToGlobal( elt_id, l, c1 ).index();

            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            {
                v[q] = v_*context.symmetricGradient( ldof, q );
            }
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::gradInterpolate(  matrix_node_type __ptsReal, grad_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    // create analysys map : id -> List of pt
    auto __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
    auto it = __loc->result_analysis_begin();
    auto it_end = __loc->result_analysis_end();

    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );

    const size_type gmc_v = vm::JACOBIAN|vm::KB|vm::POINT;
    const size_type fec_v = gmc_v|vm::GRAD;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );
    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        auto itL=it->second.begin();
        auto itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );


        //update precompute of basis functions
        __pc->update( pts );
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
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
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dn_( ContextType const & context, dn_array_type& v ) const
{
    index_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_v<index_type> )
        return;

    const size_type Q = context.xRefs().size2();

    auto const& s = M_functionspace->dof()->localToGlobalSigns( elt_id );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            //size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            size_type gdof = M_functionspace->dof()->localToGlobal( elt_id, l, c1 ).index();
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < Q; ++q )
                v[q] += context.dn( ldof, q )*(s(ldof)*v_);
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dnInterpolate( matrix_node_type __ptsReal, dn_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    CHECK( false ) << "TODO";
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

    index_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_v<index_type> )
        return;

    const size_type Q = context.xRefs().size2();

    auto const& s = M_functionspace->dof()->localToGlobalSigns( elt_id );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents:1;
        const int ncdof2 = is_product?nComponents2:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = M_functionspace->dof()->localToGlobal( elt_id, l, c1 ).index();
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );
            value_type v_ = this->globalValue( gdof );
            for ( size_type q = 0; q < Q; ++q )
            {
                v[q] += context.div( ldof, q )*(s(ldof)*v_);
            }
        }
    }
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::divInterpolate( matrix_node_type __ptsReal, div_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    // create analysys map : id -> List of pt
    auto __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
    auto it = __loc->result_analysis_begin();
    auto it_end = __loc->result_analysis_end();

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
    const size_type gmc_v = vm::JACOBIAN|vm::KB|vm::POINT;
    const size_type fec_v = gmc_v|vm::DIV|vm::FIRST_DERIVATIVE;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );
    
    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        auto itL=it->second.begin();
        auto itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
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

    index_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_v<index_type> )
        return;

    auto const& s = M_functionspace->dof()->localToGlobalSigns( elt_id );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = M_functionspace->dof()->localToGlobal( elt_id, l, c1 ).index();
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
#if 0
                v[q] += s(ldof)*v_*context.curl( ldof, q );
#else
                if ( nRealDim == 3 )
                {
                    for ( typename array_type::index i = 0; i < nRealDim; ++i )
                    {
                        v[q]( i,0 ) += s(ldof)*v_*context.curl( ldof, i, 0, q );
                    }
                }

                else if ( nRealDim == 2 )
                {
                    v[q]( 0,0 ) += s(ldof)*v_*context.curl( ldof, 0, 0, q );
                }
#endif
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
    if ( elt_id == invalid_v<index_type> )
        return;

    auto const& s = M_functionspace->dof()->localToGlobalSigns( elt_id );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = M_functionspace->dof()->localToGlobal( elt_id, l, c1 ).index();
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
                if ( nRealDim == 3 )
                {
                    v[q]( 0,0 ) += s(ldof)*v_*context.curl( ldof, comp, 0, q );
                }

                else if ( nRealDim == 2 )
                {
                    v[q]( 0,0 ) += s(ldof)*v_*context.curl( ldof, 0, 0, q );
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
    // create analysys map : id -> List of pt
    auto __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
    auto it = __loc->result_analysis_begin();
    auto it_end = __loc->result_analysis_end();

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
    const size_type gmc_v = vm::JACOBIAN|vm::KB|vm::POINT;
    const size_type fec_v = gmc_v|vm::FIRST_DERIVATIVE|vm::CURL;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        auto itL=it->second.begin();
        auto itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curl_type __curl( this->curl( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nRealDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curl( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curl( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curl( 2,0,k );
            }
        }

        else if ( nRealDim == 2 )
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
    // create analysys map : id -> List of pt
    auto __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
    auto it = __loc->result_analysis_begin();
    auto it_end = __loc->result_analysis_end();

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
    const size_type gmc_v = vm::KB|vm::POINT;
    const size_type fec_v = gmc_v|vm::FIRST_DERIVATIVE|vm::CURL;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        auto itL=it->second.begin();
        auto itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curlx_type __curlx( this->curlx( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nRealDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curlx( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curlx( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curlx( 2,0,k );
            }
        }

        else if ( nRealDim == 2 )
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
    // create analysys map : id -> List of pt
    auto __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
    auto it = __loc->result_analysis_begin();
    auto it_end = __loc->result_analysis_end();

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
    const size_type gmc_v = vm::KB|vm::POINT;
    const size_type fec_v = gmc_v|vm::FIRST_DERIVATIVE|vm::CURL;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        auto itL=it->second.begin();
        auto itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curly_type __curly( this->curly( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nRealDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curly( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curly( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curly( 2,0,k );
            }
        }

        else if ( nRealDim == 2 )
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
    // create analysys map : id -> List of pt
    auto __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
    auto it = __loc->result_analysis_begin();
    auto it_end = __loc->result_analysis_end();

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
    const size_type gmc_v = vm::KB|vm::POINT;
    const size_type fec_v = gmc_v|vm::FIRST_DERIVATIVE|vm::CURL;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        auto itL=it->second.begin();
        auto itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curlz_type __curlz( this->curlz( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nRealDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curlz( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curlz( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curlz( 2,0,k );
            }
        }

        else if ( nRealDim == 2 )
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
            size_type gdof = M_functionspace->dof()->localToGlobal( context.eId(), i, c1 ).index();
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
    // create analysys map : id -> List of pt
    auto __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
    auto it = __loc->result_analysis_begin();
    auto it_end = __loc->result_analysis_end();

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
    const size_type gmc_v = vm::JACOBIAN|vm::KB|vm::POINT;
    const size_type fec_v = gmc_v|vm::GRAD;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );
    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        auto itL=it->second.begin();
        auto itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
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

    typedef typename Localization<mesh_type>::localization_ptrtype localization_ptrtype;
    typedef typename Localization<mesh_type>::container_search_iterator_type analysis_iterator_type;
    typedef typename Localization<mesh_type>::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
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
    const size_type gmc_v = vm::KB|vm::POINT;
    const size_type fec_v = gmc_v|vm::GRAD;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );

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
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
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

    typedef typename Localization<mesh_type>::localization_ptrtype localization_ptrtype;
    typedef typename Localization<mesh_type>::container_search_iterator_type analysis_iterator_type;
    typedef typename Localization<mesh_type>::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
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
    const size_type gmc_v = vm::KB|vm::POINT;
    const size_type fec_v = gmc_v|vm::GRAD;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );
    
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
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
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
            size_type gdof = M_functionspace->dof()->localToGlobal( context.eId(), i, c1 ).index();
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

    typedef typename Localization<mesh_type>::localization_ptrtype localization_ptrtype;
    typedef typename Localization<mesh_type>::container_search_iterator_type analysis_iterator_type;
    typedef typename Localization<mesh_type>::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
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
    const size_type gmc_v = vm::JACOBIAN|vm::KB|vm::POINT|vm::HESSIAN;
    const size_type fec_v = gmc_v|vm::GRAD;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );
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
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
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


//
// Laplacian
//
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::laplacian_( ContextType const & context, id_array_type& v ) const
{
    laplacian_( context, v, mpl::int_<rank>() );
}
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::laplacian_( ContextType const & context, id_array_type& v, mpl::int_<0> ) const
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
            size_type gdof = M_functionspace->dof()->localToGlobal( context.eId(), i, c1 ).index();
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
                v[q]( 0,0 ) += v_*context.laplacian( ldof, 0, 0, q );
            }

        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::laplacian_( ContextType const & context, id_array_type& v, mpl::int_<1> ) const
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
            size_type gdof = M_functionspace->dof()->localToGlobal( context.eId(), i, c1 ).index();
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( i )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
                for ( int c2 = 0; c2 < nRealDim; ++c2 )
                    v[q]( c2,0 ) += v_*context.laplacian( ldof, c2, 0, q );
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::laplacianInterpolate( matrix_node_type __ptsReal, id_array_type& v,
                                                                          bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename Localization<mesh_type>::localization_ptrtype localization_ptrtype;
    typedef typename Localization<mesh_type>::container_search_iterator_type analysis_iterator_type;
    typedef typename Localization<mesh_type>::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_v<index_type>, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_v<index_type> );
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
    const size_type gmc_v = vm::JACOBIAN|vm::KB|vm::POINT|vm::HESSIAN;
    const size_type fec_v = gmc_v|vm::LAPLACIAN|vm::FIRST_DERIVATIVE;
    auto __c = __gm->template context<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );

    fec_ptr_t<fec_v,gmc_v> __ctx = std::make_shared<fec_t<fec_v,gmc_v>>( this->functionSpace()->fe(),
                                                                         __c,
                                                                         __pc );

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
        __c->template update<gmc_v>( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );
        
        //evaluate element for these points
        laplacian_type __lap( this->laplacian( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            v[boost::get<0>( *itL )]( 0,0 ) =  __lap( 0,0,k );
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
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::template Context<geoelement_type> gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef typename gm_type::template Context<geoelement_type> gm_context_type;
    typedef typename geoelement_type::gm1_type gm1_type;
    typedef typename gm1_type::template Context<geoelement_type> gm1_context_type;
    typedef std::shared_ptr<gm_context_type> gm_context_ptrtype;
    typedef std::shared_ptr<gm1_context_type> gm1_context_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm1_context_ptrtype> > map_gmc1_type;


    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

    // basis
    typedef typename element_type::functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context< context, fe_type, gm_type, geoelement_type> fecontext_type;
    typedef std::shared_ptr<fecontext_type> fecontext_ptrtype;
    //typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, fecontext_ptrtype> > map_gmc_type;

    // expression
    //typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;
    typedef typename expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename expression_type::template tensor<map_gmc1_type> t_expr1_type;
    typedef typename t_expr_type::shape shape;

    //
    // start
    //
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

    auto const& initElt = boost::unwrap_ref( *it );
    gm_context_ptrtype __c = this->mesh()->gm()->template context<context>( initElt,__geopc, ex.dynamicContext() );
    gm1_context_ptrtype __c1 = this->mesh()->gm1()->template context<context>( initElt,__geopc1, ex.dynamicContext() );

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
            __c->template update<context>( curElt );
            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
            tensor_expr.update( mapgmc );
            __fe->interpolate( tensor_expr, IhLoc );
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->template update<context>( curElt );
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
                __c->template update<context>( curElt );
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                tensor_expr.update( mapgmc );
                __fe->interpolate( tensor_expr, IhLoc );
            }

            else
            {
                __c1->template update<context>( curElt );
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
    typedef typename boost::unwrap_reference<typename IteratorType::value_type>::type::GeoShape range_geoshape_type;
    typedef Element<Y,Cont> element_type;
    typedef typename element_type::functionspace_type::mesh_type::face_type::GeoShape fe_geoshape_type;

    this->onImpl( r, ex, prefix, geomap_strategy, accumulate, verbose, mpl::int_<MESH_FACES>(),
                  mpl::bool_< std::is_same<range_geoshape_type,fe_geoshape_type>::value >() );
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
                                                            mpl::int_<MESH_FACES>, mpl::true_ )
{
    typedef ExprType expression_type;
    typedef Element<Y,Cont> element_type;

    // mesh element
    typedef typename element_type::functionspace_type::mesh_type::element_type geoelement_type;
    typedef typename element_type::functionspace_type::mesh_type::face_type face_type;
    //typedef typename geoelement_type::face_type face_type;

    // geometric mapping context
    typedef typename geoelement_type::gm_type gm_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    typedef typename geoelement_type::gm1_type gm1_type;
    typedef std::shared_ptr<gm1_type> gm1_ptrtype;

    typedef typename element_type::functionspace_type::fe_type fe_type;
    const size_type context = mpl::if_< mpl::or_<mpl::bool_<is_hdiv_conforming>, mpl::bool_<is_hcurl_conforming> >,
                                        mpl::int_<ExprType::context|vm::POINT|vm::JACOBIAN>,
                                        mpl::int_<ExprType::context|vm::POINT> >::type::value;


    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

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

    auto __gm = this->functionSpace()->gm();
    auto __gm1 = this->functionSpace()->gm1();


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
    auto __c = __gm->template context<context>( firstFace.element( 0 ), __geopc, __face_id, ex.dynamicContext() );
    auto __c1 = __gm1->template context<context>( firstFace.element( 0 ), __geopc1, __face_id, ex.dynamicContext() );

    auto expr_evaluator = ex.evaluator( vf::mapgmc(__c) );
    auto expr1_evaluator = ex.evaluator( vf::mapgmc(__c1) );


    bool hasMeshSupportPartial = __dof->hasMeshSupport() && __dof->meshSupport()->isPartialSupport();
    bool hasDofTableMPIExtended = __dof->buildDofTableMPIExtended();

#if 0
    size_type nbFaceDof = invalid_v<size_type>;

    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );

    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    DVLOG(2)  << "[projector::operator(MESH_FACES)] nbFaceDof = " << nbFaceDof << "\n";
#endif
    auto IhLoc = __fe->faceLocalInterpolant();
    for ( ; __face_it != __face_en; ++__face_it )
    {
        face_type const& curFace = boost::unwrap_ref(*__face_it);

        __face_id = curFace.pos_first();
        uint16_type faceConnectionId = 0;
        if ( hasMeshSupportPartial )
        {
            auto const& elt0 = curFace.element( 0 );
            if ( !__dof->meshSupport()->hasElement( elt0.id() ) || ( !hasDofTableMPIExtended && elt0.isGhostCell() ) )
            {
                if ( !curFace.isConnectedTo1() )
                    continue;
                auto const& elt1 = curFace.element( 1 );
                if ( !__dof->meshSupport()->hasElement( elt1.id() ) || ( !hasDofTableMPIExtended && elt1.isGhostCell() ) )
                    continue;
                __face_id = curFace.pos_second();
                faceConnectionId = 1;
            }
        }
        else if ( !hasDofTableMPIExtended && curFace.element( 0 ).isGhostCell() )
        {
            DCHECK( curFace.isConnectedTo1() ) << "invalid face, no other connection";
            __face_id = curFace.pos_second();
            faceConnectionId = 1;
        }

        DVLOG(2) << "[projector] FACE_ID = " << curFace.id()
                 << " element id= " << ((faceConnectionId == 0)? curFace.ad_first() : curFace.ad_second() )
                 << " pos in elt= " << ((faceConnectionId == 0)? curFace.pos_first() : curFace.pos_second() )
                 << " hasMarker: " << curFace.hasMarker() << "\n";
        DVLOG(2) << "[projector] FACE_ID = " << curFace.id() << " real pts=" << curFace.G() << "\n";


#if 0
        std::pair<size_type,size_type> range_dof( std::make_pair( this->start(),
                this->functionSpace()->nDof() ) );
        DVLOG(2)  << "[projector] dof start = " << range_dof.first << "\n";
        DVLOG(2)  << "[projector] dof range = " << range_dof.second << "\n";
#endif
        switch ( geomap_strategy )
        {
        default:
        case GeomapStrategyType::GEOMAP_OPT:
        case GeomapStrategyType::GEOMAP_HO:
        {
            __c->template update<context>( curFace.element( faceConnectionId ), __face_id );
            DVLOG(2) << "[projector] FACE_ID = " << curFace.id() << "  ref pts=" << __c->xRefs() << "\n";
            DVLOG(2) << "[projector] FACE_ID = " << curFace.id() << " real pts=" << __c->xReal() << "\n";

            expr_evaluator.update( vf::mapgmc( __c ) );
            __fe->faceInterpolate( expr_evaluator, IhLoc );
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->template update<context>( curFace.element( faceConnectionId ), __face_id );
            DVLOG(2) << "[projector] FACE_ID = " << curFace.id() << "  ref pts=" << __c1->xRefs() << "\n";
            DVLOG(2) << "[projector] FACE_ID = " << curFace.id() << " real pts=" << __c1->xReal() << "\n";

            expr1_evaluator.update( vf::mapgmc( __c1 ) );
            __fe->faceInterpolate( expr1_evaluator, IhLoc );
        }
        break;
        }
        if ( accumulate )
            this->plus_assign( curFace, faceConnectionId, IhLoc );
        else
            this->assign( curFace, faceConnectionId, IhLoc );
    } // face_it

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
                                                            mpl::int_<MESH_FACES>, mpl::false_ )
{
    typedef Element<Y,Cont> element_type;

    typedef typename boost::unwrap_reference<typename IteratorType::value_type>::type _range_face_type;
    typedef typename boost::remove_const< typename boost::remove_reference<_range_face_type>::type >::type range_face_type;
    typedef typename range_face_type::super2::template Element<range_face_type>::type range_geoelement_type;
    typedef typename range_geoelement_type::gm_type range_gm_type;
    typedef std::shared_ptr<range_gm_type> range_gm_ptrtype;
    typedef typename range_geoelement_type::permutation_type range_permutation_type;
    typedef typename range_gm_type::precompute_ptrtype range_geopc_ptrtype;
    typedef typename range_gm_type::precompute_type range_geopc_type;

    auto __face_it = r.first;
    auto const __face_en = r.second;

    if ( __face_it == __face_en )
        return;

    auto const* dof = this->functionSpace()->dof().get();
    auto const* fe = this->functionSpace()->fe().get();
    auto meshFe = this->functionSpace()->mesh();
    auto gmFe = this->functionSpace()->mesh()->gm();

    auto const& firstFace = boost::unwrap_ref(*__face_it);

    const bool feMeshIsSubmesh = meshFe->isSubMeshFrom( firstFace.mesh() );
    CHECK ( feMeshIsSubmesh ) << "only implemented for ; fe is submesh";

    auto const& eltConnectedToFirstFace = firstFace.element( 0 );
    uint16_type fid_in_element = firstFace.pos_first();


    range_gm_ptrtype gmRange = eltConnectedToFirstFace.gm();
    if ( !gmRange )
        gmRange.reset( new range_gm_type );

    typedef typename element_type::functionspace_type::basis_type::template ChangeDim<range_geoelement_type::nDim>::type fe_range_dim_type;
    fe_range_dim_type feRangeDim;

    std::vector<std::map<range_permutation_type, range_geopc_ptrtype> > geopcRange( range_geoelement_type::numTopologicalFaces );

    for ( uint16_type __f = 0; __f < range_geoelement_type::numTopologicalFaces; ++__f )
    {
        for ( range_permutation_type __p( range_permutation_type::IDENTITY );
              __p < range_permutation_type( range_permutation_type::N_PERMUTATIONS ); ++__p )
        {
            // special case with lagrange P0d : dof point in 2d/3d faces are not connected to 1d/2d element, trick use barycenter of faces
            if constexpr ( is_lagrange_polynomialset_P0d_v<typename element_type::functionspace_type::basis_type> )
            {
                auto fb = gmRange->referenceConvex().faceBarycenter( __f );
                matrix_node_type fb_convert( fb.size(), 1 );
                ublas::column(fb_convert,0) = fb;
                geopcRange[__f][__p] = range_geopc_ptrtype(  new range_geopc_type( gmRange, fb_convert ) );
            }
            else
                geopcRange[__f][__p] = range_geopc_ptrtype(  new range_geopc_type( gmRange, feRangeDim.points( __f ) ) );
        }
    }

    constexpr size_type context = (is_hdiv_conforming || is_hcurl_conforming)?ExprType::context|vm::POINT|vm::JACOBIAN:ExprType::context|vm::POINT;

    auto gmcRange = gmRange->template context<context,1>( eltConnectedToFirstFace, geopcRange, fid_in_element, ex.dynamicContext() );
    auto expr_evaluator = ex.evaluatorWithPermutation( vf::mapgmc(gmcRange) );

    // geomap context on fe (allow to get relation between geomap context on face range )
    size_type eltIdRelatedToFace = meshFe->meshToSubMesh( firstFace.id() );
    auto const& firstEltRelatedToFace =  this->mesh()->element( eltIdRelatedToFace );
    auto geopcFe = gmFe->preCompute( fe->points() );
    auto gmcFe = gmFe->template context<vm::POINT>( firstEltRelatedToFace, geopcFe );

    CHECK( gmcFe->nPoints() == gmcRange->nPoints() ) << "should be have same number of point : " << gmcFe->nPoints() << " vs " << gmcRange->nPoints();
    uint16_type nPointsGmc = gmcFe->nPoints();
    std::vector<uint16_type> mapBetweenGmc( nPointsGmc, invalid_uint16_type_value );

    auto IhLoc = fe->localInterpolant();
    for ( ; __face_it != __face_en; ++__face_it )
    {
        auto const& curFace = boost::unwrap_ref(*__face_it);
        fid_in_element = curFace.pos_first();
        uint16_type faceConnectionId = 0;
        gmcRange->template update<context>( curFace.element( faceConnectionId ), fid_in_element );
        expr_evaluator.update( vf::mapgmc( gmcRange ) );

        // get dof relation between fe and face in range
        eltIdRelatedToFace = meshFe->meshToSubMesh( curFace.id() );
        auto const& curEltRelatedToFace =  meshFe->element( eltIdRelatedToFace );
        gmcFe->template update<vm::POINT>( curEltRelatedToFace );
        double dofPtCompareTol = std::max(1e-15,curEltRelatedToFace.hMin()*1e-5);
        for ( int q = 0 ; q < nPointsGmc ; ++q )
        {
            mapBetweenGmc[q] = invalid_uint16_type_value;
            auto const& gmcPtFe = gmcFe->xReal(q);
            for ( int q2 = 0 ; q2 < nPointsGmc ; ++q2 )
            {
                auto const& gmcPtRange = gmcRange->xReal(q2);

                bool findDof = true;
                for (uint16_type d=0;d< element_type::functionspace_type::mesh_type::nRealDim;++d)
                {
                    findDof = findDof && (std::abs( gmcPtFe[d]-gmcPtRange[d] )<dofPtCompareTol);
                }
                if ( findDof )
                {
                    //std::cout << "eltId " << curEltRelatedToFace.id() << " :  " << q << "-"<<q2<<"\n";
                    mapBetweenGmc[q] = q2;
                    break;
                }
            }
            CHECK( mapBetweenGmc[q] != invalid_uint16_type_value ) << "not found dof relation";
        }

        // given permutation to tensor expression
        expr_evaluator.setPermutation( mapBetweenGmc );
        // get interpolated value by using the permtutations
        fe->interpolate( expr_evaluator, IhLoc );

        // assign value at dofs
        if ( accumulate )
            this->plus_assign( curEltRelatedToFace, IhLoc );
        else
            this->assign( curEltRelatedToFace, IhLoc );

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
                                                            mpl::int_<MESH_EDGES> )
{
    LOG(INFO) << "onImpl on Mesh Edges";
    // TODO : check that we do not use hdiv hcurl or other type of elements
    const size_type context = ExprType::context|vm::POINT;

    auto mesh = this->functionSpace()->mesh().get();
    auto const* __dof = this->functionSpace()->dof().get();
    auto const* __fe = this->functionSpace()->fe().get();

    // get [first,last) entity iterators over the range
    auto entity_it = r.first;
    auto entity_en = r.second;
    if ( entity_it == entity_en )
        return;

    auto gm = mesh->gm();
    auto const& firstEntity = boost::unwrap_ref( *entity_it );
    DVLOG(3) << "entity " << firstEntity.id() << " with hasMarker "
             << firstEntity.hasMarker() << " nb: " << std::distance(entity_it,entity_en);
    index_type eid = firstEntity.elements().begin()->first;
    uint16_type eid_in_element = firstEntity.elements().begin()->second;
    auto const& elt = mesh->element( eid );
    //auto geopc = gm->preCompute( __fe->edgePoints(eid_in_element) );
    //auto ctx = gm->template context<context,2>( elt, geopc );
    typedef Element<Y,Cont> element_type;
    typedef typename element_type::functionspace_type::mesh_type::element_type geoelement_type;
    typedef typename geoelement_type::template PermutationSubEntity<2> permutation_type;
    typedef typename geoelement_type::gm_type::precompute_ptrtype geopc_ptrtype;
    std::vector<std::map<permutation_type, geopc_ptrtype> > geopc( geoelement_type::numEdges );
    for ( uint16_type __f = 0; __f < geoelement_type::numEdges; ++__f )
    {
        for ( permutation_type __p( permutation_type::IDENTITY ); __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            geopc[__f][__p] = gm->preCompute( __fe->edgePoints(__f) );
        }
    }
    auto ctx = gm->template context<context,2>( elt, geopc, eid_in_element, ex.dynamicContext() );

    auto expr_evaluator = ex.evaluator( vf::mapgmc(ctx) );

    auto IhLoc = __fe->edgeLocalInterpolant();
    for ( ; entity_it != entity_en; ++entity_it )
    {
        auto const& curEntity = boost::unwrap_ref(*entity_it);
        //std::cout << "ok2: " << curEntity.elements().size()  << std::endl;
        eid = invalid_v<index_type>;
        for ( auto const& eltConnectedToEdge : curEntity.elements() )
        {
            index_type eltIdConnected = eltConnectedToEdge.first;
            //std::cout << "eltIdConnected: " << eltIdConnected << std::endl;
            if ( __dof->isElementDone( eltIdConnected ) )
            {
                eid = eltIdConnected;
                eid_in_element = eltConnectedToEdge.second;
                break;
            }
        }
        if ( eid == invalid_v<index_type> )
            continue;
        //std::cout << "entity " << curEntity.id() << " element " << eid << " id in element "
        //<< eid_in_element<< " with hasMarker " << curEntity.hasMarker() << std::endl;
        auto const& elt = mesh->element( eid );
        ctx->template update<context>( elt, eid_in_element );

        expr_evaluator.update( vf::mapgmc( ctx ) );
        __fe->edgeInterpolate( expr_evaluator, IhLoc );
        //std::cout << "Ihloc: " << IhLoc << " eid: " << eid << " eid_in_element:" << eid_in_element << std::endl;
        if ( accumulate )
            this->plus_assign( curEntity, IhLoc, std::make_pair(eid,eid_in_element) );
        else
            this->assign( curEntity, IhLoc, std::make_pair(eid,eid_in_element) );
    } // entity_it
    LOG(INFO) << "onImpl on Mesh Edges done";

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
                                                            mpl::int_<MESH_POINTS> )
{
    LOG(INFO) << "onImpl on Mesh Points";google::FlushLogFiles(google::GLOG_INFO);
    // TODO : check that we do not use hdiv hcurl or other type of elements
    const size_type context = ExprType::context|vm::POINT;
    DVLOG(3)  << "assembling Dirichlet conditions\n";google::FlushLogFiles(google::GLOG_INFO);
    auto mesh = this->functionSpace()->mesh().get();
    auto const* __dof = this->functionSpace()->dof().get();
    auto const* __fe = this->functionSpace()->fe().get();

    // get [first,last) point iterators over the range
    auto pt_it = r.first;
    auto pt_en = r.second;
    if ( pt_it == pt_en )
        return;

    auto gm = mesh->gm();
    auto const& firstPt = boost::unwrap_ref( *pt_it );
    DVLOG(3) << "point " << firstPt.id() << " with hasMarker " << firstPt.hasMarker() << " nb: " << std::distance(pt_it,pt_en);
    google::FlushLogFiles(google::GLOG_INFO);

    index_type eid = firstPt.elements().begin()->first;
    uint16_type ptid_in_element = firstPt.elements().begin()->second;
    auto const& elt = mesh->element( eid );
    //auto geopc = gm->preCompute( __fe->vertexPoints(ptid_in_element) );
    //auto ctx = gm->template context<context>( elt, geopc );

    typedef Element<Y,Cont> element_type;
    typedef typename element_type::functionspace_type::mesh_type::element_type geoelement_type;
    typedef typename geoelement_type::template PermutationSubEntity<geoelement_type::nDim> permutation_type;
    typedef typename geoelement_type::gm_type::precompute_ptrtype geopc_ptrtype;
    std::vector<std::map<permutation_type, geopc_ptrtype> > geopc( geoelement_type::numVertices );
    for ( uint16_type __f = 0; __f < geoelement_type::numVertices; ++__f )
    {
        for ( permutation_type __p( permutation_type::IDENTITY ); __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            geopc[__f][__p] = gm->preCompute( __fe->vertexPoints(__f) );
        }
    }
    auto ctx = gm->template context<context,geoelement_type::nDim>( elt, geopc, ptid_in_element, ex.dynamicContext() );

    //t_expr_type expr( ex, mapgmc(ctx) );
    auto expr_evaluator = ex.evaluator( vf::mapgmc(ctx) );

    size_type nbVertexDof = fe_type::nDofPerVertex;
    DVLOG(3)  << "[projector::operator(MESH_POINTS)] nbVertexDof = " << nbVertexDof << "\n";

    auto IhLoc = __fe->vertexLocalInterpolant();
    for ( ; pt_it != pt_en; ++pt_it )
    {
        auto const& curPt = boost::unwrap_ref(*pt_it);
        DVLOG(3) << "point " << curPt.id() << " with hasMarker " << curPt.hasMarker();
        eid = invalid_v<index_type>;
        for ( auto const& eltConnectedToPoint : curPt.elements() )
        {
            size_type eltIdConnected = eltConnectedToPoint.first;
            if ( __dof->isElementDone( eltIdConnected ) )
            {
                eid = eltIdConnected;
                ptid_in_element = eltConnectedToPoint.second;
                break;
            }
        }
        if ( eid == invalid_v<index_type> )
            continue;

        auto const& elt = mesh->element( eid );
        //geopc = gm->preCompute( __fe->vertexPoints(ptid_in_element) );
        //ctx->update( elt, eid, geopc, mpl::int_<0>() );
        ctx->template update<context>( elt, ptid_in_element );

        expr_evaluator.update( vf::mapgmc( ctx ) );
        __fe->vertexInterpolate( expr_evaluator, IhLoc );
        if ( accumulate )
            this->plus_assign( curPt, IhLoc, std::make_pair(eid,ptid_in_element) );
        else
            this->assign( curPt, IhLoc, std::make_pair(eid,ptid_in_element) );
    } // pt_it
    LOG(INFO) << "onImpl on Mesh Points done";
}


}


#endif
