/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2014-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<int i>
typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::template sub_element<i>::type
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::elementImpl( std::string const& name,
                                                                      bool updateOffViews )
{
    return functionspace_type::legacy_composite_element_ops_type::template buildSubElementView<i>( *this,
                                                                                                  M_containersOffProcess,
                                                                                                  name,
                                                                                                  updateOffViews );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<int i,typename ExprT>
typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::template sub_element<i>::type &
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::element( ExprT e, std::string const& name,
                                                             bool updateOffViews,
                                                             typename std::enable_if<std::is_base_of<ExprBase,ExprT>::value >::type*  )
{
    this->element<i>(name,updateOffViews).on( _range=elements(this->mesh()), _expr=e );
    return this->element<i>(name,updateOffViews);
}

//
// Element implementation
//
template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::Element()
    :
    super(),
    M_start( 0 ),
    M_ct( ComponentType::NO_COMPONENT ),
    M_ct2( ComponentType::NO_COMPONENT ),
    M_containersOffProcess( boost::none )
{
    this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::Element( Element const& __e )
    :
    super( __e ),
    M_backingGuard( __e.M_backingGuard ),
    M_functionspace( __e.M_functionspace ),
    M_name( __e.M_name ),
    M_desc( __e.M_desc ),
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
template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::Element( Element && __e )
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

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
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

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
                                                             std::string const& __name,
                                                             size_type __start,
                                                             ComponentType __ct )
    :
    Element( __functionspace, __name, __name, __start, __ct )
{
}


template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
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

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
                                                             container_type const& __c,
                                                             std::string const& __name,
                                                             size_type __start,
                                                             ComponentType __ct, ComponentType __ct2 )
    :
    Element( __functionspace, __c, __name, __name, __start, __ct, __ct2 )
{}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
                                                             size_type nActiveDof, value_type* arrayActiveDof,
                                                             size_type nGhostDof, value_type* arrayGhostDof )
    :
    super( nActiveDof, arrayActiveDof, nGhostDof, arrayGhostDof, __functionspace->dof() ),
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

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::~Element()
{
    VLOG(1) << "Element destructor...";
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::initFromSpace( functionspace_ptrtype const& __functionspace,
        container_type const& __c )
{
    M_functionspace = __functionspace;
    ( container_type )*this = __c;
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::initSubElementView( mpl::true_ )
{
    fusion::for_each( M_elements,
                      Feel::detail::LegacyCompositeInitializeElement<Element<Y,Cont>>( this ) );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::operator=( Element<Y,Cont> const& __e )
{
    if (  this != &__e )
    {
        M_functionspace = __e.M_functionspace;

        if ( __e.M_name != "unknown" )
            M_name = __e.M_name;
        M_desc = __e.M_desc;

        M_start = __e.M_start;
        M_ct = __e.M_ct;
        M_ct2 = __e.M_ct2;
        M_backingGuard = __e.M_backingGuard;
        M_containersOffProcess = __e.M_containersOffProcess;

        super::operator=( __e );
        this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );
    }

    return *this;
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::operator=( Element<Y,Cont> && __e )
{
    if (  this != &__e )
    {
        M_functionspace = std::move(__e.M_functionspace);
        M_name = std::move(__e.M_name);
        M_desc = std::move(__e.M_desc);
        M_start = __e.M_start;
        M_ct = __e.M_ct;
        M_ct2 = __e.M_ct2;
        M_backingGuard = std::move( __e.M_backingGuard );
        M_containersOffProcess = __e.M_containersOffProcess;

        super::operator=( __e );
        this->initSubElementView( mpl::bool_<functionspace_type::is_composite>() );

        __e.M_functionspace = nullptr;
        __e.M_start = 0;
        __e.M_ct = ComponentType::NO_COMPONENT;
        __e.M_ct2 = ComponentType::NO_COMPONENT;
        __e.M_backingGuard.reset();
        __e.M_containersOffProcess = boost::none;
    }

    return *this;
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename ContOtherType>
typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::operator=( Element<Y,ContOtherType> const& v )
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

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename VectorExpr>
typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::operator=( VectorExpr const& __v )
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
template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::id_type
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::operator()( node_type const& __x, bool extrapolate, bool parallel ) const
{
    this->updateGlobalValues();

    node_type __x_ref;
    size_type __cv_id = invalid_v<size_type>;
    rank_type rank = functionSpace()->mesh()->comm().rank();
    rank_type nprocs = functionSpace()->mesh()->comm().size();
    std::vector<uint8_type/*int*/> found_pt( nprocs, 0 );
    std::vector<uint8_type/*int*/> global_found_pt( nprocs, 0 );

    bool foundPoint = functionSpace()->findPoint( __x, __cv_id, __x_ref );
    if ( !foundPoint && extrapolate && __cv_id != invalid_v<size_type> )
    {
        // findPoint returns the closest candidate, but its last reference
        // coordinates need not belong to that candidate.
        typename gm_type::Inverse inverse( functionSpace()->gm(),
            this->mesh()->element( __cv_id ), this->worldComm().subWorldCommSeqPtr() );
        inverse.setXReal( __x );
        __x_ref = inverse.xRef();
        foundPoint = true;
    }
    if ( foundPoint )
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
template<typename A0, typename A1, typename A2, typename A3>
template<typename Y, typename Cont>
template<typename ContextType, typename Visitor>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::forEachLocalCoefficient(
    ContextType const& context, Visitor&& visit ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    index_type elementId = context.eId();
    auto const& contextMesh = context.gmContext()->element().mesh();
    if ( contextMesh->isSubMeshFrom( this->mesh() ) )
        elementId = contextMesh->subMeshToMesh( elementId );
    if ( contextMesh->isParentMeshOf( this->mesh() ) )
        elementId = this->mesh()->meshToSubMesh( context.eId() );
    if ( elementId == invalid_v<index_type> )
        return;

    for ( auto const& entry : M_functionspace->dof()->localDof( elementId ) )
    {
        auto const localDof = entry.first.localDof();
        visit( localDof, this->orientedLocalCoefficient( elementId, localDof ) );
    }
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y, typename Cont>
template<Feel::size_type Context, typename OutputArray, typename Evaluation>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::interpolateAtPoints(
    matrix_node_type const& points, OutputArray& values, bool conformalEval,
    matrix_node_type const& setPointsConf, Evaluation&& evaluate ) const
{
    // Synchronize before localization: a rank with no local points must also
    // participate in the coefficient exchange.
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    auto localization = this->mesh()->tool_localization();
    if ( conformalEval )
        localization->localize( points, invalid_v<index_type>, setPointsConf, std::true_type{} );
    else
        localization->localize( points, invalid_v<index_type> );
    auto it = localization->resultAnalysisBegin();
    auto const end = localization->resultAnalysisEnd();
    if ( it == end )
        return;

    constexpr bool needsMappingDerivatives =
        vm::has_mapping_hessian_v<Context> ||
        vm::has_laplacian_v<Context> ||
        ( ( is_hdiv_conforming || is_hcurl_conforming ) &&
          ( vm::has_grad_v<Context> || vm::has_first_derivative_v<Context> ) );
    constexpr size_type geometryContext = vm::POINT|vm::JACOBIAN|vm::KB|
        ( needsMappingDerivatives ? vm::HESSIAN : 0 );
    constexpr size_type finiteElementContext = vm::POINT|vm::JACOBIAN|vm::KB|Context;

    auto gm = this->functionSpace()->gm();
    auto referencePoints = matrix_node_type( boost::get<1>( it->second.front() ).size(), it->second.size() );
    auto updateReferencePoints = [&]( auto const& locatedPoints )
    {
        referencePoints.resize( referencePoints.size1(), locatedPoints.size(), false );
        size_type q = 0;
        for ( auto const& point : locatedPoints )
            ublas::column( referencePoints, q++ ) = boost::get<1>( point );
    };
    updateReferencePoints( it->second );
    auto geometricPrecompute = std::make_shared<geopc_type>( gm, referencePoints );
    auto basisPrecompute = std::make_shared<pc_type>( this->functionSpace()->fe(), referencePoints );
    auto geometry = gm->template context<geometryContext>( this->mesh()->element( it->first ), geometricPrecompute );
    auto context = std::make_shared<fec_t<finiteElementContext,geometryContext>>(
        this->functionSpace()->fe(), geometry, basisPrecompute );

    for ( ; it != end; ++it )
    {
        updateReferencePoints( it->second );
        geometricPrecompute->update( referencePoints );
        basisPrecompute->update( referencePoints );
        geometry->template update<geometryContext>( this->mesh()->element( it->first ), geometricPrecompute );
        context->update( geometry, basisPrecompute );
        auto const evaluated = evaluate( *context );

        size_type q = 0;
        for ( auto const& point : it->second )
        {
            auto& value = values[boost::get<0>( point )];
            for ( int i = 0; i < value.dimension( 0 ); ++i )
                for ( int j = 0; j < value.dimension( 1 ); ++j )
                    value( i,j ) = evaluated( i,j,q );
            ++q;
        }
    }
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<Feel::size_type Context, typename OutputArray, typename Evaluation>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::evaluateOnElement(
    geoelement_type const& element, OutputArray& values, Evaluation&& evaluate ) const
{
    constexpr size_type geometryContext = vm::POINT|vm::JACOBIAN|vm::KB|vm::HESSIAN;
    constexpr size_type finiteElementContext = vm::POINT|vm::JACOBIAN|vm::KB|Context;
    auto geometry = this->geomapPtr( element );
    auto precompute = this->pcPtr( element );
    auto context = std::make_shared<fec_t<finiteElementContext,geometryContext>>(
        this->functionSpace()->fe(), geometry, precompute );
    values.resize( boost::extents[geometry->nPoints()] );
    for ( auto& value : values )
        value.setZero();
    evaluate( *context, values );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y, typename Cont>
template<typename Context_t>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::id_( Context_t const & context, id_array_type& v ) const
{
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            v[q] += context.id( localDof, q )*coefficient;
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::idInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<0>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->id( context ); } );
}

//
// Grad
//
template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::grad_type
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::grad( node_type const& __x ) const
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

        constexpr size_type gmc_v = vm::JACOBIAN|vm::KB|vm::POINT|
            ( ( is_hdiv_conforming || is_hcurl_conforming ) ? vm::HESSIAN : 0 );
        constexpr size_type fec_v = vm::JACOBIAN|vm::KB|vm::POINT|vm::GRAD|vm::FIRST_DERIVATIVE;
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
            rank_type owner = 0;
            while ( owner < global_found_pt.size() && global_found_pt[owner] == 0 )
                ++owner;
            mpi::broadcast( functionSpace()->mesh()->comm(), g_, owner );
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
template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::grad_( ContextType const & context, grad_array_type& v ) const
{
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            for ( int c = 0; c < nComponents1; ++c )
                for ( int d = 0; d < nRealDim; ++d )
                    v[q]( c,d ) += coefficient*context.grad( localDof, c, d, q );
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename ContextType,typename EType>
    requires EType::is_vectorial
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::symmetricGradient( ContextType const & context,
                                                                       grad_array_type& v ) const
{
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            v[q] += coefficient*context.symmetricGradient( localDof, q );
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::gradInterpolate(  matrix_node_type __ptsReal, grad_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::GRAD|vm::FIRST_DERIVATIVE>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->grad( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y, typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::symm_gradInterpolate(
    matrix_node_type points, grad_array_type& values, bool conformalEval,
    matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::SYMM|vm::GRAD|vm::FIRST_DERIVATIVE>(
        points, values, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->symm_grad( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename ContextType>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::dn_( ContextType const & context, dn_array_type& v ) const
{
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            v[q] += context.dn( localDof, q )*coefficient;
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::dnInterpolate( matrix_node_type __ptsReal, dn_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    throw std::logic_error(
        "Normal-derivative interpolation requires a face context with an explicit normal; point coordinates alone are insufficient" );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::div_( ContextType const & context, div_array_type& v ) const
{
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            v[q] += context.div( localDof, q )*coefficient;
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::divInterpolate( matrix_node_type __ptsReal, div_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::DIV|vm::FIRST_DERIVATIVE>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->div( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::curl_( ContextType const & context, curl_array_type& v ) const
{
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
        {
            if constexpr ( nRealDim == 3 )
                for ( int c = 0; c < nRealDim; ++c )
                    v[q]( c,0 ) += coefficient*context.curl( localDof, c, 0, q );
            else if constexpr ( nRealDim == 2 )
                v[q]( 0,0 ) += coefficient*context.curl( localDof, 0, 0, q );
        }
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::curl_( ContextType const & context, comp_curl_array_type& v, int comp ) const
{
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
        {
            if constexpr ( nRealDim == 3 )
                v[q]( 0,0 ) += coefficient*context.curl( localDof, comp, 0, q );
            else if constexpr ( nRealDim == 2 )
                v[q]( 0,0 ) += coefficient*context.curl( localDof, 0, 0, q );
        }
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::curlInterpolate( matrix_node_type __ptsReal, curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::CURL|vm::FIRST_DERIVATIVE>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->curl( context ); } );
}


template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::curlxInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::CURL|vm::FIRST_DERIVATIVE>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->curlx( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::curlyInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::CURL|vm::FIRST_DERIVATIVE>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->curly( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::curlzInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::CURL|vm::FIRST_DERIVATIVE>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->curlz( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::d_( int N, ContextType const & context, id_array_type& v ) const
{
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            for ( int c = 0; c < nComponents1; ++c )
                v[q]( c,0 ) += coefficient*context.d( localDof, c, N, q );
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::dxInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::GRAD|vm::FIRST_DERIVATIVE>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->dx( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::dyInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::GRAD|vm::FIRST_DERIVATIVE>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->dy( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::dzInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::GRAD|vm::FIRST_DERIVATIVE>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->dz( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::hess_( ContextType const & context, hess_array_type& v ) const
{
    static_assert( rank <= 1, "Hessian evaluation requires a scalar or vector field" );
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            for ( int i = 0; i < nRealDim; ++i )
                for ( int j = 0; j < nRealDim; ++j )
                    if constexpr ( rank == 0 )
                        v[q]( i,j ) += coefficient*context.hess( localDof, i, j, q );
                    else
                        for ( int c = 0; c < nComponents1; ++c )
                            v[q]( c + nComponents1*i,j ) += coefficient*context.hess( localDof, c, i, j, q );
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::hessInterpolate( matrix_node_type __ptsReal, hess_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::HESSIAN|vm::SECOND_DERIVATIVE|vm::GRAD>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->hess( context ); } );
}


//
// Laplacian
//
template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::laplacian_( ContextType const & context, id_array_type& v ) const
{
    static_assert( rank <= 1, "Laplacian evaluation requires a scalar or vector field" );
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            for ( int c = 0; c < nComponents1; ++c )
                v[q]( c,0 ) += coefficient*context.laplacian( localDof, c, 0, q );
    } );
}


template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::laplacianInterpolate( matrix_node_type __ptsReal, id_array_type& v,
                                                                          bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::LAPLACIAN|vm::SECOND_DERIVATIVE|vm::GRAD>(
        __ptsReal, v, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->laplacian( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y, typename Cont>
template<typename ContextType>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::graddiv_( ContextType const& context, id_array_type& values ) const
{
    static_assert( is_vectorial, "graddiv requires a vector field" );
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            values[q] += coefficient * context.graddiv( localDof, q );
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y, typename Cont>
template<typename ContextType>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::curlcurl_( ContextType const& context, id_array_type& values ) const
{
    static_assert( is_vectorial, "curlcurl requires a vector field" );
    this->forEachLocalCoefficient( context, [&]( auto localDof, auto coefficient )
    {
        for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            values[q] += coefficient * context.curlcurl( localDof, q );
    } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y, typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::graddivInterpolate(
    matrix_node_type points, id_array_type& values, bool conformalEval,
    matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::GRADDIV>( points, values, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->graddiv( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y, typename Cont>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::curlcurlInterpolate(
    matrix_node_type points, id_array_type& values, bool conformalEval,
    matrix_node_type const& setPointsConf ) const
{
    this->template interpolateAtPoints<vm::CURLCURL>( points, values, conformalEval, setPointsConf,
        [this]( auto const& context ) { return this->curlcurl( context ); } );
}

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename IteratorType,  typename ExprType>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::onImpl( std::pair<IteratorType,IteratorType> const& r,
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



    // return if no elements
    if ( it == en )
        return;

    auto const& initElt = boost::unwrap_ref( *it );
    gm_context_ptrtype __c = this->mesh()->gm()->template context<context>( initElt,__geopc, ex.dynamicContext() );
    gm1_context_ptrtype __c1 = this->mesh()->gm1()->template context<context>( initElt,__geopc1, ex.dynamicContext() );

    typedef typename t_expr_type::shape shape;
    constexpr bool is_rank_ok = ( shape::M == nComponents1 &&
                                  shape::N == nComponents2 );

    BOOST_MPL_ASSERT_MSG( is_rank_ok,//mpl::bool_<is_rank_ok>::value,
                          INVALID_TENSOR_RANK,
                          ( mpl::int_<shape::M>, mpl::int_<nComponents>, shape ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    //t_expr_type tensor_expr( basis_type::isomorphism( ex ), mapgmc );
    t_expr_type tensor_expr( ex, mapgmc );

    map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

    t_expr1_type tensor_expr1( ex, mapgmc1 );

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
        if ( isP0Continuous<fe_type>::result ||
             ( is_lagrange_polynomialset_v<fe_type> &&
               std::is_same_v<typename fe_type::continuity_type, Continuous> && __fe->order() == 0 ) )
        {
            break;
        }
    }
} // onImpl (MESH_ELEMENTS)

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename IteratorType,  typename ExprType>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::onImpl( std::pair<IteratorType,IteratorType> const& r,
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

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename IteratorType,  typename ExprType>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::onImpl( std::pair<IteratorType,IteratorType> const& r,
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
    bool hasDofTableMPIExtended = __dof->hasDofTableExtended();

#if 0
    size_type nbFaceDof = invalid_v<size_type>;

    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * __fe->dofPerVertex() +
                      face_type::numEdges * __fe->dofPerEdge() +
                      face_type::numFaces * __fe->dofPerFace() );

    else
        nbFaceDof = face_type::numVertices * __fe->dofPerVertex();

    DVLOG(2)  << "[projector::operator(MESH_FACES)] nbFaceDof = " << nbFaceDof << "\n";
#endif
    auto IhLoc = __fe->faceLocalInterpolant();
    for ( ; __face_it != __face_en; ++__face_it )
    {
        face_type const& curFace = boost::unwrap_ref(*__face_it);
#if 0
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
#else
        uint16_type faceConnectionId = invalid_v<uint16_type>;
        if ( curFace.isConnectedTo1() )
        {
            if ( hasMeshSupportPartial )
            {
                bool hasElt0 = __dof->meshSupport()->hasElement( curFace.element( 0 ).id() );
                bool hasElt1 = __dof->meshSupport()->hasElement( curFace.element( 1 ).id() );
                if ( hasElt0 && hasElt1 )
                {
                    // prefer non ghost element if possible
                    if ( curFace.element( 0 ).isGhostCell() && !curFace.element( 1 ).isGhostCell() )
                        faceConnectionId = 1;
                    else
                        faceConnectionId = 0;
                }
                else if ( hasElt0 )
                    faceConnectionId = 0;
                else if ( hasElt1 )
                    faceConnectionId = 1;
            }
            else
            {
                // prefer non ghost element if possible
                if ( curFace.element( 0 ).isGhostCell() && !curFace.element( 1 ).isGhostCell() )
                    faceConnectionId = 1;
                else
                    faceConnectionId = 0;
            }
        }
        else
        {
            if ( hasMeshSupportPartial )
            {
                if ( __dof->meshSupport()->hasElement( curFace.element( 0 ).id() ) )
                    faceConnectionId = 0;
            }
            else
                faceConnectionId = 0;
        }
        if ( faceConnectionId == invalid_v<uint16_type> )
            continue;
        __face_id = faceConnectionId == 0 ? curFace.pos_first() : curFace.pos_second();
#endif

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

template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename IteratorType,  typename ExprType>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::onImpl( std::pair<IteratorType,IteratorType> const& r,
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


    const bool feMeshIsSubmesh = meshFe->isSubMeshFrom( unwrap_ref(*__face_it).mesh() );
    CHECK ( feMeshIsSubmesh ) << "only implemented for ; fe is submesh";

    size_type eltIdRelatedToFace = meshFe->meshToSubMesh( unwrap_ref(*__face_it).id() );
    while ( eltIdRelatedToFace == invalid_v<size_type> )
    {
        ++__face_it;
        if ( __face_it == __face_en )
            break;
        eltIdRelatedToFace = meshFe->meshToSubMesh( unwrap_ref( *__face_it ).id() );
    }
    if ( eltIdRelatedToFace == invalid_v<size_type> )
        return;

    auto const& firstFace = boost::unwrap_ref(*__face_it);
    auto const& eltConnectedToFirstFace = firstFace.element( 0 );
    uint16_type fid_in_element = firstFace.pos_first();


    range_gm_ptrtype gmRange = eltConnectedToFirstFace.gm();
    if ( !gmRange )
    {
        if constexpr ( range_geoelement_type::is_order_dynamic )
            gmRange = std::make_shared<range_gm_type>( RuntimeOrder{ eltConnectedToFirstFace.order() } );
        else
            gmRange = std::make_shared<range_gm_type>();
    }

    typedef typename element_type::functionspace_type::basis_type::template ChangeDim<range_geoelement_type::nDim>::type fe_range_dim_type;
    auto feRangeDim = [&]() {
        if constexpr ( std::is_constructible_v<fe_range_dim_type, RuntimeOrder> )
            return fe_range_dim_type( RuntimeOrder{ fe->order() } );
        else
            return fe_range_dim_type();
    }();
    bool const isP0d = is_lagrange_polynomialset_v<typename element_type::functionspace_type::basis_type> &&
                       !element_type::functionspace_type::basis_type::continuity_type::is_continuous &&
                       fe->order() == 0;

    std::vector<std::map<range_permutation_type, range_geopc_ptrtype> > geopcRange( range_geoelement_type::numTopologicalFaces );

    for ( uint16_type __f = 0; __f < range_geoelement_type::numTopologicalFaces; ++__f )
    {
        for ( range_permutation_type __p( range_permutation_type::IDENTITY );
              __p < range_permutation_type( range_permutation_type::N_PERMUTATIONS ); ++__p )
        {
            // special case with lagrange P0d : dof point in 2d/3d faces are not connected to 1d/2d element, trick use barycenter of faces
            if ( isP0d )
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
    auto gmcRangeConnection1 = gmRange->template context<context,1>( eltConnectedToFirstFace, geopcRange, fid_in_element, ex.dynamicContext() );
    auto exprOneSideEvaluator = ex.evaluatorWithPermutation( vf::mapgmc(gmcRange) );
    auto exprTwoSideEvaluator = ex.evaluatorWithPermutation( vf::mapgmc(gmcRange,gmcRangeConnection1) );

    // geomap context on fe (allow to get relation between geomap context on face range )
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
        // get dof relation between fe and face in range
        eltIdRelatedToFace = meshFe->meshToSubMesh( curFace.id() );
        auto const& curEltRelatedToFace =  meshFe->element( eltIdRelatedToFace );
        gmcFe->template update<vm::POINT>( curEltRelatedToFace );
        double dofPtCompareTol = std::max(1e-15,curEltRelatedToFace.hMin()*1e-5);

        bool useTwoConnections = curFace.isConnectedTo1();
        uint16_type faceConnection0Id = 0, faceConnection1Id = 1;
        gmcRange->template update<context>( curFace.element( faceConnection0Id ), faceConnection0Id == 0 ? curFace.pos_first() : curFace.pos_second() );
        if ( useTwoConnections)
        {
            gmcRangeConnection1->template update<context>( curFace.element( faceConnection1Id ), faceConnection1Id == 0 ? curFace.pos_first() : curFace.pos_second() );
            exprTwoSideEvaluator.update( vf::mapgmc( gmcRange,gmcRangeConnection1 ) );
        }
        else
        {
            exprOneSideEvaluator.update( vf::mapgmc( gmcRange ) );
        }

        // compute points mapping between fe and face range
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

        if ( useTwoConnections)
        {
            // given permutation to tensor expression
            exprTwoSideEvaluator.setPermutation( mapBetweenGmc );
            // get interpolated value by using the permtutations
            fe->interpolate( exprTwoSideEvaluator, IhLoc );
        }
        else
        {
            // given permutation to tensor expression
            exprOneSideEvaluator.setPermutation( mapBetweenGmc );
            // get interpolated value by using the permtutations
            fe->interpolate( exprOneSideEvaluator, IhLoc );
        }

        // assign value at dofs
        if ( accumulate )
            this->plus_assign( curEltRelatedToFace, IhLoc );
        else
            this->assign( curEltRelatedToFace, IhLoc );

    }
}



template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename IteratorType,  typename ExprType>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::onImpl( std::pair<IteratorType,IteratorType> const& r,
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


template<typename A0, typename A1, typename A2, typename A3>
template<typename Y,  typename Cont>
template<typename IteratorType,  typename ExprType>
void
FunctionSpace<A0, A1, A2, A3>::Element<Y,Cont>::onImpl( std::pair<IteratorType,IteratorType> const& r,
                                                            ExprType const& ex,
                                                            std::string const& prefix,
                                                            GeomapStrategyType geomap_strategy,
                                                            bool accumulate,
                                                            bool verbose,
                                                            mpl::int_<MESH_POINTS> )
{
    LOG(INFO) << "onImpl on Mesh Points";
    Logger::flush();
    // TODO : check that we do not use hdiv hcurl or other type of elements
    const size_type context = ExprType::context|vm::POINT;
    DVLOG(3)  << "assembling Dirichlet conditions\n";
    Logger::flush();

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
    #if !defined(FEELPP_HAS_SPDLOG)

    google::FlushLogFiles(google::GLOG_INFO);

    #else

    Logger::flush();

    #endif

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

    size_type nbVertexDof = __fe->dofPerVertex();
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
