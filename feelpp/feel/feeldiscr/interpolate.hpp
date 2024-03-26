/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-01-31

  Copyright (C) 2008-2011 Universite Joseph Fourier (Grenoble I)

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
   \file interpolate.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-01-31
 */
#ifndef __interpolate_H
#define __interpolate_H 1

#include <feel/feelmesh/intersect.hpp>
#include <feel/feelvf/operators.hpp>

namespace Feel
{
enum { INTERPOLATE_DIFFERENT_MESH=0, INTERPOLATE_SAME_MESH = 1 };


template<typename SpaceType, typename FunctionType>
bool
interpolate_copy( std::shared_ptr<SpaceType> const& space,
                  FunctionType const& f,
                  typename SpaceType::element_type& interp )
{
    if constexpr (
#if 0
        std::is_same_v<typename SpaceType::mesh_type, typename FunctionType::functionspace_type::mesh_type> && 
        std::is_same_v<typename SpaceType::basis_type, typename FunctionType::functionspace_type::basis_type>
#else
        std::is_same_v<SpaceType, typename FunctionType::functionspace_type> 
#endif
                  )
    {
        if ( space == f.functionSpace() )
        {
            DVLOG(2) << "[interpolate_copy] Same mesh and same space\n";
            auto n = interp.name();
            interp = f;
            interp.setName( n );
            return true;
        }
    }
    return false;
}

template<typename InterpType, typename SizeT>
void
interpolate_sync( InterpType & interp, bool hasMeshSupportPartialDomain, std::set<SizeT> const& dofUsedWithPartialMeshSupport )
{
    static const bool interp_is_vector = is_std_vector_v<InterpType>;
    if constexpr ( !interp_is_vector )
    {
        if ( hasMeshSupportPartialDomain )
            sync( interp, "=", dofUsedWithPartialMeshSupport );
        else
            sync( interp, "=" );
    }
    else
    {
        for ( auto & interp_c1 : interp )
            for ( auto & interp_c2 : interp_c1 )
                interpolate_sync( unwrap_ptr( interp_c2 ), hasMeshSupportPartialDomain, dofUsedWithPartialMeshSupport );
    }
}
/**
 * Given a space \p space using a lagrange basis, compute the
 * interpolation \p interp of \p f belonging to another function
 * space.
 *
 * <pre>
 * FunctionSpace<mesh_type,fusion::vector<Whatever> > Yh;
 * FunctionSpace<mesh_type,fusion::vector<Whatever> >::element_type f;
 * FunctionSpace<mesh_type,fusion::vector<Lagrange<> > > Xh;
 * FunctionSpace<mesh_type,fusion::vector<Lagrange<> > >::element_type u;
 * interpolate( Xh, f, u );
 * </pre>
 */
template<typename SpaceType, typename FunctionType, typename InterpType>
void
interpolate( std::shared_ptr<SpaceType> const& space,
             FunctionType const& f,
             /*typename SpaceType::element_type*/  InterpType & interp )
{
    static const bool interp_is_vector = is_std_vector_v<InterpType>;
    typedef typename SpaceType::value_type value_type;
    typedef boost::multi_array<value_type,3> array_type;
    typedef typename SpaceType::element_type interp_element_type;

    typedef typename SpaceType::mesh_type mesh_type;
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename mesh_type::element_iterator mesh_element_iterator;
    using size_type = typename SpaceType::size_type;

    typedef typename FunctionType::functionspace_type::mesh_type domain_mesh_type;
    typedef typename domain_mesh_type::element_type domain_geoelement_type;
    typedef typename domain_mesh_type::element_iterator domain_mesh_element_iterator;
    // geometric mapping context
    typedef typename mesh_type::gm_type gm_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    static const size_type gmc_v = vm::POINT|vm::GRAD|vm::KB|vm::JACOBIAN;
    typedef typename gm_type::template Context<geoelement_type> gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;

    typedef typename domain_mesh_type::gm_type domain_gm_type;
    typedef std::shared_ptr<domain_gm_type> domain_gm_ptrtype;
    typedef typename domain_gm_type::template Context<domain_geoelement_type> domain_gmc_type;
    typedef std::shared_ptr<domain_gmc_type> domain_gmc_ptrtype;

    typedef typename FunctionType::functionspace_type::fe_type f_fe_type;
    typedef typename f_fe_type::template Context<vm::POINT|vm::GRAD|vm::KB|vm::JACOBIAN, f_fe_type, domain_gm_type, domain_geoelement_type,gmc_v> f_fectx_type;
    typedef std::shared_ptr<f_fectx_type> f_fectx_ptrtype;

    // dof
    typedef typename SpaceType::dof_type dof_type;

    // basis
    typedef typename SpaceType::basis_type basis_type;

    if constexpr ( !interp_is_vector )
        {
            // if same space type and mesh  then return the function itself
            if ( interpolate_copy( space, f, interp ) )
                return;
        }

    dof_type const* __dof = space->dof().get();
    basis_type const* __basis = space->basis().get();
    gm_ptrtype __gm = space->gm();
    typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
    typedef typename gm_type::precompute_type geopc_type;
    typedef typename domain_gm_type::precompute_ptrtype domain_geopc_ptrtype;
    typedef typename domain_gm_type::precompute_type domain_geopc_type;

    geopc_ptrtype __geopc( new geopc_type( __gm, __basis->dual().points() ) );


    bool inputUseDofTableMPIExtended = f.functionSpace()->dof()->buildDofTableMPIExtended();
    bool outputUseDofTableMPIExtended = space->dof()->buildDofTableMPIExtended();
    bool upExtendedElt = ( space->mesh()->worldComm().localSize()>1 && inputUseDofTableMPIExtended && outputUseDofTableMPIExtended );

    bool applyVectorSync = !upExtendedElt && outputUseDofTableMPIExtended;

    // if same mesh but not same function space (different order)
    if ( f.functionSpace()->mesh()->isSameMesh( space->mesh() ) )
    {
        Range<typename FunctionType::functionspace_type::mesh_type,MESH_ELEMENTS> rangeElt( space->mesh() );
        bool hasMeshSupportPartialDomain = f.functionSpace()->dof()->hasMeshSupport() && f.functionSpace()->dof()->meshSupport()->isPartialSupport();
        bool hasMeshSupportPartialImage = space->dof()->hasMeshSupport() && space->dof()->meshSupport()->isPartialSupport();
        if ( hasMeshSupportPartialDomain && hasMeshSupportPartialImage )
            rangeElt = intersect( f.functionSpace()->dof()->meshSupport()->rangeElements(), space->dof()->meshSupport()->rangeElements() );
        else if ( hasMeshSupportPartialDomain )
            rangeElt = f.functionSpace()->dof()->meshSupport()->rangeElements();
        else if ( hasMeshSupportPartialImage )
            rangeElt = space->dof()->meshSupport()->rangeElements();
        else
        {
            EntityProcessType entityProcess = (upExtendedElt)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
            rangeElt = elements( f.functionSpace()->mesh(), entityProcess );
        }
        // need to sync if the function f is defined on a partial mesh support
        if ( hasMeshSupportPartialDomain )
            applyVectorSync = true;

        auto it = rangeElt.template get<1>();
        auto en = rangeElt.template get<2>();

        std::set<size_type> dofUsedWithPartialMeshSupport;
        if ( it == en )
        {
            if ( applyVectorSync )
            {
                interpolate_sync( interp, hasMeshSupportPartialDomain, dofUsedWithPartialMeshSupport );
            }
            return;
        }

        DVLOG(2) << "[interpolate] Same mesh but not same space";

        if constexpr ( is_lagrange_polynomialset_v<basis_type> && (basis_type::nOrder == 1) &&
                       is_lagrange_polynomialset_v<f_fe_type> && (f_fe_type::nOrder > 0 ) &&
                       !SpaceType::is_mortar && !FunctionType::functionspace_type::is_mortar  )
        {
            // we guess that the vertex local dofs id are the same
            DVLOG(2) << "[interpolate] optimization with Lagrange fe P1";
            for ( ; it != en; ++ it )
            {
                auto const& curElt = unwrap_ref(*it);
                auto const& f_indices = f.functionSpace()->dof()->localToGlobalIndices( curElt.id() );
                for( auto const& ldof : space->dof()->localDof( curElt.id() ) )
                {
                    index_type index = ldof.second.index();
                    dofUsedWithPartialMeshSupport.insert( index );
                    if constexpr ( !interp_is_vector )
                    {
                        uint16_type f_ldofId = ldof.first.localDof();
                        if constexpr ( f_fe_type::nComponents > 1 )
                        {
                            uint16_type comp = f_ldofId/basis_type::nLocalDof;
                            f_ldofId = space->fe()->dofParent( f_ldofId ) + comp*f_fe_type::nLocalDof;
                        }
                        DCHECK( f_ldofId < f_indices.size() ) << "something wrong " << f_ldofId << " vs " << f_indices.size();
                        index_type f_index = f_indices[f_ldofId];
                        interp( index ) = f( f_index );
                    }
                    else
                    {
                        for ( uint16_type c1=0; c1<interp.size(); ++c1 )
                        {
                            for ( uint16_type c2=0; c2<interp[c1].size(); ++c2 )
                            {
                                uint16_type newLocalDofId = ldof.first.localDof()+(c2+f_fe_type::nComponents2*c1)*f_fe_type::nLocalDof;
                                DCHECK( newLocalDofId <  f_indices.size() ) << "something wrong " << newLocalDofId << " vs " << f_indices.size();
                                index_type f_index = f_indices[newLocalDofId];
                                unwrap_ptr(interp[c1][c2])( index ) = f(f_index);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            auto __fe = space->fe();
            auto ex = idv(f);
            const size_type context = ex.context|vm::POINT|vm::KB|vm::JACOBIAN;
            auto gmc = __gm->template context<context>( unwrap_ref(*it), __geopc );
            auto expr_evaluator = ex.evaluator( vf::mapgmc(gmc) );
            if constexpr ( !interp_is_vector )
            {
                auto IhLoc = __fe->localInterpolant();
                for ( ; it != en; ++ it )
                {
                    auto const& curElt = unwrap_ref(*it);
                    gmc->template update<context>( curElt );
                    expr_evaluator.update( vf::mapgmc( gmc ) );
                    __fe->interpolate( expr_evaluator, IhLoc );

                    auto const& s = space->dof()->localToGlobalSigns( curElt.id() );
                    for( auto const& ldof : space->dof()->localDof( curElt.id() ) )
                    {
                        index_type index = ldof.second.index();
                        interp( index ) = s(ldof.first.localDof())*IhLoc( ldof.first.localDof() );
                        dofUsedWithPartialMeshSupport.insert( index );
                    }
                }
            }
            else
            {
                std::vector<std::tuple<uint16_type,uint16_type,decltype(__fe->localInterpolant())>> IhLocsByComp;
                for ( uint16_type c1=0; c1<interp.size(); ++c1 )
                    for ( uint16_type c2=0; c2<interp[c1].size(); ++c2 )
                        IhLocsByComp.push_back( std::make_tuple( c1,c2,__fe->localInterpolant()) );

                for ( ; it != en; ++ it )
                {
                    auto const& curElt = unwrap_ref(*it);
                    gmc->template update<context>( curElt );
                    expr_evaluator.update( vf::mapgmc( gmc ) );

                    for ( auto & [c1,c2,IhLoc] : IhLocsByComp )
                        for( int q = 0; q <  __fe->nLocalDof; ++q )
                            IhLoc( q ) = expr_evaluator.evalq( c1, c2, q );

                    auto const& s = space->dof()->localToGlobalSigns( curElt.id() );
                    for( auto const& ldof : space->dof()->localDof( curElt.id() ) )
                    {
                        index_type gindex = ldof.second.index();
                        uint16_type lindex = ldof.first.localDof();
                        for ( auto const& [c1,c2,IhLoc] : IhLocsByComp )
                            unwrap_ptr(interp[c1][c2])( gindex ) = s(lindex)*IhLoc( lindex );
                        dofUsedWithPartialMeshSupport.insert( gindex );
                    }
                }
            }
        }
        if ( applyVectorSync )
        {
            interpolate_sync( interp, hasMeshSupportPartialDomain, dofUsedWithPartialMeshSupport );
        }

        DVLOG(2) << "[interpolate] Same mesh but not same space done";
    } // same mesh

    else if constexpr ( !interp_is_vector ) // INTERPOLATE_DIFFERENT_MESH
    {
        EntityProcessType entityProcess = (upExtendedElt)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
        auto rangeElt = elements( f.functionSpace()->mesh(), entityProcess );
        auto it = rangeElt.template get<1>();
        auto en = rangeElt.template get<2>();

        if ( it == en )
        {
            if ( applyVectorSync )
                sync( interp, "=" );
            return;
        }


        DVLOG(2) << "[interpolate] different meshes\n";
        domain_gm_ptrtype __dgm = f.functionSpace()->gm();
        // get only one point
        typename matrix_node<value_type>::type pts( mesh_type::nDim, 1 );
        ublas::column( pts, 0 ) = ublas::column( __basis->dual().points(), 0 );

        domain_geopc_ptrtype __dgeopc( new domain_geopc_type( __dgm, pts ) );

        domain_gmc_ptrtype __c = __dgm->template context<gmc_v>( unwrap_ref( *it ), __dgeopc );
        auto pc = f.functionSpace()->fe()->preCompute( f.functionSpace()->fe(), __c->xRefs() );

        f_fectx_ptrtype fectx( new f_fectx_type( f.functionSpace()->fe(),
                               __c,
                               pc ) );
        typedef boost::multi_array<typename f_fectx_type::id_type,1> array_type;
        array_type fvalues( f.idExtents( *fectx ) );

        using index_type = typename domain_mesh_type::index_type;
        MeshInverse<domain_mesh_type> meshinv( f.functionSpace()->mesh() );

        /* initialisation of the mesh::inverse data structure */
        typename SpaceType::dof_type::dof_points_const_iterator it_dofpt = space->dof()->dofPointBegin();
        typename SpaceType::dof_type::dof_points_const_iterator en_dofpt = space->dof()->dofPointEnd();
        size_type nbpts = 0;

        for ( ; it_dofpt != en_dofpt; ++it_dofpt, ++nbpts  )
        {
            meshinv.addPointWithId( it_dofpt->second );
        }

        FEELPP_ASSERT( meshinv.nPoints() == nbpts )( meshinv.nPoints() )( nbpts ).error( "invalid number of points " );
        meshinv.distribute();

        std::vector<bool> dof_done( nbpts );
        std::fill( dof_done.begin(), dof_done.end(), false );
        std::vector<boost::tuple<index_type,uint16_type > > itab;

        size_type first_dof = space->dof()->firstDof();
        typename f_fectx_type::id_type m_id;
        for ( ; it != en; ++ it )
        {
            domain_geoelement_type const& curElt = boost::unwrap_ref(*it);
            __c->template update<gmc_v>( curElt );
            meshinv.pointsInConvex( curElt.id(), itab );

            if ( itab.size() == 0 )
                continue;

            for ( size_type i = 0; i < itab.size(); ++i )
            {
                // get dof id in target dof table
                size_type dof;
                uint16_type comp;
                boost::tie( dof, comp ) = itab[i];
#if !defined( NDEBUG )
                DVLOG(2) << "[interpolate] element : " << curElt.id() << " npts: " << itab.size() << " ptid: " << i
                              << " gdof: " << dof << " comp = " << comp << "\n";
#endif

                if ( !dof_done[dof-first_dof] )
                {
                    dof_done[dof-first_dof]=true;
                    ublas::column( pts, 0 ) = meshinv.referenceCoords().find(dof)->second;
                    __dgeopc->update( pts );
                    //std::cout << "------------------------------------------------------------\n";
                    //std::cout << "pts = " << pts << "\n";
                    __c->template update<gmc_v>( *it );
                    pc->update( __c->xRefs() );
                    fectx->update( __c, pc );
                    //typename FunctionType::pc_type pc( f.functionSpace()->fe(), __c->xRefs() );
                    //typename FunctionType::id_type interpfunc( f.id( *__c, pc ) );
                    //typename FunctionType::id_type interpfunc;

                    std::fill( fvalues.data(), fvalues.data()+fvalues.num_elements(), m_id.constant(0.) );
                    f.id( *fectx, fvalues );
                    //std::cout << "interpfunc :  " << interpfunc << "\n";

                    //for ( uint16_type comp = 0;comp < basis_type::nComponents;++comp )
                    {
                        //size_type globaldof =  basis_type::nLocalDof*comp+first_dof+dof;
                        //size_type globaldof =  first_dof+ndofcomp*comp+dof;
                        //size_type globaldof =  first_dof+dof;
                        size_type globaldof = dof;

                        // update only values on the processor
                        if ( globaldof >= interp.firstLocalIndex() &&
                                globaldof < interp.lastLocalIndex() )
                            interp( globaldof ) = fvalues[0]( comp,0 );

                        //interp( globaldof ) = interpfunc(comp,0,0);
                    }
                }
            }
        }

#if 0

        for ( size_type i = 0; i < dof_done.size(); ++i )
        {
            if ( dof_done[i] != true )
            {
                LOG(INFO) << "[interpolate] dof not treated\n";
                //FEELPP_ASSERT( dof_done[i] == true )( i ).warn ( "invalid dof, was not treated" );

                typename SpaceType::dof_type::dof_points_const_iterator it_dofpt = space->dof()->dofPointBegin();
                typename SpaceType::dof_type::dof_points_const_iterator en_dofpt = space->dof()->dofPointEnd();
                size_type nbpts = 0;

                for ( ; it_dofpt != en_dofpt; ++it_dofpt, ++nbpts  )
                {
                    //meshinv.addPointWithId( *it_dofpt );

                    // be careful with indices in parallel
                    if ( boost::get<1>( *it_dofpt ) == i )
                    {
                        LOG(INFO) << "   id :  " << boost::get<1>( *it_dofpt ) << "\n";
                        LOG(INFO) << "coord :  " << boost::get<0>( *it_dofpt ) << "\n";
                        LOG(INFO) << " comp :  " << boost::get<2>( *it_dofpt ) << "\n";

                        LOG(INFO) << "f( " << boost::get<0>( *it_dofpt ) << ")=" << f( boost::get<0>( *it_dofpt ) ) << "\n";
                    }
                }
            }
        }

#endif

    if ( applyVectorSync )
        sync( interp, "=" );

    }


    //std::cout << "interp=" << interp << "\n";
} // interpolate

}

#endif /* __interpolate_H */
