/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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

namespace Feel
{
enum { INTERPOLATE_DIFFERENT_MESH=0, INTERPOLATE_SAME_MESH = 1 };

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
template<typename SpaceType, typename FunctionType>
void
interpolate( boost::shared_ptr<SpaceType> const& space,
             FunctionType const& f,
             typename SpaceType::element_type& interp, int same_mesh = INTERPOLATE_DIFFERENT_MESH )
{
    typedef typename SpaceType::value_type value_type;
    typedef boost::multi_array<value_type,3> array_type;
    typedef typename SpaceType::element_type interp_element_type;

    typedef typename SpaceType::mesh_type mesh_type;
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename mesh_type::element_iterator mesh_element_iterator;

    typedef typename FunctionType::functionspace_type::mesh_type domain_mesh_type;
    typedef typename domain_mesh_type::element_type domain_geoelement_type;
    typedef typename domain_mesh_type::element_iterator domain_mesh_element_iterator;
    // geometric mapping context
    typedef typename mesh_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::template Context<vm::POINT|vm::GRAD|vm::KB|vm::JACOBIAN, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;

    typedef typename domain_mesh_type::gm_type domain_gm_type;
    typedef boost::shared_ptr<domain_gm_type> domain_gm_ptrtype;
    typedef typename domain_gm_type::template Context<vm::POINT|vm::GRAD|vm::KB|vm::JACOBIAN, domain_geoelement_type> domain_gmc_type;
    typedef boost::shared_ptr<domain_gmc_type> domain_gmc_ptrtype;

    typedef typename FunctionType::functionspace_type::fe_type f_fe_type;
    typedef typename f_fe_type::template Context<vm::POINT|vm::GRAD|vm::KB|vm::JACOBIAN, f_fe_type, domain_gm_type, domain_geoelement_type,domain_gmc_type::context> f_fectx_type;
    typedef boost::shared_ptr<f_fectx_type> f_fectx_ptrtype;

    // dof
    typedef typename SpaceType::dof_type dof_type;

    // basis
    typedef typename SpaceType::basis_type basis_type;


    const bool same_basis = boost::is_same<basis_type, typename FunctionType::functionspace_type::basis_type>::value;
    DVLOG(2) << "[interpolate] are the basis the same " << same_basis << "\n";
    DVLOG(2) << "[interpolate] are the meshes the same " << same_mesh << "\n";

    // if same space type and mesh  then return the function itself
    if ( same_basis && same_mesh == INTERPOLATE_SAME_MESH )
    {
        DVLOG(2) << "[interpolate] Same mesh and same space\n";
        interp = f;
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


    f.updateGlobalValues();

    //auto it = f.functionSpace()->mesh()->beginElementWithProcessId();
    //auto en = f.functionSpace()->mesh()->endElementWithProcessId();
    bool upExtendedElt = ( space->mesh()->worldComm().localSize()>1 && f.functionSpace()->dof()->buildDofTableMPIExtended() && space->dof()->buildDofTableMPIExtended() );
    EntityProcessType entityProcess = (upExtendedElt)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    auto rangeElt = elements( f.functionSpace()->mesh(), entityProcess );
    auto it = rangeElt.template get<1>();
    auto en = rangeElt.template get<2>();
    if ( it==en ) return;

    //gmc_ptrtype __c( new gmc_type( __gm, *it, __geopc ) );

    //f.id( *fectx, fvalues );

    // if same mesh but not same function space (different order)
    //if ( f.functionSpace()->mesh() == space->mesh() )
    //if ( same_mesh == INTERPOLATE_SAME_MESH )
    if ( ( MeshBase* )f.functionSpace()->mesh().get() == ( MeshBase* )space->mesh().get() )
    {
        DVLOG(2) << "[interpolate] Same mesh but not same space\n";

        domain_gm_ptrtype __dgm = f.functionSpace()->gm();
        typedef typename domain_gm_type::precompute_ptrtype domain_geopc_ptrtype;
        typedef typename domain_gm_type::precompute_type domain_geopc_type;
        domain_geopc_ptrtype __dgeopc( new domain_geopc_type( __dgm, __basis->dual().points() ) );

        domain_gmc_ptrtype __c( new domain_gmc_type( __dgm, *it, __dgeopc ) );
        auto pc = f.functionSpace()->fe()->preCompute( f.functionSpace()->fe(), __c->xRefs() );

        f_fectx_ptrtype fectx( new f_fectx_type( f.functionSpace()->fe(),
                               __c,
                               pc ) );

        typedef boost::multi_array<typename f_fectx_type::id_type,1> array_type;
        array_type fvalues( f.idExtents( *fectx ) );

        for ( ; it != en; ++ it )
        {
            domain_geoelement_type const& curElt = boost::unwrap_ref(*it);
            __c->update( curElt );
            fectx->update( __c, pc );
            std::fill( fvalues.data(), fvalues.data()+fvalues.num_elements(), f_fectx_type::id_type::Zero() );
            f.id( *fectx, fvalues );

            //std::cout << "interpfunc :  " << interpfunc << "\n";
            for ( uint16_type l = 0; l < basis_type::nLocalDof; ++l )
            {

                const int ncdof = basis_type::is_product?basis_type::nComponents:1;

                for ( uint16_type comp = 0; comp < ncdof; ++comp )
                {
                    size_type globaldof =  boost::get<0>( __dof->localToGlobal( curElt.id(),
                                                          l, comp ) );

#if 0
                    size_type globaldof_f =  boost::get<0>( f.functionSpace()->dof()->localToGlobal( curElt.id(),l, 0 ) );
                    std::cout << "elt : " << curElt.id() << "\n"
                              << "  l : " << l << "\n"
                              << " comp: " << comp << "\n"
                              << " dof: " << globaldof_f << "\n"
                              << "  value: " << f( globaldof_f ) << "\n";
#endif

                    //DVLOG(2) << "globaldof = " << globaldof << " firstldof = " << interp.firstLocalIndex() << " lastldof " << interp.lastLocalIndex() << "\n";
                    // update only values on the processor
                    if ( globaldof >= interp.firstLocalIndex() &&
                            globaldof < interp.lastLocalIndex() )
                    {
                        interp( globaldof ) = fvalues[l]( comp,0 );
                        //DVLOG(2) << "interp( " << globaldof << ")=" << interp( globaldof ) << "\n";
                        //std::cout << "interp( " << globaldof << ")=" << interp( globaldof ) << "\n";
                    }
                }
            }
        }


#if 0
        if ( space->mesh()->worldComm().localSize()>1 )
        {
            if ( f.functionSpace()->dof()->buildDofTableMPIExtended() && space->dof()->buildDofTableMPIExtended() )
            {
                std::set<size_type> eltGhostDone;
                auto face_it = f.functionSpace()->mesh()->interProcessFaces().first;
                auto const face_en = f.functionSpace()->mesh()->interProcessFaces().second;
                for ( ; face_it!=face_en ; ++face_it )
                {
                    auto const& elt0 = face_it->element0();
                    auto const& elt1 = face_it->element1();
                    const bool elt0isGhost = elt0.isGhostCell();
                    auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

                    if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;

                    __c->update( eltOffProc );
                    fectx->update( __c, pc );
                    std::fill( fvalues.data(), fvalues.data()+fvalues.num_elements(), f_fectx_type::id_type::Zero() );
                    f.id( *fectx, fvalues );

                    for ( uint16_type l = 0; l < basis_type::nLocalDof; ++l )
                    {
                        const int ncdof = basis_type::is_product?basis_type::nComponents:1;

                        for ( uint16_type comp = 0; comp < ncdof; ++comp )
                        {
                            size_type globaldof =  boost::get<0>( __dof->localToGlobal( eltOffProc.id(),l, comp ) );

                            // update only values on the processor
                            if ( globaldof >= interp.firstLocalIndex() &&
                                 globaldof < interp.lastLocalIndex() )
                            {
                                interp( globaldof ) = fvalues[l]( comp,0 );
                            }
                        }
                    }
                }
            }
        }
#endif


        DVLOG(2) << "[interpolate] Same mesh but not same space done\n";
    } // same mesh

    else // INTERPOLATE_DIFFERENT_MESH
    {
        DVLOG(2) << "[interpolate] different meshes\n";
        domain_gm_ptrtype __dgm = f.functionSpace()->gm();
        // get only one point
        typename matrix_node<value_type>::type pts( mesh_type::nDim, 1 );
        ublas::column( pts, 0 ) = ublas::column( __basis->dual().points(), 0 );

        domain_geopc_ptrtype __dgeopc( new domain_geopc_type( __dgm, pts ) );

        domain_gmc_ptrtype __c( new domain_gmc_type( __dgm, *it, __dgeopc ) );
        auto pc = f.functionSpace()->fe()->preCompute( f.functionSpace()->fe(), __c->xRefs() );

        f_fectx_ptrtype fectx( new f_fectx_type( f.functionSpace()->fe(),
                               __c,
                               pc ) );
        typedef boost::multi_array<typename f_fectx_type::id_type,1> array_type;
        array_type fvalues( f.idExtents( *fectx ) );


        typename domain_mesh_type::Inverse meshinv( f.functionSpace()->mesh() );

        /* initialisation of the mesh::inverse data structure */
        typename SpaceType::dof_type::dof_points_const_iterator it_dofpt = space->dof()->dofPointBegin();
        typename SpaceType::dof_type::dof_points_const_iterator en_dofpt = space->dof()->dofPointEnd();
        size_type nbpts = 0;

        for ( ; it_dofpt != en_dofpt; ++it_dofpt, ++nbpts  )
        {
            meshinv.addPointWithId( *it_dofpt );
        }

        FEELPP_ASSERT( meshinv.nPoints() == nbpts )( meshinv.nPoints() )( nbpts ).error( "invalid number of points " );
        meshinv.distribute();

        std::vector<bool> dof_done( nbpts );
        std::fill( dof_done.begin(), dof_done.end(), false );
        std::vector<boost::tuple<size_type,uint16_type > > itab;

        size_type first_dof = space->dof()->firstDof();

        for ( ; it != en; ++ it )
        {
            domain_geoelement_type const& curElt = boost::unwrap_ref(*it);
            __c->update( curElt );
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
                    __c->update( *it );
                    pc->update( __c->xRefs() );
                    fectx->update( __c, pc );
                    //typename FunctionType::pc_type pc( f.functionSpace()->fe(), __c->xRefs() );
                    //typename FunctionType::id_type interpfunc( f.id( *__c, pc ) );
                    //typename FunctionType::id_type interpfunc;

                    std::fill( fvalues.data(), fvalues.data()+fvalues.num_elements(), f_fectx_type::id_type::Zero() );
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
    }

    //std::cout << "interp=" << interp << "\n";
} // interpolate

}

#endif /* __interpolate_H */
