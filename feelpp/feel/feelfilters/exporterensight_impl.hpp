/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-21

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file exporterensight.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-21
 */
#ifndef __EXPORTERENSIGHT_CPP
#define __EXPORTERENSIGHT_CPP

#include <feel/feelcore/feel.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/timeset.hpp>
#include <feel/feelfilters/exporterensight.hpp>


namespace Feel
{
template<typename MeshType, int N>
ExporterEnsight<MeshType,N>::ExporterEnsight( worldcomm_ptr_t const& worldComm )
:
super( worldComm ),
M_element_type()

{
    init();
}
template<typename MeshType, int N>
ExporterEnsight<MeshType,N>::ExporterEnsight( std::string const& __p, int freq, worldcomm_ptr_t const& worldComm )
    :
    super( "ensight", __p, freq, worldComm ),
    M_element_type()
{
    init();
}
template<typename MeshType, int N>
ExporterEnsight<MeshType,N>::ExporterEnsight( po::variables_map const& vm, std::string const& exp_prefix, worldcomm_ptr_t const& worldComm )
    :
    super( vm, exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
ExporterEnsight<MeshType,N>::ExporterEnsight( std::string const& exp_prefix, worldcomm_ptr_t const& worldComm )
    :
    super( exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
ExporterEnsight<MeshType,N>::ExporterEnsight( ExporterEnsight const & __ex )
    :
    super( __ex ),
    M_element_type( __ex.M_element_type )
{
}

template<typename MeshType, int N>
ExporterEnsight<MeshType,N>::~ExporterEnsight()
{}

template<typename MeshType, int N>
void
ExporterEnsight<MeshType,N>::init()
{
    if ( mesh_type::nDim == 1 )
        if ( mesh_type::Shape == SHAPE_LINE )
            M_element_type = ( mesh_type::nOrder == 1 )?"bar2":"bar3";

    if ( mesh_type::nDim == 2 )
    {
        if ( mesh_type::Shape == SHAPE_TRIANGLE )
            M_element_type = ( mesh_type::nOrder == 1 )?"tria3":"tria6";

        else if ( mesh_type::Shape == SHAPE_QUAD )
            M_element_type = ( mesh_type::nOrder == 1 )?"quad4":"quad8";
    }

    if ( mesh_type::nDim == 3 )
    {
        if ( mesh_type::Shape == SHAPE_TETRA )
            M_element_type = ( mesh_type::nOrder == 1 )?"tetra4":"tetra10";

        else if ( mesh_type::Shape == SHAPE_HEXA )
            M_element_type = ( mesh_type::nOrder == 1 )?"hexa8":"hexa20";
    }
}
template<typename MeshType, int N>
void
ExporterEnsight<MeshType,N>::save( steps_write_on_disk_type const& stepsToWriteOnDisk ) const
{
    DVLOG(2) << "[ExporterEnsight::save] export in ensight format\n";

    tic();
    _F_writeSoSFile();
    toc("ExporterEnsight::save sos",FLAGS_v>1);

    tic();
    _F_writeCaseFile();
    toc("ExporterEnsight::save case",FLAGS_v>1);

    tic();
    _F_writeGeoFiles();
    toc("ExporterEnsight::save geo",FLAGS_v>1);

    tic();
    _F_writeVariableFiles();
    toc("ExporterEnsight::save variable",FLAGS_v>1);
}

template<typename MeshType, int N>
void
ExporterEnsight<MeshType,N>::_F_writeSoSFile() const
{
    // only on proc 0
    if ( this->worldComm().rank() == this->worldComm().masterRank() )
    {
        // first save for paraview
        {
            std::ostringstream filestr;
            filestr << this->path() << "/" << this->prefix() << "-paraview-" << this->worldComm().globalSize() << ".sos";
            std::ofstream __out( filestr.str().c_str() );
            
            if ( __out.fail() )
            {
                DVLOG(2) << "cannot open " << filestr.str()  << "\n";
                exit( 0 );
            }
            
            __out << "FORMAT:\n"
                << "type: master_server gold \n"
                << "SERVERS\n"
                << "number of servers: " << this->worldComm().globalSize() << "\n";
            
            for ( int pid = 0 ; pid < this->worldComm().globalSize(); ++pid )
            {
                
                __out << "#Server " << pid+1 << "\n"
                    << "machine id: " << mpi::environment::processor_name() << "\n"
                    << "executable: /usr/local/bin/ensight76/bin/ensight7.server\n"
                    << "data_path: " << fs::current_path().string() << "\n"
                    << "casefile: " << this->prefix() << "-" << this->worldComm().globalSize() << "_" << pid << ".case\n";
            }
            __out.close();
        }
        {
            // second save for ensight
            std::ostringstream filestr;
            filestr << this->path() << "/" << this->prefix() << "-" << this->worldComm().globalSize() << ".sos";
            std::ofstream __out( filestr.str().c_str() );
        
            if ( __out.fail() )
            {
                DVLOG(2) << "cannot open " << filestr.str()  << "\n";
                exit( 0 );
            }
            __out << "FORMAT:\n"
            << "type: master_server gold \n\n"
            << "MULTIPLE_CASEFILES\n"
            << "total number of cfiles: " << this->worldComm().globalSize() << "\n"
            << "# cfiles global path: " << fs::current_path().string() << "\n"
            << "cfiles pattern: "<<this->prefix() << "-" << this->worldComm().globalSize() << "_*.case\n"
            << "cfiles start number: 0\n"
            << "cfiles increment: 1\n\n"
            << "SERVERS\n"
            << "number of servers: "<< (this->worldComm().globalSize()/100)+1 <<" repeat\n";
            __out.close();
        }
    }
}
template<typename MeshType, int N>
void
ExporterEnsight<MeshType,N>::_F_writeCaseFile() const
{
    std::ostringstream filestr;
    filestr << this->path() << "/"
            << this->prefix() << "-"
            << this->worldComm().globalSize() << "_" << this->worldComm().globalRank() << ".case";
    std::ofstream __out( filestr.str().c_str() );

    if ( __out.fail() )
    {
        DVLOG(2) << "cannot open " << filestr.str()  << "\n";
        exit( 0 );
    }

    __out << "FORMAT:\n"
          << "type: ensight \n"
          << "GEOMETRY:\n";

    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();

    switch ( this->exporterGeometry() )
    {
    case EXPORTER_GEOMETRY_STATIC:
    {
        timeset_ptrtype __ts = *__ts_it;
        __out << "model: " << this->prefix()
              << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank() << ".geo";
    }
    break;
    default:
    case EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY:
    case EXPORTER_GEOMETRY_CHANGE:
    {
        while ( __ts_it != __ts_en )
        {
            timeset_ptrtype __ts = *__ts_it;

            __out << "model: " << __ts->index() << " " << this->prefix()
                  << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank() << ".geo***";
            if ( this->exporterGeometry() == EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY )
                __out << " change_coords_only";

            ++__ts_it;
        }
    }
    break;
    }
    __out << "\n";

    __out << "VARIABLES:" << "\n";

    __ts_it = this->beginTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;

        auto __tstp_st = __ts->beginStep();
        auto __tstp_en = __ts->endStep();

        /* protect this portion of code */
        /* if we don't have time steps */
        /* happens when we only have the mesh */
        if(__tstp_st != __tstp_en)
        {
            auto __tstp_it = boost::prior(__tstp_en);

            typename timeset_type::step_type::nodal_const_iterator __it = ( *__tstp_it )->beginNodal();
            typename timeset_type::step_type::nodal_const_iterator __end = ( *__tstp_it )->endNodal();
            while ( __it != __end )
            {
                auto const& nodalData =  __it->second;
                switch ( nodalData.first )
                {
                case FunctionSpaceType::SCALAR :  __out << "scalar";break;
                case FunctionSpaceType::VECTORIAL :  __out << "vector";break;
                case FunctionSpaceType::TENSOR2 :  __out << "tensor asym";break;
                case FunctionSpaceType::TENSOR2_SYMM : __out << "tensor symm";break;
                }
                auto const& nodalField = unwrap_ptr( nodalData.second[0][0] );
                __out << " per node: "
                    << __ts->index() << " " // << *__ts_it->beginStep() << " "
                    << nodalField.name() << " " << this->prefix() << "." << __it->first << "-" << this->worldComm().globalSize() << "_" << nodalField.worldComm().localRank() << ".***" << "\n";// important localRank !!
                ++__it;
            }

            typename timeset_type::step_type::element_const_iterator __it_el = ( *__tstp_it )->beginElement();
            typename timeset_type::step_type::element_const_iterator __end_el = ( *__tstp_it )->endElement();
            while ( __it_el != __end_el )
            {
                auto const& elementData =  __it_el->second;
                switch ( elementData.first )
                {
                case FunctionSpaceType::SCALAR :  __out << "scalar";break;
                case FunctionSpaceType::VECTORIAL :  __out << "vector";break;
                case FunctionSpaceType::TENSOR2 :  __out << "tensor asym";break;
                case FunctionSpaceType::TENSOR2_SYMM :  __out << "tensor symm";break;
                }
                auto const& elementField = unwrap_ptr( elementData.second[0][0] );
                __out << " per element: "
                    << __ts->index() << " " // << *__ts_it->beginStep() << " "
                    << elementField.name() << " " << this->prefix() << "." << __it_el->first << "-" << this->worldComm().globalSize() << "_" << elementField.worldComm().localRank() << ".***" << "\n";// important localRank !!
                ++__it_el;
            }

        }

        ++__ts_it;
    }

    __out << "TIME:\n";
    __ts_it = this->beginTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;
        typename timeset_type::step_const_iterator __its = __ts->beginStep();
        typename timeset_type::step_const_iterator __ens = __ts->endStep();

        if(__its != __ens)
        {
            __out << "time set:        " << __ts->index() << "\n"
                << "number of steps: " << __ts->numberOfSteps() << "\n"
                << "filename start number: " << ( *__its )->index() << "\n"
                << "filename increment: " << 1 << "\n"
                << "time values: ";
        }
        else
        {
            __out << "time set:        " << TS_INITIAL_INDEX << "\n"
                << "number of steps: " << 1 << "\n"
                << "filename start number: " << TS_INITIAL_INDEX << "\n"
                << "filename increment: " << 1 << "\n"
                << "time values: 1.0";
        }

        uint16_type __l = 0;

        while ( __its != __ens )
        {

            __out << ( *__its )->time() << " ";

            if ( __l++ % 10 == 0 )
                __out << "\n";

            ++__its;
        }

#if 0
        namespace lambda = boost::lambda;
        std::for_each( __ts->beginStep(), __ts->endStep(),
                       __out << lambda::bind( &timeset_type::step_type::time, *lambda::_1 ) << boost::lambda::constant( ' ' ) );
        std::for_each( __ts->beginStep(), __ts->endStep(),
                       std::cerr << lambda::bind( &timeset_type::step_type::time, *lambda::_1 ) << boost::lambda::constant( ' ' ) );
#endif
        ++__ts_it;
    }

    __out << "\n";
    __out.close();

}

template<typename MeshType, int N>
void
ExporterEnsight<MeshType,N>::_F_writeGeoFiles() const
{
    namespace lambda = boost::lambda;

    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;

        typename timeset_type::step_const_iterator __it = __ts->beginStep();
        typename timeset_type::step_const_iterator __end = __ts->endStep();
        int timeIndex = TS_INITIAL_INDEX;
        mesh_ptrtype mesh = NULL;

        /* if we do not have time steps, we try to save the mesh at least */
        if( __it == __end )
        {
            if ( __ts->hasMesh() )
            { mesh = __ts->mesh(); }
        }
        /* otherwise we save the mesh */
        else
        {
            __it = boost::prior( __end );

            /* check that step is in memory */
            if( (*__it)->isInMemory() && (*__it)->hasMesh() )
            { mesh = (*__it)->mesh(); }

            /* record step index */
            timeIndex = (*__it)->index();
        }

        /* if we were not able to get a mesh, there is a problem */
        if(!mesh)
        {
            LOG(INFO) << "GEO: Unable to get mesh data" << std::endl;
            return;
        }

        if ( this->exporterGeometry() == EXPORTER_GEOMETRY_STATIC )
        {
            std::ostringstream __geofname;
            __geofname << this->path() << "/"
                       << this->prefix()
                       << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank()
                       << ".geo";
            M_filename =  __geofname.str();
            CHECK( (*__it)->hasMesh() || __ts->hasMesh()  ) << "Invalid mesh data structure in static geometry mode\n";
            mesh->accept( const_cast<ExporterEnsight<MeshType,N>&>( *this ) );
        }
        else
        {
            std::ostringstream __geofname;

            __geofname << this->path() << "/"
                << this->prefix()
                << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank()
                << ".geo" << std::setfill( '0' ) << std::setw( 3 ) << timeIndex;

            //__writegeo( __step->mesh(), __ts->name(), __geofname.str() );
            //, __ts->name(), __geofname.str() );
            M_filename =  __geofname.str();
            mesh->accept( const_cast<ExporterEnsight<MeshType,N>&>( *this ) );
        }

        ++__ts_it;
    }
}

template<typename MeshType, int N>
void
ExporterEnsight<MeshType,N>::_F_writeVariableFiles() const
{
    namespace lambda = boost::lambda;


    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;

        typename timeset_type::step_const_iterator __it = __ts->beginStep();
        typename timeset_type::step_const_iterator __end = __ts->endStep();

        if(__it != __end)
        {
            __it = boost::prior( __end );

            while ( __it != __end )
            {
                typename timeset_type::step_ptrtype __step = *__it;;

                if ( __step->isInMemory() )
                {
                    saveFields<true>( __step, __step->beginNodal(), __step->endNodal() );
                    saveFields<false>( __step, __step->beginElement(), __step->endElement() );
                }

                ++__it;
            }
        }

        ++__ts_it;
    }
}


template<typename MeshType, int N>
template<bool IsNodal,typename Iterator>
void
ExporterEnsight<MeshType,N>::saveFields( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const
{
    //auto mit = elements(__step->mesh());
    //Feel::detail::MeshPoints<float> mp( __step->mesh().get(), this->worldComm(), mit.template get<1>(), mit.template get<2>(), false, true );
    while ( __var != en )
    {
        auto const& fieldData =  __var->second;
        auto const& field00 = unwrap_ptr( fieldData.second[0][0] );

        if ( !field00.worldComm().isActive() ) return;

        std::ostringstream __varfname;

        __varfname << this->path() << "/" << this->prefix() << "." << __var->first
                   << "-" << this->worldComm().globalSize() << "_" << field00.worldComm().localRank() // important localRank
                   << "." << std::setfill( '0' ) << std::setw( 3 ) << __step->index();
        DVLOG(2) << "[ExporterEnsight::saveFields] saving " << __varfname.str() << "...\n";
        std::fstream __out( __varfname.str().c_str(), std::ios::out | std::ios::binary );

        char buffer[ 80 ];
        strcpy( buffer, field00.name().c_str() );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );

        uint16_type nComponents = invalid_uint16_type_value, nComponents1 = invalid_uint16_type_value, nComponents2 = invalid_uint16_type_value;
        bool isTensor2Symm = false;
        if ( fieldData.first == FunctionSpaceType::SCALAR )
        {
            nComponents = 1;
            nComponents1 = 1;
            nComponents2 = 1;
        }
        else if ( fieldData.first == FunctionSpaceType::VECTORIAL )
        {
            nComponents = 3;
            nComponents1 = 3;
            nComponents2 = 1;
        }
        else if ( fieldData.first == FunctionSpaceType::TENSOR2 )
        {
            nComponents = 9;
            nComponents1 = 3;
            nComponents2 = 3;
        }
        else if ( fieldData.first == FunctionSpaceType::TENSOR2_SYMM )
        {
            nComponents = 6;
            nComponents1 = 3;
            nComponents2 = 3;
            isTensor2Symm = true;
        }

        /**
         * BE CAREFUL HERE some points in the mesh may not be present in the
         * mesh element connectivity, we really need to have an array of the
         * size of the number of points in the mesh even though if some are not
         * in the connectivity and not an array of the dimension of the function
         * space which has the "right" size.
         */
        Eigen::VectorXf m_field;
        auto const& d = field00.functionSpace()->dof().get();
        int reorder_tensor2symm[6] = { 0,3,4,1,5,2 };

        if constexpr ( IsNodal )
            {
                if (!M_mapNodalArrayToDofId )
                {
                    auto rangeElements = __step->mesh()->elementsWithProcessId();
                    auto elt_it = std::get<0>( rangeElements );
                    auto elt_en = std::get<1>( rangeElements );

                    M_cache_mp.try_emplace(  0 /*dummy value*/, __step->mesh().get(), this->worldComm(), elt_it, elt_en, false, true );
                    auto& mp = M_cache_mp.at(  0 /*dummy value*/ );

                    index_type nValuesPerComponent = mp.ids.size();
                    m_field = Eigen::VectorXf::Zero( nComponents*nValuesPerComponent );

                    M_mapNodalArrayToDofId = std::make_optional<std::vector<size_type>>( nValuesPerComponent,invalid_size_type_value );
                    auto & mapArrayToDofId = *M_mapNodalArrayToDofId;

                    for ( ; elt_it != elt_en; ++elt_it )
                    {
                        auto const& elt = unwrap_ref( *elt_it );
                        auto const& locglob_ind = d->localToGlobalIndices( elt.id() );
                        for ( uint16_type p = 0; p < __step->mesh()->numLocalVertices(); ++p )
                        {
                            size_type ptid = mp.old2new[elt.point( p ).id()]-1;
                            size_type dof_id = locglob_ind(d->localDofId(p,0));
                            mapArrayToDofId[ptid] = dof_id;
                            for ( uint16_type c1 = 0; c1 < fieldData.second.size(); ++c1 )
                            {
                                for ( uint16_type c2 = 0; c2 < fieldData.second[c1].size(); ++c2 )
                                {
                                    auto const& fieldComp = unwrap_ptr( fieldData.second[c1][c2] );
                                    uint16_type cMap = c2*nComponents1+c1;
                                    if ( isTensor2Symm )
                                        cMap = reorder_tensor2symm[Feel::detail::symmetricIndex( c1,c2, nComponents1 )];
                                    size_type global_node_id = nComponents * ptid + cMap;
                                    m_field[global_node_id] = fieldComp.globalValue( dof_id );
                                }
                            }
                        }
                    }
                }
                else
                {
                    auto const& mapArrayToDofId = *M_mapNodalArrayToDofId;
                    index_type nValuesPerComponent = mapArrayToDofId.size();
                    m_field = Eigen::VectorXf::Zero( nComponents*nValuesPerComponent );
                    for ( size_type k=0;k<nValuesPerComponent;++k )
                    {
                        size_type dof_id = mapArrayToDofId[k];
                        for ( uint16_type c1 = 0; c1 < fieldData.second.size(); ++c1 )
                        {
                            for ( uint16_type c2 = 0; c2 < fieldData.second[c1].size(); ++c2 )
                            {
                                auto const& fieldComp = unwrap_ptr( fieldData.second[c1][c2] );
                                uint16_type cMap = c2*nComponents1+c1;
                                if ( isTensor2Symm )
                                    cMap = reorder_tensor2symm[Feel::detail::symmetricIndex( c1,c2, nComponents1 )];
                                size_type global_node_id = nComponents*k + cMap;
                                m_field[global_node_id] = fieldComp.globalValue( dof_id );
                            }
                        }
                    }
                }
                __out.write( ( char * ) m_field.data(), m_field.size() * sizeof( float ) );
            }
        else
        {
            for ( auto const& [fragmentId,fragmentData] : fragmentationMarkedElements( __step->mesh() ) )
            {
                auto const& [range,mIds,fragmentName] = fragmentData;
                sprintf( buffer, "part %d",fragmentId );
                __out.write( ( char * ) & buffer, sizeof( buffer ) );
                DVLOG(2) << "part " << buffer << "\n";
                strcpy( buffer, this->elementType().c_str() );
                __out.write( ( char * ) & buffer, sizeof( buffer ) );
                DVLOG(2) << "element type " << buffer << "\n";

                auto itFindMapElementArrayToDofId = M_mapElementArrayToDofId.find(fragmentId);
                if ( itFindMapElementArrayToDofId == M_mapElementArrayToDofId.end() )
                {
                    auto rangeMarkedElements = range;
                    auto elt_m_it = boost::get<1>( rangeMarkedElements );
                    auto const elt_m_en = boost::get<2>( rangeMarkedElements );
                    index_type nValuesPerComponent = std::distance( elt_m_it, elt_m_en );
                    m_field = Eigen::VectorXf::Zero( nComponents*nValuesPerComponent );
                    auto & mapArrayToDofId =  M_mapElementArrayToDofId[fragmentId];
                    mapArrayToDofId.resize( nValuesPerComponent,invalid_size_type_value );

                    for ( index_type e=0; elt_m_it != elt_m_en; ++elt_m_it,++e )
                    {
                        auto const& elt = unwrap_ref( *elt_m_it );
                        auto const& locglob_ind = d->localToGlobalIndices( elt.id() );
                        size_type dof_id = locglob_ind(d->localDofId(0,0));
                        mapArrayToDofId[e] = dof_id;
                        for ( uint16_type c1 = 0; c1 < fieldData.second.size(); ++c1 )
                        {
                            for ( uint16_type c2 = 0; c2 < fieldData.second[c1].size(); ++c2 )
                            {
                                auto const& fieldComp = unwrap_ptr( fieldData.second[c1][c2] );
                                uint16_type cMap = c2*nComponents1+c1;
                                if ( isTensor2Symm )
                                    cMap = reorder_tensor2symm[Feel::detail::symmetricIndex( c1,c2, nComponents1 )];
                                size_type global_node_id = nComponents * e + cMap;
                                m_field(global_node_id) = fieldComp.globalValue( dof_id );
                            }
                        }
                    }
                }
                else
                {
                    auto const& mapArrayToDofId = itFindMapElementArrayToDofId->second;
                    index_type nValuesPerComponent = mapArrayToDofId.size();
                    m_field = Eigen::VectorXf::Zero( nComponents*nValuesPerComponent );
                    for ( size_type k=0;k<nValuesPerComponent;++k )
                    {
                        size_type dof_id = mapArrayToDofId[k];
                        for ( uint16_type c1 = 0; c1 < fieldData.second.size(); ++c1 )
                        {
                            for ( uint16_type c2 = 0; c2 < fieldData.second[c1].size(); ++c2 )
                            {
                                auto const& fieldComp = unwrap_ptr( fieldData.second[c1][c2] );
                                uint16_type cMap = c2*nComponents1+c1;
                                if ( isTensor2Symm )
                                    cMap = reorder_tensor2symm[Feel::detail::symmetricIndex( c1,c2, nComponents1 )];
                                size_type global_node_id = nComponents*k + cMap;
                                m_field[global_node_id] = fieldComp.globalValue( dof_id );
                            }
                        }
                    }
                }
                __out.write( ( char * ) m_field.data(), m_field.size() * sizeof( float ) );
            }
        }

        DVLOG(2) << "[ExporterEnsight::saveFields] saving " << __varfname.str() << "done\n";
        ++__var;
    }
}

template<typename MeshType, int N>
void
ExporterEnsight<MeshType,N>::visit( mesh_type* __mesh )
{
    char buffer[ 80 ];
    std::vector<int> idnode, idelem;

    std::fstream __out( M_filename.c_str(), std::ios::out | std::ios::binary );

    // get only the filename (maybe with full path)
    fs::path gp = M_filename;
    std::string theFileName = gp.filename().string();
    CHECK( theFileName.length() <= 80 ) << "the file name is too long : theFileName=" << theFileName << "\n";

    strcpy( buffer, "C Binary" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, theFileName.c_str() );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "elements" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "node id given" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "element id given" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "coordinates" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );

    auto mit = elements(__mesh);
    M_cache_mp.try_emplace( 0 /*dummy value*/,  __mesh, this->worldComm(), mit.template get<1>(), mit.template get<2>(), false, true );
    auto& mp = M_cache_mp.at(  0 /*dummy value*/ );
    size_type __nv = mp.ids.size();
    __out.write( ( char * ) &__nv, sizeof( int ) );
    LOG(INFO) << "n pts = " << __nv << " numppoints=" << __mesh->numPoints();
    __out.write( ( char * ) & mp.ids.front(), mp.ids.size() * sizeof( int ) );
    __out.write( ( char * ) mp.coords.data(), mp.coords.size() * sizeof( float ) );

    for ( auto const& [fragmentId,fragmentData] : fragmentationMarkedElements( __mesh ) )
    {
        auto const& [range,mIds,fragmentName] = fragmentData;

        sprintf( buffer, "part %d",fragmentId );
        //    strcpy( buffer, "part 1" );

        __out.write( ( char * ) & buffer, sizeof( buffer ) );
        sprintf( buffer, "Marker %d (%s)", fragmentId, fragmentName.substr(0, 32).c_str());

        __out.write( ( char * ) & buffer, sizeof( buffer ) );

        strcpy( buffer, this->elementType().c_str() );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );

        auto rangeMarkedElements = range;
        auto elt_it = boost::get<1>( rangeMarkedElements );
        auto const elt_en = boost::get<2>( rangeMarkedElements );

        //	int __ne = __mesh->numElements();
        //int __ne = p_it->second;
        int __ne = std::distance( elt_it, elt_en );

        DVLOG(2) << "num Elements to save : " << __ne << "\n";

        __out.write( ( char * ) &__ne, sizeof( int ) );

        idelem.resize( __ne );


        for ( size_type e = 0; elt_it != elt_en; ++elt_it, ++e )
        {
            idelem[e] = boost::unwrap_ref( *elt_it ).id() + 1;
        }

        __out.write( ( char * ) & idelem.front(), idelem.size() * sizeof( int ) );

        std::vector<int> eids( __mesh->numLocalVertices()*__ne );
        size_type e= 0;
        elt_it = boost::get<3>( rangeMarkedElements )->begin();
        for ( ; elt_it != elt_en; ++elt_it, ++e )
        {
            auto const& elt = boost::unwrap_ref( *elt_it );
            for ( size_type j = 0; j < __mesh->numLocalVertices(); j++ )
            {
                // ensight id start at 1
                int __id = mp.old2new[ elt.point( j ).id()];
                eids[__mesh->numLocalVertices()*e+j] = __id;
            }
        }
        __out.write( ( char * ) eids.data() , eids.size()*sizeof( int ) );

    }
}

#if 0
#if defined( FEELPP_INSTANTIATION_MODE )
//
// explicit instances
//
template class ExporterEnsight<Mesh<Simplex<1,1,1> > >;
template class ExporterEnsight<Mesh<Simplex<1,1,2> > >;
template class ExporterEnsight<Mesh<Simplex<2,1,2> > >;
template class ExporterEnsight<Mesh<Simplex<2,2,2> > >;
template class ExporterEnsight<Mesh<Simplex<2,1,3> > >;
template class ExporterEnsight<Mesh<Simplex<3,1,3> > >;

template class ExporterEnsight<Mesh<Simplex<3,2,3> > >;

template class ExporterEnsight<Mesh<Hypercube<1,1,1> > >;
template class ExporterEnsight<Mesh<Hypercube<2,1,2> > >;
template class ExporterEnsight<Mesh<Hypercube<3,1,3> > >;
template class ExporterEnsight<Mesh<Hypercube<3,2,3> > >;

template class ExporterEnsight<Mesh<Simplex<2,3,2> > >;
template class ExporterEnsight<Mesh<Hypercube<2,2> > >;
template class ExporterEnsight<Mesh<Hypercube<2,3> > >;

#endif // FEELPP_INSTANTIATION_MODE
#endif
}
#endif // __EXPORTERENSIGHT_CPP
