/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
ExporterEnsight<MeshType,N>::ExporterEnsight( WorldComm const& worldComm )
:
super( worldComm ),
M_element_type()

{
    init();
}
template<typename MeshType, int N>
ExporterEnsight<MeshType,N>::ExporterEnsight( std::string const& __p, int freq, WorldComm const& worldComm )
    :
    super( "ensight", __p, freq, worldComm ),
    M_element_type()
{
    init();
}
template<typename MeshType, int N>
ExporterEnsight<MeshType,N>::ExporterEnsight( po::variables_map const& vm, std::string const& exp_prefix, WorldComm const& worldComm )
    :
    super( vm, exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
ExporterEnsight<MeshType,N>::ExporterEnsight( std::string const& exp_prefix, WorldComm const& worldComm )
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
ExporterEnsight<MeshType,N>::save() const
{
    if ( !this->worldComm().isActive() ) return;

    //static int freq = 0;

    DVLOG(2) << "[ExporterEnsight::save] checking if frequency is ok\n";


    if ( this->cptOfSave() % this->freq()  )
    {
        this->saveTimeSet();
        return;
    }

    boost::timer ti;
    DVLOG(2) << "[ExporterEnsight::save] export in ensight format\n";

    DVLOG(2) << "[ExporterEnsight::save] export sos\n";
    _F_writeSoSFile();
    DVLOG(2) << "[ExporterEnsight::save] export sos ok, time " << ti.elapsed() << "\n";

    ti.restart();
    DVLOG(2) << "[ExporterEnsight::save] export case file\n";
    _F_writeCaseFile();
    DVLOG(2) << "[ExporterEnsight::save] export case file ok, time " << ti.elapsed() << "\n";

    ti.restart();
    DVLOG(2) << "[ExporterEnsight::save] export geo(mesh) file\n";
    _F_writeGeoFiles();
    DVLOG(2) << "[ExporterEnsight::save] export geo(mesh) file ok, time " << ti.elapsed() << "\n";

    ti.restart();
    DVLOG(2) << "[ExporterEnsight::save] export variable file\n";
    _F_writeVariableFiles();
    DVLOG(2) << "[ExporterEnsight::save] export variable files ok, time " << ti.elapsed() << "\n";

    ti.restart();
    DVLOG(2) << "[ExporterEnsight::save] export time set\n";
    this->saveTimeSet();
    DVLOG(2) << "[ExporterEnsight::save] export time set ok, time " << ti.elapsed() << "\n";
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
        __out << "model: " << __ts->name()
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

            __out << "model: " << __ts->index() << " " << __ts->name()
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

            typename timeset_type::step_type::nodal_scalar_const_iterator __it = ( *__tstp_it )->beginNodalScalar();
            typename timeset_type::step_type::nodal_scalar_const_iterator __end = ( *__tstp_it )->endNodalScalar();

            while ( __it != __end )
            {
                __out << "scalar per node: "
                    << __ts->index() << " " // << *__ts_it->beginStep() << " "
                    << __it->second.name() << " " << __it->first << "-" << this->worldComm().globalSize() << "_" << __it->second.worldComm().localRank() << ".***" << "\n";// important localRank !!
                ++__it;
            }

            typename timeset_type::step_type::nodal_vector_const_iterator __itv = ( *__tstp_it )->beginNodalVector();
            typename timeset_type::step_type::nodal_vector_const_iterator __env = ( *__tstp_it )->endNodalVector();

            while ( __itv != __env )
            {
                __out << "vector per node: "
                    << __ts->index() << " " // << *__ts_it->beginStep() << " "
                    << __itv->second.name() << " " << __itv->first << "-" << this->worldComm().globalSize() << "_" << __itv->second.worldComm().localRank() << ".***" << "\n";// important localRank !!
                ++__itv;
            }

            typename timeset_type::step_type::nodal_tensor2_const_iterator __itt = ( *__tstp_it )->beginNodalTensor2();
            typename timeset_type::step_type::nodal_tensor2_const_iterator __ent = ( *__tstp_it )->endNodalTensor2();

            while ( __itt != __ent )
            {
                __out << "tensor per node: "
                    << __ts->index() << " " // << *__ts_it->beginStep() << " "
                    << __itt->second.name() << " " << __itt->first << "-" << this->worldComm().globalSize() << "_" << __itt->second.worldComm().localRank() << ".***" << "\n"; // important localRank !!
                ++__itt;
            }

            typename timeset_type::step_type::element_scalar_const_iterator __it_el = ( *__tstp_it )->beginElementScalar();
            typename timeset_type::step_type::element_scalar_const_iterator __end_el = ( *__tstp_it )->endElementScalar();

            while ( __it_el != __end_el )
            {
                __out << "scalar per element: "
                    << __ts->index() << " " // << *__ts_it->beginStep() << " "
                    << __it_el->second.name() << " " << __it_el->first << "-" << this->worldComm().globalSize() << "_" << __it_el->second.worldComm().localRank() << ".***" << "\n";// important localRank !!
                ++__it_el;
            }

            typename timeset_type::step_type::element_vector_const_iterator __itv_el = ( *__tstp_it )->beginElementVector();
            typename timeset_type::step_type::element_vector_const_iterator __env_el = ( *__tstp_it )->endElementVector();

            while ( __itv_el != __env_el )
            {
                __out << "vector per element: "
                    << __ts->index() << " " // << *__ts_it->beginStep() << " "
                    << __itv_el->second.name() << " " << __itv_el->first << "-" << this->worldComm().globalSize() << "_" << __itv_el->second.worldComm().localRank() << ".***" << "\n"; // important localRank !!
                ++__itv_el;
            }

            typename timeset_type::step_type::element_tensor2_const_iterator __itt_el = ( *__tstp_it )->beginElementTensor2();
            typename timeset_type::step_type::element_tensor2_const_iterator __ent_el = ( *__tstp_it )->endElementTensor2();

            while ( __itt_el != __ent_el )
            {
                __out << "tensor per element: "
                    << __ts->index() << " " // << *__ts_it->beginStep() << " "
                    << __itt_el->second.name() << " " << __itt_el->first << "-" << this->worldComm().globalSize() << "_" << __itt_el->second.worldComm().localRank() << ".***" << "\n"; // important localRank !!
                ++__itt_el;
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
                       << __ts->name()
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
                << __ts->name()
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
                    saveNodal( __step, __step->beginNodalScalar(), __step->endNodalScalar() );
                    saveNodal( __step, __step->beginNodalVector(), __step->endNodalVector() );
                    saveNodal( __step, __step->beginNodalTensor2(), __step->endNodalTensor2() );

                    saveElement( __step, __step->beginElementScalar(), __step->endElementScalar() );
                    saveElement( __step, __step->beginElementVector(), __step->endElementVector() );
                    saveElement( __step, __step->beginElementTensor2(), __step->endElementTensor2() );

                }

                ++__it;
            }
        }

        ++__ts_it;
    }
}


template<typename MeshType, int N>
template<typename Iterator>
void
ExporterEnsight<MeshType,N>::saveNodal( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const
{
    auto mit = elements(__step->mesh());
    Feel::detail::MeshPoints<float> mp( __step->mesh().get(), this->worldComm(), mit.template get<1>(), mit.template get<2>(), false, true );
    while ( __var != en )
    {
        if ( !__var->second.worldComm().isActive() ) return;

        std::ostringstream __varfname;

        __varfname << this->path() << "/" << __var->first
                   << "-" << this->worldComm().globalSize() << "_" << __var->second.worldComm().localRank() // important localRank
                   << "." << std::setfill( '0' ) << std::setw( 3 ) << __step->index();
        DVLOG(2) << "[ExporterEnsight::saveNodal] saving " << __varfname.str() << "...\n";
        std::fstream __out( __varfname.str().c_str(), std::ios::out | std::ios::binary );

        char buffer[ 80 ];
        strcpy( buffer, __var->second.name().c_str() );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );

        uint16_type nComponents = __var->second.nComponents;

        if ( __var->second.is_vectorial )
            nComponents = 3;

        /**
         * BE CAREFUL HERE some points in the mesh may not be present in the
         * mesh element connectivity, we really need to have an array of the
         * size of the number of points in the mesh even though if some are not
         * in the connectivity and not an array of the dimension of the function
         * space which has the "right" size.
         */
        std::vector<float> m_field( nComponents*mp.ids.size() );
        CHECK( m_field.size()/nComponents == __var->second.localSize()/__var->second.nComponents ) << "Invalid size : " << m_field.size() << "!=" << __var->second.localSize();

        //__field.clear();
        typename mesh_type::element_const_iterator elt_it, elt_en;
        boost::tie( boost::tuples::ignore, elt_it, elt_en ) = elements( *__step->mesh() );
        size_type e = 0;

        if ( !__var->second.areGlobalValuesUpdated() )
            __var->second.updateGlobalValues();

        for ( ; elt_it != elt_en; ++elt_it )
        {
            for ( uint16_type c = 0; c < nComponents; ++c )
            {
                for ( uint16_type p = 0; p < __step->mesh()->numLocalVertices(); ++p, ++e )
                {
                    size_type ptid = mp.old2new[elt_it->point( p ).id()]-1;
                    size_type global_node_id = nComponents * ptid + c ;
#if 0
                    DCHECK( ptid < __step->mesh()->numPoints() ) << "Invalid point id " << ptid << " element: " << elt_it->id()
                                                                 << " local pt:" << p
                                                                 << " mesh numPoints: " << __step->mesh()->numPoints();
                    //DCHECK( global_node_id < __field_size ) << "Invalid dof id : " << global_node_id << " max size : " << __field_size;
#endif
                    if ( c < __var->second.nComponents )
                    {
                        size_type dof_id = boost::get<0>( __var->second.functionSpace()->dof()->localToGlobal( elt_it->id(),p, c ) );

                        m_field[global_node_id] = __var->second.globalValue( dof_id );
                    }
                    else
                        m_field[global_node_id] = 0;
                }
            }
        }

        //std::vector<float> field;
        //std::for_each( m_field.begin(), m_field.end(), [&field]( std::pair<size_type, float> const& p ) { field.push_back( p.second ); });
        __out.write( ( char * ) m_field.data(), m_field.size() * sizeof( float ) );

        DVLOG(2) << "[ExporterEnsight::saveNodal] saving " << __varfname.str() << "done\n";
        ++__var;
    }
}
template<typename MeshType, int N>
template<typename Iterator>
void
ExporterEnsight<MeshType,N>::saveElement( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const
{
    while ( __evar != __evaren )
    {
        if ( !__evar->second.worldComm().isActive() ) return;

        std::ostringstream __evarfname;

        __evarfname << this->path() << "/" << __evar->first
                    << "-" << this->worldComm().globalSize() << "_" << __evar->second.worldComm().localRank() // important localRank
                    << "." << std::setfill( '0' ) << std::setw( 3 ) << __step->index();
        DVLOG(2) << "[ExporterEnsight::saveElement] saving " << __evarfname.str() << "...\n";
        std::fstream __out( __evarfname.str().c_str(), std::ios::out | std::ios::binary );

        char buffer[ 80 ];
        strcpy( buffer, __evar->second.name().c_str() );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );

        typename mesh_type::parts_const_iterator_type p_it = __step->mesh()->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = __step->mesh()->endParts();

        for ( ; p_it != p_en; ++p_it )
        {
            sprintf( buffer, "part %d",p_it->first );
            __out.write( ( char * ) & buffer, sizeof( buffer ) );
            DVLOG(2) << "part " << buffer << "\n";
            strcpy( buffer, this->elementType().c_str() );
            __out.write( ( char * ) & buffer, sizeof( buffer ) );
            DVLOG(2) << "element type " << buffer << "\n";

            uint16_type nComponents = __evar->second.nComponents;

            if ( __evar->second.is_vectorial )
                nComponents = 3;

            size_type __field_size = nComponents*__evar->second.size()/__evar->second.nComponents;
            ublas::vector<float> __field( __field_size );
            __field.clear();
            typename mesh_type::marker_element_const_iterator elt_it;
            typename mesh_type::marker_element_const_iterator elt_en;
            boost::tie( elt_it, elt_en ) = __step->mesh()->elementsWithMarker( p_it->first,
                                                                               __evar->second.worldComm().localRank() ); // important localRank!!!!

            if ( !__evar->second.areGlobalValuesUpdated() )
                __evar->second.updateGlobalValues();

            DVLOG(2) << "[saveElement] firstLocalIndex = " << __evar->second.firstLocalIndex() << "\n";
            DVLOG(2) << "[saveElement] lastLocalIndex = " << __evar->second.lastLocalIndex() << "\n";
            DVLOG(2) << "[saveElement] field.size = " << __field_size << "\n";
            size_type e = 0;

            for ( ; elt_it != elt_en; ++elt_it, ++e )
            {
                DVLOG(2) << "pid : " << this->worldComm().globalRank()
                              << " elt_it :  " << elt_it->id()
                              << " e : " << e << "\n";

                for ( int c = 0; c < nComponents; ++c )
                {
                    size_type global_node_id = nComponents * e + c ;

                    if ( c < __evar->second.nComponents )
                    {
                        size_type dof_id = boost::get<0>( __evar->second.functionSpace()->dof()->localToGlobal( elt_it->id(),0, c ) );

                        DVLOG(2) << "c : " << c
                                      << " gdofid: " << global_node_id
                                      << " dofid : " << dof_id
                                      << " f.size : " <<  __field.size()
                                      << " e.size : " <<  __evar->second.size()
                                      << "\n";

                        __field[global_node_id] = __evar->second.globalValue( dof_id );
#if 0

                        if ( dof_id >= __evar->second.firstLocalIndex() &&
                                dof_id < __evar->second.lastLocalIndex()  )
                            __field[global_node_id] = __evar->second( dof_id );

                        else
                            __field[global_node_id] = 0;

#endif

#if 1
                        //__field[global_node_id] = __evar->second.globalValue(dof_id);
                        DVLOG(2) << "c : " << c
                                      << " gdofid: " << global_node_id
                                      << " dofid : " << dof_id
                                      << " field :  " << __field[global_node_id]
                                      << " evar: " << __evar->second.globalValue( dof_id ) << "\n";
#endif
                    }

                    else
                        __field[global_node_id] = 0;
                }
            }

            __out.write( ( char * ) __field.data().begin(), nComponents * e * sizeof( float ) );
        }

        DVLOG(2) << "[ExporterEnsight::saveElement] saving " << __evarfname.str() << "done\n";
        ++__evar;
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
    Feel::detail::MeshPoints<float> mp( __mesh, this->worldComm(), mit.template get<1>(), mit.template get<2>(), false, true );
    size_type __nv = mp.ids.size();
    __out.write( ( char * ) &__nv, sizeof( int ) );
    LOG(INFO) << "n pts = " << __nv << " numppoints=" << __mesh->numPoints();
    __out.write( ( char * ) & mp.ids.front(), mp.ids.size() * sizeof( int ) );
    __out.write( ( char * ) mp.coords.data(), mp.coords.size() * sizeof( float ) );

    typename mesh_type::parts_const_iterator_type p_it = __mesh->beginParts();
    typename mesh_type::parts_const_iterator_type p_en = __mesh->endParts();

    for ( ; p_it != p_en; ++p_it )
    {
        sprintf( buffer, "part %d",p_it->first );
        //    strcpy( buffer, "part 1" );

        __out.write( ( char * ) & buffer, sizeof( buffer ) );
        sprintf( buffer, "Material %d",p_it->first );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );

        strcpy( buffer, this->elementType().c_str() );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );

        //    typename mesh_type::element_const_iterator elt_it = __mesh->beginElement();
        //    typename mesh_type::element_const_iterator elt_en = __mesh->endElement();
        typename mesh_type::marker_element_const_iterator elt_it;// = __mesh->beginElementWithMarker(p_it->first);
        typename mesh_type::marker_element_const_iterator elt_en;// = __mesh->endElementWithMarker(p_it->first);
        boost::tie( elt_it, elt_en ) = __mesh->elementsWithMarker( p_it->first,
                                                                   __mesh->worldComm().localRank() ); // important localRank!!!!

        //	int __ne = __mesh->numElements();
        //int __ne = p_it->second;
        int __ne = std::distance( elt_it, elt_en );

        DVLOG(2) << "num Elements to save : " << __ne << "\n";

        __out.write( ( char * ) &__ne, sizeof( int ) );

        idelem.resize( __ne );


        for ( size_type e = 0; elt_it != elt_en; ++elt_it, ++e )
        {
            idelem[e] = elt_it->id() + 1;
        }

        __out.write( ( char * ) & idelem.front(), idelem.size() * sizeof( int ) );

        //	elt_it = __mesh->beginElement();
        boost::tie( elt_it, elt_en ) = __mesh->elementsWithMarker( p_it->first,
                                                                   __mesh->worldComm().localRank() ); // important localRank!!!!
        //elt_it = __mesh->beginElementWithMarker(p_it->first);
        std::vector<int> eids( __mesh->numLocalVertices()*__ne );
        size_type e= 0;
        for ( ; elt_it != elt_en; ++elt_it, ++e )
        {
            for ( size_type j = 0; j < __mesh->numLocalVertices(); j++ )
            {
                // ensight id start at 1
                int __id = mp.old2new[elt_it->point( j ).id()];
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
