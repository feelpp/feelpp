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
#ifndef __EXPORTERENSIGHTGOLD_CPP
#define __EXPORTERENSIGHTGOLD_CPP

#include <feel/feelcore/feel.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/timeset.hpp>
#include <feel/feelfilters/exporterensight.hpp>


namespace Feel
{
template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( WorldComm const& worldComm )
:
super( worldComm ),
_M_element_type()

{
    init();
}
template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( std::string const& __p, int freq, WorldComm const& worldComm )
    :
    super( "ensightgold", __p, freq, worldComm ),
    _M_element_type()
{
    init();
}
template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( po::variables_map const& vm, std::string const& exp_prefix, WorldComm const& worldComm )
    :
    super( vm, exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( ExporterEnsightGold const & __ex )
    :
    super( __ex ),
    _M_element_type( __ex._M_element_type )
{
}

template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::~ExporterEnsightGold()
{}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::init()
{
    if ( mesh_type::nDim == 1 )
        if ( mesh_type::Shape == SHAPE_LINE )
            _M_element_type = ( mesh_type::nOrder == 1 )?"bar2":"bar3";

    if ( mesh_type::nDim == 2 )
    {
        if ( mesh_type::Shape == SHAPE_TRIANGLE )
            _M_element_type = ( mesh_type::nOrder == 1 )?"tria3":"tria6";

        else if ( mesh_type::Shape == SHAPE_QUAD )
            _M_element_type = ( mesh_type::nOrder == 1 )?"quad4":"quad8";
    }

    if ( mesh_type::nDim == 3 )
    {
        if ( mesh_type::Shape == SHAPE_TETRA )
            _M_element_type = ( mesh_type::nOrder == 1 )?"tetra4":"tetra10";

        else if ( mesh_type::Shape == SHAPE_HEXA )
            _M_element_type = ( mesh_type::nOrder == 1 )?"hexa8":"hexa20";
    }
}
template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::save() const
{
    if ( !this->worldComm().isActive() ) return;

    //static int freq = 0;

    DVLOG(2) << "checking if frequency is ok\n";


    if ( this->cptOfSave() % this->freq()  )
    {
        this->saveTimeSet();
        return;
    }

    boost::timer ti;
    DVLOG(2) << "export in ensight format\n";

    DVLOG(2) << "export sos\n";
    _F_writeSoSFile();
    DVLOG(2) << "export sos ok, time " << ti.elapsed() << "\n";

    ti.restart();
    DVLOG(2) << "export case file\n";
    _F_writeCaseFile();
    DVLOG(2) << "export case file ok, time " << ti.elapsed() << "\n";

    ti.restart();
    DVLOG(2) << "export geo(mesh) file\n";
    _F_writeGeoFiles();
    DVLOG(2) << "export geo(mesh) file ok, time " << ti.elapsed() << "\n";

    ti.restart();
    DVLOG(2) << "export variable file\n";
    _F_writeVariableFiles();
    DVLOG(2) << "export variable files ok, time " << ti.elapsed() << "\n";

    ti.restart();
    DVLOG(2) << "export time set\n";
    this->saveTimeSet();
    DVLOG(2) << "export time set ok, time " << ti.elapsed() << "\n";
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::_F_writeSoSFile() const
{
    // only on proc 0
    if ( this->worldComm().rank() == this->worldComm().masterRank() )
    {
        std::ostringstream filestr;
        filestr << this->path() << "/" << this->prefix() << "-" << this->worldComm().globalSize() << ".sos";
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
    }
}
template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::_F_writeCaseFile() const
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
          << "type: ensight gold\n"
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

            if ( this->useSingleTransientFile() )
                __out << "model: " << __ts->index() << " 1 " << __ts->name()
                      << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank() << ".geo";
            else
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

        typename timeset_type::step_type::nodal_scalar_const_iterator __it = ( *__ts->rbeginStep() )->beginNodalScalar();
        typename timeset_type::step_type::nodal_scalar_const_iterator __end = ( *__ts->rbeginStep() )->endNodalScalar();

        while ( __it != __end )
        {
            if ( this->useSingleTransientFile() )
                __out << "scalar per node: "
                      << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                      << __it->second.name() << " " << __it->first << "-" << this->worldComm().globalSize() << "_" << __it->second.worldComm().localRank() << "\n";// important localRank !!
            else
                __out << "scalar per node: "
                      << __ts->index() << " " // << *__ts_it->beginStep() << " "
                      << __it->second.name() << " " << __it->first << "-" << this->worldComm().globalSize() << "_" << __it->second.worldComm().localRank() << ".***" << "\n";// important localRank !!

            ++__it;
        }

        typename timeset_type::step_type::nodal_vector_const_iterator __itv = ( *__ts->rbeginStep() )->beginNodalVector();
        typename timeset_type::step_type::nodal_vector_const_iterator __env = ( *__ts->rbeginStep() )->endNodalVector();

        while ( __itv != __env )
        {
            if ( this->useSingleTransientFile() )
                __out << "vector per node: "
                      << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                      << __itv->second.name() << " " << __itv->first << "-" << this->worldComm().globalSize() << "_" << __itv->second.worldComm().localRank() << "\n";// important localRank !!
            else
                __out << "vector per node: "
                      << __ts->index() << " " // << *__ts_it->beginStep() << " "
                      << __itv->second.name() << " " << __itv->first << "-" << this->worldComm().globalSize() << "_" << __itv->second.worldComm().localRank() << ".***" << "\n";// important localRank !!
            ++__itv;
        }

        typename timeset_type::step_type::nodal_tensor2_const_iterator __itt = ( *__ts->rbeginStep() )->beginNodalTensor2();
        typename timeset_type::step_type::nodal_tensor2_const_iterator __ent = ( *__ts->rbeginStep() )->endNodalTensor2();

        while ( __itt != __ent )
        {
            if ( this->useSingleTransientFile() )
                __out << "tensor per node: "
                      << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                      << __itt->second.name() << " " << __itt->first << "-" << this->worldComm().globalSize() << "_" << __itt->second.worldComm().localRank() << "\n"; // important localRank !!
            else
                __out << "tensor per node: "
                      << __ts->index() << " " // << *__ts_it->beginStep() << " "
                      << __itt->second.name() << " " << __itt->first << "-" << this->worldComm().globalSize() << "_" << __itt->second.worldComm().localRank() << ".***" << "\n"; // important localRank !!
            ++__itt;
        }

        typename timeset_type::step_type::element_scalar_const_iterator __it_el = ( *__ts->rbeginStep() )->beginElementScalar();
        typename timeset_type::step_type::element_scalar_const_iterator __end_el = ( *__ts->rbeginStep() )->endElementScalar();

        while ( __it_el != __end_el )
        {
            if ( this->useSingleTransientFile() )
                __out << "scalar per element: "
                      << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                      << __it_el->second.name() << " " << __it_el->first << "-" << this->worldComm().globalSize() << "_" << __it_el->second.worldComm().localRank() << "\n";// important localRank !!
            else
                __out << "scalar per element: "
                      << __ts->index() << " " // << *__ts_it->beginStep() << " "
                      << __it_el->second.name() << " " << __it_el->first << "-" << this->worldComm().globalSize() << "_" << __it_el->second.worldComm().localRank() << ".***" << "\n";// important localRank !!

            ++__it_el;
        }

        typename timeset_type::step_type::element_vector_const_iterator __itv_el = ( *__ts->rbeginStep() )->beginElementVector();
        typename timeset_type::step_type::element_vector_const_iterator __env_el = ( *__ts->rbeginStep() )->endElementVector();

        while ( __itv_el != __env_el )
        {
            if ( this->useSingleTransientFile() )
                __out << "vector per element: "
                      << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                      << __itv_el->second.name() << " " << __itv_el->first << "-" << this->worldComm().globalSize() << "_" << __itv_el->second.worldComm().localRank() << "\n"; // important localRank !!
            else
                __out << "vector per element: "
                      << __ts->index() << " " // << *__ts_it->beginStep() << " "
                      << __itv_el->second.name() << " " << __itv_el->first << "-" << this->worldComm().globalSize() << "_" << __itv_el->second.worldComm().localRank() << ".***" << "\n"; // important localRank !!
            ++__itv_el;
        }

        typename timeset_type::step_type::element_tensor2_const_iterator __itt_el = ( *__ts->rbeginStep() )->beginElementTensor2();
        typename timeset_type::step_type::element_tensor2_const_iterator __ent_el = ( *__ts->rbeginStep() )->endElementTensor2();

        while ( __itt_el != __ent_el )
        {
            __out << "tensor per element: "
                  << __ts->index() << " " // << *__ts_it->beginStep() << " "
                  << __itt_el->second.name() << " " << __itt_el->first << "-" << this->worldComm().globalSize() << "_" << __itt_el->second.worldComm().localRank() << ".***" << "\n"; // important localRank !!
            ++__itt_el;
        }

        ++__ts_it;
    }

    __out << "TIME:\n";
    __ts_it = this->beginTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;
        typename timeset_type::step_const_iterator __its = __ts->beginStep();

        __out << "time set:        " << __ts->index() << "\n"
              << "number of steps: " << __ts->numberOfSteps() << "\n"
              << "filename start number: " << ( *__its )->index() << "\n"
              << "filename increment: " << 1 << "\n"
              << "time values: ";

        uint16_type __l = 0;
        typename timeset_type::step_const_iterator __ens = __ts->endStep();

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

    if ( this->useSingleTransientFile() )
    {
        __out << "FILE\n";
        __out << "file set : 1\n";
        auto ts = *this->beginTimeSet();
        __out << "number of steps: " << ts->numberOfSteps() << "\n";
    }

    __out << "\n";

    LOG(INFO) << "rank : " << this->worldComm().globalRank();
    LOG(INFO) << "gg"  << (( Environment::numberOfProcessors() > 1 )  && ( this->worldComm().globalRank() == 0 )) << "\n";
    if ( ( Environment::numberOfProcessors() > 1 )  && ( this->worldComm().globalRank() == 0 ) )
    {
        __out << "APPENDED_CASEFILES\n"
              << "total number of cfiles: " << Environment::numberOfProcessors()-1 << "\n"
              << "cfiles global path: " << fs::current_path().string() << "\n"
              << "cfiles: ";
        for(int p = 1; p < Environment::numberOfProcessors(); ++p )
        {
            std::ostringstream filestr;
            filestr << this->prefix() << "-"
                    << this->worldComm().globalSize() << "_" << p << ".case";
            __out << filestr.str() << "\n        ";
        }
    }
    __out.close();

}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::_F_writeGeoFiles() const
{
    namespace lambda = boost::lambda;

    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;

        typename timeset_type::step_const_iterator __it = __ts->beginStep();
        typename timeset_type::step_const_iterator __end = __ts->endStep();
        __it = boost::prior( __end );

        if ( this->exporterGeometry() == EXPORTER_GEOMETRY_STATIC )
        {
            std::ostringstream __geofname;
            __geofname << this->path() << "/"
                       << __ts->name()
                       << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank()
                       << ".geo";
            _M_filename =  __geofname.str();
            CHECK( (*__it)->hasMesh() || __ts->hasMesh()  ) << "Invalid mesh data structure in static geometry mode\n";
            if ( __ts->hasMesh() )
                __ts->mesh()->accept( const_cast<ExporterEnsightGold<MeshType,N>&>( *this ) );
            if ( (*__it)->hasMesh() && !__ts->hasMesh() )
                (*__it)->mesh()->accept( const_cast<ExporterEnsightGold<MeshType,N>&>( *this ) );
        }

        while ( __it != __end )
        {
            typename timeset_type::step_ptrtype __step = *__it;


            std::ostringstream __geofname;

            if ( this->exporterGeometry() != EXPORTER_GEOMETRY_STATIC )
            {
                __geofname << this->path() << "/"
                           << __ts->name()
                           << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank()
                           << ".geo";
                if ( !this->useSingleTransientFile() )
                    __geofname << std::setfill( '0' ) << std::setw( 3 ) << __step->index();

                if ( __step->isInMemory() )
                {
                    //__writegeo( __step->mesh(), __ts->name(), __geofname.str() );
                    //, __ts->name(), __geofname.str() );
                    _M_filename =  __geofname.str();
                    __step->mesh()->accept( const_cast<ExporterEnsightGold<MeshType,N>&>( *this ) );
                }
            }
            ++__it;
        }

        ++__ts_it;
    }
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::_F_writeVariableFiles() const
{
    namespace lambda = boost::lambda;


    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;

        typename timeset_type::step_const_iterator __it = __ts->beginStep();
        typename timeset_type::step_const_iterator __end = __ts->endStep();
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

        ++__ts_it;
    }
}


template<typename MeshType, int N>
template<typename Iterator>
void
ExporterEnsightGold<MeshType,N>::saveNodal( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const
{
    while ( __var != en )
    {
        if ( !__var->second.worldComm().isActive() ) return;

        std::ostringstream __varfname;

        __varfname << this->path() << "/" << __var->first
                   << "-" << this->worldComm().globalSize() << "_" << __var->second.worldComm().localRank(); // important localRank
        if ( !this->useSingleTransientFile() )
            __varfname << "." << std::setfill( '0' ) << std::setw( 3 ) << __step->index();
        DVLOG(2) << "[ExporterEnsightGold::saveNodal] saving " << __varfname.str() << "...\n";
        std::fstream __out;
        if ( this->useSingleTransientFile() )
            __out.open( __varfname.str().c_str(), std::ios::out | std::ios::app | std::ios::binary );
        else
            __out.open( __varfname.str().c_str(), std::ios::out | std::ios::binary );

        char buffer[ 80 ];

        if ( this->useSingleTransientFile() )
        {
            strcpy(buffer,"BEGIN TIME STEP");
            __out.write((char*)&buffer,sizeof(buffer));
        }

        strcpy( buffer, __var->second.name().c_str() );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );

        typename mesh_type::parts_const_iterator_type p_it = __step->mesh()->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = __step->mesh()->endParts();

        for ( ; p_it != p_en; ++p_it )
        {
            strcpy( buffer, "part" );
            __out.write( ( char * ) & buffer, sizeof( buffer ) );
            int partid = p_it->first;
            __out.write( ( char * ) & partid, sizeof(int) );
            DVLOG(2) << "part " << buffer << "\n";
            strcpy( buffer, "coordinates" );
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
            size_type __field_size = __step->mesh()->numPoints();
            if ( __var->second.is_vectorial )
                __field_size *= 3;
            ublas::vector<float> __field( __field_size );

            __field.clear();
            //typename mesh_type::element_const_iterator elt_it, elt_en;
            //boost::tie( boost::tuples::ignore, elt_it, elt_en ) = elements( *__step->mesh() );

            typename mesh_type::marker_element_const_iterator elt_it;
            typename mesh_type::marker_element_const_iterator elt_en;
            boost::tie( elt_it, elt_en ) = __step->mesh()->elementsWithMarker( p_it->first,
                                                                               __var->second.worldComm().localRank() ); // important localRank!!!!


            size_type e = 0;

            if ( !__var->second.areGlobalValuesUpdated() )
                __var->second.updateGlobalValues();

            for ( ; elt_it != elt_en; ++elt_it )
            {
                for ( uint16_type c = 0; c < nComponents; ++c )
                {
                    for ( uint16_type p = 0; p < __step->mesh()->numLocalVertices(); ++p, ++e )
                    {
                        size_type ptid = elt_it->point( p ).id();
                        size_type global_node_id = nComponents * ptid + c ;
                        DCHECK( ptid < __step->mesh()->numPoints() ) << "Invalid point id " << ptid << " element: " << elt_it->id()
                                                                     << " local pt:" << p
                                                                     << " mesh numPoints: " << __step->mesh()->numPoints();
                        DCHECK( global_node_id < __field_size ) << "Invalid dof id : " << global_node_id << " max size : " << __field_size;

                        if ( c < __var->second.nComponents )
                        {
                            size_type dof_id = boost::get<0>( __var->second.functionSpace()->dof()->localToGlobal( elt_it->id(),p, c ) );

#if 0

                            if ( dof_id >= __var->second.firstLocalIndex() &&
                                 dof_id < __var->second.lastLocalIndex()  )
                                __field[global_node_id] = __var->second( dof_id );

                            else
                                __field[global_node_id] = 0;

#else
                            __field[global_node_id] = __var->second.globalValue( dof_id );
#endif
                        }

                        else
                            __field[global_node_id] = 0;
                    }
                }
            }

            __out.write( ( char * ) __field.data().begin(), __field.size() * sizeof( float ) );

        } // parts loop

        if ( this->useSingleTransientFile() )
        {
            strcpy(buffer,"END TIME STEP");
            __out.write((char*)&buffer,sizeof(buffer));
        }
        DVLOG(2) << "[ExporterEnsightGold::saveNodal] saving " << __varfname.str() << "done\n";
        ++__var;
    }
}
template<typename MeshType, int N>
template<typename Iterator>
void
ExporterEnsightGold<MeshType,N>::saveElement( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const
{
    while ( __evar != __evaren )
    {
        if ( !__evar->second.worldComm().isActive() ) return;

        std::ostringstream __evarfname;

        __evarfname << this->path() << "/" << __evar->first
                    << "-" << this->worldComm().globalSize() << "_" << __evar->second.worldComm().localRank() // important localRank
                    << "." << std::setfill( '0' ) << std::setw( 3 ) << __step->index();
        DVLOG(2) << "[ExporterEnsightGold::saveElement] saving " << __evarfname.str() << "...\n";
        std::fstream __out( __evarfname.str().c_str(), std::ios::out | std::ios::binary );

        char buffer[ 80 ];

        if ( this->useSingleTransientFile() )
        {
            strcpy(buffer,"BEGIN TIME STEP");
            __out.write((char*)&buffer,sizeof(buffer));
        }

        strcpy( buffer, __evar->second.name().c_str() );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );

        typename mesh_type::parts_const_iterator_type p_it = __step->mesh()->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = __step->mesh()->endParts();

        for ( ; p_it != p_en; ++p_it )
        {
            strcpy( buffer, "part" );
            __out.write( ( char * ) & buffer, sizeof( buffer ) );
            //sprintf( buffer, "%d",p_it->first );
            int partid = p_it->first;
            __out.write( ( char * ) & partid, sizeof(int) );
            DVLOG(2) << "part " << buffer << "\n";
            //strcpy( buffer, this->elementType().c_str() );
            strcpy( buffer, "coordinates" );
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
       if ( this->useSingleTransientFile() )
       {
            strcpy(buffer,"END TIME STEP");
            __out.write((char*)&buffer,sizeof(buffer));
       }


        DVLOG(2) << "[ExporterEnsightGold::saveElement] saving " << __evarfname.str() << "done\n";
        ++__evar;
    }
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::visit( mesh_type* __mesh )
{
    char buffer[ 80 ];
    std::vector<int> idnode, idelem;

    std::fstream __out;
    if ( this->useSingleTransientFile() )
        __out.open( _M_filename.c_str(), std::ios::out |  std::ios::app | std::ios::binary );
    else
        __out.open( _M_filename.c_str(), std::ios::out | std::ios::binary );


    strcpy( buffer, "C Binary" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );

    if ( this->useSingleTransientFile() )
    {
        strcpy(buffer,"BEGIN TIME STEP");
        __out.write((char*)&buffer,sizeof(buffer));
    }




    strcpy( buffer, _M_filename.c_str() );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "elements" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "node id given" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "element id given" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );

    typename mesh_type::parts_const_iterator_type p_it = __mesh->beginParts();
    typename mesh_type::parts_const_iterator_type p_en = __mesh->endParts();

    for ( ; p_it != p_en; ++p_it )
    {
        strcpy( buffer, "part" );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );
        int partid = p_it->first;
        __out.write( ( char * ) & partid, sizeof(int) );

        sprintf( buffer, "Material %d",p_it->first );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );


        strcpy( buffer, "coordinates" );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );



        size_type __nv = __mesh->numPoints();
        __out.write( ( char * ) &__nv, sizeof( int ) );


        idnode.resize( __nv );

        size_type __coords_size = 3*__nv;
        ublas::vector<float> __coords( __coords_size );
        __coords.clear();
        typename mesh_type::point_const_iterator pt_it = __mesh->beginPoint();
        typename mesh_type::point_const_iterator pt_en = __mesh->endPoint();

        for ( size_type __count = 0; pt_it != pt_en; ++pt_it, ++__count )
        {
            size_type __id = pt_it->id();
            idnode[__count] = __id+1;
            __coords[__count] = ( float ) pt_it->node()[0];

        if ( mesh_type::nRealDim >= 2 )
            __coords[__nv+__count] = ( float ) pt_it->node()[1];

        else
            __coords[__nv+__count] = float( 0 );

        if ( mesh_type::nRealDim >= 3 )
            __coords[2*__nv+__count] = float( pt_it->node()[2] );

        else
            __coords[2*__nv+__count] = float( 0 );
        }

        __out.write( ( char * ) & idnode.front(), idnode.size() * sizeof( int ) );
        __out.write( ( char * ) __coords.data().begin(), __coords.size() * sizeof( float ) );

        // elements
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

        for ( ; elt_it != elt_en; ++elt_it )
        {
            for ( size_type j = 0; j < __mesh->numLocalVertices(); j++ )
            {
                // ensight id start at 1
                int __id = elt_it->point( j ).id()+1;
                __out.write( ( char * ) & __id , sizeof( int ) );
            }
        }
    }
    if ( this->useSingleTransientFile() )
    {
        strcpy(buffer,"END TIME STEP");
        __out.write((char*)&buffer,sizeof(buffer));
    }
}

#if 0
#if defined( FEELPP_INSTANTIATION_MODE )
//
// explicit instances
//
template class ExporterEnsightGold<Mesh<Simplex<1,1,1> > >;
template class ExporterEnsightGold<Mesh<Simplex<1,1,2> > >;
template class ExporterEnsightGold<Mesh<Simplex<2,1,2> > >;
template class ExporterEnsightGold<Mesh<Simplex<2,2,2> > >;
template class ExporterEnsightGold<Mesh<Simplex<2,1,3> > >;
template class ExporterEnsightGold<Mesh<Simplex<3,1,3> > >;

template class ExporterEnsightGold<Mesh<Simplex<3,2,3> > >;

template class ExporterEnsightGold<Mesh<Hypercube<1,1,1> > >;
template class ExporterEnsightGold<Mesh<Hypercube<2,1,2> > >;
template class ExporterEnsightGold<Mesh<Hypercube<3,1,3> > >;
template class ExporterEnsightGold<Mesh<Hypercube<3,2,3> > >;

template class ExporterEnsightGold<Mesh<Simplex<2,3,2> > >;
template class ExporterEnsightGold<Mesh<Hypercube<2,2> > >;
template class ExporterEnsightGold<Mesh<Hypercube<2,3> > >;

#endif // FEELPP_INSTANTIATION_MODE
#endif
}
#endif // __EXPORTERENSIGHT_CPP
