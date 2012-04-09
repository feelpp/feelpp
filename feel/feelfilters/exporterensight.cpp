/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
ExporterEnsight<MeshType,N>::ExporterEnsight( std::string const& __p, int freq, WorldComm const& worldComm )
    :
    super( "ensight", __p, freq, worldComm ),
    _M_element_type()
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
ExporterEnsight<MeshType,N>::ExporterEnsight( ExporterEnsight const & __ex )
    :
    super( __ex ),
    _M_element_type( __ex._M_element_type )
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
ExporterEnsight<MeshType,N>::save() const
{
    if ( !this->worldComm().isActive() ) return;

    //static int freq = 0;

    Debug( 8006 ) << "[ExporterEnsight::save] checking if frequency is ok\n";


    if ( this->cptOfSave() % this->freq()  )
    {
        this->saveTimeSet();
        return;
    }

    boost::timer ti;
    Debug( 8006 ) << "[ExporterEnsight::save] export in ensight format\n";

    Debug( 8006 ) << "[ExporterEnsight::save] export sos\n";
    _F_writeSoSFile();
    Debug( 8006 ) << "[ExporterEnsight::save] export sos ok, time " << ti.elapsed() << "\n";

    ti.restart();
    Debug( 8006 ) << "[ExporterEnsight::save] export case file\n";
    _F_writeCaseFile();
    Debug( 8006 ) << "[ExporterEnsight::save] export case file ok, time " << ti.elapsed() << "\n";

    ti.restart();
    Debug( 8006 ) << "[ExporterEnsight::save] export geo(mesh) file\n";
    _F_writeGeoFiles();
    Debug( 8006 ) << "[ExporterEnsight::save] export geo(mesh) file ok, time " << ti.elapsed() << "\n";

    ti.restart();
    Debug( 8006 ) << "[ExporterEnsight::save] export variable file\n";
    _F_writeVariableFiles();
    Debug( 8006 ) << "[ExporterEnsight::save] export variable files ok, time " << ti.elapsed() << "\n";

    ti.restart();
    Debug( 8006 ) << "[ExporterEnsight::save] export time set\n";
    this->saveTimeSet();
    Debug( 8006 ) << "[ExporterEnsight::save] export time set ok, time " << ti.elapsed() << "\n";
}

template<typename MeshType, int N>
void
ExporterEnsight<MeshType,N>::_F_writeSoSFile() const
{
    // only on proc 0
    if ( this->worldComm().rank() == this->worldComm().masterRank() )
    {
        std::ostringstream filestr;
        filestr << this->path() << "/" << this->prefix() << "-" << this->worldComm().globalSize() << ".sos";
        std::ofstream __out( filestr.str().c_str() );

        if ( __out.fail() )
        {
            Debug( 3100 ) << "cannot open " << filestr.str()  << "\n";
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
ExporterEnsight<MeshType,N>::_F_writeCaseFile() const
{
    std::ostringstream filestr;
    filestr << this->path() << "/"
            << this->prefix() << "-"
            << this->worldComm().globalSize() << "_" << this->worldComm().globalRank() << ".case";
    std::ofstream __out( filestr.str().c_str() );

    if ( __out.fail() )
    {
        Debug( 3100 ) << "cannot open " << filestr.str()  << "\n";
        exit( 0 );
    }

    __out << "FORMAT:\n"
          << "type: ensight \n"
          << "GEOMETRY:\n";

    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;
        __out << "model: " << __ts->index() << " " << __ts->name()
              << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank() << ".geo***"  << "\n";
        ++__ts_it;
    }

    __out << "VARIABLES:" << "\n";

    __ts_it = this->beginTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;

        typename timeset_type::step_type::nodal_scalar_const_iterator __it = ( *__ts->rbeginStep() )->beginNodalScalar();
        typename timeset_type::step_type::nodal_scalar_const_iterator __end = ( *__ts->rbeginStep() )->endNodalScalar();

        while ( __it != __end )
        {
            __out << "scalar per node: "
                  << __ts->index() << " " // << *__ts_it->beginStep() << " "
                  << __it->second.name() << " " << __it->first << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank() << ".***" << "\n";// important localRank !!
            ++__it;
        }

        typename timeset_type::step_type::nodal_vector_const_iterator __itv = ( *__ts->rbeginStep() )->beginNodalVector();
        typename timeset_type::step_type::nodal_vector_const_iterator __env = ( *__ts->rbeginStep() )->endNodalVector();

        while ( __itv != __env )
        {
            __out << "vector per node: "
                  << __ts->index() << " " // << *__ts_it->beginStep() << " "
                  << __itv->second.name() << " " << __itv->first << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank() << ".***" << "\n";// important localRank !!
            ++__itv;
        }

        typename timeset_type::step_type::nodal_tensor2_const_iterator __itt = ( *__ts->rbeginStep() )->beginNodalTensor2();
        typename timeset_type::step_type::nodal_tensor2_const_iterator __ent = ( *__ts->rbeginStep() )->endNodalTensor2();

        while ( __itt != __ent )
        {
            __out << "tensor per node: "
                  << __ts->index() << " " // << *__ts_it->beginStep() << " "
                  << __itt->second.name() << " " << __itt->first << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank() << ".***" << "\n"; // important localRank !!
            ++__itt;
        }

        typename timeset_type::step_type::element_scalar_const_iterator __it_el = ( *__ts->rbeginStep() )->beginElementScalar();
        typename timeset_type::step_type::element_scalar_const_iterator __end_el = ( *__ts->rbeginStep() )->endElementScalar();

        while ( __it_el != __end_el )
        {
            __out << "scalar per element: "
                  << __ts->index() << " " // << *__ts_it->beginStep() << " "
                  << __it_el->second.name() << " " << __it_el->first << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank() << ".***" << "\n";// important localRank !!
            ++__it_el;
        }

        typename timeset_type::step_type::element_vector_const_iterator __itv_el = ( *__ts->rbeginStep() )->beginElementVector();
        typename timeset_type::step_type::element_vector_const_iterator __env_el = ( *__ts->rbeginStep() )->endElementVector();

        while ( __itv_el != __env_el )
        {
            __out << "vector per element: "
                  << __ts->index() << " " // << *__ts_it->beginStep() << " "
                  << __itv_el->second.name() << " " << __itv_el->first << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank() << ".***" << "\n"; // important localRank !!
            ++__itv_el;
        }

        typename timeset_type::step_type::element_tensor2_const_iterator __itt_el = ( *__ts->rbeginStep() )->beginElementTensor2();
        typename timeset_type::step_type::element_tensor2_const_iterator __ent_el = ( *__ts->rbeginStep() )->endElementTensor2();

        while ( __itt_el != __ent_el )
        {
            __out << "tensor per element: "
                  << __ts->index() << " " // << *__ts_it->beginStep() << " "
                  << __itt_el->second.name() << " " << __itt_el->first << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank() << ".***" << "\n"; // important localRank !!
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
        __it = boost::prior( __end );

        while ( __it != __end )
        {
            typename timeset_type::step_ptrtype __step = *__it;;


            std::ostringstream __geofname;

            __geofname << this->path() << "/"
                       << __ts->name()
                       << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank()
                       << ".geo" << std::setfill( '0' ) << std::setw( 3 ) << __step->index();

            if ( __step->isInMemory() )
            {
                //__writegeo( __step->mesh(), __ts->name(), __geofname.str() );
                //, __ts->name(), __geofname.str() );
                _M_filename =  __geofname.str();
                __step->mesh()->accept( const_cast<ExporterEnsight<MeshType,N>&>( *this ) );
            }

            ++__it;
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
ExporterEnsight<MeshType,N>::saveNodal( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const
{
    while ( __var != en )
    {
        if ( !__var->second.worldComm().isActive() ) return;

        std::ostringstream __varfname;

        __varfname << this->path() << "/" << __var->first
                   << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank() // important localRank
                   << "." << std::setfill( '0' ) << std::setw( 3 ) << __step->index();
        Debug( 8006 ) << "[ExporterEnsight::saveNodal] saving " << __varfname.str() << "...\n";
        std::fstream __out( __varfname.str().c_str(), std::ios::out | std::ios::binary );

        char buffer[ 80 ];
        strcpy( buffer, __var->second.name().c_str() );
        __out.write( ( char * ) & buffer, sizeof( buffer ) );

        uint16_type nComponents = __var->second.nComponents;

        if ( __var->second.is_vectorial )
            nComponents = 3;

        size_type __field_size = nComponents*__var->second.size()/__var->second.nComponents;
        ublas::vector<float> __field( __field_size );
        __field.clear();
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
                    size_type ptid = elt_it->point( p ).id();
                    size_type global_node_id = nComponents * ptid + c ;

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

        Debug( 8006 ) << "[ExporterEnsight::saveNodal] saving " << __varfname.str() << "done\n";
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
                    << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank() // important localRank
                    << "." << std::setfill( '0' ) << std::setw( 3 ) << __step->index();
        Debug( 8006 ) << "[ExporterEnsight::saveElement] saving " << __evarfname.str() << "...\n";
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
            Debug( 8006 ) << "part " << buffer << "\n";
            strcpy( buffer, this->elementType().c_str() );
            __out.write( ( char * ) & buffer, sizeof( buffer ) );
            Debug( 8006 ) << "element type " << buffer << "\n";

            uint16_type nComponents = __evar->second.nComponents;

            if ( __evar->second.is_vectorial )
                nComponents = 3;

            size_type __field_size = nComponents*__evar->second.size()/__evar->second.nComponents;
            ublas::vector<float> __field( __field_size );
            __field.clear();
            typename mesh_type::marker_element_const_iterator elt_it;
            typename mesh_type::marker_element_const_iterator elt_en;
            boost::tie( elt_it, elt_en ) = __step->mesh()->elementsWithMarker( p_it->first,
                                           this->worldComm().localRank() ); // important localRank!!!!

            if ( !__evar->second.areGlobalValuesUpdated() )
                __evar->second.updateGlobalValues();

            Debug( 8006 ) << "[saveElement] firstLocalIndex = " << __evar->second.firstLocalIndex() << "\n";
            Debug( 8006 ) << "[saveElement] lastLocalIndex = " << __evar->second.lastLocalIndex() << "\n";
            Debug( 8006 ) << "[saveElement] field.size = " << __field_size << "\n";
            size_type e = 0;

            for ( ; elt_it != elt_en; ++elt_it, ++e )
            {
                Debug( 8006 ) << "pid : " << this->worldComm().globalRank()
                              << " elt_it :  " << elt_it->id()
                              << " e : " << e << "\n";

                for ( int c = 0; c < nComponents; ++c )
                {
                    size_type global_node_id = nComponents * e + c ;

                    if ( c < __evar->second.nComponents )
                    {
                        size_type dof_id = boost::get<0>( __evar->second.functionSpace()->dof()->localToGlobal( elt_it->id(),0, c ) );

                        Debug( 8006 ) << "c : " << c
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
                        Debug( 8006 ) << "c : " << c
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

        Debug( 8006 ) << "[ExporterEnsight::saveElement] saving " << __evarfname.str() << "done\n";
        ++__evar;
    }
}

template<typename MeshType, int N>
void
ExporterEnsight<MeshType,N>::visit( mesh_type* __mesh )
{
    char buffer[ 80 ];
    std::vector<int> idnode, idelem;

    std::fstream __out( _M_filename.c_str(), std::ios::out | std::ios::binary );


    strcpy( buffer, "C Binary" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, _M_filename.c_str() );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "elements" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "node id given" );
    __out.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "element id given" );
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
        __coords[3*__count+0] = ( float ) pt_it->node()[0];

        if ( mesh_type::nRealDim >= 2 )
            __coords[3*__count+1] = ( float ) pt_it->node()[1];

        else
            __coords[3*__count+1] = float( 0 );

        if ( mesh_type::nRealDim >= 3 )
            __coords[3*__count+2] = float( pt_it->node()[2] );

        else
            __coords[3*__count+2] = float( 0 );
    }

    __out.write( ( char * ) & idnode.front(), idnode.size() * sizeof( int ) );
    __out.write( ( char * ) __coords.data().begin(), __coords.size() * sizeof( float ) );


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
                                       this->worldComm().localRank() ); // important localRank!!!!

        //	int __ne = __mesh->numElements();
        //int __ne = p_it->second;
        int __ne = std::distance( elt_it, elt_en );

        Debug( 8006 ) << "num Elements to save : " << __ne << "\n";

        __out.write( ( char * ) &__ne, sizeof( int ) );

        idelem.resize( __ne );


        for ( size_type e = 0; elt_it != elt_en; ++elt_it, ++e )
        {
            idelem[e] = elt_it->id() + 1;
        }

        __out.write( ( char * ) & idelem.front(), idelem.size() * sizeof( int ) );

        //	elt_it = __mesh->beginElement();
        boost::tie( elt_it, elt_en ) = __mesh->elementsWithMarker( p_it->first,
                                       this->worldComm().localRank() ); // important localRank!!!!
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
