/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-21

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2011 Feel++ Consortium

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
#ifndef FEELPP_EXPORTERENSIGHTGOLD_CPP
#define FEELPP_EXPORTERENSIGHTGOLD_CPP 1

#include <feel/feelcore/feel.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/timeset.hpp>
#include <feel/feelfilters/exporterensightgold.hpp>
#include <feel/feelfilters/detail/fileindex.hpp>

namespace Feel
{
template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( WorldComm const& worldComm )
:
super( worldComm ),
M_worldCommBase(worldComm),
M_element_type()
{
    init();
}
template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( std::string const& __p, int freq, WorldComm const& worldComm )
    :
    super( "ensightgold", __p, freq, worldComm ),
    M_worldCommBase(worldComm),
    M_element_type()
{
    init();
}
template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( po::variables_map const& vm, std::string const& exp_prefix, WorldComm const& worldComm )
    :
    super( vm, exp_prefix, worldComm ),
    M_worldCommBase(worldComm)
{
    init();
}

template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( std::string const& exp_prefix, WorldComm const& worldComm )
    :
    super( exp_prefix, worldComm ),
    M_worldCommBase(worldComm)
{
    init();
}

template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( ExporterEnsightGold const & __ex )
    :
    super( __ex ),
    M_worldCommBase( __ex.worldCommBase() ),
    M_element_type( __ex.M_element_type )
{
}

template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::~ExporterEnsightGold()
{}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::init()
{
    /* define ensight named constant for the different faces/elements */
    if ( mesh_type::nDim == 1 )
        if ( mesh_type::Shape == SHAPE_LINE )
        {
            M_element_type = ( mesh_type::nOrder == 1 )?"bar2":"bar3";
            M_face_type = "point";
        }

    if ( mesh_type::nDim == 2 )
    {
        if ( mesh_type::Shape == SHAPE_TRIANGLE )
            M_element_type = ( mesh_type::nOrder == 1 )?"tria3":"tria6";

        else if ( mesh_type::Shape == SHAPE_QUAD )
            M_element_type = ( mesh_type::nOrder == 1 )?"quad4":"quad8";

        M_face_type = ( mesh_type::nOrder == 1 )?"bar2":"bar3";
    }

    if ( mesh_type::nDim == 3 )
    {
        if ( mesh_type::Shape == SHAPE_TETRA )
        {
            M_element_type = ( mesh_type::nOrder == 1 )?"tetra4":"tetra10";
            M_face_type = ( mesh_type::nOrder == 1 )?"tria3":"tria6";
        }

        else if ( mesh_type::Shape == SHAPE_HEXA )
        {
            M_element_type = ( mesh_type::nOrder == 1 )?"hexa8":"hexa20";
            M_face_type = ( mesh_type::nOrder == 1 )?"quad4":"quad8";
        }
    }

    /* TODO Do a cleanup of previous stored files */
    /* to avoid conflicts */
    /* example case: where a new simulation has fewer timesteps */
    /* it can be confusing if files with higher timesteps remain */

    /* Init number of digit for maximum time step */
    M_timeExponent = 4;

    /* if we do not want to merge the results from the different processes */
    /* we isolate each process by using worldCommSeq() from Environment */
    /* Each process will be in a group seeing only itself and not every */
    /* process as worldComm() */
    if( ! boption( _name = "exporter.ensightgold.merge.markers" ) )
    {
        //this->setWorldComm(Environment::worldCommSeq());
        this->setWorldComm(this->worldComm().subWorldCommSeq());
    }
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::save() const
{
    if ( !this->worldComm().isActive() ) return;

    DVLOG(2) << "checking if frequency is ok\n";

    if ( this->cptOfSave() % this->freq()  )
    {
        this->saveTimeSet();
        return;
    }

    /* Check that we have steps to save */
    /* Ensures that we do not end up in a segfault */
    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();
    bool hasSteps = true;
    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;

        typename timeset_type::step_const_iterator __it = __ts->beginStep();
        typename timeset_type::step_const_iterator __end = __ts->endStep();

        mesh_ptrtype mesh = NULL;

        /* check if we have steps for the current dataset */
        if(__it == __end)
        {
            LOG(INFO) << "Timeset " << __ts->name() << " (" << __ts->index() << ") contains no timesteps (Consider using add() or addRegions())" << std::endl;
            hasSteps = false;

            /* if we have no time steps, we still save the geometry */
            if(__ts->hasMesh())
            {
                mesh = __ts->mesh(); 
            }
        }
        /* if we have steps */
        else
        {
            /* get the step that will be saved */
            __it = boost::prior( __end );

            /* get the mesh for the timestep to save */
            typename timeset_type::step_ptrtype __step = *__it;
            mesh = __step->mesh();
        }

        /* If we haven't yet computed the markers to be written */
        /* or we restart a simulation, we need to update the markers to be written */
        if( mesh && M_markersToWrite.size() == 0)
        {
            this->computeMarkersToWrite(mesh);
        }

        ++__ts_it;
    }

    /*
    if(!hasSteps)
    {
        return;
    }
    */

    boost::timer ti;
    DVLOG(2) << "export in ensight format\n";

    ti.restart();
    DVLOG(2) << "export geo(mesh) file\n";
    writeGeoFiles();
    DVLOG(2) << "export geo(mesh) file ok, time " << ti.elapsed() << "\n";

    LOG(INFO) << "Geo File written" << std::endl;

    ti.restart();
    DVLOG(2) << "export variable file\n";
    /* only try to write variable data when we have time steps */
    if(hasSteps)
    { writeVariableFiles(); }
    DVLOG(2) << "export variable files ok, time " << ti.elapsed() << "\n";

    ti.restart();
    DVLOG(2) << "export time set\n";
    this->saveTimeSet();
    DVLOG(2) << "export time set ok, time " << ti.elapsed() << "\n";

    ti.restart();
    DVLOG(2) << "export case file\n";
    writeCaseFile();
    DVLOG(2) << "export case file ok, time " << ti.elapsed() << "\n";

    DVLOG(2) << "export sos\n";
    writeSoSFile();
    DVLOG(2) << "export sos ok, time " << ti.elapsed() << "\n";
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeSoSFile() const
{
    // only on proc 0
    if ( this->worldComm().isMasterRank() )
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
              << "type: master_server gold \n\n";
        if( boption( _name="exporter.ensightgold.merge.markers") )
        {
            __out << "MULTIPLE_CASEFILES\n"
                << "total number of cfiles: 1 " << std::endl
                << "cfiles global path: " << fs::current_path().string() << "\n"
                << "cfiles pattern: "<<this->prefix() << ".case\n"
                << "cfiles start number: 0\n"
                << "cfiles increment: 1\n\n";
        }
        else
        {
            __out << "MULTIPLE_CASEFILES\n"
                << "total number of cfiles: " << this->worldCommBase().globalSize() << "\n"
                << "cfiles global path: " << fs::current_path().string() << "\n"
                << "cfiles pattern: "<<this->prefix() << "-" << this->worldCommBase().globalSize() << "_*.case\n"
                << "cfiles start number: 0\n"
                << "cfiles increment: 1\n\n";
        }
        __out << "SERVERS\n"
              << "number of servers: "<< (this->worldCommBase().globalSize()/100)+1 <<" repeat\n";

        //
        // save also a sos that paraview can understand, the previous format
        // does not seem to be supported by paraview
        //
        std::ostringstream filestrparaview;
        filestrparaview << this->path() << "/" << this->prefix() << "-paraview-" << this->worldComm().globalSize() << ".sos";
        std::ofstream __outparaview( filestrparaview.str().c_str() );

        __outparaview << "FORMAT:\n"
                      << "type: master_server gold \n"
                      << "SERVERS\n";

        /* set the number of servers to one */
        /* if we merged the markers in a single file */
        if( boption( _name="exporter.ensightgold.merge.markers") )
        {
            __outparaview << "number of servers: 1 " << std::endl;

            __outparaview << "#Server " << 1 << "\n"
                << "machine id: " << mpi::environment::processor_name() << "\n"
                << "executable: /usr/local/bin/ensight76/bin/ensight7.server\n"
                << "data_path: " << fs::current_path().string() << "\n"
                << "casefile: " << this->prefix() << ".case\n";
        }
        else
        {
            __outparaview << "number of servers: " << this->worldCommBase().globalSize() << "\n";

            for ( int pid = 0 ; pid < this->worldCommBase().globalSize(); ++pid )
            {
                __outparaview << "#Server " << pid+1 << "\n"
                    << "machine id: " << mpi::environment::processor_name() << "\n"
                    << "executable: /usr/local/bin/ensight76/bin/ensight7.server\n"
                    << "data_path: " << fs::current_path().string() << "\n"
                    << "casefile: " << this->prefix() << "-" << this->worldCommBase().globalSize() << "_" << pid << ".case\n";
            }
        }
    }
}
template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeCaseFile() const
{
    // only on proc 0
    if( ! boption( _name="exporter.ensightgold.merge.markers") || ( boption( _name="exporter.ensightgold.merge.markers") && this->worldComm().isMasterRank() ) )
    {
        std::ostringstream filestr;

        timeset_const_iterator __ts_it = this->beginTimeSet();
        timeset_const_iterator __ts_en = this->endTimeSet();

        filestr << this->path() << "/"
            << this->prefix();
        if( ! boption( _name="exporter.ensightgold.merge.markers") )
        { filestr << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
        filestr << ".case";

        std::ofstream __out( filestr.str().c_str() );

        if ( __out.fail() )
        {
            DVLOG(2) << "cannot open " << filestr.str()  << "\n";
            exit( 0 );
        }

        __out << "FORMAT:\n"
            << "type: ensight gold\n"
            << "GEOMETRY:\n";

        switch ( this->exporterGeometry() )
        {
            case EXPORTER_GEOMETRY_STATIC:
                {
                    timeset_ptrtype __ts = *__ts_it;
                    __out << "model: " << __ts->name();
                    if( ! boption( _name="exporter.ensightgold.merge.markers") )
                    { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                    __out << ".geo";
                }
                break;
            default:
            case EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY:
            case EXPORTER_GEOMETRY_CHANGE:
                {
                    while ( __ts_it != __ts_en )
                    {
                        timeset_ptrtype __ts = *__ts_it;

                        if( boption( _name="exporter.ensightgold.merge.timesteps") )
                        {
                            __out << "model: " << __ts->index() << " 1 " << __ts->name();
                            if( ! boption( _name="exporter.ensightgold.merge.markers") )
                            { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                            /* if we want to pack data in several files instead of one, we add an index to the filename */
                            if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
                            { __out << ".*"; }
                            __out << ".geo";
                        }
                        else
                        {
                            __out << "model: " << __ts->index() << " " << __ts->name();
                            if( ! boption( _name="exporter.ensightgold.merge.markers") )
                            { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                            __out << ".geo" << "." << std::string(M_timeExponent, '*');
                        }

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

                // treat constant per case
                auto s_it = ( *__tstp_it )->beginScalar();
                auto s_en = ( *__tstp_it )->endScalar();
                for ( ; s_it != s_en; ++s_it )
                {
                    if ( s_it->second.second )
                    {
                        // constant over time
                        __out << "constant per case: " << s_it->first << " " << s_it->second.first << "\n";
                    }
                    else
                    {
                        if ( this->worldComm().isMasterRank() )
                        {

                            __out << "constant per case file: " << __ts->index() << " " << s_it->first << " " << s_it->first << ".scl";
                            // loop over time
                            auto stepit = __ts->beginStep();
                            auto stepen = __ts->endStep();

                            std::ofstream ofs;
                            int d = std::distance( stepit, stepen );
                            LOG(INFO) << "distance = " << d;
                            if ( d > 1 )
                                ofs.open( s_it->first+".scl", std::ios::out | std::ios::app );
                            else
                                ofs.open( s_it->first+".scl", std::ios::out );

                            auto step = *boost::prior(stepen);
                            ofs << step->scalar( s_it->first ) << "\n";
                            ofs.close();
                            __out << "\n";
                        }

                    }
                }

                typename timeset_type::step_type::nodal_scalar_const_iterator __it = ( *__tstp_it )->beginNodalScalar();
                typename timeset_type::step_type::nodal_scalar_const_iterator __end = ( *__tstp_it )->endNodalScalar();

                while ( __it != __end )
                {
                    if( boption( _name="exporter.ensightgold.merge.timesteps") )
                    {
                        __out << "scalar per node: "
                            << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                            << __it->second.name() << " " << __it->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        /* if we want to pack data in several files instead of one, we add an index to the filename */
                        if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
                        { __out << ".*"; }
                        __out << ".scl" << std::endl;
                    }
                    else
                    {
                        __out << "scalar per node: "
                            << __ts->index() << " " // << *__ts_it->beginStep() << " "
                            << __it->second.name() << " " << __it->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        __out << ".scl" << "." << std::string(M_timeExponent, '*') << std::endl;
                    }
                    ++__it;
                }

                typename timeset_type::step_type::nodal_vector_const_iterator __itv = ( *__tstp_it )->beginNodalVector();
                typename timeset_type::step_type::nodal_vector_const_iterator __env = ( *__tstp_it )->endNodalVector();

                while ( __itv != __env )
                {
                    if( boption( _name="exporter.ensightgold.merge.timesteps") )
                    {
                        __out << "vector per node: "
                            << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                            << __itv->second.name() << " " << __itv->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        /* if we want to pack data in several files instead of one, we add an index to the filename */
                        if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
                        { __out << ".*"; }
                        __out << ".vec" << std::endl;
                    }
                    else
                    {
                        __out << "vector per node: "
                            << __ts->index() << " " // << *__ts_it->beginStep() << " "
                            << __itv->second.name() << " " << __itv->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        __out << ".vec" << "." << std::string(M_timeExponent, '*') << std::endl;
                    }
                    ++__itv;
                }

                typename timeset_type::step_type::nodal_tensor2_const_iterator __itt = ( *__tstp_it )->beginNodalTensor2();
                typename timeset_type::step_type::nodal_tensor2_const_iterator __ent = ( *__tstp_it )->endNodalTensor2();

                while ( __itt != __ent )
                {
                    if( boption( _name="exporter.ensightgold.merge.timesteps") )
                    {
                        __out << "tensor per node: "
                            << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                            << __itt->second.name() << " " << __itt->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        /* if we want to pack data in several files instead of one, we add an index to the filename */
                        if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
                        { __out << ".*"; }
                        __out << ".tsr" << std::endl;
                    }
                    else
                    {
                        __out << "tensor per node: "
                            << __ts->index() << " " // << *__ts_it->beginStep() << " "
                            << __itt->second.name() << " " << __itt->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        __out << ".tsr" << "." << std::string(M_timeExponent, '*') << std::endl; 
                    }
                    ++__itt;
                }

                typename timeset_type::step_type::element_scalar_const_iterator __it_el = ( *__tstp_it )->beginElementScalar();
                typename timeset_type::step_type::element_scalar_const_iterator __end_el = ( *__tstp_it )->endElementScalar();

                while ( __it_el != __end_el )
                {
                    if( boption( _name="exporter.ensightgold.merge.timesteps") )
                    {
                        __out << "scalar per element: "
                            << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                            << __it_el->second.name() << " " << __it_el->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        /* if we want to pack data in several files instead of one, we add an index to the filename */
                        if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
                        { __out << ".*"; }
                        __out << ".scl" << std::endl;
                    }
                    else
                    {
                        __out << "scalar per element: "
                            << __ts->index() << " " // << *__ts_it->beginStep() << " "
                            << __it_el->second.name() << " " << __it_el->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        __out << ".scl" << "." << std::string(M_timeExponent, '*') << std::endl;
                    }
                    ++__it_el;
                }

                typename timeset_type::step_type::element_vector_const_iterator __itv_el = ( *__tstp_it )->beginElementVector();
                typename timeset_type::step_type::element_vector_const_iterator __env_el = ( *__tstp_it )->endElementVector();

                while ( __itv_el != __env_el )
                {
                    if( boption( _name="exporter.ensightgold.merge.timesteps") )
                    {
                        __out << "vector per element: "
                            << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                            << __itv_el->second.name() << " " << __itv_el->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        /* if we want to pack data in several files instead of one, we add an index to the filename */
                        if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
                        { __out << ".*"; }
                        __out << ".vec" << std::endl;
                    }
                    else
                    {
                        __out << "vector per element: "
                            << __ts->index() << " " // << *__ts_it->beginStep() << " "
                            << __itv_el->second.name() << " " << __itv_el->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        __out << ".vec" << "." << std::string(M_timeExponent, '*') << std::endl;
                    }
                    ++__itv_el;
                }

                typename timeset_type::step_type::element_tensor2_const_iterator __itt_el = ( *__tstp_it )->beginElementTensor2();
                typename timeset_type::step_type::element_tensor2_const_iterator __ent_el = ( *__tstp_it )->endElementTensor2();

                while ( __itt_el != __ent_el )
                {
                    if( boption( _name="exporter.ensightgold.merge.timesteps") )
                    {
                        __out << "tensor per element: "
                            << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                            << __itt_el->second.name() << " " << __itt_el->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        /* if we want to pack data in several files instead of one, we add an index to the filename */
                        if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
                        { __out << ".*"; }
                        __out << ".tsr" << std::endl;
                    }
                    else
                    {
                        __out << "tensor per element: "
                            << __ts->index() << " " // << *__ts_it->beginStep() << " "
                            << __itt_el->second.name() << " " << __itt_el->first;
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().globalRank(); }
                        __out << ".tsr" << "." << std::string(M_timeExponent, '*') << std::endl;
                    }
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

            ++__ts_it;
        }
        __out << "\n";

        if( boption( _name="exporter.ensightgold.merge.timesteps") )
        {
            auto ts = *(this->beginTimeSet());

            // create several filesets if we pack data into groups in several files
            if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
            {
                int stepsPerFile = ioption( _name="exporter.ensightgold.pack.timesteps" );
                int nbFiles = ts->numberOfSteps() / stepsPerFile;
                int r = ts->numberOfSteps() % stepsPerFile;

                __out << "FILE\n";
                __out << "file set: 1\n";
                // create a file set for each step pack
                int i = 0;
                for( i = 0; i < nbFiles; i++)
                {
                    __out << "filename index: " << i + 1 << "\n";
                    __out << "number of steps: " << stepsPerFile << "\n";
                }

                // if we have remainder steps to pack
                // we create an additional file pack
                if( r != 0 )
                {
                    __out << "filename index: " << i + 1 << "\n";
                    __out << "number of steps: " << r << "\n";
                }
            }
            else
            {
                __out << "FILE\n";
                __out << "file set: 1\n";
                __out << "number of steps: " << ts->numberOfSteps() << "\n";
            }
        }

        __out << "\n";

        if ( boption( _name="exporter.ensightgold.use-sos" ) == false )
        {
            // In the following line, we substituted the Environment::numberOfProcessors
            // by the size of the worldComm passed to the exporter to ensure that we are using
            // only the data for the current processor
            if ( ( this->worldComm().globalSize() > 1 )  && ( this->worldComm().globalRank() == 0 ) )
            {
                __out << "APPENDED_CASEFILES\n"
                    << "total number of cfiles: " << this->worldComm().globalSize()-1 << "\n"
                    // no need for that
                    // << "cfiles global path: " << fs::current_path().string() << "\n"
                    << "cfiles: ";
                for(int p = 1; p < this->worldComm().globalSize(); ++p )
                {
                    std::ostringstream filestr;
                    filestr << this->prefix() << "-"
                        << this->worldComm().globalSize() << "_" << p << ".case";
                    __out << filestr.str() << "\n        ";
                }
            }
        } // use-sos

        __out.close();
    }
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeGeoFiles() const
{
    namespace lambda = boost::lambda;

    int size;
    char buffer[80];

    MPI_File fh;
    MPI_Status status;
    MPI_Info info;
    MPI_Offset offset = 0;

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

        /* static geometry */
        if ( this->exporterGeometry() == EXPORTER_GEOMETRY_STATIC )
        {
            LOG(INFO) << "GEO: Static geo mode" << std::endl;

            /* only write the geometry in the first timestep */
            if( __it == __ts->beginStep() ) 
            {
                /* generate geo filename */
                std::ostringstream __geofname;
                __geofname << this->path() << "/" << __ts->name();
                if( ! boption( _name="exporter.ensightgold.merge.markers") )
                { __geofname << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().localRank(); }
                __geofname << ".geo";
                M_filename =  __geofname.str();
                //CHECK( (*__it)->hasMesh() || __ts->hasMesh()  ) << "Invalid mesh data structure in static geometry mode\n";

                /* Open File with MPI IO */
                char * str = strdup(__geofname.str().c_str());

                /* Check if file exists and delete it, if so */
                /* (MPI IO does not have a truncate mode ) */
                if(this->worldComm().isMasterRank() && fs::exists(str))
                {
                    MPI_File_delete(str, MPI_INFO_NULL);
                }
                MPI_Barrier( this->worldComm().comm() );

                MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
                free(str);

                Feel::detail::FileIndex index;

                if( boption( _name="exporter.ensightgold.merge.timesteps" ) )
                {
                    // first read
                    index.read(fh);
                    // we position the cursor at the beginning of the file
                    MPI_File_seek_shared(fh, 0, MPI_SEEK_SET);

                    /* write C binary if we didn't find the index <=> first pass on the file */
                    if( !index.defined() )
                    {
                        if( this->worldComm().isMasterRank() )
                        { size = sizeof(buffer); }
                        else
                        { size = 0; }
                        memset(buffer, '\0', sizeof(buffer));
                        strcpy(buffer, "C Binary");
                        MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                    }

                    /* Write time step start */
                    if( this->worldComm().isMasterRank() )
                    { size = sizeof(buffer); }
                    else
                    { size = 0; }
                    memset(buffer, '\0', sizeof(buffer));
                    strcpy(buffer,"BEGIN TIME STEP");
                    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                    LOG(INFO) << "saveNodal out: " << buffer;

                    /* add the beginning of the new block to the file */
                    MPI_File_get_position_shared(fh, &offset);
                    index.add( offset );
                }

                /* Write the file */
                this->writeGeoMarkers(fh, mesh);

                if( boption( _name="exporter.ensightgold.merge.timesteps" ) )
                {
                    /* write timestep end */
                    if( this->worldComm().isMasterRank() )
                    { size = sizeof(buffer); }
                    else
                    { size = 0; }
                    memset(buffer, '\0', sizeof(buffer));
                    strcpy(buffer,"END TIME STEP");
                    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

                    // write back the file index
                    index.write( fh );
                }

                /* close file */
                MPI_File_close(&fh);
            }
        }
        /* changing geometry */
        else
        {
            LOG(INFO) << "GEO: Changing geo mode" << std::endl;
            /* Transient mode */
            if( boption( _name="exporter.ensightgold.merge.timesteps") )
            {
                /* TODO */
                /* MPI_File_Open -> f */
                /* P0 : Write "C Binary" */
                /* for each Timestep T in TS */
                /* P0 : Write "BEGIN TIME STEP" */
                /* WriteGeo(T, f) */
                /* P0 : Write "END TIME STEP" */

                std::ostringstream __geofname;

                __geofname << this->path() << "/"
                    << __ts->name();
                if( ! boption( _name="exporter.ensightgold.merge.markers") )
                { __geofname << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().localRank(); }
                /* if we want to pack data in several files instead of one */
                /* we compute an index to add to the filename */
                if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
                {
                    // timestep indices start at 1
                    __geofname << "." << ((timeIndex - TS_INITIAL_INDEX) / ioption( _name="exporter.ensightgold.pack.timesteps" ) + 1);
                }
                __geofname << ".geo";

                M_filename =  __geofname.str();

                /* Open File with MPI IO */
                char * str = strdup(M_filename.c_str());

                /* Check if file exists and delete it, if so */
                /* (MPI IO does not have a truncate mode ) */
                if(this->worldComm().isMasterRank() && __it == __ts->beginStep() && fs::exists(str))
                {
                    MPI_File_delete(str, MPI_INFO_NULL);
                }
                MPI_Barrier(this->worldComm().comm());

                MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
                free(str);

                Feel::detail::FileIndex index;

                // first read
                index.read(fh);

                /* Move to the beginning of the fie index section */
                /* to overwrite it */
                if ( index.defined() && (timeIndex - TS_INITIAL_INDEX) > 0 ) {
                    MPI_File_seek_shared(fh, index.fileblock_n_steps, MPI_SEEK_SET);
                }
                else {
                    // we position the cursor at the beginning of the file
                    MPI_File_seek_shared(fh, 0, MPI_SEEK_SET);
                }

                /* write C binary if we didn't find the index <=> first pass on the file */
                if( !index.defined() )
                {
                    if( this->worldComm().isMasterRank() )
                    { size = sizeof(buffer); }
                    else
                    { size = 0; }
                    memset(buffer, '\0', sizeof(buffer));
                    strcpy(buffer, "C Binary");
                    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                }

                /* Write time step start */
                if( this->worldComm().isMasterRank() )
                { size = sizeof(buffer); }
                else
                { size = 0; }
                memset(buffer, '\0', sizeof(buffer));
                strcpy(buffer,"BEGIN TIME STEP");
                MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                LOG(INFO) << "saveNodal out: " << buffer;

                /* add the beginning of the new block to the file */
                MPI_File_get_position_shared(fh, &offset);
                index.add( offset );

                /* write data for timestep */
                this->writeGeoMarkers(fh, mesh);

                /* write timestep end */
                if( this->worldComm().isMasterRank() )
                { size = sizeof(buffer); }
                else
                { size = 0; }
                memset(buffer, '\0', sizeof(buffer));
                strcpy(buffer,"END TIME STEP");
                MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

                // write back the file index
                index.write( fh );
            }
            /* non transient */
            else
            {
                // TODO
                /* for each Timestep T in TS */
                /* MPI_File_Open -> f */
                /* WriteGeo(T, f) */

                std::ostringstream __geofname;

                __geofname << this->path() << "/" << __ts->name();
                if( ! boption( _name="exporter.ensightgold.merge.markers") )
                { __geofname << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().localRank(); }
                __geofname << ".geo" << "." << std::setfill( '0' ) << std::setw( M_timeExponent ) << timeIndex;

                /* no need to check if step is in memory */
                /* as it is performed at the beginning of the function */
                M_filename =  __geofname.str();

                /* Open File with MPI IO */
                char * str = strdup(M_filename.c_str());

                /* Check if file exists and delete it, if so */
                /* (MPI IO does not have a truncate mode ) */
                if(this->worldComm().isMasterRank() && fs::exists(str))
                {
                    MPI_File_delete(str, MPI_INFO_NULL);
                }
                MPI_Barrier(this->worldComm().comm());

                MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
                free(str);

                /* Either write every marker in one file */
                this->writeGeoMarkers(fh, mesh);

                /* close file */
                MPI_File_close(&fh);
            }
        }

        ++__ts_it;
    }
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType, N>::computeMarkersToWrite(mesh_ptrtype mesh) const
{
    /* Protect the function call in case we have a NULL mesh */
    if(mesh == NULL)
    {
        return;
    }

    // TODO write the faces
    // Integrate them into the parts corresponding to the elements ?
    int nbmarkers = mesh->markerNames().size();

    LOG(INFO) << "nMarkers " << nbmarkers << std::endl;

    if(boption( _name="exporter.ensightgold.merge.markers" ))
    {
        // TODO find a better place for this code to be executed
        // as we have allreduce, this can take serious execution time
        /* determine for which markers, we have data to write */
#if 0
        /* Display info about markers */
        std::ostringstream ossmn;
        ossmn << this->worldComm().rank() << " markers ";
        BOOST_FOREACH( auto marker, mesh->markerNames() )
        {
            ossmn << " " << marker.second[0];
        }
        LOG(INFO) << ossmn.str() << std::endl;
#endif

#if 0
        M_markersToWrite.clear();
        std::ostringstream osselts;
        osselts << mesh->worldComm().rank();
        BOOST_FOREACH( auto marker, mesh->markerNames() )
        {
            /* Check whether at least one process has elements to write */
            int localNElts = std::distance(mesh->beginElementWithMarker(marker.second[0]), mesh->endElementWithMarker(marker.second[0]));
            int globalNElts = 0;

            osselts << " " << marker.second[0] << " (" << localNElts << ")";

            mpi::all_reduce(this->worldComm(), localNElts, globalNElts, mpi::maximum<int>());

            /* if we have at least one element for the current marker */
            /* all the processes need to parse it to avoid deadlocks with gather in MeshPoints */
            if(globalNElts)
            {
                M_markersToWrite.push_back(marker.second[0]);
            }
        }
        LOG(INFO) << osselts.str() << std::endl;
#endif

        // TODO If the number of parts per process is not the sam
        // some processes might not enter some instance of the function in the following loop
        // and inside it there is a allgather that would cause the processes to be stuck in it
        // check the number of parts
        /*
        int localNParts = std::distance(p_it, p_en);
        int globalNParts = 0;

        mpi::all_reduce(this->worldComm(), localNParts, globalNParts, mpi::maximum<int>());

        LOG(INFO) << this->worldComm().rank() << " " << localNParts << " " << globalNParts << " " << mesh->markerNames().size() << std::endl;
        */

        // TODO Removed this loop, as it was causing MPI deadlocks when the numebr of parts was different
        // from one process to the other
        /* Write elements */
        /*
        typename mesh_type::parts_const_iterator_type p_it = mesh->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = mesh->endParts();

        for ( ; p_it != p_en; ++p_it )
        {
            this->writeGeoMarkedElements(fh, mesh, p_it);
        }
        */

        /* Check whether successive part id are the same on different processes */
        /* does not seem to be the case */
        typename mesh_type::parts_const_iterator_type p_st = mesh->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = mesh->endParts();

        std::ostringstream osspi;
        osspi << this->worldComm().rank() << " partid";
        for(auto p_it = p_st ; p_it != p_en; ++p_it )
        {
            osspi << " " << p_it->first << "(" << p_it->second << ")";
        }
        LOG(INFO) << osspi.str() << std::endl;

        /* iterate over the local markers to get the different markers needed to be written */
        std::vector<int> localMarkers;
        for(auto p_it = p_st ; p_it != p_en; ++p_it )
        {
            localMarkers.push_back(p_it->first);
        }

        /* gather all the markers to be written on the different processes */
        /* to order the writing step */
        std::vector<std::vector<int> > globalMarkers;
        mpi::all_gather(this->worldComm(), localMarkers, globalMarkers);

        for(int i = 0; i < globalMarkers.size(); i++)
        {
            for(int j = 0; j < globalMarkers[i].size(); j++)
            {
                M_markersToWrite.insert(globalMarkers[i][j]);
            }
        }

        std::ostringstream osss;
        osss << this->worldComm().rank() << " parts/markers";
        for(std::set<int>::iterator it = M_markersToWrite.begin(); it != M_markersToWrite.end(); it++)
        {
            osss << " " << *it << " (" << std::distance(mesh->beginElementWithMarker(*it), mesh->endElementWithMarker(*it)) << ")";
        }
        LOG(INFO) << osss.str() << std::endl;
    }
    else
    {
        /* iterate over the markers to get the different markers needed to be written */
        /* (in this case all the markers local to the processor) */
        typename mesh_type::parts_const_iterator_type p_st = mesh->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = mesh->endParts();
        for(auto p_it = p_st ; p_it != p_en; ++p_it )
        {
            M_markersToWrite.insert(p_it->first);
        }
    }

}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeGeoMarkers(MPI_File fh, mesh_ptrtype mesh) const
{
    /* Function WriteGeoMarker(Timestep T, MPI_File f) */
        /* P0 : Write "C Binary", part, #, desc, coords */
        /* for each part Marker in T */
            /* P0 : Write "part", Write # */
            /* P* : Compute nn, build points id array, build points array */
            /* P0 : Write nn */
            /* P* : Write_shared id */
            /* P* : Write_shared points */

            /* for each face/element type */
                /* P* : compute ne, build elements id array, build elements array */
                /* P0 : Write element type */
                /* P0 : Write ne */
                /* P* : Write_shared id */
                /* P* : Write_shared elements */
                /* WriteGeoPart(P, f) */

    /* Write file header */
    this->writeGeoHeader(fh);

    /* Write faces */
    if ( boption( _name="exporter.ensightgold.save-face" ) )
    {
        for( std::pair<const std::string, std::vector<size_type> > & m : mesh->markerNames() )
        {
            this->writeGeoMarkedFaces(fh, mesh, m);
        }
    }

    /* Working with marker names instead */
    for(std::set<int>::iterator mit = M_markersToWrite.begin(); mit != M_markersToWrite.end(); mit++)
    {
        this->writeGeoMarkedElements(fh, mesh, *mit);
    }
}


template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType, N>::writeGeoHeader(MPI_File fh) const
{
    MPI_Status status;

    int size = 0;
    char buffer[80];

    DVLOG(2) << "Merging markers : " << "\n";

    /* write header */
    /* little trick to only perform collective operation (optimized) */
    /* and avoid scattering offset and reseting the shared pointer if we would write this only on master proc */
    if( this->worldComm().isMasterRank() )
    { size = sizeof(buffer); }
    else
    { size = 0; }

    // only write C Binary if we are not mergin timesteps
    // as it is oalready writtent at the beginning of the file
    if( ! boption( _name="exporter.ensightgold.merge.timesteps") )
    {
        memset(buffer, '\0', sizeof(buffer));
        strcpy(buffer, "C Binary");
        MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
        //LOG(INFO) << "wrote " << buffer << std::endl;
    }

    // get only the filename (maybe with full path)
    fs::path gp = M_filename;
    std::string theFileName = gp.filename().string();
    CHECK( theFileName.length() < 80 ) << "the file name is too long : theFileName=" << theFileName << "\n";

    memset(buffer, '\0', sizeof(buffer));
    strcpy( buffer, theFileName.c_str() );
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

    memset(buffer, '\0', sizeof(buffer));
    strcpy( buffer, "elements" );
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

    memset(buffer, '\0', sizeof(buffer));
    strcpy( buffer, "node id given" );
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

    memset(buffer, '\0', sizeof(buffer));
    strcpy( buffer, "element id given" );
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeGeoMarkedFaces(MPI_File fh, mesh_ptrtype mesh, std::pair<const std::string, std::vector<size_type> > & m) const
{
    int size;
    char buffer[80];

    MPI_Status status;

    std::vector<int32_t> idnode, idelem;

    LOG(INFO) << "Marker " << m.first << std::endl;
    /* save faces */
    if ( m.second[1] != mesh->nDim-1 )
        return;

    VLOG(1) << "writing face with marker " << m.first << " with id " << m.second[0];
    auto pairit = mesh->facesWithMarker( m.second[0], this->worldComm().localRank() );
    auto fit = pairit.first;
    auto fen = pairit.second;
    Feel::detail::MeshPoints<float> mp( mesh.get(), this->worldComm(), fit, fen, true, true, true );
    int32_t __ne = std::distance( fit, fen );
    int nverts = fit->numLocalVertices;
    DVLOG(2) << "Faces : " << __ne << "\n";

    if( this->worldComm().isMasterRank() )
    { size = sizeof(buffer); }
    else
    { size = 0; }

    memset(buffer, '\0', sizeof(buffer));
    strcpy( buffer, "part" );
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

    int32_t partid = m.second[0]; 
    if( this->worldComm().isMasterRank() )
    { size = 1; }
    else
    { size = 0; }
    MPI_File_write_ordered(fh, &partid, size, MPI_INT32_T, &status);

    if( this->worldComm().isMasterRank() )
    { size = sizeof(buffer); }
    else
    { size = 0; }
    memset(buffer, '\0', sizeof(buffer));
    strncpy(buffer, m.first.c_str(), sizeof(buffer) - 1 );
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

    if( this->worldComm().isMasterRank() )
    { size = sizeof(buffer); }
    else
    { size = 0; }
    memset(buffer, '\0', sizeof(buffer));
    strcpy(buffer, "coordinates");
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

    // write points coordinates
    fit = pairit.first;
    fen = pairit.second;

    size_type __nv = mp.ids.size();
    int32_t gnop = (int32_t)(mp.globalNumberOfPoints());
    if( this->worldComm().isMasterRank() )
    { size = 1; }
    else
    { size = 0; }
    // write number of points
    MPI_File_write_ordered(fh, &gnop, size, MPI_INT32_T, &status);
    /* write points ids */
    MPI_File_write_ordered(fh, mp.ids.data(), mp.ids.size(), MPI_INT32_T, &status );
    /* write points coordinates in the order x1 ... xn y1 ... yn z1 ... zn */
    for(int i = 0; i < 3; i++)
    {
        MPI_File_write_ordered(fh, mp.coords.data() + i * __nv, __nv, MPI_FLOAT, &status );
    }

    // write connectivity
    fit = pairit.first;
    fen = pairit.second;

    if( this->worldComm().isMasterRank() )
    { size = sizeof(buffer); }
    else
    { size = 0; }
    memset(buffer, '\0', sizeof(buffer));
    strcpy( buffer, M_face_type.c_str() );
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
    VLOG(1) << "face type " << buffer;

    if( this->worldComm().isMasterRank() )
    { size = 1; }
    else
    { size = 0; }
    MPI_File_write_ordered(fh, &__ne, size, MPI_INT32_T, &status);
    VLOG(1) << "n faces " << __ne;

    idnode.resize( __ne );
    fit = pairit.first;
    size_type e = 0;
    for ( ; fit != fen; ++fit, ++e )
    {
        idnode[e] = fit->id() + 1;
    }
    CHECK( e == idnode.size() ) << "Invalid number of face id for part " << m.first;
    MPI_File_write_ordered(fh, idnode.data(), idnode.size(), MPI_INT32_T, &status);

    idelem.resize( __ne*nverts );
    fit = pairit.first;
    e = 0;
    for( ; fit != fen; ++fit, ++e )
    {
        for ( size_type j = 0; j < nverts; j++ )
        {
            // ensight id start at 1
            idelem[e*nverts+j] = mp.old2new[fit->point( j ).id()];
        }
    }
    CHECK( e*nverts == idelem.size() ) << "Invalid number of faces " << e*nverts << " != " << idelem.size() << " in connectivity for part " << m.first;
    MPI_File_write_ordered(fh, idelem.data(), idelem.size(), MPI_INT32_T, &status);
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeGeoMarkedElements(MPI_File fh, mesh_ptrtype mesh, size_type markerid) const
{
    MPI_Status status;

    int size = 0;
    char buffer[80];

    std::vector<int32_t> idnode, idelem;

    //auto r = markedelements(mesh, part->first, EntityProcessType::ALL );
    auto r = markedelements(mesh, markerid, EntityProcessType::ALL );
    auto allelt_it = r.template get<1>();
    auto allelt_en = r.template get<2>();

    //VLOG(1) << "material : " << m << " total nb element: " << std::distance(allelt_it, allelt_en );
    //VLOG(1) << "material : " << m << " ghost nb element: " << std::distance(gelt_it, gelt_en );
    //VLOG(1) << "material : " << m << " local nb element: " << std::distance(lelt_it, lelt_en );
    Feel::detail::MeshPoints<float> mp( mesh.get(), this->worldComm(), allelt_it, allelt_en, true, true, true );
    VLOG(1) << "mesh pts size : " << mp.ids.size();

    // part
    memset(buffer, '\0', sizeof(buffer));
    strcpy( buffer, "part" );

    if( this->worldComm().isMasterRank() )
    { size = sizeof(buffer); }
    else
    { size = 0; }
    // write number of points
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

    // Was previously using p_it->first as partid
    // part id
    int32_t partid = markerid;

    if ( this->worldComm().isMasterRank() )
        LOG(INFO) << "writing part " << partid << std::endl;

    if( this->worldComm().isMasterRank() )
    { size = 1; }
    else
    { size = 0; }
    // write number of points
    MPI_File_write_ordered(fh, &partid, size, MPI_INT32_T, &status);

    // material
    memset(buffer, '\0', sizeof(buffer));
    //sprintf(buffer, "Material %d", part->first);
    sprintf(buffer, "Material %d", (int)(markerid));
    if( this->worldComm().isMasterRank() )
    { size = sizeof(buffer); }
    else
    { size = 0; }
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

    memset(buffer, '\0', sizeof(buffer));
    strcpy( buffer, "coordinates" );
    if( this->worldComm().isMasterRank() )
    { size = sizeof(buffer); }
    else
    { size = 0; }
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

    size_type __nv = mp.ids.size();
    int32_t gnop = (int32_t)(mp.globalNumberOfPoints());
    if( this->worldComm().isMasterRank() )
    { size = 1; }
    else
    { size = 0; }
    MPI_File_write_ordered(fh, &gnop, size, MPI_INT32_T, &status );
    // now we need to move with respect to the processors for the coordinates

    // TODO modify this code !
    // Integrate id translation in the MeshPoints constructor
    /* write points ids */
    std::vector<int32_t> pointids;
    for(int i = 0 ; i < mp.ids.size() ; ++i )
    {
        pointids.push_back(mp.ids.at(i) + mp.offsets_pts);
    }
    MPI_File_write_ordered(fh, pointids.data(), pointids.size(), MPI_INT32_T, &status );

    /* write points coordinates in the order x1 ... xn y1 ... yn z1 ... zn */
    for(int i = 0; i < 3; i++)
    {
        MPI_File_write_ordered(fh, mp.coords.data() + i * __nv, __nv, MPI_FLOAT, &status );
    }

    /* write element type */
    memset(buffer, '\0', sizeof(buffer));
    strcpy( buffer, this->elementType().c_str() );
    if( this->worldComm().isMasterRank() )
    { size = sizeof(buffer); }
    else
    { size = 0; }
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status );

    /* compute and write number of elements */
    //auto r2 = markedelements(mesh, part->first, EntityProcessType::LOCAL_ONLY );
    auto r2 = markedelements(mesh, markerid, EntityProcessType::LOCAL_ONLY );
    auto lelt_it = r2.template get<1>();
    auto lelt_en = r2.template get<2>();

    Feel::detail::MeshPoints<float> mpl( mesh.get(), this->worldComm(), lelt_it, lelt_en, true, true, true );
    int32_t gnole = mpl.globalNumberOfElements();
    //LOG(INFO) << "Global nb elements: " << gnole << std::endl;
    if( this->worldComm().isMasterRank() )
    { size = 1; }
    else
    { size = 0; }
    MPI_File_write_ordered(fh, &gnole, size, MPI_INT32_T, &status );
    // now we need to move with respect to the processors for the coordinates
    //MPI_File_seek(fh, mp.offsets_elts, MPI_SEEK_CUR );

    // TODO Warning: Element ids are not renumbered, fuzzy behaviours might appear like element
    // getting the same id on two different processes when number of elements tends to be equal
    // to the number of processes

    /* Write element ids */
    std::vector<int32_t> elids;
    for(auto it = lelt_it ; it != lelt_en; ++it )
    {
        auto const& elt = boost::unwrap_ref( *it );
        //LOG(INFO) << std::distance(lelt_it, lelt_en) << " " << mp.offsets_pts << " " << mp.offsets_elts << " " << elt.id() << std::endl;
        elids.push_back(elt.id() + mp.offsets_elts + 1);
    }

    MPI_File_write_ordered(fh, elids.data(), elids.size(), MPI_INT32_T, &status );

    /* Write point ids of vertices */
    /* gather points */
    std::vector<int32_t> ptids;
    for(auto it = lelt_it ; it != lelt_en; ++it )
    {
        auto const& elt = boost::unwrap_ref( *it );
        for ( size_type j = 0; j < elt.numLocalVertices; j++ )
        {
            ptids.push_back(elt.point( j ).id());
        }
    }

    /* translate the point ids using the global id system (where LOCAL and GHOST cells are taken into account) */
    mp.translatePointIds(ptids);

    MPI_File_write_ordered(fh, ptids.data(), ptids.size(), MPI_INT32_T, &status );

    /* Write ghost elements */
    if ( this->worldComm().globalSize() > 1 )
    {
        // get ghost elements
        //auto r1 = markedelements(mesh, part->first, EntityProcessType::GHOST_ONLY );
        auto r1 = markedelements(mesh, markerid, EntityProcessType::GHOST_ONLY );
        auto gelt_it = r1.template get<1>();
        auto gelt_en = r1.template get<2>();

        std::string ghost_t = "g_" + this->elementType();
        // ghosts elements
        memset(buffer, '\0', sizeof(buffer));
        strcpy( buffer, ghost_t.c_str() );
        if( this->worldComm().isMasterRank() )
        { size = sizeof(buffer); }
        else
        { size = 0; }
        MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status );

        Feel::detail::MeshPoints<float> mpg( mesh.get(), this->worldComm(), gelt_it, gelt_en, true, true, true );
        //VLOG(1) << "material : " << p_it->first << " ghost nb element: " << __ne;

        int32_t gnoge = mpg.globalNumberOfElements();
        if( this->worldComm().isMasterRank() )
        { size = 1; }
        else
        { size = 0; }
        MPI_File_write_ordered(fh, &gnoge, size, MPI_INT32_T, &status );

        /* Write elements ids */
        std::vector<int32_t> elids;
        for(auto it = gelt_it ; it != gelt_en; ++it )
        {
            auto const& elt = boost::unwrap_ref( *it );
            elids.push_back(elt.id() + mp.offsets_elts + 1);
        }

        MPI_File_write_ordered(fh, elids.data(), elids.size(), MPI_INT32_T, &status );

        /* Write point ids of vertices */
        std::vector<int32_t> ptids;
        for( auto it = gelt_it; it != gelt_en; ++it )
        {
            auto const& elt = boost::unwrap_ref( *it );
            for ( size_type j = 0; j < elt.numLocalVertices; j++ )
            {
                ptids.push_back(elt.point( j ).id());
            }
        }

        /* translate the point ids using the global id system (where LOCAL and GHOST cells are taken into account) */
        mp.translatePointIds(ptids);

        MPI_File_write_ordered(fh, ptids.data(), ptids.size(), MPI_INT32_T, &status );
    }
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeVariableFiles() const
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
                int dist = 0;

                if( (dist = std::distance(__step->beginNodalScalar(), __step->endNodalScalar())) != 0)
                { LOG(INFO) << "NodalScalar: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginNodalVector(), __step->endNodalVector())) != 0)
                { LOG(INFO) << "NodalVector: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginNodalTensor2(), __step->endNodalTensor2())) != 0)
                { LOG(INFO) << "NodalTensor2: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginElementScalar(), __step->endElementScalar())) != 0)
                { LOG(INFO) << "ElementScalar: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginElementVector(), __step->endElementVector())) != 0)
                { LOG(INFO) << "ElementVector: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginElementTensor2(), __step->endElementTensor2())) != 0)
                { LOG(INFO) << "ElementTensor2: " << dist << std::endl; }

                saveNodal( __step, (__it == __ts->beginStep()), __step->beginNodalScalar(), __step->endNodalScalar() );
                saveNodal( __step, (__it == __ts->beginStep()), __step->beginNodalVector(), __step->endNodalVector() );
                saveNodal( __step, (__it == __ts->beginStep()), __step->beginNodalTensor2(), __step->endNodalTensor2() );

                saveElement( __step, (__it == __ts->beginStep()), __step->beginElementScalar(), __step->endElementScalar() );
                saveElement( __step, (__it == __ts->beginStep()), __step->beginElementVector(), __step->endElementVector() );
                saveElement( __step, (__it == __ts->beginStep()), __step->beginElementTensor2(), __step->endElementTensor2() );

            }

            ++__it;
        }

        ++__ts_it;
    }
}


template<typename MeshType, int N>
template<typename Iterator>
void
ExporterEnsightGold<MeshType,N>::saveNodal( typename timeset_type::step_ptrtype __step, bool isFirstStep, Iterator __var, Iterator en ) const
{
    int size = 0;
    char buffer[ 80 ];

    MPI_Offset offset = 0;
    MPI_File fh;
    MPI_Status status;

    while ( __var != en )
    {
        if ( !__var->second.worldComm().isActive() ) return;

        std::ostringstream __varfname;

        auto __mesh = __step->mesh();

        __varfname << this->path() << "/" << __var->first;
        if( ! boption( _name="exporter.ensightgold.merge.markers") )
        { __varfname << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().localRank(); }
        /* if we want to pack data in several files instead of one */
        /* we compute an index to add to the filename */
        if( boption( _name = "exporter.ensightgold.merge.timesteps")
        && ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
        {
            // timestep indices start at 1
            __varfname << "." << ((__step->index() - TS_INITIAL_INDEX) / ioption( _name="exporter.ensightgold.pack.timesteps" ) + 1);
        }
        // add extension
        if(__var->second.is_scalar)
        { __varfname << ".scl"; }
        else if(__var->second.is_vectorial)
        { __varfname << ".vec"; }
        else if(__var->second.is_tensor2)
        { __varfname << ".tsr"; }
        else
        { __varfname << ".scl"; LOG(ERROR) << "Could not detect data type (scalar, vector, tensor2). Defaulted to scalar." << std::endl; }

        if( ! boption( _name="exporter.ensightgold.merge.timesteps") )
        {
            __varfname << "." << std::setfill( '0' ) << std::setw( M_timeExponent ) << __step->index();
        }
        DVLOG(2) << "[ExporterEnsightGold::saveNodal] saving " << __varfname.str() << "...\n";
        std::fstream __out;

        /* Open File with MPI IO */
        char * str = strdup(__varfname.str().c_str());

        /* Check if file exists if we are on step one and delete it if so */
        /* (MPI IO does not have a truncate mode ) */
        // std::cout << "Nodes " << this->worldComm().isMasterRank() << " " << __step->index() << " " << isFirstStep << " " << fs::exists(str) << std::endl;
        if(this->worldComm().isMasterRank() && isFirstStep && fs::exists(str))
        {
            MPI_File_delete(str, MPI_INFO_NULL);
        }
        MPI_Barrier(this->worldComm().comm());

        MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
        free(str);

        Feel::detail::FileIndex index;

        if( boption( _name = "exporter.ensightgold.merge.timesteps") )
        {
            // first read
            index.read(fh);

            /* Move to the beginning of the fie index section */
            /* to overwrite it */
            if ( index.defined() && (__step->index() - TS_INITIAL_INDEX) > 0 ) {
                MPI_File_seek_shared(fh, index.fileblock_n_steps, MPI_SEEK_SET);
            }
            else {
                // we position the cursor at the beginning of the file
                MPI_File_seek_shared(fh, 0, MPI_SEEK_SET);
            }

            /* Write time step start */
            if( this->worldComm().isMasterRank() )
            { size = sizeof(buffer); }
            else
            { size = 0; }
            memset(buffer, '\0', sizeof(buffer));
            strcpy(buffer,"BEGIN TIME STEP");
            MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

            /* add the beginning of the new block to the file */
            MPI_File_get_position_shared(fh, &offset);
            index.add( offset );
        }

        if( this->worldComm().isMasterRank() )
        { size = sizeof(buffer); }
        else
        { size = 0; }
        memset(buffer, '\0', sizeof(buffer));
        strcpy( buffer, __var->second.name().c_str() );
        MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

        /* handle faces data */
        if ( boption( _name="exporter.ensightgold.save-face" ) )
        {
            BOOST_FOREACH( auto m, __mesh->markerNames() )
            {
                if ( m.second[1] != __mesh->nDim-1 )
                    continue;
                VLOG(1) << "writing face with marker " << m.first << " with id " << m.second[0];
                auto pairit = __mesh->facesWithMarker( m.second[0], this->worldComm().localRank() );
                auto fit = pairit.first;
                auto fen = pairit.second;

                Feel::detail::MeshPoints<float> mp( __mesh.get(), this->worldComm(), fit, fen, true, true, true );
                int __ne = std::distance( fit, fen );

                int nverts = fit->numLocalVertices;
                DVLOG(2) << "Faces : " << __ne << "\n";

                if( this->worldComm().isMasterRank() )
                { size = sizeof(buffer); }
                else
                { size = 0; }
                memset(buffer, '\0', sizeof(buffer));
                strcpy( buffer, "part" );
                MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

                int32_t partid = m.second[0];
                if( this->worldComm().isMasterRank() )
                { size = 1; }
                else
                { size = 0; }
                MPI_File_write_ordered(fh, &partid, size, MPI_INT32_T, &status);

                if( this->worldComm().isMasterRank() )
                { size = sizeof(buffer); }
                else
                { size = 0; }
                memset(buffer, '\0', sizeof(buffer));
                strcpy( buffer, "coordinates" );
                MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

                // write values
                fit = pairit.first;
                fen = pairit.second;

                uint16_type nComponents = __var->second.nComponents;
                if ( __var->second.is_vectorial )
                    nComponents = 3;

                int nfaces = mp.ids.size();
                std::vector<float> field( nComponents*nfaces, 0.0 );
                for( ; fit != fen; ++fit )
                {
                    for ( uint16_type c = 0; c < nComponents; ++c )
                    {
                        for ( size_type j = 0; j < nverts; j++ )
                        {
                            size_type pid = mp.old2new[fit->point( j ).id()]-1;
                            size_type global_node_id = nfaces*c + pid ;
                            if ( c < __var->second.nComponents )
                            {
                                size_type thedof =  __var->second.start() +
                                    boost::get<0>(__var->second.functionSpace()->dof()->faceLocalToGlobal( fit->id(), j, c ));

                                field[global_node_id] = __var->second.globalValue( thedof );
                            }
                            else
                            {
                                field[global_node_id] = 0;
                            }
                        }
                    }
                }
                /* Write each component separately */
                for ( uint16_type c = 0; c < __var->second.nComponents; ++c )
                {
                    MPI_File_write_ordered(fh, field.data() + nfaces * c, nfaces, MPI_FLOAT, &status);
                }
            } // boundaries loop
        }

        /* handle elements */
        for( std::set<int>::iterator mit = M_markersToWrite.begin(); mit != M_markersToWrite.end(); mit++)
        {
            if( this->worldComm().isMasterRank() )
            { size = sizeof(buffer); }
            else
            { size = 0; }
            memset(buffer, '\0', sizeof(buffer));
            strcpy( buffer, "part" );
            MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

            if( this->worldComm().isMasterRank() )
            { size = 1; }
            else
            { size = 0; }

            int32_t partid = *mit;
            MPI_File_write_ordered(fh, &partid, size, MPI_INT32_T, &status);

            if( this->worldComm().isMasterRank() )
            { size = sizeof(buffer); }
            else
            { size = 0; }
            memset(buffer, '\0', sizeof(buffer));
            strcpy( buffer, "coordinates" );
            MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

            uint16_type nComponents = __var->second.nComponents;

            VLOG(1) << "nComponents field: " << nComponents;
            if ( __var->second.is_vectorial )
            {
                nComponents = 3;
                VLOG(1) << "nComponents field(is_vectorial): " << nComponents;
            }

            /* we get that from the local processor */
            /* We do not need the renumbered global index */
            //auto r = markedelements(__mesh,(boost::any)p_it->first,EntityProcessType::ALL);
            auto r = markedelements(__mesh, *mit, EntityProcessType::ALL);
            auto elt_it = r.template get<1>();
            auto elt_en = r.template get<2>();

            Feel::detail::MeshPoints<float> mp( __step->mesh().get(), this->worldComm(), elt_it, elt_en, true, true, true );

            /* create an array to store data per node */
            int npts = mp.ids.size();
            size_type __field_size = npts;
            if ( __var->second.is_vectorial )
                __field_size *= 3;
            ublas::vector<float> __field( __field_size, 0.0 );
            size_type e = 0;
            VLOG(1) << "field size=" << __field_size;
            if ( !__var->second.areGlobalValuesUpdated() )
                __var->second.updateGlobalValues();

            /*
            std::cout << this->worldComm().rank() << " marker=" << *mit << " nbPts:" << npts << " nComp:" << nComponents 
                      << " __evar->second.nComponents:" << __var->second.nComponents << std::endl;
            */

            /* loop on the elements */
            int index = 0;
            for ( ; elt_it != elt_en; ++elt_it )
            {
                VLOG(3) << "is ghost cell " << elt_it->get().isGhostCell();
                /* looop on the ccomponents is outside of the loop on the vertices */
                /* because we need to pack the data in the x1 x2 ... xn y1 y2 ... yn z1 z2 ... zn order */
                for ( uint16_type c = 0; c < nComponents; ++c )
                {
                    for ( uint16_type p = 0; p < __step->mesh()->numLocalVertices(); ++p, ++e )
                    {
                        size_type ptid = mp.old2new[elt_it->get().point( p ).id()]-1;
                        size_type global_node_id = mp.ids.size()*c + ptid ;
                        //LOG(INFO) << elt_it->get().point( p ).id() << " " << ptid << " " << global_node_id << std::endl;
                        DCHECK( ptid < __step->mesh()->numPoints() ) << "Invalid point id " << ptid << " element: " << elt_it->get().id()
                                                                     << " local pt:" << p
                                                                     << " mesh numPoints: " << __step->mesh()->numPoints();
                        DCHECK( global_node_id < __field_size ) << "Invalid dof id : " << global_node_id << " max size : " << __field_size;

                        if ( c < __var->second.nComponents )
                        {
                            size_type dof_id = boost::get<0>( __var->second.functionSpace()->dof()->localToGlobal( elt_it->get().id(), p, c ) );

                            __field[global_node_id] = __var->second.globalValue( dof_id );
                            //__field[npts*c + index] = __var->second.globalValue( dof_id );
                            //DVLOG(3) << "v[" << global_node_id << "]=" << __var->second.globalValue( dof_id ) << "  dof_id:" << dof_id;
                            DVLOG(3) << "v[" << (npts*c + index) << "]=" << __var->second.globalValue( dof_id ) << "  dof_id:" << dof_id;
                        }
                        else
                        {
                            __field[global_node_id] = 0.0;
                            //__field[npts*c + index] = 0.0;
                        }
                    }
                }

                /* increment index of vertex */
                index++;
            }

            /* Write each component separately */
            for ( uint16_type c = 0; c < nComponents; ++c )
            {
                MPI_File_write_ordered(fh, ((float *)(__field.data().begin())) + npts * c, npts, MPI_FLOAT, &status);
            }
        } // parts loop

        if( boption(_name="exporter.ensightgold.merge.timesteps") )
        {
            /* write timestep end */
            if( this->worldComm().isMasterRank() )
            { size = sizeof(buffer); }
            else
            { size = 0; }
            memset(buffer, '\0', sizeof(buffer));
            strcpy(buffer,"END TIME STEP");
            MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

            // write back the file index
            index.write( fh );
        }
        DVLOG(2) << "[ExporterEnsightGold::saveNodal] saving " << __varfname.str() << "done\n";

        MPI_File_close(&fh);

        ++__var;
    }
}

template<typename MeshType, int N>
template<typename Iterator>
void
ExporterEnsightGold<MeshType,N>::saveElement( typename timeset_type::step_ptrtype __step, bool isFirstStep, Iterator __evar, Iterator __evaren ) const
{
    int size;
    char buffer[ 80 ];

    MPI_File fh;
    MPI_Offset offset = 0;
    MPI_Status status;

    while ( __evar != __evaren )
    {
        if ( !__evar->second.worldComm().isActive() ) return;

        std::ostringstream __evarfname;

        __evarfname << this->path() << "/" << __evar->first;

        if( ! boption( _name="exporter.ensightgold.merge.markers") )
        { __evarfname << "-" << this->worldCommBase().globalSize() << "_" << this->worldCommBase().localRank(); }
        /* if we want to pack data in several files instead of one */
        /* we compute an index to add to the filename */
        if( boption( _name = "exporter.ensightgold.merge.timesteps")
        && ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
        {
            // timestep indices start at 1
            __evarfname << "." << ((__step->index() - TS_INITIAL_INDEX) / ioption( _name="exporter.ensightgold.pack.timesteps" ) + 1);
        }
        // add extension
        if(__evar->second.is_scalar)
        { __evarfname << ".scl"; }
        else if(__evar->second.is_vectorial)
        { __evarfname << ".vec"; }
        else if(__evar->second.is_tensor2)
        { __evarfname << ".tsr"; }
        else
        { __evarfname << ".scl"; LOG(ERROR) << "Could not detect data type (scalar, vector, tensor2). Defaulted to scalar." << std::endl; }

        if(! boption( _name="exporter.ensightgold.merge.timesteps") )
        {
            __evarfname << "." << std::setfill( '0' ) << std::setw( M_timeExponent ) << __step->index();
        }
        DVLOG(2) << "[ExporterEnsightGold::saveElement] saving " << __evarfname.str() << "...\n";

        auto __mesh = __step->mesh();

        /* Open File with MPI IO */
        char * str = strdup(__evarfname.str().c_str());

        /* Check if file exists if we are on step one and delete it if so */
        /* (MPI IO does not have a truncate mode ) */
        // std::cout << "Elt: " << this->worldComm().isMasterRank() << " " << __step->index() << " " << isFirstStep << " " << fs::exists(str) << std::endl;
        if(this->worldComm().isMasterRank() && isFirstStep && fs::exists(str))
        {
            MPI_File_delete(str, MPI_INFO_NULL);
        }
        MPI_Barrier(this->worldComm().comm());

        MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
        free(str);

        Feel::detail::FileIndex index;

        if( boption( _name = "exporter.ensightgold.merge.timesteps") )
        {
            // first read
            index.read(fh);

            /* Move to the beginning of the fie index section */
            /* to overwrite it */
            if ( index.defined() && (__step->index() - TS_INITIAL_INDEX) > 0 ) {
                MPI_File_seek_shared(fh, index.fileblock_n_steps, MPI_SEEK_SET);
            }
            else {
                // we position the cursor at the beginning of the file
                MPI_File_seek_shared(fh, 0, MPI_SEEK_SET);
            }

            /* Write time step start */
            if( this->worldComm().isMasterRank() )
            { size = sizeof(buffer); }
            else
            { size = 0; }
            memset(buffer, '\0', sizeof(buffer));
            strcpy(buffer,"BEGIN TIME STEP");
            MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

            /* add the beginning of the new block to the file */
            MPI_File_get_position_shared(fh, &offset);
            index.add( offset );
        }

        if( this->worldComm().isMasterRank() )
        { size = sizeof(buffer); }
        else
        { size = 0; }
        memset(buffer, '\0', sizeof(buffer));
        strcpy( buffer, __evar->second.name().c_str() );
        MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

        // iterate over the markers
        for( std::set<int>::iterator mit = M_markersToWrite.begin() ; mit != M_markersToWrite.end(); mit++ )
        {
            if( this->worldComm().isMasterRank() )
            { size = sizeof(buffer); }
            else
            { size = 0; }
            memset(buffer, '\0', sizeof(buffer));
            strcpy( buffer, "part" );
            MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

            int32_t partid = *mit;
            if( this->worldComm().isMasterRank() )
            { size = 1; }
            else
            { size = 0; }
            MPI_File_write_ordered(fh, &partid, size, MPI_INT32_T, &status);
            DVLOG(2) << "part " << buffer << "\n";

            if( this->worldComm().isMasterRank() )
            { size = sizeof(buffer); }
            else
            { size = 0; }
            memset(buffer, '\0', sizeof(buffer));
            strcpy( buffer, this->elementType().c_str() );
            MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
            DVLOG(2) << "element type " << buffer << "\n";

            uint16_type nComponents = __evar->second.nComponents;

            if ( __evar->second.is_vectorial )
                nComponents = 3;

            size_type __field_size = nComponents * __evar->second.size()/__evar->second.nComponents;
            ublas::vector<float> __field( __field_size );
            __field.clear();

            auto r = markedelements(__step->mesh(), *mit, EntityProcessType::LOCAL_ONLY );
            auto elt_st = r.template get<1>();
            auto elt_en = r.template get<2>();

            if ( !__evar->second.areGlobalValuesUpdated() )
                __evar->second.updateGlobalValues();

            DVLOG(2) << "[saveElement] firstLocalIndex = " << __evar->second.firstLocalIndex() << "\n";
            DVLOG(2) << "[saveElement] lastLocalIndex = " << __evar->second.lastLocalIndex() << "\n";
            DVLOG(2) << "[saveElement] field.size = " << __field_size << "\n";

            //size_type ncells = __evar->second.size()/__evar->second.nComponents;
            size_type ncells = std::distance( elt_st, elt_en );

            /*
            std::cout << this->worldComm().rank() << " marker=" << *mit << " nbElts:" << ncells << " nComp:" << nComponents 
                      << " __evar->second.nComponents:" << __evar->second.nComponents << std::endl;
            */

            for ( int c = 0; c < nComponents; ++c )
            {
                size_type e = 0;
                for ( auto elt_it = elt_st ; elt_it != elt_en; ++elt_it, ++e )
                {
                    auto const& elt = boost::unwrap_ref( *elt_it );
                    DVLOG(2) << "pid : " << this->worldComm().globalRank()
                             << " elt_it :  " << elt.id()
                             << " e : " << e << "\n";

                    size_type global_node_id = c*ncells+e ;

                    if ( c < __evar->second.nComponents )
                    {
                        size_type dof_id = boost::get<0>( __evar->second.functionSpace()->dof()->localToGlobal( elt.id(),0, c ) );

                        DVLOG(2) << "c : " << c
                                 << " gdofid: " << global_node_id
                                 << " dofid : " << dof_id
                                 << " f.size : " <<  __field.size()
                                 << " e.size : " <<  __evar->second.size()
                                 << "\n";

                        __field[global_node_id] = __evar->second.globalValue( dof_id );

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

            /* Write each component separately */
            for ( uint16_type c = 0; c < nComponents; ++c )
            {
                MPI_File_write_ordered(fh, __field.data().begin() + ncells * c, ncells, MPI_FLOAT, &status);
            }
        }

        if( boption(_name="exporter.ensightgold.merge.timesteps") )
        {
            /* write timestep end */
            if( this->worldComm().isMasterRank() )
            { size = sizeof(buffer); }
            else
            { size = 0; }
            memset(buffer, '\0', sizeof(buffer));
            strcpy(buffer,"END TIME STEP");
            MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);

            // write back the file index
            index.write( fh );
        }

        MPI_File_close(&fh);

        DVLOG(2) << "[ExporterEnsightGold::saveElement] saving " << __evarfname.str() << "done\n";
        ++__evar;
    }
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::visit( mesh_type* __mesh )
{
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
#endif // FEELPP_EXPORTERENSIGHT_CPP
