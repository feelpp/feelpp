/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( worldcomm_ptr_t const& worldComm )
:
super( worldComm ),
M_element_type()
{
    init();
}
template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( std::string const& __p, int freq, worldcomm_ptr_t const& worldComm )
    :
    super( "ensightgold", __p, freq, worldComm ),
    M_element_type()
{
    init();
}
template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( po::variables_map const& vm, std::string const& exp_prefix, worldcomm_ptr_t const& worldComm )
    :
    super( vm, exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( std::string const& exp_prefix, worldcomm_ptr_t const& worldComm )
    :
    super( exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( ExporterEnsightGold const & __ex )
    :
    super( __ex ),
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
    tic();
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

    tic();
    writeGeoFiles();
    toc("ExporterEnsightGold::save geo",FLAGS_v>1);


    /* only try to write variable data when we have time steps */
    if(hasSteps)
    {
        tic();
        writeVariableFiles();
        toc("ExporterEnsightGold::save variables",FLAGS_v>1);
    }
    tic();
    this->saveTimeSet();
    toc("ExporterEnsightGold::save timeset",FLAGS_v>1);

    tic();
    writeCaseFile();
    toc("ExporterEnsightGold::save case",FLAGS_v>1);

    tic();
    writeSoSFile();
    toc("ExporterEnsightGold::save sos",FLAGS_v>1);
    toc("ExporterEnsightGold::save", FLAGS_v > 0 );
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
                << "cfiles pattern: "<< this->prefix() << ".case\n"
                << "cfiles start number: 0\n"
                << "cfiles increment: 1\n\n";
        }
        else
        {
            __out << "MULTIPLE_CASEFILES\n"
                << "total number of cfiles: " << this->worldComm().globalSize() << "\n"
                << "cfiles global path: " << fs::current_path().string() << "\n"
                << "cfiles pattern: "<<this->prefix() << "-" << this->worldComm().globalSize() << "_*.case\n"
                << "cfiles start number: 0\n"
                << "cfiles increment: 1\n\n";
        }
        __out << "SERVERS\n"
              << "number of servers: "<< (this->worldComm().globalSize()/100)+1 <<" repeat\n";

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
            __outparaview << "number of servers: " << this->worldComm().globalSize() << "\n";

            for ( int pid = 0 ; pid < this->worldComm().globalSize(); ++pid )
            {
                __outparaview << "#Server " << pid+1 << "\n"
                    << "machine id: " << mpi::environment::processor_name() << "\n"
                    << "executable: /usr/local/bin/ensight76/bin/ensight7.server\n"
                    << "data_path: " << fs::current_path().string() << "\n"
                    << "casefile: " << this->prefix() << "-" << this->worldComm().globalSize() << "_" << pid << ".case\n";
            }
        }
    }
}
template<typename MeshType, int N>
template<typename Iterator,typename TSt>
void
ExporterEnsightGold<MeshType,N>::writeCaseFileVariables( Iterator it, Iterator end,
                                                         std::string const& loc,
                                                         std::string const& /*type*/,
                                                         std::string const& /*ext*/,
                                                         TSt const& __ts,
                                                         std::ostream& __out ) const
{
    for( ; it != end; ++it )
    {
        std::string type, ext;
        if ( it->second.first == FunctionSpaceType::SCALAR )
        {
            type = "scalar";
            ext = "scl";
        }
        else if ( it->second.first == FunctionSpaceType::VECTORIAL )
        {
            type = "vector";
            ext = "vec";
        }
        else if ( it->second.first == FunctionSpaceType::TENSOR2 )
        {
            type = "tensor asym";
            ext = "tsr";
        }
        else if ( it->second.first == FunctionSpaceType::TENSOR2_SYMM )
        {
            type = "tensor symm";
            ext = "tsrs";
        }

        if( boption( _name="exporter.ensightgold.merge.timesteps") )
        {
            __out << type << " per " << loc << ": "
                  << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                  << unwrap_ptr(it->second.second[0][0]).name() << " " << __ts->name() << "." << it->first;
            if( ! boption( _name="exporter.ensightgold.merge.markers") )
            {
                __out << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank();
            }
            /* if we want to pack data in several files instead of one, we add an index to the filename */
            if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
            { __out << ".*"; }
            __out << "." << ext << std::endl;
        }
        else
        {
            __out << type << " per " << loc << ": "
                  << __ts->index() << " " // << *__ts_it->beginStep() << " "
                  << unwrap_ptr(it->second.second[0][0]).name() << " " << __ts->name() << "." << it->first;
            if( ! boption( _name="exporter.ensightgold.merge.markers") )
            {
                __out << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank();
            }
            __out << "." << ext << "." << std::string(M_timeExponent, '*') << std::endl;
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

        while ( __ts_it != __ts_en )
        {
            timeset_ptrtype __ts = *__ts_it;

            filestr << this->path() << "/"
                << __ts->name();
            if( ! boption( _name="exporter.ensightgold.merge.markers") )
            { filestr << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank(); }
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
                    __out << "model: " << __ts->name();
                    if( ! boption( _name="exporter.ensightgold.merge.markers") )
                    { __out << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank(); }
                    __out << ".geo";
                }
                break;
                default:
                case EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY:
                case EXPORTER_GEOMETRY_CHANGE:
                {
                    if( boption( _name="exporter.ensightgold.merge.timesteps") )
                    {
                        __out << "model: " << __ts->index() << " 1 " << __ts->name();
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank(); }
                        /* if we want to pack data in several files instead of one, we add an index to the filename */
                        if( ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
                        { __out << ".*"; }
                        __out << ".geo";
                    }
                    else
                    {
                        __out << "model: " << __ts->index() << " " << __ts->name();
                        if( ! boption( _name="exporter.ensightgold.merge.markers") )
                        { __out << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank(); }
                        __out << ".geo" << "." << std::string(M_timeExponent, '*');
                    }

                    if ( this->exporterGeometry() == EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY )
                    {
                        __out << " change_coords_only";
                    }
                }
                break;
            }
            __out << "\n";

            __out << "VARIABLES:" << "\n";

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
                        __out << "constant per case: " << s_it->first << " " << __ts->name() << "." << s_it->second.first << "\n";
                    }
                    else
                    {
                        if ( this->worldComm().isMasterRank() )
                        {

                            __out << "constant per case file: " << __ts->index() << " " << s_it->first << " " << __ts->name() << "." << s_it->first << ".cst";
                            // loop over time
                            auto stepit = __ts->beginStep();
                            auto stepen = __ts->endStep();

                            std::ofstream ofs;
                            int d = std::distance( stepit, stepen );
                            LOG(INFO) << "distance = " << d;
                            if ( d > 1 )
                                ofs.open( this->path()+ "/"+ __ts->name()+"."+s_it->first+".cst", std::ios::out | std::ios::app );
                            else
                                ofs.open( this->path()+ "/" +__ts->name()+"."+s_it->first+".cst", std::ios::out );

                            auto step = *boost::prior(stepen);
                            ofs << step->scalar( s_it->first ) << "\n";
                            ofs.close();
                            __out << "\n";
                        }

                    }
                }

                writeCaseFileVariables( ( *__tstp_it )->beginNodal(),
                                        ( *__tstp_it )->endNodal(),
                                        "node", "scalar", "scl",
                                        __ts, __out );
                writeCaseFileVariables( ( *__tstp_it )->beginElement(),
                                        ( *__tstp_it )->endElement(),
                                        "element",  "scalar", "scl",
                                        __ts, __out );

#if 0
                writeCaseFileVariables( ( *__tstp_it )->beginNodalVector(),
                                        ( *__tstp_it )->endNodalVector(),
                                        "node", "vector", "vec",
                                        __ts, __out );
                writeCaseFileVariables( ( *__tstp_it )->beginNodalTensor2(),
                                        ( *__tstp_it )->endNodalTensor2(),
                                        "node", "tensor asym", "tsr",
                                        __ts, __out );
                writeCaseFileVariables( ( *__tstp_it )->beginNodalTensor2Symm(),
                                        ( *__tstp_it )->endNodalTensor2Symm(),
                                        "node", "tensor symm", "tsrs",
                                        __ts, __out );
                writeCaseFileVariables( ( *__tstp_it )->beginElementScalar(),
                                        ( *__tstp_it )->endElementScalar(),
                                        "element",  "scalar", "scl",
                                        __ts, __out );
                writeCaseFileVariables( ( *__tstp_it )->beginElementVector(),
                                        ( *__tstp_it )->endElementVector(),
                                        "element", "vector", "vec",
                                        __ts, __out );
                writeCaseFileVariables( ( *__tstp_it )->beginElementTensor2(),
                                        ( *__tstp_it )->endElementTensor2(),
                                        "element", "tensor asym", "tsr",
                                        __ts, __out );
                writeCaseFileVariables( ( *__tstp_it )->beginElementTensor2Symm(),
                                        ( *__tstp_it )->endElementTensor2Symm(),
                                        "element", "tensor symm", "tsrs",
                                        __ts, __out );
#endif

            }

            __out << "TIME:\n";

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

            uint16_type __l = 1;

            while ( __its != __ens )
            {

                __out << std::scientific << std::setprecision( 6 ) << ( *__its )->time() << " ";

                if ( __l++ % 10 == 0 )
                    __out << "\n";

                ++__its;
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
                        filestr << __ts->name() << "-"
                            << this->worldComm().globalSize() << "_" << p << ".case";
                        __out << filestr.str() << "\n        ";
                    }
                }
            } // use-sos

            __out.close();
            ++__ts_it;
        }
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
    //MPI_Info info;
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
                { __geofname << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank(); }
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

                //init MPI_Info object fromhints defined as environment variables
                //MPI_Info info = initIoInfoFromEnvVars();
                if( boption( _name="exporter.ensightgold.merge.timesteps" ) )
                    MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL , &fh );
                else
                    MPI_File_open( this->worldComm().comm(), str, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL , &fh );

                //MPI_Info_free(&info);
                free(str);
                // Initializing cursor in file
                posInFile=0;


                Feel::detail::FileIndex index;

                if( boption( _name="exporter.ensightgold.merge.timesteps" ) )
                {
                    // first read
                    index.read(fh);
                    // we position the cursor at the beginning of the file
                    // MPI_File_seek_shared(fh, 0, MPI_SEEK_SET);
                    posInFile=0;

                    /* write C binary if we didn't find the index <=> first pass on the file */
                    if( !index.defined() )
                    {
                        if( this->worldComm().isMasterRank() )
                        {
                            memset(buffer, '\0', sizeof(buffer));
                            strncpy(buffer, "C Binary", 80);
                            //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                            MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                        }
                        // all procs increment their offset
                        posInFile+=80;
                    }


                    /* Write time step start */
                    if( this->worldComm().isMasterRank() )
                    {
                        memset(buffer, '\0', sizeof(buffer));
                        strncpy(buffer, "BEGIN TIME STEP", 80);
                        //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                        MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                        LOG(INFO) << "saveNodal out: " << buffer;
                    }
                    // all procs increment their offset
                    posInFile+=80;

                    /* add the beginning of the new block to the file */
                    //MPI_File_get_position_shared(fh, &offset);
                    index.add( posInFile );
                }

                /* Write the file */
                this->writeGeoMarkers(fh, mesh);
                if( boption( _name="exporter.ensightgold.merge.timesteps" ) )
                {
                    /* write timestep end */
                    if( this->worldComm().isMasterRank() )
                    {
                        char bufferEnd[80];
                        memset(buffer, '\0', 80);
                        strcpy(buffer,"END TIME STEP");
                        //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                        MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                    }
                    posInFile+=80;
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
                { __geofname << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank(); }
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

                //init MPI_Info object from hints defined as environment variables
                //MPI_Info info = initIoInfoFromEnvVars();

                MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );

                //MPI_Info_free(&info);
                free(str);

                Feel::detail::FileIndex index;

                // first read
                index.read(fh);

                /* Move to the beginning of the fie index section */
                /* to overwrite it */
                if ( index.defined() && (timeIndex - TS_INITIAL_INDEX) > 0 ) {
                    // MPI_File_seek_shared(fh, index.fileblock_n_steps, MPI_SEEK_SET);
                    posInFile=index.fileblock_n_steps;
                }
                else {
                    // we position the cursor at the beginning of the file
                    // MPI_File_seek_shared(fh, 0, MPI_SEEK_SET);
                    posInFile=0;
                }

                /* write C binary if we didn't find the index <=> first pass on the file */
                if( !index.defined() )
                {
                    if( this->worldComm().isMasterRank() ){
                        memset(buffer, '\0', sizeof(buffer));
                        strncpy(buffer, "C Binary", 80);
                        MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                    }
                    posInFile+=80;
                }

                /* Write time step start */
                if( this->worldComm().isMasterRank() ){
                    memset(buffer, '\0', sizeof(buffer));
                    strncpy(buffer, "BEGIN TIME STEP", 80);
                    //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                    MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                    LOG(INFO) << "saveNodal out: " << buffer;
                }
                posInFile+=80;

                /* add the beginning of the new block to the file */
                index.add( posInFile );

                /* write data for timestep */
                this->writeGeoMarkers(fh, mesh);

                /* write timestep end */
                if( this->worldComm().isMasterRank() )
                {
                    memset(buffer, '\0', sizeof(buffer));
                    strncpy(buffer, "END TIME STEP", 80);
                    //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                    MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                }
                posInFile+=80;
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
                { __geofname << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank(); }
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

                //init MPI_Info object fromhints defined as environment variables
                //MPI_Info info = initIoInfoFromEnvVars();
                MPI_File_open( this->worldComm().comm(), str, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
                //MPI_Info_free(&info);

                free(str);

                // Initializing cursor in file
                posInFile=0;


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
            auto rangeMarkedElements = mesh->elementsWithMarker( *it );
            osss << " " << *it << " (" << std::distance( std::get<0>( rangeMarkedElements ),std::get<1>( rangeMarkedElements ) ) << ")";
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
        for( auto & m : mesh->markerNames() )
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

    DVLOG(2) << "Merging markers : " << "\n";

    /* write header */
    /* little trick to only perform collective operation (optimized) */
    /* and avoid scattering offset and reseting the shared pointer if we would write this only on master proc */

    // only write C Binary if we are not mergin timesteps
    // as it is oalready writtent at the beginning of the file
    if( ! boption( _name="exporter.ensightgold.merge.timesteps") )
    {
        char buffer[80];
        if( this->worldComm().isMasterRank() )
        {
            memset(buffer, '\0', sizeof(buffer));
            strncpy(buffer, "C Binary", 80);
            MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
            //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
        }
        posInFile+=80;
        //LOG(INFO) << "wrote " << buffer << std::endl;
    }

    char buffer[320];
    if( this->worldComm().isMasterRank() )
    {
        // get only the filename (maybe with full path)
        fs::path gp = M_filename;
        std::string theFileName = gp.filename().string();
        CHECK( theFileName.length() < 80 ) << "the file name is too long : theFileName=" << theFileName << "\n";
        char buffer2[320];
        memset(buffer2, '\0', sizeof(buffer2));
        strcpy( buffer2, theFileName.c_str() );
        strcpy( &buffer2[80], "elements" );
        strcpy( &buffer2[160], "node id given" );
        strcpy( &buffer2[240], "element id given" );
        //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
        MPI_File_write_at(fh, posInFile, buffer2, sizeof(buffer2), MPI_CHAR, &status);
    }
    posInFile+=320;

}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeGeoMarkedFaces(MPI_File fh, mesh_ptrtype mesh, std::pair<const std::string, std::vector<index_type> > & m) const
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
    auto rangeMarkedFaces = mesh->facesWithMarker( m.second[0], this->worldComm().localRank() );
    auto fit = std::get<0>( rangeMarkedFaces );
    auto fen = std::get<1>( rangeMarkedFaces );
    Feel::detail::MeshPoints<float> mp( mesh.get(), this->worldComm(), fit, fen, true, true, true );
    int32_t __ne = std::distance( fit, fen );
    if ( __ne == 0 ) return;
    int nverts = boost::unwrap_ref( *fit ).numLocalVertices;
    DVLOG(2) << "Faces : " << __ne << "\n";

    int sizeOfInt32_t;
    int sizeOfFloat;
    MPI_Type_size( MPI_INT32_T , &sizeOfInt32_t );
    MPI_Type_size( MPI_FLOAT , &sizeOfFloat );
    int localOffset, sumOffsets; //needed to calculate offsets for individual mpiio cursors

    if( this->worldComm().isMasterRank() )
    {
        memset(buffer, '\0', sizeof(buffer));
        strncpy(buffer, "part", 80);
        //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
        MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
    }
    posInFile+=80;

    if( this->worldComm().isMasterRank() )
    {
        int32_t partid = m.second[0] + 1;
        MPI_File_write_at(fh, posInFile, &partid, 1, MPI_INT32_T, &status);
        // MPI_File_write_ordered(fh, &partid, 1, MPI_INT32_T, &status);
    }
    posInFile += sizeOfInt32_t;

    if( this->worldComm().isMasterRank() )
    {
        char buffer2[160];
        memset(buffer2, '\0', sizeof(buffer2));
        strncpy(buffer2, m.first.c_str(), 80 - 1 );
        strcpy(&buffer2[80], "coordinates");
        // MPI_File_write_ordered(fh, buffer2, size, MPI_CHAR, &status);
        MPI_File_write_at(fh, posInFile, buffer2, sizeof(buffer2), MPI_CHAR, &status);
    }
    posInFile+=160;

    // write points coordinates
    // fit = pairit.first;
    // fen = pairit.second;
    fit = std::get<2>( rangeMarkedFaces )->begin();

    size_type __nv = mp.ids.size();
    int32_t gnop = (int32_t)(mp.globalNumberOfPoints());

    if( this->worldComm().isMasterRank() )
    {
        // write number of points
        // MPI_File_write_ordered(fh, &gnop, size, MPI_INT32_T, &status);
        MPI_File_write_at(fh, posInFile, &gnop, 1, MPI_INT32_T, &status);
    }
    posInFile += sizeOfInt32_t;

    /* write points ids */
    tic();
    // all procs writing :
    // - calculate mp.ids.size()*sizeOfInt32_t
    int ptIdWritingSize = mp.ids.size()*sizeOfInt32_t;
    localOffset = 0;
    // - calculate every proc offset "localOffset" with Mpi_Exscan(...)
    MPI_Exscan(&ptIdWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
    // - write in file on cursor : posInFile + localOffset
    MPI_File_write_at(fh, posInFile+localOffset, mp.ids.data(), mp.ids.size(),MPI_INT32_T, &status);
    // - calculate the whole offset to increment posInFile :
    sumOffsets = localOffset + ptIdWritingSize;
    MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
    posInFile += sumOffsets;
    toc("ExporterEnsightGold writeVariableFiles write ids",FLAGS_v>0);

    /* write points coordinates in the order x1 ... xn y1 ... yn z1 ... zn */
    tic();
    // All procs write :
    // - calculate every proc size of part to write
    int coordsWritingSize = __nv*sizeOfFloat;
    localOffset = 0;
    // - calculate every proc offset "localOffset" with Mpi_Exscan
    MPI_Exscan(&coordsWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
    sumOffsets = localOffset + coordsWritingSize;
    MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
    // - write in file on cursor : posInFile + localOffset + i*sumOffsets
    for(int i = 0; i < 3; i++)
    {
        MPI_File_write_at(fh, posInFile + i*sumOffsets + localOffset, \
                          mp.coords.data() + i * __nv, __nv, MPI_FLOAT, &status );
        //MPI_File_write_ordered(fh, mp.coords.data() + i * __nv, __nv, MPI_FLOAT, &status );
    }
    posInFile += 3*sumOffsets;
    toc("ExporterEnsightGold writeVariableFiles write coords",FLAGS_v>0);

    // write connectivity
    // fit = pairit.first;
    // fen = pairit.second;
    fit = std::get<2>( rangeMarkedFaces )->begin();

    if( this->worldComm().isMasterRank() )
    {
        memset(buffer, '\0', sizeof(buffer));
        strcpy( buffer, M_face_type.c_str() );
        MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
        //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
        VLOG(1) << "face type " << buffer;
    }
    posInFile +=80;


    if( this->worldComm().isMasterRank() )
    {
        MPI_File_write_at(fh, posInFile, &__ne, 1, MPI_INT32_T, &status);
        //MPI_File_write_ordered(fh, &__ne, size, MPI_INT32_T, &status);
        VLOG(1) << "n faces " << __ne;
    }
    posInFile += sizeOfInt32_t;

    idnode.resize( __ne );
    // fit = pairit.first;
    fit = std::get<2>( rangeMarkedFaces )->begin();
    size_type e = 0;
    for ( ; fit != fen; ++fit, ++e )
    {
        auto const& face = boost::unwrap_ref( *fit );
        idnode[e] = face.id() + 1;
    }
    CHECK( e == idnode.size() ) << "Invalid number of face id for part " << m.first;

    // All procs write :
    // - calculate every proc size of part to write
    int idnodeWritingSize = idnode.size()*sizeOfInt32_t;
    localOffset = 0;
    // - calculate every proc offset "localOffset" with Mpi_Exscan
    MPI_Exscan(&idnodeWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
    // - write in file on cursor : posInFile + localOffset
    MPI_File_write_at(fh, posInFile+localOffset, idnode.data(), idnode.size(), MPI_INT32_T, &status);
    //MPI_File_write_ordered(fh, idnode.data(), idnode.size(), MPI_INT32_T, &status);

    sumOffsets = localOffset + idnodeWritingSize;
    MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
    posInFile += sumOffsets; // all procs go to the end of last writing.



    idelem.resize( __ne*nverts );
    // fit = pairit.first;
    fit = std::get<2>( rangeMarkedFaces )->begin();
    e = 0;
    for( ; fit != fen; ++fit, ++e )
    {
        auto const& face = boost::unwrap_ref( *fit );
        for ( size_type j = 0; j < nverts; j++ )
        {
            // ensight id start at 1
            idelem[e*nverts+j] = mp.old2new[face.point( j ).id()];
        }
    }
    CHECK( e*nverts == idelem.size() ) << "Invalid number of faces " << e*nverts << " != " << idelem.size() << " in connectivity for part " << m.first;


    // All procs write :
    // - calculate every proc size of part to write
    int idelemWritingSize = idelem.size()*sizeOfInt32_t;
    localOffset = 0;
    // - calculate every proc offset "localOffset" with Mpi_Exscan
    MPI_Exscan(&idelemWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
    // - write in file on cursor : posInFile + localOffset

    MPI_File_write_at(fh, posInFile+localOffset, idelem.data(), idelem.size(), MPI_INT32_T, &status);
    //MPI_File_write_ordered(fh, idelem.data(), idelem.size(), MPI_INT32_T, &status);

    sumOffsets = localOffset + idelemWritingSize;
    MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
    posInFile += sumOffsets; // all procs go to the end of last writing.
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeGeoMarkedElements(MPI_File fh, mesh_ptrtype mesh, size_type markerid) const
{
    MPI_Status status;

    int size = 0;
    char buffer[80];

    int sizeOfInt32_t;
    int sizeOfFloat;
    MPI_Type_size( MPI_INT32_T , &sizeOfInt32_t );
    MPI_Type_size( MPI_FLOAT , &sizeOfFloat );
    int localOffset, sumOffsets;

    std::vector<int32_t> idnode, idelem;

    //auto r = markedelements(mesh, part->first, EntityProcessType::ALL );
    auto r = markedelements(mesh, markerid, EntityProcessType::LOCAL_ONLY );
    auto allelt_it = r.template get<1>();
    auto allelt_en = r.template get<2>();

    //VLOG(1) << "material : " << m << " total nb element: " << std::distance(allelt_it, allelt_en );
    //VLOG(1) << "material : " << m << " ghost nb element: " << std::distance(gelt_it, gelt_en );
    //VLOG(1) << "material : " << m << " local nb element: " << std::distance(lelt_it, lelt_en );
    M_cache_mp.try_emplace( markerid, mesh.get(), this->worldComm(), allelt_it, allelt_en, true, true, true );
    auto& mp = M_cache_mp.at( markerid );
    
    VLOG(1) << "mesh pts size : " << mp.ids.size();
    // part

    if( this->worldComm().isMasterRank() )
    {
        memset(buffer, '\0', sizeof(buffer));
        strncpy(buffer, "part", 80);
        // write number of points
        MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
        //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
    }
    posInFile+=80;
    // Was previously using p_it->first as partid
    // part id
    if ( this->worldComm().isMasterRank() )
    {
        int32_t partid = markerid + 1;
        LOG(INFO) << "writing part " << partid << std::endl;
        // write number of points
        MPI_File_write_at(fh, posInFile, &partid, 1, MPI_INT32_T, &status);
        // MPI_File_write_ordered(fh, &partid, size, MPI_INT32_T, &status);
    }
    posInFile+=sizeOfInt32_t;

    if( this->worldComm().isMasterRank() )
    {
        char buffer2[160];
        // material
        memset(buffer2, '\0', sizeof(buffer2));
        //sprintf(buffer, "Material %d", part->first);
        // Restrict the marker name to 32 chars to avoid buffer overflow (48 chars left for text + id)
        sprintf(buffer2, "Marker %d (%s)", (int)(markerid), mesh->markerName(markerid).substr(0, 32).c_str());
        strcpy( &buffer2[80], "coordinates" );

        MPI_File_write_at(fh, posInFile, buffer2, sizeof(buffer2), MPI_CHAR, &status);
        //MPI_File_write_ordered(fh, buffer2, size, MPI_CHAR, &status);
    }
    posInFile+=160;

    size_type __nv = mp.ids.size();
    int32_t gnop = (int32_t)(mp.globalNumberOfPoints());

    if( this->worldComm().isMasterRank() )
    {
        MPI_File_write_at(fh, posInFile, &gnop, 1, MPI_INT32_T, &status );
        //MPI_File_write_ordered(fh, &gnop, size, MPI_INT32_T, &status );
    }
    posInFile+=sizeOfInt32_t;

    // TODO modify this code !
    // Integrate id translation in the MeshPoints constructor
    /* write points ids */
    std::vector<int32_t> pointids;
    for(int i = 0 ; i < mp.ids.size() ; ++i )
    {
        pointids.push_back(mp.ids.at(i) + mp.offsets_pts);
    }

    // All procs write :
    // - calculate every proc size of part to write
    int pointidsWritingSize = pointids.size()*sizeOfInt32_t;
    localOffset = 0;
    // - calculate every proc offset "localOffset" with Mpi_Exscan
    MPI_Exscan(&pointidsWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
    // - write in file on cursor : posInFile + localOffset

    MPI_File_write_at(fh, posInFile+localOffset, pointids.data(), pointids.size(), MPI_INT32_T, &status );
    //MPI_File_write_ordered(fh, pointids.data(), pointids.size(), MPI_INT32_T, &status );

    sumOffsets = localOffset + pointidsWritingSize;
    // last proc bcasts to all
    MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
    posInFile += sumOffsets; // all procs go to the end of last writing.

    /* write points coordinates in the order x1 ... xn y1 ... yn z1 ... zn */

    // All procs write :
    // - calculate every proc size of part to write
    int coordsWritingSize = __nv*sizeOfFloat;
    localOffset = 0;
// - calculate every proc offset "localOffset" with Mpi_Exscan
    MPI_Exscan(&coordsWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());

    sumOffsets = localOffset + coordsWritingSize;
    MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());

    // - write in file on cursor : posInFile + localOffset
    for(int i = 0; i < 3; i++)
    {
        MPI_File_write_at(fh, posInFile + localOffset + i*sumOffsets, \
                          mp.coords.data() + i * __nv, __nv, MPI_FLOAT, &status );
        //MPI_File_write_ordered(fh, mp.coords.data() + i * __nv, __nv, MPI_FLOAT, &status );
    }
    posInFile += 3*sumOffsets;


    /* write element type */
    if( this->worldComm().isMasterRank() )
    {
        memset(buffer, '\0', sizeof(buffer));
        strncpy( buffer, this->elementType().c_str(), 80);
        MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
        //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status );
    }
    posInFile+=80;

    /* compute and write number of elements */
    //auto r2 = markedelements(mesh, part->first, EntityProcessType::LOCAL_ONLY );
    auto r2 = markedelements(mesh, markerid, EntityProcessType::LOCAL_ONLY );
    auto lelt_it = r2.template get<1>();
    auto lelt_en = r2.template get<2>();

    M_cache_mp.try_emplace( markerid, mesh.get(), this->worldComm(), lelt_it, lelt_en, true, true, true );
    mp = M_cache_mp.at( markerid );
    
    int32_t gnole = mp.globalNumberOfElements();
    //LOG(INFO) << "Global nb elements: " << gnole << std::endl;
    if( this->worldComm().isMasterRank() )
    {
        MPI_File_write_at(fh, posInFile, &gnole, 1, MPI_INT32_T, &status );
        //MPI_File_write_ordered(fh, &gnole, size, MPI_INT32_T, &status );
    }
    posInFile+=sizeOfInt32_t;

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

    // All procs write :
    // - calculate every proc size of part to write
    int elidsWritingSize = elids.size()*sizeOfInt32_t;
    localOffset = 0;
    // - calculate every proc offset "localOffset" with Mpi_Exscan
    MPI_Exscan(&elidsWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
    // - write in file on cursor : posInFile + localOffset

    MPI_File_write_at(fh, posInFile+localOffset, elids.data(), elids.size(), MPI_INT32_T, &status );
    //MPI_File_write_ordered(fh, elids.data(), elids.size(), MPI_INT32_T, &status );
    sumOffsets = localOffset + elidsWritingSize;
    MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
    posInFile += sumOffsets; // all procs go to the end of last writing.

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

    // All procs write :
    // - calculate every proc size of part to write
    int ptidsWritingSize = ptids.size()*sizeOfInt32_t;
    localOffset = 0;
    // - calculate every proc offset "localOffset" with Mpi_Exscan
    MPI_Exscan(&ptidsWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
    // - write in file on cursor : posInFile + localOffset

    MPI_File_write_at(fh, posInFile+localOffset, ptids.data(), ptids.size(), MPI_INT32_T, &status );
    //MPI_File_write_ordered(fh, ptids.data(), ptids.size(), MPI_INT32_T, &status );

    sumOffsets = localOffset + ptidsWritingSize;
    MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
    posInFile += sumOffsets; // all procs go to the end of last writing.


    /* Write ghost elements */
    if ( 0  ) //this->worldComm().globalSize() > 1 )
    {
        // get ghost elements
        //auto r1 = markedelements(mesh, part->first, EntityProcessType::GHOST_ONLY );
        auto r1 = markedelements(mesh, markerid, EntityProcessType::GHOST_ONLY );
        auto gelt_it = r1.template get<1>();
        auto gelt_en = r1.template get<2>();

        std::string ghost_t = "g_" + this->elementType();
        // ghosts elements
        if( this->worldComm().isMasterRank() )
        {
            memset(buffer, '\0', sizeof(buffer));
            strcpy( buffer, ghost_t.c_str() );
            //buffer = ghost_t.c_str();
            MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
            //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status );
        }
        posInFile +=80;

        Feel::detail::MeshPoints<float> mpg( mesh.get(), this->worldComm(), gelt_it, gelt_en, true, true, true );
        //VLOG(1) << "material : " << p_it->first << " ghost nb element: " << __ne;

        int32_t gnoge = mpg.globalNumberOfElements();
        if( this->worldComm().isMasterRank() )
        {
            MPI_File_write_at(fh, posInFile, &gnoge, 1, MPI_INT32_T, &status );
            //MPI_File_write_ordered(fh, &gnoge, size, MPI_INT32_T, &status );
        }
        posInFile +=sizeOfInt32_t;

        /* Write elements ids */
        std::vector<int32_t> elids;
        for(auto it = gelt_it ; it != gelt_en; ++it )
        {
            auto const& elt = boost::unwrap_ref( *it );
            elids.push_back(elt.id() + mp.offsets_elts + 1);
        }


        // All procs write :
        // - calculate every proc size of part to write
        int elidsWritingSize = elids.size()*sizeOfInt32_t;
        localOffset = 0;
        // - calculate every proc offset "localOffset" with Mpi_Exscan
        MPI_Exscan(&elidsWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
        // - write in file on cursor : posInFile + localOffset

        MPI_File_write_at(fh, posInFile+localOffset, elids.data(), elids.size(), MPI_INT32_T, &status );
        // MPI_File_write_ordered(fh, elids.data(), elids.size(), MPI_INT32_T, &status );

        sumOffsets = localOffset + elidsWritingSize;
        MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
        posInFile += sumOffsets; // all procs go to the end of last writing.

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

        // All procs write :
        // - calculate every proc size of part to write
        int ptidsWritingSize = ptids.size()*sizeOfInt32_t;
        localOffset = 0;
        // - calculate every proc offset "localOffset" with Mpi_Exscan
        MPI_Exscan(&ptidsWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
        // - write in file on cursor : posInFile + localOffset

        MPI_File_write_at(fh, posInFile+localOffset, ptids.data(), ptids.size(), MPI_INT32_T, &status );
        // MPI_File_write_ordered(fh, ptids.data(), ptids.size(), MPI_INT32_T, &status );

        sumOffsets = localOffset + ptidsWritingSize;
        MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
        posInFile += sumOffsets; // all procs go to the end of last writing.

        // TODO : ASSERT ALL CURSORS EQUAL
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

                if( (dist = std::distance(__step->beginNodal(), __step->endNodal())) != 0)
                { LOG(INFO) << "NodalScalar: " << dist << std::endl; }
#if 0
                if( (dist = std::distance(__step->beginNodalVector(), __step->endNodalVector())) != 0)
                { LOG(INFO) << "NodalVector: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginNodalTensor2(), __step->endNodalTensor2())) != 0)
                { LOG(INFO) << "NodalTensor2: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginNodalTensor2Symm(), __step->endNodalTensor2Symm())) != 0)
                { LOG(INFO) << "NodalTensor2Symm: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginElementScalar(), __step->endElementScalar())) != 0)
                { LOG(INFO) << "ElementScalar: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginElementVector(), __step->endElementVector())) != 0)
                { LOG(INFO) << "ElementVector: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginElementTensor2(), __step->endElementTensor2())) != 0)
                { LOG(INFO) << "ElementTensor2: " << dist << std::endl; }
                if( (dist = std::distance(__step->beginElementTensor2Symm(), __step->endElementTensor2Symm())) != 0)
                { LOG(INFO) << "ElementTensor2Symm: " << dist << std::endl; }
#endif
                saveNodal<true>( __ts, __step, (__it == __ts->beginStep()), __step->beginNodal(), __step->endNodal() );
                saveNodal<false>( __ts, __step, (__it == __ts->beginStep()), __step->beginElement(), __step->endElement() );
#if 0
                saveNodal( __ts, __step, (__it == __ts->beginStep()), __step->beginNodalVector(), __step->endNodalVector() );
                saveNodal( __ts, __step, (__it == __ts->beginStep()), __step->beginNodalTensor2(), __step->endNodalTensor2() );
                saveNodal( __ts, __step, (__it == __ts->beginStep()), __step->beginNodalTensor2Symm(), __step->endNodalTensor2Symm() );

                saveElement( __ts, __step, (__it == __ts->beginStep()), __step->beginElementScalar(), __step->endElementScalar() );
                saveElement( __ts, __step, (__it == __ts->beginStep()), __step->beginElementVector(), __step->endElementVector() );
                saveElement( __ts, __step, (__it == __ts->beginStep()), __step->beginElementTensor2(), __step->endElementTensor2() );
                saveElement( __ts, __step, (__it == __ts->beginStep()), __step->beginElementTensor2Symm(), __step->endElementTensor2Symm() );
#endif

            }

            ++__it;
        }

        ++__ts_it;
    }
}


template<typename MeshType, int N>
template<bool IsNodal,typename Iterator>
void
ExporterEnsightGold<MeshType,N>::saveNodal( timeset_ptrtype __ts, typename timeset_type::step_ptrtype __step, bool isFirstStep, Iterator __var, Iterator en ) const
{
    tic();
    int size = 0;
    char buffer[ 80 ];

    MPI_Offset offset = 0;
    MPI_File fh;
    MPI_Status status;

    int sizeOfInt32_t;
    int sizeOfFloat;
    MPI_Type_size( MPI_INT32_T , &sizeOfInt32_t );
    MPI_Type_size( MPI_FLOAT , &sizeOfFloat );
    int localOffset, sumOffsets;
    while ( __var != en )
    {
        tic();
        auto const& nodalData =  __var->second;
        auto const& nodalField00 = unwrap_ptr( nodalData.second[0][0] );

        if ( !nodalField00.worldComm().isActive() ) return;

        uint16_type nComponents = invalid_uint16_type_value, nComponents1 = invalid_uint16_type_value, nComponents2 = invalid_uint16_type_value;
        std::string fileExt;
        bool isTensor2Symm = false;
        if ( nodalData.first == FunctionSpaceType::SCALAR )
        {
            nComponents = 1;
            nComponents1 = 1;
            nComponents2 = 1;
            fileExt = ".scl";
        }
        else if ( nodalData.first == FunctionSpaceType::VECTORIAL )
        {
            nComponents = 3;
            nComponents1 = 3;
            nComponents2 = 1;
            fileExt = ".vec";
        }
        else if ( nodalData.first == FunctionSpaceType::TENSOR2 )
        {
            nComponents = 9;
            nComponents1 = 3;
            nComponents2 = 3;
            fileExt = ".tsr";
        }
        else if ( nodalData.first == FunctionSpaceType::TENSOR2_SYMM )
        {
            nComponents = 6;
            nComponents1 = 3;
            nComponents2 = 3;
            fileExt = ".tsrs";
            isTensor2Symm = true;
        }


        std::ostringstream __varfname;

        auto __mesh = __step->mesh();

        __varfname << this->path() << "/" << __ts->name() << "." << __var->first;
        if( ! boption( _name="exporter.ensightgold.merge.markers") )
        { __varfname << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank(); }
        /* if we want to pack data in several files instead of one */
        /* we compute an index to add to the filename */
        if( boption( _name = "exporter.ensightgold.merge.timesteps")
        && ioption( _name="exporter.ensightgold.pack.timesteps" ) > 1 )
        {
            // timestep indices start at 1
            __varfname << "." << ((__step->index() - TS_INITIAL_INDEX) / ioption( _name="exporter.ensightgold.pack.timesteps" ) + 1);
        }
        // add extension
        __varfname << fileExt;
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

        //init MPI_Info object fromhints defined as environment variables
        //MPI_Info info = initIoInfoFromEnvVars();
        if( boption( _name = "exporter.ensightgold.merge.timesteps") )
            MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL , &fh );
        else
            MPI_File_open( this->worldComm().comm(), str, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL , &fh );

        //MPI_Info_free(&info);
        free(str);

        // INIT CURSOR IN FILE
        posInFile = 0;
        Feel::detail::FileIndex index;

        if( boption( _name = "exporter.ensightgold.merge.timesteps") )
        {
            // first read
            index.read(fh);

            /* Move to the beginning of the fie index section */
            /* to overwrite it */
            if ( index.defined() && (__step->index() - TS_INITIAL_INDEX) > 0 ) {
                // MPI_File_seek_shared(fh, index.fileblock_n_steps, MPI_SEEK_SET);
                posInFile = index.fileblock_n_steps;
            }
            else {
                // we position the cursor at the beginning of the file
                // MPI_File_seek_shared(fh, 0, MPI_SEEK_SET);
                posInFile = 0;
            }

            /* Write time step start */
            if( this->worldComm().isMasterRank() )
            {
                memset(buffer, '\0', sizeof(buffer));
                strncpy(buffer, "BEGIN TIME STEP", 80);
                // buffer = "BEGIN TIME STEP";
                MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                // memset(buffer, '\0', sizeof(buffer));
                // strcpy(buffer,"BEGIN TIME STEP");
                // MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
            }
            posInFile+=80;

            /* add the beginning of the new block to the file */
            // MPI_File_get_position_shared(fh, &offset);
            index.add( posInFile );
        }

        if( this->worldComm().isMasterRank() )
        {
            memset(buffer, '\0', sizeof(buffer));
            strncpy(buffer, nodalField00.name().c_str(), 80);
            // buffer = __var->second.name().c_str();
            MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
            // MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
        }
        posInFile+=80;

        /* handle faces data */
        if ( boption( _name="exporter.ensightgold.save-face" ) )
        {
            BOOST_FOREACH( auto m, __mesh->markerNames() )
            {
                if ( m.second[1] != __mesh->nDim-1 )
                    continue;
                VLOG(1) << "writing face with marker " << m.first << " with id " << m.second[0];
                //auto pairit = __mesh->facesWithMarker( m.second[0], this->worldComm().localRank() );
                auto rangeMarkedFaces = __mesh->facesWithMarker( m.second[0], this->worldComm().localRank() );
                auto fit = std::get<0>( rangeMarkedFaces );
                auto fen = std::get<1>( rangeMarkedFaces );

                // auto fit = pairit.first;
                // auto fen = pairit.second;

                Feel::detail::MeshPoints<float> mp( __mesh.get(), this->worldComm(), fit, fen, true, true, true );
                int __ne = std::distance( fit, fen );
                if ( __ne == 0 ) continue;
                int nverts = boost::unwrap_ref( *fit ).numLocalVertices;
                DVLOG(2) << "Faces : " << __ne << "\n";

                if( this->worldComm().isMasterRank() )
                {
                    memset(buffer, '\0', sizeof(buffer));
                    strcpy( buffer, "part" );
                    //buffer = "part";
                    MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                    // memset(buffer, '\0', sizeof(buffer));
                    // strcpy( buffer, "part" );
                    // MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                }
                posInFile+=80;

                if( this->worldComm().isMasterRank() )
                {
                    int32_t partid = m.second[0] + 1;
                    MPI_File_write_at(fh, posInFile, &partid, 1, MPI_INT32_T, &status);
                    //MPI_File_write_ordered(fh, &partid, size, MPI_INT32_T, &status);
                }
                posInFile+=sizeOfInt32_t;

                if( this->worldComm().isMasterRank() )
                {
                    memset(buffer, '\0', sizeof(buffer));
                    strcpy( buffer, "coordinates" );
                    //buffer = "coordinates";
                    MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                    // memset(buffer, '\0', sizeof(buffer));
                    // strcpy( buffer, "coordinates" );
                    // MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                }
                posInFile+=80;

                // write values
                // fit = pairit.first;
                // fen = pairit.second;
                fit = std::get<2>( rangeMarkedFaces )->begin();

                int nfaces = mp.ids.size();
                std::vector<float> field( nComponents*nfaces, 0.0 );
                for( ; fit != fen; ++fit )
                {
                    auto const& face = boost::unwrap_ref( * fit );
                    for ( uint16_type j = 0; j < nverts; j++ )
                    {
                        size_type pid = mp.old2new[face.point( j ).id()]-1;
                        for ( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                        {
                            if ( c1 >= nodalData.second.size() )
                                continue;
                            for ( uint16_type c2 = 0; c2 < nComponents2; ++c2 )
                            {
                                if ( c2 >= nodalData.second[c1].size() )
                                    continue;
                                auto const& nodalField = unwrap_ptr( nodalData.second[c1][c2] );
                                uint16_type c = c2*nComponents1+c1;
                                size_type global_node_id = nfaces*c + pid ;
                                size_type thedof =  nodalField.start() +
                                    nodalField.functionSpace()->dof()->faceLocalToGlobal( face.id(), j, c ).index();
                                field[global_node_id] = nodalField.globalValue( thedof );
                            }
                        }
                    }
                }
                /* Write each component separately */

                // All procs write :
                // - calculate every proc size of part to write
                int fieldWritingSize = nfaces*sizeOfFloat;
                localOffset = 0;
                // - calculate every proc offset "localOffset" with Mpi_Exscan
                MPI_Exscan(&fieldWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
                sumOffsets = localOffset + fieldWritingSize;
                MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
                // - write in file on cursor : posInFile + localOffset
                for ( uint16_type c = 0; c < nodalField00.nComponents; ++c )
                {
                    MPI_File_write_at(fh, posInFile + localOffset + c*sumOffsets, \
                                      field.data() + nfaces * c, nfaces, MPI_FLOAT, &status);
                    // MPI_File_write_ordered(fh, field.data() + nfaces * c, nfaces, MPI_FLOAT, &status);
                }
                posInFile += nodalField00.nComponents * sumOffsets;
            } // boundaries loop
        }
        toc("saveNodal intro",FLAGS_v>0);
        /* handle elements */
        for( std::set<int>::iterator mit = M_markersToWrite.begin(); mit != M_markersToWrite.end(); mit++)
        {
            tic();
            if( this->worldComm().isMasterRank() )
            {
                memset(buffer, '\0', sizeof(buffer));
                strcpy( buffer, "part" );
                //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                int32_t partid = *mit + 1;
                MPI_File_write_at(fh, posInFile+80, &partid, 1, MPI_INT32_T, &status);
                DVLOG(2) << "part " << buffer << "\n";
                if ( IsNodal )
                {
                    strcpy(buffer, "coordinates");
                }
                else
                {
                    memset(buffer, '\0', sizeof(buffer));
                    strncpy(buffer, this->elementType().c_str(), 80);
                }
                MPI_File_write_at(fh, posInFile+80+sizeOfInt32_t, buffer, sizeof(buffer), MPI_CHAR, &status);
            }
            posInFile+=160+sizeOfInt32_t;

            /* we get that from the local processor */
            /* We do not need the renumbered global index */
            //auto r = markedelements(__mesh,(boost::any)p_it->first,EntityProcessType::ALL);
            // auto r = markedelements(__mesh, *mit, EntityProcessType::ALL);
            auto r = markedelements(__mesh, *mit );
            auto elt_it = r.template get<1>();
            auto elt_en = r.template get<2>();
            toc("saveNodal part",FLAGS_v>0);


            /* create an array to store data per node */
            int nValuesPerComponent = 0;

            //int npts = mp.ids.size();
            //size_type __field_size = npts*nComponents;

            Eigen::VectorXf __field;// = Eigen::VectorXf::Zero( __field_size );
            //VLOG(1) << "field size=" << __field_size;

            int reorder_tensor2symm[6] = { 0,3,4,1,5,2 };

            /* Loop on the elements */
            tic();
            auto const& d = nodalField00.functionSpace()->dof().get();

            if constexpr ( IsNodal )
                {
                    tic();
                    M_cache_mp.try_emplace( *mit, __step->mesh().get(), this->worldComm(), elt_it, elt_en, true, true, true );
                    auto& mp = M_cache_mp.at( *mit );
                    toc( "ExporterEnsightGold::writeVariables MeshPoints", FLAGS_v > 0 );

                    nValuesPerComponent = mp.ids.size();
                    size_type __field_size = nValuesPerComponent*nComponents;
                    VLOG(1) << "field size=" << __field_size;
                    __field = Eigen::VectorXf::Zero( __field_size );

                    int index = 0;
                    const int np = __step->mesh()->numLocalVertices();
                    for ( ; elt_it != elt_en; ++elt_it )
                    {
                        auto const& elt = elt_it->get();
                        auto const& locglob_ind = d->localToGlobalIndices( elt.id() );

                        /* looop on the ccomponents is outside of the loop on the vertices */
                        /* because we need to pack the data in the x1 x2 ... xn y1 y2 ... yn z1 z2 ... zn order */
                        for ( uint16_type p = 0; p < np; ++p )
                        {
                            size_type ptid = mp.old2new[elt.point( p ).id()]-1;
                            for ( uint16_type c1 = 0; c1 < nodalData.second.size(); ++c1 )
                            {
                                for ( uint16_type c2 = 0; c2 < nodalData.second[c1].size(); ++c2 )
                                {
                                    auto const& nodalField = unwrap_ptr( nodalData.second[c1][c2] );
                                    uint16_type cMap = c2*nComponents1+c1;
                                    if ( isTensor2Symm )
                                        cMap = reorder_tensor2symm[Feel::detail::symmetricIndex( c1,c2, nComponents1 )];

                                    size_type global_node_id = mp.ids.size()*cMap + ptid ;
                                    //LOG(INFO) << elt_it->get().point( p ).id() << " " << ptid << " " << global_node_id << std::endl;
                                    DCHECK( ptid < __step->mesh()->numPoints() ) << "Invalid point id " << ptid << " element: " << elt.id()
                                                                                 << " local pt:" << p
                                                                                 << " mesh numPoints: " << __step->mesh()->numPoints();
                                    DCHECK( global_node_id < __field_size ) << "Invalid dof id : " << global_node_id << " max size : " << __field_size;

                                    size_type dof_id = locglob_ind(d->localDofId(p,0));

                                    __field(global_node_id) = nodalField.globalValue( dof_id );
                                    DVLOG(3) << "v[" << global_node_id << "]=" << nodalField.globalValue( dof_id ) << "  dof_id:" << dof_id;
                                }
                            }
                        }

                        /* increment index of vertex */
                        index++;
                    }
                } // IsNodal
            else
            {
                nValuesPerComponent = std::distance( elt_it, elt_en );
                size_type __field_size = nValuesPerComponent*nComponents;
                VLOG(1) << "field size=" << __field_size;
                __field = Eigen::VectorXf::Zero( __field_size );
                for ( index_type e=0 ; elt_it != elt_en; ++elt_it, ++e )
                {
                    auto const& elt = elt_it->get();
                    auto const& locglob_ind = d->localToGlobalIndices( elt.id() );
                    size_type dof_id = locglob_ind(d->localDofId(0,0));
                    for ( uint16_type c1 = 0; c1 < nodalData.second.size(); ++c1 )
                    {
                        for ( uint16_type c2 = 0; c2 < nodalData.second[c1].size(); ++c2 )
                        {
                            auto const& nodalField = unwrap_ptr( nodalData.second[c1][c2] );
                            uint16_type cMap = c2*nComponents1+c1;
                            if ( isTensor2Symm )
                                cMap = reorder_tensor2symm[Feel::detail::symmetricIndex( c1,c2, nComponents1 )];

                            size_type global_node_id = cMap*nValuesPerComponent + e;
                            __field(global_node_id) = nodalField.globalValue( dof_id );
                            DVLOG(3) << "v[" << global_node_id << "]=" << nodalField.globalValue( dof_id ) << "  dof_id:" << dof_id;
                        }
                    }
                }
            }

            toc("saveNodal element loop",FLAGS_v>0);
            tic();
            /* Write each component separately */
            // All procs write :
            // - calculate every proc size of part to write
            int fieldWritingSize = nValuesPerComponent*sizeOfFloat;
            localOffset = 0;
            // - calculate every proc offset "localOffset" with Mpi_Exscan
            MPI_Exscan(&fieldWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
            sumOffsets = localOffset + fieldWritingSize;
            MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
            // - write in file on cursor : posInFile + localOffset +
            for ( uint16_type c = 0; c < nComponents; ++c )
            {
                MPI_File_write_at(fh, posInFile + localOffset + c*sumOffsets, \
                                  ((float *)(__field.data())) + nValuesPerComponent * c, nValuesPerComponent, MPI_FLOAT, &status);
                // MPI_File_write_ordered(fh, ((float *)(__field.data().begin())) + nValuesPerComponent * c, nValuesPerComponent, MPI_FLOAT, &status);
            }
            posInFile += nComponents*sumOffsets;
            toc("saveNodal write part",FLAGS_v>0);
        } // parts loop

        if( boption(_name="exporter.ensightgold.merge.timesteps") )
        {
            /* write timestep end */
            if( this->worldComm().isMasterRank() )
            {
                memset(buffer, '\0', sizeof(buffer));
                strncpy(buffer,"END TIME STEP", 80);
                MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                // MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
            }
            posInFile+=80;
            // write back the file index
            index.write( fh );
        }
        DVLOG(2) << "[ExporterEnsightGold::saveNodal] saving " << __varfname.str() << "done\n";

        MPI_File_close(&fh);

        ++__var;
    }
    toc("saveNodal", FLAGS_v>0);
}

template<typename MeshType, int N>
template<typename Iterator>
void
ExporterEnsightGold<MeshType,N>::saveElement( timeset_ptrtype __ts, typename timeset_type::step_ptrtype __step, bool isFirstStep, Iterator __evar, Iterator __evaren ) const
{
#if 0 // TODO
    int size;
    char buffer[ 80 ];

    MPI_File fh;
    MPI_Offset offset = 0;
    MPI_Status status;

    int sizeOfInt32_t;
    int  sizeOfFloat;
    MPI_Type_size( MPI_INT32_T , &sizeOfInt32_t );
    MPI_Type_size( MPI_FLOAT , &sizeOfFloat );
    int localOffset, sumOffsets;
    while ( __evar != __evaren )
    {
        if ( !__evar->second.worldComm().isActive() ) return;

        std::ostringstream __evarfname;

        __evarfname << this->path() << "/" << __ts->name() << "." << __evar->first;

        if( ! boption( _name="exporter.ensightgold.merge.markers") )
        { __evarfname << "-" << this->worldComm().globalSize() << "_" << this->worldComm().localRank(); }
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
        else if(__evar->second.is_tensor2symm)
        { __evarfname << ".tsrs"; }
        else
        { __evarfname << ".scl"; LOG(ERROR) << "Could not detect data type (scalar, vector, tensor2, tensor2symm). Defaulted to scalar." << std::endl; }

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

        //init MPI_Info object fromhints defined as environment variables
        //MPI_Info info = initIoInfoFromEnvVars();
        if( boption( _name = "exporter.ensightgold.merge.timesteps") )
            MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
        else
            MPI_File_open( this->worldComm().comm(), str, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );

        //MPI_Info_free(&info);
        free(str);

        // INIT CURSOR IN FILE
        posInFile = 0;

        Feel::detail::FileIndex index;

        if( boption( _name = "exporter.ensightgold.merge.timesteps") )
        {
            // first read
            index.read(fh);

            /* Move to the beginning of the fie index section */
            /* to overwrite it */
            if ( index.defined() && (__step->index() - TS_INITIAL_INDEX) > 0 ) {
                // MPI_File_seek_shared(fh, index.fileblock_n_steps, MPI_SEEK_SET);
                posInFile = index.fileblock_n_steps;
            }
            else {
                // we position the cursor at the beginning of the file
                //MPI_File_seek_shared(fh, 0, MPI_SEEK_SET);
                posInFile = 0;
            }

            /* Write time step start */
            if( this->worldComm().isMasterRank() )
            {
                memset(buffer, '\0', sizeof(buffer));
                strncpy(buffer, "BEGIN TIME STEP",80);
                MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                // MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
            }
            posInFile+=80;

            /* add the beginning of the new block to the file */
            // MPI_File_get_position_shared(fh, &offset);
            index.add( posInFile );
        }

        if( this->worldComm().isMasterRank() )
        {
            memset(buffer, '\0', sizeof(buffer));
            strncpy( buffer, __evar->second.name().c_str(), 80);
            // buffer = __evar->second.name().c_str();
            MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
            // MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
        }
        posInFile+=80;
        // iterate over the markers
        for( std::set<int>::iterator mit = M_markersToWrite.begin() ; mit != M_markersToWrite.end(); mit++ )
        {
            if( this->worldComm().isMasterRank() )
            {
                memset(buffer, '\0', sizeof(buffer));
                strncpy(buffer, "part", 80);
                MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                int32_t partid = *mit + 1;
                MPI_File_write_at(fh, posInFile+80, &partid, 1, MPI_INT32_T, &status);
                DVLOG(2) << "part " << buffer << "\n";
                memset(buffer, '\0', sizeof(buffer));
                strncpy(buffer, this->elementType().c_str(), 80);
                MPI_File_write_at(fh, posInFile+80+sizeOfInt32_t, buffer, sizeof(buffer), MPI_CHAR, &status);
                DVLOG(2) << "element type " << buffer << "\n";

            }
            posInFile+=160+sizeOfInt32_t;

            uint16_type nComponents = __evar->second.nComponents;

            if ( __evar->second.is_vectorial )
                nComponents = 3;
            if ( __evar->second.is_tensor2 )
                nComponents = 9;
            if ( __evar->second.is_tensor2symm )
                nComponents = 6;

            int nc = __evar->second.nComponents;
            if ( __evar->second.is_tensor2symm )
                nc = __evar->second.nComponents1*(__evar->second.nComponents1+1)/2;

            size_type __field_size = nComponents * __evar->second.size()/nc;

            Eigen::VectorXf __field = Eigen::VectorXf::Zero( __field_size );

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
            int reorder_tensor2symm[6] = { 0, 3, 1, 4, 5, 2 };
            auto const& d = __evar->second.functionSpace()->dof().get();
            tic();
            size_type e = 0;
            for ( auto elt_it = elt_st ; elt_it != elt_en; ++elt_it, ++e )
            {
                auto const& elt = boost::unwrap_ref( *elt_it );
                DVLOG(2) << "pid : " << this->worldComm().globalRank()
                         << " elt_it :  " << elt.id()
                         << " e : " << e << "\n";

                auto const& locglob_ind = d->localToGlobalIndices( elt.id() );
                
                for ( int c = 0; c < nComponents; ++c )
                {
                    uint16_type c1= c;
                    if ( __evar->second.is_tensor2symm )
                        c1=reorder_tensor2symm[c];
                    int nc = __evar->second.nComponents;
                    if ( __evar->second.is_tensor2symm )
                        nc = __evar->second.nComponents1*(__evar->second.nComponents1+1)/2;
                        
                    size_type global_node_id = c1*ncells+e ;

                    if ( c < nc)
                    {
                        size_type dof_id = locglob_ind( d->localDofId( 0, c ) );
                        
                        DVLOG(2) << "c : " << c
                                 << " gdofid: " << global_node_id
                                 << " dofid : " << dof_id
                                 << " f.size : " <<  __field.size()
                                 << " e.size : " <<  __evar->second.size()
                                 << "\n";

                        __field(global_node_id) = __evar->second.globalValue( dof_id );

#if 1
                        //__field[global_node_id] = __evar->second.globalValue(dof_id);
                        DVLOG(2) << "c : " << c
                                 << " gdofid: " << global_node_id
                                 << " dofid : " << dof_id
                                 << " field :  " << __field[global_node_id]
                                 << " evar: " << __evar->second.globalValue( dof_id ) << "\n";
#endif
                    }
                }
            }
            toc("saveElement element loop",FLAGS_v>0);
            /* Write each component separately */
            tic();
            // All procs write :
            // - calculate every proc size of part to write
            int fieldWritingSize = ncells*sizeOfFloat;
            localOffset = 0;
            // - calculate every proc offset "localOffset" with Mpi_Exscan
            MPI_Exscan(&fieldWritingSize, &localOffset, 1, MPI_INT, MPI_SUM, this->worldComm());
            sumOffsets = localOffset + fieldWritingSize;
            MPI_Bcast(&sumOffsets, 1, MPI_INT, this->worldComm().globalSize()-1, this->worldComm());
            // - write in file on cursor : posInFile + localOffset

            for ( uint16_type c = 0; c < nComponents; ++c )
            {
                MPI_File_write_at(fh, posInFile + c*sumOffsets + localOffset, \
                                  __field.data() + ncells * c, ncells, MPI_FLOAT, &status);
                // MPI_File_write_ordered(fh, __field.data().begin() + ncells * c, ncells, MPI_FLOAT, &status);
            }
            toc("saveElement write",FLAGS_v>0);
            posInFile += nComponents*sumOffsets;
        }

        if( boption(_name="exporter.ensightgold.merge.timesteps") )
        {
            /* write timestep end */
            if( this->worldComm().isMasterRank() )
            {
                memset(buffer, '\0', sizeof(buffer));
                strcpy(buffer,"END TIME STEP");
                // MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
            }
            posInFile+=80;
            // write back the file index
            index.write( fh );
        }

        MPI_File_close(&fh);

        DVLOG(2) << "[ExporterEnsightGold::saveElement] saving " << __evarfname.str() << "done\n";
        ++__evar;
    }
#endif
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
