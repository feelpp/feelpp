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
ExporterEnsightGold<MeshType,N>::ExporterEnsightGold( std::string const& exp_prefix, worldcomm_ptr_t const& worldComm )
    :
    super( exp_prefix, worldComm )
{
    init();
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



    M_nodesOrderingInElementToEnsight[ Simplex<1,2>::type() + "_1d_g2" ] = { 0,2,1 };
    M_nodesOrderingInElementToEnsight[ Simplex<2,2>::type() + "_2d_g2" ] = { 0,1,2,4,5,3 };
    M_nodesOrderingInElementToEnsight[ Simplex<3,2>::type() + "_3d_g2" ] = { 0,1,2,3,5,6,4,7,8,9 };

    M_nodesOrderingInElementToEnsight[ Hypercube<1,2>::type() + "_1d_g2" ] = { 0,2,1 };
    M_nodesOrderingInElementToEnsight[ Hypercube<2,2>::type() + "_2d_g2" ] = { 0,1,2,3,4,5,6,7 };
    M_nodesOrderingInElementToEnsight[ Hypercube<3,2>::type() + "_3d_g2" ] = { 0,1,2,3,4,5,6,7,8,9,10,11,17,12,16,18,13,19,14,15 };


    /* TODO Do a cleanup of previous stored files */
    /* to avoid conflicts */
    /* example case: where a new simulation has fewer timesteps */
    /* it can be confusing if files with higher timesteps remain */

    /* Init number of digit for maximum time step */
    M_timeExponent = 4;

    M_mergeTimeSteps = boption( _name="exporter.ensightgold.merge.timesteps" );
    M_packTimeSteps = ioption( _name="exporter.ensightgold.pack.timesteps" );
    if ( M_mergeTimeSteps && ( M_packTimeSteps == 1 ) )
        M_mergeTimeSteps = false;
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::save( steps_write_on_disk_type const& stepsToWriteOnDisk ) const
{
    tic();

    for ( auto const& [__ts,steps] : stepsToWriteOnDisk  )
    {
        if ( steps.empty() && __ts->numberOfSteps() == 0 && __ts->hasMesh() ) // save only the mesh
        {
            tic();
            int stepIndex = TS_INITIAL_INDEX;
            writeGeoFiles( __ts, __ts->mesh(), stepIndex, true );
            toc("ExporterEnsightGold::save geo",FLAGS_v>1);
        }
        else
        {
            for ( auto const&  __step : steps )
            {
                mesh_ptrtype mesh = __step->mesh();
                int stepIndex = __step->activeIndex();
                bool isFirstStep = ( stepIndex == (*__ts->beginStep())->index() );
                tic();
                writeGeoFiles( __ts, mesh,stepIndex,isFirstStep );
                toc("ExporterEnsightGold::save geo",FLAGS_v>1);
                tic();
                writeVariableFiles( __ts, __step );
                toc("ExporterEnsightGold::save variables",FLAGS_v>1);
            }
        }
    }

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
        __out << "MULTIPLE_CASEFILES\n"
              << "total number of cfiles: 1 " << std::endl
              << "cfiles global path: " << fs::current_path().string() << "\n"
              << "cfiles pattern: "<< this->prefix() << ".case\n"
              << "cfiles start number: 0\n"
              << "cfiles increment: 1\n\n";
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

        __outparaview << "number of servers: 1 " << std::endl;
        __outparaview << "#Server " << 1 << "\n"
                      << "machine id: " << mpi::environment::processor_name() << "\n"
                      << "executable: /usr/local/bin/ensight76/bin/ensight7.server\n"
                      << "data_path: " << fs::current_path().string() << "\n"
                      << "casefile: " << this->prefix() << ".case\n";
    }
}
template<typename MeshType, int N>
template<typename Iterator,typename TSt>
void
ExporterEnsightGold<MeshType,N>::writeCaseFileVariables( Iterator it, Iterator end,
                                                         std::string const& loc,
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

        if ( M_mergeTimeSteps )
        {
            __out << type << " per " << loc << ": "
                  << __ts->index() << " 1 " // << *__ts_it->beginStep() << " "
                  << unwrap_ptr(it->second.second[0][0]).name() << " " << __ts->name() << "." << it->first;
            /* if we want to pack data in several files instead of one, we add an index to the filename */
            if ( M_packTimeSteps )
                __out << ".*";
            __out << "." << ext << std::endl;
        }
        else
        {
            __out << type << " per " << loc << ": "
                  << __ts->index() << " " // << *__ts_it->beginStep() << " "
                  << unwrap_ptr(it->second.second[0][0]).name() << " " << __ts->name() << "." << it->first;
            __out << "." << ext << "." << std::string(M_timeExponent, '*') << std::endl;
        }
    }
}

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeCaseFile() const
{
    // only on proc 0
    if( this->worldComm().isMasterRank() )
    {
        std::ostringstream filestr;

        timeset_const_iterator __ts_it = this->beginTimeSet();
        timeset_const_iterator __ts_en = this->endTimeSet();

        while ( __ts_it != __ts_en )
        {
            timeset_ptrtype __ts = *__ts_it;

            filestr << this->path() << "/"
                << __ts->name();
            filestr << ".case";

            std::ofstream __out( filestr.str().c_str() );

            if ( __out.fail() )
            {
                CHECK( false ) << "cannot open " << filestr.str()  << "\n";
            }

            __out << "FORMAT:\n"
                  << "type: ensight gold\n"
                  << "GEOMETRY:\n";

            switch ( this->exporterGeometry() )
            {
                case EXPORTER_GEOMETRY_STATIC:
                {
                    __out << "model: " << __ts->name();
                    __out << ".geo";
                }
                break;
                default:
                case EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY:
                case EXPORTER_GEOMETRY_CHANGE:
                {
                    if( M_mergeTimeSteps )
                    {
                        __out << "model: " << __ts->index() << " 1 " << __ts->name();
                        /* if we want to pack data in several files instead of one, we add an index to the filename */
                        if ( M_packTimeSteps > 0 )
                            __out << ".*";
                        __out << ".geo";
                    }
                    else
                    {
                        __out << "model: " << __ts->index() << " " << __ts->name();
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
                                        "node",
                                        __ts, __out );
                writeCaseFileVariables( ( *__tstp_it )->beginElement(),
                                        ( *__tstp_it )->endElement(),
                                        "element",
                                        __ts, __out );

            }

            __out << "TIME:\n";

            if ( auto __firstActiveStep =  __ts->firstActiveStep() )
            {
                __out << "time set:        " << __ts->index() << "\n"
                    << "number of steps: " << __ts->numberOfActiveSteps() << "\n"
                    << "filename start number: " << __firstActiveStep->activeIndex() << "\n"
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
            typename timeset_type::step_const_iterator __its = __ts->beginStep();
            typename timeset_type::step_const_iterator __ens = __ts->endStep();
            for ( ; __its != __ens ; ++__its )
            {
                if ( (*__its)->isIgnored() )
                    continue;
                __out << std::scientific << std::setprecision( 6 ) << ( *__its )->time() << " ";

                if ( __l++ % 10 == 0 )
                    __out << "\n";
            }
            __out << "\n";

            if ( M_mergeTimeSteps )
            {
                auto ts = *(this->beginTimeSet());

                // create several filesets if we pack data into groups in several files
                if( M_packTimeSteps > 0 )
                {
                    int stepsPerFile = M_packTimeSteps;
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
ExporterEnsightGold<MeshType,N>::writeGeoFiles( timeset_ptrtype __ts, mesh_ptrtype mesh, int timeIndex, bool isFirstStep ) const
{
    // prepare/udate cache
    bool buildNewCache = true;
    auto itFindCache = M_cache_mp.find( __ts->name() );
    if ( itFindCache != M_cache_mp.end() )
    {
        auto mpFound = itFindCache->second;
        bool sameMesh = mesh->isSameMesh( mpFound->mesh() );
        if ( this->exporterGeometry() == EXPORTER_GEOMETRY_STATIC  || this->exporterGeometry() == EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY )
        {
            CHECK( sameMesh ) << "the mesh has changed but you use EXPORTER_GEOMETRY_STATIC or EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY";
            buildNewCache = false;
        }
        if ( this->exporterGeometry() == EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY )
            mpFound->updateNodesCoordinates();
    }
    if ( buildNewCache )
    {
        tic();
        M_cache_mp[__ts->name()] = std::make_shared<mesh_contiguous_numbering_mapping_type>( mesh.get(), false, this->meshFragmentation() );
        toc( "ExporterEnsightGold::writeGeoFiles init cache", FLAGS_v > 0 );
        // clear others caches with export of fields
        M_mapNodalArrayToDofId.clear();
        M_mapElementArrayToDofId.clear();
    }

    // mesh already write, do nothing
    if ( !isFirstStep && this->exporterGeometry() == EXPORTER_GEOMETRY_STATIC )
        return;

    //int size;
    //char buffer[80];

    MPI_File fh;
    //MPI_Status status;
    //MPI_Info info;
    //MPI_Offset offset = 0;

    auto const& mp = *M_cache_mp.find( __ts->name() )->second;

    bool writeNewGeoFile = true;

    std::ostringstream __geofname;
    __geofname << this->path() << "/" << __ts->name();

    if( M_mergeTimeSteps )
    {
        if ( M_packTimeSteps > 0 )
        {
            // timestep indices start at 1
            int startTsIndex = (__ts->numberOfSteps() == 0 )? TS_INITIAL_INDEX : (*__ts->beginStep())->index();
            __geofname << "." << ((timeIndex - startTsIndex)/M_packTimeSteps + 1);
            writeNewGeoFile = ( ( (timeIndex - startTsIndex) % M_packTimeSteps ) == 0 ) ;
        }
        else
            writeNewGeoFile = isFirstStep;
    }
    __geofname << ".geo";
    if ( !M_mergeTimeSteps && this->exporterGeometry() != EXPORTER_GEOMETRY_STATIC )
        __geofname << "." << std::setfill( '0' ) << std::setw( M_timeExponent ) << timeIndex;


    M_filename =  __geofname.str();


    /* Open File with MPI IO */
    char * str = strdup(__geofname.str().c_str());

    if ( writeNewGeoFile )
    {
        /* Check if file exists and delete it, if so */
        /* (MPI IO does not have a truncate mode ) */
        if(this->worldComm().isMasterRank() && fs::exists(str))
        {
            MPI_File_delete(str, MPI_INFO_NULL);
        }
        MPI_Barrier( this->worldComm().comm() );
    }

    if( M_mergeTimeSteps && !writeNewGeoFile )
        MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL , &fh );
    else
        MPI_File_open( this->worldComm().comm(), str, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL , &fh );

    //MPI_Info_free(&info);
    free(str);
    // Initializing cursor in file
    posInFile=0;

    bool writeHeaderBeginFile = writeNewGeoFile;
    bool writeBeginEndTimeSet = M_mergeTimeSteps;

    Feel::detail::FileIndex index( this->worldCommPtr() );
    // read previous FILE_INDEX
    if ( M_mergeTimeSteps && !writeNewGeoFile )
    {
        index.read(fh);
        if ( index.defined() )
            posInFile=index.nextFreePosFile();
    }

    /* Write the geo file */
    this->writeGeoMarkers( fh, mp, writeHeaderBeginFile, writeBeginEndTimeSet, index );

    // write new FILE_INDEX
    if ( M_mergeTimeSteps )
        index.write( fh, posInFile );

    /* close file */
    MPI_File_close(&fh);

}


template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeGeoMarkers(MPI_File fh, mesh_contiguous_numbering_mapping_type const& mp, bool writeHeaderBeginFile, bool writeBeginEndTimeSet, Feel::detail::FileIndex & index ) const
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
    MPI_Status status;

    /* Write file header */
    // only write C Binary if we are not mergin timesteps
    // as it is oalready writtent at the beginning of the file
    if( writeHeaderBeginFile )
    {
        char buffer[80];
        if( this->worldComm().isMasterRank() )
        {
            memset(buffer, '\0', sizeof(buffer));
            strncpy(buffer, "C Binary", 80);
            MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
        }
        posInFile+=80;
        //LOG(INFO) << "wrote " << buffer << std::endl;
    }

    // Write time step start
    if ( writeBeginEndTimeSet )
    {
        if( this->worldComm().isMasterRank() )
        {
            char buffer[80];
            memset(buffer, '\0', sizeof(buffer));
            strncpy(buffer, "BEGIN TIME STEP", 80);
            MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
        }
        // all procs increment their offset
        posInFile+=80;

        /* add the beginning of the new block to the file */
        index.add( posInFile );
    }
    if( this->worldComm().isMasterRank() )
    {
        //char buffer[320];
        // get only the filename (maybe with full path)
        fs::path gp = M_filename;
        std::string theFileName = gp.filename().string();
        CHECK( theFileName.length() < 80 ) << "the file name is too long : theFileName=" << theFileName << "\n";
        char buffer2[320];
        memset(buffer2, '\0', sizeof(buffer2));
        strcpy( buffer2, theFileName.c_str() );
        strcpy( &buffer2[80], "elements" );
        //strcpy( &buffer2[160], "node id given" );
        //strcpy( &buffer2[240], "element id given" );
        strcpy( &buffer2[160], "node id assign" );
        strcpy( &buffer2[240], "element id assign" );
        //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
        MPI_File_write_at(fh, posInFile, buffer2, sizeof(buffer2), MPI_CHAR, &status);
    }
    posInFile+=320;



#if 0
    /* Write faces */
    if ( boption( _name="exporter.ensightgold.save-face" ) )
    {
        auto mesh = mp.mesh();
        for( auto & m : mesh->markerNames() )
        {
            this->writeGeoMarkedFaces(fh, mesh, m);
        }
    }
#endif
    // Write geo elements parts : nodes and connectivity
    for ( auto const& [part,nameAndRangeElt] : mp.partIdToRangeElement() )
    {
        this->writeGeoMarkedElements( fh, mp, part );
    }

    // Write time step end
    if ( writeBeginEndTimeSet )
    {
        if ( this->worldComm().isMasterRank() )
        {
            char buffer[80];
            memset(buffer, '\0', 80);
            strcpy(buffer,"END TIME STEP");
            MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
        }
        posInFile+=80;
    }

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

#if 0
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
        strncpy( buffer, this->elementType().c_str(), 80-1);
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

#else
template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeGeoMarkedElements( MPI_File fh, mesh_contiguous_numbering_mapping_type const& mp, int part ) const
{
    MPI_Status status;
    char buffer[80];
    int sizeOfInt32_t;
    int sizeOfFloat;
    MPI_Type_size( MPI_INT32_T , &sizeOfInt32_t );
    MPI_Type_size( MPI_FLOAT , &sizeOfFloat );

    auto mesh = mp.mesh();

    rank_type currentPid = mesh->worldComm().localRank();

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
        int32_t partid = part + 1;
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
        //sprintf(buffer2, "Marker %d (%s)", (int)(markerid), mesh->markerName(markerid).substr(0, 32).c_str());
        std::string partName = mp.name( part );
        if( partName.empty() )
            partName = (boost::format("marker_%1%")%part).str();
        partName.resize( std::min((int)partName.size(),(int)80) );
        strcpy( buffer2, partName.c_str() );
            //sprintf(buffer2, "Marker %d (%s)", (int)(markerid), mesh->markerName(markerid).substr(0, 32).c_str());
        strcpy( &buffer2[80], "coordinates" );

        MPI_File_write_at(fh, posInFile, buffer2, sizeof(buffer2), MPI_CHAR, &status);
        //MPI_File_write_ordered(fh, buffer2, size, MPI_CHAR, &status);
    }
    posInFile+=160;

    //size_type __nv = mp.ids.size();
    int32_t gnop = (int32_t)(mp.numberOfPointAllProcess( part ) );
    //std::cout << "gnop=" << gnop << std::endl;
    if( this->worldComm().isMasterRank() )
    {
        MPI_File_write_at(fh, posInFile, &gnop, 1, MPI_INT32_T, &status );
        //MPI_File_write_ordered(fh, &gnop, size, MPI_INT32_T, &status );
    }
    posInFile+=sizeOfInt32_t;


#if 0
    std::vector<int32_t> pointids(gnop);
    for(int i = 0 ; i < gnop ; ++i )
        pointids[i] = i+1;
    MPI_File_write_at(fh, posInFile/*+localOffset*/,  pointids.data(), pointids.size(), MPI_INT32_T, &status );
    posInFile+=gnop*sizeOfInt32_t;
#endif

    int32_t nPointOnProcess = mp.numberOfPoint( part,currentPid );
    auto const& nodes_B = mp.nodes( part );
    std::vector<float> nodes( nodes_B.size() );
    CHECK( nodes_B.size() == 3*nPointOnProcess ) << "aie " <<  nodes_B.size() << " vs " <<  3*nPointOnProcess;
    for (int k=0;k<nPointOnProcess;++k )
    {
        nodes[k] = nodes_B[3*k];
        nodes[k+nPointOnProcess] = nodes_B[3*k+1];
        nodes[k+2*nPointOnProcess] = nodes_B[3*k+2];
    }

    int32_t localOffset = mp.startPointIds( part,currentPid )*sizeOfFloat;
    //MPI_File_write_at(fh, posInFile+localOffset,  nodes.data(), nodes.size(), MPI_FLOAT, &status );
    if ( nPointOnProcess > 0 )
    {
        MPI_File_write_at(fh, posInFile+localOffset,  nodes.data(), nPointOnProcess, MPI_FLOAT, &status ); // coord x
        MPI_File_write_at(fh, posInFile+gnop*sizeOfFloat+localOffset,  &nodes.data()[nPointOnProcess], nPointOnProcess, MPI_FLOAT, &status ); // coord y
        MPI_File_write_at(fh, posInFile+2*gnop*sizeOfFloat+localOffset,  &nodes.data()[2*nPointOnProcess], nPointOnProcess, MPI_FLOAT, &status ); // coord z
    }
    posInFile += gnop*3*sizeOfFloat;

    /* write element type */
    if( this->worldComm().isMasterRank() )
    {
        memset(buffer, '\0', sizeof(buffer));
        strncpy( buffer, this->elementType().c_str(), 80-1);
        MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
        //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status );
    }
    posInFile+=80;

    int32_t gnole = mp.numberOfElementAllProcess( part );//mp.globalNumberOfElements();
    //std::cout << "gnole="<<gnole<<std::endl;
    //LOG(INFO) << "Global nb elements: " << gnole << std::endl;
    if( this->worldComm().isMasterRank() )
    {
        MPI_File_write_at(fh, posInFile, &gnole, 1, MPI_INT32_T, &status );
        //MPI_File_write_ordered(fh, &gnole, size, MPI_INT32_T, &status );
    }
    posInFile+=sizeOfInt32_t;
#if 0
    std::vector<int32_t> elids(gnole);
    for (int k=0;k<gnole;++k)
        elids[k] = k+1;
    MPI_File_write_at(fh, posInFile/*+localOffset*/,  elids.data(), elids.size(), MPI_INT32_T, &status );
    posInFile+=gnole*sizeOfInt32_t;
#endif

    auto const& pointsIdsInElt_B = mp.pointIdsInElements( part );
    std::vector<int32_t> pointsIdsInElt(pointsIdsInElt_B.size());

    uint16_type nPointsUsedInElt = mesh_type::element_type::numPoints;
    if constexpr ( N == 1 )
    {
        for ( int k=0;k<pointsIdsInElt_B.size();++k)
            pointsIdsInElt[k] = pointsIdsInElt_B[k] +1 ;
    }
    else
    {
        auto itFindOrdering = M_nodesOrderingInElementToEnsight.find( (boost::format("%1%_%2%d_g%3%") %mesh_type::shape_type::type() %mesh_type::shape_type::nDim %mesh_type::shape_type::nOrder ).str() );
        CHECK( itFindOrdering != M_nodesOrderingInElementToEnsight.end() ) << "not found an ordering";
        auto const& mappingWithThisKindOfElement = itFindOrdering->second;
        nPointsUsedInElt = mappingWithThisKindOfElement.size();
        mp.updateOrderingOfPointsIdsInElt( part, currentPid, pointsIdsInElt, mappingWithThisKindOfElement, 1, nPointsUsedInElt );
    }
    localOffset = mp.startElementIds( part,currentPid )*nPointsUsedInElt*sizeOfInt32_t;
    MPI_File_write_at(fh, posInFile+localOffset, pointsIdsInElt.data(), pointsIdsInElt.size(), MPI_INT32_T, &status );
    posInFile += gnole*nPointsUsedInElt*sizeOfInt32_t;
}
#endif

template<typename MeshType, int N>
void
ExporterEnsightGold<MeshType,N>::writeVariableFiles( timeset_ptrtype __ts, step_ptrtype __step /*steps_write_on_disk_type const& stepsToWriteOnDisk*/ ) const
{
    int nNodalFields = std::distance(__step->beginNodal(), __step->endNodal());
    int nElementFields = std::distance(__step->beginElement(), __step->endElement());
    if ( (nNodalFields+nElementFields) == 0 )
        return;//continue;
    LOG(INFO) << "ExporterEnsightGold::writeVariableFiles Nodal: " << nNodalFields;
    LOG(INFO) << "ExporterEnsightGold::writeVariableFiles Element: " << nElementFields;

    // prepare some common information
    auto __firstActiveStep = __ts->firstActiveStep();
    bool isFirstStep = (__step->activeIndex() == __firstActiveStep->activeIndex());
    std::ostringstream ossFilenameStepIndex;
    bool writeNewFile = true;
    if( M_mergeTimeSteps )
    {
        if ( M_packTimeSteps > 0 )
        {
            int relativeTsIndex = __step->index() - (*__ts->beginStep())->index();
            ossFilenameStepIndex << (relativeTsIndex/M_packTimeSteps + 1);
            writeNewFile = ( ( relativeTsIndex % M_packTimeSteps ) == 0 ) ;
        }
        else
            writeNewFile = isFirstStep;
    }
    else
    {
        ossFilenameStepIndex << std::setfill( '0' ) << std::setw( M_timeExponent ) << __step->activeIndex();
    }
    std::string filenameStepIndex = ossFilenameStepIndex.str();

    saveFields<true>( __ts, __step, writeNewFile, filenameStepIndex, isFirstStep, __step->beginNodal(), __step->endNodal() );
    saveFields<false>( __ts, __step, writeNewFile, filenameStepIndex, isFirstStep, __step->beginElement(), __step->endElement() );
}


template<typename MeshType, int N>
template<bool IsNodal,typename Iterator>
void
ExporterEnsightGold<MeshType,N>::saveFields( timeset_ptrtype __ts, typename timeset_type::step_ptrtype __step, bool writeNewFile, std::string const& filenameStepIndex, bool isFirstStep, Iterator __var, Iterator en ) const
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

    const int reorder_tensor2symm[6] = { 0,3,4,1,5,2 };

    auto __mesh = __step->mesh();

    auto itFindCache = M_cache_mp.find( __ts->name() );
    CHECK( itFindCache != M_cache_mp.end() ) << "mesh numbering cache not found";
    auto const& mp = *itFindCache->second;
    CHECK( __mesh->isSameMesh( mp.mesh() ) ) << "something wrong mesh should be identical between step and cache";
    rank_type currentPid = __mesh->worldComm().localRank();

    while ( __var != en )
    {
        tic();
        auto const& fieldData =  __var->second;
        auto const& field00 = unwrap_ptr( fieldData.second[0][0] );

        if ( !field00.worldComm().isActive() ) return;

        uint16_type nComponents = invalid_uint16_type_value, nComponents1 = invalid_uint16_type_value, nComponents2 = invalid_uint16_type_value;
        std::string fileExt;
        bool isTensor2Symm = false;
        if ( fieldData.first == FunctionSpaceType::SCALAR )
        {
            nComponents = 1;
            nComponents1 = 1;
            nComponents2 = 1;
            fileExt = ".scl";
        }
        else if ( fieldData.first == FunctionSpaceType::VECTORIAL )
        {
            nComponents = 3;
            nComponents1 = 3;
            nComponents2 = 1;
            fileExt = ".vec";
        }
        else if ( fieldData.first == FunctionSpaceType::TENSOR2 )
        {
            nComponents = 9;
            nComponents1 = 3;
            nComponents2 = 3;
            fileExt = ".tsr";
        }
        else if ( fieldData.first == FunctionSpaceType::TENSOR2_SYMM )
        {
            nComponents = 6;
            nComponents1 = 3;
            nComponents2 = 3;
            fileExt = ".tsrs";
            isTensor2Symm = true;
        }


        std::ostringstream __varfname;

        __varfname << this->path() << "/" << __ts->name() << "." << __var->first;
        /* if we want to pack data in several files instead of one */
        /* we compute an index to add to the filename */
        if( M_mergeTimeSteps && M_packTimeSteps > 0 )
            __varfname << "." << filenameStepIndex;
        // add extension
        __varfname << fileExt;
        if( !M_mergeTimeSteps )
            __varfname << "." << filenameStepIndex;

        DVLOG(2) << "[ExporterEnsightGold::saveFields] saving " << __varfname.str() << "...\n";
        std::fstream __out;

        /* Open File with MPI IO */
        char * str = strdup(__varfname.str().c_str());

        /* Check if file exists if we are on step one and delete it if so */
        /* (MPI IO does not have a truncate mode ) */
        // std::cout << "Nodes " << this->worldComm().isMasterRank() << " " << __step->index() << " " << isFirstStep << " " << fs::exists(str) << std::endl;
#if 1
        if ( writeNewFile )
        {
            if( this->worldComm().isMasterRank()  && fs::exists(str))
            {
                MPI_File_delete(str, MPI_INFO_NULL);
            }
            MPI_Barrier(this->worldComm().comm());
        }
#endif

        //init MPI_Info object fromhints defined as environment variables
        //MPI_Info info = initIoInfoFromEnvVars();
        if ( M_mergeTimeSteps && !writeNewFile )
            MPI_File_open( this->worldComm().comm(), str, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL , &fh );
        else
            MPI_File_open( this->worldComm().comm(), str, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL , &fh );

        //MPI_Info_free(&info);
        free(str);

        // INIT CURSOR IN FILE
        posInFile = 0;
#if 0
        // trunc the file ( not works always, some part are not good sometimes ??)
        if ( writeNewFile )
            MPI_File_set_size( fh, posInFile );
#endif
        Feel::detail::FileIndex index( this->worldCommPtr()  );

        if( M_mergeTimeSteps )
        {
            // first read
            index.read(fh);

            /* Move to the beginning of the fie index section */
            /* to overwrite it */
            if ( index.defined() && (__step->index() - TS_INITIAL_INDEX) > 0 ) {
                // MPI_File_seek_shared(fh, index.fileblock_n_steps, MPI_SEEK_SET);
                posInFile = index.nextFreePosFile();
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
                MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
            }
            posInFile+=80;

            /* add the beginning of the new block to the file */
            // MPI_File_get_position_shared(fh, &offset);
            index.add( posInFile );
        }

        if( this->worldComm().isMasterRank() )
        {
            memset(buffer, '\0', sizeof(buffer));
            strncpy(buffer, field00.name().c_str(), 80-1);
            MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
        }
        posInFile+=80;

        toc("saveFields intro",FLAGS_v>0);
        /* handle faces data */
#if 0
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
                            if ( c1 >= fieldData.second.size() )
                                continue;
                            for ( uint16_type c2 = 0; c2 < nComponents2; ++c2 )
                            {
                                if ( c2 >= fieldData.second[c1].size() )
                                    continue;
                                auto const& fieldComp = unwrap_ptr( fieldData.second[c1][c2] );
                                uint16_type c = c2*nComponents1+c1;
                                size_type global_node_id = nfaces*c + pid ;
                                size_type thedof =  fieldComp.start() +
                                    fieldComp.functionSpace()->dof()->faceLocalToGlobal( face.id(), j, c ).index();
                                field[global_node_id] = fieldComp.globalValue( thedof );
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
                for ( uint16_type c = 0; c < field00.nComponents; ++c )
                {
                    MPI_File_write_at(fh, posInFile + localOffset + c*sumOffsets, \
                                      field.data() + nfaces * c, nfaces, MPI_FLOAT, &status);
                    // MPI_File_write_ordered(fh, field.data() + nfaces * c, nfaces, MPI_FLOAT, &status);
                }
                posInFile += field00.nComponents * sumOffsets;
            } // boundaries loop
        }
#endif

        /* handle elements */
        for ( auto const& [part,nameAndRangeElt] : mp.partIdToRangeElement() )
        {
            tic();
            if( this->worldComm().isMasterRank() )
            {
                memset(buffer, '\0', sizeof(buffer));
                strcpy( buffer, "part" );
                //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
                MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
                //int32_t partid = *mit + 1;
                int32_t partid = part + 1;
                MPI_File_write_at(fh, posInFile+80, &partid, 1, MPI_INT32_T, &status);
                DVLOG(2) << "part " << buffer << "\n";
                if ( IsNodal )
                {
                    strcpy(buffer, "coordinates");
                }
                else
                {
                    memset(buffer, '\0', sizeof(buffer));
                    strncpy(buffer, this->elementType().c_str(), 80-1);
                }
                MPI_File_write_at(fh, posInFile+80+sizeOfInt32_t, buffer, sizeof(buffer), MPI_CHAR, &status);
            }
            posInFile+=160+sizeOfInt32_t;
            toc("saveFields part",FLAGS_v>0);


            /* create an array to store data per node */
            index_type nValuesPerComponent = invalid_v<index_type>;
            index_type nValuesPerComponentAllProcess = invalid_v<index_type>;
            index_type offsetValuesPerComponent = invalid_v<index_type>;

            Eigen::VectorXf __field;

            /* Loop on the elements */
            tic();
            auto const& d = field00.functionSpace()->dof().get();

            if constexpr ( IsNodal )
                {

                    nValuesPerComponentAllProcess = mp.numberOfPointAllProcess( part );
                    offsetValuesPerComponent = mp.startPointIds( part,currentPid );

                    if ( M_mapNodalArrayToDofId.find( part ) == M_mapNodalArrayToDofId.end() )
                    {
                        nValuesPerComponent = mp.numberOfPoint( part,currentPid );

                        auto & mapArrayToDofId = M_mapNodalArrayToDofId[part];
                        mapArrayToDofId.resize( nValuesPerComponent,invalid_size_type_value );

                        size_type __field_size = nValuesPerComponent*nComponents;
                        VLOG(1) << "field size=" << __field_size;
                        __field = Eigen::VectorXf::Zero( __field_size );

                        //const int np = __step->mesh()->numLocalVertices();

                        auto const& r =  mp.rangeElement( part );
                        auto elt_it = r.begin();
                        auto elt_en = r.end();
                        for ( ; elt_it != elt_en; ++elt_it )
                        {
                            auto const& elt = unwrap_ref( *elt_it );
                            auto const& locglob_ind = d->localToGlobalIndices( elt.id() );
                            const int np = N==1? elt.numVertices : elt.numPoints;
                            for ( uint16_type p = 0; p < np; ++p )
                            {
                                index_type ptid = mp.pointIdToContiguous(part,elt.point( p ).id());
                                if ( ptid == invalid_v<typename mesh_type::index_type> ) // point not in this process
                                    continue;
                                ptid -= offsetValuesPerComponent;
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

                                        size_type global_node_id = nValuesPerComponent*cMap + ptid ;
                                        DCHECK( ptid < __step->mesh()->numPoints() ) << "Invalid point id " << ptid << " element: " << elt.id()
                                                                                     << " local pt:" << p
                                                                                     << " mesh numPoints: " << __step->mesh()->numPoints();
                                        DCHECK( global_node_id < __field_size ) << "Invalid dof id : " << global_node_id << " max size : " << __field_size;

                                        __field(global_node_id) = fieldComp.globalValue( dof_id );
                                        DVLOG(3) << "v[" << global_node_id << "]=" << fieldComp.globalValue( dof_id ) << "  dof_id:" << dof_id;
                                    }
                                }
                            }
                        }
                    }
                } // IsNodal
            else
            {
                nValuesPerComponentAllProcess = mp.numberOfElementAllProcess( part );
                offsetValuesPerComponent = mp.startElementIds( part,currentPid );

                if ( M_mapElementArrayToDofId.find( part ) == M_mapElementArrayToDofId.end() )
                {
                    nValuesPerComponent = mp.numberOfElement( part,currentPid );

                    auto & mapArrayToDofId =  M_mapElementArrayToDofId[part];
                    mapArrayToDofId.resize( nValuesPerComponent,invalid_size_type_value );

                    size_type __field_size = nValuesPerComponent*nComponents;
                    VLOG(1) << "field size=" << __field_size;
                    __field = Eigen::VectorXf::Zero( __field_size );
                    auto const& r =  mp.rangeElement( part );
                    auto elt_it = r.begin();
                    auto elt_en = r.end();
                    for ( ; elt_it != elt_en; ++elt_it )
                    {
                        auto const& elt = unwrap_ref( *elt_it );
                        index_type e = mp.elementIdToContiguous(part,elt.id());
                        e -= offsetValuesPerComponent;
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

                                size_type global_node_id = cMap*nValuesPerComponent + e;
                                __field(global_node_id) = fieldComp.globalValue( dof_id );
                                DVLOG(3) << "v[" << global_node_id << "]=" << fieldComp.globalValue( dof_id ) << "  dof_id:" << dof_id;
                            }
                        }
                    }
                }
            }

            if ( nValuesPerComponent == invalid_v<index_type> )
            {
                auto const& mapArrayToDofId = (IsNodal)? M_mapNodalArrayToDofId.find( part )->second : M_mapElementArrayToDofId.find( part )->second;
                nValuesPerComponent = mapArrayToDofId.size();
                size_type __field_size = nValuesPerComponent*nComponents;
                VLOG(1) << "field size=" << __field_size;
                __field = Eigen::VectorXf::Zero( __field_size );
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
                            size_type global_node_id = cMap*nValuesPerComponent + k;
                            __field(global_node_id) = fieldComp.globalValue( dof_id );
                        }
                    }
                }
            }

            toc("saveFields element loop",FLAGS_v>0);
            tic();
            // - write in file on cursor : posInFile + localOffset
            if ( nValuesPerComponent > 0 )
            {
                for ( uint16_type c = 0; c < nComponents; ++c )
                {
                    MPI_File_write_at(fh, posInFile + (offsetValuesPerComponent + c*nValuesPerComponentAllProcess)*sizeOfFloat,
                                      &__field.data()[c*nValuesPerComponent], nValuesPerComponent, MPI_FLOAT, &status);
                }
            }
            posInFile += nComponents*nValuesPerComponentAllProcess*sizeOfFloat;
            toc("saveFields write part",FLAGS_v>0);
        } // parts loop

        if ( M_mergeTimeSteps )
        {
            /* write timestep end */
            if( this->worldComm().isMasterRank() )
            {
                memset(buffer, '\0', sizeof(buffer));
                strncpy(buffer,"END TIME STEP", 80);
                MPI_File_write_at(fh, posInFile, buffer, sizeof(buffer), MPI_CHAR, &status);
            }
            posInFile+=80;
            // write back the file index
            index.write( fh, posInFile );
        }
        DVLOG(2) << "[ExporterEnsightGold::saveFields] saving " << __varfname.str() << "done\n";

        MPI_File_close(&fh);

        ++__var;
    }
    toc("saveFields", FLAGS_v>0);
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
