/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-11-09

  Copyright (C) 2004,2005 EPFL
  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)

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
/*!
 * \file exporterhdf5_impl.hpp
 * \brief HDF5 and XDMF exporter
 * \author Benjamin Vanthong <benjamin.vanthong@gmail.com>
 * \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
 * \date 2014-08-28
 */
#ifndef __Exporterhdf5_CPP
#define __Exporterhdf5_CPP 1

#if defined(FEELPP_HAS_HDF5)

#include <feel/feelcore/feel.hpp>

#include <feel/feelfilters/exporterhdf5.hpp>
#include <feel/feelcore/hdf5.hpp>

namespace Feel 
{
namespace fs = boost::filesystem;

template<typename MeshType, int N>
Exporterhdf5<MeshType,N>::Exporterhdf5( WorldComm const& worldComm )
:
super( worldComm ),
M_element_type()

{
    init();
}
template<typename MeshType, int N>
Exporterhdf5<MeshType,N>::Exporterhdf5( std::string const& __p, int freq, WorldComm const& worldComm )
    :
    super( "hdf5", __p, freq, worldComm ),
    M_element_type()
{
    init();
}
template<typename MeshType, int N>
Exporterhdf5<MeshType,N>::Exporterhdf5( po::variables_map const& vm, std::string const& exp_prefix, WorldComm const& worldComm )
    :
    super( vm, exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
Exporterhdf5<MeshType,N>::Exporterhdf5( std::string const& exp_prefix, WorldComm const& worldComm )
    :
    super( exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
Exporterhdf5<MeshType,N>::Exporterhdf5( Exporterhdf5 const & __ex )
    :
    super( __ex ),
    M_element_type( __ex.M_element_type )
{
}

template<typename MeshType, int N>
Exporterhdf5<MeshType,N>::~Exporterhdf5()
{}


template<typename MeshType, int N>
void
Exporterhdf5<MeshType,N>::init()
{
    if ( mesh_type::nDim == 1 )
        if ( mesh_type::Shape == SHAPE_LINE )
            M_element_type = ( mesh_type::nOrder == 1 )?"Polyline":"Edge_3";

    if ( mesh_type::nDim == 2 )
    {
        if ( mesh_type::Shape == SHAPE_TRIANGLE )
            M_element_type = ( mesh_type::nOrder == 1 )?"Triangle":"Tri_6";

        else if ( mesh_type::Shape == SHAPE_QUAD )
            M_element_type = ( mesh_type::nOrder == 1 )?"Quadrilateral":"Quad_8";
    }

    if ( mesh_type::nDim == 3 )
    {
        if ( mesh_type::Shape == SHAPE_TETRA )
            M_element_type = ( mesh_type::nOrder == 1 )?"Tetrahedron":"Tet_10";

        else if ( mesh_type::Shape == SHAPE_HEXA )
            M_element_type = ( mesh_type::nOrder == 1 )?"Hexahedron":"Hex_20";
    }
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::save () const 
{
    /* make sure to reset values from previous calls */
    M_XDMFContent.str("");

    /* build file name */
    M_fileName.str("");
    M_fileName << this->prefix();
    //std::cout << this->prefix() << std::endl;

    /* First write hdf5 files */
    writeHDF5();

    /* Then the Xdmf description file */
    /* This must be written after the HDF5 file */
    /* For the XDMF content to be built on each processor in the case we are using MPI IO */
    writeXDMF();
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeXDMF() const 
{
    int size;
    int stepIndex = TS_INITIAL_INDEX;
    std::ostringstream cbuf;
    MPI_File fh;
    MPI_Status status;

    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();

    timeset_ptrtype __ts = *__ts_it;

    while ( __ts_it != __ts_en )
    {
        __ts = *__ts_it;
        typename timeset_type::step_const_iterator __it = __ts->beginStep();
        typename timeset_type::step_const_iterator __end = __ts->endStep();

        typename timeset_type::step_ptrtype __step;

        /* If we have a timestep to save, we replace the defaulted step index by the corresponding index */
        if(__it != __end)
        {
            __it = boost::prior ( __end );
            __step = *__it;

            stepIndex = __step->index();
        }

        if( __it == __end || ( __it != __end && __step->isInMemory() ) )
        {
            /* build file name */
            cbuf << M_fileName.str() << "-" << stepIndex << ".xmf";

            /* Test master rank only once */
            bool isMasterRank = this->worldComm().isMasterRank();

            /* Open file with MPI IO */
            char * strTmp = strdup(cbuf.str().c_str());
            if(isMasterRank && fs::exists(strTmp))
            {
                MPI_File_delete(strTmp, MPI_INFO_NULL);
            }

            MPI_Barrier(this->worldComm().comm());
            MPI_File_open(this->worldComm().comm(), strTmp, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

            free(strTmp);

            /* write file header */
            cbuf.str("");
            cbuf << "<?xml version=\"1.0\" ?>" << std::endl;
            cbuf << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
            cbuf << "<Xdmf Version=\"2.0\">" << std::endl;
            cbuf << "<Domain>" << std::endl;
            cbuf << "<Grid Name=\"Simulation\" GridType=\"Collection\" CollectionType=\"Spatial\">" << std::endl;

            if( isMasterRank ) 
            { size = cbuf.str().size(); }
            else
            { size = 0; }
            strTmp = strdup(cbuf.str().c_str());
            MPI_File_write_ordered(fh, strTmp, size, MPI_CHAR, &status);
            free(strTmp);

            /* Write time info */
            /*
               timeset_const_iterator __ts_it = this->beginTimeSet ();
               timeset_const_iterator __ts_en = this->endTimeSet ();
               timeset_ptrtype __ts = *__ts_it;

               cbuf.str("");
               cbuf << "       <Grid Name=\"Simulation over time\" GridType=\"Collection\" CollectionType=\"Temporal\">" << "" << std::endl;
               cbuf << "           <Time TimeType=\"HyperSlab\">" << std::endl;
               M_xmf << "               <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\">" << std::endl;
               cbuf << "               " << (*(__ts->beginStep()))->time() << " " << __ts->timeIncrement() << " " << __ts->numberOfTotalSteps() << std::endl;
               cbuf << "               </DataItem>" << std::endl;
               cbuf << "           </Time>" << std::endl;

               if( isMasterRank ) 
               { size = cbuf.str().size(); }
               else
               { size = 0; }
               MPI_File_write_ordered(fh, cbuf.str().c_str(), size, MPI_CHAR, &status);
               */

            /* write Xdmf content */
            strTmp = strdup(M_XDMFContent.str().c_str());
            MPI_File_write_ordered(fh, strTmp, M_XDMFContent.str().size(), MPI_CHAR, &status);
            free(strTmp);

            /* write footer */
            cbuf.str("");
            cbuf << "</Grid>" << std::endl;
            cbuf << "</Domain>" << std::endl;
            cbuf << "</Xdmf>" << std::endl;

            if( isMasterRank ) 
            { size = cbuf.str().size(); }
            else
            { size = 0; }
            strTmp = strdup(cbuf.str().c_str());
            MPI_File_write_ordered(fh, strTmp, size, MPI_CHAR, &status);
            free(strTmp);

            MPI_File_close(&fh);
        }
        ++__ts_it;
    }
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeHDF5() const 
{
    int size;
    int stepIndex = TS_INITIAL_INDEX;
    std::ostringstream oss;

    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();
    
    timeset_ptrtype __ts = *__ts_it;

    while ( __ts_it != __ts_en )
    {
        __ts = *__ts_it;
        typename timeset_type::step_const_iterator __it = __ts->beginStep();
        typename timeset_type::step_const_iterator __end = __ts->endStep();
        typename timeset_type::step_ptrtype __step;

        /* If we have a timestep to save, we replace the defaulted step index by the corresponding index */
        if(__it != __end)
        {
            __step = *(boost::prior( __end ));
            stepIndex = __step->index();
        }

        oss.str("");
        oss << M_fileName.str() << "-" << stepIndex << ".h5";

        /* open the corresponding h5 file for output */
        M_HDF5.openFile(oss.str(), this->worldComm(), false);

        /* TODO find a way to incorporate time steps directly info the xml data */
        /*
           std::ostrinstream str
           str.str("");
           str << "           <Grid Name=\"" << M_fileNameStep << "\" GridType=\"Uniform\">" << std::endl;
           str << "               <Time Value=\"" << __step->time() << "\"/>" << std::endl;  
           */

        M_XDMFContent << "<Grid Name=\"" << oss.str() << "\" GridType=\"Uniform\">" << std::endl;

        /* check if we have steps for the current dataset */
        if(__it == __end)
        {
            LOG(INFO) << "Timeset " << __ts->name() << " (" << __ts->index() << ") contains no timesteps (Consider using add() or addRegions())" << std::endl;

            /* if we have a mesh, we write mesh data */
            if(__ts->hasMesh())
            {
                saveMesh(__ts->mesh(), stepIndex);
            }
        }
        else
        {
            if ( __step->isInMemory() )
            {
                saveMesh(__step->mesh(), stepIndex);

                /*
                if ( this->worldComm().globalRank() == this->worldComm().masterRank() )
                {
                    std::cout << "time                          : " << __step->time () << std::endl; 
                    std::cout << "time increment                : " << __ts->timeIncrement () << std::endl; 
                    std::cout << "numberOfSteps                 : " << __ts->numberOfSteps () << std::endl; 
                    std::cout << "numberOfTotalSteps            : " << __ts->numberOfTotalSteps () << std::endl;
                    std::cout << "file generated                : " << M_fileNameStep  << ".h5" << std::endl;
                }
                */

                saveNodal(__step, __step->beginNodalScalar(), __step->endNodalScalar() );
                saveNodal(__step, __step->beginNodalVector(), __step->endNodalVector() );
                saveNodal( __step, __step->beginNodalTensor2(), __step->endNodalTensor2() );

                saveElement(__step, __step->beginElementScalar(), __step->endElementScalar() );
                saveElement(__step, __step->beginElementVector(), __step->endElementVector() );
                saveElement( __step, __step->beginElementTensor2(), __step->endElementTensor2() );
            }
        }

        M_XDMFContent << "</Grid>" << std::endl;

        // TODO Check if it must be done here
        M_HDF5.closeFile();

        ++__ts_it;
    }
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::visit ( mesh_type* mesh) 
{
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::saveMesh(mesh_ptrtype mesh, int stepIndex) const 
{
    std::ostringstream oss;
    std::ostringstream groupName;

    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];

    typename mesh_type::parts_const_iterator_type p_st = mesh->beginParts();
    typename mesh_type::parts_const_iterator_type p_en = mesh->endParts();
    for(auto p_it = p_st ; p_it != p_en; ++p_it )
    {
        auto r = markedelements(mesh, p_it->first, EntityProcessType::LOCAL_ONLY );
        auto lelt_it = r.template get<1>();
        auto lelt_en = r.template get<2>();
        Feel::detail::MeshPoints<float> mp( mesh.get(), this->worldComm(), lelt_it, lelt_en, false, true, true, 0 );

        auto elt = boost::unwrap_ref( *lelt_it );

        // build up group name
        groupName << "/" << this->worldComm().globalRank() << "-" << p_it->first;

        currentCount[0] = mp.ids.size();
        currentCount[1] = 3;

        /* if we are using MPI IO all the processes must create the tables */
        /* for all processes of we end up in a deadlock (see HDF5 FAQ) */
        for (size_type i = 0; i < this->worldComm().globalSize (); i++)
        {
            std::ostringstream gName;
            gName << "/" << i << "-" << p_it->first;
            currentCount[0] = mp.numberOfPoints[i];
            //M_HDF5.createTable(groupName.str(), "point_coords", H5T_NATIVE_DOUBLE, currentCount, false);
            M_HDF5.createGroup(gName.str());
            M_HDF5.createTable(gName.str(), "point_coords", H5T_NATIVE_FLOAT, currentCount, true);
        }
        currentCount[0] = mp.ids.size(); 

        hsize_t currentOffset[2] = {0, 0};

        /* write the point coordinates */
        M_HDF5.write(groupName.str() + "point_coords", H5T_NATIVE_FLOAT, currentCount, currentOffset, mp.coords.data());

        /* close all the opened tables */
        for (size_type i = 0; i < this->worldComm().globalSize(); i++)
        {
            std::ostringstream gName;
            gName << "/" << i << "-" << p_it->first;
            M_HDF5.closeTable(gName.str()+"point_coords");
        }

        M_XDMFContent << "<Geometry GeometryType=\"XYZ\">" << std::endl;
        M_XDMFContent << "<DataItem Dimensions=\"" << mp.ids.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" Endian=\"Little\">" << std::endl;
        M_XDMFContent << M_fileName.str() << "-" << stepIndex << ".h5:" << groupName.str() << "/point_coords" << std::endl;

        M_XDMFContent << "</DataItem>" << std::endl;
        M_XDMFContent << "</Geometry>" << std::endl;

        currentSpaceDims [0] = mp.numberOfElements[this->worldComm().globalRank()];
        currentSpaceDims [1] = elt.numLocalVertices;

        for (size_type i = 0; i < this->worldComm().globalSize (); i++)
        {
            std::ostringstream gName;
            gName << "/" << i << "-" << p_it->first;
            currentSpaceDims [0] = mp.numberOfElements[i];

            M_HDF5.createTable (gName.str(), "element_nodes", H5T_STD_I32LE, currentSpaceDims, true);
        }
        currentSpaceDims [0] = mp.numberOfElements[this->worldComm().globalRank()];

        M_HDF5.write( groupName.str()+"element_nodes", /*H5T_NATIVE_LLONG*/ H5T_STD_I32LE /*H5T_NATIVE_B32*/, currentSpaceDims, currentOffset, mp.elem.data() );

        for (size_type i = 0; i < this->worldComm().globalSize(); i++)
        {
            std::ostringstream gName;
            gName << "/" << i << "-" << p_it->first;
            M_HDF5.closeTable(gName.str()+"element_nodes");
        }

        M_XDMFContent << "<Topology TopologyType=\"" << M_element_type << "\" NumberOfElements=\"" << mp.numberOfElements[this->worldComm().globalRank()] << "\" NodesPerElement=\"" << elt.numLocalVertices << "\">" << std::endl;
        M_XDMFContent << "<DataItem Dimensions=\"" << mp.numberOfElements[this->worldComm().globalRank()] << " " << elt.numLocalVertices << "\" NumberType=\"Int\" Precision=\"4\" Format=\"HDF\" Endian=\"Little\">" << std::endl;
        M_XDMFContent << M_fileName.str() << "-" << stepIndex << ".h5:" << groupName.str() << "/element_nodes" << std::endl;

        M_XDMFContent << "</DataItem>" << std::endl;
        M_XDMFContent << "</Topology>" << std::endl;
    }
}

template <typename MeshType, int N>
template <typename Iterator>
void Exporterhdf5<MeshType, N>::saveNodal ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const 
{   
    while ( __var != en )
    {    
        std::string attributeType = "Scalar";
        std::string solutionName = __var->first;
        uint16_type nComponents = __var->second.nComponents;

        if ( __var->second.is_scalar )
        {
            solutionName += ".scl"; 
            attributeType = "Scalar";
            nComponents = 1;
        }
        else if ( __var->second.is_vectorial )
        {
            solutionName += ".vec";
            attributeType = "Vector";
            nComponents = 3;
        }
        else if ( __var->second.is_tensor2 )
        {
            solutionName += ".tsr";
            attributeType = "Tensor";
            nComponents = 9;
        }

        solutionName += ".node";

        hsize_t currentSpacesDims [2];

        typename mesh_type::parts_const_iterator_type p_it = __step->mesh()->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = __step->mesh()->endParts();

        std::ostringstream groupName;
        groupName << "/" << this->worldComm().globalRank() << "-" << p_it->first;

        for (; p_it != p_en; p_it++) 
        {
            auto r = markedelements(__step->mesh(), p_it->first, EntityProcessType::LOCAL_ONLY);
            auto elt_it = r.template get<1>();
            auto elt_en = r.template get<2>();

            Feel::detail::MeshPoints<float> mp ( __step->mesh().get(), this->worldComm(), elt_it, elt_en, false, true, true, 0 );

            currentSpacesDims[0] = mp.ids.size();
            currentSpacesDims[1] = nComponents;

            for (size_type i = 0; i < this->worldComm().globalSize (); i++)
            {
                std::ostringstream gName;
                gName << "/" << i << "-" << p_it->first;
                currentSpacesDims[0] = mp.numberOfPoints[i];

                M_HDF5.createTable(gName.str(), solutionName.c_str() , H5T_NATIVE_FLOAT, currentSpacesDims, true);
            }
            
            currentSpacesDims[0] = mp.ids.size();

            std::vector<float> realBuffer;
            realBuffer.resize(currentSpacesDims[0] * currentSpacesDims[1], 0);

            size_type e = 0; 

            for (; elt_it != elt_en; ++elt_it )
            {
                for ( uint16_type c = 0; c < nComponents; ++c )
                {
                    for ( uint16_type p = 0; p < __step->mesh()->numLocalVertices(); ++p, ++e )
                    {
                        size_type ptid = mp.old2new[elt_it->get().point(p).id()];
                        size_type global_node_id = ptid * nComponents + c;

                        if ( c < __var->second.nComponents ) 
                        {
                            size_type dof_id = boost::get<0>( __var->second.functionSpace()->dof()->localToGlobal ( elt_it->get().id(), p, c ) );
                            realBuffer[global_node_id] = __var->second.globalValue ( dof_id );
                        }
                        else
                        {
                            realBuffer[global_node_id] = 0.0;
                        }
                    }
                }
            }

            hsize_t currentOffset[2] = {0, 0};

            M_HDF5.write( groupName.str() + solutionName, H5T_NATIVE_FLOAT, currentSpacesDims, currentOffset, &realBuffer[0] );

            for (size_type i = 0; i < this->worldComm().globalSize(); i++)
            {
                std::ostringstream gName;
                gName << "/" << i << "-" << p_it->first;
                M_HDF5.closeTable(gName.str()+solutionName);
            }

            M_XDMFContent << "<Attribute AttributeType=\""<< attributeType << "\" Name=\"" << solutionName << "\" Center=\"Node\">" << std::endl;
            M_XDMFContent << "<DataItem Dimensions=\""<< mp.ids.size() << " " << nComponents << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" Endian=\"Little\">" << std::endl;
            M_XDMFContent << M_fileName.str() << "-" << __step->index() << ".h5:" << groupName.str() << "/" << solutionName << std::endl;

            M_XDMFContent << "</DataItem>" << std::endl;    
            M_XDMFContent << "</Attribute>" << std::endl;
        }    
        ++__var;
    }
}

template<typename MeshType, int N>
template<typename Iterator>
void Exporterhdf5<MeshType, N>::saveElement ( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const 
{
    while ( __evar != __evaren ) 
    {
        std::string attributeType ("Scalar");
        std::string solutionName = __evar->first;
        uint16_type nComponents = __evar->second.nComponents;

        if ( __evar->second.is_scalar )
        {
            solutionName += ".scl"; 
            attributeType = "Scalar";
            nComponents = 1;
        }
        else if ( __evar->second.is_vectorial )
        {
            solutionName += ".vec";
            attributeType = "Vector";
            nComponents = 3;
        }
        else if ( __evar->second.is_tensor2 )
        {
            solutionName += ".tsr";
            attributeType = "Tensor";
            nComponents = 9;
        }
        solutionName += ".element";
        hsize_t currentSpacesDims [2];

        typename mesh_type::parts_const_iterator_type p_it = __step->mesh()->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = __step->mesh()->endParts();
        for (; p_it != p_en; p_it++ ) 
        {
            typename mesh_type::marker_element_const_iterator elt_st;
            typename mesh_type::marker_element_const_iterator elt_en;
            boost::tie( elt_st, elt_en ) = __step->mesh()->elementsWithMarker( p_it->first, __evar->second.worldComm().localRank() );

            std::ostringstream groupName;
            groupName << "/" << this->worldComm().globalRank() << "-" << p_it->first;

            Feel::detail::MeshPoints<float> mp ( __step->mesh().get(), this->worldComm(), elt_st, elt_en, false, true, true, 0 );

            currentSpacesDims[0] = mp.numberOfElements[this->worldComm().globalRank()];
            currentSpacesDims[1] = nComponents;

            for (size_type i = 0; i < this->worldComm().globalSize (); i++)
            {
                std::ostringstream gName;
                gName << "/" << i << "-" << p_it->first;
                currentSpacesDims[0] = mp.numberOfElements[i];

                M_HDF5.createTable(gName.str(), solutionName, H5T_NATIVE_FLOAT, currentSpacesDims, true);
            }
            currentSpacesDims[0] = mp.numberOfElements[this->worldComm().globalRank()];

            std::vector<float> realBuffer;
            realBuffer.resize(currentSpacesDims[0] * currentSpacesDims[1], 0);

            if ( !__evar->second.areGlobalValuesUpdated() )
                __evar->second.updateGlobalValues();

            size_type ncells = std::distance ( elt_st, elt_en );
            size_type e = 0;
            for ( auto elt_it = elt_st; elt_it != elt_en; ++elt_it, ++e )
            {
                for ( int c = 0; c < nComponents; ++c )
                {
                    size_type global_node_id = e * nComponents + c;
                    if ( c < __evar->second.nComponents )
                    {
                        size_type dof_id = boost::get<0>( __evar->second.functionSpace()->dof()->localToGlobal( elt_it->id(), 0, c ) );
                        realBuffer[global_node_id] = __evar->second.globalValue ( dof_id );
                    }
                    else 
                    {
                        realBuffer[global_node_id] = 0;
                    }
                }
            }

            hsize_t currentOffset[2] = {0, 0};

            M_HDF5.write(groupName.str() + solutionName, H5T_NATIVE_FLOAT, currentSpacesDims, currentOffset, &realBuffer[0]);

            for (size_type i = 0; i < this->worldComm().globalSize(); i++)
            {
                std::ostringstream gName;
                gName << "/" << i << "-" << p_it->first;
                M_HDF5.closeTable(gName.str() + solutionName);
            }

            M_XDMFContent << "<Attribute AttributeType=\"" << attributeType << "\" Name=\"" << solutionName << "\" Center=\"Cell\">" << std::endl;
            M_XDMFContent << "<DataItem Dimensions=\""<< mp.numberOfElements[this->worldComm().globalRank()] << " " << nComponents << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" Endian=\"Little\">" << std::endl;
            M_XDMFContent << M_fileName.str() << "-" << __step->index() << ".h5:" << groupName.str() << "/" << solutionName << std::endl;

            M_XDMFContent << "</DataItem>" << std::endl;    
            M_XDMFContent << "</Attribute>" << std::endl;
        }   
        __evar++;
    }
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeStats() const 
{
    size_type numStats = 3 ;
    hsize_t currentSpacesDims[2] = {1, numStats} ;

    M_HDF5.createTable( "stats", H5T_STD_U32LE, currentSpacesDims ) ;

    std::vector<size_type> uintBuffer;
    uintBuffer.resize (numStats) ;

    uintBuffer[0] = 0; //M_maxNumPoints
    uintBuffer[1] = 0; //M_maxNumElements
    uintBuffer[2] = 0; //M_elementNodes

    if ( this->worldComm().globalRank() == this->worldComm().masterRank() )
    {
        /*
        std::cout << "nombre de Points              : " << M_maxNumPoints << std::endl ;
        std::cout << "M_numMaxElements              : " << M_maxNumElements << std::endl ;
        std::cout << "nombre de Points par element  : " << M_elementNodes << std::endl ;
        std::cout << "M_numParts                    : " << M_numParts << std::endl ;
        */
        std::cout << "mesh_type::nRealDim           : " << mesh_type::nRealDim << std::endl ;
        std::cout << "fileName                      : " << M_fileName.str() << ".xmf" << std::endl ;
    }

    hsize_t currentOffset [2] = {0, 0} ;
    M_HDF5.write("stats", H5T_NATIVE_LLONG, currentSpacesDims, currentOffset, &uintBuffer[0]) ;
    M_HDF5.closeTable("stats") ;
}

}
#endif /* FEELL_HAS_HDF5 */
#endif /* __Exporterhdf5_CPP */

