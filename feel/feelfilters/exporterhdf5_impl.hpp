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
 * \author VANTHONG Benjamin <benjamin.vanthong@gmail.com>
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
    /*
    if ( this->worldComm().globalRank() == this->worldComm().masterRank() )
        std::cout << "exporter.merge                : " << (boption (_name = "exporter.merge" ) ? "true" : "false")  << std::endl;
    MPI_Barrier( this->worldComm().comm() );
    if ( boption ( _name = "exporter.merge" ) )
        writeMerge();
    else 
        write ();
    */

    /* make sure to reset values from previous calls */
    M_XDMFContent.str("");

    /* build file name */
    M_fileName.str("");
    M_fileName << this->prefix();
    std::cout << this->prefix() << std::endl;
    if( ! boption( _name = "exporter.merge" ) )
    {
        M_fileName << "-" << this->worldComm().globalSize() << "_" << this->worldComm().globalRank();
    }

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
    std::ostringstream cbuf;
    MPI_File fh;
    MPI_Status status;

    // sequential worldcomm
    WorldComm const & wcs = this->worldComm().subWorldCommSeq();

    /* build file name */
    cbuf << M_fileName.str() << ".xmf";

    /* Test master rank only once */
    /* whether we used a normal worldcomm or a sequentialized one */
    bool isMasterRank = false;
    if( boption( _name = "exporter.merge" ) )
    {
        isMasterRank = this->worldComm().isMasterRank();
    }
    else
    {
        isMasterRank = wcs.isMasterRank();
    }

    /* Open file with MPI IO */
    char * strTmp = strdup(cbuf.str().c_str());
    if(isMasterRank && fs::exists(strTmp))
    {
        MPI_File_delete(strTmp, MPI_INFO_NULL);
    }
    if( boption( _name = "exporter.merge" ) )
    {
        MPI_Barrier(this->worldComm().comm());
        MPI_File_open(this->worldComm().comm(), strTmp, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    }
    else
    {
        MPI_Barrier(wcs.comm());
        MPI_File_open(wcs.comm(), strTmp, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    }
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

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeHDF5() const 
{
    int size;
    std::ostringstream oss;

    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();
    
    timeset_ptrtype __ts = *__ts_it;

    while ( __ts_it != __ts_en )
    {
        __ts = *__ts_it;
        typename timeset_type::step_const_iterator __it = __ts->beginStep();
        typename timeset_type::step_const_iterator __end = __ts->endStep();
        __it = boost::prior ( __end );

        while ( __it != __end )
        {
            typename timeset_type::step_ptrtype __step = * __it;

            if ( __step->isInMemory() )
            {
                oss << M_fileName.str() << "-" << __step->index() << ".h5";

                /* open the corresponding h5 file for output */
                if( boption( _name = "exporter.merge" ) )
                {
                    M_HDF5.openFile(oss.str(), this->worldComm(), false);
                }
                else
                {
                    M_HDF5.openFile(oss.str(), this->worldComm().subWorldCommSeq(), false);
                }

                /* TODO find a way to incorporate time steps directly info the xml data */
                /*
                M_str.str("");
                M_str << "           <Grid Name=\"" << M_fileNameStep << "\" GridType=\"Uniform\">" << std::endl;
                M_str << "               <Time Value=\"" << __step->time() << "\"/>" << std::endl;  
                */

                M_XDMFContent << "<Grid Name=\"" << oss.str() << "\" GridType=\"Uniform\">" << std::endl;

                writePoints(__step);
                writeElements(__step);
                
                if ( this->worldComm().globalRank() == this->worldComm().masterRank() )
                {
                    //std::cout << "time                          : " << __step->time () << std::endl; 
                    //std::cout << "time increment                : " << __ts->timeIncrement () << std::endl; 
                    //std::cout << "numberOfSteps                 : " << __ts->numberOfSteps () << std::endl; 
                    //std::cout << "numberOfTotalSteps            : " << __ts->numberOfTotalSteps () << std::endl;
                    //std::cout << "file generated                : " << M_fileNameStep  << ".h5" << std::endl;
                }

                saveNodal(__step, __step->beginNodalScalar(), __step->endNodalScalar() );
                saveNodal(__step, __step->beginNodalVector(), __step->endNodalVector() );
                saveNodal( __step, __step->beginNodalTensor2(), __step->endNodalTensor2() );

                saveElement(__step, __step->beginElementScalar(), __step->endElementScalar() );
                saveElement(__step, __step->beginElementVector(), __step->endElementVector() );
                saveElement( __step, __step->beginElementTensor2(), __step->endElementTensor2() );

                M_XDMFContent << "</Grid>" << std::endl;

                // TODO Check if it must be done here
                M_HDF5.closeFile();
            }
            ++__it;
        }
        ++__ts_it;
    }
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::visit ( mesh_type* mesh) 
{
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writePoints(typename timeset_type::step_ptrtype __step) const 
{
    std::ostringstream oss;
    auto pt_it = __step->mesh()->beginPointWithProcessId();
    auto const pt_en = __step->mesh()->endPointWithProcessId();
    M_maxNumPoints= std::distance (pt_it, pt_en);

    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];

    currentSpaceDims[0] = 1;
    currentSpaceDims[1] = M_maxNumPoints;

    currentCount[0] = M_maxNumPoints;
    currentCount[1] = 3;

    /* if we are using MPI IO all the processes must create the tables */
    /* for all processes of we end up in a deadlock (see HDF5 FAQ) */
    if( boption( _name = "exporter.merge" ) )
    {
        oss.str("");
        oss << this->worldComm().globalRank();
        std::vector <size_type> globalNumPoint;
        globalNumPoint.resize (this->worldComm().globalSize(), 0);
        boost::mpi::all_gather(this->worldComm(), M_maxNumPoints, globalNumPoint);
        for (size_type i = 0; i < this->worldComm().globalSize (); i++)
        {
            std::ostringstream filestr;
            filestr << i;
            currentCount[0] = globalNumPoint[i];
            M_HDF5.createTable(filestr.str(), "point_coords", H5T_IEEE_F64BE, currentCount, false);
        }
        currentCount[0] = M_maxNumPoints;
    }
    else
    {
        M_HDF5.createTable ("point_coords", H5T_IEEE_F64BE, currentCount);
        M_HDF5.createTable ("point_ids", H5T_STD_U32BE, currentSpaceDims);
    }

    M_uintBuffer.resize (currentSpaceDims[0]*currentSpaceDims[1], 0);
    M_realBuffer.resize (currentCount[0]*currentCount[1], 0);

    for (size_type i = 0; i < M_maxNumPoints; i++ , pt_it++) 
    {
        M_uintBuffer[i] = pt_it->id ();

        M_realBuffer[3*i] = pt_it->node()[0];
        if (mesh_type::nRealDim >= 2)
            M_realBuffer[3*i + 1] = pt_it->node()[1];
        if (mesh_type::nRealDim >= 3)
            M_realBuffer[3*i + 2] = pt_it->node()[2];
    }

    bubbleSort (&M_uintBuffer[0], &M_realBuffer[0], M_maxNumPoints);

    for (size_type i = 0; i < M_maxNumPoints; i ++) 
        M_newPointId[M_uintBuffer[i]] = i;

    hsize_t currentOffset[2] = {0, 0};

    /* write the point coordinates */
    if( boption( _name = "exporter.merge" ) )
    {
        M_HDF5.write(oss.str() + "point_coords", H5T_NATIVE_DOUBLE, currentCount, currentOffset, &M_realBuffer[0]);

        /* close all the opened tables */
        for (size_type i = 0; i < this->worldComm().globalSize(); i++)
        {
            std::ostringstream filestr1;
            filestr1 << i;
            M_HDF5.closeTable(filestr1.str()+"point_coords");
        }
    }
    else
    {
        M_HDF5.write("point_coords", H5T_NATIVE_DOUBLE, currentCount, currentOffset, &M_realBuffer[0]);
        M_HDF5.write("point_ids", H5T_NATIVE_LLONG, currentSpaceDims, currentOffset , &M_uintBuffer[0]);

        /* close all the opened tables */
        M_HDF5.closeTable("point_coords");
        M_HDF5.closeTable("point_ids");    
    }

    M_XDMFContent << "<Geometry GeometryType=\"XYZ\">" << std::endl;
    M_XDMFContent << "<DataItem Dimensions=\"" << M_maxNumPoints << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">" << std::endl;
    if( boption( _name = "exporter.merge" ) )
    {
        M_XDMFContent << M_fileName.str() << "-" << __step->index() << ".h5:/" << this->worldComm().rank() << "/point_coords" << std::endl;
    }
    else
    {
        M_XDMFContent << M_fileName.str() << "-" << __step->index() << ".h5:/point_coords" << std::endl;
    }
    M_XDMFContent << "</DataItem>" << std::endl;
    M_XDMFContent << "</Geometry>" << std::endl;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeElements(typename timeset_type::step_ptrtype __step) const 
{
    std::ostringstream oss;

    typename mesh_type::parts_const_iterator_type p_it = __step->mesh()->beginParts();
    typename mesh_type::parts_const_iterator_type p_en = __step->mesh()->endParts();
    M_numParts = std::distance (p_it, p_en);
    M_maxNumElements = 0;
    for (int i = 0; i < M_numParts; i++, p_it ++) 
    {
        auto elt_it = __step->mesh()->beginElementWithMarker (p_it->first);
        auto elt_en = __step->mesh()->endElementWithMarker (p_it->first);
        M_maxNumElements += std::distance (elt_it, elt_en);
    }
    M_elementNodes = __step->mesh()-> numLocalVertices ();

    hsize_t currentSpacesDims[2];
    hsize_t currentSpacesDims2[2];

    currentSpacesDims [0] = M_maxNumElements;
    currentSpacesDims [1] = M_elementNodes;

    currentSpacesDims2 [0] = 1;
    currentSpacesDims2 [1] = M_maxNumElements;

    if( boption( _name = "exporter.merge" ) )
    {
        std::ostringstream filestr0;
        filestr0 << this->worldComm().globalRank();
        std::vector <size_type> globalNumElements;
        globalNumElements.resize (this->worldComm().globalSize(), 0);
        boost::mpi::all_gather(this->worldComm(), M_maxNumElements, globalNumElements);

        for (size_type i = 0; i < this->worldComm().globalSize (); i++)
        {
            std::ostringstream filestr;
            filestr << i;
            currentSpacesDims [0] = globalNumElements[i];

            M_HDF5.createTable (filestr.str(), "element_nodes", H5T_STD_U32BE, currentSpacesDims, true);
        }
        currentSpacesDims [0] = M_maxNumElements;
    }
    else
    {
        M_HDF5.createTable ("element_ids", H5T_STD_U32BE, currentSpacesDims2);
        M_HDF5.createTable ("element_nodes", H5T_STD_U32BE, currentSpacesDims);
    }

    M_uintBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], 0);
    std::vector<size_type> idsBuffer;
    idsBuffer.resize (currentSpacesDims2[1], 0);

    M_uintBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], 0);

    size_type k = 0;
    size_type i = 0;
    for (p_it = __step->mesh()->beginParts (); k < M_numParts;  p_it++ , k++)
    {
        auto elt_it = __step->mesh()->beginElementWithMarker(p_it->first);
        auto elt_en = __step->mesh()->endElementWithMarker(p_it->first);
        for (; elt_it != elt_en; ++elt_it , i ++)
        {
            idsBuffer[i] = elt_it->id ();
            for ( size_type j = 0; j < M_elementNodes; j ++ )
                M_uintBuffer[j + M_elementNodes*i] = M_newPointId[elt_it->point(j).id()] ; 
        }
    }


    hsize_t currentOffset[2] = {0, 0};

    if( boption( _name = "exporter.merge" ) )
    {
        std::ostringstream filestr0;
        filestr0 << this->worldComm().globalRank();
        M_HDF5.write( filestr0.str()+"element_nodes", H5T_NATIVE_LLONG, currentSpacesDims, currentOffset, &M_uintBuffer[0] );

        for (size_type i = 0; i < this->worldComm().globalSize(); i++)
        {
            std::ostringstream filestr1;
            filestr1 << i;
            M_HDF5.closeTable(filestr1.str()+"element_nodes");
        }
    }
    else
    {
        M_HDF5.write ( "element_ids", H5T_NATIVE_LLONG, currentSpacesDims2, currentOffset, &idsBuffer[0] );
        M_HDF5.write ( "element_nodes", H5T_NATIVE_LLONG, currentSpacesDims, currentOffset, &M_uintBuffer[0] );

        M_HDF5.closeTable ("element_ids");
        M_HDF5.closeTable ("element_nodes");
    }

    M_XDMFContent << "<Topology TopologyType=\"" << M_element_type << "\" NumberOfElements=\"" << M_maxNumElements << "\" NodesPerElement=\"" << M_elementNodes << "\">" << std::endl;
    M_XDMFContent << "<DataItem Dimensions=\"" <<M_maxNumElements << " " << M_elementNodes << "\" NumberType=\"Int\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">" << std::endl;
    if( boption( _name = "exporter.merge" ) )
    {
        M_XDMFContent << M_fileName.str() << "-" << __step->index() << ".h5:/" << this->worldComm().rank() << "/element_nodes" << std::endl;
    }
    else
    {
        M_XDMFContent << M_fileName << "-" << __step->index() << ".h5:/element_nodes" << std::endl;
    }
    M_XDMFContent << "</DataItem>" << std::endl;
    M_XDMFContent << "</Topology>" << std::endl;

}

template <typename MeshType, int N>
template <typename Iterator>
void Exporterhdf5<MeshType, N>::saveNodal ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const 
{   
    while ( __var != en )
    {    
        std::string attributeType ("Scalar");
        std::string solutionName = __var->first;
        uint16_type nComponents = __var -> second.nComponents;

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

//        if ( this->worldComm().globalRank() == this->worldComm().masterRank() )
//            std::cout << "solution name                 : " << solutionName << std::endl;

        hsize_t currentSpacesDims [2];

        currentSpacesDims [0] = nComponents;
        currentSpacesDims [1] = M_maxNumPoints;

        if( boption( _name = "exporter.merge") )
        {
            std::ostringstream filestr0;
            filestr0 << this->worldComm().globalRank();

            std::vector <size_type> globalNumPoint;
            globalNumPoint.resize (this->worldComm().globalSize(), 0);
            // TODO remove this communication, already done in writePoints
            boost::mpi::all_gather(this->worldComm(), M_maxNumPoints, globalNumPoint);

            for (size_type i = 0; i < this->worldComm().globalSize (); i++)
            {
                std::ostringstream filestr;
                filestr << i;
                currentSpacesDims [1] = globalNumPoint[i];

                M_HDF5.createTable (filestr.str(), solutionName.c_str() , H5T_IEEE_F64BE, currentSpacesDims, true);
            }

            currentSpacesDims [1] = M_maxNumPoints;
        }
        else
        {
            M_HDF5.createTable (solutionName.c_str(), H5T_IEEE_F64BE, currentSpacesDims);
        }

        M_realBuffer.resize (M_maxNumPoints*nComponents, 0);

        typename mesh_type::parts_const_iterator_type p_it = __step->mesh()->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = __step->mesh()->endParts();
        M_numParts = std::distance (p_it, p_en);

        for (; p_it != p_en; p_it++) 
        {

            auto r = markedelements(__step->mesh(), p_it->first, EntityProcessType::LOCAL_ONLY);
            auto elt_it = r.template get<1>();
            auto elt_en = r.template get<2>();

            Feel::detail::MeshPoints<float> mp ( __step->mesh().get(), this->worldComm(), elt_it, elt_en, true, true, true );

            size_type e = 0; 

            for (; elt_it != elt_en; ++elt_it )
            {
                for ( uint16_type c = 0; c < nComponents; ++c )
                {
                    for ( uint16_type p = 0; p < __step->mesh()->numLocalVertices(); ++p, ++e )
                    {
                        size_type  ptid = M_newPointId[elt_it->get().point(p).id()];
                        size_type global_node_id = mp.ids.size()*c + ptid;
                        if ( c < __var->second.nComponents ) 
                        {
                            size_type dof_id = boost::get<0>( __var->second.functionSpace()->dof()->localToGlobal ( elt_it->get().id(), p, c ) );
                            M_realBuffer[global_node_id] = __var->second.globalValue ( dof_id );
                        }
                        else
                        {
                            M_realBuffer[global_node_id] = 0.0;
                        }
                    }
                }
            }
        }    
        hsize_t currentOffset[2] = {0, 0};

        if( boption( _name = "exporter.merge" ) )
        {
            std::ostringstream filestr0;
            filestr0 << this->worldComm().globalRank ();

            M_HDF5.write ( filestr0.str() + solutionName, H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, &M_realBuffer[0] );

            for (size_type i = 0; i < this->worldComm().globalSize(); i++)
            {
                std::ostringstream filestr1;
                filestr1 << i;
                M_HDF5.closeTable(filestr1.str()+solutionName);
            }
        }
        else
        {
            M_HDF5.write ( solutionName.c_str(), H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, &M_realBuffer[0] );

            M_HDF5.closeTable (solutionName.c_str());    
        }

        M_XDMFContent << "<Attribute AttributeType=\""<< attributeType << "\" Name=\"" << solutionName << "\" Center=\"Node\">" << std::endl;
        M_XDMFContent << "<DataItem Dimensions=\""<< nComponents << " " << M_maxNumPoints << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">" << std::endl;
        if( boption( _name = "exporter.merge" ) )
        {
            M_XDMFContent << M_fileName.str() << "-" << __step->index() << ".h5:/" << this->worldComm().rank() << "/" << solutionName << std::endl;
        }
        else
        {
            M_XDMFContent << M_fileName << "-" << __step->index() << ".h5:/" << solutionName << std::endl;
        }
        M_XDMFContent << "</DataItem>" << std::endl;    
        M_XDMFContent << "</Attribute>" << std::endl;
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

        currentSpacesDims [0] = nComponents;
        currentSpacesDims [1] = M_maxNumElements;

        if( boption( _name = "exporter.merge" ) )
        {
            std::ostringstream filestr0;
            filestr0 << this->worldComm().globalRank();
            std::vector <size_type> globalNumElements;
            globalNumElements.resize (this->worldComm().globalSize(), 0);
            boost::mpi::all_gather(this->worldComm(), M_maxNumElements, globalNumElements);

            for (size_type i = 0; i < this->worldComm().globalSize (); i++)
            {
                std::ostringstream filestr;
                filestr << i;
                currentSpacesDims [1] = globalNumElements[i];

                M_HDF5.createTable (filestr.str(), solutionName, H5T_IEEE_F64BE, currentSpacesDims, true);
            }
            currentSpacesDims [1] = M_maxNumElements;
        }
        else
        {
            M_HDF5.createTable (solutionName.c_str(), H5T_IEEE_F64BE, currentSpacesDims);
        }

        M_realBuffer.resize (M_maxNumElements*nComponents, 0);

        typename mesh_type::parts_const_iterator_type p_it = __step->mesh()->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = __step->mesh()->endParts();
        for (; p_it != p_en; p_it++ ) 
        {
            typename mesh_type::marker_element_const_iterator elt_st;
            typename mesh_type::marker_element_const_iterator elt_en;
            boost::tie( elt_st, elt_en ) = __step->mesh()->elementsWithMarker( p_it->first, __evar->second.worldComm().localRank() );

            if ( !__evar->second.areGlobalValuesUpdated() )
                __evar->second.updateGlobalValues();

            size_type ncells = std::distance ( elt_st, elt_en );
            for ( int c = 0; c < nComponents; ++c )
            {
                size_type e = 0;
                for ( auto elt_it = elt_st; elt_it != elt_en; ++elt_it, ++e )
                {
                    size_type global_node_id = c*ncells+e;
                    if ( c < __evar->second.nComponents )
                    {
                        size_type dof_id = boost::get<0>( __evar->second.functionSpace()->dof()->localToGlobal( elt_it->id(), 0, c ) );
                        M_realBuffer[global_node_id] = __evar->second.globalValue ( dof_id );
                    }
                    else 
                        M_realBuffer[global_node_id] = 0;
                }
            }
        }   
        hsize_t currentOffset[2] = {0, 0};

        if( boption( _name = "exporter.merge" ) )
        {
            std::ostringstream filestr0;
            filestr0 << this->worldComm().globalRank();
            
            M_HDF5.write ( filestr0.str() + solutionName, H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, &M_realBuffer[0] );

            for (size_type i = 0; i < this->worldComm().globalSize(); i++)
            {
                std::ostringstream filestr1;
                filestr1 << i;
                M_HDF5.closeTable(filestr1.str()+solutionName);
            }
        }
        else
        {
            M_HDF5.write ( solutionName.c_str(), H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, &M_realBuffer[0] );

            M_HDF5.closeTable (solutionName.c_str());    
        }

        M_XDMFContent << "<Attribute AttributeType=\"" << attributeType << "\" Name=\"" << solutionName << "\" Center=\"Cell\">" << std::endl;
        M_XDMFContent << "<DataItem Dimensions=\""<< nComponents <<" "<< M_maxNumElements << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">" << std::endl;
        if( boption( _name = "exporter.merge" ) )
        {
            M_XDMFContent << M_fileName.str() << "-" << __step->index() << ".h5:/" << this->worldComm().rank() << "/" << solutionName << std::endl;
        }
        else
        {
            M_XDMFContent << M_fileName.str() << "-" << __step->index() << ".h5:/" << solutionName << std::endl;
        }
        M_XDMFContent << "</DataItem>" << std::endl;    
        M_XDMFContent << "</Attribute>" << std::endl;
        __evar++;
    }
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeStats () const 
{
    size_type numStats = 3 ;
    hsize_t currentSpacesDims[2] = {1, numStats} ;

    M_HDF5.createTable ( "stats", H5T_STD_U32BE, currentSpacesDims ) ;
    M_uintBuffer.resize (numStats) ;

    M_uintBuffer[0] = M_maxNumPoints ;
    M_uintBuffer[1] = M_maxNumElements ;
    M_uintBuffer[2] = M_elementNodes ;

    if ( this->worldComm().globalRank() == this->worldComm().masterRank() )
    {
        std::cout << "nombre de Points              : " << M_maxNumPoints << std::endl ;
        std::cout << "M_numMaxElements              : " << M_maxNumElements << std::endl ;
        std::cout << "nombre de Points par element  : " << M_elementNodes << std::endl ;
        std::cout << "M_numParts                    : " << M_numParts << std::endl ;
        std::cout << "mesh_type::nRealDim           : " << mesh_type::nRealDim << std::endl ;
        //std::cout << "fileNameStep                  : " << M_fileNameStep << ".h5" << std::endl ;
        std::cout << "fileName                      : " << M_fileName << ".xmf" << std::endl ;
    }

    hsize_t currentOffset [2] = {0, 0} ;
    M_HDF5.write ("stats", H5T_NATIVE_LLONG, currentSpacesDims, currentOffset, &M_uintBuffer[0]) ;
    M_HDF5.closeTable ("stats") ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::bubbleSort (size_type * ids, value_type * coords, size_type n) const 
{
    size_type int_tmp;
    value_type float_tmp;
    bool swapped = false;
    do 
    {
        swapped = false;
        for (size_type j = 0; j < n-1; j ++) 
        {
            if (ids[j] > ids[j+1]) 
            {
                int_tmp = ids[j];
                ids[j] = ids[j+1];
                ids[j+1] = int_tmp;

                float_tmp = coords[3*j];
                coords[3*j] = coords [3*(j+1)];
                coords[3*(j+1)] = float_tmp;

                float_tmp = coords[3*j+1];
                coords[3*j+1] = coords [3*(j+1)+1];
                coords[3*(j+1)+1] = float_tmp;

                float_tmp = coords[3*j+2];
                coords[3*j+2] = coords [3*(j+1)+2];
                coords[3*(j+1)+2] = float_tmp;

                swapped = true;
            }
        }
        n = n -1;
    }
    while (swapped);
}

}
#endif /* FEELL_HAS_HDF5 */
#endif /* __Exporterhdf5_CPP */

