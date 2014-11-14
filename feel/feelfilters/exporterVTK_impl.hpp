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
   \file exporterVTK_impl.cpp
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \date 2014-11-13
*/
#ifndef __EXPORTERVTK_CPP
#define __EXPORTERVTK_CPP 1

#include <feel/feelfilters/exporterVTK.hpp>

namespace Feel
{
template<typename MeshType, int N>
ExporterVTK<MeshType,N>::ExporterVTK( WorldComm const& worldComm )
:
super( worldComm ),
M_element_type()
{
    init();
}
template<typename MeshType, int N>
ExporterVTK<MeshType,N>::ExporterVTK( std::string const& __p, int freq,
                                        WorldComm const& worldComm )
    :
    super( "vtk", __p, freq, worldComm )
{
    init();
}
template<typename MeshType, int N>
ExporterVTK<MeshType,N>::ExporterVTK( po::variables_map const& vm, std::string const& exp_prefix,
                                        WorldComm const& worldComm )
    :
    super( vm, exp_prefix, worldComm )
{
    init();
}
template<typename MeshType, int N>
ExporterVTK<MeshType,N>::ExporterVTK( ExporterVTK const & __ex )
    :
    super( __ex )
{}
template<typename MeshType, int N>
ExporterVTK<MeshType,N>::~ExporterVTK()
{}

template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::init()
{
    if ( mesh_type::nDim == 1 )
        if ( mesh_type::Shape == SHAPE_LINE )
        {
            M_element_type = ( mesh_type::nOrder == 1 )? VTK_LINE : VTK_POLY_LINE;
            M_face_type = VTK_VERTEX;
        }

    if ( mesh_type::nDim == 2 )
    {
        if ( mesh_type::Shape == SHAPE_TRIANGLE )
            M_element_type = ( mesh_type::nOrder == 1 )? VTK_TRIANGLE : VTK_QUADRATIC_TRIANGLE;

        else if ( mesh_type::Shape == SHAPE_QUAD )
            M_element_type = ( mesh_type::nOrder == 1 )? VTK_QUAD : VTK_QUADRATIC_QUAD;

        M_face_type = ( mesh_type::nOrder == 1 )? VTK_LINE : VTK_POLY_LINE;
    }

    if ( mesh_type::nDim == 3 )
    {
        if ( mesh_type::Shape == SHAPE_TETRA )
        {
            M_element_type = ( mesh_type::nOrder == 1 )? VTK_TETRA: VTK_QUADRATIC_TETRA;
            M_face_type = ( mesh_type::nOrder == 1 )? VTK_TRIANGLE : VTK_QUADRATIC_TRIANGLE;
        }

        else if ( mesh_type::Shape == SHAPE_HEXA )
        {
            M_element_type = ( mesh_type::nOrder == 1 )? VTK_HEXAHEDRON : VTK_QUADRATIC_HEXAHEDRON;
            M_face_type = ( mesh_type::nOrder == 1 )? VTK_QUAD : VTK_QUADRATIC_QUAD;
        }
    }

    //std::cout << "Faces: " << M_face_type << "; Elements: " << M_element_type << std::endl;
}

template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::save() const
{
    DVLOG(2) << "[ExporterVTK] checking if frequency is ok\n";

    if ( this->cptOfSave() % this->freq()  )
    {
        this->saveTimeSet();
        return;
    }

    DVLOG(2) << "[ExporterVTK] frequency is ok\n";

    DVLOG(2) << "[ExporterVTK] save()...\n";

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
            typename timeset_type::step_ptrtype __step = *__it;

            if ( __step->isInMemory() )
            {
                /* instanciante data object */
                vtkSmartPointer<vtkout_type> out = vtkSmartPointer<vtkout_type>::New();

                /* write data into vtk object */
                saveMesh(__step, out);
                saveNodeData( __step, __step->beginNodalScalar(), __step->endNodalScalar(), out );
                saveNodeData( __step, __step->beginNodalVector(), __step->endNodalVector(), out );
                saveNodeData( __step, __step->beginNodalTensor2(), __step->endNodalTensor2(), out );
                saveElementData( __step, __step->beginElementScalar(), __step->endElementScalar(), out );
                saveElementData( __step, __step->beginElementVector(), __step->endElementVector(), out );
                saveElementData( __step, __step->beginElementTensor2(), __step->endElementTensor2(), out );

                /* Build file name */
                std::ostringstream fname;
                fname   << __ts->name()  //<< this->prefix() //this->path()
                    << "-" << __step->index()
                    << "-" << this->worldComm().size() << "_" << this->worldComm().rank()
                    << ".vtk";

                /*
                out->GetOutputInformation(0).Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_NUMBER_OF_PIECES(), this->worldComm().globalSize());
                out->GetOutputInformation(0).Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_PIECE_NUMBER(), this->worldComm().rank());
                */

#if 0
                vtkSmartPointer<vtkMultiBlockDataSet> mbds = vtkSmartPointer<vtkMultiBlockDataSet>::New();
                mbds->SetNumberOfBlocks( this->worldComm().globalSize() );

                std::ostringstream oss;

                for( unsigned int block = 0 ; block < mbds->GetNumberOfBlocks(); ++block )
                {
                    if( block == this->worldComm().rank() )
                    {
                        oss.str("");
                        oss << "P" << this->worldComm().rank();
                        mbds->SetBlock( block, out );
                        mbds->GetMetaData(block)->Set(vtkCompositeDataSet::NAME(), oss.str().c_str() );
                        mbds->GetMetaData(block)->Set(vtkCompositeDataSet::DATA_TIME_STEP(), __step->index());
                    } // END if we own the block
                    else
                    {
                        mbds->SetBlock( block, NULL );
                    } // END else we don't own the block
                } // END for all blocks

                /* Build vtk objects while reusing the current mpi communicator */
                /* before version 5.10, we cannot initialize a MPIController with an external MPI_Comm */
                #if VTK_MAJOR_VERSION <= 5 && VTK_MINOR_VERSION < 10
                    //#warning parallel export not supported with this version of VTK
                #else  
#if 1
                MPI_Comm comm = this->worldComm().comm();
                vtkMPICommunicatorOpaqueComm * mpicoc = new vtkMPICommunicatorOpaqueComm(&comm);
                vtkSmartPointer<vtkMPICommunicator> mpicomm = vtkSmartPointer<vtkMPICommunicator>::New();
                    mpicomm->InitializeExternal(mpicoc);
                vtkSmartPointer<vtkMPIController> mpictrl = vtkSmartPointer<vtkMPIController>::New();
                    mpictrl->SetCommunicator(mpicomm);
#endif
                fname.str("");
                fname  << __ts->name()  //<< this->prefix() //this->path()
                    << "-" << __step->index()
                    << ".vtm";

                vtkSmartPointer<vtkXMLPMultiBlockDataWriter> xmlpw = vtkSmartPointer<vtkXMLPMultiBlockDataWriter>::New();
                    xmlpw->SetController(mpictrl);
                    xmlpw->SetFileName(fname.str().c_str());
#if VTK_MAJOR_VERSION <= 5
                    xmlpw->SetInput(mbds);
#else
                    xmlpw->SetInputData(mbds);
#endif
                    xmlpw->Update();

                //delete mpicomm;
                #endif

#else
#if 0
                vtkSmartPointer<vtkoutwriter_type> dw = vtkSmartPointer<vtkoutwriter_type>::New();
#if VTK_MAJOR_VERSION <= 5
                dw->SetInput(out);
#else
                dw->SetInputData(out);
#endif
                dw->SetFileName(fname.str().c_str());
                dw->Update();
#else
                std::ostringstream mfname;
                mfname   << __ts->name()  //<< this->prefix() //this->path()
                    << "-" << __step->index()
                    << "-" << this->worldComm().size() << "_" << this->worldComm().rank()
                    << ".vtm";
                vtkSmartPointer<vtkMultiBlockDataSet> mbds = vtkSmartPointer<vtkMultiBlockDataSet>::New();
                    mbds->SetNumberOfBlocks( 1 );

                    int blockID = 0;
                    mbds->SetBlock(blockID, out);
                    mbds->GetMetaData(blockID)->Set(vtkCompositeDataSet::NAME(), "marker 1");
                    //mbds->GetMetaData(blockID)->Set(vtkCompositeDataSet::DATA_TIME_STEP(), __step->index());

                vtkXMLMultiBlockDataWriter* mbw = vtkXMLMultiBlockDataWriter::New();
                    mbw->SetFileName(mfname.str().c_str());
#if VTK_MAJOR_VERSION <= 5
                    mbw->SetInput(mbds);
#else
                    mbw->SetInputData(mbds);
#endif
                    mbw->Write();
#endif
#endif
            }
            __it++;
        }
        __ts_it++;
    }

    DVLOG(2) << "[ExporterVTK] saving done\n";

    this->saveTimeSet();
}

template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::saveMesh( typename timeset_type::step_ptrtype step, vtkSmartPointer<vtkout_type> out ) const
{
    /* get local elements */
    //auto r = markedelements(step->mesh(), markerid, EntityProcessType::LOCAL_ONLY );
    auto r = elements(step->mesh());
    auto elt_it = r.template get<1>();
    auto elt_en = r.template get<2>();

    /* Gather points and elements into vectors + info about data */
    Feel::detail::MeshPoints<float> mp( step->mesh().get(), this->worldComm(), elt_it, elt_en, false, true, true, 0 );

    /* Add points to data structure */
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetDataTypeToFloat();
        points->SetNumberOfPoints(mp.ids.size());

    for ( int i = 0; i < mp.ids.size() ; i++ )
    {
        points->SetPoint( (vtkIdType)(mp.ids[i]), (float *)(mp.coords.data()) + i * mesh_type::element_type::numPoints );
    } 

    out->SetPoints(points);

    /* Add cells to data structure */
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkelement_type> cell;

    for( int i = 0; i < mp.elem.size(); i+=mesh_type::element_type::numPoints )
    {
        cell = vtkSmartPointer<vtkelement_type>::New();
        for( int p=0; p < mesh_type::element_type::numPoints; ++p )
        {
            cell->GetPointIds()->SetId(p, mp.elem[i + p]);
        }
        cells->InsertNextCell(cell);
    }

    out->SetCells(M_element_type, cells);
}

template<typename MeshType, int N>
template<typename Iterator>
void
ExporterVTK<MeshType,N>::saveNodeData( typename timeset_type::step_ptrtype step, Iterator __var, Iterator en, vtkSmartPointer<vtkout_type> out ) const
{
    while ( __var != en )
    {
        if ( !__var->second.worldComm().isActive() ) return;

        /* handle faces data */
#if 0
        if ( option( _name="exporter.ensightgold.save-face" ).template as<bool>() )
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
#endif
        /* handle elements */
        uint16_type nComponents = __var->second.nComponents;

        VLOG(1) << "nComponents field: " << nComponents;
        if ( __var->second.is_vectorial )
        {
            nComponents = 3;
            VLOG(1) << "nComponents field(is_vectorial): " << nComponents;
        }

        /* we get that from the local processor */
        /* We do not need the renumbered global index */
        //auto r = markedelements(__mesh, M_markersToWrite[i], EntityProcessType::ALL);
        auto r = elements( step->mesh() );
        auto elt_it = r.template get<1>();
        auto elt_en = r.template get<2>();

        Feel::detail::MeshPoints<float> mp( step->mesh().get(), this->worldComm(), elt_it, elt_en, false, true, true, 0 );
        
        // previous implementation
        //size_type __field_size = mp.ids.size();
        //int nelts = std::distance(elt_it, elt_en);
        int npts = mp.ids.size();
        size_type __field_size = npts;
        if ( __var->second.is_vectorial )
            __field_size *= 3;
        std::vector<float> __field( __field_size, 0.0 );
        size_type e = 0;
        VLOG(1) << "field size=" << __field_size;
        if ( !__var->second.areGlobalValuesUpdated() )
            __var->second.updateGlobalValues();

        vtkSmartPointer<vtkFloatArray> da = vtkSmartPointer<vtkFloatArray>::New();
        da->SetName(__var->first.c_str());

        /* set array parameters */
        /* no need for preallocation if we are using Insert* methods */
        da->SetNumberOfComponents(nComponents);
        da->SetNumberOfTuples(npts);

        /*
           std::cout << this->worldComm().rank() << " nbPts:" << npts << " nComp:" << nComponents 
           << " __var->second.nComponents:" << __var->second.nComponents << std::endl;
        */

        /*
           std::cout << this->worldComm().rank() << " marker=" << *mit << " nbPts:" << npts << " nComp:" << nComponents 
           << " __evar->second.nComponents:" << __var->second.nComponents << std::endl;
           */

        /* loop on the elements */
        int index = 0;
        for ( ; elt_it != elt_en; ++elt_it )
        {
            VLOG(3) << "is ghost cell " << elt_it->isGhostCell();
            /* looop on the ccomponents is outside of the loop on the vertices */
            for ( uint16_type c = 0; c < nComponents; ++c )
            {
                for ( uint16_type p = 0; p < step->mesh()->numLocalVertices(); ++p, ++e )
                {
                    size_type ptid = mp.old2new[elt_it->point( p ).id()];
                    size_type global_node_id = ptid * nComponents + c;
                    //size_type global_node_id = mp.ids.size()*c + ptid ;
                    //LOG(INFO) << elt_it->get().point( p ).id() << " " << ptid << " " << global_node_id << std::endl;
                    DCHECK( ptid < step->mesh()->numPoints() ) << "Invalid point id " << ptid << " element: " << elt_it->id()
                        << " local pt:" << p
                        << " mesh numPoints: " << step->mesh()->numPoints();
                    DCHECK( global_node_id < __field_size ) << "Invalid dof id : " << global_node_id << " max size : " << __field_size;

                    if ( c < __var->second.nComponents )
                    {
                        size_type dof_id = boost::get<0>( __var->second.functionSpace()->dof()->localToGlobal( elt_it->id(), p, c ) );

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

        /* insert data into array */
        for(int i = 0; i < npts; i++)
        {
            da->SetTuple(mp.ids[i], __field.data() + i * nComponents);
        }

        /* add data array into the vtk object */
        out->GetPointData()->AddArray(da);

        DVLOG(2) << "[ExporterVTK::saveNodal] saving " << __var->first << "done\n";

        ++__var;
    }
}

template<typename MeshType, int N>
template<typename Iterator>
void
ExporterVTK<MeshType,N>::saveElementData( typename timeset_type::step_ptrtype step, Iterator __evar, Iterator __evaren, vtkSmartPointer<vtkout_type> out ) const
{
    while ( __evar != __evaren )
    {
        if ( !__evar->second.worldComm().isActive() ) return;

        auto mesh = step->mesh();

        auto r = elements( step->mesh() );
        auto elt_st = r.template get<1>();
        auto elt_en = r.template get<2>();

        if ( !__evar->second.areGlobalValuesUpdated() )
            __evar->second.updateGlobalValues();

        //size_type ncells = __evar->second.size()/__evar->second.nComponents;
        size_type ncells = std::distance( elt_st, elt_en );

        uint16_type nComponents = __evar->second.nComponents;

        if ( __evar->second.is_vectorial )
            nComponents = 3;

        size_type __field_size = nComponents * ncells;
        std::vector<float> __field( __field_size );
        __field.clear();

        DVLOG(2) << "[saveElement] firstLocalIndex = " << __evar->second.firstLocalIndex() << "\n";
        DVLOG(2) << "[saveElement] lastLocalIndex = " << __evar->second.lastLocalIndex() << "\n";
        DVLOG(2) << "[saveElement] field.size = " << __field_size << "\n";

        vtkSmartPointer<vtkFloatArray> da = vtkSmartPointer<vtkFloatArray>::New();
        da->SetName(__evar->first.c_str());

        /* set array parameters */
        da->SetNumberOfComponents(nComponents);
        da->SetNumberOfTuples(ncells);

        /*
           std::cout << this->worldComm().rank() << " nbElts:" << ncells << " nComp:" << nComponents 
           << " __evar->second.nComponents:" << __evar->second.nComponents << std::endl;
        */
        /*
           std::cout << this->worldComm().rank() << " marker=" << *mit << " nbElts:" << ncells << " nComp:" << nComponents 
           << " __evar->second.nComponents:" << __evar->second.nComponents << std::endl;
           */

        float * array = new float[nComponents];

        size_type e = 0;
        for ( auto elt_it = elt_st ; elt_it != elt_en; ++elt_it, ++e )
        {
            auto const& elt = boost::unwrap_ref( *elt_it );
            DVLOG(2) << "pid : " << this->worldComm().globalRank()
                << " elt_it :  " << elt.id()
                << " e : " << e << "\n";

            for ( int c = 0; c < nComponents; ++c )
            {
                size_type global_node_id = e*ncells+c ;

                if ( c < __evar->second.nComponents )
                {
                    size_type dof_id = boost::get<0>( __evar->second.functionSpace()->dof()->localToGlobal( elt.id(),0, c ) );

                    DVLOG(2) << "c : " << c
                        << " gdofid: " << global_node_id
                        << " dofid : " << dof_id
                        << " f.size : " <<  __field.size()
                        << " e.size : " <<  __evar->second.size()
                        << "\n";

                    //__field[global_node_id] = __evar->second.globalValue( dof_id );
                    array[c] =  __evar->second.globalValue( dof_id );

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
                {
                    //__field[global_node_id] = 0;
                    array[c] = 0.0;
                }
            }
            da->SetTuple(e, array);
        }

        delete[] array;

        /* add data array into the vtk object */
        out->GetCellData()->AddArray(da);

        DVLOG(2) << "[ExporterVTK::saveElement] saving " << __evar->first << "done\n";
        ++__evar;
    }
}


template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::visit( mesh_type* )
{
}

#if 0
#if defined( FEELPP_INSTANTIATION_MODE )

//
// explicit instances
//

# define DIMS BOOST_PP_TUPLE_TO_LIST(3,(1,2,3))
# define ORDERS BOOST_PP_TUPLE_TO_LIST(5,(1,2,3,4,5))
# define ORDERS_FUN_GMSH BOOST_PP_TUPLE_TO_LIST(5,(1,2,3,4,5))

// exporter VTK
# define FACTORY(LDIM,LORDER,ORDERFUN) template class ExporterVTK<Mesh<Simplex<LDIM,LORDER,LDIM> >, ORDERFUN >;
# define FACTORY_OP(_, GDO) FACTORY GDO

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_OP, 3, ( DIMS, ORDERS, ORDERS_FUN_GMSH ) )

#endif // FEELPP_INSTANTIATION_MODE
#endif
}
#endif // __EXPORTERVTK_CPP
