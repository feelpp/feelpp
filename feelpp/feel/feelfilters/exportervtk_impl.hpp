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
   \file exporterVTK_impl.cpp
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \date 2014-11-13
*/
#ifndef __EXPORTERVTK_CPP
#define __EXPORTERVTK_CPP 1

namespace Feel
{
template<typename MeshType, int N>
ExporterVTK<MeshType,N>::ExporterVTK( worldcomm_ptr_t const& worldComm )
:
super( worldComm ),
M_element_type()
{
    init();
}
template<typename MeshType, int N>
ExporterVTK<MeshType,N>::ExporterVTK( std::string const& __p, int freq,
                                        worldcomm_ptr_t const& worldComm )
    :
    super( "vtk", __p, freq, worldComm )
{
    init();
}
template<typename MeshType, int N>
ExporterVTK<MeshType,N>::ExporterVTK( po::variables_map const& vm, std::string const& exp_prefix,
                                        worldcomm_ptr_t const& worldComm )
    :
    super( vm, exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
ExporterVTK<MeshType,N>::ExporterVTK( std::string const& exp_prefix,
                                      worldcomm_ptr_t const& worldComm )
    :
    super( exp_prefix, worldComm )
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
{
#if VTK_MAJOR_VERSION >= 6 && defined(VTK_HAS_PARALLEL)
#if defined(FEELPP_VTK_INSITU_ENABLED)
    /* end up in situ simulation */
    if(boption( _name="exporter.vtk.insitu.enable" ))
    {
        inSituProcessor->Finalize();
        //inSituProcessor->Delete();
    }
#endif

    /* clean memory */
    delete this->opaqueComm;
#endif
}

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

    M_inSituEnable = false;
    M_inSituSave = false;
#if VTK_MAJOR_VERSION >= 6 && defined(VTK_HAS_PARALLEL)
    /* before version 5.10, we cannot initialize a MPIController with an external MPI_Comm */
    this->lComm = this->worldComm().comm();
    /* initialize the VTK communicator from the current MPI communicator */
    this->opaqueComm = new vtkMPICommunicatorOpaqueComm(&(this->lComm));

#if defined(FEELPP_VTK_INSITU_ENABLED)
    /* initialize in-situ visualization if needed */
    M_inSituEnable = boption( _name="exporter.vtk.insitu.enable" );
    M_inSituSave = boption( _name="exporter.vtk.insitu.save" );

    if(boption( _name="exporter.vtk.insitu.enable" ))
    {
        if(inSituProcessor == NULL)
        {
            inSituProcessor = vtkSmartPointer<vtkCPProcessor>::New();
            inSituProcessor->Initialize(*(this->opaqueComm));
            //inSituProcessor->DebugOn();
        }
        else
        {
            inSituProcessor->RemoveAllPipelines();
        }

        /* specify a user script */
        std::string pyscript = soption( _name="exporter.vtk.insitu.pyscript" );
        if(pyscript != "")
        {
            vtkSmartPointer<vtkCPPythonScriptPipeline> pipeline = vtkSmartPointer<vtkCPPythonScriptPipeline>::New();

            pipeline->Initialize(pyscript.c_str());
            inSituProcessor->AddPipeline(pipeline.GetPointer());
        }
        /* else revert to a basic VTK pipeline */
        else
        {
            vtkSmartPointer<vtkBaseInsituPipeline> pipeline = vtkSmartPointer<vtkBaseInsituPipeline>::New();
            pipeline->Initialize();
            inSituProcessor->AddPipeline(pipeline.GetPointer());
        }
    }
#endif
#endif

    //std::cout << "Faces: " << M_face_type << "; Elements: " << M_element_type << std::endl;
}

#ifdef FEELPP_HAS_LIBXML2
template<typename MeshType, int N>
int ExporterVTK<MeshType,N>::writeTimePVD(std::string xmlFilename, double timestep, std::string dataFilename, int partNo) const
{

    /*
       <!-- Sample xml code for PVD format -->
       <?xml version="1.0"?>
       <VTKFile type="Collection" version="0.1">
           <Collection>
               <DataSet timestep="0" group="" part="0" file="ts_0.vtm"/>
               <DataSet timestep="0.5" group="" part="0" file="ts_1.vtm"/>
           </Collection>
       </VTKFile>
   */

    int retcode = 0;

    xmlDocPtr doc = NULL;
    xmlNodePtr root = NULL, node1 = NULL, node2 = NULL;
    std::ostringstream oss;

    /* only do this on the master rank */
    if(this->worldComm().isMasterRank())
    {
        /* First Step: Find the Collection node */
        /* Either in a newly created file or in an existing file */
        /* check if the time file already exists */
        /* if so we update its data by adding a new DataSet node */
        if(boost::filesystem::exists(xmlFilename))
        {
            doc = xmlReadFile(xmlFilename.c_str(), NULL, 0);
            if (doc == NULL) {
                //std::cerr << "Failed to parse %s" << std::endl;
                return 1;
            }
            root = xmlDocGetRootElement(doc);

            /* check that we have VTKFile as a first entry */
            if(xmlStrncmp(root->name, BAD_CAST "VTKFile", 7) == 0 && root->children != NULL)
            {
                /* get first child */
                node1 = root->children;
                if(xmlStrncmp(node1->name, BAD_CAST "Collection", 10) != 0)
                { node1 = NULL; retcode = 1; }
            } 
            /* mark this as an error */
            else
            { retcode = 1; }
        }
        /* if the file does not already exists we create it */
        else
        {
            /* create a new document */
            doc = xmlNewDoc(BAD_CAST "1.0");
            root = xmlNewNode(NULL, BAD_CAST "VTKFile");
            xmlSetProp(root, BAD_CAST "type", BAD_CAST "Collection");
            xmlSetProp(root, BAD_CAST "version", BAD_CAST "0.1");
            xmlDocSetRootElement(doc, root);

            node1 = xmlNewNode(NULL, BAD_CAST "Collection");
            xmlAddChild(root, node1);
        }

        /* Second step */
        /* Create a new dataset entry to add to the Collection node */
        if(node1)
        {
            node2 = xmlNewNode(NULL, BAD_CAST "DataSet");
            xmlAddChild(node1, node2);

            oss.str("");
            oss << timestep;
            xmlSetProp(node2, BAD_CAST "timestep", BAD_CAST oss.str().c_str());
            xmlSetProp(node2, BAD_CAST "group", BAD_CAST "");
            oss.str("");
            oss << partNo;
            xmlSetProp(node2, BAD_CAST "part", BAD_CAST oss.str().c_str());
            xmlSetProp(node2, BAD_CAST "file", BAD_CAST dataFilename.c_str());

            /*
            xmlChar * mem = NULL;
            int size = 0;
            xmlDocDumpFormatMemory(doc, &mem, &size, 1);
            std::cout << mem << std::endl;
            xmlFree(mem);
            */

            //std::cout << "Writing file " << xmlFilename << std::endl;
            FILE * f = fopen(xmlFilename.c_str(), "w");
            xmlDocDump(f, doc);
            fclose(f);
        }

        xmlFreeDoc(doc);
    }

    return retcode;
}
#endif

template<typename MeshType, int N>
vtkSmartPointer<vtkMultiBlockDataSet>
ExporterVTK<MeshType,N>::buildMultiBlockDataSet( double time, std::map<std::string,vtkSmartPointer<vtkout_type>> const& outs ) const
{
    vtkSmartPointer<vtkMultiBlockDataSet> mbds = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    mbds->SetNumberOfBlocks( outs.size() );

    int block = 0;
    for ( auto const& [blockName,out] : outs )
    {
        mbds->SetBlock( block, out );
        mbds->GetMetaData(block)->Set(vtkCompositeDataSet::NAME(), blockName.c_str() );
        mbds->GetMetaData(block)->Set(vtkDataObject::DATA_TIME_STEP(), time);
        ++block;
    }
    return mbds;
}

template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::write( int stepIndex, std::string filename, vtkSmartPointer<vtkMultiBlockDataSet> mbds ) const
{
    /*
       out->getoutputinformation(0).set(vtk.vtkstreamingdemanddrivenpipeline.update_number_of_pieces(), this->worldcomm().globalsize());
       out->getoutputinformation(0).set(vtk.vtkstreamingdemanddrivenpipeline.update_piece_number(), this->worldcomm().rank());
       */

    /* InitializeExternal is only supported from 5.10+, */
    /* but lets aim for the latest major version 6 to reduce the complexity */
#if VTK_MAJOR_VERSION >= 6 && defined(VTK_HAS_PARALLEL)
    /* Build vtk objects while reusing the current mpi communicator */
    /* before version 5.10, we cannot initialize a MPIController with an external MPI_Comm */
    vtkSmartPointer<vtkMPICommunicator> mpicomm = vtkSmartPointer<vtkMPICommunicator>::New();
        mpicomm->InitializeExternal(this->opaqueComm);
    vtkSmartPointer<vtkMPIController> mpictrl = vtkSmartPointer<vtkMPIController>::New();
        mpictrl->SetCommunicator(mpicomm);

    vtkSmartPointer<vtkXMLPMultiBlockDataWriter> xmlpw = vtkSmartPointer<vtkXMLPMultiBlockDataWriter>::New();
        xmlpw->SetController(mpictrl);
        //xmlpw->SetTimeStep(stepIndex - TS_INITIAL_INDEX);
        xmlpw->SetFileName(filename.c_str());
        xmlpw->SetInputData(mbds);
        /* only write the meta file on the first processor */
        if( this->worldComm().isMasterRank() )
        { xmlpw->SetWriteMetaFile(1); }
        else
        { xmlpw->SetWriteMetaFile(0); }
        xmlpw->Update();
#else
    vtkSmartPointer<vtkXMLMultiBlockDataWriter> mbw = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
        //mbw->SetTimeStep(stepIndex - TS_INITIAL_INDEX);
        mbw->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION <= 5
        mbw->SetInput(mbds);
#else
        mbw->SetInputData(mbds);
#endif
        mbw->Update();
#endif

    /* basic exporter code */
#if 0
    vtkSmartPointer<vtkoutwriter_type> dw = vtkSmartPointer<vtkoutwriter_type>::New();
#if VTK_MAJOR_VERSION <= 5
        dw->SetInput(out);
#else
        dw->SetInputData(out);
#endif
        dw->SetFileName(filename.c_str());
        dw->Update();
#endif
}

template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::save( steps_write_on_disk_type const& stepsToWriteOnDisk ) const
{
    DVLOG(2) << "[ExporterVTK] save()...\n";

    for ( auto const& [__ts,steps] : stepsToWriteOnDisk  )
    {
        std::map<std::string,vtkSmartPointer<vtkout_type>> outs;
        if ( steps.empty() && __ts->numberOfSteps() == 0 && __ts->hasMesh() ) // save only the mesh
        {
            this->saveMesh( __ts,__ts->mesh(), outs );
        }
        else
        {
            for ( auto const&  __step : steps )
            {
                int stepIndex = __step->index();
                double time = __step->time();

                /* write data into vtk object */
                this->saveMesh( __ts,__step->mesh(), outs );
                this->saveFields( __ts, __step, outs );

#if VTK_MAJOR_VERSION >= 6
                for ( auto & [partname,out] :outs )
                    out->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), time);
#endif

                /* Build a multi block dataset based on gathered data */
                vtkSmartPointer<vtkMultiBlockDataSet> mbds = this->buildMultiBlockDataSet( time, outs );

#if defined(FEELPP_VTK_INSITU_ENABLED)
                if( M_inSituEnable && inSituProcessor->GetNumberOfPipelines() > 0)
                    this->updateInSituProcessor( mbds, stepIndex, time );
#endif

                if ( !M_inSituEnable || M_inSituSave )
                    this->saveData( mbds, stepIndex, time );

            }
        }
    }

    DVLOG(2) << "[ExporterVTK] saving done\n";
}


template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::saveData( vtkSmartPointer<vtkMultiBlockDataSet> mbds, int stepIndex, double time ) const
{
    /* InitializeExternal is only supported from 5.10+, */
    /* but lets aim for the latest major version 6 to reduce the complexity */
    std::ostringstream fname;
    //fname.str("");
    fname << this->path() << "/" << this->prefix()  //<< this->prefix() //this->path()
          << "-" << (stepIndex - TS_INITIAL_INDEX);
#if VTK_MAJOR_VERSION < 6 || !defined(VTK_HAS_PARALLEL)
    fname << "-" << this->worldComm().size() << "_" << this->worldComm().rank();
#endif
    fname << ".vtm";

    /* write VTK files */
    this->write(stepIndex, fname.str(), mbds);

    /* write additional file for handling time steps */
    /* only write on master rank */
    if(this->worldComm().isMasterRank())
    {
        /* rebuild the filename for the pvd file */
        /* This removes the path used for exporting */
        std::ostringstream lfile;
        lfile << this->prefix()  //<< this->prefix() //this->path()
              << "-" << (stepIndex - TS_INITIAL_INDEX);
#if VTK_MAJOR_VERSION < 6 || !defined(VTK_HAS_PARALLEL)
        lfile << "-" << this->worldComm().size() << "_" << this->worldComm().rank();
#endif
        lfile << ".vtm";
                
        /* check if we are on the initial timestep */
        /* if so, we delete the previous pvd file */
        /* otherwise we would append dataset to already existing data */
        std::string pvdFilename = this->path() + "/" + this->prefix() + ".pvd";
        if( (stepIndex - TS_INITIAL_INDEX) == 0 && fs::exists(pvdFilename.c_str()))
        {
            fs::remove(pvdFilename.c_str()); 
        }
#ifdef FEELPP_HAS_LIBXML2
#if VTK_MAJOR_VERSION < 6 || !defined(VTK_HAS_PARALLEL)
        /* when we are not writing data with parallel filters */
        /* we provide the info about the different parts from with */
        /* a dataset is built: the different file names and the part id */
        std::ostringstream oss;
        for(int i = 0; i < this->worldComm().size(); i++)
        {
            oss.str("");
            oss << this->prefix() << "-" << (stepIndex - TS_INITIAL_INDEX)
                << "-" << this->worldComm().size() << "_" << i
                << ".vtm";
            this->writeTimePVD(pvdFilename, time, oss.str(), i);
        }
#else
        /* When writing in parallel, we only write one entry in the pvd file */
        this->writeTimePVD(pvdFilename, time, lfile.str());
#endif
#endif
    }
}
template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::saveMesh( timeset_ptrtype __ts, mesh_ptrtype mesh, std::map<std::string,vtkSmartPointer<vtkout_type>> & outs ) const
{
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
        M_cache_mp[__ts->name()] = std::make_shared<mesh_contiguous_numbering_mapping_type>( mesh.get(), true );
    auto const& mp = *M_cache_mp.find( __ts->name() )->second;

    rank_type currentPid = mesh->worldComm().localRank();

    for ( auto const& [part,nameAndRangeElt] : mp.partIdToRangeElement() )
    {
        vtkSmartPointer<vtkout_type> out = vtkSmartPointer<vtkout_type>::New();

        //! add points to data structure
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetDataTypeToFloat();
        index_type nPointOnProcess = mp.numberOfPoint( part,currentPid );
        points->SetNumberOfPoints( nPointOnProcess );
        index_type startPtId = mp.startPointIds( part,currentPid );
        auto const& nodes = mp.nodes( part ).data();
        for ( index_type i = 0; i < nPointOnProcess ; i++ )
            points->SetPoint( (vtkIdType)(i), nodes + i*3 );

        out->SetPoints(points);

        index_type nElt = mp.numberOfElement( part,currentPid );
        out->Allocate( nElt, nElt );

        //! add cell to data structure
        vtkSmartPointer<vtkelement_type> cell = vtkSmartPointer<vtkelement_type>::New();

        auto const& pointsIdsInElt = mp.pointIdsInElements( part );
        for( index_type i = 0; i < nElt; ++i )
        {
            for( uint16_type p=0; p < mesh_type::element_type::numPoints; ++p )
                cell->GetPointIds()->SetId( p, pointsIdsInElt[i*mesh_type::element_type::numPoints+p] - startPtId );
            out->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
        }

        outs[ mp.name( part ) ] = out;
    }

}

template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::saveFields( timeset_ptrtype __ts, typename timeset_type::step_ptrtype step, std::map<std::string,vtkSmartPointer<vtkout_type>> & outs ) const
{
    auto __mesh = step->mesh();
    auto itFindCache = M_cache_mp.find( __ts->name() );
    CHECK( itFindCache != M_cache_mp.end() ) << "mesh numbering cache not found";
    auto const& mp = *itFindCache->second;
    CHECK( __mesh->isSameMesh( mp.mesh() ) ) << "something wrong mesh should be identical between step and cache";

    for ( auto const& [part,nameAndRangeElt] : mp.partIdToRangeElement() )
    {
        this->saveFields<true>( step, mp, part, step->beginNodal(), step->endNodal(), outs[ mp.name( part ) ] );
        this->saveFields<false>( step, mp, part, step->beginElement(), step->endElement(), outs[ mp.name( part ) ] );
    }

}
template<typename MeshType, int N>
template<bool IsNodal,typename Iterator>
void
ExporterVTK<MeshType,N>::saveFields( typename timeset_type::step_ptrtype step, mesh_contiguous_numbering_mapping_type const& mp, int part, Iterator __var, Iterator en, vtkSmartPointer<vtkout_type> out ) const
{
    rank_type currentPid = mp.mesh()->worldComm().localRank();

    while ( __var != en )
    {
        auto const& fieldData =  __var->second;
        auto const& field00 = unwrap_ptr( fieldData.second[0][0] );
        if ( !field00.worldComm().isActive() ) return;

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

        VLOG(1) << "nComponents field: " << nComponents;

        /* we get that from the local processor */
        /* We do not need the renumbered global index */
        //auto r = markedelements(__mesh, M_markersToWrite[i], EntityProcessType::ALL);
        //auto r = elements( step->mesh() );
        auto const& r =  mp.rangeElement( part );
        auto elt_it = r.begin();
        auto elt_en = r.end();

        vtkSmartPointer<vtkFloatArray> da = vtkSmartPointer<vtkFloatArray>::New();
        da->SetName(__var->first.c_str());
        /* set array parameters */
        /* no need for preallocation if we are using Insert* methods */
        da->SetNumberOfComponents(nComponents);

        float * array = new float[nComponents];
        for ( int k = 0 ; k<nComponents ; ++k )
            array[k] = 0;


        auto const& d = field00.functionSpace()->dof().get();
        //int reorder_tensor2symm[6] = { 0,3,4,1,2,5 };
        int reorder_tensor2symm[6] = { 0,3,5,1,4,2 };

        if constexpr ( IsNodal )
            {
                index_type nValuesPerComponent = mp.numberOfPoint( part,currentPid );
                index_type startPointId = mp.startPointIds( part,currentPid );
                da->SetNumberOfTuples(nValuesPerComponent);

                /* loop on the elements */
                for ( ; elt_it != elt_en; ++elt_it )
                {
                    auto const& elt = unwrap_ref( *elt_it );
                    auto const& locglob_ind = d->localToGlobalIndices( elt.id() );
                    for ( uint16_type p = 0; p < step->mesh()->numLocalVertices(); ++p )
                    {
                        index_type ptid = mp.pointIdToContiguous(part,elt.point( p ).id());
                        DCHECK( ptid != invalid_v<typename mesh_type::index_type> ) << "point not in this process";
                        ptid -= startPointId;

                        size_type dof_id = locglob_ind(d->localDofId(p,0));
                        for ( uint16_type c1 = 0; c1 < fieldData.second.size(); ++c1 )
                        {
                            for ( uint16_type c2 = 0; c2 < fieldData.second[c1].size(); ++c2 )
                            {
                                auto const& fieldComp = unwrap_ptr( fieldData.second[c1][c2] );
                                uint16_type cMap = c2*nComponents1+c1;
                                if ( isTensor2Symm )
                                    cMap = reorder_tensor2symm[Feel::detail::symmetricIndex( c1,c2, nComponents1 )];
                                array[cMap] = fieldComp.globalValue( dof_id );
                            }
                        }
                        da->SetTuple(ptid, array);
                    }
                }

                /* add data array into the vtk object */
                out->GetPointData()->AddArray(da);

                /* Set the first scalar/vector/tensor data, we process as active */
                if ( fieldData.first == FunctionSpaceType::SCALAR && !(out->GetPointData()->GetScalars()) )
                    out->GetPointData()->SetActiveScalars(da->GetName());
                else if ( fieldData.first == FunctionSpaceType::VECTORIAL && !(out->GetPointData()->GetVectors()) )
                    out->GetPointData()->SetActiveVectors(da->GetName());
                else if ( fieldData.first == FunctionSpaceType::TENSOR2 && !(out->GetPointData()->GetTensors()) )
                    out->GetPointData()->SetActiveTensors(da->GetName());
                else if ( fieldData.first == FunctionSpaceType::TENSOR2_SYMM && !(out->GetPointData()->GetTensors()) )
                    out->GetPointData()->SetActiveTensors(da->GetName());

            }
        else
        {
            index_type nValuesPerComponent = std::distance(elt_it,elt_en);
            da->SetNumberOfTuples( nValuesPerComponent);
            index_type startEltId = mp.startElementIds( part,currentPid );
            for ( ; elt_it != elt_en; ++elt_it )
            {
                auto const& elt = unwrap_ref( *elt_it );
                auto const& locglob_ind = d->localToGlobalIndices( elt.id() );
                size_type dof_id = locglob_ind(d->localDofId(0,0));
                index_type e = mp.elementIdToContiguous(part,elt.id());
                e -= startEltId;
                for ( uint16_type c1 = 0; c1 < fieldData.second.size(); ++c1 )
                {
                    for ( uint16_type c2 = 0; c2 < fieldData.second[c1].size(); ++c2 )
                    {
                        auto const& fieldComp = unwrap_ptr( fieldData.second[c1][c2] );
                        uint16_type cMap = c2*nComponents1+c1;
                        if ( isTensor2Symm )
                            cMap = reorder_tensor2symm[Feel::detail::symmetricIndex( c1,c2, nComponents1 )];
                        size_type global_node_id = nComponents * e + cMap;
                        array[cMap] = fieldComp.globalValue( dof_id );
                    }
                }
                da->SetTuple(e, array);
            }
            /* add data array into the vtk object */
            out->GetCellData()->AddArray(da);

            /* Set the first scalar/vector/tensor data, we process as active */
            if ( fieldData.first == FunctionSpaceType::SCALAR && !(out->GetCellData()->GetScalars()) )
                out->GetCellData()->SetActiveScalars(da->GetName());
            else if ( fieldData.first == FunctionSpaceType::VECTORIAL && !(out->GetCellData()->GetVectors()) )
                out->GetCellData()->SetActiveVectors(da->GetName());
            else if ( fieldData.first == FunctionSpaceType::TENSOR2 && !(out->GetCellData()->GetTensors()) )
                out->GetCellData()->SetActiveTensors(da->GetName());
            else if ( fieldData.first == FunctionSpaceType::TENSOR2_SYMM && !(out->GetCellData()->GetTensors()) )
                out->GetCellData()->SetActiveTensors(da->GetName());
        }

        delete[] array;

        DVLOG(2) << "[ExporterVTK::saveNodal] saving " << __var->first << "done\n";

        ++__var;
    }
}

template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::visit( mesh_type* )
{
}

#if defined(FEELPP_VTK_INSITU_ENABLED) 
template<typename MeshType, int N>
void
ExporterVTK<MeshType,N>::updateInSituProcessor( vtkSmartPointer<vtkMultiBlockDataSet> mbds, int stepIndex, double time ) const
{
    //std::cout << "Processing timestep In-Situ " << (__step->index() - TS_INITIAL_INDEX) << " " << __step->time() << std::endl;
    vtkSmartPointer<vtkCPDataDescription> dataDescription = vtkSmartPointer<vtkCPDataDescription>::New();
    dataDescription->AddInput("input");
    dataDescription->SetTimeData(time, stepIndex /*- TS_INITIAL_INDEX*/);

    vtkStdString sh = soption( _name="exporter.vtk.insitu.hostname");

    vtkSmartPointer<vtkStringArray> hname = vtkSmartPointer<vtkStringArray>::New();
    hname->SetName( "hostname" );
    hname->InsertNextValue(sh);
    vtkSmartPointer<vtkIntArray> port = vtkSmartPointer<vtkIntArray>::New();
    port->SetName( "port" );
    port->InsertNextValue( ioption( _name="exporter.vtk.insitu.port") );

    vtkSmartPointer<vtkFieldData> fdata = vtkSmartPointer<vtkFieldData>::New();
    fdata->AddArray(hname);
    fdata->AddArray(port);

    dataDescription->SetUserData(fdata);

    if(inSituProcessor->RequestDataDescription(dataDescription.GetPointer()) != 0)
    {
        dataDescription->GetInputDescriptionByName("input")->SetGrid(mbds);
        //dataDescription->SetForceOutput(true);
        //std::cout << "CoProcess " << inSituProcessor->CoProcess(dataDescription.GetPointer())<< std::endl;
        inSituProcessor->CoProcess(dataDescription.GetPointer());
    }
}
#endif

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
