
namespace Feel
{


    template <typename MeshType>
    Eigen::VectorXd
    ShadingMask<MeshType>::get_element_normal(Eigen::VectorXd p1,Eigen::VectorXd p2,Eigen::VectorXd p3)
    {
        Eigen::VectorXd element_normal = ((p3-p1).head<3>()).cross((p2-p1).head<3>());
        element_normal.normalize();
        return(element_normal);
    }


    template <typename MeshType>
    auto
    ShadingMask<MeshType>::commonComputePartCTRL(int NumOption,matrix_node_type const& element_points,int n_rays_thread, int id_thread)
    {
        
        std::vector<double> random_direction(dim);
        matrixSize = M_azimuthSize * M_altitudeSize;
        //NumOption1
        Eigen::MatrixXd SM_table   (M_azimuthSize,M_altitudeSize); SM_table.setZero(); 
        Eigen::MatrixXd Angle_table(M_azimuthSize,M_altitudeSize); Angle_table.setZero(); 
            
        //NumOption2
        Eigen::VectorXd SM_vector    (matrixSize); SM_vector.setZero();
        Eigen::VectorXd Angle_vector (matrixSize); Angle_vector.setZero();
            
        int index_altitude;
        int index_azimuth;
        int initial_index_rays = n_rays_thread * id_thread ;

            for(int j=0;j<n_rays_thread;j++)
            {
                Eigen::VectorXd random_origin= get_random_point(element_points);
                Eigen::VectorXd rand_dir(dim),origin(3);
                bool inward_ray=false;
                Eigen::VectorXd p1(dim),p2(dim),p3(dim);
                for(int i=0;i<dim;i++)
                {
                    p1(i)=column(element_points, 0)[i]; p2(i)=column(element_points, 1)[i]; p3(i)=column(element_points, 2)[i];
                }
                Eigen::VectorXd element_normal=get_element_normal(p1,p2,p3);

                if(dim==3)
                {
                    for(int i=0;i<dim;i++)
                    { 
                        origin(i) = random_origin[i];
                    }

                    // Choose the direction randomly among the latitude and azimuth
                    random_direction = std::get<0>(M_raysdirections[initial_index_rays + j]);
                    index_azimuth    = std::get<1>(M_raysdirections[initial_index_rays + j]);
                    index_altitude   = std::get<2>(M_raysdirections[initial_index_rays + j]);                              
                    for(int i=0;i<dim;i++)              { rand_dir(i) = random_direction[i]; }
                    if(rand_dir.dot(element_normal)>=0) { inward_ray=true; }
                }

                
                //BVHRay<mesh_type::nRealDim> ray( origin, rand_dir, 1e-8 );
                BVHRay<mesh_type::nRealDim> ray( origin, rand_dir);

                int closer_intersection_element = -1;
                if(inward_ray)
                {
                    closer_intersection_element = 1;
                }
                else
                {
                    if (NumOption==1) { 
                        for(auto& [building_name,bvh_building_tree] : M_bvh_tree_vector)
                        {
                            //auto rayIntersectionResult =  bvh_building_tree->intersect(ray) ;
                            auto rayIntersectionResult =  bvh_building_tree->intersect(_ray=ray) ;
                            if ( !rayIntersectionResult.empty() ) closer_intersection_element = 1;
                            if (closer_intersection_element >=0 ) break;
                        }
                    }
                    
                    if (NumOption==2) { 
                        //auto rayIntersectionResult =  M_bvh->intersect(ray) ;
                        auto rayIntersectionResult =  M_bvh->intersect(_ray=ray) ;
                        if ( !rayIntersectionResult.empty() ) closer_intersection_element = 1;    
                    } 
                                
                }
                            
                // If there is an intersection, increase the shading mask table entry by 1 and augment the angle table by 1 as well
                if (NumOption==1) { 
                    if ( closer_intersection_element >=0 )
                    {
                        SM_table(index_azimuth,index_altitude)++; Angle_table(index_azimuth,index_altitude)++;
                    }
                    else
                    {
                        Angle_table(index_azimuth,index_altitude)++;
                    }
                }

                            
                if (NumOption==2) { 
                    int vector_entry = index_azimuth + M_azimuthSize*index_altitude;
                    if ( closer_intersection_element >=0 )
                    {
                        SM_vector(vector_entry)++; Angle_vector(vector_entry)++;
                    }
                    else
                    {
                        Angle_vector(vector_entry)++;
                    }
                }
            }

            auto Result=std::make_pair(SM_table,Angle_table); 

            //if (NumOption==1) { Result=std::make_pair(SM_table,Angle_table);  }
            //if (NumOption==2) { Result=std::make_pair(SM_vector,Angle_vector); }
        return Result; 
    }




/*
   template <typename MeshType>
    auto
    ShadingMask<MeshType>::commonComputePartCTRL_GPU(int NumOption,matrix_node_type const& element_points,int n_rays_thread, int id_thread)
    {
        Eigen::MatrixXd SM_table(M_azimuthSize,M_altitudeSize);
        SM_table.setZero();

        Eigen::MatrixXd Angle_table(M_azimuthSize,M_altitudeSize);
        Angle_table.setZero();

        std::vector<int> indices_altitude(n_rays);
        std::vector<int> indices_azimuth(n_rays);
        std::vector<BVHRay> DeviceRays(n_rays);
        std::vector<int> inward_rays(n_rays, 0);
        int index_altitude, index_azimuth;

        for(int i=0;i<n_rays;i++)
        {

            // Construct the ray emitting from a random point of the element
            auto random_origin = get_random_point(el.second.vertices());

            Eigen::VectorXd rand_dir(dim);
            Eigen::VectorXd p1(dim),p2(dim),p3(dim),origin(3);
            bool inward_ray=false;
                            
            if(dim==3)
            {
                for(int i=0;i<dim;i++)
                {
                    p1(i)=column(el.second.vertices(), 0)[i];
                    p2(i)=column(el.second.vertices(), 1)[i];
                    p3(i)=column(el.second.vertices(), 2)[i];
                    origin(i) = random_origin[i];
                }
                
                auto element_normal = ((p3-p1).head<3>()).cross((p2-p1).head<3>());
                element_normal.normalize();

                // Choose the direction randomly among the latitude and azimuth
                getRandomDirectionSM(random_direction,M_gen,M_gen2,index_azimuth,index_altitude);
                for(int i=0;i<dim;i++)
                {
                    rand_dir(i) = random_direction[i];
                }
                                
                if(rand_dir.dot(element_normal)>=0)
                {
                    inward_rays[i]=1; // true
                }
            }

            BVHRay ray(origin,rand_dir);

            // add the ray to the list of rays
            DeviceRays[i] = ray;
            indices_altitude[i] = index_altitude;
            indices_azimuth[i] = index_azimuth;
        }

        int maxintersection_ray = 10 ; // defines the size of the resulting vectors from the GPU,
        int buildingnumber = M_bvh_tree_vector.size(); // they are of fixed size
        int MaxIntersections = n_rays * maxintersection_ray;

        std::vector<std::vector<int>> results(buildingnumber, std::vector<int>(MaxIntersections, -1));

        int index = 0;
        // Compute the intersection of the rays with each building's bvh tree
        for(auto& [building_name,bvh_building_tree] : M_bvh_tree_vector)
        {
            results[index] = GPUraySearchWrapper(DeviceRays, &bvh_building_tree);
            index++;
        }

        for (int i = 0 ; i < n_rays ; i++)
        {
            if (inward_rays[i] == 1)
            {
                for (int index = 0 ; index < buildingnumber ; index++)
                {
                    results[index][i] = 1;
                }
            }
        }

                        // increase the Shading Tables at the correct indices
        for(int i=0;i<n_rays;i++)
        {
            for(int index = 0 ; index < buildingnumber ; index++)
            {
                if(results[index][i] >= 0)
                {
                    SM_table(indices_azimuth[i],indices_altitude[i]) += 1;
                    Angle_table(indices_azimuth[i],indices_altitude[i]) += results[index][i];
                }
                else
                {
                    Angle_table(indices_azimuth[i],indices_altitude[i]) ++;
                }
            }
        }
        SM_table_marker += SM_table;
        Angle_table_marker += Angle_table;

        auto Result=std::make_pair(SM_table,Angle_table); 

        return Result; 
    }
*/



 // Compute shading masks for one building only

    template <typename MeshType>
    void 
    ShadingMask<MeshType>::computeMasksOneBuildingCTRL(std::string building_name)
    {
        int dim = M_submeshes[building_name]->realDimension();
        std::vector<double> random_direction(dim);
        std::cout << "Size=" << M_submeshes[building_name]->markerNames().size();
        std::cout << " Submeshes markers=" << M_submeshes[building_name]->markerNames() << std::endl;
        //std::cout << "Size=" << M_submeshes.size() << std::endl;
        
        // Loop over the markers of the building
        Eigen::MatrixXd SM_table_marker(M_azimuthSize,M_altitudeSize);
        Eigen::MatrixXd Angle_table_marker(M_azimuthSize,M_altitudeSize);

        for(auto  [marker,marker_id] : M_submeshes[building_name]->markerNames())
        {
            SM_table_marker.setZero();
            Angle_table_marker.setZero();
            auto ray_submesh = createSubmesh(_mesh=M_submeshes[building_name],_range=markedelements(M_submeshes[building_name],marker));

            // Launch Nrays from each triangle of each marker
            for(auto const &el : ray_submesh->elements() ) // from each element of the submesh, launch M_Nrays randomly oriented
            {
                auto rays_from_element = [&](int NumOption,matrix_node_type const& element_points,int n_rays_thread, int id_thread){
                    return commonComputePartCTRL(NumOption,element_points,n_rays_thread,id_thread);
                };

                // Execute the lambda function on multiple threads using
                // std::async and std::future to collect the results
                std::vector<int> n_rays_thread;
                n_rays_thread.push_back(M_Nrays - (M_Nthreads-1) * (int)(M_Nrays / M_Nthreads));
                for(int t= 1; t < M_Nthreads; ++t){ n_rays_thread.push_back( M_Nrays / M_Nthreads); }
                int NumOption=1;
                matrix_node_type const& element_points=el.second.vertices();



                auto MyAlgo000=[&](const int& k) {  
                    auto values_Pair=commonComputePartCTRL(NumOption,element_points,n_rays_thread[k],k);
                    SM_table_marker +=values_Pair.first;
                    Angle_table_marker += values_Pair.second;
                return true;};


            //std::cout<<"\n[INFO: numModeTaskUsed="<<numModeTaskUsed<<"]\n";

            //BEGIN: TaskDispatch part
            /*
                if (numModeTaskUsed==1)
                {
                    TasksDispatch Fg1; 
                    Fg1.init(numTypeThread,M_Nthreads,QSaveTypeThreadDotON);
                    Fg1.run(MyAlgo000);
                }
            */
            //END: TaskDispatch part


            

             //BEGIN: Task part

             
                if (numModeTaskUsed==2)
                {
                    Task::Task TsK(M_Nthreads,numTypeThread);
                    TsK.setSave( QSaveTypeThreadDotON );
                    TsK.setInfo( false );
                    for ( int k = 0; k < M_Nthreads; k++ )
                    {
                        auto const& idk = k;
                        TsK.add( _param(idk), _tasks = MyAlgo000 );
                    }
                    TsK.run();
                    TsK.close();
                }


                if (numModeTaskUsed==3)
                {
                    Task::Task TsK(M_Nthreads,numTypeThread);
                    TsK.setSave( QSaveTypeThreadDotON );
                    TsK.setInfo( false );
                    TsK.add_loop(MyAlgo000,M_Nthreads);
                    TsK.run();
                    TsK.close();
                }
            
            
            //END: Task part


            }
            // Divide the shading mask by the corresponding value of the angle table
            // If an angle combination has not been selected, suppose there is no shadow
            auto shadingMatrix = SM_table_marker.array().binaryExpr( Angle_table_marker.array() , [](auto x, auto y) { return y==0 ? 0 : x/y; });

            // Shading mask value 0 means that the surface is not shadowed, value 1 it is fully shadowed
            // Save the shading mask table to a csv file
            if (M_saveMasks) 
            { 
                saveShadingMask("SM_Matrix_",building_name,marker,shadingMatrix.matrix());
                if (QSaveControlFiles) { 
                    saveShadingMask("SM_Matrix_CTRL_",building_name,marker,shadingMatrix.matrix());
                 }              
            }
        }
    }


template <typename MeshType>
    void 
    ShadingMask<MeshType>::saveMetadataInfoPart()
    {
        //BEGIN:SAVE META INFO
        auto timeComputation = toc("Shading masks computed using raytracing");
        M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
        M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
        M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;

        auto end_computation = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end_computation);
        M_metadataJson["shadingMask"]["Timestamp"]["End"] = strtok(std::ctime(&end_time),"\n");

        auto timeAllDuration = toc("Time All Duration");
        M_metadataJson["shadingMask"]["Timestamp"]["Time All Duration"] = timeAllDuration;  

        saveMetadata("shadingmask_metadata_"+std::to_string(M_Nthreads)+"_"+std::to_string(M_Nrays)); 
        std::cout<<"[INFO: Metadata Saved]\n";
        //END:SAVE META INFO
    }



template <typename MeshType>
    void 
    ShadingMask<MeshType>::computeMasksSubPartList()
    {
        // [INFO]: refactoring OK for this parts
        tic();
        int nbObjects = 0;
        std::string nameFile,buildingName;
        std::vector<std::string> listObjects;
        nl::json const&  markersVolume = j_["Buildings"]["list"].get<std::vector<std::string>>(); nbObjects=markersVolume.size(); 

        std::cout << "[computeMasksSubPartList] : nbObjects="<<nbObjects<< std::endl;
        std::cout << "[computeMasksSubPartList] : M_Nthreads =" <<M_Nthreads<< std::endl;

        for(int idx = 0 ; idx <nbObjects; ++idx)
        {
            buildingName=markersVolume[idx];
            computeMasksOneBuildingCTRL(buildingName);
        }
        
        //BEGIN:SAVE META INFO
        saveMetadataInfoPart();
        //END:SAVE META INFO

    }





 template <typename MeshType>
    void 
    ShadingMask<MeshType>:: computeMasksSubPartSurfaceVolumes(int numOp)
    {
        // [INFO]: refactoring OK for this parts
        tic();
        int nbObjects = 0;
        std::string nameFile,buildingName;
        std::vector<std::string> listObjects;
        
        if (numOp==2) { nameFile=Environment::expand(j_["Buildings"]["fileSurfaces"].get<std::string>()); }
        if (numOp==3) { nameFile=Environment::expand(j_["Buildings"]["fileVolumes"].get<std::string>());  }
 
        if ((numOp==2) || (numOp==3)) { listObjects=GetListNameObjects(nameFile); nbObjects=listObjects.size(); }

        std::cout << "[computeMasksSubPartSurfaceVolumes] : nbObjects="<<nbObjects<< std::endl;
        std::cout << "[computeMasksSubPartSurfaceVolumes] : M_Nthreads =" <<M_Nthreads<< std::endl;

        /**
        if (numTypeThread==2) {
            SpRuntime runtime(M_Nthreads);
            for(int idx = 0 ; idx <nbObjects; ++idx)
            {
                runtime.task(SpRead(idx),
                        [&](const int & k) -> bool {
                            if ((numOp==2)|| (numOp==3)) { buildingName=listObjects[k]; }
                            computeMasksOneBuildingCTRL(buildingName);
                            return true;
                        }
                ).setTaskName("Op("+std::to_string(idx)+")");
                usleep(10);
                std::atomic_int counter(0);
            }
            runtime.waitAllTasks();
            runtime.stopAllThreads();
        }
        else 
        */
        {
            for(int idx = 0 ; idx <nbObjects; ++idx)
            {
                if ((numOp==2)|| (numOp==3)) { buildingName=listObjects[idx];   }
                computeMasksOneBuildingCTRL(buildingName);
            }
        }
        
        //BEGIN:SAVE META INFO
        saveMetadataInfoPart();
        //END:SAVE META INFO
    }


template <typename MeshType>
    bool
    ShadingMask<MeshType>::computePartMarker(
        std::vector<std::string> marker_list_thread, 
        int id_thread, int start_index)
    {
        std::vector<double> random_direction(dim);
        int index_altitude;
        int index_azimuth;
                            
        int initial_index_marker;
        int i_marker = 0;

        int len_marker_list_thread = marker_list_thread.size();

        int vector_entry;

        for( auto const& marker : marker_list_thread)
        {
            auto faces_with_marker = M_listMarkerFaceEntity[marker];
                                
            initial_index_marker = start_index + i_marker;

            auto initial_index_SM = SM_tables_Alpha.begin() +  initial_index_marker * matrixSize;
            auto initial_index_Angles = Angle_tables_Alpha.begin() +  initial_index_marker * matrixSize;

                                // Extract a view from the vectors SM_tables and Angle_tables
            auto SM_vector = Eigen::Map<Eigen::VectorXd>( &(*initial_index_SM), matrixSize);
            auto Angle_vector = Eigen::Map<Eigen::VectorXd>( &(*initial_index_Angles), matrixSize);


            for(auto const& face : faces_with_marker)
            {
                for(int j=0;j<M_Nrays;j++)
                {
                    // Construct the ray emitting from a random point of the element
                    auto random_origin = get_random_point(face.vertices());

                    Eigen::VectorXd rand_dir(3);
                    Eigen::VectorXd p1(3),p2(3),p3(3),origin(3);
                    bool inward_ray=false;
                                        
                     for(int i=0;i<3;i++)
                    {
                        p1(i)=column(face.vertices(), 0)[i];
                        p2(i)=column(face.vertices(), 1)[i];
                        p3(i)=column(face.vertices(), 2)[i];
                        origin(i) = random_origin[i];
                    }
                    
                    Eigen::VectorXd element_normal=get_element_normal(p1,p2,p3);

                    // Choose the direction randomly among the latitude and azimuth
                                        
                    random_direction = std::get<0>(M_raysdirections[j]);
                    index_azimuth = std::get<1>(M_raysdirections[j]);
                    index_altitude = std::get<2>(M_raysdirections[j]);
                    for(int i=0;i<3;i++)
                    {
                        rand_dir(i) = random_direction[i];
                    }
                    if(rand_dir.dot(element_normal)>=0)
                    {
                        inward_ray=true;
                    }

                    //BVHRay<mesh_type::nRealDim> ray( origin, rand_dir, 1e-8 );
                    BVHRay<mesh_type::nRealDim> ray( origin, rand_dir);

                    int closer_intersection_element = -1;
                    if(inward_ray)
                    {
                        closer_intersection_element = 1;
                    }
                    else
                    {
                        //auto rayIntersectionResult =  M_bvh->intersect(ray) ;
                        auto rayIntersectionResult =  M_bvh->intersect(_ray=ray) ;
                        if ( !rayIntersectionResult.empty() )
                        closer_intersection_element = 1;                                
                    }
                    
                    // Compute the index associated to the entry to modify
                    // The vector is constituted of M_altitudeSize blocks of M_azimuthSize stacked onto each other
                    vector_entry = index_azimuth + M_azimuthSize*index_altitude;

                    // If there is an intersection, increase the shading mask table entry by 1 and augment the angle table by 1 as well
                    if ( closer_intersection_element >=0 )
                    {
                        SM_vector(vector_entry)++;
                        Angle_vector(vector_entry)++;
                    }
                    else
                    {
                        Angle_vector(vector_entry)++;
                    }
                }
             }

            i_marker += 1;
            // std::cout << "I_marker " << i_marker << " thread number " << id_thread << " marker " << marker << std::endl;
        }
        return true;
    }


template <typename MeshType>
    void 
    ShadingMask<MeshType>::computeMasksSubPartMarkersCTRL()
    {        
            std::vector<double> random_direction(3);            
            std::map<std::string, int> markerLineMap;
            matrixSize = M_azimuthSize * M_altitudeSize;

            SM_tables_Alpha.assign(M_listFaceMarkers.size() * matrixSize, 0);
            Angle_tables_Alpha.assign(M_listFaceMarkers.size() * matrixSize, 0);

            std::cout << std::endl;
            std::cout << "M_Nthreads="<<M_Nthreads<< std::endl;
            std::cout << "Allocated SM_tables and Angle_tables of size " << M_listFaceMarkers.size() * matrixSize << std::endl;
            std::cout << " - M_listFaceMarkers.size " << M_listFaceMarkers.size()<< std::endl;
            std::cout << " - matrixSize " << matrixSize << std::endl;
            
                    auto multithreading_over_markers = [&](std::vector<std::string> marker_list_thread, int id_thread, int start_index) {
                        computePartMarker(marker_list_thread,id_thread,start_index);
                        return true;
                    };

                    
                    tic();
                    // Store the index of M_listFaceMarkers where each thread will stop its computations
                std::vector<int> marker_threads_list_length;

                marker_threads_list_length.push_back(M_listFaceMarkers.size() - (M_Nthreads-1) * (int)(M_listFaceMarkers.size() / M_Nthreads));
                for(int t= 1; t < M_Nthreads; ++t){
                    marker_threads_list_length.push_back( marker_threads_list_length[t-1] + M_listFaceMarkers.size() / M_Nthreads);
                }

                    std::cout << std::endl;
                    std::cout << "Size of marker list per thread " << marker_threads_list_length << std::endl;
                    std::cout << std::endl;
                    std::cout << "---marker_threads_list_length   = " <<marker_threads_list_length.size()<< std::endl;
                    std::cout << "---M_listFaceMarkers            = " <<M_listFaceMarkers.size()<< std::endl;
                    std::cout << "---M_listFaceMarkers/M_Nthreads = " <<M_listFaceMarkers.size()/M_Nthreads<< std::endl;
                    std::cout << std::endl;

                
                std::vector<std::vector<std::string>> marker_thread_lists(marker_threads_list_length.size());    
                // Store the index of M_listFaceMarkers from where each thread will start its computations
                std::vector<int> start_index_list(marker_threads_list_length.size());
                start_index_list[0]=0;
                int n0 = 0;
                int t = 0;                 
                for(auto n : marker_threads_list_length)
                {
                    for(int i=n0; i< n; i++) { marker_thread_lists[t].push_back(M_listFaceMarkers[i]);}
                    n0 = n;
                    t += 1;
                    start_index_list[t]=n;
                }

                
                auto MyAlgo000=[&](const int& k) {  
                    multithreading_over_markers(marker_thread_lists[k],k,start_index_list[k]);
                return true;};

                //std::cout<<"\n[INFO: numModeTaskUsed="<<numModeTaskUsed<<"]\n";
                std::cout<<"\n[INFO: M_Nthreads="<<M_Nthreads<<"]\n";

                //BEGIN: TaskDispatch part
                /*
                if (numModeTaskUsed==1)
                {
                    TasksDispatch Fg1; 
                    Fg1.init(numTypeThread,M_Nthreads,QSaveTypeThreadDotON);
                    Fg1.run(MyAlgo000);
                }
                */
                
                //END: TaskDispatch part

                //BEGIN: Task part

                
                if (numModeTaskUsed==2)
                {
                    Task::Task TsK(M_Nthreads,numTypeThread);
                    TsK.setSave( QSaveTypeThreadDotON );
                    TsK.setInfo( false );
                    for ( int k = 0; k < M_Nthreads; k++ )
                    {
                        auto const& idk = k;
                        TsK.add( _param(idk), _tasks = MyAlgo000 );
                    }
                    TsK.run();
                    TsK.close();
                }
                

                if (numModeTaskUsed==3)
                {
                    Task::Task TsK(M_Nthreads,numTypeThread);
                    TsK.setSave( QSaveTypeThreadDotON );
                    TsK.setInfo( false );
                    TsK.add_loop(MyAlgo000,M_Nthreads);
                    TsK.run();
                    TsK.close();
                }

                
                //END: Task part


                auto timeComputation = toc("Shading masks computed using raytracing");
                M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
                M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
                M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;



            // Divide the shading mask by the corresponding value of the angle table
            // If an angle combination has not been selected, suppose there is no shadow
            std::transform(SM_tables_Alpha.begin(),SM_tables_Alpha.end(),Angle_tables_Alpha.begin(),SM_tables_Alpha.begin(),std::divides<double>());  
            computeSaveMasks(SM_tables_Alpha);      
    }


// Compute shading masks for the buildings in the json file
template <typename MeshType>
    void 
    ShadingMask<MeshType>::computeMasksMaster()
    {
        if (numTypeThread==2) { std::cout<<"\n[INFO] : Resolution with Specx...\n"; }
        else { std::cout<<"\n[INFO] : Resolution with std::async...\n"; }

        if ( j_["/Buildings"_json_pointer].contains("list") )         { std::cout << "\n[INFO: STEP2 >COMPUTE MASKS MASTER LIST]" << std::endl;     computeMasksSubPartList();  }
        if ( j_["/Buildings"_json_pointer].contains("fileVolumes"))   { std::cout << "\n[INFO: STEP2 >COMPUTE MASKS MASTER VOLUMES]" << std::endl;  computeMasksSubPartSurfaceVolumes(3);  }
        if ( j_["/Buildings"_json_pointer].contains("fileSurfaces") ) { std::cout << "\n[INFO: STEP2 >COMPUTE MASKS MASTER SURFACES]" << std::endl; computeMasksSubPartSurfaceVolumes(2);  }
        if ( j_["/Buildings"_json_pointer].contains("fileFaces") ||  j_["/Buildings"_json_pointer].contains("aggregatedMarkers") ) 
        {            
             if (M_mthreadtype == "markers") { std::cout << "\n[INFO:STEP2 >COMPUTE MASKS MASTER AGGREGATE MARKERS]" << std::endl; computeMasksSubPartMarkersCTRL();  }
        }
    }



 template <typename MeshType>
    void 
    ShadingMask<MeshType>::computeSaveMasks(std::vector<double> SM_tables)
    {
        // [INFO]: refactoring OK for this parts
        std::cout << "\n[INFO: STEP3 >SAVE]\n";
        if(M_saveMasks)
        {
            tic();
            int numOp = 0;
            std::string building_name, marker="";
            if (j_["/Buildings"_json_pointer].contains("fileFaces"))         { numOp=1; }
            if (j_["/Buildings"_json_pointer].contains("aggregatedMarkers")) { numOp=2; }
            // a csv containing the face markers is provided, or they are computed using aggregated markers

            for(int i=0; i< M_listFaceMarkers.size(); i++)
            {
                if (numOp==1) { building_name = std::to_string(i); marker = std::to_string(i); }
                if (numOp==2) { building_name = M_listFaceMarkers[i]; }
                auto initial_index_SM = SM_tables.begin() +  i * matrixSize;
                auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),M_azimuthSize, M_altitudeSize);   
                saveShadingMask("SM_Matrix_",building_name,marker,shadingMatrix.matrix());
                if (QSaveControlFiles) { saveShadingMask("SM_Matrix_CTRL_",building_name,marker,shadingMatrix.matrix()); }
            }
  
            auto timeCSVsaving = toc("Mask CSV saved");
            M_metadataJson["shadingMask"]["Timer"]["MaskCSVsaving"] = timeCSVsaving;     
        }
        auto end_computation = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end_computation);
        M_metadataJson["shadingMask"]["Timestamp"]["End"] = strtok(std::ctime(&end_time),"\n");
        std::time_t end_time_all = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        
        auto timeAllDuration = toc("Time All Duration");
        M_metadataJson["shadingMask"]["Timestamp"]["Time All Duration"] = timeAllDuration;  

        saveMetadata("shadingmask_metadata_"+std::to_string(M_Nthreads)+"_"+std::to_string(M_Nrays)); 
        std::cout<<"[INFO: Metadata Saved]\n";
    }

}