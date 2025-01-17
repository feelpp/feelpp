namespace Feel
{

    template <typename MeshType>
    void
    ShadingMask<MeshType>::saveShadingMask(std::string prefix_name,std::string building_name, std::string marker_name, const Eigen::Ref<const Eigen::MatrixXd>& M)
    {
        int i_long=0;
        int j_alt=0;

        Eigen::MatrixXd matrix_sm(M_azimuthSize,M_altitudeSize);

        // Save the matrix into a CSV file, inside the shadingMasks subfolder of the results folder
        std::ofstream matrix_file;
        std::string shadingMaskFolder = (fs::path(Environment::appRepository())/("shadingMasks")).string();
        if (!fs::exists(shadingMaskFolder))
            fs::create_directory(shadingMaskFolder);
        std::string matrix_filename = shadingMaskFolder+"/"+prefix_name+building_name+"_"+marker_name+".csv";
        matrix_file.open(matrix_filename,std::ios_base::trunc);
        for(int i=0; i<M_azimuthSize; i++)
        {
            for(int j=0; j<M_altitudeSize-1; j++)
            {
                matrix_file << M(i,j) << ",";
            }
            matrix_file << M(i,M_altitudeSize-1) << "\n";
        }
        matrix_file.close();
    }



    template <typename MeshType>
    void
    ShadingMask<MeshType>::saveMetadata(std::string name)
    {
        std::string folder = (fs::path(Environment::appRepository())).string();
        if (!fs::exists(folder)) { fs::create_directory(folder); }
        std::string json_filename = folder+"/"+name+".json";
        std::ofstream o(json_filename);
        o << std::setw(6) << M_metadataJson << std::endl;
    }
}
