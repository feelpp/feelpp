#include <fmt/ostream.h>
#include "rht.hpp"

// Readers for view factor data

namespace Feel
{
    // Save and store view factor matrix
    template<int Dim, int Order>
    void RHT<Dim,Order>::saveVF(const std::string& cavity_name, const Eigen::Ref<const Eigen::MatrixXd>& M)
    {
        //Store VF matrix in M_matrix_vf_map
        auto markers_list_vf = M_markers_map[cavity_name];
        auto n_markers = markers_list_vf.size();
        
        // Direct copy: input matrix M already has correct dimensions and ordering
        Eigen::MatrixXd matrix_vf = M;
        M_matrix_vf_map.insert(std::make_pair(cavity_name, matrix_vf));

        // Save the matrix into a CSV file
        std::ofstream matrix_file;
        std::string matrix_filename = fmt::format("VF_Matrix_{}.csv", cavity_name);
        matrix_file.open(matrix_filename, std::ios_base::out);
        
        // Header row
        matrix_file << fmt::format("X,{}\n", fmt::join(markers_list_vf, ","));
        
        // Data rows
        for(size_t i = 0; i < n_markers; ++i) {
            std::vector<double> row_data(matrix_vf.row(i).data(), 
                                          matrix_vf.row(i).data() + n_markers);
            matrix_file << fmt::format("{},{}\n", markers_list_vf[i], fmt::join(row_data, ","));
        }
        matrix_file.close();
        LOG(INFO) << fmt::format("View factor matrix saved to: {}", matrix_filename);
    }

    // Load view factor matrix
    template<int Dim, int Order>
    void RHT<Dim,Order>::loadVF(const std::string& cavity_name, const std::string& filename)
    {
        std::fstream f;
        f.open(Environment::expand( filename));
        if(!f.is_open())
        {
            LOG(ERROR) << fmt::format("Failed to open file: {}", filename);
            throw std::runtime_error(fmt::format("Cannot open file: {}", filename));
        }
        std::vector<std::string> markers;
        std::string line, entry, temp;
        
        // Read first row of the matrix to get the order and number of the markers
        //f >> temp;
        std::getline(f, line,'\n');
        std::stringstream s(line); 
        int number_markers=0;
        while (std::getline(s, entry, ','))
        {
            number_markers++;
        }
        number_markers -= 1; //first column contains no name

        // Read the VF matrix CSV and markers
        // The first column contains always the marker names
        bool is_marker=true;
        int row_mat=0;
        int col_mat=0;
        Eigen::MatrixXd matrix_vf(number_markers,number_markers);
        while (std::getline(f, line,'\n')) 
        {     
            std::stringstream s(line);        
            while (std::getline(s, entry, ',')) {    
                if(is_marker)
                {
                    markers.push_back(entry);
                    is_marker=false;
                }
                else
                {
                    matrix_vf(row_mat,col_mat) = stod(entry);
                    col_mat+=1;   
                }                
            }
            row_mat+=1;
            col_mat=0;
            is_marker=true;

        }
        M_markers_map.insert(std::make_pair(cavity_name,markers));
        M_matrix_vf_map.insert(std::make_pair(cavity_name,matrix_vf));
    }
    template<int Dim, int Order>
    void RHT<Dim,Order>::computeVF(const std::string& cavity_name, const std::string& filename)
    {
        auto jsonfile = removeComments( readFromFile( Environment::expand( filename ) ) );
        std::istringstream astr( jsonfile );
        json json_vf = json::parse( astr );    

        auto markers = json_vf["viewfactor"]["markers"];    
        M_markers_map.insert(std::make_pair(cavity_name,markers));

        if(json_vf["viewfactor"]["type"]=="UnobstructedPlanar")
        {
            UnobstructedPlanarViewFactor<mesh_t> upvf( M_mesh, json_vf );                
            upvf.compute();
            LOG(DEBUG) << fmt::format("View factors:\n{}", upvf.viewFactors());
            saveVF(cavity_name,upvf.viewFactors());
        }
        else if(json_vf["viewfactor"]["type"]=="Raytracing")
        {
            LOG(WARNING) << "Raytracing not implemented yet";
        }
    }
    // Compute the view factor matrix; for the moment, only unobstructed view factor computation
    // via numerical integration
        template<int Dim, int Order>
        void RHT<Dim,Order>::computeVF_and_save()
        {
            if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
            {
                // Loop over cavities: 
                for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
                {
                    auto vf_status = value["viewfactors"]["status"];                
                    if(vf_status=="load")
                    {
                        auto vf_filename = value["viewfactors"]["filename"];
                        LOG(INFO) << fmt::format("Loading view factors for '{}' from {}", bc, vf_filename.dump());
                        loadVF(bc,vf_filename);
                    }
                    else if(vf_status=="compute")
                    {
                        auto vf_filename = value["viewfactors"]["filename"];
                        LOG(INFO) << fmt::format("Computing view factors for '{}' from {}", bc, vf_filename.dump());
                        computeVF(bc,vf_filename);
                    }
                    else
                    {
                        LOG(ERROR) << fmt::format("Invalid view factor status: {}", vf_status.dump());
                    }
                }
            }        

        } // end RHT<Dim,Order>::computeVF_and_save
}