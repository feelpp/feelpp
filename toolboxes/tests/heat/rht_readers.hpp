#include <fmt/ostream.h>
#include "rht.hpp"

// Readers for view factor data

namespace Feel
{
    // Save and store view factor matrix
    template<int Dim, int Order>
    void RHT<Dim,Order>::saveVF(std::string cavity_name,const Eigen::Ref<const Eigen::MatrixXd>& M)
    {
        //Store VF matrix in M_matrix_vf_map
        int i_mark=0;
        int j_mark=0;
        auto markers_list_vf=M_markers_map[cavity_name];
        auto n_markers=markers_list_vf.size();
        Eigen::MatrixXd matrix_vf(n_markers,n_markers);
        for(auto & marker1 : markers_list_vf)
        {
            // Compute the first index to place view factors in the vf matrix
            int index1;
            auto it = std::find(markers_list_vf.begin(),markers_list_vf.end(),marker1);
            if( it != markers_list_vf.end())
            {
                index1= it -  markers_list_vf.begin();
            }
            else
            {
                std::cout << fmt::format("Marker {} not in the list",marker1)<<std::endl;
                break;
            }

            for(auto & marker2 :  markers_list_vf)
            {
                // Compute the second index to place view factors in the vf matrix
                int index2;
                auto it = std::find(markers_list_vf.begin(),markers_list_vf.end(),marker2);
                if( it != markers_list_vf.end())
                {
                    index2= it -  markers_list_vf.begin();
                }
                else
                {
                    std::cout << fmt::format("Marker {} not in the list",marker1)<<std::endl;
                    break;
                }                   
                // Dispatch the newly computed view factors to the matrix from before and save this matrix
                matrix_vf(index1,index2)=M(i_mark,j_mark);

                j_mark+=1;
            }
            j_mark=0;
            i_mark+=1;
        }
        M_matrix_vf_map.insert(std::make_pair(cavity_name,matrix_vf));

        // Save the matrix into a CSV file
        std::ofstream matrix_file;
        std::string matrix_filename="VF_Matrix_"+cavity_name+".csv";
        matrix_file.open(matrix_filename,std::ios_base::out);
        for(int i=0; i<markers_list_vf.size(); i++)
        {
            if(i==0)
            {
                matrix_file << " X,";
                for(int j=0; j<markers_list_vf.size()-1; j++)
                {
                    matrix_file << markers_list_vf[j] << ",";
                }
                matrix_file << markers_list_vf[markers_list_vf.size()-1] << "\n";
            }
            if(markers_list_vf.size()==1)
                matrix_file << markers_list_vf[i] << ",";
            else
            {
                for(int j=0; j<markers_list_vf.size()-1; j++)
                {
                    if(j==0)
                    {                    
                        matrix_file << markers_list_vf[i] << ",";                                 
                    }
                    matrix_file << matrix_vf(i,j) << ",";
                }        
            }
            matrix_file << matrix_vf(i,markers_list_vf.size()-1) << "\n";
        }
        matrix_file.close();
        LOG(INFO) << fmt::format("View factor matrix has been saved and stored") << std::endl;
    }

    // Load view factor matrix
    template<int Dim, int Order>
    void RHT<Dim,Order>::loadVF(std::string cavity_name, std::string filename)
    {
        std::fstream f;
        f.open(Environment::expand( filename));
        if(!f.is_open())
        {
            std::cout << "file not opened" << std::endl;
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
    void RHT<Dim,Order>::computeVF(std::string cavity_name,std::string filename)
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
            std::cout << upvf.viewFactors() << std::endl;
            saveVF(cavity_name,upvf.viewFactors());
        }
        else if(json_vf["viewfactor"]["type"]=="Raytracing")
        {
            std::cout << "Raytracing not implemented at the moment" <<std::endl;
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
                        fmt::format("Loading view factors for {} from {}",bc,vf_filename.dump());
                        loadVF(bc,vf_filename);
                    }
                    else if(vf_status=="compute")
                    {
                        auto vf_filename = value["viewfactors"]["filename"];
                        fmt::format( "Computing view factors for {} from {}", bc, vf_filename.dump() );
                        computeVF(bc,vf_filename);
                    }
                    else
                    {
                        std::cout << "View factor status not correct" << std::endl;
                    }
                }
            }        

        } // end RHT<Dim,Order>::computeVF_and_save
}