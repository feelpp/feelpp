/**
 * @file rht.cpp
 * @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 * @brief Radiative heat transfer
 * @version 0.1
 * @date 2022-07-22
 *
 * @copyright Copyright (c) 2022 Université de Strasbourg
 *
 */
#include <iostream>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelviewfactor/unobstructedplanarviewfactor.hpp>
#include <assert.h> 

namespace Feel
{      
    template<int Dim,int Order=1>
    class RHT  
    {
    public:
        using mesh_t = Mesh<Simplex<Dim>>;
        using mesh_ptr_t = std::shared_ptr<mesh_t>;
        using mesh_trace_t = typename mesh_t::trace_mesh_type;
        using mesh_trace_ptr_t = typename mesh_t::trace_mesh_ptrtype;
        using spacevect_t = Pchv_type<mesh_t, Order>;
        using spacevect_ptr_t = Pchv_ptrtype<mesh_t, Order>;
        using space_t = Pch_type<mesh_t, Order>;
        using space_ptr_t = Pch_ptrtype<mesh_t, Order>;
        using spacedisc_t = Pdh_type<mesh_trace_t, Order-1>;
        using spacedisc_ptr_t = Pdh_ptrtype<mesh_trace_t, Order-1>;
        using spacedisc_interm_ptr_t = Pdh_ptrtype<mesh_trace_t, Order>;
        static constexpr int nDim = Dim;
        using scalar_t =  double;
        using coord_t = Eigen::Matrix<double, Dim, 1>;
        using tensor_t = Eigen::Matrix<double, Dim, Dim>;
        using elementvect_t = typename spacevect_t::element_type;
        using elementvect_ptr_t = typename spacevect_t::element_ptrtype;
        using element_t = typename space_t::element_type;
        using element_ptr_t = typename space_t::element_ptrtype;
        using elementdisc_ptr_t = typename spacedisc_t::element_ptrtype;

        RHT(nl::json specs)
        {
            // Assign the json structures to the members of the class
            this->specs=specs;            
        }    

        void init();
        void solveHeatEquation(element_ptr_t T, element_ptr_t q );

        typedef Backend<double> backend_type;
        typedef std::shared_ptr<backend_type> backend_ptrtype;

        /*matrix*/
        typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
        typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
        typedef typename backend_type::vector_type vector_type;
        typedef typename backend_type::vector_ptrtype vector_ptrtype;

        typedef Bdf<space_t>  bdf_type;
        typedef std::shared_ptr<bdf_type> bdf_ptrtype;

        typedef Exporter<mesh_t,1> exporter_type;
        typedef std::shared_ptr <exporter_type> exporter_ptrtype;

        void rebuild_D(element_ptr_t T);
        void rebuild_N(element_ptr_t q);
        void build_M();
        void execute();
        void computeVF_and_save();
        void saveVF(std::string cavity_name,const Eigen::Ref<const Eigen::MatrixXd>&  M);
        void loadVF(std::string cavity_name,std::string filename);
        void computeVF(std::string cavity_name,std::string filename);
        void checkResults();

        struct TandQ
        {
            element_ptr_t T_;
            element_ptr_t q_;

            element_ptr_t T(){return T_;}
            element_ptr_t q(){return q_;}

            // void setT(element_ptr_t& T){T_=boost::unwrap_ref(T);}
            // void setq(element_ptr_t& q){q_=boost::unwrap_ref(q);}
            void setT(element_ptr_t& T){T_=T;}
            void setq(element_ptr_t& q){q_=q;}
        };

        void solvePicardIteration();
        void exportHeat();
        void initHeatEquation();        
        nl::json specs;
        nl::json j_viewfactor;
        space_ptr_t M_Xh;
        spacevect_ptr_t M_Xhvec;
        spacedisc_ptr_t M_Xhd0;
        spacedisc_interm_ptr_t M_Xhd1;
        std::map< std::string, spacedisc_ptr_t > M_Xhd0_map;
        std::map< std::string, spacedisc_interm_ptr_t > M_Xhd1_map;
        mesh_ptr_t M_mesh;
        std::map< std::string, mesh_trace_ptr_t> M_surface_submesh_map;
        std::map< std::string, std::vector<std::string> > M_markers_map;   
        bdf_ptrtype M_bdf;
        TandQ M_currentTempAndFlux;        

        sparse_matrix_ptrtype M_M, M_M_block; // matrix for problem on multiple irradiating surfaces
        std::map< std::string, Eigen::MatrixXd > M_matrix_vf_map; // matrix storing view factors
        vector_ptrtype M_D, M_N, M_D_block,M_N_block; // vectors for problem on multiple irradiating surfaces

        sparse_matrix_ptrtype M_a,M_at; // matrices for heat transfer PDE
        vector_ptrtype M_l,M_lt; // right-hand sides for heat transfer PDE

        exporter_ptrtype M_e; // BDF exporter
    };

// Add options to the main file 
    inline Feel::po::options_description
    makeOptions()
    {
        Feel::po::options_description options( "rht options" );
        options.add_options()
            ( "specs", Feel::po::value<std::string>(),"json spec file for rht" )
            ( "steady", Feel::po::value<bool>()->default_value( 1 ),"if 1: steady else unsteady" )
            ("deactivate-exporters",Feel::po::value<bool>()->default_value( false ),"deactivate exporters");
        options.add( Feel::backend_options("interp1"));
        options.add( Feel::backend_options("interp2"));
        options.add( Feel::backend_options("heatEq"));
        options.add( Feel::backend_options("Modesteq"));

        return options.add( Feel::feel_options() );
    
    } // end makeOptions

// Build the time-invariant part of the heat equation
    template<int Dim, int Order>
    void RHT<Dim,Order>::initHeatEquation()
    {
        auto a = form2( _test = M_Xh, _trial = M_Xh,_matrix=M_a,_backend=backend(_name="heatEq") );    
        auto l = form1( _test = M_Xh,_vector=M_l );

        auto u = M_Xh->element();
        auto v = M_Xh->element();

        M_e=exporter(_mesh=M_mesh,_name="HeatEqn");

        M_bdf->initialize( unwrap_ptr(M_currentTempAndFlux.T()) );
        M_bdf->start();    
        
        // For each material, integrate rho*C*idt(u)*id(v)/dt + \int_MAT k_MAT * grad(u) * grad(v)
        for ( auto [key, material] : specs["/Models/heat/materials"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "material {}", material );
            std::string mat = fmt::format( "/Materials/{}/k", material.get<std::string>() );
            auto k = specs[nl::json::json_pointer( mat )].get<std::string>();
            std::string matRho = fmt::format( "/Materials/{}/rho", material.get<std::string>() );
            std::string matCp = fmt::format( "/Materials/{}/Cp", material.get<std::string>() );
            auto Rho = specs[nl::json::json_pointer( matRho )].get<std::string>();
            auto Cp = specs[nl::json::json_pointer( matCp )].get<std::string>();

            a += integrate( _range = markedelements( M_mesh, material.get<std::string>() ), 
                    _expr =  inner( expr( k ) * gradt( u ) , grad( v ) ) );
            a += integrate(_range = markedelements( M_mesh, material.get<std::string>() ), 
                    _expr = M_bdf->polyDerivCoefficient( 0 ) *expr( Rho ) * expr( Cp ) * idt( u ) * id( v ) );
        }  

        // BC Robin -grad(u)*N = h*(T-Text)
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "convective_heat_flux" ) )
        {
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/convective_heat_flux"_json_pointer].items() )
            {
                LOG( INFO ) << fmt::format( "convective_heat_flux {}: {}", bc, value.dump() );
                auto h = value["h"].get<std::string>();
                auto Text = expr(value["Text"].get<std::string>());

                Text.setParameterValues({{"Tref_C",specs["/Parameters/Tref_C"_json_pointer].get<double>()}});

                a += integrate( _range = markedfaces( M_mesh, bc ),
                        _expr = expr( h ) * id( v ) * idt( u ) );
                l += integrate( _range = markedfaces( M_mesh, bc ),
                        _expr =  expr( h ) * Text * id( v ) );
            }
        }

        // ===== BC RHT =====
        // Blackbody radiative condition (no view factors)
        // set the term sigma * epsilon * T_0^4
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_blackbody_heat_flux" ) )
        {
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_blackbody_heat_flux"_json_pointer].items() )
            {
                LOG( INFO ) << fmt::format( "radiative_blackbody_heat_flux {}: {}", bc, value.dump() );
                            
                auto sigma = expr(value["sigma"].get<std::string>());                
                auto Tref = expr(value["Tref"].get<std::string>());
                Tref.setParameterValues({{"Tref_C",specs["/Parameters/Tref_C"_json_pointer].get<double>()}});
                sigma.setParameterValues({{"sigma",specs["/Parameters/sigma"_json_pointer].get<double>()}});

                auto Tref4 = Tref*Tref*Tref*Tref;

                // Recover the emissivity epsilon of the coating material associated to the emitting face
                for ( std::string mark :  value.at("markers") )
                {
                    auto done = 0;
                    for ( auto& [key, coat] : specs["/Coating"_json_pointer].items() )
                    {
                        for ( auto markcoat : coat.at("markers") )
                        {                        
                            if ( mark == markcoat )
                            {
                                std::cout << mark << " "<<markcoat <<std::endl;
                                auto epsilon = coat["epsilon"].get<std::string>();

                                l += integrate( _range = markedfaces( M_mesh, mark ), 
                                        _expr =  sigma * expr( epsilon ) * Tref4 * id( v ) );
                                done = 1;
                                break;
                            }
                        }
                        if ( done ){ break; }
                    }
                }
            }
        }   
    } // end RHT<Dim,Order>::initHeatEquation

// Export the temperature and the radiative flux
    template<int Dim, int Order>
    void RHT<Dim,Order>::exportHeat()
    {               
        M_e->step(M_bdf->time())->add( "T", idv(M_currentTempAndFlux.T()) );
        M_e->step(M_bdf->time())->add( "q", idv(M_currentTempAndFlux.q()) );
        M_e->save();
    } // end RHT<Dim,Order>::exportHeat()

// Solve the heat transfer problem
    template<int Dim, int Order>
    void RHT<Dim,Order>::solveHeatEquation( element_ptr_t T, element_ptr_t q )
    {
        // Create the linear and bilinear forms
        auto at = form2(_trial=M_Xh, _test= M_Xh, _matrix=M_at,_backend=backend(_name="heatEq",_rebuild=true));
        auto lt = form1(_test=M_Xh, _vector=M_lt);

        auto a = form2(_trial=M_Xh, _test= M_Xh, _matrix=M_a,_backend=backend(_name="heatEq"));
        auto l = form1(_test=M_Xh, _vector=M_l);

        // Start from the time-independent forms computed in initHeatEquation
        at = a;
        lt = l;

        auto u = M_Xh->element();
        auto v = M_Xh->element();

        auto T_sol = M_Xh->elementPtr();

        // Save T_{n-1} in another variable to avoid its modification 
        *T_sol = *T;

        // BC Neumann
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "flux" ) )
        {
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/flux"_json_pointer].items() )
            {
                LOG( INFO ) << fmt::format( "flux {}: {}", bc, value.dump() );
                auto flux = expr(value["expr"].get<std::string>());

                lt += integrate( _range = markedfaces( M_mesh, bc ),
                        _expr = - flux  * id( v ) );
            }
        }

        // Update of RHS for heat equations starting from T with the term rho*C*T_{n-1}/dt (at order 1)
        for ( auto [key, material] : specs["/Models/heat/materials"_json_pointer].items() )
        {
            std::string matRho = fmt::format( "/Materials/{}/rho", material.get<std::string>() );
            std::string matCp = fmt::format( "/Materials/{}/Cp", material.get<std::string>() );
            auto Rho = specs[nl::json::json_pointer( matRho )].get<std::string>();
            auto Cp = specs[nl::json::json_pointer( matCp )].get<std::string>();

            lt += integrate( _range = markedelements( M_mesh, material.get<std::string>() ), 
                             //_expr = expr( Rho ) * expr( Cp ) * (((idv(M_bdf->polyDeriv()) + (cst(M_bdf->polyDerivCoefficient( 1 ))*( -idv(M_bdf->unknowns()[0]) + idv(T))) ))) * id( v ) );                   
                             //_expr = expr( Rho ) * expr( Cp ) *  (cst(M_bdf->polyDerivCoefficient( 1 ))*idv(T))* id( v ) );            
                             _expr = expr( Rho ) * expr( Cp ) * idv(M_bdf->polyDeriv())  * id( v ) );                          
        }

        // Blackbody radiative condition
        // Integrating sigma*epsilon*T_{n-1}^3*idt(u)*id(v)
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_blackbody_heat_flux" ) )
        {
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_blackbody_heat_flux"_json_pointer].items() )
            {                
                auto sigma = expr(value["sigma"].get<std::string>());

                sigma.setParameterValues({{"sigma",specs["/Parameters/sigma"_json_pointer].get<double>()}});

                // Recover the emissivity epsilon of the coating material associated to the marked face
                for ( std::string mark :  M_markers_map[bc] )
                {
                    auto done = 0;
                    for ( auto& [key, coat] : specs["/Coating"_json_pointer].items() )
                    {
                        for ( auto markcoat : coat.at("markers") )
                        {
                            if ( mark == markcoat )
                            {
                                auto epsilon = coat["epsilon"].get<std::string>();                            
                                auto idvu3 = idv( T ) * idv( T ) * idv( T );

                                at += integrate( _range = markedfaces( M_mesh, mark ),
                                        _expr =  sigma * expr( epsilon ) * idvu3 * idt( u ) * id( v ));
                                done = 1;
                                break;
                            }
                        }
                        if ( done ){ break; }
                    }    
                }
            }
        }

        // Radiative enclosure
        // Integrating directly the radiative flux obtained via the solution of the problem coming from Modest book
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
        {
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
            {                
                auto sigma = expr(value["sigma"].get<std::string>());

                sigma.setParameterValues({{"sigma",specs["/Parameters/sigma"_json_pointer].get<double>()}});
                
                for ( std::string mark :  M_markers_map[bc] )
                {
                    auto done = 0;
                    for ( auto& [key, coat] : specs["/Coating"_json_pointer].items() )
                    {
                        for ( auto markcoat : coat.at("markers") )
                        {
                            if ( mark == markcoat )
                            {                                
#if 1                                
                                lt += integrate( _range = markedfaces( M_mesh, mark ),
                                _expr = -(idv(q)) * id( v ) );                                       
#endif
                                done = 1;
                                break;
                            }
                        }
                        if ( done ){ break; }
                    }   
                }
            }
        }
        for ( auto& [bc, value] : specs["/BoundaryConditions/heat/temperature"_json_pointer].items())
        {
            auto dirichletBc = expr(value["expr"].get<std::string>());
            at+=on( _range=markedfaces(M_mesh,bc), _rhs=lt, _element=*T_sol, _expr=dirichletBc );
        }
        // Solution of the heat transfer problem
        at.solve( _rhs = lt, _solution = T_sol,_name="heatEq",_rebuild=true );

        M_currentTempAndFlux.setT(T_sol);

    } // end RHT<Dim,Order>::solveHeatEquation

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
                    fmt::format("Loading view factors for {} from {}",bc,vf_filename);
                    loadVF(bc,vf_filename);
                }
                else if(vf_status=="compute")
                {
                    auto vf_filename = value["viewfactors"]["filename"];
                    fmt::format("Computing view factors for {} from {}",bc,vf_filename);                    
                    computeVF(bc,vf_filename);
                }
                else
                {
                    std::cout << "View factor status not correct" << std::endl;
                }
            }
        }
    } // end RHT<Dim,Order>::computeVF_and_save


// Initialization function: charging mesh, defining fem spaces, view factors computation
    template<int Dim, int Order>
    void RHT<Dim,Order>::init()
    {
        // Initialize mesh, spaces, marker list and view factor matrix
        M_mesh = loadMesh( _mesh = new mesh_t, _filename = this->specs["/Meshes/heat/Import/filename"_json_pointer].template get<std::string>() );
        M_Xh = space_t::New(_mesh=M_mesh);             
        M_Xhvec = spacevect_t::New(_mesh=M_mesh);  
        tic();
        this->computeVF_and_save();        
        toc("Computation of view factors");  
              
        int n_cavities = M_markers_map.size();
        BlocksBaseVector<double> myblock_D(n_cavities);
        BlocksBaseVector<double> myblock_N(n_cavities);
        BlocksBaseSparseMatrix<double> myblock_M(n_cavities,n_cavities);

        int i =0;
        for(const auto& [cavity_name,markers]: M_markers_map)
        {            
            auto cavity_submesh = createSubmesh(_mesh=M_mesh,_range=markedfaces(M_mesh,markers),_update=0);
            M_surface_submesh_map.insert(std::make_pair(cavity_name,cavity_submesh));
            if(M_mesh->isParentMeshOf(cavity_submesh))
                LOG(INFO) << fmt::format("M_mesh is parent mesh of {} submesh",cavity_name)<< std::endl;
            
            // Create a discontinuous space, linear per face, over the view-factor markers
            M_Xhd0 = Pdh<Order-1>(cavity_submesh,true);
            M_Xhd1 = Pdh<Order>(cavity_submesh,true);
            M_Xhd0_map.insert(std::make_pair(cavity_name,M_Xhd0));
            M_Xhd1_map.insert(std::make_pair(cavity_name,M_Xhd1));

            M_D = backend()->newVector(M_Xhd0);
            M_N = backend()->newVector(M_Xhd0); 
            auto mat_graph = stencil( _test=M_Xhd0,_trial=M_Xhd0)->graph();
            M_M = backend()->newMatrix(0,0,0,0,mat_graph);  
            myblock_D(i,0)=M_D;
            myblock_N(i,0)=M_N; 
            myblock_M(i,i)=M_M;
            
            i++;
        }       
        M_M_block = backend()->newBlockMatrix(_block=myblock_M);
        M_N_block = backend()->newBlockVector(_block=myblock_N);
        M_D_block = backend()->newBlockVector(_block=myblock_D);

        // Create the bdf structure
        M_bdf = bdf( _space = M_Xh );
        
        auto T = M_Xh->elementPtr();
        auto q = M_Xh->elementPtr();

        // Initialize temperature for the heat equation
        auto T_init=expr(specs["/InitialConditions/heat/temperature/Expression/Tini/expr"_json_pointer].get<std::string>());
        T_init.setParameterValues({{"Tinit_C",specs["/Parameters/Tinit_C"_json_pointer].get<double>()}});
        for(auto & [bc,value]: specs["/InitialConditions/heat/temperature/Expression"_json_pointer].items())
        {            
            for(std::string mark : value.at("markers") )
            {
                T->on(_range=markedelements(M_mesh,mark),_expr=T_init);
            }
        }
        if(!boption(_name="deactivate-exporters"))
        {
            auto e_opInterp = exporter(_mesh=M_mesh,_name="init");
            e_opInterp->addRegions();
            e_opInterp->add("T",idv(T));        
            e_opInterp->save();
        }
        // Initialize the temperature and flux in the data structure needed for the Picard loop
        M_currentTempAndFlux.setT(T);
        M_currentTempAndFlux.setq(q);      

        // Initialize matrix and right-hand sides of the heat transfer problem
        M_lt = backend()->newVector(M_Xh);
        M_l = backend()->newVector(M_Xh); 
        auto cg_graph = stencil( _test=M_Xh,_trial=M_Xh)->graph();
        M_a = backend()->newMatrix(0,0,0,0,cg_graph);  
        M_at = backend()->newMatrix(0,0,0,0,cg_graph);  

        LOG(INFO) << fmt::format("Init routine finished")<< std::endl;
    
    }  // end RHT<Dim,Order>::init

// Build the M matrix for radiative heat transfer with several
// emitting surfaces
    template<int Dim, int Order>
    void RHT<Dim,Order>::build_M()
    {       
        auto v = M_Xhd0->element();
        auto u = M_Xhd0->element();      

        // Insert terms of the form 1/epsilon on the diagonal of the emitting surfaces
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
        {
            int i_cavity = 0;
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
            {
                for ( std::string mark :  M_markers_map[bc] )
                {
                    auto done = 0;
                    // Loop over the "Coating" section where epsilon is stored
                    for ( auto& [key, coat] : specs["/Coating"_json_pointer].items() )
                    {                        
                        for ( auto markcoat : coat.at("markers") )
                        {
                            if ( mark == markcoat )
                            {                    
                                auto epsilon = coat["epsilon"].get<std::string>();
                                
                                form2(_trial=M_Xhd0_map[bc],_test=M_Xhd0_map[bc],_matrix=M_M_block,
                                        _rowstart=i_cavity, _colstart=i_cavity) += 
                                            integrate( _range= markedfaces(M_mesh,mark),
                                                 _expr = cst(1.)/expr(epsilon)*inner(idt(u),id(v)) );       

                                done = 1;
                                break;
                            }
                        }
                        if ( done ){ break; }
                    }
                }
                i_cavity++;
            }            
        }
    } // RHT<Dim,Order>::build_M

// Rebuild the N vector for radiative heat transfer with several
// emitting surfaces
    template<int Dim, int Order>
    void RHT<Dim,Order>::rebuild_N(  element_ptr_t flux )
    {
        M_N_block->zero();   
        int i_cavity=0;
        // Insert terms of the form  \int_bdry (1/epsilon-1)*q(x')dF(x,x') in the vector N
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
        {
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
            {   
                auto v = M_Xhd0_map[bc]->element();
                auto u = M_Xhd0_map[bc]->element();    
                
                auto flux_interp = M_Xhd0_map[bc]->element();
                auto flux_interp_interm = M_Xhd1_map[bc]->element();
                auto vf_field = M_Xhd0_map[bc]->element();                
                // Interpolate the continous flux onto the discontinuous space
                auto opI_contToDisc = opInterpolation(_domainSpace=M_Xh,
                                            _imageSpace=M_Xhd1_map[bc],
                                            _backend=backend(_name="interp1",_rebuild=true),
                                            _type=InterpolationConforme());
                auto opI_1to0 = opInterpolation(_domainSpace=M_Xhd1_map[bc],
                                            _imageSpace=M_Xhd0_map[bc],
                                            _backend=backend(_name="interp2",_rebuild=true),
                                            _type=InterpolationConforme());     

                opI_contToDisc->apply(*flux,flux_interp_interm);                                                       
                opI_1to0->apply( flux_interp_interm, flux_interp );
                int i_mark=0; 
                int j_mark=0;        
                auto cavity_markers = M_markers_map[bc];   
                for ( std::string mark :  cavity_markers ) // Loop on index i
                {
                    //std::cout << "For marker " << mark << std::endl;                    
                    for ( std::string mark2 :  cavity_markers ) // Loop on index j
                    {
                        // Compute the index relative to marker mark
                        auto it = std::find(cavity_markers.begin(),cavity_markers.end(),mark);                    
                        i_mark = it - cavity_markers.begin();
                        // Compute the index relative to marker mark2
                        auto jt = std::find(cavity_markers.begin(),cavity_markers.end(),mark2);                    
                        j_mark = jt - cavity_markers.begin();
                        //std::cout << fmt::format("we add a contribution from index {} whose view factor is {}",j_mark,M_vf_matrix(i_mark,j_mark))<< std::endl;
                        auto done = 0;
                        for ( auto& [key, coat] : specs["/Coating"_json_pointer].items() )
                        {                        
                            for ( auto markcoat : coat.at("markers") )
                            {
                                vf_field.zero();
                                if ( mark2 == markcoat )
                                {                    
                                    auto epsilon_mark2 = coat["epsilon"].get<std::string>();
                                    //std::cout << fmt::format("we add a contribution from marker {} whose emissivity is {}",markcoat,epsilon_mark2)<< std::endl;
                                    vf_field.on(_range=markedfaces(M_mesh,mark2),_expr=cst(M_matrix_vf_map[bc](i_mark,j_mark)));
                                    
                                    auto q_scal_vf = integrate(_range=markedfaces(M_mesh,mark2), _expr=inner(idv(flux),idv(vf_field))).evaluate()(0,0);
                                    auto measure_mark2 = integrate( _range=markedfaces(M_mesh,mark2),_expr= cst(1.) ).evaluate()(0,0);
                                    
                                    form1(_test=M_Xhd0_map[bc],_vector=M_N_block,
                                        _rowstart=i_cavity) += integrate( _range=markedfaces(M_mesh,mark),
                                        _expr = (cst(1.)/expr(epsilon_mark2)-cst(1.))/cst(measure_mark2) * cst(q_scal_vf) * id(v) );

                                    done = 1;


                                }
                            }
                            if ( done ){ break; }
                        }  
                    }                  
                }
                i_cavity++;
            }
        }
    } // end RHT<Dim,Order>::rebuild_N


// Rebuild the D right-hand side for radiative heat transfer with several
// emitting surfaces
    template<int Dim, int Order>
    void RHT<Dim,Order>::rebuild_D(  element_ptr_t T )
    {           
        // M_D->zero();
        M_D_block->zero();  
        int i_cavity=0;      
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
        {
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
            {      
                auto v = M_Xhd0_map[bc]->element();
                auto u = M_Xhd0_map[bc]->element();    
                
                auto T4_interp = M_Xhd0_map[bc]->element();
                auto T4_interp_interm = M_Xhd1_map[bc]->element();
                auto vf_field = M_Xhd0_map[bc]->element();

                // Interpolate the continous temperature onto the discontinuous space    
                auto opI_interm = opInterpolation(_domainSpace=M_Xh,
                                            _imageSpace=M_Xhd1_map[bc],                                    
                                            _backend=backend(_name="interp1",_rebuild=true),
                                            _type=InterpolationConforme());
                auto opI_1to0 = opInterpolation(_domainSpace=M_Xhd1_map[bc],
                                            _imageSpace=M_Xhd0_map[bc],                                    
                                            _backend=backend(_name="interp2",_rebuild=true),
                                            _type=InterpolationConforme());                                    
                auto field = M_Xh->element();
                field.on(_range=elements(M_mesh),_expr=idv(T)*idv(T)*idv(T)*idv(T));
                opI_interm->apply(field,T4_interp_interm);
                opI_1to0->apply( T4_interp_interm,T4_interp);

                //  Loop over the emissivities and build the D right-hand side
                // Integrate term of the form sigma*T^4 - \int_bdry (sigma*T^4(x') dF(x,x'))
                int i_mark=0; 
                int j_mark=0;           
                auto sigma = expr(value["sigma"].get<std::string>());

                sigma.setParameterValues({{"sigma",specs["/Parameters/sigma"_json_pointer].get<double>()}});
                
                auto cavity_markers = M_markers_map[bc];   

                int i_marker_to_bb=0;
                for ( std::string mark :  cavity_markers )
                {
                    //std::cout << "For marker " << mark << std::endl;
                    for ( std::string mark2 : cavity_markers ) // Loop on index j
                    {
                        // Compute the index relative to marker mark
                        auto it = std::find(cavity_markers.begin(),cavity_markers.end(),mark);                    
                        i_mark = it - cavity_markers.begin();
                        // Compute the index relative to marker mark2
                        auto jt = std::find(cavity_markers.begin(),cavity_markers.end(),mark2);                    
                        j_mark = jt - cavity_markers.begin();
                        //std::cout << fmt::format("we add a contribution from index {} whose view factor is {}",j_mark,M_vf_matrix(i_mark,j_mark))<< std::endl;
                        
                        vf_field.zero();                                                       
                        //std::cout << fmt::format("we add a contribution from marker {} whose emissivity is {}",markcoat,epsilon_mark2)<< std::endl;
                        vf_field.on(_range=markedfaces(M_mesh,mark2),_expr=cst(M_matrix_vf_map[bc](i_mark,j_mark)));                      
                    
        
                        auto T4_scal_vf = integrate(_range=markedfaces(M_mesh,mark2), _expr=inner(idv(T4_interp),idv(vf_field))).evaluate()(0,0);                                                        
                        
                        auto measure_mark2 = integrate( _range=markedfaces(M_mesh,mark2),_expr= cst(1.) ).evaluate()(0,0);
                        
                        form1(_test=M_Xhd0_map[bc],_vector=M_D_block,
                                        _rowstart=i_cavity) += integrate( _range=markedfaces(M_mesh,mark),
                                        _expr= cst(-1)/cst(measure_mark2)*sigma*cst(T4_scal_vf) *id(v) );  
                                            
                    }

                    // If the cavity is open, and it is assumed that a black body of fixed temperature
                    // T_ref exchances heat with the cavity, and additional term is added: 
                    // - the view factor is computed using the reciprocity formula FijAi = FjiAj
                    // - the contribution is of the form \sigma T^4 * Fij
                    if(value["enclosure"]=="open")
                    {
                        //std::cout << "Open enclosure " << bc << std::endl;
                        double vf_marker_to_bb = 1 - M_matrix_vf_map[bc].row(i_marker_to_bb).sum(); //Fij
                        double vf_bb_to_marker = 1 - M_matrix_vf_map[bc].col(i_marker_to_bb).sum(); //Fji
                        auto measure_mark1 = integrate( _range=markedfaces(M_mesh,mark),_expr= cst(1.) ).evaluate()(0,0); //Ai
                        double area_bb = measure_mark1 * vf_marker_to_bb / vf_bb_to_marker;
                        auto T_bb = value["Tref"].get<double>();
                        auto T_bb4 = cst(T_bb)*cst(T_bb)*cst(T_bb)*cst(T_bb);
                        form1(_test=M_Xhd0_map[bc],_vector=M_D_block,
                                        _rowstart=i_cavity) += integrate( _range=markedfaces(M_mesh,mark),
                                        _expr= cst(-1)*sigma*T_bb4*cst(vf_bb_to_marker) *id(v) ); 
                    }
                    
                    form1(_test=M_Xhd0_map[bc],_vector=M_D_block,
                                        _rowstart=i_cavity) += integrate( _range=markedfaces(M_mesh,mark),
                                        _expr= inner(sigma*idv(T4_interp), id(v) ));
                    i_marker_to_bb++;
                }
                i_cavity++;
                if(!boption(_name="deactivate-exporters"))
                {
                    auto e_opInterp = exporter(_mesh=M_mesh,_name="functD");
                    e_opInterp->addRegions();
                    e_opInterp->add("T4interp"+bc,idv(T4_interp)); 
                    e_opInterp->save();
                }
            }
        }        
         
    } // end RHT<Dim,Order>::rebuild_D

// Solve the double Picard loop for temperature and radiative flux convergence
// by solving the heat equation and the multiple-surface equation in Modest's book
    template<int Dim, int Order>
    void RHT<Dim,Order>::solvePicardIteration()
    {        
        // Initialize the tolerances, norm values and initial temperatures and fluxes
        double tol_T = 1e-2;
        double tol_q = 1e-2;
        double norm_T_Tnew = 1;
        double norm_q_qnew = 1;
        auto T = M_Xh->elementPtr(); 
        auto q = M_Xh->elementPtr(); 
        *T = *M_bdf->unknowns()[0]; // at the beginning of the loop, the temperature it the one at the previous time step
        *q = *M_currentTempAndFlux.q();// at the beginning of the loop, the flux it the one at the previous time step
        
        int n_cavities = M_markers_map.size();
        BlocksBaseVector<double> myblock_L(n_cavities);
        BlocksBaseVector<double> myblock_sol(n_cavities);
        BlocksBaseSparseMatrix<double> myblock_M(n_cavities,n_cavities);

        vector_ptrtype L;
        vector_ptrtype L_block,sol_block;

        int i =0;
        for(const auto& [cavity_name,Xhd0_space]: M_Xhd0_map)
        {
            L = backend()->newVector(Xhd0_space);           
            myblock_L(i,0)=L;           
            myblock_sol(i,0)=Xhd0_space->elementPtr();
            i++;
        }

        // Initial guess for the flux
        // Important for the correct initialization of the problem
        this->rebuild_N(M_currentTempAndFlux.q()); 
        this->rebuild_D(T);
                
        L_block=backend()->newBlockVector(_block=myblock_L);
        L_block->add(*M_N_block);
        L_block->add(*M_D_block);

        auto q_disc=M_Xhd0->element();
        auto q_disc_interm=M_Xhd1->element();

        sol_block = backend()->newBlockVector(_block=myblock_sol);
        backend(_rebuild=true)->solve( _matrix=M_M_block, _rhs=L_block, _solution=sol_block );
        myblock_sol.localize(sol_block);

        // Solve the Modest equation with q from previous iteration and new temperature
        
        int i_cavity=0;
        for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
        {
            auto q_interp = M_Xh->element();
            auto q_disc=M_Xhd0_map[bc]->element();
            auto q_disc_interm=M_Xhd1_map[bc]->element();
            // Interpolate the discontinuous solution onto the continuous space
            auto opI1 = opInterpolation(_domainSpace=M_Xhd0_map[bc],
                                    _imageSpace=M_Xhd1_map[bc],
                                    _type=InterpolationConforme());
            auto opI2 = opInterpolation(_domainSpace=M_Xhd1_map[bc],
                                    _imageSpace=M_Xh,
                                    _range=markedfaces(M_mesh,M_markers_map[bc]),
                                    _type=InterpolationConforme());  
            q_disc.container()= (*myblock_sol(i_cavity,0));            
            opI1->apply( q_disc, q_disc_interm );
            opI2->apply( q_disc_interm,q_interp );
            q->on(_range = markedfaces(M_mesh,M_markers_map[bc]), _expr=idv(q_interp));
            i_cavity++;
        }

        // Initial guess for temperature and flux 
        // to be inserted in the Picard loop
        // Results of this function are stored in the structure M_currentTempAndFlux
        this->solveHeatEquation(M_bdf->unknowns()[0],q);
        
        auto T_new=M_Xh->elementPtr();
        auto q_new_k_l=M_Xh->elementPtr();
        auto T_relaxed=M_Xh->elementPtr();
        // 
        *T_new = *M_currentTempAndFlux.T();
       
        auto norm_T_init = normL2(_range=elements(M_mesh),_expr = idv(T));


        auto e_opInterp = exporter(_mesh=M_mesh,_name="Tloop");    
        
        // Construct the emission term with view factors for the Modest equation
        this->rebuild_N(M_currentTempAndFlux.q());         

        std::cout  << "inside Picard iteration" << std::endl;
        double iter_T=0;
        if(boption(_name="deactivate-exporters"))
        {
            M_M_block->printMatlab("M_block.m");
            M_N_block->printMatlab("N_block.m");
            M_D_block->printMatlab("D_block.m");
        }
        // Since Modest and heat equations are coupled via temperature and radiative flux
        // a double Picard loop is proposed to have the convergence of the two quantities
        while(norm_T_Tnew>tol_T)
        {
            this->rebuild_D(T_new);
            
            L_block=backend()->newBlockVector(_block=myblock_L);
            L_block->add(*M_N_block);
            L_block->add(*M_D_block);

            sol_block = backend()->newBlockVector(_block=myblock_sol);
            backend(_rebuild=true)->solve( _matrix=M_M_block, _rhs=L_block, _solution=sol_block );
            myblock_sol.localize(sol_block);

            auto q_disc=M_Xhd0->element();
            auto q_disc_interm=M_Xhd1->element();            

            // Interpolate the discontinuous solution onto the continuous space
            int i_cavity=0;
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
            {       
                auto q_interp = M_Xh->element();         
                auto q_disc=M_Xhd0_map[bc]->element();
                auto q_disc_interm=M_Xhd1_map[bc]->element();
                auto opI1 = opInterpolation(_domainSpace=M_Xhd0_map[bc],
                                        _imageSpace=M_Xhd1_map[bc],
                                        _type=InterpolationConforme());
                auto opI2 = opInterpolation(_domainSpace=M_Xhd1_map[bc],
                                        _imageSpace=M_Xh,
                                        _range=markedfaces(M_mesh,M_markers_map[bc]),
                                        _type=InterpolationConforme());      
                q_disc.container() = *myblock_sol(i_cavity,0);
                opI1->apply( q_disc , q_disc_interm );
                opI2->apply( q_disc_interm, q_interp );
                q->on(_range = markedfaces(M_mesh,M_markers_map[bc]), _expr=idv(q_interp));
                i_cavity++;
            }
                    
            auto e_opInterp2 = exporter(_mesh=M_mesh,_name="insideLoop");
            double iter=0;
            auto norm_q_init = 0.0;
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
            {
                norm_q_init += normL2(_range=markedfaces(M_mesh,M_markers_map[bc]),
                                _expr = idv(q));
            }
            double norm_q_qnew = 1;
            // Loop over the radiative flux for a fixed temperature field
            while(norm_q_qnew > tol_q)
            {
                // Rebuild the N term of the Modest equation
                this->rebuild_N(q);                 

                L_block=backend()->newBlockVector(_block=myblock_L);
                L_block->add(*M_N_block);
                L_block->add(*M_D_block);
           
                auto q_new_k_l_disc=M_Xhd0->element();       
                auto q_new_k_l_disc_interm=M_Xhd1->element();       
                tic();                

                sol_block = backend()->newBlockVector(_block=myblock_sol);
                backend(_rebuild=true)->solve( _matrix=M_M_block, _rhs=L_block, _solution=sol_block );
                myblock_sol.localize(sol_block);
                toc("Solve radiative heat transfer on surfaces");
                
                int i_cavity=0;
                for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
                {         
                    auto q_interp = M_Xh->element();       
                    auto q_new_k_l_disc=M_Xhd0_map[bc]->element();       
                    auto q_new_k_l_disc_interm=M_Xhd1_map[bc]->element(); 
                    auto opI1 = opInterpolation(_domainSpace=M_Xhd0_map[bc],
                                            _imageSpace=M_Xhd1_map[bc],
                                            _type=InterpolationConforme());
                    auto opI2 = opInterpolation(_domainSpace=M_Xhd1_map[bc],
                                            _imageSpace=M_Xh,
                                            _range=markedfaces(M_mesh,M_markers_map[bc]),
                                            _type=InterpolationConforme());
                    q_new_k_l_disc.container() = (*myblock_sol(i_cavity,0));
                    opI1->apply( q_new_k_l_disc, q_new_k_l_disc_interm );                    
                    opI2->apply( q_new_k_l_disc_interm, q_interp );
                    q_new_k_l->on(_range = markedfaces(M_mesh,M_markers_map[bc]), _expr=idv(q_interp));
                    i_cavity++;
                }
                // Check if the difference between two iterates of the Picard loop is small enough
                norm_q_qnew=0.0;
                for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
                {   
                    norm_q_qnew += normL2(_range=markedfaces(M_mesh,M_markers_map[bc]),
                                _expr = (idv(q_new_k_l)-idv(q)));
                }
                norm_q_qnew /=norm_q_init;

                *q=*q_new_k_l;        
                std::cout << "norm_q_qnew" << norm_q_qnew << std::endl;    

                M_M_block->printMatlab("M_1_block.m");
                M_N_block->printMatlab("N_1_block.m");
                M_D_block->printMatlab("D_1_block.m");            
               
                iter++;
            }
            // Solve the heat trasfer equation with the newly computed radiative heat flux from Modest equation
            tic();
            this->solveHeatEquation(T_new,q);  
            toc("Solve heat PDE");
            *T = *T_new;

            *T_new = *M_currentTempAndFlux.T();

            // Check if two iterates of the solution of the heat transfer equations are close enough
            norm_T_Tnew = normL2(_range=elements(M_mesh),
                            _expr = (idv(T_new)-idv(T)));

            norm_T_Tnew /= norm_T_init;

            std::cout << "norm_T_Tnew" << norm_T_Tnew << std::endl;
            if(!boption(_name="deactivate-exporters"))
            {
                e_opInterp->step(iter_T)->add("oldT",idv(T));
                e_opInterp->step(iter_T)->add("newT",idv(T_new));       
                e_opInterp->save();
            }
            iter_T++; 
        }
        M_currentTempAndFlux.setq(q);
        // The Picard loop has finished; save the temperature as T_{n} to proceed to next iteration
        std::cout << "end Picard iteration" << std::endl;
        auto newT = unwrap_ptr(M_currentTempAndFlux.T());
        M_bdf->next(newT);        

    } // end RHT<Dim,Order>::solvePicardIteration()


    // Routine solving the heat transfer problem
    template<int Dim, int Order>
    void RHT<Dim,Order>::execute()
    {        
        // Initialize spaces, mesh, matrices
        this->init();
        
        // Initialize the terms of the heat transfer problem that are time-independent
        this->initHeatEquation();

        // Compute the diagonal matrix to solve the problem with multiple radiating surfaces
        this->build_M();
        //this->build_M_block();

        // Solve the time dependent heat transfer problem
        for( ;!M_bdf->isFinished(); )
        {
            std::cout << fmt::format("It's time {} of the resolution of the heat equation",M_bdf->time()) <<std::endl;

            // Solve the heat transfer problem at each time instant by solving a Picard loop to ensure convergence
            // of temperature and fluxes
            tic();
            this->solvePicardIteration();
            toc("solve Picard iteration");
            // Export the temperature, and flux on the radiating surfaces
            if(!boption(_name="deactivate-exporters"))
                this->exportHeat();
        }
        this->checkResults();
    } // end RHT<Dim,Order>::execute()

    // Routine comparing the results with literature
    template<int Dim, int Order>
    void RHT<Dim,Order>::checkResults()
    {
        for(auto [key,struc] : specs["/Checker"_json_pointer].items())
        {
            std::cout << "key" << key << std::endl;
            if(struc["type"].get<std::string>()=="average")
            {
                auto markers = struc["markers"].get<std::vector<std::string>>();
                auto value = struc["exact_value"].get<double>();
                auto tol = struc["rel_tolerance"].get<double>();
                
                if(struc["quantity"]=="flux")
                {
                    // average normal radiative flux
                    auto average_q = integrate(_range=markedfaces(M_mesh,markers),_expr=idv(M_currentTempAndFlux.q())).evaluate()(0,0);
                    auto measure_markers = integrate(_range=markedfaces(M_mesh,markers),_expr=cst(1.)).evaluate()(0,0);
                    average_q /=measure_markers;

                    auto rel_difference = (average_q-value)/value;

                    assert(math::abs(rel_difference) <= tol);
                    std::cout << "average_q" << average_q << std::endl;
                    std::cout << "markers" << markers << std::endl;
                    std::cout << "value" << value << std::endl;
                    std::cout << "tol" << tol << std::endl;
                    std::cout << "rel_difference" << rel_difference << std::endl;
                }
                else if (struc["quantity"]=="temperature")
                {
                    // average temperature on surface
                    auto average_T = integrate(_range=markedfaces(M_mesh,markers),_expr=idv(M_currentTempAndFlux.T())).evaluate()(0,0);
                    auto measure_markers = integrate(_range=markedfaces(M_mesh,markers),_expr=cst(1.)).evaluate()(0,0);
                    average_T /=measure_markers;

                    auto rel_difference = (average_T-value)/value;

                    assert(math::abs(rel_difference) <= tol);
                    std::cout << "average_T" << average_T << std::endl;
                    std::cout << "markers" << markers << std::endl;
                    std::cout << "value" << value << std::endl;
                    std::cout << "tol" << tol << std::endl;
                    std::cout << "rel_difference" << rel_difference << std::endl;

                }
                else
                {
                    std::cout << fmt::format("Quantity average checker not supported ");
                }
            }
            else
            {
                std::cout << fmt::format("Non average checker not supported ");
            }
        }
        // Case in which an average value is available        
        // Case in which a profile is compared
    }


} // namespace Feel



int main( int argc, char** argv )
{
    using namespace Feel;
    try
    {
        Environment env( _argc = argc, _argv = argv,
                         _desc = makeOptions(),
                         _about = about( _name = fmt::format( "rht-{}dp{}", FEELPP_DIM, FEELPP_ORDER ),
                                         _author = "Feel++ Consortium",
                                         _email = "feelpp@cemosis.fr" ) );

        // Read the json file associated to the heat transfer problem
        auto jsonfile = removeComments( readFromFile( Environment::expand( soption( "specs" ) ) ) );
        std::istringstream istr( jsonfile );
        json specs = json::parse( istr );    

        // Read the json file associated to viewfactor computation ( it must contain only the markers for which view 
        // factors are computed )
        // jsonfile = removeComments( readFromFile( Environment::expand( soption( "viewfactor" ) ) ) );
        // std::istringstream astr( jsonfile );
        // json viewfactor = json::parse( astr );   

        // Instantiate the class for the solution of the heat transfer problem
        RHT<FEELPP_DIM, FEELPP_ORDER> rht(specs);

        // Solve the heat transfer problem
        rht.execute();

    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
} // end main()