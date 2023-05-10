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
#include <math.h>
// #include <feel/feeldiscr/context.hpp>

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
        void solveHeatEquationNonLinear(element_ptr_t T );

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
        void executeNonLinear();
        void computeVF_and_save();
        void saveVF(std::string cavity_name,const Eigen::Ref<const Eigen::MatrixXd>&  M);
        void loadVF(std::string cavity_name,std::string filename);
        void computeVF(std::string cavity_name,std::string filename);
        void checkResults();
        void computeMatrixD(mesh_trace_ptr_t cavity_submesh, Eigen::MatrixXd& matrixD);
        void addCavityGraph(std::shared_ptr<GraphCSR> modified_graph_with_cavities,mesh_trace_ptr_t cavity_submesh);

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
        spacedisc_ptr_t M_Xhd0;
        spacedisc_interm_ptr_t M_Xhd1;
        std::map< std::string, spacedisc_ptr_t > M_Xhd0_map;
        std::map< std::string, spacedisc_interm_ptr_t > M_Xhd1_map;
        mesh_ptr_t M_mesh;
        std::map< std::string, mesh_trace_ptr_t> M_surface_submesh_map;
        std::map< std::string, Eigen::VectorXd> M_surface_submesh_local_areas;
        std::map< std::string, std::vector<std::string> > M_markers_map;   
        bdf_ptrtype M_bdf;
        TandQ M_currentTempAndFlux;        

        sparse_matrix_ptrtype M_M, M_M_block; // matrix for problem on multiple irradiating surfaces
        std::map< std::string, Eigen::MatrixXd > M_matrix_vf_map; // matrix storing view factors
        vector_ptrtype M_D, M_N, M_D_block,M_N_block; // vectors for problem on multiple irradiating surfaces
        std::map< std::string, Eigen::MatrixXd > M_cavity_coupling_mat;

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
        // if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "convective_heat_flux" ) )
        // {
        //     for ( auto& [bc, value] : specs["/BoundaryConditions/heat/convective_heat_flux"_json_pointer].items() )
        //     {
        //         LOG( INFO ) << fmt::format( "convective_heat_flux {}: {}", bc, value.dump() );
        //         auto h = value["h"].get<std::string>();
        //         auto Text = expr(value["Text"].get<std::string>());

        //         Text.setParameterValues({{"Tref_C",specs["/Parameters/Tref_C"_json_pointer].get<double>()}});

        //         a += integrate( _range = markedfaces( M_mesh, bc ),
        //                 _expr = expr( h ) * id( v ) * idt( u ) );
        //         l += integrate( _range = markedfaces( M_mesh, bc ),
        //                 _expr =  expr( h ) * Text * id( v ) );
        //     }
        // }

        // // ===== BC RHT =====
        // // Blackbody radiative condition (no view factors)
        // // set the term sigma * epsilon * T_0^4
        // if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_blackbody_heat_flux" ) )
        // {
        //     for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_blackbody_heat_flux"_json_pointer].items() )
        //     {
        //         LOG( INFO ) << fmt::format( "radiative_blackbody_heat_flux {}: {}", bc, value.dump() );
                            
        //         auto sigma = expr(value["sigma"].get<std::string>());                
        //         auto Tref = expr(value["Tref"].get<std::string>());
        //         Tref.setParameterValues({{"Tref_C",specs["/Parameters/Tref_C"_json_pointer].get<double>()}});
        //         sigma.setParameterValues({{"sigma",specs["/Parameters/sigma"_json_pointer].get<double>()}});

        //         auto Tref4 = Tref*Tref*Tref*Tref;

        //         // Recover the emissivity epsilon of the coating material associated to the emitting face
        //         for ( std::string mark :  value.at("markers") )
        //         {
        //             auto done = 0;
        //             for ( auto& [key, coat] : specs["/Coating"_json_pointer].items() )
        //             {
        //                 for ( auto markcoat : coat.at("markers") )
        //                 {                        
        //                     if ( mark == markcoat )
        //                     {
        //                         std::cout << mark << " "<<markcoat <<std::endl;
        //                         auto epsilon = coat["epsilon"].get<std::string>();

        //                         l += integrate( _range = markedfaces( M_mesh, mark ), 
        //                                 _expr =  sigma * expr( epsilon ) * Tref4 * id( v ) );
        //                         done = 1;
        //                         break;
        //                     }
        //                 }
        //                 if ( done ){ break; }
        //             }
        //         }
        //     }
        // }   
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
    
    // Compute inverse matrices for the Jacobian of the radiative condition
    template<int Dim, int Order>
    void RHT<Dim,Order>::computeMatrixD(mesh_trace_ptr_t mesh, Eigen::MatrixXd& matrixD )
    {
        auto size_mat = nelements(elements(mesh));
        Eigen::MatrixXd detailedVFmatrix(size_mat,size_mat);
        matrixD.resize(size_mat,size_mat);
        Eigen::MatrixXd matrixG(size_mat,size_mat), inverseMatrixG(size_mat,size_mat);
        Eigen::VectorXd areas(size_mat), emissivity(size_mat);
        detailedVFmatrix.setZero();
        inverseMatrixG.setZero(); 
        areas.setZero(); 
        emissivity.setZero();
        matrixD.setZero();
        std::map<int,double > markerIntToEmissivity;

        // Create a map markerInt -> emissivity of the marker
        for ( auto& [key, coat] : specs["/Coating"_json_pointer].items() )
        {                        
            for ( std::string markcoat : coat.at("markers") )
            {
                markerIntToEmissivity.insert( std::make_pair(mesh->markerName(markcoat), std::stod( coat["epsilon"].get<std::string>() )) ) ;

                if(std::stod( coat["epsilon"].get<std::string>()) == 1 )
                {
                    std::cout << "There are black facets in the cavity. Abort" << std::endl;
                    return;
                }

            }
        }
        std::cout << markerIntToEmissivity << std::endl;
        for ( auto& [bc, value] : this->specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
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

                auto jsonfile = removeComments( readFromFile( Environment::expand( vf_filename ) ) );
                std::istringstream astr( jsonfile );
                json json_vf = json::parse( astr );    


                if(json_vf["viewfactor"]["type"]=="UnobstructedPlanar")
                {
                    UnobstructedPlanarViewFactor<mesh_t> upvf( this->M_mesh, json_vf , true);              
                    upvf.compute();    
                    detailedVFmatrix = upvf.viewFactors();          
                    areas = upvf.areas();    
                    M_surface_submesh_local_areas.insert(std::make_pair(bc,areas));
                }
            }
            else
            {
                std::cout << "View factor status not correct" << std::endl;
            }
        }             
        matrixG = detailedVFmatrix.array().rowwise() * areas.transpose().array();
        matrixG *= -1;
        // for(int i=0; i < size_mat ; i++)
        int i=0;
        for(auto elem : elements(mesh))
        {
            auto const& elt = boost::unwrap_ref( elem );
            matrixG(i,i) += areas[i] / ( 1 - markerIntToEmissivity[ elt.marker().value() ]);
            i++;
        }
        inverseMatrixG = matrixG.inverse();        


        // Compute the matrix Dab containing the facet-facet radiative interaction terms
        i=0;
        for(auto elem_a : elements(mesh))
        {
            auto const& elt_a = boost::unwrap_ref( elem_a );
            int j=0;
            // if facet a is black-body like -> emissivity = 1
            auto emissivity_a = markerIntToEmissivity[ elt_a.marker().value()];
            if( emissivity_a == 1 )
            {
                std::cout << "Case of a non gray is not treated" << std::endl;
                for(auto elem_b : elements(mesh))
                {
                    auto const& elt_b = boost::unwrap_ref( elem_b );

                    auto emissivity_b = markerIntToEmissivity[ elt_b.marker().value()];
                    // if facet b is black-body like -> emissivity = 1
                    if( emissivity_b == 1 )
                    {
                        
                    }
                    else // if facet b is gray-body like -> emissivity < 1
                    {

                    }
                    j++;
                }
            }
            else // if facet a is gray-body like -> emissivity < 1
            {
                for(auto elem_b : elements(mesh))
                {
                    auto const& elt_b = boost::unwrap_ref( elem_b );

                    auto emissivity_b = markerIntToEmissivity[ elt_b.marker().value()];
                    // if facet b is black-body like -> emissivity = 1
                    if( emissivity_b == 1 )
                    {
                        std::cout << "Case of b non gray is not treated" << std::endl;
                    }
                    else // if facet b is gray-body like -> emissivity < 1
                    {
                        matrixD(i,j) = inverseMatrixG(i,j) * emissivity_a / (1-emissivity_a) * areas[i] * emissivity_b / (1-emissivity_b) * areas[j];                        
                    }

                    j++;
                }
            }
            i++;
        }       

    }

// Add the appropriate lines in the graph of the Jacobian matrix of the radiative heat boundary conditions
template<int Dim, int Order>
void RHT<Dim,Order>::addCavityGraph(std::shared_ptr<GraphCSR> modified_graph_with_cavities, mesh_trace_ptr_t cavity_submesh)
{
    // For each DOF in the cavity, one must insert all the dofs of the cavity as potential 
    // non-zero entries. For this reason, first a loop is performed to collect all 
    // global numbering of the cavity dofs, and then all of these are inserted on the corresponding
    // lines in the graph

    // Collect all cavity dofs 
    std::vector<GraphCSR::size_type> cavity_dofs;
    for(auto elem: elements(cavity_submesh))
    {
        auto const& elt = boost::unwrap_ref( elem );

        const int dof_s = M_Xh->dof()->getIndicesSize();
        std::vector<GraphCSR::size_type> element_dof_local(dof_s);
        for( auto const& ldof : M_Xh->dof()->faceLocalDof( elt.id() ) )
        {
            cavity_dofs.push_back(ldof.index() );
        }        
    }    

    // Delete doubles by transforming into set
    std::set<int> s( cavity_dofs.begin(), cavity_dofs.end() );
    cavity_dofs.assign( s.begin(), s.end() );

    // Create a line in the matrix graph with non-zero entries where cavity dofs are 
    for(auto il1 : cavity_dofs)
    {

        const GraphCSR::size_type ig1 = M_Xh->dof()->mapGlobalProcessToGlobalCluster()[il1];
        auto theproc = M_Xh->dof()->procOnGlobalCluster( ig1 );

        {
            GraphCSR::row_type& row = modified_graph_with_cavities->row( ig1 );
            row.get<0>() = theproc;
            row.get<1>() = il1;
            row.get<2>().insert( cavity_dofs.begin(), cavity_dofs.end() );
        }
    
    }
}



// Initialization function: charging mesh, defining fem spaces, view factors computation
    template<int Dim, int Order>
    void RHT<Dim,Order>::init()
    {
        // Initialize mesh, spaces, marker list and view factor matrix
        M_mesh = loadMesh( _mesh = new mesh_t, _filename = this->specs["/Meshes/heat/Import/filename"_json_pointer].template get<std::string>() );
        M_Xh = space_t::New(_mesh=M_mesh); 
        tic();
        this->computeVF_and_save();        
        toc("Computation of view factors");  
              
        int n_cavities = M_markers_map.size();        

        int i =0;

        auto modified_graph_with_cavities = stencil( _test=M_Xh,_trial=M_Xh,_close=false)->graph();
        
        
        for(const auto& [cavity_name,markers]: M_markers_map)
        {            
            auto cavity_submesh = createSubmesh(_mesh=M_mesh,_range=markedfaces(M_mesh,markers),_update=0,_view=1);
            M_surface_submesh_map.insert(std::make_pair(cavity_name,cavity_submesh));
            if(M_mesh->isParentMeshOf(cavity_submesh))
                LOG(INFO) << fmt::format("M_mesh is parent mesh of {} submesh",cavity_name)<< std::endl;

            M_Xhd0 = Pdh<Order-1>(cavity_submesh,true);
            M_Xhd0_map.insert(std::make_pair(cavity_name,M_Xhd0));

            M_Xhd1 = Pdh<Order>(cavity_submesh,true);
            M_Xhd1_map.insert(std::make_pair(cavity_name,M_Xhd1));

            // Compute the matrix coupling all the dofs in the cavity
            Eigen::MatrixXd cavityCouplingMatrix;            
            computeMatrixD(cavity_submesh,cavityCouplingMatrix);            
            M_cavity_coupling_mat.insert(std::make_pair(cavity_name,cavityCouplingMatrix));

            // Add the appropriate lines in the matrix graph, accomodating the presence of
            // radiative boundary conditions for the Jacobian of the radiative condition
            addCavityGraph(modified_graph_with_cavities,cavity_submesh);
        }       
        modified_graph_with_cavities->close();
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
        
        
        M_a = backend()->newMatrix(0,0,0,0,modified_graph_with_cavities);  
        M_at = backend()->newMatrix(0,0,0,0,modified_graph_with_cavities);  

        M_a->printMatlab("M_a.m");

        LOG(INFO) << fmt::format("Init routine finished")<< std::endl;
    
    }  // end RHT<Dim,Order>::init

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

     // Routine solving the heat transfer problem
    template<int Dim, int Order>
    void RHT<Dim,Order>::executeNonLinear()
    {        
        // Initialize spaces, mesh, matrices
        this->init();
        
        // Initialize the terms of the heat transfer problem that are time-independent
        this->initHeatEquation();

        // Solve the time dependent heat transfer problem
        for( ;!M_bdf->isFinished(); )
        {
            std::cout << fmt::format("It's time {} of the resolution of the heat equation",M_bdf->time()) <<std::endl;

            // Solve the heat transfer problem at each time instant by solving a Picard loop to ensure convergence
            // of temperature and fluxes
            tic();
            this->solveHeatEquationNonLinear(M_bdf->unknowns()[0]);
            toc("solve non linear heat equation");
            // Export the temperature, and flux on the radiating surfaces
            if(!boption(_name="deactivate-exporters"))
                this->exportHeat();
        }
        this->checkResults();
    } // end RHT<Dim,Order>::executeNonLinear()

    template<int Dim, int Order>
    void RHT<Dim,Order>::solveHeatEquationNonLinear(element_ptr_t T )
    {
        
        auto T_sol = M_Xh->elementPtr();
        auto T_vec = *T;
        *T_sol = *T;

        auto Res = backend()->newVector( M_Xh );
        auto Jac = backend()->newMatrix( 0,0,0,0, M_a->graph() );

        typename space_t::basis_type::points_type p(mesh_t::nDim,1);
        auto basispc = M_Xh->basis()->preCompute( M_Xh->basis(), p );   



        auto update_jacobian = [=]( const vector_ptrtype& T_vec, sparse_matrix_ptrtype& at_mat ) {

            auto at = form2(_trial=M_Xh, _test= M_Xh, _matrix=at_mat,_backend=backend(_name="heatEq"));
            auto a = form2(_trial=M_Xh, _test= M_Xh, _matrix=M_a,_backend=backend(_name="heatEq"));
            at = a;
            // Start from the time-independent forms computed in initHeatEquation            
            auto u = M_Xh->element(T_vec);                        

            // Radiative enclosure
            if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
            {
                for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
                {                
                    auto opI = opInterpolation(_domainSpace=M_Xh,
                                            _imageSpace=M_Xhd1_map[bc],                                    
                                            _backend=backend(_name="interp1",_rebuild=true),
                                            _type=InterpolationConforme());                    

                    auto fieldTemperature = M_Xh->element();
                    fieldTemperature.on(_range=elements(M_mesh),_expr=idv(u));
                    auto interpTemperatureHd1 = M_Xhd1_map[bc]->element();
                    opI->apply(fieldTemperature,interpTemperatureHd1);

                    // Elementwise mean of the temperature on the surface of the cavity
                    // The temperature is constant on each surface element
                    auto TmeanElementwise = interpTemperatureHd1.ewiseMean(M_Xhd0_map[bc]);                    

                    auto sigma = expr(value["sigma"].get<std::string>());

                    sigma.setParameterValues({{"sigma",specs["/Parameters/sigma"_json_pointer].get<double>()}});

                    auto sigma_double = specs["/Parameters/sigma"_json_pointer].get<double>();

                    // Fill the Jacobian matrix with the radiative heat transfer contributions:
                    // For each pair of cavity elements, compute the local coupling matrices and 
                    // dispatch them into the global jacobian matrix
                    auto current_range = markedfaces(this->M_mesh,value["markers"].get<std::vector<std::string>>());//elements(M_surface_submesh_map[bc]);
                    int i = 0;

                    auto the_im1 = im(this->M_mesh,3);//im( M_surface_submesh_map[bc],3);
                    auto current_pts1 = [&the_im1]( auto f, auto p )
                    {
                        //return the_im1.fpoints( f, p.value() );
                        return the_im1.points();
                    }; 

                    auto current_ctx = context( _element = boost::unwrap_ref( *begin( current_range ) ), _type = on_facets_t(), _geomap = this->M_mesh->gm()/*M_surface_submesh_map[bc]->gm()*/,_pointset = current_pts1 );
                    auto remote_ctx  = context( _element = boost::unwrap_ref( *begin( current_range ) ), _type = on_facets_t(), _geomap = this->M_mesh->gm()/*M_surface_submesh_map[bc]->gm()*/,_pointset = current_pts1 );

                    for(auto face_a : current_range)
                    {
                        auto const& f_a = boost::unwrap_ref( face_a );
                        int j=0;
                        current_ctx->template update<vm::POINT|vm::NORMAL|vm::JACOBIAN>(f_a.element0(),f_a.idInElement0());

                        // std::cout << "index_loop " << i << " face id heat solver " <<f_a.id() << std::endl;
                        for(auto face_b : current_range )
                        {
                            auto const& f_b = boost::unwrap_ref( face_b );
                            remote_ctx->template update<vm::POINT|vm::NORMAL|vm::JACOBIAN>(f_b.element0(),f_b.idInElement0());
                            if(f_a == f_b)
                            {
                                double Ta_average = 0;
                                // Compute the M^{aa} local matrix
                                double Sa = M_surface_submesh_local_areas[bc][i];
                                for( auto const& ldof : M_Xhd0_map[bc]->dof()->localDof( M_Xhd0_map[bc]->mesh()->meshToSubMesh(f_a.id()) ) )
                                {
                                    // std::cout << TmeanElementwise[  ldof.second.index() ] << std::endl;
                                    Ta_average += TmeanElementwise[  ldof.second.index() ];                                
                                }                                
                                double sum_Dab = M_cavity_coupling_mat[bc].row(i).sum() - M_cavity_coupling_mat[bc](i,j);
                                std::cout << std::scientific << "coeff_Mab " << M_cavity_coupling_mat[bc].row(i).sum() << " " <<  M_cavity_coupling_mat[bc](i,j) << std::endl;
                                double coeff_Maa = 4 * sigma_double * math::pow(Ta_average,3) * sum_Dab / math::pow(Sa,2);
                                coeff_Maa *= current_ctx->J(0) * remote_ctx->J(0); // multiply by jacobian geometric transformation Ja * Ja (for the two moment integrals)
                                std::cout << std::scientific << "coeff_Mab " << sigma_double << " " << math::pow(Ta_average,3) << " " << math::pow(Sa,2) << " " << sum_Dab << std::endl;
                                std::cout << std::scientific << "coeff_Mab2 " << current_ctx->J(0) * remote_ctx->J(0) << std::endl;


                                for( auto const& ldof : M_Xh->dof()->faceLocalDof( f_a.id() ) )
                                {
                                    // Find the values of the global dof index of the face
                                    index_type thedof = ldof.index();
                                    thedof = at.dofIdToContainerIdTrial()[ thedof ];
                                    for( auto const& ldof2 : M_Xh->dof()->faceLocalDof( f_b.id() ))
                                    {
                                        index_type thedof2 = ldof2.index();
                                        thedof2 = at.dofIdToContainerIdTrial()[ thedof2 ];
                                        double first_moment_a1 = basispc->firstMoment(ldof.localDof());
                                        double first_moment_a2 = basispc->firstMoment(ldof2.localDof());
                                        // std::cout << std::scientific << "dof,dof2 " << thedof << ", " << thedof2 << " contribution coeff_Maa" << coeff_Maa  << std::endl;
                                        // std::cout << std::scientific << "dof,dof2 " << thedof << ", " << thedof2 << " contribution first_moment_a1" <<  first_moment_a1  << std::endl;
                                        // std::cout << std::scientific << "dof,dof2 " << thedof << ", " << thedof2 << " contribution first_moment_a2" <<  first_moment_a2 << std::endl;
                                        
                                        at.set(thedof,thedof2, coeff_Maa * first_moment_a1 * first_moment_a2);
                                    }
                                }
                            }
                            else
                            {
                                // Compute the M^{ab} local matrix
                                double Tb_average = 0;
                                double Sa = M_surface_submesh_local_areas[bc][i];
                                double Sb = M_surface_submesh_local_areas[bc][j];
                                for( auto const& ldof :M_Xhd0_map[bc]->dof()->localDof( M_Xhd0_map[bc]->mesh()->meshToSubMesh(f_b.id()) ))
                                    Tb_average += TmeanElementwise[ ldof.second.index() ];
                                double Dab = M_cavity_coupling_mat[bc](i,j);

                                double coeff_Mab = - 4 * sigma_double * math::pow(Tb_average,3) * Dab / (Sa*Sb);                                
                                coeff_Mab *= current_ctx->J(0) * remote_ctx->J(0); // multiply by jacobian geometric transformation Ja * Ja (for the two moment integrals)

                                for( auto const& ldof : M_Xh->dof()->faceLocalDof( f_a.id() ) )
                                {
                                    // Find the values of the global dof index of the face
                                    index_type thedof = ldof.index();
                                    thedof = at.dofIdToContainerIdTrial()[ thedof ];
                                    for( auto const& ldof2 : M_Xh->dof()->faceLocalDof( f_b.id() ))
                                    {
                                        index_type thedof2 = ldof2.index();
                                        thedof2 = at.dofIdToContainerIdTrial()[ thedof2 ];
                                        double first_moment_a1 = basispc->firstMoment(ldof.localDof());
                                        double first_moment_a2 = basispc->firstMoment(ldof2.localDof());
                                        // std::cout << std::scientific << "dof,dof2 " << thedof << ", " << thedof2 << " contribution coeff_Mab" << coeff_Mab  << std::endl;
                                        // std::cout << std::scientific << "dof,dof2 " << thedof << ", " << thedof2 << " contribution first_moment_a1" <<  first_moment_a1  << std::endl;
                                        // std::cout << std::scientific << "dof,dof2 " << thedof << ", " << thedof2 << " contribution first_moment_a2" <<  first_moment_a2 << std::endl;

                                        at.set(thedof,thedof2, coeff_Mab * first_moment_a1 * first_moment_a2);
                                    }
                                }
                            }
                            j++;
                        }
                        i++;
                    }                 
                }
            } 
            at.close();
            auto res = backend()->newVector( M_Xh );
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/temperature"_json_pointer].items())
            {
                auto dirichletBc = expr(value["expr"].get<std::string>());
                at += on( _range=markedfaces(M_Xh->mesh(),bc), _element=u, _rhs=res, _expr=cst(0)*dirichletBc );                
            }                    
        };
        auto update_residual = [=]( const vector_ptrtype& T_vec, vector_ptrtype& lt_vec ) {

            ////////////////////////////////////
            // Add linear terms to the residual 
            ////////////////////////////////////
            
            auto lt = form1(_test=M_Xh, _vector=lt_vec);
            auto l = form1(_test=M_Xh, _vector=M_l);
            lt.zero();
            lt += l;             
            // std::cout << fmt::format("Norm lt vec initial {}\n",lt.vector().l2Norm());  

            // Start from the time-independent forms computed in initHeatEquation
            
            auto u = M_Xh->element(T_vec);
            auto v = M_Xh->element();

            // BC Neumann
            if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "flux" ) )
            {
                for ( auto& [bc, value] : specs["/BoundaryConditions/heat/flux"_json_pointer].items() )
                {
                    LOG( INFO ) << fmt::format( "flux {}: {}", bc, value.dump() );
                    auto flux = expr(value["expr"].get<std::string>());

                    lt += integrate( _range = markedfaces( M_mesh, bc ),
                            _expr = flux  * id( v ) );
                }
            }
            // std::cout << fmt::format("Norm lt vec initial {}\n",lt.vector().l2Norm());
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

                lt += integrate( _range = markedelements( M_mesh, material.get<std::string>() ), 
                        _expr =  inner( expr( k ) * gradv( u ) , grad( v ) ) );
                lt += integrate(_range = markedelements( M_mesh, material.get<std::string>() ), 
                        _expr = M_bdf->polyDerivCoefficient( 0 ) *expr( Rho ) * expr( Cp ) * idv( u ) * id( v ) );
            }

            // Update of RHS for heat equations starting from T with the term rho*C*T_{n-1}/dt (at order 1)
            for ( auto [key, material] : specs["/Models/heat/materials"_json_pointer].items() )
            {
                std::string matRho = fmt::format( "/Materials/{}/rho", material.get<std::string>() );
                std::string matCp = fmt::format( "/Materials/{}/Cp", material.get<std::string>() );
                auto Rho = specs[nl::json::json_pointer( matRho )].get<std::string>();
                auto Cp = specs[nl::json::json_pointer( matCp )].get<std::string>();

                lt += integrate( _range = markedelements( M_mesh, material.get<std::string>() ),                                                     
                                _expr = - expr( Rho ) * expr( Cp ) * idv(M_bdf->polyDeriv())  * id( v ) );                          
            }
            
            // Radiative enclosure
            // Integrating directly the radiative flux obtained via the solution of the problem coming from Modest book

             
            if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
            {
                for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
                {                
                    auto opI = opInterpolation(_domainSpace=M_Xh,
                                            _imageSpace=M_Xhd1_map[bc],                                    
                                            _backend=backend(_name="interp1",_rebuild=true),
                                            _type=InterpolationConforme());                    

                    auto fieldTemperature = M_Xh->element();
                    fieldTemperature.on(_range=elements(M_mesh),_expr=idv(u));
                    auto interpTemperatureHd1 = M_Xhd1_map[bc]->element();
                    opI->apply(fieldTemperature,interpTemperatureHd1);

                    auto the_im = im( this->M_mesh,3);//M_surface_submesh_map[bc], 3);
                    auto current_pts1 = [&the_im]( auto f, auto p )
                    {
                        return the_im.points();
                    }; 


                    // Elementwise mean of the temperature on the surface of the cavity
                    // The temperature is constant on each surface element
                    auto TmeanElementwise = interpTemperatureHd1.ewiseMean(M_Xhd0_map[bc]);
                    
                    auto sigma = expr(value["sigma"].get<std::string>());

                    sigma.setParameterValues({{"sigma",specs["/Parameters/sigma"_json_pointer].get<double>()}});
                    auto sigma_double = specs["/Parameters/sigma"_json_pointer].get<double>();
                    auto current_range =  markedfaces(this->M_mesh,value["markers"].get<std::vector<std::string>>());//elements(M_surface_submesh_map[bc]);

                    int i = 0;
                    auto current_ctx = context( _element = boost::unwrap_ref( *begin( current_range ) ), _type = on_facets_t(), _geomap = this->M_mesh->gm()/*M_surface_submesh_map[bc]->gm()*/,_pointset = current_pts1 );

                    for(auto face_a : current_range )
                    {                        
                        auto const& f_a = boost::unwrap_ref( face_a );
                        current_ctx->template update<vm::POINT|vm::NORMAL|vm::JACOBIAN>(f_a.element0(),f_a.idInElement0());
                        int j=0;
                        double Sa = M_surface_submesh_local_areas[bc][i];
                        double Ta_average = 0;
                        for( auto const& ldof :M_Xhd0_map[bc]->dof()->localDof( M_Xhd0_map[bc]->mesh()->meshToSubMesh(f_a.id()) ))
                        {
                            Ta_average += TmeanElementwise[ldof.second.index() ];
                        }
                        double coeff_residual =0;
                        for(auto face_b : current_range )
                        {
                            double Tb_average=0;
                            auto const& f_b = boost::unwrap_ref( face_b );                            
                            if(f_b != f_a)
                            {
                                // Compute the local radiative flux contributions                               
                                double Dab = M_cavity_coupling_mat[bc](i,j);
                                for( auto const& ldof :M_Xhd0_map[bc]->dof()->localDof( M_Xhd0_map[bc]->mesh()->meshToSubMesh(f_b.id()) ))
                                {
                                    Tb_average += TmeanElementwise[  ldof.second.index() ];
                                }
                                coeff_residual += sigma_double * (math::pow(Tb_average,4)-math::pow(Ta_average,4)) * Dab ;
                                //std::cout << std::scientific << "coeffresidual " << " f_a.id() " << f_a.id() << "meshtoSubmesh" << M_Xhd0_map[bc]->mesh()->meshToSubMesh(f_a.id())<<" f_b.id() " << f_b.id()  << sigma_double << " " << math::pow(Tb_average,4) << " " << math::pow(Ta_average,4) << " " << Dab << std::endl;                            
                            }
                            j++;
                        }
                        // Insert the contribution on the residual side
                        for( auto const& ldof : M_Xh->dof()->faceLocalDof( f_a.id() ) )
                        {
                            // Find the values of the global dof index of the face
                            index_type thedof = ldof.index();
                            thedof = l.dofIdToContainerId()[ thedof ];
                            double first_moment_a1 = basispc->firstMoment(ldof.localDof());
                            // std::cout << std::scientific << "dof,dof2 " << thedof << " contribution coeff_residual" << coeff_residual  << std::endl;
                            // std::cout << std::scientific << "dof,dof2 " << thedof << " contribution first_moment_a1" <<  first_moment_a1  << std::endl;
                            lt.set(thedof,  coeff_residual  *current_ctx->J(0)*  first_moment_a1 /Sa);                            
                        }
                        i++;
                    }

                }                
                                    // If the cavity is open, and it is assumed that a black body of fixed temperature
                                    // T_ref exchances heat with the cavity, and additional term is added: 
                                    // - the view factor is computed using the reciprocity formula FijAi = FjiAj
                                    // - the contribution is of the form \sigma T^4 * Fij
                                    // if(value["enclosure"]=="open")
                                    // {
                                    //     //std::cout << "Open enclosure " << bc << std::endl;
                                    //     double vf_marker_to_bb = 1 - M_matrix_vf_map[bc].row(i_marker_to_bb).sum(); //Fij
                                    //     double vf_bb_to_marker = 1 - M_matrix_vf_map[bc].col(i_marker_to_bb).sum(); //Fji
                                    //     auto measure_mark1 = integrate( _range=markedfaces(M_mesh,mark),_expr= cst(1.) ).evaluate()(0,0); //Ai
                                    //     double area_bb = measure_mark1 * vf_marker_to_bb / vf_bb_to_marker;
                                    //     auto T_bb = value["Tref"].get<double>();
                                    //     auto T_bb4 = cst(T_bb)*cst(T_bb)*cst(T_bb)*cst(T_bb);

                                        
                                    //     lt += integrate( _range = markedfaces( M_mesh, mark ),
                                    //         _expr =  sigma * expr( epsilon ) * T_bb4 * cst(vf_bb_to_marker) * id( v ));

                                    // }
                                    
                //                     i_marker_to_bb++;

                //                     done = 1;
                //                     break;
                //                 }
                //             }
                //             if ( done ){ break; }
                //         }    
                //     }
                // }
            }            
            // if(specs["/BoundaryConditions/heat"_json_pointer].contains( "solar_radiation" ) )
            // {
            //     for ( auto& [meteo_station_name, value] : specs["/BoundaryConditions/heat/solar_radiation"_json_pointer].items() )
            //     {
                    
            //         auto current_Tsky = accessTsky(meteo_station_name,M_bdf->time(), M_bdf->timeStep());
            //         current_Tsky += 273.15; // transform in Kelvin
            //         //std::cout <<  fmt::format("Current Tsky  {} meteo station {}", current_Tsky,meteo_station_name);
            //         auto buildings_specs_string = fmt::format("/BoundaryConditions/heat/solar_radiation/{}/buildings",meteo_station_name);
                    
            //         for(auto const& [building_name,building_structure] : specs[json::json_pointer(buildings_specs_string)].items() )
            //         {
            //             for(std::string surface_name : building_structure.at("surfaces"))
            //             {
                            
            //                 auto current_solarRad = accessSolarRadiation(meteo_station_name, building_name, surface_name, M_bdf->time(), M_bdf->timeStep());
            //                 //std::cout <<  fmt::format("Current solar radiation {} surface {}", current_solarRad,surface_name);
            //                 // Add solar radiation as Neumann boundary condition on the exposed surfaces
            //                 lt += integrate( _range = markedfaces( M_mesh, surface_name ),
            //                                 _expr =  cst(current_solarRad) * id( v ));
                            
            //                 // Add black-body radiation from sky temperature on both 
            //                 auto sigma = expr(specs["/Parameters/sigma"_json_pointer].get<double>());   
                            
            //                 auto Tsky4 = cst(current_Tsky) * cst(current_Tsky) * cst(current_Tsky) * cst(current_Tsky);  

            //                 auto idvu3 = idv( T ) * idv( T ) * idv( T );

            //                 lt += integrate( _range = markedfaces( M_mesh, surface_name ), 
            //                                 _expr = sigma * (Tsky4 -idvu3 * idv( T ))  * id( v ) );                                                                                                                                            
                            
            //             }
            //         }
            //     }
            // }
            lt_vec->close();
            auto temp = M_Xh->element();
            temp = *lt_vec;
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/temperature"_json_pointer].items())
            {
                auto dirichletBc = expr(value["expr"].get<std::string>());
                temp.on( _range=markedfaces(M_mesh,bc),_expr=cst(0)*dirichletBc );
            }

            *lt_vec = temp;

            // std::cout << "T_vec" << *T_vec << std::endl;

        };
        
        for ( auto& [bc, value] : specs["/BoundaryConditions/heat/temperature"_json_pointer].items())
        {
            auto dirichletBc = expr(value["expr"].get<std::string>());
            T_sol->on( _range=markedfaces(M_mesh,bc),_expr=dirichletBc );
        }

        auto e = exporter(_mesh=M_Xh->mesh(),_name="beforeSol");
        e->add("T_sol",*T_sol);
        e->save();

        backend(_name="heatEq" )->nlSolver()->jacobian = update_jacobian;
        backend(_name="heatEq" )->nlSolver()->residual = update_residual;
        backend(_name="heatEq")->nlSolve( _solution=T_sol,_jacobian=Jac,_residual=Res );

        auto e1 = exporter(_mesh=M_Xh->mesh(),_name="afterSol");
        e1->add("T_sol",*T_sol);
        e1->save();

        M_currentTempAndFlux.setT(T_sol);
        
        auto newT = unwrap_ptr(M_currentTempAndFlux.T());
        M_bdf->next(newT);    

        

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
        // rht.execute();
        rht.executeNonLinear();

    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
} // end main()