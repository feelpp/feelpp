
#include "rht_readers.hpp"

namespace Feel
{      


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

        // BC Robin -k*grad(u)*N = h*(T-Text)
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "convective_heat_flux" ) )
        {
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/convective_heat_flux"_json_pointer].items() )
            {
                LOG( INFO ) << fmt::format( "convective_heat_flux {}: {}", bc, value.dump() );
                auto h = value["h"].get<std::string>();
                // auto Text = expr(value["Text"].get<std::string>());

                // Text.setParameterValues({{"Tref_C",specs["/Parameters/Tref_C"_json_pointer].get<double>()}});

                a += integrate( _range = markedfaces( M_mesh, bc ),
                        _expr = expr( h ) * id( v ) * idt( u ) );
                // l += integrate( _range = markedfaces( M_mesh, bc ),
                //         _expr =  expr( h ) * Text * id( v ) );
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
        auto e_opInterp2 = exporter(_mesh=M_mesh,_name="solveHeat");

        e_opInterp2->add("Tinit",idv(T));

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

                                lt += integrate( _range = markedfaces( M_mesh, mark ),
                                        _expr = - sigma * expr( epsilon ) * idvu3 * idv( T ) * id( v ));

                                at += integrate( _range = markedfaces( M_mesh, mark ),
                                        _expr = cst(4) * sigma * expr( epsilon ) * idvu3 * idt( u ) * id( v ));
                                done = 1;
                                break;
                            }
                        }
                        if ( done ){ break; }
                    }    
                }
            }
        }

        if(specs["/BoundaryConditions/heat"_json_pointer].contains( "solar_radiation" ) )
        {
            for ( auto& [meteo_station_name, value] : specs["/BoundaryConditions/heat/solar_radiation"_json_pointer].items() )
            {
                
                auto current_Tsky = accessTsky(meteo_station_name,M_bdf->time(), M_bdf->timeStep());
                current_Tsky += 273.15; // transform in Kelvin
                //std::cout <<  fmt::format("Current Tsky  {} meteo station {}", current_Tsky,meteo_station_name);
                auto buildings_specs_string = fmt::format("/BoundaryConditions/heat/solar_radiation/{}/buildings",meteo_station_name);
                
                for(auto const& [building_name,building_structure] : specs[json::json_pointer(buildings_specs_string)].items() )
                {
                    for(std::string surface_name : building_structure.at("surfaces"))
                    {
                        
                        auto current_solarRad = accessSolarRadiation(meteo_station_name, building_name, surface_name, M_bdf->time(), M_bdf->timeStep());
                        //std::cout <<  fmt::format("Current solar radiation {} surface {}", current_solarRad,surface_name);
                        // Add solar radiation as Neumann boundary condition on the exposed surfaces
                        lt += integrate( _range = markedfaces( M_mesh, surface_name ),
                                        _expr =  cst(current_solarRad) * id( v ));
                        
                        // Add black-body radiation from sky temperature on both 
                        auto sigma = expr(specs["/Parameters/sigma"_json_pointer].get<double>());   
                        
                        auto Tsky4 = cst(current_Tsky) * cst(current_Tsky) * cst(current_Tsky) * cst(current_Tsky);  

                        auto idvu3 = idv( T ) * idv( T ) * idv( T );

                        lt += integrate( _range = markedfaces( M_mesh, surface_name ), 
                                        _expr = - sigma * (Tsky4 -idvu3 * idv( T ))  * id( v ) );                        
                        
                        
                        at += integrate( _range = markedfaces( M_mesh, surface_name ),
                                _expr =  cst(4) * sigma * idvu3 * idt( u ) * id( v ));
                        
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

        e_opInterp2->add("Tsol",idv(T_sol));
        e_opInterp2->save();

        M_currentTempAndFlux.setT(T_sol);

    } // end RHT<Dim,Order>::solveHeatEquation


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
        if(n_cavities != 0)
        {
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
        }
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

        read_Tsky();
        read_solarRadiation();
        read_Text();

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
        // this->solveHeatEquationNonLinear(M_bdf->unknowns()[0]);  
        
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
            // this->solveHeatEquationNonLinear(T_new);  
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

    // Solve non linear Heat equation with radiative BC
    template<int Dim, int Order>
    void RHT<Dim,Order>::solveHeatEquationNonLinear(element_ptr_t T )
    {
        // Recover the initial guess of the nonlinear iteration
        auto T_sol = M_Xh->elementPtr();
        auto T_vec = *T;
        *T_sol = *T;

        auto e_opInterp2 = exporter(_mesh=M_mesh,_name="solveHeatNonLin");

        e_opInterp2->add("Tinit",idv(T));

        auto Res = backend()->newVector( M_Xh );
        auto Jac = backend()->newMatrix( _test=M_Xh, _trial=M_Xh );


        auto update_jacobian  = [=]( const vector_ptrtype& T_vec, sparse_matrix_ptrtype& at_mat ) {

                            auto at = form2(_trial=M_Xh, _test= M_Xh, _matrix=at_mat,_backend=backend(_name="heatEq"));
                            auto a = form2(_trial=M_Xh, _test= M_Xh, _matrix=M_a,_backend=backend(_name="heatEq"));
                            
                            // Start from the time-independent forms computed in initHeatEquation
                            at = a;
                            
                            auto u = M_Xh->element();
                            auto v = M_Xh->element();
                            u=*T_vec;

                            if(specs["/BoundaryConditions/heat"_json_pointer].contains( "solar_radiation" ) )
                            {
                                for ( auto& [meteo_station_name, value] : specs["/BoundaryConditions/heat/solar_radiation"_json_pointer].items() )
                                {
                                    
                                    auto current_Tsky = accessTsky(meteo_station_name, M_bdf->time(), M_bdf->timeStep());
                                    current_Tsky += 273.15; // transform in Kelvin
                                    // std::cout <<  fmt::format("Current Tsky  {} meteo station {}", current_Tsky,meteo_station_name);
                                    auto buildings_specs_string = fmt::format("/BoundaryConditions/heat/solar_radiation/{}/buildings",meteo_station_name);
                                    
                                    for(auto const& [building_name,building_structure] : specs[json::json_pointer(buildings_specs_string)].items() )
                                    {
                                        for(std::string surface_name : building_structure["OpaqueSurfaces"]["surfaces"])
                                        {
                                                                                                
                                            // Add black-body radiation from sky temperature on both 
                                            auto sigma = expr(specs["/Parameters/sigma"_json_pointer].get<double>());                                                                                                       

                                            auto idvu3 = idv( u ) * idv( u ) * idv( u );
                                            
                                            at += integrate( _range = markedfaces( M_mesh, surface_name ),
                                                    _expr = cst(4) * sigma * idvu3 * idt( u ) * id( v ));
                                            
                                        }
                                        for(std::string surface_name : building_structure["TransparentSurfaces"]["surfaces"])
                                        {
                                            // Add U-value term for the transparent surfaces
                                            auto Uvalue = expr(building_structure["TransparentSurfaces"]["Uvalue"].get<double>()); 

                                            at += integrate( _range = markedfaces( M_mesh, surface_name ),
                                                    _expr = Uvalue * idt( u ) * id( v ));
                                            
                                        }
                                    }
                                }
                            }

                    // Radiative cavity
                            if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
                            {
                                for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
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
                                                    auto idvu3 = idv( u ) * idv( u ) * idv( u );                                                            
                                                    // Term 4 sigma * epsilon * T^3 * increment
                                                    at += integrate( _range = markedfaces( M_mesh, mark ),
                                                            _expr = cst(4) *  sigma * expr( epsilon ) * idvu3 * idt( u ) * id( v ));
                                                    // Term - 4 sigma * epsilon \int T^3 dFij * increment

                                                    // auto cavity_markers = M_markers_map[bc];   

                                                    int i_marker_to_bb=0;
                                                    // int i_mark=0; 
                                                    // int j_mark=0; 
                                                    // auto vf_field = M_Xhd0_map[bc]->element();
                                                    // for ( std::string mark2 : cavity_markers ) // Loop on index j
                                                    // {
                                                    //     // Compute the index relative to marker mark
                                                    //     auto it = std::find(cavity_markers.begin(),cavity_markers.end(),mark);                    
                                                    //     i_mark = it - cavity_markers.begin();
                                                    //     // Compute the index relative to marker mark2
                                                    //     auto jt = std::find(cavity_markers.begin(),cavity_markers.end(),mark2);                    
                                                    //     j_mark = jt - cavity_markers.begin();
                                                    //     //std::cout << fmt::format("we add a contribution from index {} whose view factor is {}",j_mark,M_vf_matrix(i_mark,j_mark))<< std::endl;
                                                        
                                                    //     vf_field.zero();                                                       
                                                    //     //std::cout << fmt::format("we add a contribution from marker {} whose emissivity is {}",markcoat,epsilon_mark2)<< std::endl;
                                                    //     vf_field.on(_range=markedfaces(M_mesh,mark2),_expr=cst(M_matrix_vf_map[bc](i_mark,j_mark)));

                                                    //     at += integrate( _range = markedfaces( M_mesh, mark ),
                                                    //         _expr = cst(4) *  sigma * cst(M_matrix_vf_map[bc](i_mark,j_mark)) * expr( epsilon ) * idvu3 * idt( u ) * id( v ));
                                                        
                                                    //     // auto T4_scal_vf = integrate(_range=markedfaces(M_mesh,mark2), _expr=inner(idv(T4_interp),idv(vf_field))).evaluate()(0,0);
                                                    //     auto T3jacobian_scal_vf = integrate(_range=markedfaces(M_mesh,mark2), _expr=inner(cst(4) *  sigma * expr( epsilon ) * idvu3,idv(vf_field))).evaluate()(0,0);

                                                    //     auto measure_mark2 = integrate( _range=markedfaces(M_mesh,mark2),_expr= cst(1.) ).evaluate()(0,0);
                                                        
                                                    //     at += integrate( _range=markedfaces(M_mesh,mark),
                                                    //                     _expr= cst(-1)/cst(measure_mark2)*cst(T3jacobian_scal_vf) * idt(u) *id(v) );  
                                                    //     // auto done2=0;
                                                    //     // for ( auto& [key2, coat2] : specs["/Coating"_json_pointer].items() )
                                                    //     // {
                                                    //     //     for ( auto markcoat2 : coat2.at("markers") )
                                                    //     //     {
                                                    //     //         if ( mark2 == markcoat2 )
                                                    //     //         {
                                                    //     //             auto epsilon_mark2 = coat2["epsilon"].get<std::string>(); 

                                                    //     //             at += integrate( _range=markedfaces(M_mesh,mark),
                                                    //     //                             _expr = expr(epsilon) * (cst(1. )/expr(epsilon_mark2)-cst(1.))/cst(measure_mark2) * cst(400) * gradt(u) *N() * id(v) );
                                                    //     //         done2 = 1;
                                                    //     //         break;
                                                    //     //         }
                                                    //     //     }
                                                    //     //     if ( done ){ break; }
                                                    //     // }
                                                                            
                                                    // }

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
                                                        auto T_bb3 = cst(T_bb)*cst(T_bb)*cst(T_bb);

                                                        
                                                        at += integrate( _range = markedfaces( M_mesh, mark ),
                                                            _expr = cst(4) *  sigma * expr( epsilon ) * T_bb3 * cst(vf_bb_to_marker) * idt( u ) * id( v ));

                                                    }
                                                    
                                                    i_marker_to_bb++;

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
                                auto RR = backend()->newVector( M_Xh );
                                auto dirichletBc = expr(value["expr"].get<std::string>());
                                at+=on( _range=markedfaces(M_mesh,bc), _rhs=RR, _element=u, _expr=cst(0)*dirichletBc );
                            }
        
        
        };

        auto update_residual = [=]( const vector_ptrtype& T_vec, vector_ptrtype& lt_vec ) {

                            auto lt = form1(_test=M_Xh, _vector=lt_vec);
                            auto l = form1(_test=M_Xh, _vector=M_l);
                            // Start from the time-independent forms computed in initHeatEquation

                            lt.zero();
                            lt += l;
                        
                            auto u = M_Xh->element();
                            auto v = M_Xh->element();

                            u=*T_vec;
                        
                            // for ( auto& [meteo_station_name_Text, value] : specs["/BoundaryConditions/heat/meteo"_json_pointer].items() )
                            // {
                            //     auto current_Text = accessText(meteo_station_name_Text, M_bdf->time(), M_bdf->timeStep());
                            //     auto Text = expr(current_Text + 273.15); // transform in Kelvin

                            // For each material, integrate rho*C*d/dt idt(u)*id(v + \int_MAT k_MAT * grad(u) * grad(v)
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
                                    lt += integrate( _range = markedelements( M_mesh, material.get<std::string>() ),                                                     
                                                        _expr = - expr( Rho ) * expr( Cp ) * idv(M_bdf->polyDeriv())  * id( v ) );    
                                }  

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

                            // BC Robin -k*grad(u)*N = h*(T-Text)
                                // if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "convective_heat_flux" ) )
                                // {
                                //     for ( auto& [bc, value] : specs["/BoundaryConditions/heat/convective_heat_flux"_json_pointer].items() )
                                //     {
                                //         LOG( INFO ) << fmt::format( "convective_heat_flux {}: {}", bc, value.dump() );
                                //         auto h = value["h"].get<std::string>();
                                //         // auto Text = expr(value["Text"].get<std::string>());

                                //         //Text.setParameterValues({{"Tref_C",specs["/Parameters/Tref_C"_json_pointer].get<double>()}});

                                //         lt += integrate( _range = markedfaces( M_mesh, bc ),
                                //                 _expr = expr( h ) * (idv( u ) - Text) * id( v )  );
                                //     }
                                // }
                            // Solar radiation
                                // if(specs["/BoundaryConditions/heat"_json_pointer].contains( "solar_radiation" ) )
                                // {
                                //     for ( auto& [meteo_station_name, value] : specs["/BoundaryConditions/heat/solar_radiation"_json_pointer].items() )
                                //     {
                                        
                                //         auto current_Tsky = accessTsky(meteo_station_name, M_bdf->time(), M_bdf->timeStep());
                                //         current_Tsky += 273.15; // transform in Kelvin
                                        
                                //         // std::cout <<  fmt::format("Current Tsky residual  {} meteo station {}", current_Tsky, meteo_station_name);

                                //         auto buildings_specs_string = fmt::format("/BoundaryConditions/heat/solar_radiation/{}/buildings",meteo_station_name);
                                        
                                //         for(auto const& [building_name,building_structure] : specs[json::json_pointer(buildings_specs_string)].items() )
                                //         {
                                //             for(std::string surface_name : building_structure["OpaqueSurfaces"]["surfaces"])
                                //             {
                                                
                                //                 auto current_solarRad = accessSolarRadiation(meteo_station_name, building_name, surface_name, M_bdf->time(), M_bdf->timeStep());
                                //                 // std::cout <<  fmt::format("Current solar radiation {} surface {}", current_solarRad,surface_name);
                                //                 auto absorptivity = building_structure["OpaqueSurfaces"]["absorptivity"].get<double>();
                                //                 // Add solar radiation as Neumann boundary condition on the exposed surfaces
                                //                 lt += integrate( _range = markedfaces( M_mesh, surface_name ),
                                //                                 _expr = - expr(absorptivity) * cst(current_solarRad) * id( v ));
                                                
                                //                 // Add black-body radiation from sky temperature on both 
                                //                 auto sigma = expr(specs["/Parameters/sigma"_json_pointer].get<double>());   
                                                
                                //                 auto Tsky4 = cst(current_Tsky) * cst(current_Tsky) * cst(current_Tsky) * cst(current_Tsky);  

                                //                 auto idvu3 = idv( T ) * idv( T ) * idv( T );

                                //                 lt += integrate( _range = markedfaces( M_mesh, surface_name ), 
                                //                                 _expr = sigma * (idvu3 * idv( T ) - Tsky4)  * id( v ) );                                                                                                                                            
                                                
                                //                 // auto mean_surf = integrate(_range = markedfaces( M_mesh, surface_name ), _expr = cst(1.0) ).evaluate()(0,0);
                                //                 // auto mean_surf_temp = integrate(_range = markedfaces( M_mesh, surface_name ), _expr = idv( T ) ).evaluate()(0,0);
                                //                 // auto mean_surf_tsky = integrate(_range = markedfaces( M_mesh, surface_name ), _expr = cst(current_Tsky) ).evaluate()(0,0);
                                //                 // auto mean_surf_t4 = integrate(_range = markedfaces( M_mesh, surface_name ), _expr = sigma * (idvu3 * idv( T ) - Tsky4) ).evaluate()(0,0);

                                //                 // std::cout << fmt::format("Surface {} external temp {} mean temperature {} mean temperature sky {}  mean temperature power {}",surface_name, current_Text, mean_surf_temp/mean_surf, mean_surf_tsky/mean_surf, mean_surf_t4/mean_surf) << std::endl;
                                //             }

                                //             for(std::string surface_name : building_structure["TransparentSurfaces"]["surfaces"])
                                //             {
                                                
                                //                 auto current_solarRad = accessSolarRadiation(meteo_station_name, building_name, surface_name, M_bdf->time(), M_bdf->timeStep());                                            

                                //                 // Read SHGC term for the transparent surfaces
                                //                 auto SHGC = expr(building_structure["TransparentSurfaces"]["SHGC"].get<double>());

                                //                 std::cout <<  fmt::format("Current solar radiation {} surface {}", current_solarRad,surface_name);
                                //                 // Add solar radiation as Neumann boundary condition on the exposed surfaces
                                //                 lt += integrate( _range = markedfaces( M_mesh, surface_name ),
                                //                                 _expr = - SHGC * cst(current_solarRad) * id( v ));

                                //                 // Add U-value term for the transparent surfaces
                                //                 auto Uvalue = expr(building_structure["TransparentSurfaces"]["Uvalue"].get<double>());                                             

                                //                 lt += integrate( _range = markedfaces( M_mesh, surface_name ), 
                                //                                 _expr = Uvalue * (idv( T ) - Text)  * id( v ) );                                                                                                                                            
                                                
                                //             }
                                //         }
                                //     }
                                // }

                         // Integrating directly the radiative flux obtained via the solution of the problem coming from Modest book
                        // Radiative enclosure
                                if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
                                {
                                    for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
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
                                                        auto idvu4 = idv( u ) * idv( u ) * idv( u )* idv( u ) ;                                                            
                                                        // Term sigma * epsilon * T^4 
                                                        lt += integrate( _range = markedfaces( M_mesh, mark ),
                                                                _expr = sigma * expr( epsilon ) * idvu4 *id( v ));
                                                        
                                                        // Term - sigma * epsilon \int T^4 dFij 

                                                        auto cavity_markers = M_markers_map[bc];   

                                                        int i_marker_to_bb=0;
                                                        int i_mark=0; 
                                                        int j_mark=0; 
                                                        auto vf_field = M_Xhd0_map[bc]->element();
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
                                                            
                                                            auto T4residual_scal_vf = integrate(_range=markedfaces(M_mesh,mark2), _expr=inner(sigma * expr( epsilon ) * idvu4,idv(vf_field))).evaluate()(0,0);

                                                            auto measure_mark2 = integrate( _range=markedfaces(M_mesh,mark2),_expr= cst(1.) ).evaluate()(0,0);
                                                            auto measure_mark1 = integrate( _range=markedfaces(M_mesh,mark),_expr= cst(1.) ).evaluate()(0,0);
                                                            
                                                            lt += integrate( _range=markedfaces(M_mesh,mark),
                                                                            _expr= cst(-1)*cst(T4residual_scal_vf)/cst(measure_mark2) *id(v) );  

                                                            auto done2=0;
                                                            for ( auto& [key2, coat2] : specs["/Coating"_json_pointer].items() )
                                                            {
                                                                for ( auto markcoat2 : coat2.at("markers") )
                                                                {
                                                                    if ( mark2 == markcoat2 )
                                                                    {
                                                                        auto epsilon_mark2 = coat2["epsilon"].get<std::string>(); 

                                                                        auto flux_other_surface = integrate(_range=markedfaces(M_mesh,mark2), _expr= cst(400) * gradv(u) * idv(vf_field) * N() ).evaluate()(0,0);

                                                                        lt += integrate( _range=markedfaces(M_mesh,mark),
                                                                                        _expr = cst(-1)*expr(epsilon)* (cst(1. ) /expr(epsilon_mark2)-cst(1.))/cst(measure_mark2) * cst(flux_other_surface) * id(v) );
                                                                        
                                                                        done2 = 1;
                                                                        break;
                                                                    }
                                                                }
                                                                if ( done ){ break; }
                                                            }
                                                                                
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

                                                            
                                                            lt += integrate( _range = markedfaces( M_mesh, mark ),
                                                                _expr =  sigma * expr( epsilon ) * T_bb4 * cst(vf_bb_to_marker) * id( v ));

                                                        }
                                                        
                                                        i_marker_to_bb++;

                                                        done = 1;
                                                        break;
                                                    }
                                                }
                                                if ( done ){ break; }
                                            }    
                                        }
                                    }
                                }
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

        };
                                

        // Lambda function creating and updating the Jacobian matrix
        // auto update_jacobian2 = [=]( const vector_ptrtype& T_vec, sparse_matrix_ptrtype& at_mat ) {
                                  
        //                         auto at = form2(_trial=M_Xh, _test= M_Xh, _matrix=at_mat,_backend=backend(_name="heatEq"));
        //                         auto a = form2(_trial=M_Xh, _test= M_Xh, _matrix=M_a,_backend=backend(_name="heatEq"));
                                
        //                         // Start from the time-independent forms computed in initHeatEquation
        //                         at = a;
                                
        //                         auto u = M_Xh->element();
        //                         auto v = M_Xh->element();
        //                         u=*T_vec;
        //                         auto int_temp = integrate(_range=elements(M_mesh),_expr=idv(u)).evaluate()(0,0);

        //                 // Blackbody radiative condition
        //                         // Integrating 4*sigma*epsilon*T_{n-1}^3*idt(u)*id(v)
        //                         if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_blackbody_heat_flux" ) )
        //                         {
        //                             for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_blackbody_heat_flux"_json_pointer].items() )
        //                             {                
        //                                 auto sigma = expr(value["sigma"].get<std::string>());

        //                                 sigma.setParameterValues({{"sigma",specs["/Parameters/sigma"_json_pointer].get<double>()}});

        //                                 // Recover the emissivity epsilon of the coating material associated to the marked face
        //                                 for ( std::string mark :  M_markers_map[bc] )
        //                                 {
        //                                     auto done = 0;
        //                                     for ( auto& [key, coat] : specs["/Coating"_json_pointer].items() )
        //                                     {
        //                                         for ( auto markcoat : coat.at("markers") )
        //                                         {
        //                                             if ( mark == markcoat )
        //                                             {
        //                                                 auto epsilon = coat["epsilon"].get<std::string>();                            
        //                                                 auto idvu3 = idv( u ) * idv( u ) * idv( u );                                                            

        //                                                 at += integrate( _range = markedfaces( M_mesh, mark ),
        //                                                         _expr = cst(4) *  sigma * expr( epsilon ) * idvu3 * idt( u ) * id( v ));
        //                                                 done = 1;
        //                                                 break;
        //                                             }
        //                                         }
        //                                         if ( done ){ break; }
        //                                     }    
        //                                 }
        //                             }
        //                         }

        //                 // Radiative enclosure
        //                         if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
        //                         {
        //                             for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
        //                             {                
        //                                 auto sigma = expr(value["sigma"].get<std::string>());

        //                                 sigma.setParameterValues({{"sigma",specs["/Parameters/sigma"_json_pointer].get<double>()}});
                                        
        //                                 // Recover the emissivity epsilon of the coating material associated to the marked face
        //                                 for ( std::string mark :  M_markers_map[bc] )
        //                                 {
        //                                     auto done = 0;
        //                                     for ( auto& [key, coat] : specs["/Coating"_json_pointer].items() )
        //                                     {
        //                                         for ( auto markcoat : coat.at("markers") )
        //                                         {
        //                                             if ( mark == markcoat )
        //                                             {
        //                                                 auto epsilon = coat["epsilon"].get<std::string>();                            
        //                                                 auto idvu3 = idv( u ) * idv( u ) * idv( u );                                                            
        //                                                 // Term 4 sigma * epsilon * T^3 * increment
        //                                                 at += integrate( _range = markedfaces( M_mesh, mark ),
        //                                                         _expr = cst(4) *  sigma * expr( epsilon ) * idvu3 * idt( u ) * id( v ));
        //                                                 // Term - 4 sigma * epsilon \int T^3 dFij * increment

        //                                                 auto cavity_markers = M_markers_map[bc];   

        //                                                 int i_marker_to_bb=0;
        //                                                 int i_mark=0; 
        //                                                 int j_mark=0; 
        //                                                 auto vf_field = M_Xhd0_map[bc]->element();
        //                                                 for ( std::string mark2 : cavity_markers ) // Loop on index j
        //                                                 {
        //                                                     // Compute the index relative to marker mark
        //                                                     auto it = std::find(cavity_markers.begin(),cavity_markers.end(),mark);                    
        //                                                     i_mark = it - cavity_markers.begin();
        //                                                     // Compute the index relative to marker mark2
        //                                                     auto jt = std::find(cavity_markers.begin(),cavity_markers.end(),mark2);                    
        //                                                     j_mark = jt - cavity_markers.begin();
        //                                                     //std::cout << fmt::format("we add a contribution from index {} whose view factor is {}",j_mark,M_vf_matrix(i_mark,j_mark))<< std::endl;
                                                            
        //                                                     vf_field.zero();                                                       
        //                                                     //std::cout << fmt::format("we add a contribution from marker {} whose emissivity is {}",markcoat,epsilon_mark2)<< std::endl;
        //                                                     vf_field.on(_range=markedfaces(M_mesh,mark2),_expr=cst(M_matrix_vf_map[bc](i_mark,j_mark)));
                                                            
        //                                                     // auto T4_scal_vf = integrate(_range=markedfaces(M_mesh,mark2), _expr=inner(idv(T4_interp),idv(vf_field))).evaluate()(0,0);
        //                                                     auto T3jacobian_scal_vf = integrate(_range=markedfaces(M_mesh,mark2), _expr=inner(cst(4) *  sigma * expr( epsilon ) * idvu3,idv(vf_field))).evaluate()(0,0);

        //                                                     auto measure_mark2 = integrate( _range=markedfaces(M_mesh,mark2),_expr= cst(1.) ).evaluate()(0,0);
                                                            
        //                                                     at += integrate( _range=markedfaces(M_mesh,mark),
        //                                                                     _expr= cst(-1)/cst(measure_mark2)*cst(T3jacobian_scal_vf) * idt(u) *id(v) );  
        //                                                     auto done2=0;
        //                                                     for ( auto& [key2, coat2] : specs["/Coating"_json_pointer].items() )
        //                                                     {
        //                                                         for ( auto markcoat2 : coat2.at("markers") )
        //                                                         {
        //                                                             if ( mark2 == markcoat2 )
        //                                                             {
        //                                                                 auto epsilon_mark2 = coat2["epsilon"].get<std::string>(); 

        //                                                                 at += integrate( _range=markedfaces(M_mesh,mark),
        //                                                                                 _expr = expr(epsilon) * (cst(1. )/expr(epsilon_mark2)-cst(1.))/cst(measure_mark2) * cst(400) * gradt(u) *N() * id(v) );
        //                                                             done2 = 1;
        //                                                             break;
        //                                                             }
        //                                                         }
        //                                                         if ( done ){ break; }
        //                                                     }
                                                                                
        //                                                 }

        //                                                 // If the cavity is open, and it is assumed that a black body of fixed temperature
        //                                                 // T_ref exchances heat with the cavity, and additional term is added: 
        //                                                 // - the view factor is computed using the reciprocity formula FijAi = FjiAj
        //                                                 // - the contribution is of the form \sigma T^4 * Fij
        //                                                 if(value["enclosure"]=="open")
        //                                                 {
        //                                                     //std::cout << "Open enclosure " << bc << std::endl;
        //                                                     double vf_marker_to_bb = 1 - M_matrix_vf_map[bc].row(i_marker_to_bb).sum(); //Fij
        //                                                     double vf_bb_to_marker = 1 - M_matrix_vf_map[bc].col(i_marker_to_bb).sum(); //Fji
        //                                                     auto measure_mark1 = integrate( _range=markedfaces(M_mesh,mark),_expr= cst(1.) ).evaluate()(0,0); //Ai
        //                                                     double area_bb = measure_mark1 * vf_marker_to_bb / vf_bb_to_marker;
        //                                                     auto T_bb = value["Tref"].get<double>();
        //                                                     auto T_bb3 = cst(T_bb)*cst(T_bb)*cst(T_bb);

                                                            
        //                                                     at += integrate( _range = markedfaces( M_mesh, mark ),
        //                                                         _expr = cst(4) *  sigma * expr( epsilon ) * T_bb3 * cst(vf_bb_to_marker) * idt( u ) * id( v ));

        //                                                 }
                                                        
        //                                                 i_marker_to_bb++;

        //                                                 done = 1;
        //                                                 break;
        //                                             }
        //                                         }
        //                                         if ( done ){ break; }
        //                                     }    
        //                                 }
        //                             }
        //                         }

        //                 // Solar radiation boundary condition
        //                         if(specs["/BoundaryConditions/heat"_json_pointer].contains( "solar_radiation" ) )
        //                         {
        //                             for ( auto& [meteo_station_name, value] : specs["/BoundaryConditions/heat/solar_radiation"_json_pointer].items() )
        //                             {
                                        
        //                                 auto current_Tsky = accessTsky(meteo_station_name, M_bdf->time(), M_bdf->timeStep());
        //                                 current_Tsky += 273.15; // transform in Kelvin
        //                                 // std::cout <<  fmt::format("Current Tsky  {} meteo station {}", current_Tsky,meteo_station_name);
        //                                 auto buildings_specs_string = fmt::format("/BoundaryConditions/heat/solar_radiation/{}/buildings",meteo_station_name);
                                        
        //                                 for(auto const& [building_name,building_structure] : specs[json::json_pointer(buildings_specs_string)].items() )
        //                                 {
        //                                     for(std::string surface_name : building_structure.at("surfaces"))
        //                                     {
                                                                                                    
        //                                         // Add black-body radiation from sky temperature on both 
        //                                         auto sigma = expr(specs["/Parameters/sigma"_json_pointer].get<double>());                                                                                                       

        //                                         auto idvu3 = idv( u ) * idv( u ) * idv( u );
                                                
        //                                         // at += integrate( _range = markedfaces( M_mesh, surface_name ),
        //                                         //         _expr = - cst(4) * sigma * idvu3 * idt( u ) * id( v ));
                                                
        //                                     }
        //                                 }
        //                             }
        //                         }
        //                         for ( auto& [bc, value] : specs["/BoundaryConditions/heat/temperature"_json_pointer].items())
        //                         {
        //                             auto RR = backend()->newVector( M_Xh );
        //                             auto dirichletBc = expr(value["expr"].get<std::string>());
        //                             at+=on( _range=markedfaces(M_mesh,bc), _rhs=RR, _element=u, _expr=cst(0)*dirichletBc );
        //                         }
        //                     };

        // Lambda function creating and updating the residual vector
        // auto update_residual2 = [=]( const vector_ptrtype& T_vec, vector_ptrtype& lt_vec ) {
                                
        //                         auto lt = form1(_test=M_Xh, _vector=lt_vec);
        //                         auto l = form1(_test=M_Xh, _vector=M_l);
        //                         // Start from the time-independent forms computed in initHeatEquation

        //                         lt.zero();
        //                         lt += l;
                            
        //                         auto u = M_Xh->element();
        //                         auto v = M_Xh->element();

        //                         u=*T_vec;

        //                 // BC Neumann
        //                         if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "flux" ) )
        //                         {
        //                             for ( auto& [bc, value] : specs["/BoundaryConditions/heat/flux"_json_pointer].items() )
        //                             {
        //                                 LOG( INFO ) << fmt::format( "flux {}: {}", bc, value.dump() );
        //                                 auto flux = expr(value["expr"].get<std::string>());

        //                                 lt += integrate( _range = markedfaces( M_mesh, bc ),
        //                                         _expr = flux  * id( v ) );
        //                             }
        //                         }

        //                 // For each material, integrate rho*C*idt(u)*id(v)/dt + \int_MAT k_MAT * grad(u) * grad(v)
        //                         for ( auto [key, material] : specs["/Models/heat/materials"_json_pointer].items() )
        //                         {
        //                             LOG( INFO ) << fmt::format( "material {}", material );
        //                             std::string mat = fmt::format( "/Materials/{}/k", material.get<std::string>() );
        //                             auto k = specs[nl::json::json_pointer( mat )].get<std::string>();
        //                             std::string matRho = fmt::format( "/Materials/{}/rho", material.get<std::string>() );
        //                             std::string matCp = fmt::format( "/Materials/{}/Cp", material.get<std::string>() );
        //                             auto Rho = specs[nl::json::json_pointer( matRho )].get<std::string>();
        //                             auto Cp = specs[nl::json::json_pointer( matCp )].get<std::string>();

        //                             lt += integrate( _range = markedelements( M_mesh, material.get<std::string>() ), 
        //                                     _expr =  inner( expr( k ) * gradv( u ) , grad( v ) ) );
        //                             lt += integrate(_range = markedelements( M_mesh, material.get<std::string>() ), 
        //                                     _expr = M_bdf->polyDerivCoefficient( 0 ) *expr( Rho ) * expr( Cp ) * idv( u ) * id( v ) );
        //                         }

        //                 // Update of RHS for heat equations starting from T with the term rho*C*T_{n-1}/dt (at order 1)
        //                         for ( auto [key, material] : specs["/Models/heat/materials"_json_pointer].items() )
        //                         {
        //                             std::string matRho = fmt::format( "/Materials/{}/rho", material.get<std::string>() );
        //                             std::string matCp = fmt::format( "/Materials/{}/Cp", material.get<std::string>() );
        //                             auto Rho = specs[nl::json::json_pointer( matRho )].get<std::string>();
        //                             auto Cp = specs[nl::json::json_pointer( matCp )].get<std::string>();

        //                             lt += integrate( _range = markedelements( M_mesh, material.get<std::string>() ),                                                     
        //                                             _expr = - expr( Rho ) * expr( Cp ) * idv(M_bdf->polyDeriv())  * id( v ) );                          
        //                         }
                                
                                
        //                 // Integrating directly the radiative flux obtained via the solution of the problem coming from Modest book
        //                 // Radiative enclosure
        //                         if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
        //                         {
        //                             for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
        //                             {                
        //                                 auto sigma = expr(value["sigma"].get<std::string>());

        //                                 sigma.setParameterValues({{"sigma",specs["/Parameters/sigma"_json_pointer].get<double>()}});
                                        
        //                                 // Recover the emissivity epsilon of the coating material associated to the marked face
        //                                 for ( std::string mark :  M_markers_map[bc] )
        //                                 {
        //                                     auto done = 0;
        //                                     for ( auto& [key, coat] : specs["/Coating"_json_pointer].items() )
        //                                     {
        //                                         for ( auto markcoat : coat.at("markers") )
        //                                         {
        //                                             if ( mark == markcoat )
        //                                             {
        //                                                 auto epsilon = coat["epsilon"].get<std::string>();                            
        //                                                 auto idvu4 = idv( u ) * idv( u ) * idv( u )* idv( u ) ;                                                            
        //                                                 // Term sigma * epsilon * T^4 
        //                                                 lt += integrate( _range = markedfaces( M_mesh, mark ),
        //                                                         _expr = sigma * expr( epsilon ) * idvu4 *id( v ));
        //                                                 // Term - sigma * epsilon \int T^4 dFij 

        //                                                 auto cavity_markers = M_markers_map[bc];   

        //                                                 int i_marker_to_bb=0;
        //                                                 int i_mark=0; 
        //                                                 int j_mark=0; 
        //                                                 auto vf_field = M_Xhd0_map[bc]->element();
        //                                                 for ( std::string mark2 : cavity_markers ) // Loop on index j
        //                                                 {
        //                                                     // Compute the index relative to marker mark
        //                                                     auto it = std::find(cavity_markers.begin(),cavity_markers.end(),mark);                    
        //                                                     i_mark = it - cavity_markers.begin();
        //                                                     // Compute the index relative to marker mark2
        //                                                     auto jt = std::find(cavity_markers.begin(),cavity_markers.end(),mark2);                    
        //                                                     j_mark = jt - cavity_markers.begin();
        //                                                     //std::cout << fmt::format("we add a contribution from index {} whose view factor is {}",j_mark,M_vf_matrix(i_mark,j_mark))<< std::endl;
                                                            
        //                                                     vf_field.zero();                                                       
        //                                                     //std::cout << fmt::format("we add a contribution from marker {} whose emissivity is {}",markcoat,epsilon_mark2)<< std::endl;
        //                                                     vf_field.on(_range=markedfaces(M_mesh,mark2),_expr=cst(M_matrix_vf_map[bc](i_mark,j_mark)));
                                                            
        //                                                     auto T4residual_scal_vf = integrate(_range=markedfaces(M_mesh,mark2), _expr=inner(sigma * expr( epsilon ) * idvu4,idv(vf_field))).evaluate()(0,0);

        //                                                     auto measure_mark2 = integrate( _range=markedfaces(M_mesh,mark2),_expr= cst(1.) ).evaluate()(0,0);
                                                            
        //                                                     lt += integrate( _range=markedfaces(M_mesh,mark),
        //                                                                     _expr= cst(-1)/cst(measure_mark2)*cst(T4residual_scal_vf) *id(v) );  

        //                                                     auto done2=0;
        //                                                     for ( auto& [key2, coat2] : specs["/Coating"_json_pointer].items() )
        //                                                     {
        //                                                         for ( auto markcoat2 : coat2.at("markers") )
        //                                                         {
        //                                                             if ( mark2 == markcoat2 )
        //                                                             {
        //                                                                 auto epsilon_mark2 = coat2["epsilon"].get<std::string>(); 

        //                                                                 lt += integrate( _range=markedfaces(M_mesh,mark),
        //                                                                                 _expr = expr(epsilon) * (cst(1. ) /expr(epsilon_mark2)-cst(1.))/cst(measure_mark2) * cst(400) * gradv(u) *N() * id(v) );
        //                                                             done2 = 1;
        //                                                             break;
        //                                                             }
        //                                                         }
        //                                                         if ( done ){ break; }
        //                                                     }
                                                                                
        //                                                 }

        //                                                 // If the cavity is open, and it is assumed that a black body of fixed temperature
        //                                                 // T_ref exchances heat with the cavity, and additional term is added: 
        //                                                 // - the view factor is computed using the reciprocity formula FijAi = FjiAj
        //                                                 // - the contribution is of the form \sigma T^4 * Fij
        //                                                 if(value["enclosure"]=="open")
        //                                                 {
        //                                                     //std::cout << "Open enclosure " << bc << std::endl;
        //                                                     double vf_marker_to_bb = 1 - M_matrix_vf_map[bc].row(i_marker_to_bb).sum(); //Fij
        //                                                     double vf_bb_to_marker = 1 - M_matrix_vf_map[bc].col(i_marker_to_bb).sum(); //Fji
        //                                                     auto measure_mark1 = integrate( _range=markedfaces(M_mesh,mark),_expr= cst(1.) ).evaluate()(0,0); //Ai
        //                                                     double area_bb = measure_mark1 * vf_marker_to_bb / vf_bb_to_marker;
        //                                                     auto T_bb = value["Tref"].get<double>();
        //                                                     auto T_bb4 = cst(T_bb)*cst(T_bb)*cst(T_bb)*cst(T_bb);

                                                            
        //                                                     lt += integrate( _range = markedfaces( M_mesh, mark ),
        //                                                         _expr =  sigma * expr( epsilon ) * T_bb4 * cst(vf_bb_to_marker) * id( v ));

        //                                                 }
                                                        
        //                                                 i_marker_to_bb++;

        //                                                 done = 1;
        //                                                 break;
        //                                             }
        //                                         }
        //                                         if ( done ){ break; }
        //                                     }    
        //                                 }
        //                             }
        //                         }

        //                 // BC Robin -k*grad(u)*N = h*(T-Text)
        //                         // if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "convective_heat_flux" ) )
        //                         // {
        //                         //     for ( auto& [bc, value] : specs["/BoundaryConditions/heat/convective_heat_flux"_json_pointer].items() )
        //                         //     {
        //                         //         LOG( INFO ) << fmt::format( "convective_heat_flux {}: {}", bc, value.dump() );
        //                         //         auto h = value["h"].get<std::string>();
        //                         //         auto Text = expr(value["Text"].get<std::string>());

        //                         //         Text.setParameterValues({{"Tref_C",specs["/Parameters/Tref_C"_json_pointer].get<double>()}});

        //                         //         lt += integrate( _range = markedfaces( M_mesh, bc ),
        //                         //                 _expr = expr( h ) * id( v ) * idt( u ) );
        //                         //         lt += integrate( _range = markedfaces( M_mesh, bc ),
        //                         //                 _expr = - expr( h ) * Text * id( v ) );
        //                         //     }
        //                         // }

        //                 // Solar radiation
        //                         if(specs["/BoundaryConditions/heat"_json_pointer].contains( "solar_radiation" ) )
        //                         {
        //                             for ( auto& [meteo_station_name, value] : specs["/BoundaryConditions/heat/solar_radiation"_json_pointer].items() )
        //                             {
                                        
        //                                 auto current_Tsky = accessTsky(meteo_station_name, M_bdf->time(), M_bdf->timeStep());
        //                                 current_Tsky += 273.15; // transform in Kelvin
                                        
        //                                 // std::cout <<  fmt::format("Current Tsky residual  {} meteo station {}", current_Tsky, meteo_station_name);

        //                                 auto buildings_specs_string = fmt::format("/BoundaryConditions/heat/solar_radiation/{}/buildings",meteo_station_name);
                                        
        //                                 for(auto const& [building_name,building_structure] : specs[json::json_pointer(buildings_specs_string)].items() )
        //                                 {
        //                                     for(std::string surface_name : building_structure.at("surfaces"))
        //                                     {
                                                
        //                                         auto current_solarRad = accessSolarRadiation(meteo_station_name, building_name, surface_name, M_bdf->time(), M_bdf->timeStep());
        //                                         std::cout <<  fmt::format("Current solar radiation {} surface {}", current_solarRad,surface_name);
        //                                         // Add solar radiation as Neumann boundary condition on the exposed surfaces
        //                                         lt += integrate( _range = markedfaces( M_mesh, surface_name ),
        //                                                         _expr = - cst(current_solarRad) * id( v ));
                                                
        //                                         // Add black-body radiation from sky temperature on both 
        //                                         auto sigma = expr(specs["/Parameters/sigma"_json_pointer].get<double>());   
                                                
        //                                         auto Tsky4 = cst(current_Tsky) * cst(current_Tsky) * cst(current_Tsky) * cst(current_Tsky);  

        //                                         auto idvu3 = idv( T ) * idv( T ) * idv( T );

        //                                         // lt += integrate( _range = markedfaces( M_mesh, surface_name ), 
        //                                         //                 _expr = - sigma * (idvu3 * idv( T ) - Tsky4)  * id( v ) );                                                                                                                                            
                                                
        //                                     }
        //                                 }
        //                             }
        //                         }
        //                         lt_vec->close();

        //                         auto temp = M_Xh->element();
        //                         temp = *lt_vec;

        //                         for ( auto& [bc, value] : specs["/BoundaryConditions/heat/temperature"_json_pointer].items())
        //                         {
        //                             auto dirichletBc = expr(value["expr"].get<std::string>());
        //                             temp.on( _range=markedfaces(M_mesh,bc),_expr=cst(0)*dirichletBc );
        //                         }

        //                         *lt_vec = temp;

        //                         // std::cout << fmt::format("Norm lt vec initial {}\n",lt.vector().l2Norm());
        //                     };
        
        // Impose Dirichlet boundary confitions on the initial guess
        for ( auto& [bc, value] : specs["/BoundaryConditions/heat/temperature"_json_pointer].items())
        {
            auto dirichletBc = expr(value["expr"].get<std::string>());
            T_sol->on( _range=markedfaces(M_mesh,bc),_expr=dirichletBc );
        }

        // Solve the non-linear problem
        backend(_name="heatEq" )->nlSolver()->jacobian = update_jacobian;
        backend(_name="heatEq" )->nlSolver()->residual = update_residual;
        backend(_name="heatEq")->nlSolve( _solution=T_sol,_jacobian=Jac,_residual=Res );

        e_opInterp2->add("Tsol",idv(T_sol));
        e_opInterp2->save();

        M_currentTempAndFlux.setT(T_sol);
        // Update the solution in the BDF data structure
        auto newT = unwrap_ptr(M_currentTempAndFlux.T());
        M_bdf->next(newT);    

    } // end RHT<Dim,Order>::solveHeatEquationNonLinear()

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
                else if (struc["quantity"]=="flux-from-temperature")
                {
                    auto T_val = idv(M_currentTempAndFlux.T());
                    auto T4 = T_val * T_val * T_val * T_val;
                    auto average_q = integrate(_range=markedfaces(M_mesh,markers),_expr=  cst(-400) * gradv(M_currentTempAndFlux.T()) * N() ).evaluate()(0,0);
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
        // rht.execute();
        // std::cout << "Execute non linear code " << std::endl;
        rht.executeNonLinear();

    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
} // end main()