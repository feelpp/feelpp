#include <optional>
#include <unordered_map>
#include <fmt/ostream.h>

#include "rht_readers.hpp"

namespace Feel
{

/// @brief Initialize the time-invariant part of the heat equation
/// 
/// Assembles bilinear and linear forms for the transient heat equation:
///   ρ C ∂T/∂t - ∇·(k∇T) = q_rad + q_flux
/// 
/// This method constructs:
/// - Bilinear form: thermal conductivity and capacitance terms
/// - Linear form: radiative boundary conditions (blackbody emission)
/// 
/// Boundary conditions handled:
/// - Radiative blackbody heat flux: σ ε T_ref⁴ (without view factors)
/// 
/// @note This is called once during initialization; time-dependent and
///       nonlinear terms are added later in solveHeatEquationNonLinear()
    template<int Dim, int Order>
    void RHT<Dim,Order>::initHeatEquation()
    {
        auto a = form2( _test = M_Xh, _trial = M_Xh,_matrix=M_a,_backend=backend(_name="heatEq") );
        auto l = form1( _test = M_Xh,_vector=M_l );

        auto u = M_Xh->element();
        auto v = M_Xh->element();

        M_e=exporter(_mesh=M_mesh,_name="HeatEqn");

        M_bdf->initialize( unwrap_ptr(M_currentTemp.T()) );
        M_bdf->start();

        // For each material, integrate rho*C*idt(u)*id(v)/dt + \int_MAT k_MAT * grad(u) * grad(v)
        for ( auto [key, material] : specs["/Models/heat/materials"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "material {}", material.dump() );
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

            M_conductivity->on( _range = markedelements(M_mesh, material.get<std::string>() ),
                               _expr = expr(k) );
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
                    if (auto epsilon_opt = getCoatingEpsilon(mark)) {
                        auto epsilon = *epsilon_opt;
                        LOG(DEBUG) << fmt::format("Matched marker '{}' with coating emissivity", mark);

                        l += integrate( _range = markedfaces( M_mesh, mark ),
                                _expr =  sigma * expr( epsilon ) * Tref4 * id( v ) );
                    } else {
                        LOG(WARNING) << fmt::format("No coating found for marker: {}", mark);
                    }
                }
            }
        }
    } // end RHT<Dim,Order>::initHeatEquation

/// @brief Export current temperature field for visualization
/// 
/// Writes the temperature solution to the exporter for post-processing
/// in visualization tools (ParaView, etc.)
    template<int Dim, int Order>
    void RHT<Dim,Order>::exportHeat()
    {
        M_e->step(M_bdf->time())->add( "T", idv(M_currentTemp.T()) );
        M_e->save();
    } // end RHT<Dim,Order>::exportHeat()

/// @brief Retrieve coating emissivity for a boundary marker
/// 
/// Searches the JSON specifications for a coating material associated
/// with the given boundary marker and returns its emissivity value.
/// 
/// @param marker Boundary marker name to look up
/// @return Emissivity expression string if found, std::nullopt otherwise
/// 
/// @note This helper eliminates O(n³) nested loops that were previously
///       used for coating lookup in multiple locations
    template<int Dim, int Order>
    std::optional<std::string> RHT<Dim,Order>::getCoatingEpsilon(const std::string& marker) const
    {
        for (const auto& [key, coat] : specs["/Coating"_json_pointer].items()) {
            auto markers_it = coat.find("markers");
            if (markers_it != coat.end()) {
                auto markers = markers_it->get<std::vector<std::string>>();
                if (std::find(markers.begin(), markers.end(), marker) != markers.end()) {
                    return coat["epsilon"].get<std::string>();
                }
            }
        }
        return std::nullopt;
    }

/// @brief Initialize mesh, finite element spaces, and view factor matrices
/// 
/// Performs complete problem setup:
/// 1. Load mesh from JSON specifications
/// 2. Create continuous (M_Xh) and discontinuous (M_Xhd0) finite element spaces
/// 3. Compute or load view factor matrices for radiative enclosures
/// 4. Initialize BDF (Backward Differentiation Formula) time-stepping
/// 5. Set initial temperature field from JSON specifications
/// 
/// View factors F_ij represent the fraction of radiation leaving surface i
/// that is intercepted by surface j, accounting for geometry and occlusion.
/// 
/// @note View factor computation can be expensive; results are cached to file
    template<int Dim, int Order>
    void RHT<Dim,Order>::init()
    {
        // Initialize mesh, spaces, marker list and view factor matrix
        M_mesh = loadMesh( _mesh = new mesh_t, _filename = this->specs["/Meshes/heat/Import/filename"_json_pointer].template get<std::string>() );
        M_Xh = space_t::New(_mesh=M_mesh);
        M_Xhds0 = spacedisc_surf_t::New(_mesh=M_mesh);
        M_conductivity = M_Xhds0->elementPtr();
        tic();
        this->computeVF_and_save();
        toc("Computation of view factors");

        int n_cavities = M_markers_map.size();
        if(n_cavities != 0)
        {
            for(const auto& [cavity_name,markers]: M_markers_map)
            {
                auto cavity_submesh = createSubmesh(_mesh=M_mesh,_range=markedfaces(M_mesh,markers),_update=0);

                // Create a discontinuous space, linear per face, over the view-factor markers
                M_Xhd0 = Pdh<Order-1>(cavity_submesh);
                M_Xhd0_map.insert(std::make_pair(cavity_name,M_Xhd0));
            }
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
        // Initialize the temperature and flux in the data structure needed for the Picard loop
        M_currentTemp.setT(T);

        // Initialize matrix and right-hand sides of the heat transfer problem
        M_lt = backend()->newVector(M_Xh);
        M_l = backend()->newVector(M_Xh);
        auto cg_graph = stencil( _test=M_Xh,_trial=M_Xh)->graph();
        M_a = backend()->newMatrix(0,0,0,0,cg_graph);
        M_at = backend()->newMatrix(0,0,0,0,cg_graph);

        LOG(INFO) << fmt::format("Init routine finished")<< std::endl;

    }  // end RHT<Dim,Order>::init

/// @brief Solve nonlinear heat equation with radiative boundary conditions
/// 
/// Solves the coupled conduction-radiation problem using Newton-Raphson iteration:
///   ρ C ∂T/∂t - ∇·(k∇T) = q_rad(T⁴)
/// 
/// Physics:
/// - Stefan-Boltzmann radiation: q = ε σ T⁴
/// - Net radiosity in enclosure: q_i = ε σ T_i⁴ - ε Σ_j F_ij σ T_j⁴ - Σ_j (1/ε_j - 1) F_ij q_j
/// - Linearization: ∂(T⁴)/∂T ≈ 4T³ (tangent at current iterate)
/// 
/// @param T Initial temperature guess (from BDF time-stepping or previous iterate)
/// 
/// The method constructs Jacobian and residual lambdas that:
/// - Handle view factor integration over cavity surfaces
/// - Apply correct linearization for Newton convergence
/// - Impose Dirichlet boundary conditions
/// 
/// @note This is called at each time step; the nonlinear solver iterates
///       until ||residual|| < tolerance
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


        auto update_jacobian  = [=,this]( const vector_ptrtype& T_vec, sparse_matrix_ptrtype& at_mat ) {

                            auto at = form2(_trial=M_Xh, _test= M_Xh, _matrix=at_mat,_backend=backend(_name="heatEq"));
                            auto a = form2(_trial=M_Xh, _test= M_Xh, _matrix=M_a,_backend=backend(_name="heatEq"));

                            // Start from the time-independent forms computed in initHeatEquation
                            at = a;

                            auto u = M_Xh->element();
                            auto v = M_Xh->element();
                            u=*T_vec;

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
                                        if (auto epsilon_opt = getCoatingEpsilon(mark)) {
                                            auto epsilon = *epsilon_opt;
                                            auto idvu3 = idv( u ) * idv( u ) * idv( u );
                                            // Term 4 sigma * epsilon * T^3 * increment (linearization of d(T⁴)/dT = 4T³)
                                            at += integrate( _range = markedfaces( M_mesh, mark ),
                                                    _expr = cst(STEFAN_BOLTZMANN_DERIVATIVE_COEFF) * sigma * expr( epsilon ) * idvu3 * idt( u ) * id( v ));
                                            // Term - 4 sigma * epsilon \int T^3 dFij * increment

                                            // auto cavity_markers = M_markers_map[bc];

                                            // int i_marker_to_bb=0;

                                            // // If the cavity is open, and it is assumed that a black body of fixed temperature
                                            // // T_ref exchances heat with the cavity, and additional term is added:
                                            // // - the view factor is computed using the reciprocity formula FijAi = FjiAj
                                            // // - the contribution is of the form \sigma T^4 * Fij
                                            // if(value["enclosure"]=="open")
                                            // {
                                            //     //std::cout << "Open enclosure " << bc << std::endl;
                                            //     double vf_marker_to_bb = 1 - M_matrix_vf_map[bc].row(i_marker_to_bb).sum(); //Fij
                                            //     double vf_bb_to_marker = 1 - M_matrix_vf_map[bc].col(i_marker_to_bb).sum(); //Fji
                                            //     auto measure_mark1 = integrate( _range=markedfaces(M_mesh,mark),_expr= cst(1.) ).evaluate()(0,0); //Ai
                                            //     double area_bb = measure_mark1 * vf_marker_to_bb / vf_bb_to_marker;
                                            //     auto T_bb = value["Tref"].get<double>();
                                            //     auto T_bb3 = cst(T_bb)*cst(T_bb)*cst(T_bb);


                                            //     at += integrate( _range = markedfaces( M_mesh, mark ),
                                            //         _expr = cst(4) *  sigma * expr( epsilon ) * T_bb3 * cst(vf_bb_to_marker) * idt( u ) * id( v ));

                                            // }

                                            // i_marker_to_bb++;
                                        } else {
                                            LOG(WARNING) << fmt::format("No coating found for marker: {}", mark);
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

        auto update_residual = [=,this]( const vector_ptrtype& T_vec, vector_ptrtype& lt_vec ) {

                            auto lt = form1(_test=M_Xh, _vector=lt_vec);
                            auto l = form1(_test=M_Xh, _vector=M_l);

                            // Start from the time-independent forms computed in initHeatEquation
                            lt.zero();
                            lt += l;

                            auto u = M_Xh->element();
                            auto v = M_Xh->element();

                            u=*T_vec;

                            // For each material, integrate rho*C*d/dt idt(u)*id(v + \int_MAT k_MAT * grad(u) * grad(v)
                                for ( auto [key, material] : specs["/Models/heat/materials"_json_pointer].items() )
                                {
                                    LOG( INFO ) << fmt::format( "material {}", material.dump() );
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
                                            if (auto epsilon_opt = getCoatingEpsilon(mark)) {
                                                auto epsilon = *epsilon_opt;
                                                auto idvu4 = idv( u ) * idv( u ) * idv( u )* idv( u ) ;
                                                // Term sigma * epsilon * T^4
                                                lt += integrate( _range = markedfaces( M_mesh, mark ),
                                                                _expr = sigma * expr( epsilon ) * idvu4 *id( v ));

                                                        auto cavity_markers = M_markers_map[bc];

                                                        // Build marker index map once for O(1) lookups
                                                        std::unordered_map<std::string, size_t> marker_to_index;
                                                        for (size_t idx = 0; idx < cavity_markers.size(); ++idx) {
                                                            marker_to_index[cavity_markers[idx]] = idx;
                                                        }

                                                        int i_marker_to_bb = 0;
                                                        size_t i_mark = marker_to_index[mark];
                                                        auto vf_field = M_Xhd0_map[bc]->element();
                                                        
                                                        for ( std::string mark2 : cavity_markers ) // Loop on index j
                                                        {
                                                            size_t j_mark = marker_to_index[mark2];

                                                            vf_field.zero();

                                                            vf_field.on(_range=markedfaces(M_mesh,mark2),_expr=cst(M_matrix_vf_map[bc](j_mark,i_mark)));

                                                            auto T4residual_scal_vf = integrate(_range=markedfaces(M_mesh,mark2), _expr=inner(sigma * expr( epsilon ) * idvu4,idv(vf_field))).evaluate()(0,0);

                                                            auto measure_mark1 = integrate( _range=markedfaces(M_mesh,mark),_expr= cst(1.) ).evaluate()(0,0);

                                                            // Term - sigma * epsilon \int T^4 dFij (negative for outgoing radiation)

                                                            lt += integrate( _range=markedfaces(M_mesh,mark),
                                                                            _expr= cst(-1.0) * cst(T4residual_scal_vf) / cst(measure_mark1) * id(v) );

                                                            if (auto epsilon_mark2_opt = getCoatingEpsilon(mark2)) {
                                                                auto epsilon_mark2 = *epsilon_mark2_opt;

                                                                auto flux_other_surface = integrate(_range=markedfaces(M_mesh,mark2), _expr= idv(M_conductivity) * gradv(u) * idv(vf_field) * N() ).evaluate()(0,0);

                                                                // Term - \int sigma * (1/epsilon-1) * q dFij (reflection term for non-black surfaces)

                                                                lt += integrate( _range=markedfaces(M_mesh,mark),
                                                                                _expr = cst(-1.0) * expr(epsilon) * (cst(1.0) / expr(epsilon_mark2) - cst(1.0)) / cst(measure_mark1) * cst(flux_other_surface) * id(v) );
                                                            } else {
                                                                LOG(WARNING) << fmt::format("No coating found for marker: {}", mark2);
                                                            }

                                                        }

                                                        // If the cavity is open, and it is assumed that a black body of fixed temperature
                                                        // T_ref exchanges heat with the cavity, an additional term is added:
                                                        // - the view factor is computed using the reciprocity formula FijAi = FjiAj
                                                        // - the contribution is of the form \sigma T^4 * Fij
                                                        if(value["enclosure"]=="open")
                                                        {
                                                            double vf_marker_to_bb = 1.0 - M_matrix_vf_map[bc].row(i_marker_to_bb).sum(); //Fij
                                                            double vf_bb_to_marker = 1.0 - M_matrix_vf_map[bc].col(i_marker_to_bb).sum(); //Fji
                                                            auto measure_mark1 = integrate( _range=markedfaces(M_mesh,mark),_expr= cst(1.0) ).evaluate()(0,0); //Ai
                                                            double area_bb = measure_mark1 * vf_marker_to_bb / vf_bb_to_marker;
                                                            auto T_bb = value["Tref"].get<double>();
                                                            auto T_bb4 = cst(T_bb) * cst(T_bb) * cst(T_bb) * cst(T_bb);

                                                            // Negative sign for incoming radiation from black body
                                                            lt += integrate( _range = markedfaces( M_mesh, mark ),
                                                                _expr = cst(-1.0) * sigma * expr( epsilon ) * T_bb4 * cst(vf_bb_to_marker) / cst(measure_mark1) * id( v ));
                                                        }

                                                        i_marker_to_bb++;
                                                    } else {
                                                        LOG(WARNING) << fmt::format("No coating found for marker: {}", mark);
                                                    }
                                                } // end for mark
                                            } // end for bc radiative_enclosure_heat_flux
                                        } // end if radiative_enclosure_heat_flux

                            auto temp = M_Xh->element();
                            temp = *lt_vec;

                            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/temperature"_json_pointer].items())
                            {
                                auto dirichletBc = expr(value["expr"].get<std::string>());
                                temp.on( _range=markedfaces(M_mesh,bc),_expr=cst(0)*dirichletBc );
                            }

                            *lt_vec = temp;

        };

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

        M_currentTemp.setT(T_sol);
        // Update the solution in the BDF data structure
        auto newT = unwrap_ptr(M_currentTemp.T());
        M_bdf->next(newT);

    } // end RHT<Dim,Order>::solveHeatEquationNonLinear()

/// @brief Execute time-dependent radiative heat transfer simulation
/// 
/// Main driver routine that:
/// 1. Initializes mesh, spaces, and view factors via init()
/// 2. Assembles time-invariant forms via initHeatEquation()
/// 3. Time-steps using BDF (Backward Differentiation Formula)
/// 4. Solves nonlinear problem at each time step
/// 5. Exports results for visualization
/// 
/// Time-stepping loop continues until BDF indicates completion
/// (steady-state or final time reached).
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
            LOG(INFO) << fmt::format("Solving heat equation at time t = {:.6e}", M_bdf->time());

            // Solve the heat transfer problem at each time instant by solving a Picard loop to ensure convergence
            // of temperature and fluxes
            tic();
            this->solveHeatEquationNonLinear(M_bdf->unknowns()[0]);
            toc("solve non linear heat equation");
            // Export the temperature, and flux on the radiating surfaces
            if(!boption(_name="deactivate-exporters"))
                this->exportHeat();
        }
    } // end RHT<Dim,Order>::executeNonLinear()

/// @brief Validate simulation results against reference values
/// 
/// Compares computed quantities with expected values from literature
/// or analytical solutions specified in JSON configuration.
/// 
/// Supported checks:
/// - Average temperature over specified boundary markers
/// - Average heat flux computed from temperature gradient
/// 
/// Uses Boost.Test CHECK() macro to report pass/fail with relative
/// tolerance defined in JSON specifications.
/// 
/// @note This is the test validation step; failures indicate either
///       numerical issues or incorrect problem setup
    template<int Dim, int Order>
    void RHT<Dim,Order>::checkResults()
    {
        for(auto [key,struc] : specs["/Checker"_json_pointer].items())
        {
            LOG(INFO) << fmt::format("Checker key: {}", key);
            if(struc["type"].get<std::string>()=="average")
            {
                auto markers = struc["markers"].get<std::vector<std::string>>();
                auto value = struc["exact_value"].get<double>();
                auto tol = struc["rel_tolerance"].get<double>();

                if (struc["quantity"]=="temperature")
                {
                    // average temperature on surface
                    auto average_T = integrate(_range=markedfaces(M_mesh,markers),_expr=idv(M_currentTemp.T())).evaluate()(0,0);
                    auto measure_markers = integrate(_range=markedfaces(M_mesh,markers),_expr=cst(1.)).evaluate()(0,0);
                    average_T /=measure_markers;

                    auto rel_difference = (average_T-value)/value;

                    CHECK( math::abs(rel_difference) < tol );

                    LOG(INFO) << fmt::format("Average temperature: {:.6e}", average_T);
                    LOG(INFO) << fmt::format("Markers: [{}]", fmt::join(markers, ", "));
                    LOG(INFO) << fmt::format("Expected value: {:.6e}", value);
                    LOG(INFO) << fmt::format("Tolerance: {:.6e}", tol);
                    LOG(INFO) << fmt::format("Relative difference: {:.6e}", rel_difference);

                }
                else if (struc["quantity"]=="flux-from-temperature")
                {
                    auto T_val = idv(M_currentTemp.T());
                    auto T4 = T_val * T_val * T_val * T_val;
                    // Negative sign: heat flux convention (outward normal)
                    auto average_q = integrate(_range=markedfaces(M_mesh,markers),_expr= cst(-1.0) * idv(M_conductivity) * gradv(M_currentTemp.T()) * N() ).evaluate()(0,0);
                    auto measure_markers = integrate(_range=markedfaces(M_mesh,markers),_expr=cst(1.)).evaluate()(0,0);
                    average_q /=measure_markers;

                    auto rel_difference = (average_q-value)/value;

                    CHECK( math::abs(rel_difference) < tol );

                    LOG(INFO) << fmt::format("Average flux: {:.6e}", average_q);
                    LOG(INFO) << fmt::format("Markers: [{}]", fmt::join(markers, ", "));
                    LOG(INFO) << fmt::format("Expected value: {:.6e}", value);
                    LOG(INFO) << fmt::format("Tolerance: {:.6e}", tol);
                    LOG(INFO) << fmt::format("Relative difference: {:.6e}", rel_difference);
                }
                else
                {
                    LOG(WARNING) << fmt::format("Quantity '{}' not supported for average checker", struc["quantity"].get<std::string>());
                }
            }
            else
            {
                LOG(WARNING) << fmt::format("Checker type '{}' not supported (only 'average' is implemented)", struc["type"].get<std::string>());
            }
        }

    }


} // namespace Feel
