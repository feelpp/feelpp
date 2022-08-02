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
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "rht options" );
    options.add_options()

        // mesh parameters
        ( "specs", Feel::po::value<std::string>(),
          "json spec file for rht" );

        ( "steady", Feel::po::value<std::bool>()->defaut_value( 1 ),
          "if 1: steady else unsteady" );

    return options.add( Feel::feel_options() );
}

template <int Dim = 2, int Order = 1>
int runHeat( nl::json const& specs )
{
    auto mesh = loadMesh( _mesh = new Mesh<Simplex<Dim>>, _filename = specs["/Meshes/heat/Import/filename"_json_pointer].get<std::string>() );
    auto Xh = Pch<Order>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();

    auto a = form2( _test = Xh, _trial = Xh );
    auto at = form2( _test = Xh, _trial = Xh );
    auto l = form1( _test = Xh );
    auto lt = form1( _test = Xh );

    auto M_bdf = bdf( _space = Xh );

    M_bdf->start();

    // from now if the option "steady" is set to True then M_ bdf-setSteady will set time-step=time-final
    if ( boption("steady") )
        M_bdf->setSteady();

    for ( auto [key, material] : specs["/Models/heat/materials"_json_pointer].items() )
    {
        LOG( INFO ) << fmt::format( "material {}", material );
        std::string mat = fmt::format( "/Materials/{}/k", material.get<std::string>() );
        auto k = specs[nl::json::json_pointer( mat )].get<std::string>();

        a += integrate( _range = markedelements( mesh, material.get<std::string>() ), 
                _expr = M_bdf->polyDerivCoefficient( 0 )  * gradt( v ) * trans( expr( k ) * grad( u ) ) );
    }

    // BC Neumann
    if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "flux" ) )
    {
        for ( auto& [bc, value] : specs["/BoundaryConditions/heat/flux"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "flux {}: {}", bc, value.dump() );
            auto flux = value["expr"].get<std::string>();

            l += integrate( _range = markedfaces( mesh, bc ),
                    _expr = M_bdf->polyDerivCoefficient( 0 ) * expr( flux ) * id( v ) );
        }
    }

    // BC Robin
    if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "convective_heat_flux" ) )
    {
        for ( auto& [bc, value] : specs["/BoundaryConditions/heat/convective_heat_flux"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "convective_heat_flux {}: {}", bc, value.dump() );
            auto h = value["h"].get<std::string>();
            auto Text = value["Text"].get<std::string>();

            a += integrate( _range = markedfaces( mesh, bc ),
                    _expr = M_bdf->polyDerivCoefficient( 0 ) * expr( h ) * id( v ) * idt( u ) );
            l += integrate( _range = markedfaces( mesh, bc ),
                    _expr = M_bdf->polyDerivCoefficient( 0 ) * expr( h ) * expr( Text ) * id( v ) );
        }
    }

    // ===== BC RHT =====
    // Blackbody radiative condition
    if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_blackbody_heat_flux" ) )
    {
        for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_blackbody_heat_flux"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "radiative_blackbody_heat_flux {}: {}", bc, value.dump() );
            auto sigma = value["sigma"].get<std::string>();
            auto epsilon = value["epsilon"].get<std::string>();
            auto Tref = value["Tref"].get<std::string>();

            auto Tref4 = expr( Tref ) * expr( Tref ) * expr( Tref ) * expr( Tref );

            l += integrate( _range = markedfaces( mesh, bc ), 
                    _expr = M_bdf->polyDerivCoefficient( 0 ) * expr( sigma ) * expr( epsilon ) * Tref4 * id( v ) );
        }
    }

    // Radiative enclosures
    if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
    {
        for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
        {
            // Closed enclosure
            if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/enclosure"_json_pointer] == "close" )
            {
                LOG( INFO ) << fmt::format( "radiative_closed_enclosure_heat_flux {}: {}", bc, value.dump() );

                if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}"_json_pointer].contains( "sigma" ) )
                {
                    auto sigma = value["sigma"].get<std::string>();
                }
                else
                {
                    // Message: default value of sigma will be used, sigma = 5.67e-8 W.m^-2.K^-4
                    auto sigma = 5.67e-8;
                }
                
                for (  )
                {
                    // Load OR compute view factors
                    if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}"_json_pointer].contains( "viewfactors" ) )
                    {
                        if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors"_json_pointer].contains( "status" ) )
                        {
                            if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors/status"_json_pointer] == "load" )
                            {
                                if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors"_json_pointer].contains( "filename" ) )
                                {
                                    // file where view factors are saved
                                    // Message: view factors are loading from filename
                                    goto loadVF;
                                }
                                else
                                {
                                    // return error: need filename
                                }
                            }
                            if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors/status"_json_pointer] == "compute" )
                            {
                                if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors"_json_pointer].contains( "filename" ) )
                                {
                                    // file where view factors will be saved
                                    // Message: view factors will be saved in file
                                }
                                else
                                {
                                    // filename == None
                                    // Message: view factors will NOT be saved
                                }
                                if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors"_json_pointer].contains( "type" ) )
                                {
                                    if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors/type"_json_pointer].contains( "raytracing" ) )
                                    {
                                        // Message: compute using MC method
                                        goto compute_MC;
                                    }
                                    else if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors/type"_json_pointer].contains( "unobstructed" ) )
                                    {
                                        // Message: compute using integral method
                                        goto compute_integral;
                                    }
                                    else
                                    {
                                        // Message default case: compute using MC method
                                        goto compute_MC;
                                    }
                                }
                                else
                                {
                                    // Message default case: compute using MC method
                                    goto compute_MC;
                                }
                            }
                        }
                        else
                        {
                            if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors"_json_pointer].contains( "filename" ) )
                            {
                                // file where view factors will be saved
                                // Message: view factors will be saved in file
                            }
                            else
                            {
                                // filename == None
                                // Message: view factors will NOT be saved
                            }
                            // Message default case: compute using MC method
                            goto compute_MC;
                        }
                    }
                    else
                    {
                        // Message default case: compute using MC method
                        goto compute_MC;
                    }

                    compute_MC:
                        // Compute view factors with MC raytracing method
                    
                    compute_integral:
                        // Compute view factors with integral method

                    loadVF:
                        // Load view factors

                    // Recover the emissivity epsilon of material linked to the marked face
                    auto epsilon = value["epsilon"].get<std::string>();

                    // Build left hand side matrix LHSm
                    // delta_ij/epsilon_j - (1/epsilon_j -1)*F_ij

                    // Build right hand side matrix RHSm

                    // Compute inverse matrix LHSm : invLHSm, and M1 = invLHSm * RHSm

                }
            }
            // Opened enclosure
            else if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/enclosure"_json_pointer] == "open" )
            {
                LOG( INFO ) << fmt::format( "radiative_opened_enclosure_heat_flux {}: {}", bc, value.dump() );
                auto sigma = value["sigma"].get<std::string>();
                auto Tref = value["Tref"].get<std::string>();

                for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
                {

                    if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors"_json_pointer].contains( "status" ) )
                        {
                            if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors/status"_json_pointer] == "load" )
                            {
                                if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors"_json_pointer].contains( "filename" ) )
                                {
                                    // file where view factors are saved
                                    // Message: view factors are loading from filename
                                    goto loadVF;
                                }
                                else
                                {
                                    // return error: need filename
                                }
                            }
                            if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors/status"_json_pointer] == "compute" )
                            {
                                if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors"_json_pointer].contains( "filename" ) )
                                {
                                    // file where view factors will be saved
                                    // Message: view factors will be saved in file
                                }
                                else
                                {
                                    // filename == None
                                    // Message: view factors will NOT be saved
                                }
                                if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors"_json_pointer].contains( "type" ) )
                                {
                                    if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors/type"_json_pointer].contains( "raytracing" ) )
                                    {
                                        // Message: compute using MC method
                                        goto compute_MC;
                                    }
                                    else if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors/type"_json_pointer].contains( "unobstructed" ) )
                                    {
                                        // Message: compute using integral method
                                        goto compute_integral;
                                    }
                                    else
                                    {
                                        // Message default case: compute using MC method
                                        goto compute_MC;
                                    }
                                }
                                else
                                {
                                    // Message default case: compute using MC method
                                    goto compute_MC;
                                }
                            }
                        }
                        else
                        {
                            if ( specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux/{}/viewfactors"_json_pointer].contains( "filename" ) )
                            {
                                // file where view factors will be saved
                                // Message: view factors will be saved in file
                            }
                            else
                            {
                                // filename == None
                                // Message: view factors will NOT be saved
                            }
                            // Message default case: compute using MC method
                            goto compute_MC;
                        }
                    }
                    else
                    {
                        // Message default case: compute using MC method
                        goto compute_MC;
                    }

                    compute_MC:
                        // Compute view factors with MC raytracing method
                    
                    compute_integral:
                        // Compute view factors with integral method

                    loadVF:
                        // Load view factors
                    
                    // Recover the emissivity epsilon of material linked to the marked face
                    auto epsilon = value["epsilon"].get<std::string>();
                    
                    // Build left hand side matrix LHSm
                    // delta_ij/epsilon_j - (1/epsilon_j -1)F_ij

                    // Build right hand side matrix RHSm
                    // sigma*(delta_ij - F_ij)

                    // Compute inverse matrix LHSm : invLHSm, and M1 = invLHSm * RHSm

                    // Build vecG = f_i(tk) * sigma * Tref**4

                    // Store invLHSm
                    
                }
            }
            else
            {
                // Return error: condition doesn't exist
            }
        }
    }

    M_bdf->initialize( u );

    if ( boption("steady") )
        std::cout << "\n***** Steady state *****" << std::endl;
    else
        std::cout << "The step is  " << M_bdf->timeStep() << "\n"
                  << "The initial time is " << M_bdf->timeInitial() << "\n"
                  << "The final time is " << M_bdf->timeFinal() << "\n"
                  << "BDF order :  " << M_bdf->timeOrder() << "\n" << std::endl;

    // Solve
    for ( M_bdf->start(); M_bdf->isFinished()==false; M_bdf->next(u) )
    {
        lt.zero();
        at.zero();
        at += a;
        lt += l;

        for ( auto [key, material] : specs["/Models/heat/materials"_json_pointer].items() )
        {
            std::string matRho = fmt::format( "/Materials/{}/rho", material.get<std::string>() );
            std::string matCp = fmt::format( "/Materials/{}/Cp", material.get<std::string>() );
            auto Rho = specs[nl::json::json_pointer( matRho )].get<std::string>();
            auto Cp = specs[nl::json::json_pointer( matCp )].get<std::string>();

            lt += integrate( _range = markedelements( mesh, material.get<std::string>() ), 
                    _expr = expr( Rho ) * expr( Cp ) * idv( M_bdf->polyDeriv() ) * id( v ) );
        }

        // Blackbody radiative condition
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_blackbody_heat_flux" ) )
        {
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_blackbody_heat_flux"_json_pointer].items() )
            {
                auto sigma = value["sigma"].get<std::string>();

                // Recover the emissivity epsilon of material linked to the marked face
                auto epsilon = value["epsilon"].get<std::string>();

                auto idvu3 = idv( M_bdf->polyDeriv() ) * idv( M_bdf->polyDeriv() ) * idv( M_bdf->polyDeriv() );

                at += integrate( _range = markedfaces( mesh, bc ),
                        _expr = M_bdf->polyDerivCoefficient( 0 ) * expr( sigma ) * expr( epsilon ) * idv3 * idt( u ) * id( v ) );
            }
        }
        // Radiative enclosure
        if ( specs["/BoundaryConditions/heat"_json_pointer].contains( "radiative_enclosure_heat_flux" ) )
        {
            for ( auto& [bc, value] : specs["/BoundaryConditions/heat/radiative_enclosure_heat_flux"_json_pointer].items() )
            {
                auto sigma = value["sigma"].get<std::string>();
                
                // Recover the emissivity epsilon of material linked to the marked face
                auto epsilon = value["epsilon"].get<std::string>();

                // invLHSm: inverse matrix LHSm
                // M1: invLHSm * RHSm

                // Build F = (u**3(tk)) then compute M1 * F in idtu3
                auto idtu3 = idt( u ) * idt( u ) * idt( u );

                // Add to LHS bilinear form a
                at += integrate( _range = markedfaces( mesh, bc ),
                        _expr = M_bdf->polyDerivCoefficient( 0 ) * expr( sigma ) * expr( epsilon ) * idtu3 * idt( u ) * id( v ) );
            }
        }

        at.solve( _rhs = lt, _solution = u );
    }

    // compute outputs
    auto m = mean( _range = elements( mesh ), _expr = idv( u ) );
    auto m_root = mean( _range = markedfaces( mesh, "Gamma_root" ), _expr = idv( u ) );
    if ( Environment::isMasterRank() )
    {
        std::cout << fmt::format( "- mean value: {}", m ) << std::endl;
        std::cout << fmt::format( "-  min value: {}", u.min() ) << std::endl;
        std::cout << fmt::format( "-  max value: {}", u.max() ) << std::endl;
        std::cout << fmt::format( "-  max deviation: {}", u.max() - u.min() ) << std::endl;
        std::cout << fmt::format( "-  mean root: {}", m_root ) << std::endl;
    }
    // Export
    auto e = exporter( _mesh = mesh );
    e->addRegions();
    e->add( "T", v );
    e->save();

    return 0;
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
        auto jsonfile = removeComments( readFromFile( Environment::expand( soption( "specs" ) ) ) );
        std::istringstream istr( jsonfile );
        json specs = json::parse( istr );
        runHeat<FEELPP_DIM, FEELPP_ORDER>( specs );
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
}