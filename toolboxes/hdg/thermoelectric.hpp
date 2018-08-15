#include <feel/feelmodels/hdg/mixedpoisson.hpp>

using namespace Feel;
using namespace Feel::FeelModels;

inline
po::options_description
makeThermoElectricHDGOptions( std::string prefix = "thermoelectric",
                              std::string prefixThermo = "thermo",
                              std::string prefixElectro = "electro" )
{
    po::options_description teOptions( "ThermoElectric HDG options");
    teOptions.add_options()
        ( prefixvm( prefix, "tolerance").c_str(), po::value<double>()->default_value( 1e-5 ),
          "picard tolerance" )
        ( prefixvm( prefix, "itmax").c_str(), po::value<int>()->default_value( 20 ),
          "picard max iteration")
        ;
    teOptions.add( makeMixedPoissonOptions( prefixThermo ) );
    teOptions.add( makeMixedPoissonOptions( prefixElectro ) );
    return teOptions;
}

template<int Dim, int OrderT, int OrderV>
class ThermoElectricHDG
{
public:
    using thermo_type = FeelModels::MixedPoisson<Dim, OrderT>;
    using thermo_ptrtype = std::shared_ptr<thermo_type>;
    using temp_type = typename thermo_type::Wh_element_t;
    using temp_ptrtype = typename thermo_type::Wh_element_ptr_t;
    using tempflux_type = typename thermo_type::Vh_element_t;
    using tempflux_ptrtype = typename thermo_type::Vh_element_ptr_t;

    using electro_type = FeelModels::MixedPoisson<Dim, OrderV>;
    using electro_ptrtype = std::shared_ptr<electro_type>;
    using potential_type = typename electro_type::Wh_element_t;
    using potential_ptrtype = typename electro_type::Wh_element_ptr_t;
    using current_type = typename electro_type::Vh_element_t;
    using current_ptrtype = typename electro_type::Vh_element_ptr_t;

    using mesh_type = typename thermo_type::mesh_type;
    using mesh_ptrtype = typename thermo_type::mesh_ptrtype;

private:
    std::string M_prefix;
    mesh_ptrtype M_mesh;

    std::string M_prefixThermo;
    thermo_ptrtype M_thermo;
    temp_type M_temperature;
    tempflux_type M_tempflux;

    std::string M_prefixElectro;
    electro_ptrtype M_electro;
    potential_type M_potential;
    current_type M_current;

public:
    ThermoElectricHDG( std::string prefix = "thermoelectric",
                       std::string prefixThermo = "thermo",
                       std::string prefixElectro = "electro" );
    void run();
};

template<int Dim, int OrderT, int OrderV>
ThermoElectricHDG<Dim, OrderT, OrderV>::ThermoElectricHDG( std::string prefix,
                                                           std::string prefixThermo,
                                                           std::string prefixElectro )
    :
    M_prefix(prefix),
    M_prefixThermo(prefixThermo),
    M_prefixElectro(prefixElectro)
{
    tic();
    M_thermo = thermo_type::New(M_prefixThermo);
    M_electro = electro_type::New(M_prefixElectro);
    M_mesh = loadMesh( new mesh_type );
    M_thermo->init( M_mesh );
    M_electro->init( M_mesh );
    toc("init");
}

template<int Dim, int OrderT, int OrderV>
void
ThermoElectricHDG<Dim, OrderT, OrderV>::run()
{
    auto electroMat = M_electro->modelProperties().materials();
    auto thermoMat = M_thermo->modelProperties().materials();

    tic();
    M_electro->assembleCstPart();
    toc("assembleCstElectro");
    tic();
    M_thermo->assembleCstPart();
    toc("assembleCstThermo");

    // first iteration
    M_electro->copyCstPart();
    M_electro->updateConductivityTerm( false );
    M_electro->assembleRhsBoundaryCond();
    M_electro->solve();
    M_potential = M_electro->potentialField();
    M_current = M_electro->fluxField();

    M_thermo->copyCstPart();
    M_thermo->updateConductivityTerm( false );
    M_thermo->assembleRhsBoundaryCond();
    for( auto const& pairMat : electroMat )
    {
        std::string marker = pairMat.first;
        auto mat = pairMat.second;
        auto cond = mat.getScalar(soption(prefixvm(M_prefixElectro,"conductivity_json")));
        auto rhs = inner(idv(M_current))/cond;
        M_thermo->assemblePotentialRHS( rhs, marker);
    }
    M_thermo->solve();
    M_temperature = M_thermo->potentialField();
    M_tempflux = M_thermo->fluxField();

    potential_type oldPotential = M_electro->potentialSpace()->element();
    current_type oldCurrent = M_electro->fluxSpace()->element();
    temp_type oldTemperature = M_thermo->potentialSpace()->element();
    tempflux_type oldTempflux = M_thermo->fluxSpace()->element();
    double incrV = normL2( elements(M_mesh), idv(M_potential) - idv(oldPotential) );
    double incrT = normL2( elements(M_mesh), idv(M_temperature) - idv(oldTemperature) );

    double tol = doption( prefixvm( M_prefix, "tolerance") );
    int i = 0, itmax = ioption( prefixvm( M_prefix, "itmax") );

    // Picard loop
    while( i++ < itmax && ( incrV > tol || incrT > tol ) )
    {
        tic();

        oldPotential = M_potential;
        oldCurrent = M_current;
        oldTemperature = M_temperature;
        oldTempflux = M_tempflux;

        tic();
        M_electro->copyCstPart();
        for( auto const& pairMat : electroMat )
        {
            std::string marker = pairMat.first;
            auto mat = pairMat.second;
            double alpha = mat.getDouble("alpha");
            double T0 = mat.getDouble("T0");
            double sigma0 = mat.getDouble(soption(prefixvm(M_prefixElectro,"conductivity_json")));
            auto sigma = mat.getScalar(soption(prefixvm(M_prefixElectro,"conductivityNL_json")),
                                       {"T"}, {idv(M_temperature)},
                                       {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
            M_electro->updateConductivityTerm( sigma, marker);
        }
        M_electro->assembleRhsBoundaryCond();
        toc("assembleElectro");
        tic();
        M_electro->solve();
        toc("solveElectro");
        M_potential = M_electro->potentialField();
        M_current = M_electro->fluxField();

        tic();
        M_thermo->copyCstPart();
        for( auto const& pairMat : thermoMat )
        {
            std::string marker = pairMat.first;
            auto mat = pairMat.second;
            double alpha = mat.getDouble("alpha");
            double T0 = mat.getDouble("T0");
            double k0 = mat.getDouble(soption(prefixvm(M_prefixThermo,"conductivity_json")));
            auto k = mat.getScalar(soption(prefixvm(M_prefixThermo,"conductivity_json")),
                                   {"T"}, {idv(M_temperature)},
                                   {{"k0",k0},{"T0",T0},{"alpha",alpha}});
            M_thermo->updateConductivityTerm( k, marker);
        }
        M_thermo->assembleRhsBoundaryCond();
        for( auto const& pairMat : electroMat )
        {
            std::string marker = pairMat.first;
            auto mat = pairMat.second;
            double alpha = mat.getDouble("alpha");
            double T0 = mat.getDouble("T0");
            double sigma0 = mat.getDouble(soption(prefixvm(M_prefixElectro,"conductivity_json")));
            auto sigma = mat.getScalar(soption(prefixvm(M_prefixElectro,"conductivityNL_json")),
                                       {"T"}, {idv(M_temperature)},
                                       {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
            auto rhs = inner(idv(M_current))/sigma;
            M_thermo->assemblePotentialRHS( rhs, marker);
        }
        toc("assembleThermo");
        tic();
        M_thermo->solve();
        toc("solveThermo");
        M_temperature = M_thermo->potentialField();
        M_tempflux = M_thermo->fluxField();

        incrV = normL2( elements(M_mesh), idv(M_potential) - idv(oldPotential) );
        incrT = normL2( elements(M_mesh), idv(M_temperature) - idv(oldTemperature) );

        toc("loop");
    } // Picard loop

    auto e = exporter( M_mesh );
    e->add("potential", M_potential);
    e->add("current", M_current);
    e->add("temperature", M_temperature);
    e->add("tempflux", M_tempflux);
    e->save();
}
