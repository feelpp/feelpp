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
    tempflux_type M_tempFlux;

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
    M_electro = electro_type::New(M_prefixElectro);
    M_thermo = thermo_type::New(M_prefixThermo);
    M_mesh = loadMesh( new mesh_type );
    M_electro->init( M_mesh );
    M_thermo->init( M_mesh );
    toc("init");
}

template<int Dim, int OrderT, int OrderV>
void
ThermoElectricHDG<Dim, OrderT, OrderV>::run()
{
    auto electroMat = M_electro->modelProperties().materials();
    auto thermoMat = M_thermo->modelProperties().materials();

    auto sigmaF = M_electro->potentialSpace()->element();
    auto kF = M_thermo->potentialSpace()->element();

    tic();
    M_electro->assembleCstPart();
    toc("assembleCstElectro");
    tic();
    M_thermo->assembleCstPart();
    toc("assembleCstThermo");

    // first iteration
    M_electro->updateConductivityTerm( false );
    M_electro->assembleRhsBoundaryCond();
    M_electro->assembleRHS();
    M_electro->solve();
    M_potential = M_electro->potentialField();
    M_current = M_electro->fluxField();
    double normP = normL2(elements(M_mesh), idv(M_potential) );
    double normC = normL2(elements(M_mesh), idv(M_current) );
    Feel::cout << "norm potential = " << normP << "\tnorm current = " << normC << std::endl;

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
    M_tempFlux = M_thermo->fluxField();
    double normT = normL2(elements(M_mesh), idv(M_temperature) );
    double normQ = normL2(elements(M_mesh), idv(M_tempFlux) );
    Feel::cout << "norm temperature = " << normT << "\tnorm tempFlux = " << normQ << std::endl;

    for( auto const& pairMat : electroMat )
    {
        std::string marker = pairMat.first;
        auto mat = pairMat.second;
        auto cond = mat.getScalar(soption(prefixvm(M_prefixElectro,"conductivity_json")));
        sigmaF.on( _range=markedelements(M_electro->mesh(), marker), _expr=cond );
    }
    for( auto const& pairMat : thermoMat )
    {
        std::string marker = pairMat.first;
        auto mat = pairMat.second;
        auto cond = mat.getScalar(soption(prefixvm(M_prefixThermo,"conductivity_json")));
        kF.on( _range=markedelements(M_thermo->mesh(), marker), _expr=cond );
    }

    potential_type oldPotential = M_electro->potentialSpace()->element();
    current_type oldCurrent = M_electro->fluxSpace()->element();
    temp_type oldTemperature = M_thermo->potentialSpace()->element();
    tempflux_type oldTempFlux = M_thermo->fluxSpace()->element();
    double incrV = 1;
    double incrT = 1;
    double incrJ = 1;
    double incrQ = 1;

    double tol = doption( prefixvm( M_prefix, "tolerance") );
    int i = 0, itmax = ioption( prefixvm( M_prefix, "itmax") );

    auto eV = exporter( _mesh=M_mesh, _name="electro" );
    auto eT = exporter( _mesh=M_mesh, _name="thermo" );

    eV->step(i)->add("potential", M_potential);
    eV->step(i)->add("current", M_current);
    eV->step(i)->add("sigma", sigmaF);
    eT->step(i)->add("temperature", M_temperature);
    eT->step(i)->add("tempFlux", M_tempFlux);
    eT->step(i)->add("k", kF);

    // Picard loop
    while( i++ < itmax && ( incrV > tol || incrT > tol ) )
    {
        tic();

        oldPotential = M_potential;
        oldCurrent = M_current;
        oldTemperature = M_temperature;
        oldTempFlux = M_tempFlux;

        tic();
        M_electro->setZero();
        M_electro->assembleCstPart();

        M_electro->assembleRHS();
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
            sigmaF.on( _range=markedelements(M_electro->mesh(), marker), _expr=sigma );
        }
        M_electro->assembleRhsBoundaryCond();
        toc("assembleElectro");
        tic();
        M_electro->solve();
        toc("solveElectro");
        M_potential = M_electro->potentialField();
        M_current = M_electro->fluxField();

        tic();
        M_thermo->setZero();
        M_thermo->assembleCstPart();

        for( auto const& pairMat : thermoMat )
        {
            std::string marker = pairMat.first;
            auto mat = pairMat.second;
            double alpha = mat.getDouble("alpha");
            double T0 = mat.getDouble("T0");
            double k0 = mat.getDouble(soption(prefixvm(M_prefixThermo,"conductivity_json")));
            auto k = mat.getScalar(soption(prefixvm(M_prefixThermo,"conductivityNL_json")),
                                   {"T"}, {idv(M_temperature)},
                                   {{"k0",k0},{"T0",T0},{"alpha",alpha}});
            M_thermo->updateConductivityTerm( k, marker);
            kF.on( _range=markedelements(M_thermo->mesh(), marker), _expr=k );
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
        M_tempFlux = M_thermo->fluxField();

        tic();
        eV->step(i)->add("potential", M_potential);
        eV->step(i)->add("current", M_current);
        eV->step(i)->add("sigma", sigmaF);
        eT->step(i)->add("temperature", M_temperature);
        eT->step(i)->add("tempFlux", M_tempFlux);
        eT->step(i)->add("k", kF);
        eV->save();
        eT->save();
        toc("export");

        tic();
        double normV = normL2(elements(M_mesh), idv(oldPotential) );
        double normT = normL2(elements(M_mesh), idv(oldTemperature) );
        double errV = normL2( elements(M_mesh), idv(M_potential) - idv(oldPotential) );
        double errT = normL2( elements(M_mesh), idv(M_temperature) - idv(oldTemperature) );
        double normJ = normL2(elements(M_mesh), idv(oldCurrent) );
        double normQ = normL2(elements(M_mesh), idv(oldTempFlux) );
        double errJ = normL2( elements(M_mesh), idv(M_current) - idv(oldCurrent) );
        double errQ = normL2( elements(M_mesh), idv(M_tempFlux) - idv(oldTempFlux) );
        incrV = errV/normV;
        incrT = errT/normT;
        incrJ = errJ/normJ;
        incrQ = errQ/normQ;

        Feel::cout << "Picard[" << i << "] increment norm V: " << incrV
                   << " increment norm J: " << incrJ << std::endl;
        Feel::cout << "Picard[" << i << "] increment norm T: " << incrT
                   << " increment norm Q: " << incrQ << std::endl;
        toc("increment");

        toc("loop");
    } // Picard loop

}
