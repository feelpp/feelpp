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
        ( prefixvm( prefix, "input-marker").c_str(), po::value<std::string>()->default_value("V1"),
          "marker for input (where V=0)")
        ( prefixvm( prefix, "output-marker").c_str(), po::value<std::string>()->default_value("V0"),
          "marker for output (where V!=0)")
        ( prefixvm( prefix, "continuation").c_str(), po::value<bool>()->default_value(false),
          "use continuation for Picard" )
        ( prefixvm( prefix, "continuation.steps").c_str(), po::value<int>()->default_value(1),
          "number of steps to use for continuation for Picard" )
        ( prefixvm( prefix, "load-initial-guess").c_str(), po::value<bool>()->default_value(false), "load initial guess" )
        ( prefixvm( prefix, "load-initial-guess.path").c_str(), po::value<std::string>()->default_value(""), "path to load initial guess" )
        ;
    teOptions.add( makeMixedPoissonOptions( prefixThermo ) );
    teOptions.add( makeMixedPoissonOptions( prefixElectro ) );
    return teOptions;
}

template<int Dim, int OrderT, int OrderV, int OrderG>
class ThermoElectricHDG
{
public:
    using thermo_type = FeelModels::MixedPoisson<Dim, OrderT, OrderG>;
    using thermo_ptrtype = std::shared_ptr<thermo_type>;
    using temp_type = typename thermo_type::Wh_element_t;
    using temp_ptrtype = typename thermo_type::Wh_element_ptr_t;
    using tempflux_type = typename thermo_type::Vh_element_t;
    using tempflux_ptrtype = typename thermo_type::Vh_element_ptr_t;

    using electro_type = FeelModels::MixedPoisson<Dim, OrderV, OrderG>;
    using electro_ptrtype = std::shared_ptr<electro_type>;
    using potential_type = typename electro_type::Wh_element_t;
    using potential_ptrtype = typename electro_type::Wh_element_ptr_t;
    using current_type = typename electro_type::Vh_element_t;
    using current_ptrtype = typename electro_type::Vh_element_ptr_t;

    using property_type = typename electro_type::Wh_element_t;

    using mesh_type = typename thermo_type::mesh_type;
    using mesh_ptrtype = typename thermo_type::mesh_ptrtype;

    using init_guess_space_type = FunctionSpace<mesh_type, bases<Lagrange<OrderT, Scalar, Continuous, PointSetFekete> > >;
    using init_guess_space_ptrtype = std::shared_ptr<init_guess_space_type>;
    using init_guess_type = typename init_guess_space_type::element_type;

private:
    std::string M_prefix;
    mesh_ptrtype M_mesh;

    std::string M_prefixThermo;
    thermo_ptrtype M_thermo;
    temp_type M_temperature;
    tempflux_type M_tempflux;
    std::string M_kKey;
    std::string M_kKeyNL;
    std::string M_temperatureKey;
    std::string M_heatFluxKey;

    std::string M_prefixElectro;
    electro_ptrtype M_electro;
    potential_type M_potential;
    current_type M_current;
    std::string M_sigmaKey;
    std::string M_sigmaKeyNL;
    std::string M_potentialKey;
    std::string M_currentKey;

    init_guess_space_ptrtype M_initGuessSpace;
    init_guess_type M_initGuess;

    property_type M_joule;
    property_type M_sigma;
    property_type M_k;


    int M_itMax;
    double M_tolerance;
    bool M_useContinuation;
    int M_continuationSteps;
    std::string M_inputMarker;
    std::string M_outputMarker;

public:
    ThermoElectricHDG( std::string prefix = "thermoelectric",
                       std::string prefixThermo = "thermo",
                       std::string prefixElectro = "electro" );
    void init( mesh_ptrtype const& mesh = nullptr );
    void solve();
    void exportResults();

    mesh_ptrtype mesh() const { return M_mesh; }
    thermo_ptrtype thermo() const { return M_thermo; }
    electro_ptrtype electro() const { return M_electro; }
    ModelProperties modelProperties() const { return M_electro->modelProperties(); }
};

template<int Dim, int OrderT, int OrderV, int OrderG>
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::ThermoElectricHDG( std::string prefix,
                                                           std::string prefixThermo,
                                                           std::string prefixElectro )
    :
    M_prefix(prefix),
    M_prefixThermo("hdg.poisson."+prefixThermo),
    M_kKey(soption(prefixvm(M_prefixThermo,"conductivity_json"))),
    M_kKeyNL(soption(prefixvm(M_prefixThermo,"conductivityNL_json"))),
    M_prefixElectro("hdg.poisson."+prefixElectro),
    M_sigmaKey(soption(prefixvm(M_prefixElectro,"conductivity_json"))),
    M_sigmaKeyNL(soption(prefixvm(M_prefixElectro,"conductivityNL_json"))),
    M_itMax(ioption(prefixvm( M_prefix, "itmax"))),
    M_tolerance(doption(prefixvm( M_prefix, "tolerance"))),
    M_useContinuation(boption(prefixvm( M_prefix, "continuation"))),
    M_continuationSteps(ioption(prefixvm( M_prefix, "continuation.steps"))),
    M_inputMarker(soption(prefixvm( M_prefix, "input-marker"))),
    M_outputMarker(soption(prefixvm( M_prefix, "output-marker")))
{
    tic();
    M_thermo = thermo_type::New(M_prefixThermo,MixedPoissonPhysics::Heat);
    M_temperatureKey = M_thermo->potentialKey();
    M_heatFluxKey = M_thermo->fluxKey();
    M_electro = electro_type::New(M_prefixElectro,MixedPoissonPhysics::Electric);
    M_potentialKey = M_electro->potentialKey();
    M_currentKey = M_electro->fluxKey();
    toc("construct");
}

template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::init( mesh_ptrtype const& mesh )
{
    tic();
    if( mesh )
        M_mesh = mesh;
    else
        M_mesh = loadMesh( new mesh_type );
    M_thermo->init( M_mesh );
    M_electro->init( M_mesh );
    toc("init");
}

template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::solve()
{
    auto electroProp = M_electro->modelProperties();
    auto electroMat = electroProp.materials();
    auto thermoProp = M_thermo->modelProperties();
    auto thermoMat = M_thermo->modelProperties().materials();
    double V = 0, Vinit = 0, I = 0, Iinit = 0;
    bool conditionIbc = false;
    if( M_useContinuation )
    {
        try
        {
            auto bndCnd = electroProp.boundaryConditions2().at(M_potentialKey).at("Dirichlet").at(M_outputMarker);
            Vinit = std::stod(bndCnd.expression());
            if( Vinit == 0 )
                throw std::out_of_range("no dirichlet condition on potential on "+M_outputMarker);
            else
                V = Vinit/M_continuationSteps;

        }
        catch( const std::out_of_range& oor)
        {
            try
            {
                auto bndCnd = electroProp.boundaryConditions2().at(M_currentKey).at("Integral").at(M_outputMarker);
                Iinit = std::stod(bndCnd.expression());
                if( Iinit != 0 )
                {
                    conditionIbc = true;
                    I = Iinit/M_continuationSteps;
                }
                else
                    throw std::out_of_range("no integral condition on flux on "+M_outputMarker);
            }
            catch( const std::out_of_range& oor)
            {
                Feel::cout << tc::red << "marker for continuation " << M_outputMarker
                           << " not found in Dirichlet or Ibc. Continuation not used !"
                           << tc::reset << std::endl;
                M_useContinuation = false;
            }
        }
    }

    M_joule = M_electro->potentialSpace()->element();
    M_k = M_electro->potentialSpace()->element();
    M_sigma = M_electro->potentialSpace()->element();

    tic();
    M_electro->assembleCstPart();
    toc("assembleCstElectro");

    // first iteration
#ifndef USE_SAME_MAT
    Feel::cout << tc::red << "WARNING copy cst part called !!!" << tc::reset << std::endl;
    M_electro->copyCstPart();
#endif
    if( boption("thermoelectric.load-initial-guess") )
    {
        M_initGuessSpace = std::make_shared<init_guess_space_type>(M_mesh);
        M_initGuess = M_initGuessSpace->element();
        M_initGuess.load(_path=soption("thermoelectric.load-initial-guess.path"));
        for( auto const& pairMat : electroMat )
        {
            std::string marker = pairMat.first;
            auto mat = pairMat.second;
            double alpha = mat.getDouble("alpha");
            double T0 = mat.getDouble("T0");
            double sigma0 = mat.getDouble(M_sigmaKey);
            auto sigma = mat.getScalar(M_sigmaKeyNL,
                                       {"T"}, {idv(M_initGuess)},
                                       {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
            M_electro->updateConductivityTerm( sigma, marker);
            M_sigma += vf::project( _space=M_electro->potentialSpace(),
                                       _range=markedelements(M_mesh,marker),
                                       _expr=sigma );
        }
    }
    else
    {
        M_electro->updateConductivityTerm( false );
    }
    M_electro->assembleRhsBoundaryCond();
    if( M_useContinuation )
    {
        if( conditionIbc )
        {
            M_electro->assembleRhsIBC(0, M_outputMarker, I-Iinit);
        }
        else
        {
            auto g = expr(std::to_string(V-Vinit));
            M_electro->assembleRhsDirichlet(g,M_outputMarker);
        }
    }
    M_electro->assemblePotentialRHS( cst(0.), "");
    tic();
    M_electro->solve();
    toc("solve electro");
    M_potential = M_electro->potentialField();
    M_current = M_electro->fluxField();


    tic();
    M_thermo->assembleCstPart();
    toc("assembleCstThermo");

#ifndef USE_SAME_MAT
    Feel::cout << tc::red << "WARNING copy cst part called !!!" << tc::reset << std::endl;
    M_thermo->copyCstPart();
#endif
    for( auto const& pairMat : thermoMat )
    {
        std::string marker = pairMat.first;
        auto mat = pairMat.second;
        double alpha = mat.getDouble("alpha");
        double T0 = mat.getDouble("T0");
        auto sigma0 = mat.getDouble(M_sigmaKey);
        double k0 = 0, L = 0;
        if( mat.hasProperty("Lorentz") )
            L = mat.getDouble("Lorentz");
        if( mat.hasProperty(M_kKey) )
            k0 = mat.getDouble(M_kKey);
        auto k = mat.getScalar(M_kKeyNL,
                               {{"T",T0},{"k0",k0},{"T0",T0},{"alpha",alpha},{"Lorentz",L},{"sigma0",sigma0}});
        M_thermo->updateConductivityTerm( k, marker);
        M_k += vf::project( _space=M_electro->potentialSpace(),
                            _range=markedelements(M_mesh,marker),
                            _expr=k );
    }
    M_thermo->assembleRhsBoundaryCond();
    for( auto const& pairMat : electroMat )
    {
        std::string marker = pairMat.first;
        auto mat = pairMat.second;
        if( boption("thermoelectric.load-initial-guess") )
        {
            double alpha = mat.getDouble("alpha");
            double T0 = mat.getDouble("T0");
            double sigma0 = mat.getDouble(M_sigmaKey);
            auto sigma = mat.getScalar(M_sigmaKeyNL,
                                       {"T"}, {idv(M_initGuess)},
                                       {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
            auto rhs = inner(idv(M_current),idv(M_current))/sigma;
            M_joule += vf::project( _space=M_electro->potentialSpace(),
                                  _range=markedelements(M_mesh,marker),
                                  _expr=rhs );
            M_thermo->assemblePotentialRHS( rhs, marker);
        }
        else
        {
            auto sigma = mat.getDouble(M_sigmaKey);
            auto rhs = inner(idv(M_current),idv(M_current))/sigma;
            M_joule += vf::project( _space=M_electro->potentialSpace(),
                                    _range=markedelements(M_mesh,marker),
                                    _expr=rhs );
            M_thermo->assemblePotentialRHS( rhs, marker);
        }
    }
    tic();
    M_thermo->solve();
    toc("solve thermo");
    M_temperature = M_thermo->potentialField();
    M_tempflux = M_thermo->fluxField();

    potential_type oldPotential = M_electro->potentialSpace()->element();
    current_type oldCurrent = M_electro->fluxSpace()->element();
    temp_type oldTemperature = M_thermo->potentialSpace()->element();
    tempflux_type oldTempflux = M_thermo->fluxSpace()->element();
    // auto rangeV = M_electro->potentialSpace()->dof()->meshSupport()->rangeElements();
    // auto rangeT = M_thermo->potentialSpace()->dof()->meshSupport()->rangeElements();
    // double normV = normL2( rangeV, idv(M_potential) );
    // double normT = normL2( rangeT, idv(M_temperature) );
    // double incrV = normL2( rangeV, idv(M_potential) - idv(oldPotential) );
    // double incrT = normL2( rangeT, idv(M_temperature) - idv(oldTemperature) );
    double normV = normL2( elements(M_mesh), idv(M_potential) );
    double normT = normL2( elements(M_mesh), idv(M_temperature) );
    double incrV = normL2( elements(M_mesh), idv(M_potential) - idv(oldPotential) );
    double incrT = normL2( elements(M_mesh), idv(M_temperature) - idv(oldTemperature) );
    double relV = incrV/normV;
    double relT = incrT/normT;

    double tol = doption( prefixvm( M_prefix, "tolerance") );
    int itmax = ioption( prefixvm( M_prefix, "itmax") );


    int steps = M_useContinuation ? M_continuationSteps : 1;
    for( int c = 1; c <= steps; ++c)
    {
        if( M_useContinuation )
        {
            if( conditionIbc )
            {
                I = c*Iinit/M_continuationSteps;
                Feel::cout << "Start Picard with I =  " << I << " (" << itmax
                           << " iterations max and tolerance = " << tol << ")" << std::endl;
            }
            else
            {
                V = c*Vinit/M_continuationSteps;
                Feel::cout << "Start Picard with V =  " << V << " (" << itmax
                           << " iterations max and tolerance = " << tol << ")" << std::endl;
            }
        }
        else
            Feel::cout << "Start Picard  (" << itmax
                       << " iterations max and tolerance = " << tol << ")" << std::endl;

        int i = 0;
        relV = 1;
        relT = 1;

        // Picard loop
        while( i++ < itmax && ( relV > tol || relT > tol ) )
        {
            tic();

            oldPotential = M_potential;
            oldCurrent = M_current;
            oldTemperature = M_temperature;
            oldTempflux = M_tempflux;
            M_sigma.zero();
            M_k.zero();
            M_joule.zero();

            tic();
#ifdef USE_SAME_MAT
            M_electro->setMatricesAndVectorToZero();
            M_electro->assembleCstPart();
#else
            M_electro->copyCstPart();
#endif
            for( auto const& pairMat : electroMat )
            {
                std::string marker = pairMat.first;
                auto mat = pairMat.second;
                double alpha = mat.getDouble("alpha");
                double T0 = mat.getDouble("T0");
                double sigma0 = mat.getDouble(M_sigmaKey);
                auto sigma = mat.getScalar(M_sigmaKeyNL,
                                           {"T"}, {idv(M_temperature)},
                                           {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
                M_electro->updateConductivityTerm( sigma, marker);
                M_sigma += vf::project( _space=M_electro->potentialSpace(),
                                        _range=markedelements(M_mesh,marker),
                                        _expr=sigma );
            }
            M_electro->assembleRhsBoundaryCond();
            if( M_useContinuation )
            {
                if( conditionIbc )
                {
                    M_electro->assembleRhsIBC(0, M_outputMarker, I-Iinit);
                }
                else
                {
                    auto g = expr(std::to_string(V-Vinit));
                    M_electro->assembleRhsDirichlet(g,M_outputMarker);
                }
            }
            M_electro->assemblePotentialRHS( cst(0.), "");
            toc("assembleElectro");
            tic();
            M_electro->solve();
            toc("solve electro");
            M_potential = M_electro->potentialField();
            M_current = M_electro->fluxField();

            tic();
#ifdef USE_SAME_MAT
            M_thermo->setMatricesAndVectorToZero();
            M_thermo->assembleCstPart();
#else
            M_thermo->copyCstPart();
#endif
            for( auto const& pairMat : thermoMat )
            {
                std::string marker = pairMat.first;
                auto mat = pairMat.second;
                double alpha = mat.getDouble("alpha");
                double T0 = mat.getDouble("T0");
                auto sigma0 = mat.getDouble(M_sigmaKey);
                double k0 = 0, L = 0;
                if( mat.hasProperty("Lorentz") )
                    L = mat.getDouble("Lorentz");
                if( mat.hasProperty(M_kKey) )
                    k0 = mat.getDouble(M_kKey);
                auto k = mat.getScalar(M_kKeyNL,
                                       {"T"}, {idv(M_temperature)},
                                       {{"k0",k0},{"T0",T0},{"alpha",alpha},
                                        {"Lorentz",L},{"sigma0",sigma0}});
                M_thermo->updateConductivityTerm( k, marker);
                M_k += vf::project( _space=M_electro->potentialSpace(),
                                    _range=markedelements(M_mesh,marker),
                                    _expr=k );
            }
            M_thermo->assembleRhsBoundaryCond();
            for( auto const& pairMat : electroMat )
            {
                std::string marker = pairMat.first;
                auto mat = pairMat.second;
                double alpha = mat.getDouble("alpha");
                double T0 = mat.getDouble("T0");
                double sigma0 = mat.getDouble(M_sigmaKey);
                auto sigma = mat.getScalar(M_sigmaKeyNL,
                                           {"T"}, {idv(M_temperature)},
                                           {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
                auto rhs = inner(idv(M_current),idv(M_current))/sigma;
                M_joule += vf::project( _space=M_electro->potentialSpace(),
                                        _range=markedelements(M_mesh,marker),
                                        _expr=rhs);
                M_thermo->assemblePotentialRHS( rhs, marker);
            }
            toc("assembleThermo");
            tic();
            M_thermo->solve();
            toc("solve thermo");
            M_temperature = M_thermo->potentialField();
            M_tempflux = M_thermo->fluxField();

            incrV = normL2( elements(M_mesh), idv(M_potential) - idv(oldPotential) );
            incrT = normL2( elements(M_mesh), idv(M_temperature) - idv(oldTemperature) );
            normV = normL2( elements(M_mesh), idv(oldPotential) );
            normT = normL2( elements(M_mesh), idv(oldTemperature) );
            relV = incrV/normV;
            relT = incrT/normT;

            Feel::cout << "iteration " << i
                       << "\n\tincrement on V = " << relV << " (absolute = " << incrV << ")"
                       << "\n\tincrement on T = " << relT << " (absolute = " << incrT << ")" << std::endl;

            toc("loop");
        } // Picard loop
    }
}

template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::exportResults()
{
    M_electro->exportResults();
    M_thermo->exportResults();
    auto e = exporter( M_mesh );
    e->add("potential", M_potential);
    e->add("current", M_current);
    e->add("temperature", M_temperature);
    e->add("tempflux", M_tempflux);
    e->add("joule", M_joule);
    e->add("sigma", M_sigma);
    e->add("k", M_k);
    e->save();

    double meanT = mean(_range=elements(M_mesh), _expr=idv(M_temperature))(0,0);
    auto mima = minmax(_range=elements(M_mesh), _pset=_Q<>(), _expr=idv(M_temperature));
    double miT = mima.min();
    double maT = mima.max();
    auto power = integrate(_range=elements(M_mesh), _expr=idv(M_joule)).evaluate()(0,0);
    auto int0 = integrate(_range=markedfaces(M_mesh,M_outputMarker),_expr=inner(idv(M_current),N())).evaluate()(0,0);
    auto int1 = integrate(_range=markedfaces(M_mesh,M_inputMarker),_expr=inner(idv(M_current),N())).evaluate()(0,0);
    double meanP = mean(_range=markedfaces(M_mesh,M_outputMarker), _expr=idv(M_potential))(0,0);
    auto mimaP = minmax(_range=markedfaces(M_mesh,M_outputMarker), _pset=_Q<>(), _expr=idv(M_potential));
    double miP = mimaP.min();
    double maP = mimaP.max();

    Feel::cout << std::setw(25) << "tempMax";
    Feel::cout << std::setw(25) << "tempMoy";
    Feel::cout << std::setw(25) << "tempMin";
    Feel::cout << std::setw(25) << "intensite0";
    Feel::cout << std::setw(25) << "intensite1";
    Feel::cout << std::setw(25) << "power";
    Feel::cout << std::setw(25) << "potMax";
    Feel::cout << std::setw(25) << "potMoy";
    Feel::cout << std::setw(25) << "potMin";
    Feel::cout << std::endl;
    Feel::cout << std::setw(25) << maT;
    Feel::cout << std::setw(25) << meanT;
    Feel::cout << std::setw(25) << miT;
    Feel::cout << std::setw(25) << int0;
    Feel::cout << std::setw(25) << int1;
    Feel::cout << std::setw(25) << power;
    Feel::cout << std::setw(25) << maP;
    Feel::cout << std::setw(25) << meanP;
    Feel::cout << std::setw(25) << miP;
    Feel::cout << std::endl;
    if( Environment::isMasterRank() )
    {
        std::ofstream file ( "measures.csv" );
        if( file )
        {
            file << std::setw(25) << "tempMax";
            file << std::setw(25) << "tempMoy";
            file << std::setw(25) << "tempMin";
            file << std::setw(25) << "intensite0";
            file << std::setw(25) << "intensite1";
            file << std::setw(25) << "power";
            file << std::setw(25) << "potMax";
            file << std::setw(25) << "potMoy";
            file << std::setw(25) << "potMin";
            file << std::endl;
            file << std::setw(25) << maT;
            file << std::setw(25) << meanT;
            file << std::setw(25) << miT;
            file << std::setw(25) << int0;
            file << std::setw(25) << int1;
            file << std::setw(25) << power;
            file << std::setw(25) << maP;
            file << std::setw(25) << meanP;
            file << std::setw(25) << miP;
            file << std::endl;
            file.close();
        }
    }
}
